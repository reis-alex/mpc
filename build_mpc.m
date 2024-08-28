%% Build an MPC controller using CasADi
% Output: CasADi solver
% Input: user-defined options
%
%   opt.n_states      : number of states
%   opt.n_controls    : number of control inputs
%   opt.N             : lenght of prediction horizon
% 
% Modeling:
% 
%   opt.model.type  : either 'linear' or 'nonlinear';
%   opt.model.A     : for 'model.type' = 'linear'
%   opt.model.B     : for 'model.type' = 'linear
%   opt.model.fun   : if 'model.type' = 'nonlinear', it must be a function handle with independent
%                     variables ordered as @(x,u)
%
% Cost functions (stage and terminal):
% 
% opt.stage.cost_function.function: if a given stage cost function is given, it must be
%                                   a function handle with variables named and ordered as (x,u,extra),
%                                   where extra can be extra parameters used in the function
% 
% If no function is given, then V = x'*Q*x + u'*R*u is used, with
% 
% opt.stage.cost_function.Q/R = Q/R numerical matrices
% opt.stage.cost_function.extra_parameters      : must contain the same variables
%                                                 as used in "extra" in opt.stage.cost_function.function.
% opt.terminal.cost_function.function           : if a given terminal cost function is given, it must be a 
%                                                 function handle ordered as @(x,extra), where "extra" is 
%                                                 an be extra parameters used in the function
% opt.terminal.cost_function.extra_parameters   : must contain the same variables
%                                                 as used in "extra" in opt.terminal.cost_function.function
% opt.terminal.set.A/b                          : describes the (polyhedral) terminal set (Ax(N)<=b)
% 
% Parameters:
% 
% opt.extra_parameters.constrained              : list of parameters that are constrained
% opt.extra_parameters.constraints.lower/upper  : lower/upper bound of such parameters
% opt.extra_parameters.input                    : parameters that are input to the optimization
% opt.constraints.states.lower/upper            : if constraints for states are given
%                                                 point-wise, upper/lower
% opt.constraints.control.upper/lower           : if constraints for controls are
%                                                 given point-wisely, upper/lower
% opt.constraints.polyhedral.A/b                : describes the (polyhedral) constraint set (Ax(k)<=b)
%
% Solver selection:
%
% opt.solver = 'ipopt' or 'qpoases'             : for general NLP and QP, respectively

%%
function [solver,args] = build_mpc(opt)
import casadi.* 

%generate states and control vectors
states = [];
for i = 1:opt.n_states
   states = [states; SX.sym(['x' int2str(i)])];
end

controls = [];
for i = 1:opt.n_controls
    controls = [controls; SX.sym(['u' int2str(i)])];
end

% generate model either linear or nonlinear
switch opt.model.type
    case 'linear'
        model = opt.model.A*states + opt.model.B*controls;
    case 'nonlinear'
        model = opt.model.fun(states);
end

% if continuous model, select how to integration
if isfield(opt,'continuous_model')
    % set integration scheme
    switch opt.continuous_model.integration
        case 'euler'
            integ = @(x,h,xplus) x + h*xplus;
        case 'RK4'
            % to be finished
    end
else
    % if not continuous model, then just use difference equation
    integ = @(x,h,xplus) xplus;
    opt.dt = 0;
end

% model function f(x,u), and optimization vectors
f = Function('f',{states,controls},{model}); 

U           = SX.sym('U',opt.n_controls,opt.N);         % Decision variables (controls)
Init_states = SX.sym('Init',opt.n_states); 
X           = SX.sym('X',opt.n_states,(opt.N+1));       % A vector that represents the states over the optimization problem.


% Build prediction and cost function, initialize prediction with
% measurement
obj         = 0;
g           = [];
g           = [g; X(:,1)-Init_states];                 

% Determine the cost function to be used: simple quadratic or a custom
if isfield(opt.stage.cost_function,'function')
    stagecost_fun    = opt.stage.cost_function.function;
else
    stagecost_fun    = @(x,u,extra) x'*opt.cost_function.Q*x + u'*opt.cost_function.R*u + 0*extra;
end

% if the cost function requires extra variables (as optimization variables or not)
if isfield(opt.stage.cost_function,'extra_parameters')
    extra_parameters = opt.stage.cost_function.extra_parameters;
else
    extra_parameters = 0;
end

% prediction loop
for k = 1:opt.N
    obj             = obj + stagecost_fun(X(:,k),U(:,k),extra_parameters);
    state_next      = X(:,k+1);
    f_value         = f(X(:,k),U(:,k));
    st_next_euler   = integ(X(:,k),opt.dt,f_value);
    g               = [g; state_next-st_next_euler];    
end

% if constraints in states are polyhedral
if isfield(opt.constraints,'polyhedral')
    for i = 1:opt.N
        g = [g; opt.constraints.polyhedral.A*X(:,i)-opt.constraints.polyhedral.b];
    end
end

% if terminal costs and constraints
if isfield(opt,'terminal')
    obj = obj + opt.terminal.cost_function.function(X(:,end),opt.terminal.cost_function.extra_parameters);
        if isfield(opt.terminal,'extra_parameters')
                g   = [g; opt.terminal.set.A*vertcat(X(:,end),opt.terminal.extra_parameters)-opt.terminal.set.b];
        else
                g   = [g; opt.terminal.set.A*X(:,end)-opt.terminal.set.b];
        end
end

% make the decision variable one column vector
OPT_variables   = [reshape(X(:,1:end-1),opt.n_states*opt.N,1);
                    reshape(U,opt.n_controls*opt.N,1);];
                
% add terminal constraint variables to the list
if isfield(opt,'terminal')
    OPT_variables = [OPT_variables;
                     reshape(X(:,end),opt.n_states,1)];
end

% add extra variables to the list
if isfield(opt,'extra_parameters') && isfield(opt.extra_parameters,'constrained')
    OPT_variables   = [OPT_variables; 
                       opt.extra_parameters.constrained];
end

% define external parameters and problem structure

if isfield(opt,'extra_parameters') && isfield(opt.extra_parameters,'input')
    Param = [Init_states; opt.extra_parameters.input];
else
    Param = [Init_states];
end
OPC   = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', Param);

% define equality constraints 
args = struct;

% equality for x(k+1)-x(k)
args.lbg(1:opt.n_states*(opt.N+1)) = 0;
args.ubg(1:opt.n_states*(opt.N+1)) = 0;

% if constraints on states are polyhedral
if isfield(opt.constraints,'polyhedral')
    args.lbg(length(args.ubg)+1:length(args.ubg)+opt.N*length(opt.constraints.polyhedral.b)) = -inf;
    args.ubg(length(args.ubg)+1:length(args.ubg)+opt.N*length(opt.constraints.polyhedral.b)) = 0;
end

% if terminal constraint, inequality for Xf.A*(x(N),xa,ua)<=Xf.b
if isfield(opt,'terminal')
    args.lbg(length(args.ubg)+1:length(args.ubg)+length(opt.terminal.set.b)) = -inf;
    args.ubg(length(args.ubg)+1:length(args.ubg)+length(opt.terminal.set.b)) = 0; 
end


% inequality constraints
% bounds for the states variables
for k = 1:opt.n_states
    if isfield(opt.constraints,'polyhedral')
        args.lbx(k:opt.n_states:opt.n_states*(opt.N),1) = -inf;
        args.ubx(k:opt.n_states:opt.n_states*(opt.N),1) = +inf;
    else
        args.lbx(k:opt.n_states:opt.n_states*(opt.N),1) = opt.constraints.states.lower(k);
        args.ubx(k:opt.n_states:opt.n_states*(opt.N),1) = opt.constraints.states.upper(k);
    end
end

% bounds for variables (controls)
for k = 1:opt.n_controls
    args.lbx(opt.n_states*opt.N+k:opt.n_controls:opt.n_states*opt.N+opt.n_controls*opt.N,1) = opt.constraints.control.lower(k);    
    args.ubx(opt.n_states*opt.N+k:opt.n_controls:opt.n_states*opt.N+opt.n_controls*opt.N,1) = opt.constraints.control.upper(k);     
end

% bounds for variables (X(N))
if isfield(opt,'terminal')
    if isfield(opt.constraints,'polyhedral')
        args.ubx(length(args.ubx)+1:length(args.ubx)+opt.n_states) = +inf;
        args.lbx(length(args.lbx)+1:length(args.lbx)+opt.n_states) = -inf;
    else
        args.ubx(length(args.ubx)+1:length(args.ubx)+opt.n_states) = opt.constraints.states.upper;
        args.lbx(length(args.lbx)+1:length(args.lbx)+opt.n_states) = opt.constraints.states.lower;
    end
    
end

% bounds for extra variables
if isfield(opt,'extra_parameters') && isfield(opt.extra_parameters,'constrained')
    args.ubx(length(args.ubx)+1:length(args.ubx)+length(opt.extra_parameters.constrained)) = opt.extra_parameters.constraints.upper;
    args.lbx(length(args.lbx)+1:length(args.lbx)+length(opt.extra_parameters.constrained)) = opt.extra_parameters.constraints.lower;
end

if (length(OPT_variables)~=length(args.lbx))
    error('MPC error: Number of variables and respective bounds are different')
end


% generate solver
opts                        = struct;
switch opt.solver
    case 'ipopt'
        opts.ipopt.max_iter         = 200;
        opts.ipopt.print_level      = 0;
        opts.print_time             = 0;
        opts.ipopt.acceptable_tol   = 1e-8;
        opts.ipopt.acceptable_obj_change_tol = 1e-8;
        solver = nlpsol('solver', 'ipopt', OPC,opts);
    case 'qpoases'
        options.terminationTolerance = 1e-5;
        options.boundTolerance = 1e-5;
        options.printLevel = 'none';
        options.error_on_fail = 0;
        solver = qpsol('solver','qpoases', OPC,options);
end
