%% Build an MPC controller using CasADi
% Output: CasADi solver
% Input: user-defined options
%
%   opt.n_states      : number of states
%   opt.n_controls    : number of control inputs
%   opt.N             : lenght of prediction horizon
% 
% * Modeling:
% 
%   opt.model.type     : either 'linear' or 'nonlinear';
%   opt.model.A        : for 'model.type' = 'linear'
%   opt.model.B        : for 'model.type' = 'linear
%   opt.model.function : if 'model.type' = 'nonlinear', it must be a function
%                        handle with independent variables ordered as @(x,u)
%
% * Cost functions (stage and terminal): all costs are grouped within
% opt.costs
%
%   opt.costs.stage.function  : if a given stage cost function is given, it 
%                               must be a function handle with variables named 
%                               and ordered as (x,u,extra), where extra can be 
%                               extra parameters used in the function.
% 
% If no specific function is given, then V = x'*Q*x + u'*R*u is used, with
% opt.stage.cost_function.Q/R = Q/R being numerical matrices of proper
% dimensions
% 
%   opt.costs.stage.parameters    : must contain the same variables
%                                   as used in "extra" in opt.costs.stage.function.
% 
%   opt.costs.terminal.function   : if a given terminal cost function is given, it must be a 
%                                   function handle ordered as @(x,extra), where "extra" are 
%                                   extra parameters used in the function
%   opt.costs.terminal.parameters : must contain the same variables
%                                   as used in "extra" in opt.terminal.cost_function.function
% 
% * Constraints : all are grouped within opt.constraints 
% 
%   opt.constraints.states        : variable-wise bounds. Must contain fields
%                                   "upper" and "lower" with proper dimensions.
%   opt.constraints.polyhedral    : describes the (polyhedral) state
%                                   constraint set. It must contain fields set.A and set.b.
% 
% If fields "states" or "polyhedral" are given, the states are considered
% unconstrained.
% 
%   opt.constraints.control       : variable-wise bounds. Must contain fields
%                                   "upper" and "lower" with proper dimensions.
% 
%   opt.constraints.terminal.set  : describes the (polyhedral) terminal set
%                                   (Ax(N)<=b). It must contain fields set.A
%                                   and set.b.
% 
%	opt.constraints.terminal.parameters : list of all variables that are
%                                         terminally constrained.
% 
%   opt.constraints.parameters.variables : list of all extra variables that
%                                          are constrained somehow.
% 
% * Parameters:
% 
%   opt.input.vector : list of parameters taken as input to the optimization
%                       problem.
% 
% * Solver selection:
%
%   opt.solver = 'ipopt' or 'qpoases'    : for general NLP and QP, respectively

%%
function [solver,args] = build_mpc(opt)
import casadi.* 

% check whether general constraints require inputs (vector or matrices)
vararg_i = 1;
if isfield(opt.input,'general_constraints') && isfield(opt.input.general_constraints,'matrix')
    input_matrix = SX.sym('input_matrix',opt.input.general_constraints.matrix.dim(1),opt.input.general_constraints.matrix.dim(2));
    vararg{vararg_i} = input_matrix; %opt.constraints.general.matrix = input_matrix;
    vararg_i = vararg_i + 1;
end

if isfield(opt.input,'general_constraints') && isfield(opt.input.general_constraints,'vector')
    input_vector = SX.sym('input_vector',opt.input.general_constraints.vector.dim);
    vararg{vararg_i} = input_vector;%opt.constraints.general.vector = input_vector;
end

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
        model = opt.model.function(states);
end

%--- something that checks opt.n_states and opt.n_controls and the size
% of A,B or opt.model.function 

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
if isfield(opt.costs.stage,'function')
    stagecost_fun    = opt.costs.stage.function;
else
    stagecost_fun    = @(x,u,extra) x'*opt.costs.stage.Q*x + u'*opt.costs.stage.R*u + 0*extra;
end

% if the cost function requires extra variables (as optimization variables or not)
if isfield(opt.costs.stage,'parameters')
    stage_parameters = opt.costs.stage.parameters;
else
    stage_parameters = 0;
end

% prediction loop
for k = 1:opt.N
    obj             = obj + stagecost_fun(X(:,k),U(:,k),stage_parameters);
    state_next      = X(:,k+1);
    f_value         = f(X(:,k),U(:,k));
    st_next_euler   = integ(X(:,k),opt.dt,f_value);
    g               = [g; state_next-st_next_euler];    
end

% if constraints in states are polyhedral
if isfield(opt,'constraints') && isfield(opt.constraints,'polyhedral')
    for i = 1:opt.N
        g = [g; opt.constraints.polyhedral.A*X(:,i)-opt.constraints.polyhedral.b];
    end
end

% if terminal costs and constraints A CORRIGIR
if isfield(opt,'constraints') && isfield(opt.constraints,'terminal') && isfield(opt.constraints.terminal,'parameters')
    g   = [g; opt.constraints.terminal.set.A*vertcat(X(:,end),opt.constraints.terminal.parameters)-opt.constraints.terminal.set.b];
    obj = obj + opt.costs.terminal.function(X(:,end),opt.costs.terminal.parameters);
else
    g   = [g; opt.constraints.terminal.set.A*X(:,end)-opt.constraints.terminal.set.b];
    obj = obj + opt.costs.terminal.function(X(:,end));
end


% if there are general constraints
if isfield(opt,'constraints') && isfield(opt.constraints,'general')
    for i = 1:opt.N
        g = [g; opt.constraints.general.function(X(:,i),vararg{:})];
    end
end
%% Define vector of decision variables
% make the decision variable one column vector
OPT_variables   = [reshape(X(:,1:end-1),opt.n_states*opt.N,1);
                    reshape(U,opt.n_controls*opt.N,1);];
                
% add terminal constraint variables to the list
if isfield(opt,'constraints') && isfield(opt.constraints,'terminal')
    OPT_variables = [OPT_variables;
                     reshape(X(:,end),opt.n_states,1)];
end

% add extra variables to the list
if isfield(opt,'constraints') && isfield(opt.constraints,'parameters')
    OPT_variables   = [OPT_variables; 
                       opt.constraints.parameters.variables];
end

%% define external parameters and problem structure

Param = [Init_states];
if isfield(opt,'input')
    if isfield(opt.input,'vector')
        Param = [Param; opt.input.vector];
    end
    
    if isfield(opt.input.general_constraints,'matrix')
        aux = [];
        for jj = 1:opt.input.general_constraints.matrix.dim(2)
            aux = [aux; input_matrix(:,jj)];
        end
        Param = [Param; aux];
    end
    
    if isfield(opt.input.general_constraints,'vector')
        Param = [Param; input_vector];
    end
end

OPC   = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', Param);

%% Constraints
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
if isfield(opt.constraints,'terminal') && isfield(opt.constraints.terminal,'set') 
    args.lbg(length(args.ubg)+1:length(args.ubg)+length(opt.constraints.terminal.set.b)) = -inf;
    args.ubg(length(args.ubg)+1:length(args.ubg)+length(opt.constraints.terminal.set.b)) = 0; 
end

% if general constraints
if isfield(opt.constraints,'general')
    args.lbg(length(args.ubg)+1:length(args.ubg)+opt.N*opt.constraints.general.dim) = -inf;
    args.ubg(length(args.ubg)+1:length(args.ubg)+opt.N*opt.constraints.general.dim) = 0;
end

%% inequality constraints
% bounds for the states variables
for k = 1:opt.n_states
    if isfield(opt.constraints,'polyhedral')
        args.lbx(k:opt.n_states:opt.n_states*(opt.N),1) = -inf;
        args.ubx(k:opt.n_states:opt.n_states*(opt.N),1) = +inf;
    else
        if isfield(opt.constraints,'states')
            args.lbx(k:opt.n_states:opt.n_states*(opt.N),1) = opt.constraints.states.lower(k);
            args.ubx(k:opt.n_states:opt.n_states*(opt.N),1) = opt.constraints.states.upper(k);
        else
            args.lbx(k:opt.n_states:opt.n_states*(opt.N),1) = -inf;
            args.ubx(k:opt.n_states:opt.n_states*(opt.N),1) = +inf;
        end
    end
end


% bounds for variables (controls)
if isfield(opt.constraints,'control')
    for k = 1:opt.n_controls
        args.lbx(opt.n_states*opt.N+k:opt.n_controls:opt.n_states*opt.N+opt.n_controls*opt.N,1) = opt.constraints.control.lower(k);
        args.ubx(opt.n_states*opt.N+k:opt.n_controls:opt.n_states*opt.N+opt.n_controls*opt.N,1) = opt.constraints.control.upper(k);
    end
else
    for k = 1:opt.n_controls
        args.lbx(opt.n_states*opt.N+k:opt.n_controls:opt.n_states*opt.N+opt.n_controls*opt.N,1) = -inf;
        args.ubx(opt.n_states*opt.N+k:opt.n_controls:opt.n_states*opt.N+opt.n_controls*opt.N,1) = +inf;
    end
end

% bounds for variables (X(N))
if isfield(opt,'constraints') && isfield(opt.constraints,'terminal')
    if isfield(opt.constraints.terminal,'set')
        args.ubx(length(args.ubx)+1:length(args.ubx)+opt.n_states) = +inf;
        args.lbx(length(args.lbx)+1:length(args.lbx)+opt.n_states) = -inf;
    else
        args.ubx(length(args.ubx)+1:length(args.ubx)+opt.n_states) = opt.constraints.states.upper;
        args.lbx(length(args.lbx)+1:length(args.lbx)+opt.n_states) = opt.constraints.states.lower;
    end
    
end

% bounds for extra variables
if isfield(opt,'constraints') && isfield(opt.constraints,'parameters')
    if isfield(opt.constraints.parameters,'lower')
        args.ubx(length(args.ubx)+1:length(args.ubx)+length(opt.constraints.parameters.variables)) = opt.constraints.parameters.upper;
        args.lbx(length(args.lbx)+1:length(args.lbx)+length(opt.constraints.parameters.variables)) = opt.constraints.parameters.lower;
    else
        args.ubx(length(args.ubx)+1:length(args.ubx)+length(opt.constraints.parameters.variables)) = +inf;
        args.lbx(length(args.lbx)+1:length(args.lbx)+length(opt.constraints.parameters.variables)) = -inf;
    end
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
        solver = nlpsol('solver', 'ipopt',OPC,opts);
    case 'qpoases'
        options.terminationTolerance = 1e-5;
        options.boundTolerance = 1e-5;
        options.printLevel = 'none';
        options.error_on_fail = 0;
        solver = qpsol('solver','qpoases',OPC,options);
end
