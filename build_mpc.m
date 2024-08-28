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

% Parameters
Init_states = SX.sym('Init',opt.n_states);

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
Ref_states  = SX.sym('Ref',opt.n_states);               % parameters (which include the initial state and the reference state)
X           = SX.sym('X',opt.n_states,(opt.N+1));       % A vector that represents the states over the optimization problem.


% Build prediction and cost function
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

% if terminal costs and constraints
if isfield(opt,'terminal')
    obj = obj + opt.terminal.cost_function.function(X(:,end),opt.terminal.cost_function.extra_parameters);
    if isfield(opt.extra_parameters,'constrained')
            g   = [g; opt.terminal.set.A*vertcat(X(:,end),opt.extra_parameters.constrained)-opt.terminal.set.b];
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
if isfield(opt.extra_parameters,'constrained')
    OPT_variables   = [OPT_variables; 
                       opt.extra_parameters.constrained];
end

% define external parameters and solver 

if isfield(opt.extra_parameters,'input')
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

% if terminal constraint, inequality for Xf.A*(x(N),xa,ua)<=Xf.b
if isfield(opt,'terminal')
    args.lbg(length(args.ubg)+1:length(args.ubg)+length(opt.terminal.set.b)) = -inf;
    args.ubg(length(args.ubg)+1:length(args.ubg)+length(opt.terminal.set.b)) = 0; 
end

% inequality constraints
% bounds for the states variables
for k = 1:opt.n_states
    args.lbx(k:opt.n_states:opt.n_states*(opt.N),1) = opt.constraints.states.lower(k); 
    args.ubx(k:opt.n_states:opt.n_states*(opt.N),1) = opt.constraints.states.upper(k);      
end

% bounds for variables (controls)
for k = 1:opt.n_controls
    args.lbx(opt.n_states*opt.N+k:opt.n_controls:opt.n_states*opt.N+opt.n_controls*opt.N,1) = opt.constraints.control.lower(k);    
    args.ubx(opt.n_states*opt.N+k:opt.n_controls:opt.n_states*opt.N+opt.n_controls*opt.N,1) = opt.constraints.control.upper(k);     
end

% bounds for variables (X(N))
if isfield(opt,'terminal')
    args.ubx(length(args.ubx)+1:length(args.ubx)+opt.n_states) = opt.constraints.states.upper;
    args.lbx(length(args.lbx)+1:length(args.lbx)+opt.n_states) = opt.constraints.states.lower;
end

% bounds for extra variables
if isfield(opt.extra_parameters,'constrained')
    args.ubx(length(args.ubx)+1:length(args.ubx)+length(opt.extra_parameters.constrained)) = opt.extra_parameters.constraints.upper;
    args.lbx(length(args.lbx)+1:length(args.lbx)+length(opt.extra_parameters.constrained)) = opt.extra_parameters.constraints.lower;
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
