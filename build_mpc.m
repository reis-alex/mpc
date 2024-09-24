%% Build an MPC controller using CasADi

function [solver,args] = build_mpc(opt)
import casadi.* 
 

%% Gather all parameters (decision variables, inputs - everything that becomes SX.sym)

if isfield(opt,'parameters')
    parameters{1} = []; 
        for i = [find(opt.parameters.dim(:,2)==1)]'
            parameters{i} = SX.sym(opt.parameters.name{i},opt.parameters.dim(i,1),1);
        end

        for i = [find(opt.parameters.dim(:,2)~=1)]'
            parameters{i} = reshape(SX.sym(opt.parameters.name{i},opt.parameters.dim(i,1),opt.parameters.dim(i,2)),opt.parameters.dim(i,2)*opt.parameters.dim(i,1),1);
        end
else
    parameters{1} = [];
end

% sort all parameters and their destinations
% parameters for stage cost
if isfield(opt.costs.stage,'parameters')
    parameters_stc = [];
    for i = 1:length(opt.costs.stage.parameters)
        parameters_stc = [parameters_stc; parameters{find(strcmp(opt.costs.stage.parameters{i},opt.parameters.name))}];
    end
end

% parameters for terminal cost
if isfield(opt.costs,'terminal') && isfield(opt.costs.terminal,'parameters')
    parameters_trc = []; 
    for i = 1:length(opt.costs.terminal.parameters)
        parameters_trc = [parameters_trc; parameters{find(strcmp(opt.costs.terminal.parameters{i},opt.parameters.name))}];
    end
end

% parameters for terminal constraints
if isfield(opt.constraints,'terminal') && isfield(opt.constraints.terminal,'parameters')
    parameters_tc = []; 
    for i = 1:length(opt.constraints.terminal.parameters)
        parameters_tc = [parameters_tc; parameters{find(strcmp(opt.constraints.terminal.parameters{i},opt.parameters.name))}];
    end
end

% parameters for general constraints
if isfield(opt.constraints,'parameters')
    parameters_const = [];
    for i = 1:length(opt.constraints.parameters)
        parameters_const = [parameters_const; parameters{find(strcmp(opt.constraints.parameters{i},opt.parameters.name))}];
    end
end

% parameters for inputs
if isfield(opt,'input')
    parameters_input = [];
    for i = 1:length(opt.input.vector)
        parameters_input = [parameters_input; parameters{find(strcmp(opt.input.vector{i},opt.parameters.name))}];
    end
    if isfield(opt.input,'matrix')
        for i = 1:length(opt.input.matrix)
            parameters_input = [parameters_input; parameters{find(strcmp(opt.input.matrix{i},opt.parameters.name))}];
        end
    end
end

% check whether general constraints require inputs (vector or matrices)

if isfield(opt,'input')
vararg{1} = [] ;
vararg_sc{1} = [];    
% if there are any input to general constraints
vararg_i = 1;
    if  isfield(opt.input,'general_constraints') && isfield(opt.input.general_constraints,'matrix')
        input_matrix = SX.sym('input_matrix',opt.input.general_constraints.matrix.dim(1),opt.input.general_constraints.matrix.dim(2));
        vararg{vararg_i} = input_matrix;
        vararg_i = vararg_i + 1;
    end

    if  isfield(opt.input,'general_constraints') && isfield(opt.input.general_constraints,'vector')
        input_vector = SX.sym('input_vector',opt.input.general_constraints.vector.dim);
        vararg{vararg_i} = input_vector;
    end

% if there are any input to stage costs
vararg_i = 1;
    if isfield(opt.input,'stage_costs') && isfield(opt.input.stage_costs,'matrix')
        input_matrix = SX.sym('input_matrix',opt.input.stage_costs.matrix.dim(1),opt.input.stage_costs.matrix.dim(2));
        vararg_sc{vararg_i} = input_matrix; 
        vararg_i = vararg_i + 1;
    end
    
    if  isfield(opt.input,'stage_costs') && isfield(opt.input.stage_costs,'vector')
        input_vector = SX.sym('input_vector',opt.input.stage_costs.vector.dim);
        vararg_sc{vararg_i} = input_vector;
    end
    
end

if isfield(opt.constraints,'general')
    opt.constraints.general.dim  = length(opt.constraints.general.function(zeros(opt.n_states,1),vararg{:}));
end

%% Modeling: generate states and control vectors, model function, integration
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
        model = opt.model.function(states,controls);
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
            integ = @(x,h,varargin) x + (h/6)*(varargin{1}+2*varargin{2}+2*varargin{3}+varargin{4});
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


%% Build prediction and cost function, initialize prediction

obj         = 0;
g           = [];
g           = [g; X(:,1)-Init_states];                 

% Determine the cost function to be used: simple quadratic or a custom
if isfield(opt.costs.stage,'function')
    stagecost_fun    = opt.costs.stage.function;
else
    stagecost_fun    = @(x,u,extra) x'*opt.costs.stage.Q*x + u'*opt.costs.stage.R*u + 0*extra;
end


%% Prediction loop, integration schemes
% prediction loop
if isfield(opt,'continuous_model')
    switch opt.continuous_model.integration
        case 'RK4'
            for k = 1:opt.N
                obj             = obj + stagecost_fun(X(:,k),U(:,k),parameters_stc);
                k1 = f(X(:,k),U(:,k));
                k2 = f(X(:,k)+opt.dt/2*k1,U(:,k));
                k3 = f(X(:,k)+opt.dt/2*k2,U(:,k));
                k4 = f(X(:,k)+opt.dt*k3,U(:,k));
                st_next_RK4   = integ(X(:,k),opt.dt,k1,k2,k3,k4);
                g               = [g; X(:,k+1)-st_next_RK4]; 
            end
        case 'euler'
            for k = 1:opt.N
                obj             = obj + stagecost_fun(X(:,k),U(:,k),parameters_stc);
                f_value         = f(X(:,k),U(:,k));
                st_next_euler   = integ(X(:,k),opt.dt,f_value);
                g               = [g; X(:,k+1)-st_next_euler]; 
            end
    end
else % if the model is already discrete
    for k = 1:opt.N
        obj             = obj + stagecost_fun(X(:,k),U(:,k),parameters_stc);%,vararg_sc{:});
        f_value         = f(X(:,k),U(:,k));
        st_next         = integ(X(:,k),opt.dt,f_value);
        g               = [g; X(:,k+1)-st_next];
    end
end

%%
% if constraints in states are polyhedral
if isfield(opt,'constraints') && isfield(opt.constraints,'polyhedral')
    for i = 1:opt.N
        g = [g; opt.constraints.polyhedral.A*X(:,i)-opt.constraints.polyhedral.b];
    end
end

% if terminal costs and constraints
if isfield(opt,'constraints') && isfield(opt.constraints,'terminal') && isfield(opt.constraints.terminal,'parameters')
    g   = [g; opt.constraints.terminal.set.A*vertcat(X(:,end),parameters_tc)-opt.constraints.terminal.set.b];
    obj = obj + opt.costs.terminal.function(X(:,end),parameters_trc);
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
                       parameters_tc];
end

%% define external parameters and problem structure

Param = [Init_states];
if isfield(opt,'input')
    if isfield(opt.input,'vector')
        Param = [Param; parameters_input];
    end
% if there are any input to general constraints
    if isfield(opt.input,'general_constraints') && isfield(opt.input.general_constraints,'matrix')
        aux = [];
        for jj = 1:opt.input.general_constraints.matrix.dim(2)
            aux = [aux; input_matrix(:,jj)];
        end
        Param = [Param; aux];
    end
    
    if isfield(opt.input,'general_constraints') && isfield(opt.input.general_constraints,'vector')
        Param = [Param; input_vector];
    end
    
% if there are any input to stage cost
    if isfield(opt.input,'stage_costs') && isfield(opt.input.stage_costs,'matrix')
        aux = [];
        for jj = 1:opt.input.stage_costs.matrix.dim(2)
            aux = [aux; input_matrix(:,jj)];
        end
        Param = [Param; aux];
    end
    
    if isfield(opt.input,'stage_costs') && isfield(opt.input.stage_costs,'vector')
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
        args.ubx(length(args.ubx)+1:length(args.ubx)+length(parameters_const)) = opt.constraints.parameters.upper;
        args.lbx(length(args.lbx)+1:length(args.lbx)+length(parameters_const)) = opt.constraints.parameters.lower;
    else
        args.ubx(length(args.ubx)+1:length(args.ubx)+length(parameters_const)) = +inf;
        args.lbx(length(args.lbx)+1:length(args.lbx)+length(parameters_const)) = -inf;
    end
end

if (length(OPT_variables)~=length(args.lbx))
    error('MPC error: Number of variables and respective bounds are different')
end


%% generate solver
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
        options.terminationTolerance = 1e-8;
        options.boundTolerance = 1e-8;
        options.printLevel = 'none';
        options.error_on_fail = 0;
        solver = qpsol('solver','qpoases',OPC,options);
end

