
%% Build an MPC controller using CasADi

function [solver,args] = build_mpc2(opt)
import casadi.* 
 

%% Gather all parameters (decision variables, inputs - everything that becomes SX.sym)

if isfield(opt,'parameters')
    [rows_p col_p] = size(opt.parameters.dim);
    if rows_p~=length(opt.parameters.name )
        error('MPC error: mismatching number of parameters and corresponding dimension (opt.parameters.dim)')
    end
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
        parameters_stc{i} = parameters{find(strcmp(opt.costs.stage.parameters{i},opt.parameters.name))};
    end
end

% parameters for general cost
if isfield(opt.costs,'general') && isfield(opt.costs.general,'parameters')
    parameters_gnc = [];
    for i = 1:length(opt.costs.general.parameters)
        parameters_gnc{i} = parameters{find(strcmp(opt.costs.general.parameters{i},opt.parameters.name))};
    end
end

% parameters for terminal cost
if isfield(opt.costs,'terminal') && isfield(opt.costs.terminal,'parameters')
    parameters_trc = []; 
    for i = 1:length(opt.costs.terminal.parameters)
        parameters_trc{i} = parameters{find(strcmp(opt.costs.terminal.parameters{i},opt.parameters.name))};
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
if isfield(opt.constraints,'general') && isfield(opt.constraints.general,'parameters')
    parameters_gc = []; 
    for i = 1:length(opt.constraints.general.parameters)
        parameters_gc{i} = parameters{find(strcmp(opt.constraints.general.parameters{i},opt.parameters.name))};
    end
end

% list parameters that are decision variables
if isfield(opt.constraints,'parameters')
    parameters_const = [];
    for i = 1:length(opt.constraints.parameters.name)
        parameters_const = [parameters_const; parameters{find(strcmp(opt.parameters.name,opt.constraints.parameters.name{i}))}];
    end
end

% parameters for inputs
if isfield(opt,'input')
    parameters_input = [];
    for i = 1:length(opt.input.vector)
          parameters_input{i} = parameters{find(strcmp(opt.input.vector{i},opt.parameters.name))};
    end
    if isfield(opt.input,'matrix')
        for j = i+1:i+length(opt.input.matrix)
            parameters_input{j} = parameters{find(strcmp(opt.input.matrix{j-i},opt.parameters.name))};
        end
    end
end

%% Modeling: generate states and control vectors, model function, integration
states = [];
inputs = [];
if ~isfield(opt.model,'states')
    for i = 1:opt.n_states
        states = [states; SX.sym(['x' int2str(i)])];
    end
    
    for i = 1:opt.n_controls
        inputs = [inputs; SX.sym(['u' int2str(i)])];
    end
    model = opt.model.function(states,inputs);
    f = Function('f',{states,inputs},{model}); 
else
    states = opt.model.states;
    inputs = opt.model.controls;
    f = Function('f',{states,inputs},{opt.model.function}); 
end


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


%% Build prediction and cost function, initialize prediction

if isfield(opt.model,'sizeperturbation')
    aux = parameters(strcmp(opt.parameters.name,'d'));
    U           = vertcat(SX.sym('U',opt.n_controls,opt.N),reshape(aux{:},1,opt.N));
else
    U           = SX.sym('U',opt.n_controls,opt.N);         % Decision variables (controls)
    nb_moves = opt.N;
    ubase = U;
end

if isfield(opt,'input_parametrization')
    ubase           = SX.sym('U',opt.n_controls,opt.input_parametrization.nb_moves); 
    U = opt.input_parametrization.function(ubase);
    nb_moves = opt.input_parametrization.nb_moves;
end


Init_states = SX.sym('Init',opt.n_states); 
X           = SX.sym('X',opt.n_states,(opt.N));       % A vector that represents the states over the optimization problem.
X = horzcat(Init_states,X);
obj         = 0;
g           = [];               

% Determine the cost function to be used: simple quadratic or a custom
if isfield(opt.costs.stage,'function')
    stagecost_fun    = opt.costs.stage.function;
else
    stagecost_fun    = @(x,u,extra) x'*opt.costs.stage.Q*x + u'*opt.costs.stage.R*u + 0*extra;
end

% Determine if cost function parameters are static or iterated upon
if isfield(opt.costs.stage,'sort_parameter')
    stc_fixed = opt.costs.stage.sort_parameter.fixed;
    stc_var = opt.costs.stage.sort_parameter.var;
    sort_param_stc = @(k) vertcat(parameters_stc{stc_fixed}, parameters_stc{stc_var(k)});
else
    if isfield(opt.costs.stage,'parameters')
        sort_param_stc = @(k) vertcat(parameters_stc{:});
    else
        sort_param_stc = @(k) 0;
    end
end

%% Prediction loop, integration schemes
% prediction loop

if isfield(opt,'continuous_model')
    switch opt.continuous_model.integration
        case 'RK4'
            for k = 1:opt.N
                obj             = obj + stagecost_fun(X(:,k),U(1:opt.n_controls,k),sort_param_stc(k));
                k1 = f(X(:,k),U(:,k));
                k2 = f(X(:,k)+opt.dt/2*k1,U(:,k));
                k3 = f(X(:,k)+opt.dt/2*k2,U(:,k));
                k4 = f(X(:,k)+opt.dt*k3,U(:,k));
                next_RK4 = X(:,k)+(opt.dt/6)*(k1+2*k2+2*k3+k4);
                g             = [g; X(:,k+1)-next_RK4]; 
            end
        case 'euler'
                for k = 1:opt.N
                obj             = obj + stagecost_fun(X(:,k),U(:,k),sort_param_stc(k));
                g               = [g; X(:,k+1)-(X(:,k)+opt.dt*f(X(:,k),U(:,k)))]; 
                end
    end
else % if the model is already discrete
                for k = 1:opt.N
                    obj             = obj + stagecost_fun(X(:,k),U(:,k),sort_param_stc(k));
                    f_value         = f(X(:,k),U(:,k));
                    st_next         = integ(X(:,k),opt.dt,f_value);
                    g               = [g; X(:,k+1)-st_next];
                end
end

% if terminal cost
if isfield(opt.costs,'terminal')
    if isfield(opt.costs,'terminal') && isfield(opt.costs.terminal,'parameters')
        obj = obj + opt.costs.terminal.function(X(:,end),parameters_trc{:});
    else
        obj = obj + opt.costs.terminal.function(X(:,end));
    end
end

% if extra cost
if isfield(opt.costs,'general')
    obj = obj + opt.costs.general.function(X,parameters_gnc{:});
end

if sum(size(obj))~=2
    error('Cost function is not a scalar.')
end
%% Define vector of decision variables and its constraints
               
% add terminal constraint variables to the list
% if isfield(opt,'constraints') && isfield(opt.constraints,'terminal')
%     OPT_variables = [reshape(X(:,2:end),opt.n_states*(opt.N),1);
%                     reshape(U,opt.n_controls*opt.N,1)];
% else
    OPT_variables   = [reshape(X(:,2:end),opt.n_states*(opt.N),1);
                       reshape(ubase,opt.n_controls*nb_moves,1);];
% end

% add extra (constrained) variables to the list
if isfield(opt,'constraints') && isfield(opt.constraints,'parameters')
    OPT_variables   = [OPT_variables; 
                       parameters_const];
end

% inequalities for the decision variables
args = struct;
args.lbx = [];
args.ubx = [];

% bounds for variables (states)
if isfield(opt.constraints,'polyhedral') || ~isfield(opt,'constraints')
    args.lbx(length(args.lbx)+1:opt.n_states*(opt.N),1) = -inf;
    args.ubx(length(args.ubx)+1:opt.n_states*(opt.N),1) = +inf;
else
    args.lbx(length(args.lbx)+1:length(args.lbx)+opt.n_states*(opt.N),1) = repmat(opt.constraints.states.lower,opt.N,1);
    args.ubx(length(args.ubx)+1:length(args.ubx)+opt.n_states*(opt.N),1) = repmat(opt.constraints.states.upper,opt.N,1);
end

% bounds for variables (controls)
if isfield(opt.constraints,'control')
        args.lbx(length(args.lbx)+1:length(args.lbx)+opt.n_controls*nb_moves,1) = repmat(opt.constraints.control.lower,nb_moves,1);
        args.ubx(length(args.ubx)+1:length(args.ubx)+opt.n_controls*nb_moves,1) = repmat(opt.constraints.control.upper,nb_moves,1);
else
        args.lbx(length(args.lbx)+1:length(args.lbx)+opt.n_controls*nb_moves,1) = -inf;
        args.ubx(length(args.ubx)+1:length(args.ubx)+opt.n_controls*nb_moves,1) = +inf;
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


%% Define functional constraints (lbg<=g(x)<=ubg)

% if constraints in states are polyhedral
if isfield(opt,'constraints') && isfield(opt.constraints,'polyhedral')
        for i = 1:opt.N+1
            g = [g; opt.constraints.polyhedral.A*X(:,i)-opt.constraints.polyhedral.b];
        end
end

% if terminal constraints
if isfield(opt,'constraints') && isfield(opt.constraints,'terminal') && isfield(opt.constraints.terminal,'parameters')
    g   = [g; opt.constraints.terminal.set.A*vertcat(X(:,end),parameters_tc)-opt.constraints.terminal.set.b];
else
    if isfield(opt,'constraints') && isfield(opt.constraints,'terminal')
        g   = [g; opt.constraints.terminal.set.A*X(:,end)-opt.constraints.terminal.set.b];
    end
end

% if there are general constraints
if isfield(opt,'constraints') && isfield(opt.constraints,'general')
    
    for jj = 1:length(opt.constraints.general.function)
        size_gc{jj} = 0;
%         switch opt.constraints.general.elements{jj}
%             case 'N'
%                 for i = 1:opt.N
%                     g = [g; opt.constraints.general.function{jj}(X(:,i),parameters_gc{:})];
%                     size_gc{jj} = size_gc{jj} + length(opt.constraints.general.function{jj}(X(:,i),parameters_gc{:}));
%                 end
%             case 'end'
%                 g = [g; opt.constraints.general.function{jj}(X(:,end),vertcat(parameters_gc{:}))]; 
%                 size_gc{jj} = size_gc{jj} +  length(opt.constraints.general.function{jj}(X(:,end),vertcat(parameters_gc{:})));
%             case 'all'
                g = [g; opt.constraints.general.function{jj}(X,U,vertcat(parameters_gc{:}))]; 
                size_gc{jj} = size_gc{jj} +  length(opt.constraints.general.function{jj}(X,U,vertcat(parameters_gc{:})));
%         end
    end
end

%% define external parameters and problem structure

Param = [Init_states];
if isfield(opt,'input')
    if isfield(opt.input,'vector')
        Param = [Param; vertcat(parameters_input{:})];
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

% equality for x(k+1)-x(k)
args.lbg(1:opt.n_states*(opt.N)) = 0.0;
args.ubg(1:opt.n_states*(opt.N)) = 0.0;

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
    for j = 1:length(opt.constraints.general.function)
        switch opt.constraints.general.type{j}
            case 'inequality'
                bound_gc = -inf;
            case 'equality'
                bound_gc = 0;
        end
        args.lbg(length(args.ubg)+1:length(args.ubg)+size_gc{j}) = bound_gc;
        args.ubg(length(args.ubg)+1:length(args.ubg)+size_gc{j}) = 0;
    end
end

if (length(g)~=length(args.lbg))
    error('MPC error: Number of constraints(g) and respective bounds (lbg or ubg) are different')
end

% get position of the first control in the whole decision vector
names = {};
aux = [];
for i = 1:length(OPT_variables)
    aux = OPT_variables(i);
    names{i} = aux.name;
end
args.vars{1} = OPT_variables;
args.vars{2} = Param;
args.vars{3} = find(strcmp(names,'U_0')):find(strcmp(names,'U_0'))+opt.n_controls-1;
args.debug = Function('g_fun', {Param, OPT_variables}, {g});
%% generate solver
opts                        = struct;
switch opt.solver
    case 'ipopt'
        opts.ipopt.max_iter         = 400;
        opts.ipopt.print_level      = 0;
        opts.print_time             = 0;
        opts.ipopt.tol = 1e-5;
        opts.ipopt.acceptable_tol = 1e-5;
        opts.ipopt.acceptable_obj_change_tol = 1e-5;
        opts.ipopt.acceptable_constr_viol_tol = 1e-5;
        opts.ipopt.honor_original_bounds = 'yes';
        opts.ipopt.bound_relax_factor = 0;
%         opts.ipopt.diverging_iterates_tol = 1e3;
%         opts.ipopt.expect_infeasible_problem = 'yes';
        opts.ipopt.nlp_scaling_method = 'gradient-based';
%         opts.ipopt.constr_viol_tol = 1e-10;  % Reduce constraint violation tolerance
%         opts.ipopt.tol = 1e-8;  % Reduce overall solution tolerance
%         opts.ipopt.dual_inf_tol = 1e-6;  % Reduce dual infeasibility tolerance
%         opts.ipopt.compl_inf_tol = 1e-8;  % Improve complementarity conditions
%         opts.ipopt.bound_push = 1e-10;
%         opts.ipopt.bound_frac = 1e-10;
        solver = nlpsol('solver','ipopt',OPC,opts);
    case 'qpoases'
        options.terminationTolerance = 1e-8;
        options.boundTolerance = 1e-8;
        options.printLevel = 'none';
        options.error_on_fail = 0;
        solver = qpsol('solver','qpoases',OPC,options);
end

%% Plot feasible region (initial conditions)
if isfield(opt,'feasibility') && opt.feasibility == 1

    d = SX.sym('d',3,1);
    obj2   = -d'*Param(1:3);
    OPC2   = struct('f', obj2, ...
        'x', vertcat(Param(1:3),OPT_variables), ...
        'g', g, ...
        'p', vertcat(Param(4:end),d));
    solver2 = nlpsol('solver','ipopt',OPC2,opts);

lbx2 = vertcat(opt.constraints.states.lower(1:3), args.lbx);
ubx2 = vertcat(opt.constraints.states.upper(1:3), args.ubx);
x0 = zeros(length(vertcat(Param(1:3),OPT_variables)),1);
for i = 1:50
    rd = abs(rand(3,1));
    dn = -1+2*rd;
    sol = solver2('x0', x0, ...
                  'lbx', lbx2,...
                  'ubx', ubx2,...
                  'lbg', args.lbg,...
                  'ubg', args.ubg,...
                  'p',vertcat(0.1,0,0,1,0,0,0.5,0,0,dn));
    x0 = full(sol.x);
    vertices(:,i) = full(sol.x(1:3));
end
feas = Polyhedron('V',vertices([1 2],:)');
end