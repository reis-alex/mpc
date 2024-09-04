clear all
close all
clc

%% define state equation 

%masses
m1 = 10; %kg
m2 = 2; %kg

%stiffness coef
k1 = 100; %N/m
k2 = 100; %N/m

% damping coef
a1 = 0.2; %N.s/m
a2 = 0.2; %N.s/m


% State space model

A=[[            0,     1,      0,     0];
    [-(k1 + k2)/m1, -a2/m1,  k2/m1,  a2/m1];
    [            0,     0,      0,     1];
    [        k2/m2,  a2/m2, -k2/m2, -a2/m2]];

B = [ 0 k1/m1 0 0].';
C = [ 0,0,    1,0];
D=0;

[p,~] = size(C);
[n,m] = size(B);

Q = eye(n);
R = eye(m);
[K,P] = dlqr(A,B,Q,R); K=-K;

x = sdpvar(n,1); u = sdpvar(m,1);
xbound = 5; ubound = 0.2;
Xc = Polyhedron('A',vertcat(eye(n),-eye(n)),'b',xbound*ones(2*n,1));
Uc = Polyhedron('A',vertcat(eye(m),-eye(m)),'b',ubound*ones(2*m,1));
Z  = Xc*Uc; Z.minHRep();


ZMatrix = [eye(n) zeros(n) zeros(n,m); K -K eye(m);...
          zeros(n) eye(n) zeros(n,m); zeros(m,n) zeros(m,n) eye(m)];
ZMatrix = Polyhedron('A',blkdiag(Z.A,Z.A)*ZMatrix,'b',vertcat(Z.b,0.99*Z.b));

Aw = [A+B*K -B*K B; zeros(n) eye(n) zeros(n,m); zeros(m,n) zeros(m,n) eye(m)];

Aw =  LTISystem('A',Aw);
Omega = Aw.invariantSet();
Omega = Omega.intersect(ZMatrix);
Omega.minHRep();
T = 100*P;
%% Define MPC
opt.N           = 10;
opt.n_controls  = m;
opt.n_states    = n;
opt.model.type	= 'linear';
opt.model.A     = A;
opt.model.B     = B;

xs = SX.sym('Xs',opt.n_states);
us = SX.sym('Us',opt.n_controls);


% Define costs
opt.costs.stage.function = @(x,u,param) (x-param(1:opt.n_states))'*Q*(x-param(1:opt.n_states)) + ...
                                         (u-param(opt.n_states+1:end))'*R*(u-param(opt.n_states+1:end));
opt.costs.stage.parameters = [xs;us];

ref = SX.sym('Ref',opt.n_states);
opt.costs.terminal.function = @(x,param) (x-param(1:opt.n_states))'*P*(x-param(1:opt.n_states)) + ...
                                           (param(1:opt.n_states)-param(opt.n_states+1:end))'*T*(param(1:opt.n_states)-param(opt.n_states+1:end));
opt.costs.terminal.parameters = [xs;ref];

%% Define constraints
% terminal constraints
opt.constraints.terminal.set = Omega;
opt.constraints.terminal.parameters = [xs;us];

% control and state constraints
opt.constraints.polyhedral = Xc;
opt.constraints.control.upper = [0.2 0.2];
opt.constraints.control.lower = [-0.2 -0.2];

% constraint on parameters
opt.constraints.parameters.variables = [xs;us];

%% Define inputs to optimization
opt.input.vector = ref;

%% Define the solver and generate it
opt.solver = 'ipopt';
[solver,args] = build_mpc(opt);

%% Simulation loop
tmax = 200;

x0 = [-5;0;-5;0];            % initial condition.
xsimu(:,1) = x0;            % xsimu contains the history of states
u0 = zeros(m,opt.N);            % two control inputs for each robot
X0 = zeros(opt.n_states*opt.N,1);      % initialization of the states decision variables

% Intialize MPC
u = [];
args.x0 = [X0;reshape(u0',2*opt.N,1);zeros(opt.n_states,1);zeros(opt.n_states,1);zeros(opt.n_controls,1)]; 

for t = 1:tmax

    if t<=50
        yref = [5;-5];
    end
    if t>50 && t<=100
        yref = [-5;0];
    end
    if t>100 && t<=150
        yref = [2;1];
    end
    if t>150 && t<=200
        yref = [-5;-5];
    end
    refsimu(:,t) = yref;

    xs = pinv([A-eye(n) B; C zeros(p,m)])*[zeros(n,1);yref];

    % set the values of the parameters vector
    args.p = [xsimu(:,t);xs(1:opt.n_states)];                                              
    
    % initial value of the optimization variables

    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);

    % get control sequence from MPC
    for i = 0:opt.N
        usol(:,i+1) = full(sol.x(opt.n_states*opt.N+1+i*opt.n_controls:opt.n_states*opt.N+i*opt.n_controls+2))';
    end
    u(:,t) = usol(:,1);

    % get artificial reference
    ya(:,t) = C*reshape(full(sol.x(opt.N*opt.n_states+opt.N*opt.n_controls+opt.n_states+1:opt.N*opt.n_states+opt.N*opt.n_controls+2*opt.n_states)),opt.n_states,1);

    xsimu(:,t+1) = A*xsimu(:,t) + B*u(:,t);
    y(:,t) = C*xsimu(:,t);
%     X0 = full(sol.x(1:N*n_states));
    args.x0 = full(sol.x); 
end