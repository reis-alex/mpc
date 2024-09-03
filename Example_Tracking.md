# Linear Tracking MPC [(Limon et al., 2008)](https://www.sciencedirect.com/science/article/abs/pii/S0005109808001106)

This example implements a tracking MPC scheme for linear systems based on artificial references. 

Consider a simple discrete-time double-integrator:

$$
x_{k+1}= \begin{bmatrix} 1 & 1 & 0 & 0\\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 1 \end{bmatrix}x_k + \begin{bmatrix}0 & 0 \\ 1 & 0 \\ 0 & 0 \\ 0 & 1\end{bmatrix}u_k
$$


```matlab
%% Define system, constraint and invariant sets

A = [1 1; 0 1];
A = blkdiag(A,A);
B = [0 0; 1 0; 0 0; 0 1];
C = [1 0 0 0; 0 0 1 0];

[p,~] = size(C);
[n,m] = size(B);

Q = eye(n);
R = eye(m);
[K,P] = dlqr(A,B,Q,R); K=-K;

x = sdpvar(n,1); u = sdpvar(m,1);
xbound = 5; ubound = 0.1;
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
%%
opt.N           = 10;
opt.n_controls  = m;
opt.n_states    = n;
opt.model.type	= 'linear';
opt.model.A     = A;
opt.model.B     = B;

xs = SX.sym('Xs',opt.n_states);
us = SX.sym('Us',opt.n_controls);
radius = 16;

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
tic
[solver,args] = build_mpc(opt);
toc
```