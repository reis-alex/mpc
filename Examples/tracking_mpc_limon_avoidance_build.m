clear all
close all
clc
import casadi.*

addpath(genpath('D:\adosreis\Documents\Optimization'))
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
opt.model.function	= @(x,u) A*x+B*u;
radius = 16;

% Define parameters
opt.parameters.name = {'Xs','Us','Ref','Matrix','Vector'};
opt.parameters.dim = [opt.n_states 1; opt.n_controls 1; opt.n_states 1; 10 2; 10 1];

% Define costs

opt.costs.stage.parameters = {'Xs','Us'};
opt.costs.stage.function = @(x,u,varargin) (x-varargin{1})'*Q*(x-varargin{1}) + ...
                                           (u-varargin{2})'*R*(u-varargin{2}) + ...
                                            1000*max([0, -((C*x)'*eye(2)*(C*x)-radius)])^2;

opt.costs.terminal.parameters = {'Xs','Us','Ref'};
opt.costs.terminal.function = @(x,varargin) (x-varargin{1})'*P*(x-varargin{1}) + ...
                                            (varargin{1}-varargin{3})'*T*(varargin{1}-varargin{3}) + ...
                                            + 1000*max([0, -((C*varargin{1})'*eye(2)*(C*varargin{1})-radius)])^2;


%% Define constraints
% terminal constraints
opt.constraints.terminal.set = Omega;
opt.constraints.terminal.parameters = {'Xs','Us'};

% control and state constraints
% opt.constraints.polyhedral = Xc;
opt.constraints.states.upper = xbound*ones(opt.n_states,1);
opt.constraints.states.lower = -xbound*ones(opt.n_states,1);
opt.constraints.control.upper = [0.2 0.2]';
opt.constraints.control.lower = [-0.2 -0.2]';

% constraint on parameters
opt.constraints.parameters = {'Xs','Us'};

%% Define inputs to optimization
opt.input.vector = {'Ref'};

%% Define the solver and generate it
opt.solver = 'ipopt';
tic
[solver,args] = build_mpc(opt);
toc
%% Simulation loop
clear i j k
t0 = 0;
tmax = 200;
t(1) = t0;

x0 = [-5;0;-5;0];            % initial condition.
xsimu(:,1) = x0;            % xsimu contains the history of states
u0 = zeros(m,opt.N);            % two control inputs for each robot
X0 = zeros(opt.n_states*opt.N,1);      % initialization of the states decision variables

% Start MPC
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
    args.p = [xsimu(:,t);xs(1:opt.n_states);];                                              
    
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

%% Plot

figure(1)
stairs(1:tmax,u(1,:))
hold on
stairs(1:tmax,u(2,:))

figure(2)
hold on
plot(Xc.projection(1:2),'color','w')
xplot = sdpvar(2,1);
plot(xplot'*eye(2)*xplot<=radius)
plot(xsimu(1,:),xsimu(3,:),'-ok')

figure(3)
t = 0:tmax-1;
hold on
stairs(t,refsimu(1,:),'k')
stairs(t,refsimu(2,:),'--k')
stairs(t,ya(1,:),'r')
stairs(t,ya(2,:),'--r')
stairs(t,y(1,:),'b')
stairs(t,y(2,:),'--b')

figure(4)

hold on
stairs(0:tmax-1,xsimu(1,1:end-1),'b')