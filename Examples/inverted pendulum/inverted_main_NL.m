clear 
close all
clc
import casadi.*

%% DEFINE MODEL Physics PARAMETERS
m = 1;
b = 1;
g = 9.81;
l = 1;
I = m*l^2;
%% Compute invariants


% State space model

A=[[         0,     1];
    [(m*g*l/I),    -b]];

B = [ 0 1 ].';
C = [ 1,0 ];

D=0;
Ts = 1;

sys = c2d(ss(A,B,C,D),Ts);
[A,B,C,D] = ssdata(sys);

[p,~] = size(C);
[n,m] = size(B);

Q = 100*eye(n);
R = eye(m);
[K,P] = dlqr(A,B,Q,R); K=-K;

xbound = 10; ubound = 4;
Xc = Polyhedron('A',vertcat(eye(n),-eye(n)),'b',xbound*ones(2*n,1));
Uc = Polyhedron('A',vertcat(eye(m),-eye(m)),'b',ubound*ones(2*m,1));
Z  = Xc*Uc; Z.minHRep();

ZMatrix = [eye(n) zeros(n) zeros(n,m); K -K eye(m);...
          zeros(n) eye(n) zeros(n,m); zeros(m,n) zeros(m,n) eye(m)];
ZMatrix = Polyhedron('A',blkdiag(Z.A,Z.A)*ZMatrix,'b',vertcat(Z.b,0.99*Z.b));

Aw = [A+B*K -B*K B; zeros(n) eye(n) zeros(n,m); zeros(m,n) zeros(m,n) eye(m)];

Aw =  LTISystem('A',Aw,'Ts',Ts);
Omega = Aw.invariantSet();
Omega = Omega.intersect(ZMatrix);
Omega.minHRep();
T = 100*P;

%% define state equation 

theta_ddot = @(x,u)((u - b*x(2) - (m*g*l/I)*sin(x(1))));
                                       
dx_dt = @(x,u)([x(2); u - b*x(2) - (m*g*l/I)*sin(x(1))]);

%% controler set up 

opt.N           = 30;
opt.n_controls  = m;
opt.n_states    = n;
opt.model.type	= 'nonlinear';


% Define parameters
opt.parameters.name = {'Xs','Us','Ref'};
opt.parameters.dim = [opt.n_states 1; opt.n_controls 1; opt.n_states 1];

%% define SS and Cost function
opt.model.function     =  @(x,u)([x(2); u - b*x(2) - (m*g*l/I)*sin(x(1))]);


opt.costs.stage.parameters = {'Xs','Us'};
opt.costs.stage.function = @(x,u,varargin) (x-varargin{:}(1:2))'*Q*(x-varargin{:}(1:2)) + ...
                                           (u-varargin{:}(3))'*R*(u-varargin{:}(3)) + ...
                                           +0;%+ 1000000*max(norm(m*x(2)^2)-10,0)^2;                                      
                                         
opt.costs.terminal.parameters = {'Xs','Us','Ref'};
opt.costs.terminal.function = @(x,varargin) (x-varargin{:}(1:2))'*P*(x-varargin{:}(1:2)) + ...
                                            (varargin{:}(1:2)-varargin{:}(4:5))'*T*(varargin{:}(1:2)-varargin{:}(4:5)) + ...
                                           +0;%+ 1000000*max(norm(m*x(2)^2)-10,0)^2;

%% Define constraints
% terminal constraints
opt.constraints.terminal.set = Omega;
opt.constraints.terminal.parameters = {'Xs','Us'};


% control and state constraints
opt.constraints.polyhedral = Xc;


opt.constraints.control.upper = [ubound ubound];
opt.constraints.control.lower = -[ubound ubound];

% constraint on parameters
opt.constraints.parameters = {'Xs','Us'};

%% Define inputs to optimization
opt.input.vector = {'Ref'};
opt.continuous_model.integration = 'RK4';
opt.dt = Ts;
%% Define the solver and generate it
opt.solver = 'ipopt';
[solver,args] = build_mpc(opt);

%% Simulation loop
tmax = 200;

x0 = [0;0];                     % initial condition.
xsimu(:,1) = x0;                    % xsimu contains the history of states
u0 = zeros(m,opt.N);                % two control inputs for each robot
X0 = zeros(opt.n_states*opt.N,1);   % initialization of the states decision variables

% Intialize MPC
u = [];
args.x0 = [X0;reshape(u0',opt.N,1);zeros(opt.n_states,1);zeros(opt.n_states,1);zeros(opt.n_controls,1)]; 


for t = 1:tmax

    %yref = 0.1*sin(1*t);

    yref = [pi/6];
    refsimu(:,t) = yref;
    
    A=[[        0,     1];
       [(m*g*l/I),    -b]];

    B = [ 0 1].';
    
    xs = pinv([A-eye(n) B; [1 0], zeros(1,m)])*[zeros(n,1);yref];
    
    % set the values of the parameters vector
    args.p = [xsimu(:,t);xs(1:opt.n_states)];                                              
    
    % initial value of the optimization variables

    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    cost_out(:,t) = full(sol.f);

    % get control sequence from MPC
    u(:,t) = full(sol.x(opt.n_states*opt.N+1));
    
    % get artificial reference
    ya(:,t) = C*reshape(full(sol.x(opt.N*opt.n_states+opt.N*opt.n_controls+opt.n_states+1:opt.N*opt.n_states+opt.N*opt.n_controls+2*opt.n_states)),opt.n_states,1);
    
    xsimu(:,t+1) = [xsimu(2,t); (u(:,t) - b*xsimu(2,t) - (m*g*l/I)*sin(xsimu(1,t)))];
    y(:,t) = C*xsimu(:,t);

    constraint(:,t) = m*xsimu(1,t)^2;

    args.x0 = full(sol.x); 
end

figure
plot(u, 'g')
figure
plot(xsimu(1,:), 'g')
hold on
stairs(constraint,'--b')

figure
stairs(1:tmax+1,xsimu(1,:),'r');
hold on
stairs(1:tmax,refsimu,'--b');
