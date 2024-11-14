clear all; close all; clc;
addpath(genpath('C:\Users\Alex Reis\Documents\MATLAB\urdf2casadi-matlab-master'))
import casadi.*


robot_path = 'C:\Users\Alex Reis\Desktop\GitHub\mpc\Examples\inverted_pendulum\inverted_pend.urdf';
robot = importrobot(robot_path);
robot.DataFormat = 'row';
robot_acceleration = urdf2casadi.Dynamics.symbolicForwardDynamics(robot_path,0);

% define state variables
q  = SX.sym('q', 1,1);
qd = SX.sym('qd',1,1);
tau = SX.sym('tau',1,1); 

opt.N = 10;  
opt.dt = 0.1;
opt.n_controls  = 1;
opt.n_states    = 2;
opt.model.function = [[qd]; robot_acceleration([q],[qd],[0 0 -10],[tau])]; % a simple double integrator
opt.model.states   =  [q;qd];
opt.model.controls = [tau];
opt.continuous_model.integration = 'euler';

% Define parameters
opt.parameters.name = {'Ref'};
opt.parameters.dim = [opt.n_states 1];


% Define costs
Q = 1000000*eye(opt.n_states);
R = 0.1*eye(opt.n_controls);

opt.costs.stage.parameters = {'Ref'};
opt.costs.stage.function = @(x,u,varargin) (x-varargin{:}(1:2))'*Q*(x-varargin{:}(1:2)) + ...
                                           (u)'*R*(u);
                  
% control and state constraints
xbound = pi;
opt.constraints.states.upper  = xbound*ones(opt.n_states,1);
opt.constraints.states.lower  = -xbound*ones(opt.n_states,1);
opt.constraints.control.upper = 10*ones(opt.n_controls,1);
opt.constraints.control.lower = -10*ones(opt.n_controls,1);
opt.constraints.general.function = @(x,varargin) x(:,end)-varargin{:};
opt.constraints.general.parameters = {'Ref'};
opt.constraints.general.type = 'equality';
opt.constraints.general.elements = 'end';

% Define inputs to optimization
opt.input.vector = {'Ref'};

% Define the solver and generate it
opt.solver = 'ipopt';
tic
[solver,args] = build_mpc(opt);
toc

%% Time simulation
tmax = 200;
x0 = [-30*pi/180;0];                % initial condition.
xsimu(:,1) = x0;                    % xsimu contains the history of states
u0 = zeros(opt.n_controls,opt.N);                % two control inputs for each robot
X0 = zeros(opt.n_states*opt.N,1);   % initialization of the states decision variables
args.x0 = [X0;reshape(u0',opt.N,1);zeros(opt.n_states,1)]; 

for k = 1:tmax

    args.p = [xsimu(:,k);-160*pi/180;0];     
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    u(:,k) = full(sol.x(opt.n_states*(opt.N+1)+1));
    aux2 = robot_acceleration(xsimu(1,k),xsimu(2,k),[0 0 -10],[u(:,k)]);
    xsimu(:,k+1) = xsimu(:,k) + opt.dt*[xsimu(2,k); aux2.full()];
    args.x0 = full(sol.x); 
end

%% Plot

for k = 1:length(xsimu)
    show(robot,xsimu(1,k)');
    axis auto
    view([90 90 90])
    drawnow
end