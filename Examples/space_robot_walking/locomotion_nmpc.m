import casadi.*
clc
clear
close all
%% Import robot 
%--- URDF filename ---%

filename = 'VISPA_crawling.urdf';
SCname = 'pulsar_PL_SP_CMG.urdf';
%filename= 'SC_3DoF.urdf';
%filenameSC = 'Astrolabe_simple';
%--- Create robot model ---%
[robot,robot_keys] = urdf2robot_flex_visu(filename);
n_q = robot.n_q;

%--- Create spacecraft model ---%
[SC,SC_keys] = urdf2robot_flex_visu(SCname);
n_w = SC.n_q; %number of wheels
%n_q = 12;
% define state variables

q  = SX.sym('q', n_q,1);
qdot = SX.sym('qd',n_q,1);
tau_q = SX.sym('tau_q',n_q,1); 

% torso Euler angles and torso frame angular vel in XYZ convention
Eul_ang = SX.sym('theta', 3,1);
phi = Eul_ang(1);
theta = Eul_ang(2);
psi = Eul_ang(3);

%theta0 = [phi; theta; psi];
Rt = Angles123_DCM(Eul_ang);

rt = SX.sym('rt', 3,1);

qt = [Eul_ang; rt];
rt_dot = SX.sym('rt_dot', 3,1);
omega_tt = SX.sym('omega_tt', 3,1);

qt_dot = [omega_tt ; rt_dot];

gen_vel = [qt_dot; qdot];
%% forward dynamics robot
%--- Kinematics ---%
%Kinematics
[RJ,RL,rJ,rL,e,g]=Kinematics_sym(Rt,rt,q,robot);
%Diferential Kinematics
[Bij,Bi0,Pt,pm]=DiffKinematics_sym(Rt,rt,rL,e,g,robot);
%Velocities
[tt,tq]=Velocities_sym(Bij,Bi0,Pt,pm,qt_dot,qdot,robot);
% %Jacobian of the last link
% [J0n, Jmn]=Jacob(rL(1:3,end),r0,rL,P0,pm,robot.n_links_joints,robot);
N    = NOC_sym(rt,rL,Pt,pm,robot);
Ndot = NOCdot_sym(rt,tt,rL,tq,Pt,pm,robot);
%--- Dynamics ---%
%Inertias in inertial frames
[It,Iq]=I_I_sym(Rt,RL,robot);

%Generalized Inertia matrix
H = GIM_NOC_sym(N,It,Iq,robot); 
%Generalized Convective Inertia matrix
C = CIM_NOC_sym(N,Ndot,tt,tq,It,Iq,robot);
% H = [Ht, Htq; Htq.', H];
% C = [Ct, Ctq; Ctq  , Cm];
Ht = H(1:6,1:6);
Htq = H(1:6, 7:end);
% %Jacobian of the last link of last chain (it is docked
[Jte, Jqe]=Jacob_sym(rL(1:3,end),rt,rL,Pt,pm,robot.n_links_joints,robot);

%% systems dynamics
invH = simplify(inv(H)); %solve(M, SX.eye(M.size1())); % check inv in casadi in python api 
Qv = C*gen_vel;
real_u = [SX.zeros(6,1); tau_q];

epsilon_ddot = invH*(real_u-Qv);
%
epsilon_ddot = simplify(epsilon_ddot);

%% MPC

opt.N = 10;  
opt.dt = 0.1;
opt.n_controls  = n_q;
opt.n_states    = 2*(6+n_q);
opt.model.function = [[qt_dot; qdot]; epsilon_ddot];
opt.model.states   =  [qt; q; qt_dot; qdot];
opt.model.controls = [tau_q];
opt.continuous_model.integration = 'euler';

% Define parameters
opt.parameters.name = {'Ref'};
opt.parameters.dim = [opt.n_states 1];

% control and state constraints
Theta = pi/4*ones(3,1);
R = inf*ones(3,1);
Q = pi*ones(n_q,1);
Thetad = 1*ones(3,1);
V = 5*ones(3,1);
Qd = 1*ones(n_q,1);

state_constraints = [Theta; R; Q; Thetad; V; Qd]; 
opt.constraints.states.upper  =  state_constraints;
opt.constraints.states.lower  = -state_constraints;

control_constraints = ones(opt.n_controls,1);
opt.constraints.control.upper =  control_constraints;
opt.constraints.control.lower = -control_constraints;

% end point tracking constraint but overwritten by saturation inequality
% constraint ! 
opt.constraints.general.function = @(x,varargin) x(:,end)-varargin{:};

opt.constraints.general.parameters  = {'Ref'};
opt.constraints.general.type        = 'equality';
opt.constraints.general.elements 	= 'end';

% x(6+robot.n_q + 1:6+robot.n_q + 6,:) is xt_dot
% x(6+robot.n_q+7:6+robot.n_q+ 6+robot.n_q,:)) is q_dot
h_wheel_max = ones(6,1);
opt.constraints.general.function = @(x, varargin) (Ht*x(6+robot.n_q + 1:6+robot.n_q + 6,:)+ Htq*x(6+robot.n_q+7:6+robot.n_q+ 6+robot.n_q,:)) - h_wheel_max;
opt.constraints.general.parameters  = {'x'};
opt.constraints.general.type        = 'inequality';
opt.constraints.general.elements 	= 'N';

% Define costs
Qc = 100*eye(opt.n_states);
Rc = 0.1*eye(opt.n_controls);
%K_om = 10e6;

opt.costs.stage.parameters = {'Ref'};
opt.costs.stage.function = @(x,u,varargin) (x-varargin{:}(1:opt.n_states))'*Qc*(x-varargin{:}(1:opt.n_states)) + ...
                                           u'*Rc*u;
                                       
% Define inputs to optimization
opt.input.vector = {'Ref'};

% Define the solver and generate it
opt.solver = 'ipopt';
tic
[solver,args] = build_mpc(opt);
toc

%% simulation dummy
tmax = 200;


%--- System_formating---%

[robot_sim]= robot_slx_format_walking(robot);
[SC_sim]= robot_slx_format_walking(SC);

robot_sim.id_EE = robot_sim.n_q;                       % Indice of the end-effector joint/link that is docked

%--- Initial conditions sat---%
qw_init = zeros(SC.n_q, 1);
Rs     = eye(3);                           % Sat's intial rotation from Rs to Rine
rs     = zeros(3,1);                       % Sat's intial position in Rine
xs0     = [0;0;0;zeros(3,1)];             % Sar's intial state=>[eulerXYZ;rt]  
q_w0     = qw_init;                           % wheels's initial value  
qs_dot0     = zeros(6,1);                   % Sats's initial velocity in Rine

qw_dot_init     = q_w0*0;                         % wheels' initial velocity in Rine                 
qw = qw_init; 
qs_dot = qs_dot0;
qw_dot = qw_dot_init;

%--- Initial conditions RObot---%

q_r=-[0 0 pi/2 0 -pi/2+pi/3 0]';
q_l= [0.4 -0.2 -0.98 3.35 0.36 0]';

q_init = [q_r; q_l];
Rt     = eye(3);                           % Base's intial rotation from Rt to Rine
rt     = zeros(3,1);                       % Base's intial position in Rine
xt0     = [0;0;0;zeros(3,1)];             % Base's intial state=>[eulerXYZ;rt]  
%q0     = q_init;                           % Actuators's initial value  
qt_dot0     = zeros(6,1);                   % Base's initial velocity in Rine
q_dot_init     = q_init*0;                         % Actuators' initial velocity in Rine                 
q = q_init; 
qt_dot = qt_dot0;
q_dot = q_dot_init;

xsimu(:,1) = x0;                    % xsimu contains the history of states
u0 = zeros(opt.n_controls,opt.N);                % two control inputs for each robot
X0 = zeros(opt.n_states*opt.N,1);   % initialization of the states decision variables for each shooting point
args.x0 = [X0;reshape(u0',opt.N,1);zeros(opt.n_states,1)]; 

args.x0 = [X0;reshape(u0',2*opt.N,1);zeros(opt.n_states,1);zeros(opt.n_states,1);zeros(opt.n_controls,1)]; 
sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
u(:,t) = full(sol.x(opt.n_states*opt.N+1));

