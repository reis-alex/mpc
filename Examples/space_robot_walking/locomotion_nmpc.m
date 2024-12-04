import casadi.*
clc
clear
close all
%% Import robot
%--- URDF filename ---%

%filename = 'VISPA_crawling.urdf';
SCname = 'pulsar_PL_SP_CMG.urdf';
filename= 'SC_3DoF.urdf';
%filenameSC = 'Astrolabe_simple';
%--- Create robot model ---%
[robot,robot_keys] = urdf2robot_flex_visu(filename);
n_q = robot.n_q;

%--- Create spacecraft model ---%
[SC,SC_keys] = urdf2robot_flex_visu(SCname);
n_w = SC.n_q; %number of wheels
%n_q = 12;
% define state variables

q_sym  = SX.sym('q', n_q,1);
qdot_sym = SX.sym('qd',n_q,1);
tau_q_sym = SX.sym('tau_q',n_q,1);

% torso Euler angles and torso frame angular vel in XYZ convention
Eul_ang_sym = SX.sym('theta', 3,1);
phi = Eul_ang_sym(1);
theta = Eul_ang_sym(2);
psi = Eul_ang_sym(3);

%theta0 = [phi; theta; psi];
Rt_sym = Angles123_DCM(Eul_ang_sym);

rt_sym = SX.sym('rt', 3,1);

qt_sym = [Eul_ang_sym; rt_sym];
rt_dot_sym = SX.sym('rt_dot', 3,1);

omega_tt_sym = SX.sym('omega_tt', 3,1);

qt_dot_omega_sym = [omega_tt_sym ; rt_dot_sym];
T_plate = [1 cos(phi)*tan(theta)  sin(phi)*tan(theta);
           0 cos(phi)            -sin(phi);
           0 sin(phi)/cos(theta)  cos(phi)/cos(theta)];

qt_dot_sym = [ T_plate*omega_tt_sym; rt_dot_sym];
%gen_vel = [qt_dot_omega_sym; qdot_sym];
%% forward dynamics robot
%--- Kinematics ---%
%Kinematics
[RJ,RL,rJ,rL,e,g]=Kinematics_sym(Rt_sym,rt_sym,q_sym,robot);
%Diferential Kinematics
[Bij,Bi0,Pt,pm]=DiffKinematics_sym(Rt_sym,rt_sym,rL,e,g,robot);
%Velocities
[t_t,tq]=Velocities_sym(Bij,Bi0,Pt,pm,qt_dot_omega_sym,qdot_sym,robot);
% %Jacobian of the last link
% [J0n, Jmn]=Jacob(rL(1:3,end),r0,rL,P0,pm,robot.n_links_joints,robot);
N    = NOC_sym(rt_sym,rL,Pt,pm,robot);
Ndot = NOCdot_sym(rt_sym,t_t,rL,tq,Pt,pm,robot);
%--- Dynamics ---%
%Inertias in inertial frames
[It,Iq]=I_I_sym(Rt_sym,RL,robot);

%Generalized Inertia matrix
H = GIM_NOC_sym(N,It,Iq,robot);
%Generalized Convective Inertia matrix
C = CIM_NOC_sym(N,Ndot,t_t,tq,It,Iq,robot);
% H = [Ht, Htq; Htq.', H];
% C = [Ct, Ctq; Ctq  , Cm];
Ht = H(1:6,1:6);
Htq = H(1:6, 7:end);
% %Jacobian of the last link of last chain (it is docked
[Jte, Jqe]=Jacob_sym(rL(1:3,end),rt_sym,rL,Pt,pm,robot.n_links_joints,robot);

% %% systems dynamics
% invH = simplify(inv(H)); %solve(M, SX.eye(M.size1())); % check inv in casadi in python api
% Qv = C*gen_vel;
real_u = [SX.zeros(6,1); tau_q_sym];
%
% epsilon_ddot = invH*(real_u-Qv);
% %
% epsilon_ddot = simplify(epsilon_ddot);


% Constraint build

end_e =4;% effector body index in robot tree

[Jt_ee, Jq_ee]=Jacob_sym(rL(1:3,end_e),rt_sym,rL,Pt,pm,end_e,robot);
J_star = [Jt_ee, Jq_ee];
[Jtdot_ee, Jqdot_ee]   = Jacobdot_sym(rL(1:3,end_e),tq(:,end_e),rt_sym,t_t,rL,tq,Pt,pm,end_e,robot);
J_star_dot = [Jtdot_ee, Jqdot_ee] ;

%% FD with constraint
alpha = 15; %baumgarte velocity stabilizer term1
gamma_holom = -J_star_dot*[qt_dot_omega_sym;qdot_sym] -2*alpha*J_star*[qt_dot_omega_sym;qdot_sym];

nc = 1;
H_aug = simplify([H, -J_star.'; J_star, zeros(6*nc) ]);

full_acc_lambda = simplify(inv(H_aug)*([C*[qt_dot_omega_sym;qdot_sym] + real_u; gamma_holom]));

epsilon_ddot = full_acc_lambda(1:6+robot.n_q);

%% MPC

opt.N = 20;
opt.dt = 0.02;
opt.n_controls  = n_q;
opt.n_states    = 2*(6+n_q);
opt.model.function = [[qt_dot_sym; qdot_sym]; epsilon_ddot];
opt.model.states   = [qt_sym; q_sym; qt_dot_omega_sym; qdot_sym];
opt.model.controls = [tau_q_sym];
opt.continuous_model.integration = 'euler';

% Define parameters
opt.parameters.name = {'Ref'};
opt.parameters.dim = [1, opt.n_states];

% control and state constraints
Theta = pi/6*ones(3,1);
R = inf*ones(3,1);
Q = pi/4*ones(n_q,1);
Thetad = 0.5*ones(3,1);
V = 5*ones(3,1);
Qd = 0.5*ones(n_q,1);

state_constraints = [Theta; R; Q; Thetad; V; Qd];
opt.constraints.states.upper  =  state_constraints;
opt.constraints.states.lower  = -state_constraints;

control_constraints = 50*ones(opt.n_controls,1);
opt.constraints.control.upper =  control_constraints;
opt.constraints.control.lower = -control_constraints;

% end point tracking constraint but overwritten by saturation inequality
% constraint !
opt.constraints.general.parameters  = {'Ref'};
opt.constraints.general.function{1} = @(x,varargin) x(:,end)-varargin{1};
opt.constraints.general.type{1}        = 'equality';
opt.constraints.general.elements{1} 	= 'end';

% x(6+robot.n_q + 1:6+robot.n_q + 6,:) is xt_dot
% x(6+robot.n_q+7:6+robot.n_q+ 6+robot.n_q,:)) is q_dot

h_wheel_max = 100*ones(6,1);

fun = Function('fc',{[qt_sym; q_sym; qt_dot_omega_sym; qdot_sym]},{abs(Ht*opt.model.states(6+robot.n_q + 1:6+robot.n_q + 6,:)+ Htq*opt.model.states(6+robot.n_q+7:6+robot.n_q+ 6+robot.n_q,:)) - h_wheel_max});
opt.constraints.general.function{2} = @(x,varargin) fun(x);
opt.constraints.general.type{2}     = 'inequality';
opt.constraints.general.elements{2}	= 'N';

% Define weights

Weight.xt = 10*[10 10 10 10 10 10 ];
Weight.q= ones(1, robot.n_q);
Weight.xtdot = [10 10 10 0.1  0.1  0.1 ];
Weight.qdot = 0.1*ones(1, robot.n_q);

Qc  = diag([Weight.xt Weight.q Weight.xtdot Weight.qdot]);
%Qc = 100*diag()
Rc = 0.1*eye(opt.n_controls);
%K_om = 10e6;

opt.costs.stage.parameters = {'Ref'};
opt.costs.stage.function = @(x,u,varargin) (x-varargin{:})'*Qc*(x-varargin{:}) + ...
    u'*Rc*u;

% Define inputs to optimization
opt.input.vector = {'Ref'};

% Define the solver and generate it
opt.solver = 'ipopt';
tic
[solver,args] = build_mpc(opt);
toc

%% simulation dummy

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
q_init = pi/4*ones(robot.n_q, 1);
%q_init = [q_r; q_l];
Rt     = eye(3);                            % Base's intial rotation from Rt to Rine
%rt    = zeros(3,1);                        % Base's intial position in Rine
rt    = [50; 50; 0];
xt0     = [0;0;0;zeros(3,1)];               % Base's intial state=>[eulerXYZ;rt]
%q0     = q_init;                           % Actuators's initial value
qt_dot0     = zeros(6,1);                   % Base's initial velocity in Rine
q_dot_init     = q_init*0;                  % Actuators' initial velocity in Rine
q_sym = q_init;
qt_dot_omega_sym = qt_dot0;
q_dot = q_dot_init;

% xsimu contains the history of states
u0 = zeros(opt.n_controls,opt.N);                % two control inputs for each robot
X0 = zeros(opt.n_states*opt.N,1);   % initialization of the states decision variables for each shooting point
%args.x0 = [X0;reshape(u0',opt.N,1);zeros(opt.n_states,1)];

xsimu(:,1) = [xt0; q_init; qt_dot0; q_dot_init];
args.x0 = [X0; reshape(u0',opt.N*opt.n_controls,1);zeros(opt.n_states,1)];

ref_theta = zeros(3,1); ref_rt = zeros(3,1); ref_omega = zeros(3,1); ref_q = 0.5*ones(robot.n_q,1); ref_rtdot = zeros(3,1); ref_qdot = zeros(n_q,1);
yref = [ref_theta; ref_rt; ref_q; ref_omega; ref_rtdot; ref_qdot];

args.p = [X0(:,1); yref];
%% double check sol.x  = x0 for first solution
tmax = 10;
dt = opt.dt;
time_vec = 0:dt:tmax;
disp("simulation starts !")
tic
for i = 1:length(time_vec)-1

    if i < 200
        ref_rt = [ 0; 0; 0.5*sin(dt*i)]; % linear position on x-axis reference
    else
        ref_rt = [ 0; 0; 0.5*sin(dt*200)];
    end
    ref_q = 0.4*pi*ones(robot.n_q,1);

    yref = [ref_theta; ref_rt; ref_q; ref_omega; ref_rtdot; ref_qdot];

    refsimu(:,i) = yref;
    args.p = [xsimu(:,i); refsimu(:,i)];
    [usol,xpred,rest,cost,solx] = solve_mpc(opt,solver,args);

    %% simulation here
    xsimu(:,i+1) = xpred(:,1);

    args.x0 = solx;
end
toc
refsimu(:,i+1) = refsimu(:,i);
[h_robot] = get_constraint_value_momentum(xsimu, robot, time_vec);

plot_solution(xsimu,time_vec,refsimu,robot)

figure 
plot(time_vec, abs(h_robot(1:3, :)));