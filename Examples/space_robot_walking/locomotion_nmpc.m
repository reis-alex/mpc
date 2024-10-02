import casadi.*
clc
%% Import robot 
%--- URDF filename ---%

filename= 'SC_3DoF.urdf';
%filenameSC = 'Astrolabe_simple';
%--- Create robot model ---%
[robot,robot_keys] = urdf2robot(filename);
n_q = robot.n_q;


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

%--- Kinematics ---%
%Kinematics
[RJ,RL,rJ,rL,e,g]=Kinematics(Rt,rt,q,robot);
%Diferential Kinematics
[Bij,Bi0,Pt,pm]=DiffKinematics(Rt,rt,rL,e,g,robot);
%Velocities
[tt,tq]=Velocities(Bij,Bi0,Pt,pm,qt_dot,qdot,robot);
% %Jacobian of the last link
% [J0n, Jmn]=Jacob(rL(1:3,end),r0,rL,P0,pm,robot.n_links_joints,robot);
N    = NOC(rt,rL,Pt,pm,robot);
Ndot = NOCdot(rt,tt,rL,tq,Pt,pm,robot);
%--- Dynamics ---%
%Inertias in inertial frames
[It,Iq]=I_I(Rt,RL,robot);

%Generalized Inertia matrix
H = GIM_NOC(N,It,Iq,robot); 
%Generalized Convective Inertia matrix
C = CIM_NOC(N,Ndot,tt,tq,It,Iq,robot);
% H = [Ht, Htq; Htq.', H];
% C = [Ct, Ctq; Ctq  , Cm];

% %Jacobian of the last link of last chain (it is docked
 [Jte, Jqe]=Jacob(rL(1:3,end),rt,rL,Pt,pm,robot.n_links_joints,robot);

%% systems dynamics


invH = simplify(inv(H)); %solve(M, SX.eye(M.size1())); % check inv in casadi in python api 
Qv = C*gen_vel;
real_u = [SX.zeros(6,1); tau_q];

epsilon_ddot = invH*(real_u-Qv);
%
epsilon_ddot = simplify(epsilon_ddot);



opt.N = 10;  
opt.dt = 0.1;
opt.n_controls  = 1;
opt.n_states    = 2;
opt.model.function = [[qt_dot; qdot]; epsilon_ddot]; % a simple double integrator
opt.model.states   =  [qt; q; qt_dot; qdot];
opt.model.controls = [tau_q];
opt.continuous_model.integration = 'euler';




