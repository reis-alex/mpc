import casadi.*

robot_path = 'inverted_pend.urdf';
robot = importrobot(robot_path);
robot_acceleration = urdf2casadi.Dynamics.symbolicForwardDynamics(robot_path,0);

% define state variables
q  = SX.sym('q', 1,1);
qd = SX.sym('qd',1,1);
tau = SX.sym('tau',1,1); 

opt.N = 12;  
opt.dt = 0.1;
opt.n_controls  = 1;
opt.n_states    = 2;
opt.model.function = [[qd]; robot_acceleration([q],[qd],[0 0 -10],[tau])]; % a simple double integrator
opt.model.states   =  [q;qd];
opt.model.controls = [tau];
opt.continuous_model.integration = 'euler';