# MPC and robotics: using an URDF to generate a prediction model

## Introduction

This examples illustrates how to generate an MPC solver by parsing a prediction model obtained from an URDF file. This code is based on [URDF2Casadi](https://github.com/robotology/urdf2casadi-matlab) from the library [Robotology](https://github.com/robotology).

The example is a simple inverted pendulum.

```matlab
import casadi.*

% define state variables
q  = SX.sym('q');
qd = SX.sym('qd');
tau = SX.sym('tau'); 

robot_path = '';
robot = importrobot(robot_path);
robot_acceleration = urdf2casadi.Dynamics.symbolicForwardDynamics(robot_path,0);

opt. N = 12;  
opt.dt = 0.1;
opt.n_controls  = 1;
opt.n_states    = 2;
opt.model.function = [[qd]; robot_acceleration([q1],[qd],[0 0 -10],[tau])]; % a simple double integrator
opt.model.states   =  [q;qd];
opt.model.controls = [tau];
opt.continuous_model.integration = 'euler';
```