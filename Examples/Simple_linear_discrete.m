clear all;
close all;
import casadi.*
opt.N = 10;
opt.n_states   	= 2;
opt.n_controls 	= 1;

A = [1 1; 0 1];
B = [0;1];

opt.model.function  	= @(x,u) A*x+B*u;

opt.constraints.states.upper =  [1; 1]; 
opt.constraints.states.lower = -[1; 1]; 
opt.constraints.control.upper =  [0.1];
opt.constraints.control.lower = -[0.1];

Q = eye(opt.n_states);
R = eye(opt.n_controls);
opt.costs.stage.function = @(x,u,varargin) x'*Q*x + u'*R*u;


opt.solver = 'ipopt';
[solver,args_mpc] = build_mpc(opt);

xsimu = [10;0];
x0 = zeros(opt.N*opt.n_states+opt.N*opt.n_controls,1);

for t = 1:100
    p = xsimu(:,t);
    MPC_sol = solver('x0', x0, 'lbx', args_mpc.lbx, 'ubx', args_mpc.ubx,'lbg', args_mpc.lbg, 'ubg', args_mpc.ubg,'p',p);
    aux = full(MPC_sol.x(opt.n_states*(opt.N)+opt.n_states+1:opt.n_states*(opt.N)+opt.N*opt.n_controls))';
    u(:,t) = aux(1);
    xsimu(:,t+1) = A*xsimu(:,t)+ B*u(:,t);
end

figure
hold on
stairs(xsimu(1,:),'b')
stairs(xsimu(2,:),'r')

figure
stairs(u(1,:))