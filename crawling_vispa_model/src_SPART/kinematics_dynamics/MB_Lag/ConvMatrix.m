function [c0, c0q, cq0, cq] = ConvMatrix(R0, r0, q0_dot, rL, t0, tL, RL, rJ, RJ, k, g, dTJ, dTL, dTJ0,dTL0, I0,Im, Bij,Bi0, P0, pm, q, q_dot, q0, robot)
%CONVMATRIX Summary of this function goes here
%   Computes Coriolis matrix from Lagrange analytical expression - Valid

nq = robot.n_q;

c0  = zeros(6,"like", R0);
c0q = zeros(6, robot.n_q,"like", R0);
cq0 = zeros(robot.n_q, 6,"like", R0);
cq  = zeros(robot.n_q,"like", R0);

grad_H0_q0 = zeros(6,6,6,"like", R0);
grad_H0q_q0 = zeros(6,nq,6,"like", R0);
grad_Hq_q0 = zeros(nq,nq,6,"like", R0);
grad_H0_q = zeros(6,6,nq,"like", R0);
grad_H0q_q= zeros(6,nq,nq,"like", R0);
grad_Hq_q = zeros(nq,nq,nq,"like", R0);

[~, ~, ~, M0, M0q, ~] = M_lagrange(I0, Im, rL, Bi0, P0, r0, pm, robot);

%[H0dot, H0mdot, Hdot] = H_dot_Lagrange(Im,Bi0,P0,R0,pm,r0,q0_dot,rL,tL,RL,robot);
[~,~,~,H0dot, H0qdot, Hdot] = M_dot_Lagrange(Im,Bi0,P0,R0,pm,r0,t0,rL,tL,RL,M0, M0q, robot);

[grad_H0_q0, grad_H0q_q0, grad_Hq_q0, grad_H0_q, grad_H0q_q, grad_Hq_q] ...
    = dH_dq_gradients(I0,Im,R0,r0,q0_dot,RJ,RL,rJ,rL,k,g,dTJ,dTL,dTJ0,dTL0,Bij,Bi0,P0,pm,q,q0,robot);



for i =1:6
    c0(i,:)  = -0.5*(q0_dot'*grad_H0_q0(:,:,i) + q_dot'*grad_H0q_q0(:,:,i)');
    c0q(i,:) = -0.5*(q0_dot'*grad_H0q_q0(:,:,i)+ q_dot'*grad_Hq_q0(:,:,i));
end

for i = 1:nq
    cq0(i,:) = -0.5*(q0_dot'*grad_H0_q(:,:,i)  + q_dot'*grad_H0q_q(:,:,i)');
    cq(i,:)  = -0.5*(q0_dot'*grad_H0q_q(:,:,i) + q_dot'*grad_Hq_q(:,:,i));
end

c0  = c0  + H0dot;
c0q = c0q + H0qdot;
cq0 = cq0 + H0qdot';
cq  = cq  + Hdot;

end
