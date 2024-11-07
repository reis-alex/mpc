function [ dI_q] = dI_dq(r0,R0,RJ,RL,rJ,rL,q, robot)
%DI_DQ Summary of this function goes here
%   Detailed explanation goes here

dI_q = zeros(3,3,robot.n_links_joints,robot.n_q, 'like', R0);
dRi_dql = zeros(3, 'like', R0);
RL_part = zeros(3, 'like', R0);
for l = 1:robot.n_q %l th derivative is a list of Inertia matrixes dI_dql
    
    for i=1:(robot.n_links_joints)
          [dRi_dq,~,~,~,~,~] = Kinematics_diffq(R0,r0,q,l,i,robot);%Part_Diffq_Kine(i,l, R0, r0,RJ,RL,rJ,rL, q, robot);
          dRi_dql = robot.con.branch(i,l)*dRi_dq(:,:,i);
          dI_q(:,:,i,l)  = dRi_dql*robot.links(i).inertia*RL(1:3,1:3,i)'...
                         + RL(1:3,1:3,i)*robot.links(i).inertia*dRi_dql';
    end
end

end


