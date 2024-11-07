function [M0_dot, M0m_dot, Mm_dot, H0dot, H0mdot, Hmdot] = M_dot_Lagrange(Im,Bi0,P0,R0,pm,r0,t0,rL,tL,RL,M0, M0m, robot)
%M_DOT_LAGRANGE Summary of this function goes here
%   Detailed explanation goes here

n = robot.n_links_joints;
n_q = robot.n_q;

[I_0d, I_Idot] = I_I_dot(R0,t0,RL,tL,robot);
Bi0_dot = zeros(6,6,robot.n_links_joints, "like", R0);

for i=1:n
    Bi0_dot(1:6,1:6,i)=[zeros(3,3), zeros(3,3); SkewSym(t0(4:6)-tL(4:6,i)), zeros(3,3)];
end

   
M0_dot = [I_0d          zeros(3);
          zeros(3)      zeros(3)];
%  M0_dot = zeros(6, "like", R0);
for i = 1:n
    mi = robot.links(i).mass;
    M0_dot  = M0_dot + Bi0_dot(:,:,i).'*[Im(:,:,i) zeros(3); ...
                                        zeros(3)  eye(3)*mi]*Bi0(:,:,i)...
                     + Bi0(:,:,i).'*[I_Idot(:,:,i) zeros(3); ...
                                    zeros(3)      zeros(3)]*Bi0(:,:,i)...
                     + Bi0(:,:,i).'*[Im(:,:,i) zeros(3); ...
                                    zeros(3)      eye(3)*mi]*Bi0_dot(:,:,i);

end
    
M0m_dot = zeros(6, n_q);

for i = 1:n
    
    [~, Jmi] = Jacob(rL(1:3,i),r0,rL,P0,pm,i,robot);
    [~, Jmi_dot] = Jacobdot(rL(1:3,i),tL(:,i),r0,t0,rL,tL,P0,pm,i,robot);
    mi = robot.links(i).mass;
    M0m_dot  = M0m_dot + Bi0_dot(:,:,i).'*[Im(:,:,i), zeros(3); ...
                                          zeros(3),  eye(3)*mi]*Jmi...
                               ...
                       + Bi0(:,:,i).'*[I_Idot(:,:,i), zeros(3); ...
                                      zeros(3),      zeros(3)]*Jmi...
                               ...
                       + Bi0(:,:,i).'*[Im(:,:,i), zeros(3); ...
                                      zeros(3) ,    eye(3)*mi]*Jmi_dot;
end                  

Mm_dot = zeros(n_q);

for i = 1:n

    [~, Jmi] = Jacob(rL(1:3,i),r0,rL,P0,pm,i,robot);
    [~, Jmi_dot] = Jacobdot(rL(1:3,i),tL(:,i),r0,t0,rL,tL,P0,pm,i,robot);
    mi = robot.links(i).mass;
    Mm_dot  = Mm_dot + Jmi_dot.'*[Im(:,:,i) zeros(3); ...
                                 zeros(3)  eye(3)*mi]*Jmi...
                     + Jmi.'*[I_Idot(:,:,i) zeros(3); ...
                             zeros(3)      zeros(3)]*Jmi...
                     + Jmi.'*[Im(:,:,i) zeros(3); ...
                             zeros(3)      eye(3)*mi]*Jmi_dot;
end

P0_d = P0_dot(R0,t0,robot);

% Hdot =[P0'*M0_dot*P0 + P0_d'*M0*P0 + P0'*M0*P0_d,  P0'*M0m_dot + P0_d'*M0m;
%        M0m_dot'*P0   + M0m'*P0_d                ,  Mm_dot];

H0dot  = M0_dot ; P0.'*M0_dot*P0 + P0_d.'*M0*P0 + P0.'*M0*P0_d;
H0mdot = M0m_dot ; P0.'*M0m_dot + P0_d.'*M0m;
Hmdot  = Mm_dot;
end

