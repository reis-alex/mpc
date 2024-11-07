function [H0, H0q, Hq, M0, M0q, Mq] = M_lagrange(I0, Im, rL, Bi0, P0, r0, pm, robot)
%M_LAGRANGE Summary of this function goes here
%   Detailed explanation goes here
n_q=robot.n_q;
n=robot.n_links_joints;

%--- H matrix ---%

%Base-link inertia matrix
M0_tilde=[I0,zeros(3,3);zeros(3,3),robot.base_link.mass*eye(3)];
M0 = M0_tilde;
M0q=zeros(6,n_q,'like',M0_tilde);
Mq=zeros(n_q,n_q,'like',M0_tilde);

for i = 1:n

    mi = robot.links(i).mass;
    [~, Jmi] = Jacob(rL(1:3,i),r0,rL,P0,pm,i,robot);

    M0  = M0 + Bi0(:,:,i).'*[Im(:,:,i) zeros(3); ...
        zeros(3)  eye(3)*mi]*Bi0(:,:,i);

    M0q = M0q + Bi0(:,:,i).'*[Im(:,:,i), zeros(3); ...
        zeros(3),  eye(3)*mi]*Jmi;

    Mq  = Mq + Jmi.'*[Im(:,:,i) zeros(3); ...
        zeros(3)  eye(3)*mi]*Jmi;
end
% %Pre-allocate Hq
% Mq=zeros(n_q,n_q,'like',M0_tilde);
% for i = 1:n
%
%     [~, Jmi] = Jacob(rL(1:3,i),r0,rL,P0,pm,i,robot);
%     mi = robot.links(i).mass;
%     Mq  = Mq + Jmi'*[Im(:,:,i) zeros(3); ...
%                      zeros(3)  eye(3)*mi]*Jmi;
% end
% %Pre-allocate H0q
% M0q=zeros(6,n_q,'like',M0_tilde);
% for i = 1:n
%
%     [~, Jmi] = Jacob(rL(1:3,i),r0,rL,P0,pm,i,robot);
%
%     mi = robot.links(i).mass;
%     M0q = M0q + Bi0(:,:,i)'*[Im(:,:,i), zeros(3); ...
%                              zeros(3),  eye(3)*mi]*Jmi;
% end
H0 = P0.'*M0*P0;
H0q = P0.'*M0q ;
Hq = Mq;
end

