function [J0,Jq]=jacob_torso(qm,ind_a,ind_b,robot)
% Computes the geometric Jacobian of the torso in the robot body frame
[~,~,~,rL,e,g]=Kinematics(eye(3),[0 0 0]',qm,robot);
[~,~,P0,pm]=DiffKinematics(eye(3),[0 0 0]',rL,e,g,robot);
%Rp=RL(:,:,ind);
ra=rL(:,ind_a);
rb=rL(:,ind_b);
[J0a, Jma]=Jacob(ra,[0 0 0]',rL,P0,pm,ind_a,robot);
[J0b, Jmb]=Jacob(rb,[0 0 0]',rL,P0,pm,ind_b,robot);
J0=[J0a;J0b];
Jq=[Jma;Jmb];