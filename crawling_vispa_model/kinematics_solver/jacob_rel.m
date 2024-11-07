function [Jp,Rp,rp]=jacob_rel(qm,ind,robot,range)
% Computes the geometric Jacobian of ind-th `rL` in the robot body frame
[~,RL,~,rL,e,g]=Kinematics(eye(3),[0 0 0]',qm,robot);
[~,~,P0,pm]=DiffKinematics(eye(3),[0 0 0]',rL,e,g,robot);
Rp=RL(:,:,ind);
rp=rL(:,ind);
[~, Jm]=Jacob(rp,[0 0 0]',rL,P0,pm,ind,robot);

Jp=[Rp*Jm(1:3,range);Jm(4:6,range)];