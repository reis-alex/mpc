function [Rp,rp]=mgd_rel(qm,ind,robot)
% Computes the forward kinematic of ind-th `rL` in the robot body frame

[~,RL,~,rL,~,~]=Kinematics(eye(3),[0 0 0]',qm,robot);

Rp=RL(:,:,ind);
rp=rL(:,ind);
