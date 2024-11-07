function [I0_dot,I_I_dot] = I_I_dot(R0, t0, RL,tL,robot)
%I_I_DOT Summary of this function goes here
%   Detailed explanation goes here
I0_dot = SkewSym(t0(1:3))*R0*robot.base_link.inertia*R0.'+R0*robot.base_link.inertia*R0.'*SkewSym(t0(1:3)).';

omega = tL(1:3,:);
I_I_dot= zeros(3,3,robot.n_links_joints, 'like', R0);

for i=1:(robot.n_links_joints)
    I_I_dot(:,:,i) = SkewSym(omega(:,i))*RL(:,:,i)*robot.links(i).inertia*RL(:,:,i).'...
                   + RL(:,:,i)*robot.links(i).inertia*RL(:,:,i).'*SkewSym(omega(:,i)).';
end

