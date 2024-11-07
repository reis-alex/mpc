function [J0i_dq, Jmj_dq]= Jacob_dq(r0,R0,P0,RJ,RL,rJ,rL, k, g, dTJ, dTL,  pm, ind_i, ind_l, dBi0_q, robot)
%JACOB_DQ Summary of this function goes here
%   Detailed explanation goes here
%Pre-allocate
J0i_dq = zeros(6,6, 'like', R0);
Jmj_dq = zeros(6,robot.n_q, 'like', R0);
% dBij_dql  = zeros(6, 'like', R0);
 dpm_dl = zeros(6,1, 'like', R0);
% dg_dq = zeros(3,1, 'like', R0);
% k = zeros(3,1, 'like', R0);
% g = zeros(3,1, 'like', R0);
% dT_q  = zeros(4, 'like', R0);

%% TEST with symbols
[RJ_ql,RL_ql,rJ_ql,rL_ql,k_ql,g_ql]  = dT_parser(dTJ, dTL, robot, ind_l, R0);

J0i_dq = dBi0_q(:,:,ind_i)*P0;
rp = rL(1:3,ind_i);
%Manipulator Jacobian
%Iterate through all "previous" joints
joints_num=0;
for j=1:ind_i
    %If joint is not fixed
    if robot.joints(j).type~=0

        if robot.con.branch(ind_i,j)==1
          
            dBij_dl = [zeros(3,6); SkewSym(rL_ql(1:3,j)-rL_ql(1:3,ind_i)), zeros(3)];

            if robot.joints(j).type==1
                %Revolute joint
                dpm_dl=[k_ql(1:3,j);cross(k_ql(1:3,j),g(1:3,j))+cross(k(1:3,j),g_ql(1:3,j))];
            elseif robot.joints(j).type==2
                %Prismatic joint
                dpm_dl=[zeros(3,1);k_ql(1:3,j)];
            elseif robot.joints(j).type==0
                %Fixed joint
                dpm_dl=zeros(6,1);
            end

            Jmj_dq(1:6,robot.joints(j).q_id) = dBij_dl*pm(1:6,j)+ [eye(3),zeros(3,3);SkewSym(rL(1:3,j)-rp),eye(3)]*dpm_dl;
        else
            Jmj_dq(1:6,robot.joints(j).q_id)=zeros(6,1);
        end
        joints_num=joints_num+1;
    end
end

end

