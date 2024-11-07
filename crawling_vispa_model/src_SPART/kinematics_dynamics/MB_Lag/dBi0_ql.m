function [dBi0_ql] = dBi0_ql(r0,R0,RJ,RL,rJ,rL,dTL,dTJ,q,robot)
%DBI0_QL Summary of this function goes here
%dBi0 6x6xnxn
%   Detailed explanation goes here

dBi0_ql = zeros(6,6,robot.n_links_joints, robot.n_q, 'like', R0);
drdq = zeros(3,1, 'like', R0);
for l = 1:robot.n_q
    [RJ_ql,RL_ql,rJ_ql,rL_ql,e_ql,g_ql]  = dT_parser(dTJ, dTL, robot, l, R0);
    for j= 1:robot.n_links_joints
        if  robot.con.branch(j,l) == 1
        
        drdq = rL_ql(:,j);
        dBi0_ql(:,:,j,l) =  [zeros(3)       zeros(3);
                             SkewSym(-drdq)  zeros(3)];
        end
    end
end

