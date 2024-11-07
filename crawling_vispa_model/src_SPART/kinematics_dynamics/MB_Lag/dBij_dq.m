function [dBij_ql] = dBij_dq(ind_i,ind_j,ind_q_l,R0,drdq,robot)
%UNTITLED Summary of this function goes here
%DBI0_QL Summary of this function goes here
%dBi0 6x6
%   Detailed explanation goes here

dBij_ql = zeros(6,6, 'like', R0);
drdqi   = zeros(3,1, 'like', R0);
drdqj   = zeros(3,1, 'like', R0);
if robot.con.branch(ind_i,ind_j)==1
    
    drdqi = robot.con.branch(ind_i,ind_q_l)*drdq(:,ind_i);
    drdqj = robot.con.branch(ind_j,ind_q_l)*drdq(:,ind_j);
    drdql = drdqj - drdqi;
    dBij_ql(:,:) =  [zeros(3)       zeros(3);
                     SkewSym(drdql)  zeros(3)];

end

end