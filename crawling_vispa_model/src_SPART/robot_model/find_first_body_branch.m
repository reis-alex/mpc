function [ind_init] = find_first_body_branch(j,robot)
%FIND_FIRST_BODY_BRANCH Summary of this function goes here
%   Detailed explanation goes here
ind_init = 152;
child_con_base = find(robot.con.child_base);
for i = 1:length(child_con_base)
    if robot.con.branch(j,child_con_base(i)) == 1
        ind_init = child_con_base(i);
        break
    end
end

