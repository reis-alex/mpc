function [Bij_dot] = Bij_dot(tL,i,j)
%Return time derivative of Bij propagation matrix
%   t torseur cinematique [omega_i, rdot_i] concatenés (6xnq)

        if robot.con.branch(i,j)==1
            %Links are in the same branch
            Bij_dot(1:6,1:6,i,j)=[zeros(3,3), zeros(3,3); SkewSym(tL(4:6,j)-tL(4:6,i)), zeros(3,3)];
        else
            %Links are not in the same branch
            Bij_dot(1:6,1:6,i,j)=zeros(6,6);
        end


end

