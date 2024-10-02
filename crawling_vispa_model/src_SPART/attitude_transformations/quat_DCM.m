function DCM = quat_DCM(q) %#codegen
%Provides the Direction Cosine Matrix (DCM) from a quaterionion (q)
%
% q = [q1;q2;q3;q4] -> With q4 being the scalar part of the quaternion.


DCM = [ 1-2*(q(2)^2+q(3)^2), 2*(q(1)*q(2)+q(3)*q(4)), 2*(q(1)*q(3)-q(2)*q(4));
        2*(q(2)*q(1)-q(3)*q(4)), 1-2*(q(1)^2+q(3)^2), 2*(q(2)*q(3)+q(1)*q(4));
        2*(q(3)*q(1)+q(2)*q(4)), 2*(q(3)*q(2)-q(1)*q(4)), 1-2*(q(1)^2+q(2)^2)];
    
end