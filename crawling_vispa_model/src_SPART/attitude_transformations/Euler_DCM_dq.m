function DCM = Euler_DCM_dq(e, alpha)
%Provides the Direction Cosine Matrix (DCM) from a Euler axis e=[e(1),e(2),e(3)]
%and angle alpha.

%Create quaternion
% q=0.5*[  e*cos(alpha/2);
%           -sin(alpha/2)];

%Convert quaternion to DCM joint derivative
%DCM = quat_DCM(q);
% DCM =[                   0,  -e(2)*(cos(alpha) - 1), -e(3)*(cos(alpha)-1);
%     -e(2)*(cos(alpha) - 1), 2*e(1)*(cos(alpha) - 1),           sin(alpha);
%     -e(3)*(cos(alpha) - 1),             -sin(alpha), 2*e(1)*(cos(alpha)- 1)];

DCM =[           -sin(alpha)*(e(2)^2 + e(3)^2),  e(3)*cos(alpha) + e(1)*e(2)*sin(alpha),   e(1)*e(3)*sin(alpha) - e(2)*cos(alpha);
        e(1)*e(2)*sin(alpha) - e(3)*cos(alpha),           -sin(alpha)*(e(1)^2 + e(3)^2),   e(1)*cos(alpha) + e(2)*e(3)*sin(alpha);
        e(2)*cos(alpha) + e(1)*e(3)*sin(alpha),  e(2)*e(3)*sin(alpha) - e(1)*cos(alpha),           -sin(alpha)*(e(1)^2 + e(2)^2)];
end

