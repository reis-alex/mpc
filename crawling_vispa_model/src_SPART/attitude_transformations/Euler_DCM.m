function DCM = Euler_DCM(e,alpha) %#codegen
%Provides the Direction Cosine Matrix (DCM) from a Euler axis e=[e1,e2,e3]
%and angle alpha.

%Create quaternion
q=[ e*sin(alpha/2);
    cos(alpha/2)];

%Convert quaternion to DCM
DCM = quat_DCM(q);

end