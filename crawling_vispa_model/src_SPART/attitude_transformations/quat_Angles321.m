function Angles = quat_Angles321(q) %#codegen
%Convert a quaternion to Euler angles (321 sequence), x-phi, y-theta, z-psi.

%Quaternion to DCM
DCM = quat_DCM(q);

%DCM to Angles (using 321 sequence)
Angles = DCM_Angles321(DCM);
    
end