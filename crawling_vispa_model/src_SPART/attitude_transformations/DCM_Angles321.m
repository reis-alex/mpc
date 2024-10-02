function Angles = DCM_Angles321(DCM) %#codegen
%Convert the DCM to Euler angles (321 sequence), x-phi, y-theta, z-psi.
%
% This solution generates angle theta between -pi/2 and pi/2,
% and angles phi and psi between -pi and pi.

Angles=[ atan2(DCM(2,3),DCM(3,3));
        asin(-DCM(1,3));
        atan2(DCM(1,2),DCM(1,1))];


end