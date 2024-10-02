function DCM = Angles123_DCM(Angles) %#codegen
% Convert the Euler angles (123 sequence), x-phi, y-theta, z-psi to DCM.
%
% DCM = Angles123_DCM(Angles)
%
% :parameters: 
%   * Angles -- Euler angles [x-phi, y-theta, z-psi] -- [3x1].
%
% :return: 
%   * DCM -- Direction Cosine Matrix -- [3x3].
%
% See also: :func:`src.attitude_transformations.Angles321_DCM` and :func:`src.attitude_transformations.DCM_Angles321`.

%Define angles
phi = Angles(1); %Rotation along x
theta = Angles(2); %Rotation along y
psi = Angles(3); %Rotation along z

%Create individual DCM
C1=[1 0 0;0 cos(phi) sin(phi);0 -sin(phi) cos(phi)]; 
C2=[cos(theta) 0 -sin(theta);0 1 0;sin(theta) 0 cos(theta)]; 
C3=[cos(psi) sin(psi) 0;-sin(psi) cos(psi) 0;0 0 1];

%Compute DCM for 123 sequence
DCM = C3*C2*C1;

end