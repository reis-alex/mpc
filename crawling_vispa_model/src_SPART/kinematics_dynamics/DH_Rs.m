function [R,s] = DH_Rs(DH,qm,type) %#codegen
% This function computes the rotation matrix R and the translation vector s
% between two joints given their Denavit-Hartenber (DH) parameters.
%
% Inputs:
%   DH -> Denavit-Hartenberg parameters.
%       DH.d -> Distance between joint origins along the joint z-axis.
%       DH.theta -> Rotation between x-axis along the joint z-axis.
%       DH.alpha -> Rotation between z-axis along the x-axis.
%       DH.a -> Distance between the common normal between the z-axis.
%   qm -> Joint variable.
%   type -> type==0 for revolute joint or type==1 for prismatic joints.
%
% Outputs:
%   R -> Rotation 3x3 matrix.
%   s -> Translation 3x1 vector.

%=== LICENSE ===%

%=== CODE ===%

%Assign the d and theta variable depending on joint type
if type==0
    %Revolute joint
    theta=DH.theta+qm;
    d=DH.d;
else
    %Prismatic joint
    theta=DH.theta;
    d=DH.d+qm;
end
    

%Rotation matrix
R = [
    cos(theta), -sin(theta)*cos(DH.alpha), sin(theta)*sin(DH.alpha);
    sin(theta), cos(theta)*cos(DH.alpha),  -cos(theta)*sin(DH.alpha);
    0,              sin(DH.alpha),                 cos(DH.alpha)
    ];      

%Translation vector
s = [DH.a*cos(theta), DH.a*sin(theta), d]';

end