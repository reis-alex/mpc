function [dBij_dq0,dBi0_dq0,dP0_dq0,dpm_dq0]=DiffKinematics_dq0(dTL0, dTJ0, R0,dT0_dq0, k,g,ind_q_0,robot)
% Computes the differential kinematics partial derivative of l-th joint of the multibody system.
%
% [Bij,Bi0,P0,pm]=DiffKinematics(R0,r0,rL,e,g,robot)
% 
% :parameters:
%   * dTJ -- Joints CCS 4x4 homogenous transformation matrices partial
%     derivatives over q with respect to inertial CCS -- as a [4x4xnxn_q]
%   * dTL -- Links CCS 4x4 homogenous transformation matrices partial
%     derivatives over q with respect to inertial CCS -- as a [4x4xnxn_q]
%   * R0 -- Rotation matrix from the base-link CCS to the inertial CCS -- [3x3].
%   * r0 -- Position of the base-link center-of-mass with respect to the origin of the inertial frame, projected in the inertial CCS -- [3x1].
%   * rL -- Positions of the links, projected in the inertial CCS -- as a [3xn] matrix.
%   * e -- Joint rotation/sliding axes, projected in the inertial CCS -- as a [3xn] matrix.
%   * g -- Vector from the origin of the ith joint CCS to the origin of the ith link CCS, projected in the inertial CCS -- as a [3xn] matrix.
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return all partial derivatives over q0 of:
%   * Bij -- Twist-propagation matrix (for manipulator i>0 and j>0) -- as a [6x6xn] matrix.
%   * Bi0 -- Twist-propagation matrix (for i>0 and j=0) -- as a [6x6xn] matrix.
%   * P0 -- Base-link twist-propagation "vector" -- as a [6x6] matrix.
%   * pm -- Manipulator twist-propagation "vector" -- as a [6xn] matrix.
%
% Use :func:`src.kinematics_dynamics.Kinematics` to compute ``rL,e``, and ``g``.
%
% See also: :func:`src.kinematics_dynamics.Kinematics` and :func:`src.kinematics_dynamics.Jacob`. 

%{  
    LICENSE

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}

%=== CODE ===%

%--- Number of links  ---%
n=robot.n_links_joints;

%--- Twist-propagation matrix ---%
[RJ_q0,RL_q0,rJ_q0,rL_q0,k_q0,g_q0]  = dT_parser(dTJ0, dTL0, robot, ind_q_0, R0);

% P0

dP0_dq0=[dT0_dq0(1:3,1:3,ind_q_0),zeros(3,3); zeros(3,3), zeros(3,3)];


%Pre-allocate Bij
dBij_dq0=zeros(6,6,n,n,'like',R0);


%Compute Bij
for j=1:n
    for i=1:n
        if robot.con.branch(i,j)==1
            %Links are in the same branch
            dBij_dq0(1:6,1:6,i,j)=[zeros(3,3), zeros(3,3); SkewSym(rL_q0(1:3,j)-rL_q0(1:3,i)), zeros(3,3)];
        else
            %Links are not in the same branch
            dBij_dq0(1:6,1:6,i,j)=zeros(6,6);
        end
    end
end

%Pre-allocate Bi0
dBi0_dq0=zeros(6,6,n,'like',R0);
dr0_dq0 = dT0_dq0(1:3,4,ind_q_0);
%Compute Bi0
for i=1:n

    dBi0_dq0(1:6,1:6,i)=[zeros(3,3), zeros(3,3); SkewSym(dr0_dq0-rL_q0(1:3,i)), zeros(3,3)];
end

%--- Twist-propagation "vector" ---%

%Pre-allocate pm
dpm_dq0=zeros(6,n,'like',R0);

%Base-link


%Forward recursion to obtain the twist-propagation "vector"
for i=1:n
    if robot.joints(i).type==1
        %Revolute joint
        dpm_dq0(1:6,i)=[k_q0(1:3,i);cross(k_q0(1:3,i),g(1:3,i))+ cross(k(1:3,i),g_q0(1:3,i))];
    elseif robot.joints(i).type==2
        %Prismatic joint
        dpm_dq0(1:6,i)=[zeros(3,1);k_q0(1:3,i)];
    elseif robot.joints(i).type==0
        %Fixed joint
        dpm_dq0(1:6,i)=zeros(6,1);
    end
end

end