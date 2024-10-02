function [H0, H0m, Hm] = GIM(M0_tilde,Mm_tilde,Bij,Bi0,P0,pm,robot)
import casadi.*
% Computes the Generalized Inertia Matrix (GIM) H of the multibody vehicle.
%
% This function uses a recursive algorithm.
%
% [H0, H0m, Hm] = GIM(M0_tilde,Mm_tilde,Bij,Bi0,P0,pm,robot)
%
% :parameters: 
%   * M0_tilde -- Base-link mass composite body matrix -- as a [6x6] matrix .
%   * Mm_tilde -- Manipulator mass composite body matrix -- as a [6x6xn] matrix.
%   * Bij -- Twist-propagation matrix (for manipulator i>0 and j>0) -- as a [6x6xn] matrix.
%   * Bi0 -- Twist-propagation matrix (for i>0 and j=0) -- as a [6x6xn] matrix.
%   * P0 -- Base-link twist-propagation "vector" -- as a [6x6] matrix.
%   * pm -- Manipulator twist-propagation "vector" -- as a [6xn] matrix.
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return: 
%   * H0 -- Base-link inertia matrix -- as a [6x6] matrix.
%   * H0m -- Base-link -- manipulator coupling inertia matrix -- as a [6xn_q] matrix.
%   * Hm -- Manipulator inertia matrix -- as a [n_qxn_q] matrix.
%   
% To obtain the full generalized inertia matrix H:
%
% .. code-block:: matlab
%   
%   %Compute H
%   [H0, H0m, Hm] = GIM(M0_tilde,Mm_tilde,Bij,Bi0,P0,pm,robot);
%   H=[H0,H0m;H0m';Hm];
%
% See also: :func:`src.kinematics_dynamics.CIM`.


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

%--- Number of links and Joints ---%
n_q=robot.n_q;
n=robot.n_links_joints;

%--- H matrix ---%

%Base-link inertia matrix
H0 = P0'*M0_tilde*P0;

%Pre-allocate Hm
Hm=SX.zeros(n_q,n_q);

%Manipulator inertia matrix Hm
for j=1:n
    for i=j:n
        if robot.joints(i).type~=0 && robot.joints(j).type~=0
            Hm(robot.joints(i).q_id,robot.joints(j).q_id)=pm(1:6,i)'*Mm_tilde{i}*Bij{i}{1,j,1}*pm(1:6,j);
            Hm(robot.joints(j).q_id,robot.joints(i).q_id)=Hm(robot.joints(i).q_id,robot.joints(j).q_id);
        end
    end
end

%Pre-allocate H0m
H0m=zeros(6,n_q,'like',M0_tilde);

%Coupling inertia matrix
for i=1:n
    if robot.joints(i).type~=0
        H0m(1:6,robot.joints(i).q_id)=(pm(1:6,i)'*Mm_tilde{i}*Bi0{i}*P0)';
    end
end

end