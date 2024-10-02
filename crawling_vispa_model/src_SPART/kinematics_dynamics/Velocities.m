function [t0,tL]=Velocities(Bij,Bi0,P0,pm,u0,um,robot)
import casadi.*
% Computes the operational-space velocities of the multibody system.
%
% [t0,tL]=Velocities(Bij,Bi0,P0,pm,u0,um,robot)
%
% :parameters:
%   * Bij -- Twist-propagation matrix (for manipulator i>0 and j>0) -- as a [6x6xn] matrix.
%   * Bi0 -- Twist-propagation matrix (for i>0 and j=0) -- as a [6x6xn] matrix.
%   * P0 -- Base-link twist-propagation "vector" -- as a [6x6] matrix.
%   * pm -- Manipulator twist-propagation "vector" -- as a [6xn] matrix.
%   * u0 -- Base-link velocities [\omega,rdot]. The angular velocity is projected in the body-fixed CCS, while the linear velocity is projected in the inertial CCS -- [6x1].
%   * um -- Joint velocities -- [n_qx1].
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return:
%   * t0 -- Base-link twist [\omega,rdot], projected in the inertial CCS -- as a [6x1] matrix.
%   * tL -- Manipulator twist [\omega,rdot], projected in the inertial CCS -- as a [6xn] matrix.
%
% Use :func:`src.kinematics_dynamics.DiffKinematics` to compute ``Bij``, ``Bi0``, ``P0``, and ``pm``.
%
% See also: :func:`src.kinematics_dynamics.Jacob`


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

%--- Number of links and joints ---%
n=robot.n_links_joints;

%Pre-allocate
tL=SX.zeros(6,n);

%Base-link
t0=P0*u0;

%Forward recursion to obtain the twist
for i=1:n
    
    if robot.joints(i).parent_link==0
        %First link
        tL(1:6,i)=Bi0{i}*t0;
    else
        %Rest of the links
        tL(1:6,i)=Bij{i}{1,i-1,1}*tL(1:6,i-1);%Bij(1:6,1:6,i,i-1)*tL(1:6,i-1);
    end
    %Add joint contribution
    if robot.joints(i).type~=0
        tL(:,i)=tL(:,i)+pm(:,i)*um(robot.joints(i).q_id);
    end
    
end

end