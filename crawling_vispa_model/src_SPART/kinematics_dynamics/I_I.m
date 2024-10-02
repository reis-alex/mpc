function [I0,Im]=I_I(R0,RL,robot)

import casadi.*
% Projects the link inertias in the inertial CCS.
%
% [I0,Im]=I_I(R0,RL,robot)
%
% :parameters: 
%   * R0 -- Rotation matrix from the base-link CCS to the inertial CCS -- [3x3].
%   * RL -- Links CCS 3x3 rotation matrices with respect to the inertial CCS -- as a [3x3xn] matrix.
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return: 
%   * I0 -- Base-link inertia matrix, projected in the inertial CCS -- as a [3x3] matrix.
%   * Im -- Links inertia matrices, projected in the inertial CCS -- as a [3x3xn] matrix.
%
% See also: :func:`src.kinematics_dynamics.MCB`. 

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

%Base-link inertia
I0 = R0*robot.base_link.inertia*R0';

%Pre-allocate inertias
Im=SX.sym('Im',3,3,robot.n_links_joints);

%Inertias of the links
for i=1:(robot.n_links_joints)
    Im{i}=RL{i}*robot.links(i).inertia*RL{i}';
end

end
