function [I0_dq0,Im_dq0]=dI_dq0(dTJ0,dTL0, R0, dR0_dq0, RL,robot)
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
n_base = 6;
%Pre-allocate inertias
I0_dq0 = zeros(3,3,n_base, "like", R0);
Im_dq0=zeros(3,3,robot.n_links_joints,n_base,'like',R0);

% dT0_dq0 = dT0_diffq0(q0);
% dR0_dq0 = dT0_dq0(1:3,1:3,:);
for coord_i = 1:n_base
    
    %Base-link inertia
    I0_dq0(:,:,coord_i) = dR0_dq0(:,:,coord_i)*robot.base_link.inertia*R0.'...
                          + R0*robot.base_link.inertia*dR0_dq0(:,:,coord_i).';
    [~,dRL_dq0,~,~,~,~]  = dT_parser(dTJ0, dTL0, robot, coord_i, R0);
    %Inertias of the links
    for i=1:(robot.n_links_joints)
        Im_dq0(1:3,1:3,i,coord_i)=dRL_dq0(:,:,i)*robot.links(i).inertia*RL(1:3,1:3,i).'...
                                  +RL(1:3,1:3,i)*robot.links(i).inertia*dRL_dq0(:,:,i).' ;
    end

end

