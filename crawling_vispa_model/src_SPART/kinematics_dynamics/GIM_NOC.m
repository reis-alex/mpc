function [H]=GIM_NOC(N,I0,Im,robot)
% Computes the Generalized Inertia Matrix (GIM) H of the multibody vehicle.
%
% This function is LESS EFFICENT that the :func:`src.kinematics_dynamics.GIM` function, as it doesn't use the recursive algorithms. 
%
% [H]=GIM_NOC(N,I0,Im,robot)
%
% :parameters: 
%   * N -- Natural Orthogonal Complement (NOC) matrix -- a [(6+6*n)x(6+n_q)] matrix.
%   * I0 -- Base-link inertia matrix, projected in the inertial CCS -- as a [3x3] matrix.
%   * Im -- Links inertia matrices, projected in the inertial CCS -- as a [3x3xn] matrix.
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return: 
%
%   * H -- Generalized Inertia Matrix -- as a [(6+n_q)x(6+n_q)] matrix.
%
% See also: :func:`src.kinematics_dynamics.GIM_NOC`.

%=== LICENSE ===%

%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU Lesser General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%=== CODE ===%

%--- generalized mass matrix M ---%

%Pre-allocate M
M=zeros(6+6*robot.n_links_joints,6+6*robot.n_links_joints,'like',N);

%Base contribution
M(1:6,1:6)=[I0(1:3,1:3),zeros(3,3);zeros(3,3),robot.base_link.mass*eye(3)];
%Manipulator Contribution
for i=1:robot.n_links_joints
    M(6+6*i-5:6+6*i,6+6*i-5:6+6*i)=[Im{i},zeros(3,3);zeros(3,3),robot.links(i).mass*eye(3)];
end

%--- Generalized Inertia Matrix H ---%
H=N'*M*N;


end