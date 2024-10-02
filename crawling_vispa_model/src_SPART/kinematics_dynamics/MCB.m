function [M0_tilde,Mm_tilde]=MCB(I0,Im,Bij,Bi0,robot)

import casadi.*
% Computes the Mass Composite Body Matrix (MCB) of the multibody system.
%
% [M0_tilde,Mm_tilde]=MCB(I0,Im,Bij,Bi0,robot)
%
% :parameters: 
%   * I0 -- Base-link inertia matrix, projected in the inertial CCS -- as a [3x3] matrix.
%   * Im -- Links inertia matrices, projected in the inertial CCS -- as a [3x3xn] matrix.
%   * Bij -- Twist-propagation matrix (for manipulator i>0 and j>0) -- as a [6x6xn] matrix.
%   * Bi0 -- Twist-propagation matrix (for i>0 and j=0) -- as a [6x6xn] matrix.
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return: 
%   * M0_tilde -- Base-link mass composite body matrix -- as a [6x6] matrix .
%   * Mm_tilde -- Manipulator mass composite body matrix -- as a [6x6xn] matrix.
%
% See also: :func:`src.kinematics_dynamics.I_I`. 

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

%Number of links and Joints
n=robot.n_links_joints;

%Pre-allocate
Mm_tilde=SX.sym('Mm_t',6,6,n);

%Backwards recursion
for i=n:-1:1
    %Initialize M tilde
    Mm_tilde{i}=[Im{i},zeros(3,3);zeros(3,3),robot.links(i).mass*eye(3)];
    %Add children contributions
    child=find(robot.con.child(:,i))';
    for j=1:length(child)
        %Mm_tilde{i}=Mm_tilde{i}+Bij(1:6,1:6,child(j),i)'*Mm_tilde(1:6,1:6,child(j))*Bij(1:6,1:6,child(j),i);
        Mm_tilde{i}=Mm_tilde{i}+Bij{i}{1,child(j),1}'*Mm_tilde{child(j)}*Bij{i}{1,child(j),1};
    end
end

%Base-link M tilde
M0_tilde=[I0,zeros(3,3);zeros(3,3),robot.base_link.mass*eye(3)];
%Add children contributions
child=find(robot.con.child_base)';
for j=1:length(child)
    M0_tilde=M0_tilde+Bi0{child(j)}'*Mm_tilde{child(j)}*Bi0{child(j)};
end


end