function [C]=CIM_NOC(N,Ndot,t0,tL,I0,Im,robot)
import casadi.*
% Computes the generalized Convective Inertia Matrix C of the multibody system.
%
% This function is LESS EFFICENT that the :func:`src.kinematics_dynamics.CIM` function, as it doesn't use the recursive algorithms. 
%
% [C]=CIM_NOC(N,Ndot,t0,tL,I0,Im,robot)
%
% :parameters: 
%   * N -- Natural Orthogonal Complement (NOC) matrix -- a [(6+6*n)x(6+n_q)] matrix.
%   * Ndot -- Natural Orthogonal Complement (NOC) matrix time-derivative -- as a [(6+6*n)x(6+n_q)] matrix.
%   * t0 -- Base-link twist [\omega,rdot], projected in the inertial CCS -- as a [6x1] matrix.
%   * tL -- Manipulator twist [\omega,rdot], projected in the inertial CCS -- as a [6xn] matrix.
%   * I0 -- Base-link inertia matrix, projected in the inertial CCS -- as a [3x3] matrix.
%   * Im -- Links inertia matrices, projected in the inertial CCS -- as a [3x3xn] matrix.
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return: 
%   C -> Generalized Convective Inertia Matrix -- as a [(6+n_q)x(6+n_q)] matrix.
%
%  See also: :func:`src.kinematics_dynamics.CIM_NOC`.


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

%=== Code ===%

%--- Generalized mass matrix ---%

%Pre-allocate M
M=zeros(6+6*robot.n_links_joints,6+6*robot.n_links_joints,'like',N);

%Base-link contribution
M(1:6,1:6)=[I0(1:3,1:3),zeros(3,3);zeros(3,3),robot.base_link.mass*eye(3)];
%Manipulator contribution
for i=1:robot.n_links_joints
    M(6+6*i-5:6+6*i,6+6*i-5:6+6*i)=[Im{i},zeros(3,3);zeros(3,3),robot.links(i).mass*eye(3)];
end

%--- Omega ---%
%Base-link Omega
Omega0=[SkewSym(t0(1:3)), zeros(3,3);
    zeros(3,3), SkewSym(t0(1:3))];

%Pre-allocate Omega
Omega=SX.sym('Omega',6,6,robot.n_links_joints);

%Compute Omega
for j=1:i
    Omega{j}=[SkewSym(tL(1:3,j)), zeros(3,3);
        zeros(3,3), SkewSym(tL(1:3,j))];
end

%--- Mdot ---%
%Pre-allocate
Mdot=SX.zeros(6+6*robot.n_links_joints,6+6*robot.n_links_joints);

%Base-link Mdot
Mdot(1:6,1:6)=[Omega0(1:3,1:3)*I0, zeros(3,3); zeros(3,3), zeros(3,3)];
%Compute Mdot
for i=1:robot.n_links_joints
    Mdot(6+6*i-5:6+6*i,6+6*i-5:6+6*i)=[Omega{i}(1:3,1:3)*Im{i}, zeros(3,3); zeros(3,3), zeros(3,3)];
end

%--- Convective Inertia Matrix ---%
C=N'*(M*Ndot+Mdot*N);

end