function [Ndot] = NOCdot(r0,t0,rL,tL,P0,pm,robot)
% Computes the Natural Orthogonal Complement (NOC) matrix time-derivative.
%
% [Ndot] = NOCdot(r0,t0,rL,tL,P0,pm,robot)
%
% :parameters: 
%   * r0 -- Position of the base-link center-of-mass with respect to the origin of the inertial frame, projected in the inertial CCS -- [3x1].
%   * t0 -- Base-link twist [\omega,rdot], projected in the inertial CCS -- as a [6x1] matrix.
%   * rL -- Positions of the links, projected in the inertial CCS -- as a [3xn] matrix.
%   * tL -- Manipulator twist [\omega,rdot], projected in the inertial CCS -- as a [6xn] matrix.
%   * P0 -- Base-link twist-propagation "vector" -- as a [6x6] matrix.
%   * pm -- Manipulator twist-propagation "vector" -- as a [6xn] matrix.
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return: 
%   * Ndot -- Natural Orthogonal Complement (NOC) matrix time-derivative -- as a [(6+6*n)x(6+n_q)] matrix.
%
% Examples:
%
%   To compute the operational-space accelerations of all links:
%
% .. code-block:: matlab
%	
%   %Compute NOC
%   [N] = NOC(r0,rL,P0,pm,robot)
%   %Compute NOC time-derivative
%   [Ndot] = NOCdot(r0,t0,rL,tL,P0,pm,robot)
%   %Twist time-derivatives of all the links
%   tdot=N*[u0dot;umdot]+Ndot*[u0;um];
%   %Twist time-derivative of the base-link
%   t0dot=tdot(1:6,1);
%   %Twist time-derivative of the ith link
%   i=2;
%   tidot=tdot(6*i:6+6*i,1);
%
% See also: :func:`src.kinematics_dynamics.Jacobdot` and :func:`src.kinematics_dynamics.NOC`. 


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

%=== Code ===%

%Pre-allocate NOCdot
Ndot=zeros(6+6*robot.n_links_joints,6+robot.n_q,'like',r0);

%Compute the NOC time derivative matrix by using the Jacobians time derivative

%Base-link contribution
Omega0=[SkewSym(t0(1:3)), zeros(3,3);
        zeros(3,3), zeros(3,3)];
Ndot(1:6,1:6+robot.n_q)=[Omega0*P0,zeros(6,robot.n_q)];

%Manipulator contribution
for i=1:robot.n_links_joints
    [J0doti, Jmdoti]=Jacobdot(rL(1:3,i),tL(1:6,i),r0,t0,rL,tL,P0,pm,i,robot);
    Ndot(6+6*i-5:6+6*i,1:6+robot.n_q)=[J0doti,Jmdoti];
end


end