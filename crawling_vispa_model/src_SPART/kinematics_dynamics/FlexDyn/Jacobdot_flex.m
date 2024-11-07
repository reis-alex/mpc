function [Jf0dot, Jfmdot] = Jacobdot_flex(rJ,r0,t0,rL,tL,P0,pm,Jf0, Jfm,u0,um,robot)
% Computes the geometric Jacobian time-derivative of the connection points of the n_f flexible links.
%
%  [Jf0dot, Jfmdot] = Jacobdot_flex(rJ,r0,t0,rL,tL,P0,pm,Jf0, Jfm,u0,um,robot)
%
% :parameters:
%   * rJ -- Positions of the joints, projected in the inertial CCS -- as a [3xn] matrix.
%   * r0 -- Position of the base-link center-of-mass with respect to the origin of the inertial frame, projected in the inertial CCS -- [3x1].
%   * t0 -- Base-link twist [\omega,rdot], projected in the inertial CCS -- as a [6x1] matrix.
%   * rL -- Positions of the links, projected in the inertial CCS -- as a [3xn] matrix.
%   * tL -- Manipulator twist [\omega,rdot], projected in the inertial CCS -- as a [6xn] matrix.
%   * P0 -- Base-link twist-propagation "vector" -- as a [6x6] matrix.
%   * pm -- Manipulator twist-propagation "vector" -- as a [6xn] matrix.
%   * Jf0 -- Base-link connection points Jacobian -- as a [6x6xn_f] matrix.
%   * Jfm -- Manipulator connection points Jacobian -- as a [6xn_qxn_f] matrix.
%   * u0 -- Base-link velocities [\omega,rdot]. The angular velocity is projected in the body-fixed CCS, while the linear velocity is projected in the inertial CCS -- [6x1].
%   * um -- Joint velocities -- [n_qx1].
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return: 
%   * Jf0dot -- Base-link connection points Jacobian time-derivative -- as a [6x6xn_f] matrix.
%   * Jfmdot -- Manipulator connection points Jacobian time-derivative -- as a [6xn_qxn_f] matrix.
%
% See also: :func:`src.kinematics_dynamics.Jacob_flex`.

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
n=robot.n_links_flex;

Jf0dot=zeros(6,6,n,'like',r0);
Jfmdot=zeros(6,robot.n_q,n,'like',r0);

for ind_l=1:n
    %The joint is the connection point of the flexible link
    rp=rL(:,robot.flexs.joint_id(ind_l));
    %Twist of the connection point
    tp=Jf0(:,:,ind_l)*u0+Jfm(:,:,ind_l)*um;
    
   [Jf0dot(:,:,ind_l), Jfmdot(:,:,ind_l)]=...
       Jacobdot(rp,tp,r0,t0,rL,tL,P0,pm,robot.flexs.link_id(ind_l),robot);

end

