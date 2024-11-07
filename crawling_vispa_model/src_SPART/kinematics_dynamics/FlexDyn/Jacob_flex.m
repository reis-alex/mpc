function [Jf0, Jfm] = Jacob_flex(rJ,r0,rL,P0,pm,robot)
% Computes the geometric Jacobian of the connection points of the n_f flexible links.
%
%  [Jf0, Jfm] = Jacob_flex(rJ,r0,rL,P0,pm,robot)
%
% :parameters:
%   * rJ -- Positions of the joints, projected in the inertial CCS -- as a [3xn] matrix.
%   * r0 -- Position of the base-link, projected in the inertial CCS -- [3x1].
%   * rL -- Positions of the links, projected in the inertial CCS -- as a [3xn] matrix.
%   * P0 -- Base-link twist-propagation "vector" -- as a [6x6] matrix.
%   * pm -- Manipulator twist-propagation "vector" -- as a [6xn] matrix.
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return: 
%   * Jf0 -- Base-link connection points Jacobian -- as a [6x6xn_f] matrix.
%   * Jfm -- Manipulator connection points Jacobian -- as a [6xn_qxn_f] matrix.
%
% See also: :func:`src.kinematics_dynamics.Jacob`.

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

Jf0=zeros(6,6,n,'like',rJ);
Jfm=zeros(6,robot.n_q,n,'like',rJ);

for ind_l=1:n
    %The joint is the connection point of the flexible link
    rp=rL(:,robot.flexs.joint_id(ind_l));
    %Compute Jacobians of the connection point
    [Jf0(:,:,ind_l), Jfm(:,:,ind_l)]=Jacob(rp,r0,rL,P0,pm,robot.flexs.link_id(ind_l),robot);

end

