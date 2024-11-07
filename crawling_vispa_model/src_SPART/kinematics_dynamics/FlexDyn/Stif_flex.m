function [Kf] = Stif_flex(robot)
% Computes the stiffness matrix of the flexible DOF of the multi-body
% systems
%
% [Kf] = Stif_flex(robot)
%
% :parameters:
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return:
%   * Kf -- Stiffness matrix of the flexible links -- as a [n_mxn_m] matrix.
%

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
n_f=robot.n_links_flex;

%--- Number of DOF  ---%
n_m=robot.flexs.nb_mode;

%Kf
Kf=zeros(n_m,n_m,'like',robot.flexs.pulse);

for ind_l=1:n_f
    %Indice of the flexible DOF
    if ind_l==n_f
        ind_DOF=robot.flexs.link_mode_id(ind_l):robot.flexs.nb_mode;
    else
        ind_DOF=robot.flexs.link_mode_id(ind_l):robot.flexs.link_mode_id(ind_l+1)-1;
    end
    %Indice of the flexible links
    ind_link=robot.flexs.link_id(ind_l);
    
    Kf(ind_DOF,ind_DOF)= diag(robot.flexs.pulse(1,ind_DOF).*robot.flexs.pulse(1,ind_DOF));
    
end