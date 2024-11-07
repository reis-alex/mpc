function [Hf0, Hfq, Hf] = GIM_flex(RL,Jf0, Jfm,robot)
% Computes the flex components of the Generalized Inertia Matrix (GIM) H of the multibody vehicle.
%
%
% [Hf0, Hfm, Hf] = GIM_flex(RL,Jf0, Jfm,robot)
%
% :parameters:
%   * RL -- Links CCS 3x3 rotation matrices with respect to the inertial CCS -- as a [3x3xn] matrix.
%   * Jf0 -- Base-link Jacobian time-derivative -- as a [6x6xn_f] matrix.
%   * Jfm -- Manipulator Jacobian time-derivative -- as a [6xn_qxn_f] matrix.
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return:
%   * Hf0 -- Flexible DOF -- Base-link -- coupling inertia matrix -- as a [n_mx6] matrix.
%   * Hfm -- Flexible DOF -- manipulator coupling inertia matrix -- as a [n_mxn_q] matrix.
%   * Hf -- Flexible DOF inertia matrix -- as a [n_mxn_m] matrix.
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

%--- Number of links and Joints ---%
n_q=robot.n_q;

%--- Number of links  ---%
n_f=robot.n_links_flex;

%--- Number of DOF  ---%
n_m=robot.flexs.nb_mode;


%--- H matrix ---%
%Pre-allocate Hf0
Hf0=zeros(n_m,6,'like',Jf0);

%Pre-allocate Hfm
Hfq=zeros(n_m,n_q,'like',Jfm);

%Hf
Hf=eye(n_m,n_m,'like',Jf0);


% loop on the flexible links
for ind_l=1:n_f
    %Indice of the flexible DOF
    if ind_l==n_f
        ind_DOF=robot.flexs.link_mode_id(ind_l):robot.flexs.nb_mode;
    else
        ind_DOF=robot.flexs.link_mode_id(ind_l):robot.flexs.link_mode_id(ind_l+1)-1;
    end
    %Indice of the flexible links
    ind_link=robot.flexs.link_id(ind_l);
    
    %rotation matrix from Inertia CCS to the link CCS
%     R_CCS_link=[RL(:,:,ind_link)' zeros(3,3,'like',RL); zeros(3,3,'like',RL) RL(:,:,ind_link)'];
%     
%     Hf0(ind_DOF,:)=robot.flexs.Li(ind_DOF,[4:6 1:3])*R_CCS_link*Jf0(:,:,ind_l);
%     Hfq(ind_DOF,:)=robot.flexs.Li(ind_DOF,[4:6 1:3])*R_CCS_link*Jfm(:,:,ind_l);
     Hf0(ind_DOF,:)=robot.flexs.Li(ind_DOF,[4:6 1:3])*Jf0(:,:,ind_l);
     Hfq(ind_DOF,:)=robot.flexs.Li(ind_DOF,[4:6 1:3])*Jfm(:,:,ind_l);
    
end