function [C0_eta, C0q_eta, Cq0_eta, Cq_eta,  C0f, Cf0, Cfq, Cqf, Cf] = CIM_flex(r0,R0,RJ,RL,rJ,rL,k,g, dTJ, dTL,dTJ0,dTL0, Bij,Bi0,P0,pm,q,q0, q_dot,q0_dot, etadot,Jf0dot, Jfqdot,robot)
% Computes the flex components of the Generalized Convective Inertia Matrix C of the multibody system.
%
%
% [Cf0, Cfm, Cf] = CIM_flex(RL,Jf0dot, Jfqdot,robot)
%
% :parameters:
%   * RL -- Links CCS 3x3 rotation matrices with respect to the inertial CCS -- as a [3x3xn] matrix.
%   * Jf0dot -- Base-link connection points Jacobian time-derivative -- as a [6x6xn_f] matrix.
%   * Jfqdot -- Manipulator connection points Jacobian time-derivative -- as a [6xn_qxn_f] matrix.
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return:
%   * Cf0 -- Flexible DOF -- Base-link -- coupling convective inertia matrix -- as a [n_mx6] matrix.
%   * Cfq -- Flexible DOF -- manipulator coupling  convective inertia matrix -- as a [n_mxn_q] matrix.
%   * Cf -- Flexible DOF convective inertia matrix -- as a [n_mxn_m] matrix.
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

%--- C matrix ---%
%Pre-allocate Cf0 and C0f

C0_eta  = zeros(6,"like", Jf0dot);
C0q_eta = zeros(6, robot.n_q,"like", Jf0dot);
Cq0_eta = zeros(robot.n_q, 6,"like", Jf0dot);
Cq_eta  = zeros(robot.n_q,"like", Jf0dot);

%--- C matrix ---%
%Pre-allocate Cf0 and C0f
C0f=zeros(6,n_m,'like',Jf0dot);
Cf0=zeros(n_m,6,'like',Jf0dot);

%Pre-allocate Cfq
Cqf=zeros(n_q,n_m,'like',Jfqdot);
Cfq=zeros(n_m,n_q,'like',Jfqdot);

%Cf
Cf=zeros(n_m,n_m,'like',Jf0dot);

% loop on the flexible links for Hflex_dot
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
    R_CCS_link=[RL(:,:,ind_link)' zeros(3,3,'like',RL); zeros(3,3,'like',RL) RL(:,:,ind_link)'];
    
    Cf0(ind_DOF,:)=robot.flexs.Li(ind_DOF,[4:6 1:3])*Jf0dot(:,:,ind_l);
    C0f(:,ind_DOF)=Cf0(ind_DOF,:)';

    Cfq(ind_DOF,:)=robot.flexs.Li(ind_DOF,[4:6 1:3])*Jfqdot(:,:,ind_l);
    Cqf(:,ind_DOF)=Cfq(ind_DOF,:)';
    Cf(ind_DOF,ind_DOF)= 2*diag(robot.flexs.pulse(1,ind_DOF).*robot.flexs.damp(1,ind_DOF));
    
end
[dHf0_dq0, dHf0_dq, dHfq_dq0, dHfq_dq, dHf_dq0, dHf_dq] = dH_dq_flex(r0,R0,P0,RJ,RL,rJ,rL,k,g, dTJ, dTL,dTJ0,dTL0, Bij,Bi0,pm,q,q0,Jf0dot, Jfqdot,robot);
% Full CC matrix computation with mass matrix derivatives

for i =1:6
    
    C0_eta(i,:)  = -0.5*(etadot'*dHf0_dq0(:,:,i)); % eq 37a
    C0q_eta(i,:) = -0.5*(etadot'*dHfq_dq0(:,:,i)); % eq 37b

    C0f(i,:)     = C0f(i,:) - 0.5*(q0_dot'*dHf0_dq0(:,:,i)' + q_dot'*dHfq_dq0(:,:,i)'); % eq 37c
end

for i = 1:n_q
    Cq0_eta(i,:) = -0.5*(etadot'*dHf0_dq(:,:,i));  % eq 37d
    Cq_eta(i,:)  = -0.5*(etadot'*dHfq_dq(:,:,i));  % eq 37e
    Cqf(i,:)     =  Cqf(i,:)  -0.5*(q0_dot'*dHf0_dq(:,:,i)'+ q_dot'*dHfq_dq(:,:,i)'); %eq 37f
end