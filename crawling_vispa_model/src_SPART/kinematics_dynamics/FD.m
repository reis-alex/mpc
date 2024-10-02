function [u0dot,umdot] = FD(tau0,taum,wF0,wFm,t0,tm,P0,pm,I0,Im,Bij,Bi0,u0,um,robot)
% This function solves the forward dynamics (FD) problem (it obtains the
% acceleration from  forces).
%
% [u0dot,umdot] = FD(tau0,taum,wF0,wFm,t0,tm,P0,pm,I0,Im,Bij,Bi0,u0,um,robot)
%
% :parameters: 
%   * tau0 -- Base-link forces [n,f]. The torque n is projected in the body-fixed CCS, while the force f is projected in the inertial CCS -- [6x1].
%   * taum -- Joint forces/torques -- as a [n_qx1] matrix.
%   * wF0 -- Wrench acting on the base-link center-of-mass [n,f], projected in the inertial CCS -- as a [6x1] matrix.
%   * wFm -- Wrench acting on the links center-of-mass  [n,f], projected in the inertial CCS -- as a [6xn] matrix.
%   * t0 -- Base-link twist [\omega,rdot], projected in the inertial CCS -- as a [6x1] matrix.
%   * tL -- Manipulator twist [\omega,rdot], projected in the inertial CCS -- as a [6xn] matrix.
%   * P0 -- Base-link twist-propagation "vector" -- as a [6x6] matrix.
%   * pm -- Manipulator twist-propagation "vector" -- as a [6xn] matrix.
%   * I0 -- Base-link inertia matrix, projected in the inertial CCS -- as a [3x3] matrix.
%   * Im -- Links inertia matrices, projected in the inertial CCS -- as a [3x3xn] matrix.
%   * Bij -- Twist-propagation matrix (for manipulator i>0 and j>0) -- as a [6x6xn] matrix.
%   * Bi0 -- Twist-propagation matrix (for i>0 and j=0) -- as a [6x6xn] matrix.
%   * u0 -- Base-link velocities [\omega,rdot]. The angular velocity is projected in the body-fixed CCS, while the linear velocity is projected in the inertial CCS -- [6x1].
%   * um -- Joint velocities -- [n_qx1].
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return: 
%   * u0dot -- Base-link accelerations [\omegadot,rddot]. The angular acceleration is projected in a body-fixed CCS, while the linear acceleration is projected in the inertial CCS -- [6x1].
%   * umdot -- Manipulator joint accelerations -- [n_qx1].
%
% See also: :func:`src.kinematics_dynamics.ID` and :func:`src.kinematics_dynamics.I_I`. 

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
n=robot.n_links_joints;
n_q=robot.n_q;

%---- Inverse Dynamics with 0 accelerations ---%
%Recompute Accelerations with u0dot=umdot=0
[t0dot,tmdot]=Accelerations(t0,tm,P0,pm,Bi0,Bij,u0,um,zeros(6,1),zeros(n_q,1),robot);
%Use the inverse dynamics
[tau0_0ddot,taum_0ddot] = ID(wF0,wFm,t0,tm,t0dot,tmdot,P0,pm,I0,Im,Bij,Bi0,robot);

%--- Forward Dynamics ---%

%Initialize solution
phi0=tau0-tau0_0ddot;
phi=taum-taum_0ddot;

%--- M hat, psi hat and psi  ---%
%Pre-allocate
M_hat=zeros(6,6,n,'like',tau0);
psi_hat=zeros(6,n,'like',tau0);
psi=zeros(6,n,'like',tau0);

%Backwards recursion
for i=n:-1:1
    %Initialize
    M_hat(1:6,1:6,i)=[Im(1:3,1:3,i),zeros(3,3);zeros(3,3),robot.links(i).mass*eye(3)];
    %Add children contributions
    for j=find(robot.con.child(:,i))'
        M_hatii=M_hat(1:6,1:6,j)-psi_hat(1:6,j)*psi(1:6,j)';
        M_hat(1:6,1:6,i)=M_hat(1:6,1:6,i)+Bij(1:6,1:6,j,i)'*M_hatii*Bij(1:6,1:6,j,i);
    end
    if robot.joints(i).type==0
        psi_hat(1:6,i)=zeros(6,1);
        psi(1:6,i)=zeros(6,1);
    else
        psi_hat(1:6,i)=M_hat(1:6,1:6,i)*pm(1:6,i);
        psi(1:6,i)=psi_hat(1:6,i)/(pm(1:6,i)'*psi_hat(1:6,i));
    end
end
%Base-link
M_hat0=[I0,zeros(3,3);zeros(3,3),robot.base_link.mass*eye(3)];
%Add children contributions
for j=find(robot.con.child_base)'
    M_hat0ii=M_hat(1:6,1:6,j)-psi_hat(1:6,j)*psi(1:6,j)';
    M_hat0=M_hat0+Bi0(1:6,1:6,j)'*M_hat0ii*Bi0(1:6,1:6,j);
end
psi_hat0=M_hat0*P0;

%--- eta ---%
%Pre-allocate and initialize
eta=zeros(6,n,'like',P0);
phi_hat=zeros(6,n,'like',P0);
phi_tilde=zeros(n_q,'like',P0);

%Backwards recursion
for i=n:-1:1
    %Initialize
    eta(1:6,i)=zeros(6,1);
    %Add children contributions
    for j=find(robot.con.child(:,i))'
        eta(1:6,i)=eta(1:6,i)+Bij(1:6,1:6,j,i)'*(psi(1:6,j)*phi_hat(j)+eta(1:6,j));
    end
    phi_hat(i)=-pm(1:6,i)'*eta(1:6,i);
    if robot.joints(i).type~=0
        phi_hat(i)=phi_hat(i)+phi(robot.joints(i).q_id);
        phi_tilde(robot.joints(i).q_id)=phi_hat(i)/(pm(1:6,i)'*psi_hat(1:6,i));
    end
end
%Base-link
eta0=zeros(6,1);
%Add children contributions
for j=find(robot.con.child_base)'
    eta0=eta0+Bi0(1:6,1:6,j)'*(psi(1:6,j)*phi_hat(j)+eta(1:6,j));
end
phi_hat0=phi0-P0'*eta0;
phi_tilde0=(P0'*psi_hat0)\phi_hat0;


%--- Base-link acceleration ---%
u0dot=phi_tilde0;

%--- Manipulator acceleration (and mu) ---%

%Pre-allocate
mu=zeros(6,n,'like',P0);
umdot=zeros(n_q,1,'like',P0);


%Forward recursion
for i=1:n
    
    if robot.joints(i).parent_link==0
        %First joint
        mu(1:6,i)=Bi0(1:6,1:6,i)*(P0*u0dot);
    else
        %Rest of the links
        if robot.joints(robot.joints(i).parent_link).type~=0
            mu_aux=(pm(1:6,robot.joints(robot.joints(i).parent_link).id)*umdot(robot.joints(i-1).q_id)+mu(1:6,robot.joints(robot.joints(i).parent_link).id));
        else
            mu_aux=mu(1:6,robot.joints(robot.joints(i).parent_link).id);
        end
        mu(1:6,i)=Bij(1:6,1:6,i,i-1)*mu_aux;
    end
    
    %Initialize
    if robot.joints(i).type~=0
        umdot(robot.joints(i).q_id,1)=phi_tilde(robot.joints(i).q_id)-psi(1:6,i)'*mu(1:6,i);
    end
end

end