function [C0, C0m, Cm0, Cm] = CIM(t0,tL,I0,Im,M0_tilde,Mm_tilde,Bij,Bi0,P0,pm,robot)

import casadi.*
% Computes the Generalized Convective Inertia Matrix C of the multibody system.
%
% :parameters: 
%   * t0 -- Base-link twist [\omega,rdot], projected in the inertial CCS -- as a [6x1] matrix.
%   * tL -- Manipulator twist [\omega,rdot], projected in the inertial CCS -- as a [6xn] matrix.
%   * I0 -- Base-link inertia matrix, projected in the inertial CCS -- as a [3x3] matrix.
%   * Im -- Links inertia matrices, projected in the inertial CCS -- as a [3x3xn] matrix.
%   * M0_tilde -- Base-link mass composite body matrix -- as a [6x6] matrix .
%   * Mm_tilde -- Manipulator mass composite body matrix -- as a [6x6xn] matrix.
%   * Bij -- Twist-propagation matrix (for manipulator i>0 and j>0) -- as a [6x6xn] matrix.
%   * Bi0 -- Twist-propagation matrix (for i>0 and j=0) -- as a [6x6xn] matrix.
%   * P0 -- Base-link twist-propagation "vector" -- as a [6x6] matrix.
%   * pm -- Manipulator twist-propagation "vector" -- as a [6xn] matrix.
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return: 
%   * C0 -> Base-link convective inertia matrix -- as a [6x6] matrix.
%   * C0m -> Base-link - manipulator coupling convective inertia matrix -- as a [6xn_q] matrix.
%   * Cm0 -> Manipulator - base-link coupling convective inertia matrix -- as a [n_qx6] matrix.
%   * Cm -> Manipulator convective inertia matrix -- as a [n_qxn_q] matrix.
%
% To obtain the full convective inertia matrix C:
%
% .. code-block:: matlab
%   
%   %Compute the Convective Inertia Matrix C
%   [C0, C0m, Cm0, Cm] = CIM(t0,tL,I0,Im,M0_tilde,Mm_tilde,Bij,Bi0,P0,pm,robot)
%   C=[C0,C0m;Cm0,Cm];
%
% See also: :func:`src.kinematics_dynamics.GIM`.

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
n=robot.n_links_joints;

%--- Omega ---%
%Base-link Omega
Omega0=[SkewSym(t0(1:3)), zeros(3,3);
    zeros(3,3), zeros(3,3)];

%Pre-allocate Omega
Omega=SX.sym('Omega_',6,6,n);

%Compute Omega
for i=1:n
    Omega{i}=[SkewSym(tL(1:3,i)), zeros(3,3);
        zeros(3,3), SkewSym(tL(1:3,i))];
end

%--- Mdot ---%
%Base-link Mdot
Mdot0=[Omega0(1:3,1:3)*I0, zeros(3,3); zeros(3,3), zeros(3,3)];

%Pre-allocate
Mdot=SX.sym('M_dot',6,6,n);

%Compute Mdot
for i=1:n
    Mdot{i}=[Omega{i}*Im{i}, zeros(3,3); zeros(3,3), zeros(3,3)];
end

%--- Mdot tilde ---%
%Pre-Allocate
Mdot_tilde=SX.sym('M_dot_ti',6,6,n);

%Backwards recursion
for i=n:-1:1
    %Initialize
    Mdot_tilde{i}=Mdot{i};
    %Add children contributions
    child=find(robot.con.child(:,i))';
    for j=1:length(child)
        Mdot_tilde{i}=Mdot_tilde{i}+Mdot_tilde{child(j)};
    end
end
%Base-link
Mdot0_tilde = Mdot0;
%Add children contributions
child=find(robot.con.child_base)';
for j=1:length(child)
    Mdot0_tilde = Mdot0_tilde+Mdot_tilde{child(j)};
end


%--- Bdot ---%

%Pre-allocate Bdotij
Bdotij=zeros(6,6,n,n,'like',t0);

%Compute Bdotij
for j=1:n
    for i=1:n
        if robot.con.branch(i,j)==1
            %Links are in the same branch
            Bdotij(1:6,1:6,i,j)=[zeros(3,3), zeros(3,3); SkewSym(tL(4:6,j)-tL(4:6,i)), zeros(3,3)];
        else
            %Links are not in the same branch
            Bdotij(1:6,1:6,i,j)=zeros(6,6);
        end
    end
end


%--- Hij tilde ---%
%Pre-allocate Hij_tilde
Hij_tilde=zeros(6,6,n,n,'like',t0);

%Hij_tilde
for i=n:-1:1
    for j=n:-1:1
        Hij_tilde(1:6,1:6,i,j)=Mm_tilde(1:6,1:6,i)*Bdotij(1:6,1:6,i,j);
        %Add children contributions
        child=find(robot.con.child(:,i))';
        for k=1:length(child)
            Hij_tilde(1:6,1:6,i,j)=Hij_tilde(1:6,1:6,i,j)+Bij(1:6,1:6,child(k),i)'*Hij_tilde(1:6,1:6,child(k),i);
        end
    end
end

%Pre-allocate Hi0_tilde and H0j_tilde
Hi0_tilde=zeros(6,6,n,'like',t0);

%Hi0_tilde
for i=n:-1:1
    Bdot=[zeros(3,3), zeros(3,3); SkewSym(t0(4:6)-tL(4:6,i)), zeros(3,3)];
    Hi0_tilde(1:6,1:6,i)=Mm_tilde(1:6,1:6,i)*Bdot;
    %Add children contributions
    child=find(robot.con.child(:,i))';
    for k=1:length(child)
        Hi0_tilde(1:6,1:6,i)=Hi0_tilde(1:6,1:6,i)+Bij(1:6,1:6,child(k),i)'*Hij_tilde(1:6,1:6,child(k),i);
    end
end

%--- C Matrix ---%
%Pre-allocate
Cm=zeros(n_q,n_q,'like',t0);
C0m=zeros(6,n_q,'like',t0);
Cm0=zeros(n_q,6,'like',t0);

%Cm Matrix
for j=1:n
    for i=1:n
        %Joints must not be fixed and links on the same branch
        if (robot.joints(i).type~=0 && robot.joints(j).type~=0) && (robot.con.branch(i,j)==1 || robot.con.branch(j,i)==1)
            %Compute Cm matrix
            if i<=j
                %Add children contributions
                child_con=zeros(6,6);
                child=find(robot.con.child(:,j))';
                for k=1:length(child)
                    child_con=child_con+Bij(1:6,1:6,child(k),i)'*Hij_tilde(1:6,1:6,child(k),j);
                end
                Cm(robot.joints(i).q_id,robot.joints(j).q_id)=pm(1:6,i)'*(Bij(1:6,1:6,j,i)'*Mm_tilde(1:6,1:6,j)*Omega(1:6,1:6,j)+child_con+Mdot_tilde(1:6,1:6,j))*pm(1:6,j);
            else
                Cm(robot.joints(i).q_id,robot.joints(j).q_id)=pm(1:6,i)'*(Mm_tilde(1:6,1:6,i)*Bij(1:6,1:6,i,j)*Omega(1:6,1:6,j)+Hij_tilde(1:6,1:6,i,j)+Mdot_tilde(1:6,1:6,i))*pm(1:6,j);
            end
        end
    end
end
%C0 matrix
%Add children contributions
child_con=zeros(6,6);
child=find(robot.con.child_base)';
for k=1:length(child)
    child_con=child_con+Bi0(1:6,1:6,child(k))'*Hi0_tilde(1:6,1:6,child(k));
end
C0 = P0'*(M0_tilde*Omega0+child_con+Mdot0_tilde)*P0;
%C0m
for j=1:n
    if  robot.joints(j).type~=0
        if j==n
            C0m(1:6,robot.joints(j).q_id)=P0'*(Bi0(1:6,1:6,j)'*Mm_tilde(1:6,1:6,j)*Omega(1:6,1:6,j)+Mdot_tilde(1:6,1:6,j))*pm(1:6,j);
        else
            %Add children contributions
            child_con=zeros(6,6);
            child=find(robot.con.child(:,j))';
            for k=1:length(child)
                child_con=child_con+Bi0(1:6,1:6,child(k))'*Hij_tilde(1:6,1:6,child(k),j);
            end
            C0m(1:6,robot.joints(j).q_id)=P0'*(Bi0(1:6,1:6,j)'*Mm_tilde(1:6,1:6,j)*Omega(1:6,1:6,j)+child_con+Mdot_tilde(1:6,1:6,j))*pm(1:6,j);
        end
    end
end
%Cm0
for i=1:n
    if  robot.joints(i).type~=0  
        Cm0(robot.joints(i).q_id,1:6)=pm(1:6,i)'*(Mm_tilde(1:6,1:6,i)*Bi0(1:6,1:6,i)*Omega0+Hi0_tilde(1:6,1:6,i)+Mdot_tilde(1:6,1:6,i))*P0;
    end
end

end

