function [tau0,taum] = ID(wF0,wFm,t0,tL,t0dot,tLdot,P0,pm,I0,Im,Bij,Bi0,robot)
% This function solves the inverse dynamics (ID) problem (it obtains the
% generalized forces from the accelerations) for a manipulator.
%
% [tau0,taum] = ID(wF0,wFm,t0,tL,t0dot,tLdot,P0,pm,I0,Im,Bij,Bi0,robot)
% 
% :parameters: 
%   * wF0 -- Wrench acting on the base-link center-of-mass [n,f], projected in the inertial CCS -- as a [6x1] matrix.
%   * wFm -- Wrench acting on the links center-of-mass  [n,f], projected in the inertial CCS -- as a [6xn] matrix.
%   * t0 -- Base-link twist [\omega,rdot], projected in the inertial CCS -- as a [6x1] matrix.
%   * tL -- Manipulator twist [\omega,rdot], projected in the inertial CCS -- as a [6xn] matrix.
%   * t0dot -- Base-link twist-rate vector \omegadot,rddot], projected in inertial frame -- as a [6x1] matrix.
%   * tLdot -- Manipulator twist-rate vector \omegadot,rddot], projected in inertial frame -- as a [6xn] matrix.
%   * P0 -- Base-link twist-propagation "vector" -- as a [6x6] matrix.
%   * pm -- Manipulator twist-propagation "vector" -- as a [6xn] matrix.
%   * I0 -- Base-link inertia matrix, projected in the inertial CCS -- as a [3x3] matrix.
%   * Im -- Links inertia matrices, projected in the inertial CCS -- as a [3x3xn] matrix.
%   * Bij -- Twist-propagation matrix (for manipulator i>0 and j>0) -- as a [6x6xn] matrix.
%   * Bi0 -- Twist-propagation matrix (for i>0 and j=0) -- as a [6x6xn] matrix.
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return: 
%   * tau0 -- Base-link forces [n,f]. The torque n is projected in the body-fixed CCS, while the force f is projected in the inertial CCS -- [6x1].
%   * taum -- Joint forces/torques -- as a [n_qx1] matrix.
%
% See also: :func:`src.kinematics_dynamics.Floating_ID` and :func:`src.kinematics_dynamics.FD`. 


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

%--- Mdot ---%
%base-link Mdot
Mdot0=[SkewSym(t0(1:3))*I0, zeros(3,3); zeros(3,3), zeros(3,3)];

%Pre-allocate
Mdot=zeros(6,6,n,'like',wF0);

%Manipulator Mdot
for i=1:n
    Mdot(1:6,1:6,i)=[SkewSym(tL(1:3,i))*Im(1:3,1:3,i), zeros(3,3); zeros(3,3), zeros(3,3)];
end

%--- Forces ---%
%Base-link
wq0=[I0,zeros(3,3);zeros(3,3),robot.base_link.mass*eye(3)]*t0dot+Mdot0*t0-wF0;

%Pre-allocate
wq=zeros(6,n,'like',wF0);

%Manipulator
for i=1:n
    wq(1:6,i)=[Im(1:3,1:3,i),zeros(3,3);zeros(3,3),robot.links(i).mass*eye(3)]*tLdot(1:6,i)+Mdot(1:6,1:6,i)*tL(1:6,i)-wFm(1:6,i);
end

%Pre-allocate
wq_tilde=zeros(6,n,'like',wF0);

%Backwards recursion
for i=n:-1:1
    %Initialize wq_tilde
    wq_tilde(1:6,i)=wq(1:6,i);
    %Add children contributions
    for j=find(robot.con.child(:,i))'
        wq_tilde(1:6,i)=wq_tilde(1:6,i)+Bij(1:6,1:6,j,i)'*wq_tilde(1:6,j);
    end
end
%Base-link
wq_tilde0=wq0;
%Add children contributions
for j=find(robot.con.child_base)'
    wq_tilde0=wq_tilde0+Bi0(1:6,1:6,j)'*wq_tilde(1:6,j);
end

%---- Joint forces ---%
%Base-link
tau0=P0'*wq_tilde0;

%Pre-allocate
taum=zeros(robot.n_q,1,'like',wF0);

%Manipulator joint forces.
for i=1:n
    if robot.joints(i).type~=0
        taum(robot.joints(i).q_id,1)=pm(1:6,i)'*wq_tilde(1:6,i);
    end
end

end