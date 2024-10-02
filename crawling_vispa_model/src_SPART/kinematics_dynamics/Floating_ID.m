function [taum,u0dot] = Floating_ID(wF0,wFm,Mm_tilde,H0,t0,tm,P0,pm,I0,Im,Bij,Bi0,u0,um,umdot,robot)
% This function solves the inverse dynamics problem (it obtains the
% generalized forces from the accelerations) for a manipulator with 
% a floating base.
%
% [taum,u0dot] = Floating_ID(wF0,wFm,Mm_tilde,H0,t0,tm,P0,pm,I0,Im,Bij,Bi0,u0,um,umdot,robot)
%
% :parameters: 
%   * wF0 -- Wrench acting on the base-link center-of-mass [n,f], projected in the inertial CCS -- as a [6x1] matrix.
%   * wFm -- Wrench acting on the links center-of-mass  [n,f], projected in the inertial CCS -- as a [6xn] matrix.
%   * M0_tilde -- Base-link mass composite body matrix -- as a [6x6] matrix .
%   * Mm_tilde -- Manipulator mass composite body matrix -- as a [6x6xn] matrix.
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
%   * umdot -- Manipulator joint accelerations -- [n_qx1].
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return: 
%   * tau0 -- Base-link forces [n,f]. The torque n is projected in the body-fixed CCS, while the force f is projected in the inertial CCS -- [6x1].
%   * taum -- Joint forces/torques -- as a [n_qx1] matrix.
%
% See also: :func:`src.kinematics_dynamics.sID` and :func:`src.kinematics_dynamics.FD`. 

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

%Number of links and Joints
n=robot.n_links_joints;

%Recompute Accelerations with qddot=0
[t0dot,tmdot]=Accelerations(t0,tm,P0,pm,Bi0,Bij,u0,um,zeros(6,1),umdot,robot);

%Use the inverse dynamics
[tau0_0ddot,taum] = ID(wF0,wFm,t0,tm,t0dot,tmdot,P0,pm,I0,Im,Bij,Bi0,robot);

%Pre-allocate kappa
kappa=zeros(n,6,'like',wF0);

%Compute Kappa
for i=1:n
    kappa(i,1:6)=pm(1:6,i)'*Mm_tilde(1:6,1:6,i)*Bi0(1:6,1:6,i);
end

%Compute base-link acceleration
u0dot=-H0\tau0_0ddot;

%Update joint forces
for i=1:n
    if robot.joints(i).type~=0
        taum(robot.joints(i).q_id)=kappa(i,1:6)*P0*u0dot+taum(robot.joints(i).q_id);
    end
end

end
