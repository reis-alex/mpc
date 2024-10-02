function [t0dot,tLdot]=Accelerations(t0,tL,P0,pm,Bi0,Bij,u0,um,u0dot,umdot,robot)
% Computes the operational-space accelerations (twist-rate) of the multibody system.
%
% [t0dot,tLdot]=Accelerations(t0,tL,P0,pm,Bi0,Bij,u0,um,u0dot,umdot,robot)
%
% :parameters: 
%   * t0 -- Base-link twist [\omega,rdot], projected in the inertial CCS -- as a [6x1] matrix.
%   * tL -- Manipulator twist [\omega,rdot], projected in the inertial CCS -- as a [6xn] matrix.
%   * Bij -- Twist-propagation matrix (for manipulator i>0 and j>0) -- as a [6x6xn] matrix.
%   * Bi0 -- Twist-propagation matrix (for i>0 and j=0) -- as a [6x6xn] matrix.
%   * P0 -- Base-link twist-propagation "vector" -- as a [6x6] matrix.
%   * pm -- Manipulator twist-propagation "vector" -- as a [6xn] matrix.
%   * u0 -- Base-link velocities [\omega,rdot]. The angular velocity is projected in the body-fixed CCS, while the linear velocity is projected in the inertial CCS -- [6x1].
%   * um -- Joint velocities -- [n_qx1].
%   * u0dot -- Base-link accelerations [\omegadot,rddot]. The angular acceleration is projected in a body-fixed CCS, while the linear acceleration is projected in the inertial CCS -- [6x1].
%   * umdot -- Manipulator joint accelerations -- [n_qx1].
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return: 
%   * t0dot -- Base-link twist-rate vector \omegadot,rddot], projected in inertial frame -- as a [6x1] matrix.
%   * tLdot -- Manipulator twist-rate vector \omegadot,rddot], projected in inertial frame -- as a [6xn] matrix.
%
% See also: :func:`src.kinematics_dynamics.Jacobdot`. 

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

%--- Omega matrices ---%
%Base-link
Omega0=[SkewSym(t0(1:3)), zeros(3,3);
    zeros(3,3), zeros(3,3)];

%Pre-allocate
Omegam=zeros(6,6,n,'like',t0);

%Compute Omega for manipulator
for i=1:n
    Omegam(1:6,1:6,i)=[SkewSym(tL(1:3,i)), zeros(3,3);
        zeros(3,3), SkewSym(tL(1:3,i))];
end

%--- Twist Rate ---%
%Base-link
t0dot = Omega0*P0*u0+P0*u0dot;

%Pre-allocate
tLdot=zeros(6,n,'like',t0);

%Forward recursion
for i=1:n
    
    if robot.joints(i).parent_link==0
        %First Link
        tLdot(1:6,i)=Bi0(1:6,1:6,i)*t0dot+[zeros(3,6);SkewSym(t0(4:6)-tL(4:6,i)),zeros(3,3)]*t0;
    else
        %Rest of the links
        tLdot(1:6,i)=Bij(1:6,1:6,i,robot.joints(i).parent_link)*tLdot(1:6,robot.joints(i).parent_link)+[zeros(3,6); SkewSym(tL(4:6,robot.joints(i).parent_link)-tL(4:6,i)), zeros(3,3)]*tL(1:6,robot.joints(i).parent_link);
    end
    
    %Add joint contribution
    if robot.joints(i).type~=0
        tLdot(1:6,i)=tLdot(1:6,i)+Omegam(1:6,1:6,i)*pm(1:6,i)*um(robot.joints(i).q_id)+pm(1:6,i)*umdot(robot.joints(i).q_id);
    end
    
    
end


end