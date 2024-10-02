function [J0dot, Jmdot]=Jacobdot(rp,tp,r0,t0,rL,tL,P0,pm,i,robot)
% Computes the geometric Jacobian time-derivative of a point `p`.
%
% [J0dot, Jmdot]=Jacobdot(rp,tp,r0,t0,rL,tL,P0,pm,i,robot)
%
% :parameters: 
%   * rp -- Position of the point of interest, projected in the inertial CCS -- [3x1].
%   * tp -- Twist of the point of interest [\omega,rdot], projected in the intertial CCS -- [6x1].
%   * r0 -- Position of the base-link center-of-mass with respect to the origin of the inertial frame, projected in the inertial CCS -- [3x1].
%   * t0 -- Base-link twist [\omega,rdot], projected in the inertial CCS -- as a [6x1] matrix.
%   * rL -- Positions of the links, projected in the inertial CCS -- as a [3xn] matrix.
%   * tL -- Manipulator twist [\omega,rdot], projected in the inertial CCS -- as a [6xn] matrix.
%   * P0 -- Base-link twist-propagation "vector" -- as a [6x6] matrix.
%   * pm -- Manipulator twist-propagation "vector" -- as a [6xn] matrix.
%   * i -- Link id where the point `p` is located -- int 0 to n.
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return: 
%   * J0dot -- Base-link Jacobian time-derivative -- as a [6x6] matrix.
%   * Jmdot -- Manipulator Jacobian time-derivative -- as a [6xn_q] matrix.
%
% Examples:
%
%   To compute the acceleration of a point ``p`` on the ith link:
%
% .. code-block:: matlab
%   
%   %Compute Jacobians
%   [J0, Jm]=Jacob(rp,r0,rL,P0,pm,i,robot);
%   Compute Jacobians time-derivatives
%   [J0dot, Jmdot]=Jacobdot(rp,tp,r0,t0,rL,tL,P0,pm,i,robot)
%   %Twist-rate of that point
%   tpdot=J0*u0dot+J0dot*u0+Jm*umdot+Jmdot*um;
%
% See also: :func:`src.kinematics_dynamics.Accelerations` and :func:`src.kinematics_dynamics.Jacob`. 

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
import casadi.*
%=== CODE ===%

%--- Omega ---%
%Base-link Omega
Omega0=[SkewSym(t0(1:3)), zeros(3,3);
    zeros(3,3), zeros(3,3)];

%Pre-allocate Omega
Omega=SX.sym('Omega',6,6,robot.n_links_joints);

%Compute Omega
for j=1:i
    Omega{j}=[SkewSym(tL(1:3,j)), zeros(3,3);
        zeros(3,3), SkewSym(tL(1:3,j))];
end


%--- Jacobian time-derivative ---%

%Base-link Jacobian
J0dot=[eye(3),zeros(3,3);SkewSym(r0-rp),eye(3)]*Omega0*P0+[zeros(3,3),zeros(3,3);SkewSym(t0(4:6)-tp(4:6)),zeros(3,3)]*P0;

%Pre-allocate
Jmdot=zeros(6,robot.n_q,'like',rp);
    
%Manipulator Jacobian
joints_num=0;
for j=1:i
    %If joint is not fixed
    if robot.joints(j).type~=0
        if robot.con.branch(i,j)==1
            Jmdot(1:6,robot.joints(j).q_id)=[eye(3),zeros(3,3);SkewSym(rL(1:3,j)-rp),eye(3)]*Omega{j}*pm(1:6,j)+[zeros(3,3),zeros(3,3);SkewSym(tL(4:6,j)-tp(4:6)),zeros(3,3)]*pm(1:6,j);
        else
            Jmdot(1:6,robot.joints(j).q_id)=zeros(6,1);
        end
        joints_num=joints_num+1;
    end
end

end