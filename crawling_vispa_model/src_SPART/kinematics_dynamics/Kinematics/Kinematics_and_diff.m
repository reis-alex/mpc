function [dTJ0,dTL0,dTJ,dTL,RJ,RL,rJ,rL,e,g]=Kinematics_and_diff(R0,r0,q0,q,robot)
% Computes the kinematics -- positions and orientations -- of the multibody
% system and geometric partial derivatives
%
% [RJ,RL,rJ,rL,e,g]=Kinematics(R0,r0,q,robot)
%
% :parameters:
%   * R0 -- Rotation matrix from the base-link CCS to the inertial CCS -- [3x3].
%   * r0 -- Position of the base-link center-of-mass with respect to the origin of the inertial frame, projected in the inertial CCS -- [3x1].
%   * q -- Displacements of the active joints -- [n_qx1].
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return:
%   * dTJ0 -- joint CCS 4x4 homogenous transformation matrices partial
%     derivatives over q0 with respect to inertial CCS -- as a [4x4xnx6]
%   * dTL0 -- Links CCS 4x4 homogenous transformation matrices partial
%     derivatives over q0 with respect to inertial CCS -- as a [4x4xnx6]
%   * dTJ -- Joints CCS 4x4 homogenous transformation matrices partial
%     derivatives over q with respect to inertial CCS -- as a [4x4xnxn_q]
%   * dTL -- Links CCS 4x4 homogenous transformation matrices partial
%     derivatives over q with respect to inertial CCS -- as a [4x4xnxn_q]
%   * RJ -- Joints CCS 3x3 rotation matrices with respect to the inertial CCS  -- as a [3x3xn] matrix.
%   * RL -- Links CCS 3x3 rotation matrices with respect to the inertial CCS -- as a [3x3xn] matrix.
%   * rJ -- Positions of the joints, projected in the inertial CCS -- as a [3xn] matrix.
%   * rL -- Positions of the links, projected in the inertial CCS -- as a [3xn] matrix.
%   * e -- Joint rotation/sliding axes, projected in the inertial CCS -- as a [3xn] matrix.
%   * g -- Vector from the origin of the ith joint CCS to the origin of the ith link CCS, projected in the inertial CCS -- as a [3xn] matrix.
%
% Remember that all the ouput magnitudes are projected in the **inertial frame**.
%
% Examples on how to retrieve the results from a specific link/joint:
%
%   To retrieve the position of the ith link: ``rL(1:3,i)``.
%
%   To retrieve the rotation matrix of the ith joint: ``RJ(1:3,1:3,i)``.
%   0
%
% See also: :func:`src.robot_model.urdf2robot` and :func:`src.robot_model.DH_Serial2robot`.

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

%--- Number of links and joints ---%
n=robot.n_links_joints;
n_q=robot.n_q;

%--- Homogeneous transformation matrices ---%

%Pre-allocate homogeneous transformations matrices
TJ=zeros(4,4,n,'like',R0);
TL=zeros(4,4,n,'like',R0);
dTJ0=zeros(4,4,n,6,'like',R0);
dTL0=zeros(4,4,n,6,'like',R0);
dTJ=zeros(4,4,n,n,'like',R0);
dTL=zeros(4,4,n,n,'like',R0);
%--- Base-link ---%
T0=[R0,r0;zeros(1,3),1];

%--- Forward kinematics recursion ---%

%Obtain the joints and links kinematics
for i=1:n

    %Get child joint
    cjoint=robot.joints(i);

    %Joint kinematics (homogeneous transformation matrix)
    if cjoint.parent_link==0
        %Parent link is the base-link
        TJ(1:4,1:4,cjoint.id)=T0*cjoint.T;
        dT0_dq0 = dT0_diffq0(q0,R0);
        for coord_i = 1:6 %derivation over base coordinates

            dT0_coordi = dT0_dq0(:,:,coord_i);%[dR0_d0i dr0_d0i; zeros(1,4)];
            dTJ0(1:4,1:4,cjoint.id,coord_i)=dT0_coordi*cjoint.T;
        end
    else
        %Parent link is not the base-link
        TJ(1:4,1:4,cjoint.id)=TL(1:4,1:4,cjoint.parent_link)*cjoint.T;
        for coord_i = 1:6 %derivation over base coordinates
            dTJ0(1:4,1:4,cjoint.id,coord_i)=dTL0(1:4,1:4,cjoint.parent_link,coord_i)*cjoint.T;
        end
    end

    %Transformation due to current joint variable
    if cjoint.type==1
        %Revolute
        T_qm=[Euler_DCM(cjoint.axis,q(cjoint.q_id))',zeros(3,1);zeros(1,3),1];
        dT_qm=[Euler_DCM_dq(cjoint.axis,q(cjoint.q_id))',zeros(3,1);zeros(1,4)];
    elseif cjoint.type==2
        %Prismatic
        T_qm=[eye(3),cjoint.axis*q(cjoint.q_id);zeros(1,3),1];
        dT_qm=[zeros(3,3) cjoint.axis;zeros(1,4)];
    else
        %Fixed
        T_qm=[eye(3),zeros(3,1);zeros(1,3),1];
        dT_qm=zeros(4,4);
    end


    %Link Kinematics (homogeneous transformation matrix)
    clink=robot.links(cjoint.child_link);
    TL(1:4,1:4,clink.id)=TJ(1:4,1:4,clink.parent_joint)*T_qm*clink.T;

    for coord_i = 1:6
        dTL0(1:4,1:4,clink.id,coord_i)=dTJ0(1:4,1:4,clink.parent_joint,coord_i)*T_qm*clink.T;
    end
    
    if cjoint.type ~= 0
        dTL(1:4,1:4,clink.id,cjoint.q_id)=TJ(1:4,1:4,clink.parent_joint)*dT_qm*clink.T; 
    end
    %Backward kinematics recursion for partial derivative relative to joint
    %variable

    parent_chain = find(robot.con.branch(cjoint.id,:));

    for ind_chain=length(parent_chain)-1:-1:1
        ind_body=parent_chain(ind_chain);
        if(robot.joints(ind_body).type~=0)
            ind_q=robot.joints(ind_body).q_id;
            if cjoint.parent_link ~= 0
                dTJ(:,:,cjoint.id,ind_q)=dTL(1:4,1:4,cjoint.parent_link,ind_q)*cjoint.T;
            else
                dTJ(:,:,cjoint.id,ind_q)=zeros(4,4, "like", R0);
            end
        end
    end

    %Backward kinematics recursion for partial derivative relative to joint
    %variable
    for ind_chain=length(parent_chain)-1:-1:1
        ind_body=parent_chain(ind_chain);
        if(robot.joints(ind_body).type~=0)
            ind_q=robot.joints(ind_body).q_id;
            dTL(1:4,1:4,clink.id,ind_q)=dTJ(1:4,1:4,clink.parent_joint,ind_q)*T_qm*clink.T;
        end
    end



end

%--- Rotation matrices, translation, position and other geometric quantities ---%

%Pre-allocate rotation matrices, translation and positions
RJ=zeros(3,3,n,'like',R0);
RL=zeros(3,3,n,'like',R0);
rJ=zeros(3,n,'like',R0);
rL=zeros(3,n,'like',R0);
%Pre-allocate rotation/sliding axis
e=zeros(3,n,'like',R0);
%Pre-allocate other geometric quantities
g=zeros(3,n,'like',R0);

%Format rotation matrices, link positions, joint axis and other geometric
%quantities

%Joint associated quantities
for i=1:n
    RJ(1:3,1:3,i)=TJ(1:3,1:3,i);
    rJ(1:3,i)=TJ(1:3,4,i);
    e(1:3,i)=RJ(1:3,1:3,i)*robot.joints(i).axis;
end
%Link associated quantities
for i=1:n
    RL(1:3,1:3,i)=TL(1:3,1:3,i);
    rL(1:3,i)=TL(1:3,4,i);
    g(1:3,i)=rL(1:3,i)-rJ(1:3,robot.links(i).parent_joint);
end


end