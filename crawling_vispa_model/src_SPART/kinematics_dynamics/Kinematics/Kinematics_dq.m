function [dTJ,dTL] = Kinematics_dq(R0,r0,RJ,RL,rJ,rL,e,g,q0,q,robot)
%KINEMATICS_DQ Summary of this function goes here
%   Detailed explanation goes here
n=robot.n_links_joints;
n_q=robot.n_q;

dTJ_temp=zeros(4,4,n,n_q,'like',R0);
dTL_temp=zeros(4,4,n,n_q,'like',R0);

dTJ=zeros(4,4,n,n_q,'like',R0);
dTL=zeros(4,4,n,n_q,'like',R0);

T0=[R0,r0;zeros(1,3),1];

for ind_l = 1:n_q
    for i=1:n

        %Get child joint
        cjoint=robot.joints(i);

        %Joint kinematics (homogeneous transformation matrix)
        if cjoint.parent_link==0
            %Parent link is the base-link
            dTJ_temp(1:4,1:4,cjoint.id,ind_l)=T0*cjoint.T;

        else
            %Parent link is not the base-link
            dTJ_temp(1:4,1:4,cjoint.id,ind_l)=dTL_temp(1:4,1:4,cjoint.parent_link,ind_l)*cjoint.T;

        end

        %%Transformation due to current joint variable
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
        if cjoint.q_id == ind_l
            dTL_temp(1:4,1:4,clink.id,ind_l)=dTJ_temp(1:4,1:4,clink.parent_joint,ind_l)*dT_qm*clink.T;
        else
            dTL_temp(1:4,1:4,clink.id,ind_l)=dTJ_temp(1:4,1:4,clink.parent_joint,ind_l)*T_qm*clink.T;  
        end
    end

    for i=1:n
        cjoint=robot.joints(i);
        clink=robot.links(cjoint.child_link);
        if  robot.con.branch(cjoint.id, ind_l)
            if cjoint.id ~= ind_l
                dTJ(1:4,1:4,cjoint.id,ind_l) = dTJ_temp(1:4,1:4,cjoint.id,ind_l);
            end
            dTL(1:4,1:4,clink.id,ind_l) = dTL_temp(1:4,1:4,clink.id,ind_l);
        else
            disp('pas dans chaine')
            disp([cjoint.q_id,ind_l])
            dTJ(1:4,1:4,cjoint.id,ind_l) = zeros(4,4,'like',R0);
            dTL(1:4,1:4,clink.id,ind_l) = zeros(4,4,'like',R0);
        end
    end
end
% dTJ_temp
% dTL_temp
% 
% for ind_l = 1:n_q
% 
% end
end