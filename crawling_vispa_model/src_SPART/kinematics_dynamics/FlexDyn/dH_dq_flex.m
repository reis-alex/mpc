function [dHf0_dq0, dHf0_dq, dHfq_dq0, dHfq_dq, dHf_dq0, dHf_dq] = dH_dq_flex(r0,R0,P0,RJ,RL,rJ,rL,k,g, dTJ, dTL,dTJ0,dTL0, Bij,Bi0,pm,q,q0,Jf0dot, Jfqdot,robot)
% Computes the flex partial derivatives of the Generalized Inertia Matrix (GIM) H of the multibody vehicle.
%

%--- Number of links and Joints ---%
n_q=robot.n_q;

%--- Number of links  ---%
n_f=robot.n_links_flex;

%--- Number of DOF  ---%
n_m=robot.flexs.nb_mode;

%init of gradients

dHf0_dq0 = zeros(n_m,6,6,'like',Jf0dot);
dHf0_dq  = zeros(n_m,6,n_q,'like',Jf0dot);

dHfq_dq0 = zeros(n_m,n_q,6,'like',Jfqdot);
dHfq_dq  = zeros(n_m,n_q,n_q,'like',Jfqdot);

dHf_dq0  = zeros(n_m,n_m,6,'like',Jf0dot);
dHf_dq   = zeros(n_m,n_m,n_q,'like',Jfqdot);


%initof partial derivations

dT0_dq0 = dT0_diffq0(q0, R0);
dBi0_q0 = dBi0_dr0();
dBi0_q  = dBi0_ql(r0,R0,RJ,RL,rJ,rL,dTL,dTJ,q,robot);

% loop on the flexible links
for ind_q0 = 1:6
    for ind_l=1:n_f
        %Indice of the flexible DOFs
        if ind_l==n_f
            ind_DOF=robot.flexs.link_mode_id(ind_l):robot.flexs.nb_mode;
        else
            ind_DOF=robot.flexs.link_mode_id(ind_l):robot.flexs.link_mode_id(ind_l+1)-1;
        end
        %Indice of the flexible links
        ind_link=robot.flexs.link_id(ind_l);
        
        %rotation matrix from Inertia CCS to the link CCS
        R_CCS_link=[RL(:,:,ind_link)' zeros(3,3,'like',RL); zeros(3,3,'like',RL) RL(:,:,ind_link)'];
        
        %partial derivative of jacobian of indl element
        [Jf0i_dq0, Jfmi_dq0] = Jacob_dq0(r0,R0,P0,q,RJ,RL,rJ,rL, k,g, dTJ0, dTL0, pm, Bij, Bi0, ind_link, ind_q0, dT0_dq0, dBi0_q0, q0, robot );
        dHf0_dq0(ind_DOF,:,ind_q0)=robot.flexs.Li(ind_DOF,[4:6 1:3])*Jf0i_dq0;
        dHfq_dq0(ind_DOF,:,ind_q0)=robot.flexs.Li(ind_DOF,[4:6 1:3])*Jfmi_dq0;

    end
end

for ind_q = 1:n_q
    for ind_l=1:n_f
        %Indice of the flexible DOFs
        if ind_l==n_f
            ind_DOF=robot.flexs.link_mode_id(ind_l):robot.flexs.nb_mode;
        else
            ind_DOF=robot.flexs.link_mode_id(ind_l):robot.flexs.link_mode_id(ind_l+1)-1;
        end
        %Indice of the flexible links
        ind_link  = robot.flexs.link_id(ind_l);
        ind_joint = robot.flexs.joint_id(ind_l);
        %rotation matrix from Inertia CCS to the link CCS
        R_CCS_link=[RL(:,:,ind_link)' zeros(3,3,'like',RL); zeros(3,3,'like',RL) RL(:,:,ind_link)'];
        
        [J0j_dq, Jmj_dq] = Jacob_dq(r0,R0,P0,RJ,RL,rJ,rL, k, g, dTJ, dTL,  pm, ind_link, robot.joints(ind_q).id, dBi0_q, robot);
        dHf0_dq(ind_DOF,:,ind_q)=robot.flexs.Li(ind_DOF,[4:6 1:3])*J0j_dq;
        dHfq_dq(ind_DOF,:,ind_q)=robot.flexs.Li(ind_DOF,[4:6 1:3])*Jmj_dq;
        
    end
end

end