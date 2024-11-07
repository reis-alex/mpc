function [grad_H0_q0, grad_H0q_q0, grad_Hq_q0, grad_H0_q, grad_H0q_q, grad_Hq_q] = dH_dq_gradients(I0,Im,R0,r0,q0_dot,RJ,RL,rJ,rL,k,g, dTJ, dTL,dTJ0,dTL0, Bij,Bi0,P0,pm,q,q0,robot)
%DH_DQ_GRADIENTS Summary of this function goes here
%   Detailed explanation goes here
nb_base_DOF = 6;
n       = robot.n_links_joints;
n_q      = robot.n_q;


%% derivation over base-link variables (6 : 3 of XYZ rotation, 3 of translation)

grad_H0_q0 = zeros(6,6,6, 'like', R0); %contains dM0/dX0 for the 6 link_base variables
grad_H0q_q0 = zeros(6,n_q,6, 'like', R0); %contains dM0q/dX0 for the 6 link_base variables
grad_Hq_q0 = zeros(n_q,n_q,6, 'like', R0); %contains dMq/dX0 for the 6 link_base variables


dT0_dq0 = dT0_diffq0(q0, R0);
dR0_q0 = dT0_dq0(1:3,1:3,:);
%dR0_q0 = dR0_d0(R0);
[I0_dq0,Im_dq0]=dI_dq0(dTJ0,dTL0, R0, dR0_q0, RL,robot);

%dBi0_q0 = dBi0_dr0();

%[H0, H0m, Hm, ~, ~, ~] = M_lagrange(I0, Im, rL, Bi0, P0, r0, pm, robot);
for coord_i = 1:nb_base_DOF
    [dBij_dq0,dBi0_q0,dP0_dq0,dpm_dq0] = DiffKinematics_dq0(dTL0, dTJ0, R0, dT0_dq0, k,g,coord_i,robot);
    dP0_dq0 = [dT0_dq0(1:3,1:3,coord_i) zeros(3);
               zeros(3)            zeros(3)];

    grad_H0_q0_tilde = P0'*[I0_dq0(:,:,coord_i), zeros(3); zeros(3,6)]*P0...
                     + dP0_dq0'*[I0,zeros(3,3);zeros(3,3),robot.base_link.mass*eye(3)]*P0 ...
                     + P0'*[I0,zeros(3,3);zeros(3,3),robot.base_link.mass*eye(3)]*dP0_dq0;

    grad_H0_q0(:,:,coord_i) = grad_H0_q0_tilde;

    for ind_j=1:n
        dIm_d0 = Im_dq0(:,:,ind_j,coord_i);%diff(Im(:,:,ind_j), q0(ind_i));
        
        [J0j, Jmj] = Jacob(rL(1:3,ind_j),r0,rL,P0,pm,ind_j,robot);
        [J0j_dq0i, Jmj_dq0i] = Jacob_dq0(r0,R0,P0,q,RJ,RL,rJ,rL, k,g, dTJ0, dTL0, pm, Bij, Bi0, ind_j, coord_i, dT0_dq0, dBi0_q0(:,:,ind_j), q0, robot );
        
        grad_H0_q0(:,:,coord_i) = grad_H0_q0(:,:,coord_i) ...
                           + J0j_dq0i'*[Im(:,:,ind_j) zeros(3);
                                        zeros(3)  eye(3)*robot.links(ind_j).mass]*J0j ...
                           + J0j'*[dIm_d0 zeros(3);
                                   zeros(3)       zeros(3)]*J0j ...    
                           + J0j'*[Im(:,:,ind_j) zeros(3);
                                   zeros(3)      eye(3)*robot.links(ind_j).mass]*J0j_dq0i;

        grad_H0q_q0(:,:,coord_i) = grad_H0q_q0(:,:,coord_i) ...
                           + J0j_dq0i'*[Im(:,:,ind_j) zeros(3);
                                        zeros(3)      eye(3)*robot.links(ind_j).mass]*Jmj ...
                           + J0j'*[dIm_d0     zeros(3);
                                   zeros(3)   zeros(3)]*Jmj...    
                           + J0j'*[Im(:,:,ind_j) zeros(3);
                                  zeros(3)   eye(3)*robot.links(ind_j).mass]*Jmj_dq0i;

        grad_Hq_q0(:,:,coord_i) = grad_Hq_q0(:,:,coord_i) ...
                           + Jmj_dq0i'*[Im(:,:,ind_j) zeros(3);
                                        zeros(3)      eye(3)*robot.links(ind_j).mass]*Jmj ...
                           + Jmj'*[dIm_d0    zeros(3);
                                   zeros(3)  zeros(3)]*Jmj...    
                           + Jmj'*[Im(:,:,ind_j) zeros(3);
                                   zeros(3)  eye(3)*robot.links(ind_j).mass]*Jmj_dq0i;

    end
end
grad_H0_q = zeros(6,6,n_q, 'like', R0); %contains dM0/dq over the joint variables
grad_H0q_q= zeros(6,n_q,n_q, 'like', R0); %contains dM0q/dq for the joint variables
grad_Hq_q = zeros(n_q,n_q,n_q, 'like', R0); %contains dMq/dq for the joint variables

%grad_Hq_q_test = zeros(n_q,n_q,n_q, 'like', R0); %contains dMq/dq for the joint variables for test

dBi0_q = zeros(6,6,robot.n_links_joints, robot.n_q, 'like', R0);


dBi0_q  = dBi0_ql(r0,R0,RJ,RL,rJ,rL,dTL,dTJ,q,robot);
%dI_q    = dI_dq(r0,R0,RJ,RL,rJ,rL,q, robot);
% 
%[~, ~, Hm, ~, ~, ~] = M_lagrange(I0, Im, rL, Bi0, P0, r0, pm, robot);
diff_I_q = zeros(3,3,robot.n_links_joints, robot.n_q, "like", R0);

for ind_k = 1: robot.n_links_joints
for l = 1:robot.n_q

[RJ_ql,RL_ql,rJ_ql,rL_ql,e_ql,g_ql]  = dT_parser(dTJ, dTL, robot, l, R0);
        
        dRj_dql = RL_ql(:,:,ind_k);
        diff_I_q(:,:,ind_k,l) = dRj_dql*robot.links(ind_k).inertia*RL(1:3,1:3,ind_k)'...
                 + RL(1:3,1:3,ind_k)*robot.links(ind_k).inertia*dRj_dql';
end
end

for ind_i_q = 1:n_q

    for ind_j=1:n
        
        [J0j, Jmj] = Jacob(rL(1:3,ind_j),r0,rL,P0,pm,ind_j,robot);
        [J0j_qi, Jmj_qi] = Jacob_dq(r0,R0,P0,RJ,RL,rJ,rL, k, g, dTJ, dTL,  pm, ind_j, robot.joints(ind_i_q).id, dBi0_q, robot);
        
 
        grad_H0_q(:,:,ind_i_q) = grad_H0_q(:,:,ind_i_q) ...
                         + J0j_qi'*[Im(:,:,ind_j) zeros(3);
                                    zeros(3)  eye(3)*robot.links(ind_j).mass]*J0j ...
                         + J0j'*[diff_I_q(:,:,ind_j,ind_i_q ) zeros(3);
                                 zeros(3)      zeros(3)]*J0j...  
                         + J0j'*[Im(:,:,ind_j) zeros(3);
                                 zeros(3)  eye(3)*robot.links(ind_j).mass]*J0j_qi;


        grad_H0q_q(:,:,ind_i_q) = grad_H0q_q(:,:,ind_i_q) ...
                          + J0j_qi'*[Im(:,:,ind_j) zeros(3);
                                     zeros(3)  eye(3)*robot.links(ind_j).mass]*Jmj ...
                          + J0j'*[diff_I_q(:,:,ind_j,ind_i_q ) zeros(3);
                                  zeros(3)      zeros(3)]*Jmj...    
                          + J0j'*[Im(:,:,ind_j) zeros(3);
                                  zeros(3)  eye(3)*robot.links(ind_j).mass]*Jmj_qi;

        
        grad_Hq_q(:,:,ind_i_q) = grad_Hq_q(:,:,ind_i_q) ...
                         + Jmj_qi'*[Im(:,:,ind_j) zeros(3);
                                    zeros(3)  eye(3)*robot.links(ind_j).mass]*Jmj...
                         + Jmj'*[diff_I_q(:,:,ind_j,ind_i_q ) zeros(3);
                                 zeros(3)      zeros(3)]*Jmj...    
                         + Jmj'*[Im(:,:,ind_j) zeros(3);
                                 zeros(3)  eye(3)*robot.links(ind_j).mass]*Jmj_qi;


    end
   
end

end

