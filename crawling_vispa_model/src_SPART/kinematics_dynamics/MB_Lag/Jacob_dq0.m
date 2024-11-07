function [J0i_dq0,Jmi_dq0] = Jacob_dq0(r0,R0,P0,q,RJ,RL,rJ,rL, k,g, dTJ0, dTL0, pm, Bij, Bi0, ind_i, coord_i, dT0_dq0, dBi0_dq0, q0, robot )
% Computes the geometric Jacobian base-derivative of a point `i` of coord_i-th base term derivative.
% returns J0i_dq0 -- Base Link Jacobian partial-derivative -- as a [6x6] matrix
% returns Jmi_dq0 -- Manipulator Jacobian partial-derivative -- as a [6xn_q] matrix

% :parameters:
%   * rp -- Position of the point of interest, projected in the inertial CCS -- [3x1].
%   * tp -- Twist of the point of interest [\omega,rdot], projected in the intertial CCS -- [6x1].
%   * r0 -- Position of the base-link center-of-mass with respect to the origin of the inertial frame, projected in the inertial CCS -- [3x1].
%   * t0 -- Base-link twist [\omega,rdot], projected in the inertial CCS -- as a [6x1] matrix.
%   * Bi0 -- Twist-propagation matrix (for i>0 and j=0) -- as a [6x6xn] matrix.
%   * pm -- Manipulator twist-propagation "vector" -- as a [6xn] matrix.
%   * i -- Link id where the point `p` is located -- int 0 to n.
%   * coord_i -- q0 coord to derivate on.
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).

%Omega0 = q0dot(1:3);


%Pre-allocate
Jmi_dq0= zeros(6,robot.n_q, 'like', R0);
J0i_dq0= zeros(6,6, 'like', R0);
dpm_dq0 = zeros(6,1, 'like', R0);
dk_dq0  = zeros(3,1, 'like', R0);
e       = zeros(3,1, 'like', R0);
RJ_part = zeros(3, 'like', R0);

[RJ_q0,RL_q0,rJ_q0,rL_q0,k_q0,g_q0]  = dT_parser(dTJ0, dTL0, robot, coord_i, R0);

%Manipulator Jacobian
joints_num=0;

dP0_dq0i = [dT0_dq0(1:3,1:3,coord_i) zeros(3);
    zeros(3)            zeros(3)];
%Bi0_dq0 = [zeros(3,6); SkewSym(dT0_dq0(1:3,4,coord_i)-rL_q0(1:3,i)), zeros(3)];
J0i_dq0 = [eye(3), zeros(3,3); SkewSym(r0-rL(1:3,ind_i)), eye(3)]*dP0_dq0i+dBi0_dq0*P0;

for j=1:ind_i
    %If joint is not fixed
    if robot.joints(j).type~=0
        %eq 56 annexe
        dk_dq0 = k_q0(:,j);
        dg_dq0 = g_q0(:,j);
        %eq 54 annexe
        if robot.con.branch(ind_i,j)==1 %&& i>j
            if robot.joints(j).type==1 %revolute
                dpm_dq0 = [dk_dq0; cross(dk_dq0, g(:,j)) + cross(k(:,j), dg_dq0)];%dispersion possible ici
            elseif robot.joints(j).type==2 %prismatic
                dpm_dq0 = [zeros(3,1); dk_dq0 ];
            end
            %eq 53 annexe
            col_jacob = [eye(3), zeros(3,3); SkewSym(rL(1:3,j)-rL(1:3,ind_i)), eye(3)]*dpm_dq0 + ...
                [zeros(3,6); SkewSym(rL_q0(1:3,j)-rL_q0(1:3,ind_i)), zeros(3)]*pm(:,j);

            Jmi_dq0(1:6,robot.joints(j).q_id) = col_jacob;
        end
        joints_num=joints_num+1;
    end
end
% end

end

