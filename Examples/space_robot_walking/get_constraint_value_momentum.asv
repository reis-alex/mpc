function [h_robot] = get_constraint_value_momentum(xsimu, robot, time_vec)
%GET_CONSTRAINT_VALUE_MOMENTUM Summary of this function goes here
%   Detailed explanation goes here



for i = 1:length(time_vec)

    xt  = xsimu(1:6, i);
    rt = xt(4:6);
    Rt = Angles321_DCM(xt(1:3));
    
    phi = xt(1);
    theta = xt(2);
    psi = xt(3);

    T_plate = [1 cos(phi)*tan(theta)  sin(phi)*tan(theta);
    0 cos(phi)            -sin(phi);
    0 sin(phi)/cos(theta)  cos(phi)/cos(theta)];
   
    q = xsimu(7: 6+robot.n_q, i);
    qt_dot_raw = xsimu(6+robot.n_q+1 :6+robot.n_q+6 ,i);
    
    euler_rate = qt_dot_raw(1:3);
    omega_t = T_plate\euler_rate;
    qt_dot = [omega_t; qt_dot_raw(4:6)];
    qdot = xsimu(12+robot.n_q+1 :12+2*robot.n_q,i);

    %Kinematics
    [RJ,RL,rJ,rL,e,g]=Kinematics(Rt,rt,q,robot);
    %Diferential Kinematics
    [Bij,Bi0,Pt,pm]=DiffKinematics(Rt,rt,rL,e,g,robot);
    %Velocities
    [t_t,tq]=Velocities(Bij,Bi0,Pt,pm,qt_dot,qdot,robot);
    % %Jacobian of the last link
    % [J0n, Jmn]=Jacob(rL(1:3,end),r0,rL,P0,pm,robot.n_links_joints,robot);
    N    = NOC(rt,rL,Pt,pm,robot);
    Ndot = NOCdot(rt,t_t,rL,tq,Pt,pm,robot);
    %--- Dynamics ---%
    %Inertias in inertial frames
    [It,Iq]=I_I(Rt,RL,robot);

    %Generalized Inertia matrix
    H = GIM_NOC(N,It,Iq,robot);
    %Generalized Convective Inertia matrix
    C = CIM_NOC(N,Ndot,t_t,tq,It,Iq,robot);
    % H = [Ht, Htq; Htq.', H];
    % C = [Ct, Ctq; Ctq  , Cm];
    Ht = H(1:6,1:6);
    Htq = H(1:6, 7:end);

    h_robot(1:6,i) = Ht*qt_dot + Htq*qdot;
end

end

