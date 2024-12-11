function x_ee = get_end_effector_pose(xsimu, robot, time_vec)
ee_ID = 6;
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

    x_ee(:,i) = rL(:,ee_ID);
end
end