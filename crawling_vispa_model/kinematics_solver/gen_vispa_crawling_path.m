%constrained dynamics URDF Tutorial

%--- Clean and clear ---%
clc
close all
clear

%--- URDF filename ---%

robot_name = 'VISPA_crawling.urdf';
SCname = 'pulsar_PL_SP_CMG.urdf';
%--- Create robot model ---%
[robot,robot_keys] = urdf2robot_flex_visu(robot_name);
%--- Create spacecraft model ---%
[SC,SC_keys] = urdf2robot_flex_visu(SCname);

%[display_rb,display_sc]=display_sc_robot(robot,SC,[0 -1 0;-1 0 0;0 0 1], [1 0 3]',zeros(12,1));
%q_l=[0 0 pi/2 0 -pi/2+pi/3 0]';
q_r=-[0 0 pi/2 0 -pi/2+pi/3 0]';
q_l= [0.4 -0.2 -0.98 3.35 0.36 0]';
%q_r=-q_l;
Rt=eye(3);
gen_path=false;
if ~gen_path
    load('path.mat');
    qstep_list=qpath(:,1);
    [RJ,RL,rJ,rL,e,g]=Kinematics(Rt,rtpath(:,1),qstep_list,robot);
    [r_com] = Center_of_Mass(rtpath(:,1),rL,robot);
    [~,ra]=mgd_rel(qstep_list,6,robot);
    [~,rb]=mgd_rel(qstep_list,12,robot);
    disp_handle=display_robot(robot,qpath(:,1),eye(3),rtpath(:,1));
    %disp_handle=display_robot(robot,rtpath(:,1),eye(3),qpath(:,1));
    %step_nb=ceil(max(abs(diff(qpath')'))/0.1);
    step_nb=ones(1,35);
    xlim([-1.5 0.5]);
    ylim([-2 2]);
    zlim([-1 1]);
    for it=2:size(qpath,2)
        step_q=(qpath(:,it)-qpath(:,it-1))/step_nb(it-1);
        step_r=(rtpath(:,it)-rtpath(:,it-1))/step_nb(it-1);
        for ind=1:step_nb(it-1)
            q_step=qpath(:,it-1)+step_q*ind;
            rt_step=rtpath(:,it-1)+step_r*ind;
            [~,ra_s]=mgd_rel(q_step,6,robot);
            [~,rb_s]=mgd_rel(q_step,12,robot);
            [RJ,RL,rJ,rL,e,g] = Kinematics(Rt,rt_step,q_step,robot);
            [r_com_s] = Center_of_Mass(rt_step,rL,robot);
            r_com=[r_com r_com_s];
            ra=[ra ra_s];
            rb=[rb rb_s];
            qstep_list=[qstep_list q_step];
            update_display_robot(disp_handle,robot,Rt,rt_step,q_step);
            pause(0.2);
        end
    end
    plot3(ra(1,:),ra(2,:),ra(3,:),'-or','LineWidth',2)
    plot3(rb(1,:),rb(2,:),rb(3,:),'-og','LineWidth',2)
    plot3(r_com(1,:),r_com(2,:),r_com(3,:),'-ob','LineWidth',2)
else
    qpath=[];
    disp_handle=display_robot(robot,eye(3),zeros(3,1),[q_l;q_r]);
    %% Pre step pose
    %solve left leg
    rd_l1=[-1;1;-0.4];
    rd_l2=[-1.2;1;-0.4];
    % [Rp,rp]=mgd_rel([q_l;q_r],6,robot);
    % Rd_l=[0 0 -1;Rp(2:3,1:2)/norm(Rp(2:3,1:2)) zeros(2,1)];
    % [q_init,e]=iterative_ik(Rd_l,rd_l,[q_l;q_r],robot,6,1:6);
    [q_mid,e]=iterative_ik_pos(rd_l1,rd_l2,[q_l;q_r],robot,6,1:5);
    fprintf('Left leg %f %f %f %f %f %f\n',e(1),e(2),e(3),e(4),e(5),e(6));
    %update_display_robot(disp_handle,robot,Rt,[0 -1 0]',q_mid);
    %solve right leg
    rd_r1=[-1;0;-0.4];
    rd_r2=[-1.2;0;-0.4];
    %[Rp,rp]=mgd_rel(q_init,12,robot);
    %Rd_r=[0 0 -1;Rp(2:3,1:2)/norm(Rp(2:3,1:2)) zeros(2,1)];
    %[q_init,e]=iterative_ik(Rd_r,rd_r,q_init,robot,12,7:12);
    [q_mid,e]=iterative_ik_pos(rd_r1,rd_r2,q_mid,robot,12,7:11);
    fprintf('Right leg %f %f %f %f %f %f\n',e(1),e(2),e(3),e(4),e(5),e(6));
    update_display_robot(disp_handle,robot,Rt,[0;-1;0],q_mid);
    qpath=[qpath q_mid];
    pause;
    %% Step right
    rd_r1=[-0.9;0;-0.4];
    rd_r2=[-1.1;0;-0.4];
    [q_mid,e]=iterative_ik_pos(rd_r1,rd_r2,q_mid,robot,12,7:11);
    update_display_robot(disp_handle,robot,Rt,[0;-1;0],q_mid);
    qpath=[qpath q_mid];
    pause(0.2);
    for dy=0.1:0.1:1
        rd_r1=[-0.9;-dy;-0.4];
        rd_r2=[-1.1;-dy;-0.4];
        [q_mid,e]=iterative_ik_pos(rd_r1,rd_r2,q_mid,robot,12,7:11);
        update_display_robot(disp_handle,robot,Rt,[0;-1;0],q_mid);
        qpath=[qpath q_mid];
        pause(0.2);
    end
    rd_r1=[-1;-1;-0.4];
    rd_r2=[-1.2;-1;-0.4];
    [q_mid,e]=iterative_ik_pos(rd_r1,rd_r2,q_mid,robot,12,7:11);
    %fprintf('Right leg %f %f %f %f %f %f\n',e(1),e(2),e(3),e(4),e(5),e(6));
    update_display_robot(disp_handle,robot,Rt,[0;-1;0],q_mid);
    qpath=[qpath q_mid];
    pause;
    %% Step left by symmetry
    q_lp=[-qpath(7:12,13).*ones(1,13) -qpath(7:12,13:-1:1)];
    q_rp=[qpath(7:12,:) qpath(7:12,end).*ones(1,13)];
    qpathb=[q_lp;q_rp];
    %% Torso to right
    qt=qpathb(:,end);
    qtr=[];
    for dy=0.1:0.1:1
        rd_r1=[-1;-1+dy;-0.4];
        rd_r2=[-1.2;-1+dy;-0.4];
        [qt,e]=iterative_ik_pos(rd_r1,rd_r2,qt,robot,12,7:11);
        qtr=[qtr qt(7:12,1)];
    end
    qtlr=[-qtr(:,10:-1:1) qtr];
    rtpath=zeros(3,36);
    rtpath(2,27:end)=[-0.1:-0.1:-1];
end