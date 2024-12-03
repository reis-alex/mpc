function verbose = plot_solution(xsimu, time_vec, refsimu, robot)
%PLOT_SOLUTION Summary of this function goes here
%   Detailed explanation goes here

size_lab  = 15;
size_tl   = 13;
size_box  = 11;
size_line = 2;
%% q
refsimu_q = refsimu(7:6+robot.n_q,:);
q = xsimu(7:6+robot.n_q,:);


figure(1)
tl = 'q';
plot(time_vec(:),q(:,:),'linewidth',1);
title(tl,'Interpreter','latex','fontsize', size_tl);
xlabel('time (s)','Interpreter','latex','fontsize', size_lab);
ylabel('$\mathbf{q}$ (rad)','Interpreter','latex','fontsize', size_lab);
hold on
plot(time_vec(:), refsimu_q, '--', 'linewidth',size_line);
legend('$\mathbf{q}$','$\mathbf{q_d}$','Interpreter','latex','fontsize', size_lab);
legend('$\mathbf{q1}$','$\mathbf{q_d1}$','Interpreter','latex','fontsize', size_lab);
grid on
set(gca,'fontsize',size_lab)
set(gca,'fontweight','bold')

qdot = xsimu(12+robot.n_q+1 :12+2*robot.n_q,:);
%% q_dot
figure(2)
tl = 'qdot';
plot(time_vec(:),qdot(:,:),'linewidth',size_line)
title(tl,'Interpreter','latex','fontsize', size_tl);
xlabel('time (s)','Interpreter','latex','fontsize', size_lab);
ylabel('$\mathbf{\dot q}$ (rad/s)','Interpreter','latex','fontsize', size_lab);
legend('$ \dot q_1 $','$ \dot q_2 $', '$ \dot q_3 $', '$ \dot q_4 $', 'Interpreter','latex','fontsize', size_lab)
grid on
set(gca,'fontsize',12)
set(gca,'fontweight','bold')

%% xt angles then positions
refsimu_xt = refsimu(1:6,:);
xt = xsimu(1:6,:);

figure(3)
tl = 'angles torso';
plot(time_vec(:),xt(1:3,:),'linewidth',1)
hold on
plot(time_vec(:),refsimu_xt(1:3,:) ,'--','linewidth',size_line)
title(tl,'Interpreter','latex','fontsize', size_tl);
xlabel('time (s)','Interpreter','latex','fontsize', size_lab);
ylabel('$\mathbf{\theta_t}$','Interpreter','latex','fontsize', size_lab);
grid on
set(gca,'fontsize',12)
set(gca,'fontweight','bold')


figure(4)
tl = 'position torso';
title(tl,'Interpreter','latex','fontsize', size_tl);
plot(time_vec(:),xt(4:6,:),'linewidth',1)
xlabel('time (s)','Interpreter','latex','fontsize', size_lab);
ylabel('$\mathbf{r_t}$','Interpreter','latex','fontsize', size_lab);
%legend('$ x_t $','$ y_t $', '$ z_t $','Interpreter','latex','fontsize', size_lab);
hold on
plot(time_vec(:),refsimu_xt(4:6,:) ,'--','linewidth',size_line);
legend('$ x_t $','$ y_t $', '$ z_t $','$ x_t ^d $','$ y_t ^d $', '$ z_t ^d $','Interpreter','latex','fontsize', size_lab);
grid on
set(gca,'fontsize',12)
set(gca,'fontweight','bold')


%% xt_dot

xt_dot = xsimu(6+robot.n_q+1 :6+robot.n_q+6 ,:);

figure(5)
tl = 'Euler angles torso rates';
plot(time_vec(:),xt_dot(1:3,:),'linewidth',size_line)
title(tl,'Interpreter','latex','fontsize', size_tl);
xlabel('time (s)','Interpreter','latex','fontsize', size_lab);
ylabel('$\mathbf{\dot {\theta_t}}$','Interpreter','latex','fontsize', size_lab);
legend('$ \dot \phi $','$ \dot \theta $', '$ \dot \psi $','Interpreter','latex','fontsize', size_lab);
grid on
set(gca,'fontsize',12)
set(gca,'fontweight','bold')


figure(6)
tl = 'linear torso speed';
plot(time_vec(:),xt_dot(4:6,:),'linewidth',size_line)
title(tl,'Interpreter','latex','fontsize', size_tl);
xlabel('time (s)','Interpreter','latex','fontsize', size_lab);
ylabel('$\mathbf{ \dot r_t }$','Interpreter','latex','fontsize', size_lab);
hold on
grid on
set(gca,'fontsize',12)
set(gca,'fontweight','bold')


figure(47)
tl = 'position torso 3D';
title(tl,'Interpreter','latex','fontsize', size_tl);
plot3(refsimu_xt(4,:),refsimu_xt(5,:),refsimu_xt(6,:), '--','linewidth',size_line);
hold on
plot3(xt(4,:),xt(5,:),xt(6,:) ,'linewidth',size_line);
grid on
set(gca,'fontsize',12)
set(gca,'fontweight','bold')

verbose = 1;
end

