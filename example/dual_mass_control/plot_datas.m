function [out] = plot_datas(t,x)
%PLOT_DATAS Summary of this function goes here

% plot the solution
figure(1)
subplot(4,1,1)
plot(t,x(:,1),'b','linewidth',1)
xlabel('t')
ylabel('x(t)')
grid on
set(gca,'fontsize',12)
set(gca,'fontweight','bold')
subplot(4,1,2)
plot(t,x(:,2),'b','linewidth',1)
xlabel('t')
ylabel("x'(t)")
grid on
set(gca,'fontsize',12)
set(gca,'fontweight','bold')
subplot(4,1,3)
plot(t,x(:,3),'b','linewidth',1)
xlabel('t')
ylabel('y(t)')
grid on
set(gca,'fontsize',12)
set(gca,'fontweight','bold')
subplot(4,1,4)
plot(t,x(:,4),'b','linewidth',1)
xlabel('t')
ylabel("y'(t)")
xlabel('time')
grid on
set(gca,'fontsize',12)
set(gca,'fontweight','bold')

% Below, plot the phase portraits of the trajectories
figure(2)
plot(x(:,1),x(:,2),'b','linewidth',1)
xlabel('x(t)')
ylabel("x'(t)")
grid on
set(gca,'fontsize',12)
set(gca,'fontweight','bold')
figure(3)
plot(x(:,3),x(:,4),'b','linewidth',1)
xlabel('y(t)')
ylabel("y'(t)")
grid on
set(gca,'fontsize',12)
set(gca,'fontweight','bold')


end

