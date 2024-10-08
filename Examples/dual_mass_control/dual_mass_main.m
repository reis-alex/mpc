%% example with spring-damper two-mass system control using MPC and momentum transfer limitation cost

clear
clear global
clc
close all
global m1 m2 k1 k2 a2

%% parameters definitions

%masses
m1 = 10; %kg
m2 = 2; %kg

%stiffness coef
k1 = 100; %N/m
k2 = 100; %N/m

% damping coef
a1 = 0.2; %N.s/m
a2 = 0.2; %N.s/m


% State space model

A=[[            0,     1,      0,     0];
    [-(k1 + k2)/m1, -a2/m1,  k2/m1,  a2/m1];
    [            0,     0,      0,     1];
    [        k2/m2,  a2/m2, -k2/m2, -a2/m2]];

B = [ 0 k1/m1 0 0].';
C = [ 0,0,    1,0];
D=0;

sys=ss(A,B,C,D);
step(sys,40)
xlabel('t')
ylabel('y(t)')
grid on
set(gca,'fontsize',12)
set(gca,'fontweight','bold')
%% Solve numerically using ODE
clear
clear global
global m1 m2 k1 k2 a2

%masses
m1 = 10; %kg
m2 = 2; %kg

%stiffness coef
k1 = 100; %N/m
k2 = 100; %N/m

% damping coef
a1 = 0.2; %N.s/m
a2 = 0.2; %N.s/m


x0=[0;0;0;0];
[t,x]=ode45(@two_mass_model,0:0.01:40,x0);

plot_datas(t,x);

function dx = two_mass_model(t,x)
global m1 m2 k1 k2 a2

x=[x(1);x(2);x(3);x(4)];
u=1;
dx1=x(2);
dx2=-((k1+k2)/m1)*x(1)-(a2/m1)*x(2)+(k2/m1)*x(3)+(a2/m1)*x(4)+(k1/m1)*u;
dx3=x(4);
dx4=(k2/m2)*x(1)+(a2/m2)*x(2)-(k2/m2)*x(3)-(a2/m2)*x(4);
dx=[dx1;dx2;dx3;dx4];
end