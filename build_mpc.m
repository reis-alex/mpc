clear all
close all
clc
import casadi.*

addpath(genpath('D:\adosreis\Documents\Optimization'))
addpath(genpath('D:\adosreis\Desktop\researchcodes\MPC & RG\'))
%% Define system, constraint and invariant sets

A = [1 1; 0 1];
A = blkdiag(A,A);
B = [0 0; 1 0; 0 0; 0 1];
C = [1 0 0 0; 0 0 1 0];

[p,~] = size(C);
[n,m] = size(B);

Q = eye(n);
R = eye(m);
[K,P] = dlqr(A,B,Q,R); K=-K;

xbound = 10; ubound = 0.2;
Xc = Polyhedron('A',vertcat(eye(n),-eye(n)),'b',xbound*ones(2*n,1));
Uc = Polyhedron('A',vertcat(eye(m),-eye(m)),'b',ubound*ones(2*m,1));
Z  = Xc*Uc; Z.minHRep();


ZMatrix = [eye(n) zeros(n) zeros(n,m); K -K eye(m);...
          zeros(n) eye(n) zeros(n,m); zeros(m,n) zeros(m,n) eye(m)];
ZMatrix = Polyhedron('A',blkdiag(Z.A,Z.A)*ZMatrix,'b',vertcat(Z.b,0.99*Z.b));

Aw = [A+B*K -B*K B; zeros(n) eye(n) zeros(n,m); zeros(m,n) zeros(m,n) eye(m)];

Aw =  LTISystem('A',Aw);
Omega = Aw.invariantSet();
Omega = Omega.intersect(ZMatrix);
Omega.minHRep();

T = 100*P;
%%
opt.N           = 5;
opt.n_controls  = m;
opt.n_states    = n;
opt.model.type	= 'linear';
opt.model.A     = A;
opt.model.B     = B;

maxrows = 50;

opt.parameters.name = {'Xs','Us','Ref','VoroA','Vorob'};
opt.parameters.dim = [opt.n_states 0; opt.n_controls 0; opt.n_states 0; maxrows 2; maxrows 0];

opt.costs.stage.parameters = {'Xs','Us','VoroA','Vorob'};
opt.costs.stage.function = @(x,u,varargin) (x-varargin{:}(1:4))'*Q*(x-varargin{:}(1:4)) + ...
                                           (u-varargin{:}(5:6))'*R*(u-varargin{:}(5:6)) + ...
                                            100*max(0,sum([varargin{:}(7:6+maxrows) varargin{:}(6+maxrows+1:6+maxrows*2)]*x([1 3])-varargin{:}(6+maxrows*2+1:end)));
 
opt.costs.terminal.parameters = {'Xs','Us','Ref'};
opt.costs.terminal.function = @(x,varargin) (x-varargin{:}(1:4))'*P*(x-varargin{:}(1:4)) + ...
                                            (varargin{:}(1:4)-varargin{:}(7:10))'*T*(varargin{:}(1:4)-varargin{:}(7:10));

%% Define constraints
% terminal constraints
opt.constraints.terminal.set = Omega;
opt.constraints.terminal.parameters = {'Xs','Us'};

% control and state constraints
opt.constraints.polyhedral = Xc;
opt.constraints.control.upper = ubound*ones(m,1);
opt.constraints.control.lower = -ubound*ones(m,1);

% constraint on parameters
opt.constraints.parameters = {'Xs','Us'};

%% define inputs
opt.input.vector = {'Ref','Vorob'};
opt.input.matrix = {'VoroA'};

%% solver
opt.solver      = 'ipopt';
[solver,args]   = build_mpc(opt);

tmax            = 120;
N_robots        = 3;

robot(:,1,1)    = [-4.1;0;1;0];
robot(:,1,2)    = [-4.15;0;0;0];
robot(:,1,3)    = [-4.1;0;-1;0];

u0 = zeros(opt.n_controls,opt.N);
X0 = zeros(opt.n_states*opt.N,1);

colors = {'r','g','b'};

n_interp = 20;

for i = 1:N_robots
    args.x0(:,i) = [X0;reshape(u0',opt.n_controls*opt.N,1);zeros(opt.n_states,1);zeros(opt.n_states,1);zeros(opt.n_controls,1)];
end

%%
for t = 1:tmax
    
    if t<=20
        yref(:,t,1)     = [robot(1,t,1)+0.5;1];
        yref(:,t,2)     = [robot(1,t,2)+0.5;0+1.1*sind(40*t)];
        yref(:,t,3)     = [robot(1,t,3)+0.5;-1];
    else
        yref(:,t,1)     = [robot(1,t,1)+0.5;-1];
        yref(:,t,2)     = [robot(1,t,2)+0.5;0+1.1*sind(20*t)];
        yref(:,t,3)     = [robot(1,t,3)+0.5;1];
    end
    
    for j = 1:N_robots

        clear vrn_cllt
        % Get a first (unconstrained) prediction for each robot
        if t == 1
            for jj = 1:N_robots
                xs                      = pinv([A-eye(n) B; C zeros(p,m)])*[zeros(n,1);yref(:,1,jj)];
                refsimu2(:,t,jj)        = xs(1:opt.n_states);
                voro_A                  = vertcat(Xc.A(:,[1 3]),zeros(maxrows-length(Xc.A(:,[1 3])),2));
                voro_b                  = vertcat(Xc.b,zeros(maxrows-length(Xc.b),1));
                args.p                  = [robot(:,t,jj); refsimu2(:,t,jj);voro_b;reshape(voro_A,maxrows*2,1)];  
                sol                     = solver('x0', args.x0(:,j), 'lbx', args.lbx, 'ubx', args.ubx,'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
                for i = 1:opt.N
                    xpred_initial(:,i,jj) = full(sol.x(1+i*opt.n_states:(i+1)*opt.n_states));
                end            
            end
        end

        % use predictions to get voronoi; compute control
        if t == 1
            actual          = xpred_initial([1 3],2:end-1,j);
            others          = xpred_initial([1 3],2:end-1,[find([1:N_robots]-j)]);
            xs              = pinv([A-eye(n) B; C zeros(p,m)])*[zeros(n,1);yref(:,1,j)];
            refsimu(:,t,j)  = xs(1:opt.n_states);
        else
            actual          = xpred([1 3],2:end-1,j,t-1);
            others          = xpred([1 3],2:end-1,[find([1:N_robots]-j)],t-1);
            xs              = pinv([A-eye(n) B; C zeros(p,m)])*[zeros(n,1);yref(:,t,j)];
            refsimu(:,t,j)  = xs(1:opt.n_states);
        end

        [vorocells,voropoly] = voronoi_extended(actual,others,n_interp,maxrows);
        H = PolyUnion(voropoly); Hull = H.convexHull;
        xplot = sdpvar(2,1);
        Hull = poly_from_plot(Hull.A*xplot<=Hull.b,xplot,maxrows);
%         Hull = intersect(Hull,Omega.projection([1 3]));
        voro_A = Hull.A;  voro_A = vertcat(voro_A,zeros(maxrows-length(voro_A),2));
        voro_b = Hull.b;  voro_b = vertcat(voro_b,zeros(maxrows-length(voro_b),1));
        args.p = [robot(:,t,j); refsimu(:,t,j);voro_b;reshape(voro_A,maxrows*2,1)];  
        tic
        sol = solver('x0', args.x0(:,j), 'lbx', args.lbx, 'ubx', args.ubx,'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
        tsol(t) = toc;

        % get control sequence from MPC
        for i = 0:opt.N
            usol(:,i+1) = full(sol.x(opt.n_states*opt.N+1+i*opt.n_controls:opt.n_states*opt.N+i*opt.n_controls+2))';
        end

        for i = 1:opt.N
            xpred(:,i,j,t) = full(sol.x(1+i*opt.n_states:(i+1)*opt.n_states));
        end
%         xpred(:,t,j) = full(sol.x(opt.n_states+1:2*opt.n_states));

        u(:,t,j) = usol(:,1);
    
        % get artificial reference
        ya(:,t,j) = C*reshape(full(sol.x(opt.N*opt.n_states+opt.N*opt.n_controls+opt.n_states+1:opt.N*opt.n_states+opt.N*opt.n_controls+2*opt.n_states)),opt.n_states,1);
    
        robot(:,t+1,j) = A*robot(:,t,j) + B*u(:,t,j);
        y(:,t,j) = C*robot(:,t,j);
    
        args.x0(:,j) = full(sol.x); 
        clear args.p usol Hull voropoly voro_A voro_b H xs actual others
        t
    end
end


%% Plot
close all
colors = {'r','g','b'};

dist_robot = zeros(4,1);
for j = 1:length(u)
    dist_robot(1,j) = norm(robot([1 3],j,1)-robot([1 3],j,2),2);
    dist_robot(2,j) = norm(robot([1 3],j,1)-robot([1 3],j,3),2);
    dist_robot(3,j) = norm(robot([1 3],j,2)-robot([1 3],j,3),2);
end

figure
hold on
for j = 1:3
    plot(1:length(u),dist_robot(j,:),['-o' colors{j}])
end

figure(5)
hold on
plot(Xc.projection(1:2),'color','w')
plot(refsimu(1,:,1),refsimu(3,:,1),'--r')
plot(refsimu(1,:,2),refsimu(3,:,2),'--g')
plot(refsimu(1,:,3),refsimu(3,:,3),'--b')

for t = 1:tmax
    figure(5)
    plot(robot(1,1:t,1),robot(3,1:t,1),['-o' colors{1}])
    plot(robot(1,1:t,2),robot(3,1:t,2),['-o' colors{2}])
    plot(robot(1,1:t,3),robot(3,1:t,3),['-o' colors{3}])
%     axis([-5 11 -3 3])
    drawnow
    pause(0.1)
%     frame = getframe(5);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if t == 1
%          imwrite(imind,cm,'N5.gif','gif', 'Loopcount',inf,'DelayTime',0.1);
%     else
%          imwrite(imind,cm,'N5.gif','gif','WriteMode','append','DelayTime',0.1);
%     end

end

figure
hold on
for i = 1:3
    subplot(121); hold on
    stairs(1:tmax,u(1,:,i),['-' colors{i}])
    subplot(122); hold on
    stairs(1:tmax,u(2,:,i),['-' colors{i}])
end


