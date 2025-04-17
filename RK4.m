function x = RK4(X,U,f,opt)
M = 1;
h = opt.dt/M;
xtemp = X(:,1);
for i = 1:M
    k1 = f(xtemp,U(:,1));
    k2 = f(xtemp+h/2*k1,U(:,1));
    k3 = f(xtemp+h/2*k2,U(:,1));
    k4 = f(xtemp+h*k3,U(:,1));
    xtemp   = xtemp+ (h/6)*(k1+2*k2+2*k3+k4);
end
x = xtemp;