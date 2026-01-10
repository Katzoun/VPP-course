clear;clc;clf;
a = 1; b = 1; q = 5; x0 = 1; T = 10;

xi = (x0*exp(a*T))/(exp(-a*T)/q  - ...
    (b^2 / (2*a))*(exp(-a*T) - exp(a*T)));

x_opt = @(t) x0*exp(a*t)+(xi*b^2/(2*a))*(exp(-a*t) - exp(a*t));
u_opt = @(t) -xi*b*exp(-a*t);

t = linspace(0,T,100);

sysfun = @(x,y) p09_sys(x,y,a,b);
syscon = @(x,y) p09_cond(x,y,x0,q);
solinit = bvpinit(linspace(0,T,50),p09_guess)

sol = bvp4c(sysfun,syscon,solinit)

subplot(2,1,1)
plot(t,x_opt(t),'k-'); hold on;
plot(sol.x,sol.y(1,:),'ro')

subplot(2,1,2)
plot(t,u_opt(t),'k-');hold on;
plot(sol.x,-b*sol.y(2,:),'ro')