clc;close;clear;
%%
g        = 9.8; 
L        = 1.0;
tinitial = 0;
tfinal   = 5.0;
X0       = [L*cos(0*pi/180),L*sin(0*pi/180),0.0,0.0]';
deltat   = 0.05;
tspan    = tinitial:deltat:tfinal;
params.g = g;
params.L = L;
dyn_pendulum = @(t,x) pendulumfixedlengthdyn(t,x,params);
options = odeset('RelTol',1e-12,'AbsTol',1e-10);
[T,X]    = ode45(dyn_pendulum,tspan,X0,options);
dlmwrite('myFile.txt',X','delimiter','\t');
