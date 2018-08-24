clear
clc
t0 = 0;
tf = 2;
x0 = 0;
xf = 2;
v0 = 0;
vf = 0;
[x,v,a,t] = PartB(x0,v0,t0,xf,vf,tf);
figure
plot(t,x)
title('x');
xlabel('t');
ylabel('x');
figure
plot(t,v)
title('v')
xlabel('t');
ylabel('v');
figure
plot(t,a);
title('a');
xlabel('t');
ylabel('a')