clear
clc
close all
x0 = 0;
xf = 2;
t0 = 0;
tf = 2;
v0 = 0;
vf = 0;
x(1) = x0;
v(1) = v0;
detT = 0.01;
t = t0:detT:tf;
for i = 2:length(t)-1
    [ x1,v1] = PartC(x0,v0,xf,vf,tf);
    x(i) = x1;
    v(i) = v1;
    x0 = x1;
    v0 = v1;
    tf = tf - detT;
end
x(i+1) = xf;
v(i+1) = vf;
figure
plot(t,x);
title('x')
xlabel('t')
figure
plot(t,v);
title('v')
xlabel('t')
axis([0,2,0,1.5])
figure
plot(t(2:end),diff(v)./detT)
title('a')
xlabel('t')