function [x,v,a,t] = PartB(x0,v0,t0,xf,vf,tf)
% x0 initial movement
% v0 initial velocity
% t0 initlal time
% xf final movement
% xf final velocity
% tf final time
A = [1,t0,t0.^2,t0.^3;
    0,1,2.*t0,3.*t0.^2;
    1,tf,tf.^2,tf.^3;
    0,1,2.*tf,3.*tf.^2];
b = [x0;v0;xf;vf]
C = inv(A)*b;
C= C(end:-1:1);
t = 0:0.01:tf;
x = polyval(C,t);
v = polyval(polyder(C),t);
a = polyval(polyder(polyder(C)),t);
