clear;
close all;

dt = 0.001;
T = 0:dt:1;
alpha_z = 25;
beta_z = 6;
alpha_x = 8;
y0 = 0;
x0 =1;
z0 = 0;
yd0 = z0;
N = 10;
g = 1;
c = [ 1.0000 0.6294 0.3962 0.2494 0.1569 0.0988 0.0622 0.0391 0.0246 0.0155];
sigma_2rd = [ 41.6667 16.3934 6.5359 2.5840 1.0235 0.4054 0.1606 0.0636 0.0252 0.0252]/1000;

base_fun = zeros( N, size(T,2) );
x = zeros(size(T));
y = zeros(size(T));
yd = zeros(size(T));
ydd = zeros(size(T));

w0 = [ 0 0 0 0 0 0 0 0 0 0];
w1 = [1000 0 0 0 0 1000 0 0 0 1000];
w2 = [-1000 0 0 0 0 -1000 0 0 0 1000];
w3 = [ -200 -200 0 0 0 0 200 0 0 0];

% for question e 
data = load('imitation.data');
demo_y = data(:, 1);
demo_yd = data(:, 2);
demo_ydd = data(:, 3);
g = demo_y(end);
y0 = demo_y(1);
yd0 = demo_yd(1);
target_f = zeros( size(T,2), 1 );
Phi = zeros( size(T,2), N );

x(1) = x0;
xd = -alpha_x * x(1);
for t = 2:size(T,2)
    x(t) = x(t-1) + xd * dt;
    xd = -alpha_x * x(t);
end

for t = 1:size(T,2)
    target_f(t) = demo_ydd(t) -  alpha_z * ( beta_z * ( g-demo_y(t) ) - demo_yd(t) );
    
     for i=1:N
        base_fun(i, t) = exp( - ( x(t)-c(i) )*( x(t)-c(i) ) / ( 2 * sigma_2rd(i) )  );
     end
    for i=1:N
        Phi(t, i) = base_fun(i, t) / sum(base_fun(:, t)) * x(t) * (g-y0);
    end
end

w_reg = inv( Phi' * Phi ) * Phi' * target_f;
%%%%%%%%%%%%%%%%%%%%%%

% w0 - w3 for question c and d
% w = w0;
% w = w1;
%  w = w2;
%  w = w3;

% for question e
w = w_reg';
% -221.3576 -415.9221 -682.9674 -837.9235 -577.8125  178.3409  926.3712  891.5268  214.8619   29.8935
disp(w);


% init for t = 0
y(1) = y0;
yd(1) = yd0;
x(1) = x0;
for i=1:N
        base_fun(i, 1) = exp( - ( x(1)-c(i) )*( x(1)-c(i) ) / ( 2 * sigma_2rd(i) )  );
end
f =dot(base_fun(:, 1), w) / sum(base_fun(:, 1)) * x(1) * (g-y0);
ydd(1) = alpha_z * ( beta_z * ( g-y(1) ) - yd(1) ) + f;
xd = -alpha_x * x(1);

iter = 2;
for t = T(2:end)
    yd(iter) = yd(iter-1) + ydd(iter-1) * dt;
    y(iter) = y(iter-1) + yd(iter-1) * dt;
    x(iter) = x(iter-1) + xd * dt;
    
    for i=1:N
        base_fun(i, iter) = exp( - ( x(iter)-c(i) )*( x(iter)-c(i) ) / ( 2 * sigma_2rd(i) )  );
    end
    f =dot(base_fun(:, iter), w) / sum(base_fun(:, iter)) * x(iter) * (g-y0);
    ydd(iter) = alpha_z * ( beta_z * ( g-y(iter) ) - yd(iter) ) + f;
    
    xd = -alpha_x * x(iter);
    
    iter =iter + 1;
end

figure(1);
hold on;
% plot(T, y);
% for question e
h1 = plot(T, y);
h2 = plot(T, demo_y, 'r');
legend([h1,h2], 'y', 'y\_demo');
title('y');

figure(2);
hold on;
% plot(T, yd);
% for question e
h1 = plot(T, yd);
h2 = plot(T, demo_yd, 'r');
legend([h1,h2], 'yd', 'yd\_demo');
title('yd');

figure(3);
hold on;
% plot(T, ydd);
% for question e
h1 = plot(T, ydd);
h2 = plot(T, demo_ydd, 'r');
legend([h1,h2], 'ydd', 'ydd\_demo');
title('ydd');

% % for question c
% figure(4);
% plot(T, x);
% title('x');
% %  
% figure(5); 
% hold on;
% for i = 1:N
%     plot(T, base_fun(i,:) );
% end
% title('10 Phi functions');




