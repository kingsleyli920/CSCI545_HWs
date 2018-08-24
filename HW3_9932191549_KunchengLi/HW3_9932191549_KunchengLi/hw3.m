% Question 1
close all; 
clear;

% a
fs = 100;
fc = 5;
[b,a] = butter(2,fc/(fs/2));
fprintf('b1=%f, b2=%f, b3=%f\n', b(1), b(2), b(3));
fprintf('a1=%f, a2=%f, a3=%f\n', a(1), a(2), a(3));

% b
filename = 'noisy.data';
delimiter=' ';
data=importdata(filename,delimiter);
yn=data(:,1);
xn=data(:,2);
un=data(:,3);
xf = filter(b, a, yn);
figure(1);
subplot(2,1,1)
plot(xn);
ylabel('xn','rotation',0);
subplot(2,1,2)
plot(xf);
ylabel('xf','rotation',0);
figure(2);
p1 = plot(xn,'b');
hold on;
p2 = plot(xf, 'r');
legend([p1,p2], 'Unfiltered values','Filtered values');
hold off;

disp("delay of Filtered values and Unfiltered values is: " + finddelay(xn,xf));

% c
length = 1000;
A = 0.9;
B = 0.5;
C = 1;
R = 1; % the covariance of the observation noise;
Q = 0.01;% the covariance of the process noise;
P = zeros(1, length);
K = zeros(1, length);
x_post = zeros(1, length);
x_pri = 0;  
p_pri = Q; % Assume P(0) = 1;
for t = 1:1000
    K(t) = p_pri * C' / ( C * p_pri * C' + R ); % Kn = pPriori * C' / (C * pPriori * C' + R);
    x_post(t) = x_pri + K(t) * ( yn(t) - C * x_pri ); % xPosteriori = xPriori + K * (yn - C * xPriori);
    P(t) = ( 1 - K(t) * C ) * p_pri; % Pn = (1 - K(n) * C) * pPriori;
    x_pri = A * x_post(t) + B * un(t); % xPriori = A * xPosteriori_n + B * un;   
    p_pri = A * P(t) * A' + Q; %pPriori = A * Pn * A' + Q;
end

figure(3);
subplot(2,1,1);
plot(xn);
ylabel('xn','rotation',0);
subplot(2,1,2);
plot(x_post);
ylabel('xf','rotation',0);

figure(4);
plot(K);
title('Gain K');

figure(5);
plot(P);
title('Posterior Covariance Matrix P');

figure(6);
p3 = plot(xn,'b');
hold on;
p4 = plot(x_post, 'r');
legend([p3,p4], 'Unfiltered values','Filtered values');
hold off;

disp("delay of Filtered values and Unfiltered values is: " + finddelay(xn,x_post));


