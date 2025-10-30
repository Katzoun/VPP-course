%% Ex8 LQR
clear; clc; close all;

dt = 0.1;       % step
gamma = 0.1;   
N = 50;         

%state vector xk = [x; vx; y; vy; g]
A = [ 1, dt, 0,  0,  0;
      0, 1 - dt*gamma, 0,  0,  0;
      0, 0, 1, dt,  0;
      0, 0, 0, 1 - dt*gamma, -dt;
      0, 0, 0,  0,  1];

B = [ 0, 0;
      dt, 0;
      0, 0;
      0, dt;
      0, 0 ];

nx = size(A,1);
nu = size(B,2);
lambda = 0.005;
sigma = diag([0.005, 0.002, 0.005, 0.001, 0]); % covariance matrix of noise


Q = zeros(nx);              % Q_k 
QN = (1 - lambda)*eye(nx);  % Q_N cena za konecny stav
R = lambda * eye(nu);       % R_k ceny za rizeni

K = cell(N+1,1);
L = cell(N,1);

K{N+1} = QN;

for k = N:-1:1
    K{k} = A'*(K{k+1} - K{k+1}*B*inv(B'*K{k+1}*B+ R)*B'*K{k+1})*A + Q;
    L{k} = - inv(B'*K{k+1}*B + R)*B'*K{k+1}*A;
end

x0 = [2; 0; 2; 0; 9.8];  % pocatecni stav [x, vx, y, vy, g]


% optimal value
J = x0'*K{1}*x0;
for k = 2:N
    J = J + trace(K{k}*sigma);
end

fprintf('optimal value of J = %.6f\n', J);


%% simulation

xs = zeros(nx, N+1);
us = zeros(nu, N);
xs(:,1) = x0;

for k = 1:N
    us(:,k) = L{k} * xs(:,k);
    xs(:,k+1) = A*xs(:,k) + B*us(:,k);
end


%% visualization
t = 0:dt:N*dt;
figure;
subplot(3,1,1);
plot(t, xs(1,:), 'LineWidth', 1.5); hold on;
plot(t, xs(3,:), 'LineWidth', 1.5);
xlabel('time [s]');
ylabel('position [m]');
legend('x','y');
grid on
title('Drone Position');

subplot(3,1,2);
plot(t, xs(2,:), 'LineWidth', 1.5); hold on;
plot(t, xs(4,:), 'LineWidth', 1.5);
xlabel('time [s]');
ylabel('velocity [m/s]');
legend('v_x','v_y');
grid on
title('Velocities');

subplot(3,1,3);
stairs(t(1:end-1), us(1,:), 'LineWidth', 1.5); hold on;
stairs(t(1:end-1), us(2,:), 'LineWidth', 1.5);
xlabel('time [s]');
ylabel('control');
legend('u_x','u_y');
grid on
title('Optimal Control');

sgtitle('Ex8');

figure(2)
plot(xs(1,:),xs(3,:), 'LineWidth', 1.5);
xlabel('position in x [m]');
ylabel('position in y [m]');
grid on
title('Drone Position');

