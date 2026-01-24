%% LQR via HJB (Riccati ODE) – simple example
% State: x = [x1; x2], dynamics: x1_dot = x2, x2_dot = u
% Cost: integral (x1^2 + x2^2 + u^2) dt

clear; clc; close all;

%% System data
A = [0 1;
     0 0];
B = [0;
     1];

Q = eye(2);   % state cost
R = 1;        % input cost (scalar)

% Analytic infinite-horizon solution from theory
K_analytic = [sqrt(3) 1;
              1       sqrt(3)];

%% Finite-horizon Riccati ODE (HJB)
% V(t,x) = x'K(t)x  =>  -Kdot = A'K + KA - KBR^-1B'K + Q
% We integrate backward from t = T to t = 0.

T  = 10;          % finite horizon (large enough to approximate infinite)
QT = eye(2);      % terminal weight K(T)

K_T_vec = QT(:);  % vectorized terminal condition

% Solve Riccati ODE backward in time
[tK, K_vec] = ode45(@(t,K) riccati_ode(t, K, A, B, Q, R), [T 0], K_T_vec);

% Take K(0) from the last point
K0 = reshape(K_vec(end,:).', 2, 2);

disp('K(0) from Riccati ODE:');
disp(K0);

disp('Analytic K*:');
disp(K_analytic);

disp('Difference K0 - K_analytic:');
disp(K0 - K_analytic);

%% Closed-loop simulation with u = -R^-1 B'' K0 x
Acl = A - B*(R\B')*K0;   % A_cl = A - B R^-1 B' K0
x0  = [-2; 10];          % initial state

Tsim = 10;               % simulation horizon
[tx, x] = ode45(@(t,x) Acl*x, [0 Tsim], x0);

% Control along the trajectory
u = zeros(size(tx));
for k = 1:length(tx)
    u(k) = - (R\B')*K0*x(k,:).';
end

%% Plots
figure;
subplot(2,1,1);
plot(tx, x(:,1), 'LineWidth', 1.2); hold on;
plot(tx, x(:,2), 'LineWidth', 1.2);
grid on;
xlabel('t [s]');
ylabel('states');
legend('x_1','x_2','Location','best');
title('Closed-loop states');

subplot(2,1,2);
plot(tx, u, 'LineWidth', 1.2);
grid on;
xlabel('t [s]');
ylabel('u(t)');
title('Optimal control');

%% --- Local functions ---

function dKdt = riccati_ode(~, Kvec, A, B, Q, R)
    % Riccati ODE for K(t) in vectorized form
    n = size(A,1);
    K = reshape(Kvec, n, n);

    % -Kdot = A'K + KA - K B R^-1 B' K + Q
    term = A'*K + K*A - K*B*(R\B')*K + Q;
    dK   = -term;

    dKdt = dK(:);
end
