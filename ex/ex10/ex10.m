%% LQR via HJB / Riccati
% x = [x1; x2], x1_dot = x2, x2_dot = u
% J = ∫ (x1^2 + x2^2 + u^2) dt

clear; clc;

%% System and cost
A = [0 1;
     0 0];
B = [0;
     1];
Q = eye(2);
R = 1;

x0 = [-2; 10];   % initial state

%% Analytic infinite-horizon solution (known from theory)
K_analytic = [sqrt(3) 1;
              1       sqrt(3)];

Jstar_analytic = x0' * K_analytic * x0;

%% Riccati ODE (finite horizon approximation)
T  = 10;        % horizon
QT = eye(2);    % terminal cost

K_T_vec = QT(:);

% solve backward from T to 0
[tK, K_vec] = ode45(@(t,K) riccati_ode(t, K, A, B, Q, R), [T 0], K_T_vec);

K0 = reshape(K_vec(end,:).', 2, 2);   % K(0) from ODE

Jstar_num = x0' * K0 * x0;

disp('K(0) from Riccati:');
disp(K0);
disp('K analytic:');
disp(K_analytic);
disp('J* analytic:');
disp(Jstar_analytic);
disp('J* numeric:');
disp(Jstar_num);

%% Closed-loop simulation with K0
Acl = A - B*(R\B')*K0;

Tsim = 10;
[t, x] = ode45(@(t,x) Acl*x, [0 Tsim], x0);

u = zeros(size(t));
for k = 1:length(t)
    u(k) = -(R\B')*K0*x(k,:).';
end

%% Plots
figure;
subplot(2,1,1);
plot(t, x(:,1), 'LineWidth', 1.1); hold on;
plot(t, x(:,2), 'LineWidth', 1.1);
grid on;
xlabel('t');
ylabel('x');
legend('x_1','x_2');

subplot(2,1,2);
plot(t, u, 'LineWidth', 1.1);
grid on;
xlabel('t');
ylabel('u');

%% ---- local function ----
function dKdt = riccati_ode(~, Kvec, A, B, Q, R)
    % Riccati ODE in vector form: -Kdot = A'K + KA - KBR^-1B'K + Q
    n = size(A,1);
    K = reshape(Kvec, n, n);

    term = A'*K + K*A - K*B*(R\B')*K + Q;
    dK   = -term;

    dKdt = dK(:);
end
