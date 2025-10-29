%% Ex8 – Linear Quadratic Regulator (LQR)
clear; clc; close all;

%% --- Parametry systému ---
dt = 0.1;       % časový krok
gamma = 0.1;    % tlumení
N = 50;         % horizont

% Stavový vektor: [x; vx; y; vy; g]
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

% Rozměry
nx = size(A,1);
nu = size(B,2);

%% --- Váhy LQR ---
lambda = 0.005;
Q = zeros(nx);         % Q_k
QN = (1 - lambda)*eye(nx); % Q_N
R = lambda * eye(nu);  % R_k

%% --- Inicializace Riccatiho rekurze ---
P = cell(N+1,1);
K = cell(N,1);

P{N+1} = QN;   % terminální podmínka

%% --- Backward Riccati recursion ---
for k = N:-1:1
    P_next = P{k+1};
    K{k} = (R + B'*P_next*B) \ (B'*P_next*A);
    P{k} = Q + A'*P_next*A - A'*P_next*B*K{k};
end

%% --- Simulace optimálního řízení ---
x0 = [5; 0; 10; 0; 9.81];  % počáteční stav [x, vx, y, vy, g]
x = zeros(nx, N+1);
u = zeros(nu, N);
x(:,1) = x0;

for k = 1:N
    u(:,k) = -K{k} * x(:,k);
    x(:,k+1) = A*x(:,k) + B*u(:,k);
end

%% --- Výpočet celkové hodnoty kritéria J ---
J = 0;
for k = 1:N
    J = J + x(:,k)'*Q*x(:,k) + u(:,k)'*R*u(:,k);
end
J = J + x(:,N+1)'*QN*x(:,N+1);

%% --- Výsledky ---
fprintf('Celková hodnota kritéria J = %.6f\n', J);

%% --- Vizualizace ---
t = 0:dt:N*dt;
figure;
subplot(3,1,1);
plot(t, x(1,:), 'LineWidth', 1.5); hold on;
plot(t, x(3,:), 'LineWidth', 1.5);
xlabel('čas [s]');
ylabel('poloha [m]');
legend('x','y');
title('Pozice dronu');

subplot(3,1,2);
plot(t, x(2,:), 'LineWidth', 1.5); hold on;
plot(t, x(4,:), 'LineWidth', 1.5);
xlabel('čas [s]');
ylabel('rychlost [m/s]');
legend('v_x','v_y');
title('Rychlosti');

subplot(3,1,3);
stairs(t(1:end-1), u(1,:), 'LineWidth', 1.5); hold on;
stairs(t(1:end-1), u(2,:), 'LineWidth', 1.5);
xlabel('čas [s]');
ylabel('řízení');
legend('u_x','u_y');
title('Optimální řízení (LQR)');

sgtitle('LQR Regulátor – Ex8');
