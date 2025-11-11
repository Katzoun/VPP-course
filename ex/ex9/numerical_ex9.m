%% Ex9 – Machine Repair Scheduling (Infinite Horizon)
clear; clc; close all;

%% ===== Parameters =====
alpha = 0.9;          % discount factor
S = 5;                % number of machine states
A = 3;                % number of actions

Xs = 1:S;             % {1=new, 2=good, 3=normal, 4=bad, 5=broken}
Us = 0:A-1;           % {0=do nothing, 1=repair, 2=buy new}

%% ===== Transition probabilities =====
% P(:,:,u+1): transition matrix for each action
P = zeros(S,S,A);

% Action u = 0 (do nothing)
P(:,:,1) = [0.7, 0.2, 0.1, 0.0, 0.0;
             0.0, 0.7, 0.2, 0.1, 0.0;
             0.0, 0.0, 0.7, 0.2, 0.1;
             0.0, 0.0, 0.0, 0.7, 0.3;
             0.0, 0.0, 0.0, 0.0, 1.0];

% Action u = 1 (repair)
P(:,:,2) = [0.1, 0.9, 0.0, 0.0, 0.0;
             0.0, 0.1, 0.9, 0.0, 0.0;
             0.0, 0.0, 0.1, 0.9, 0.0;
             0.0, 0.0, 0.0, 0.1, 0.9;
             0.0, 0.0, 0.0, 0.0, 1.0];

% Action u = 2 (buy new)
P(:,:,3) = [1, 0, 0, 0, 0;
             1, 0, 0, 0, 0;
             1, 0, 0, 0, 0;
             1, 0, 0, 0, 0;
             1, 0, 0, 0, 0];

%% ===== Rewards =====
gp = [1.0, 0.9, 0.7, 0.5, 0.0];     % profit from machine state
ga = [0.0, -0.4, -1.0];             % cost of actions
g_k = @(x,u) -(gp(x) + ga(u+1));       % total immediate reward (minimization)

%% ===== Value Iteration =====
tol = 1e-6;
J = zeros(S,1);     % initial value function
mu = zeros(S,1);    % optimal policy
iter = 0;

while true
    J_prev = J;
    for i = 1:S
        x = Xs(i);
        Q = zeros(A,1);
        for u = 1:A
            Q(u) = g_k(x,Us(u)) + alpha * (P(i,:,u) * J_prev);
        end
        [J(i), mu(i)] = min(Q);
    end
    iter = iter + 1;
    if max(abs(J - J_prev)) < tol
        break;
    end
end

fprintf('Value iteration converged in %d iterations.\n', iter);
disp('Optimal value function J*:');
disp(J);
disp('Optimal policy μ* (1=do nothing, 2=repair, 3=buy new):');
disp(mu');

%% ===== Policy Iteration =====
mu_pi = ones(S,1);   % initial policy (do nothing)
stable = false;
iter_pi = 0;

I = eye(S);

while ~stable
    iter_pi = iter_pi + 1;
    
    % --- Policy Evaluation ---
    A_pi = I;
    b_pi = zeros(S,1);
    for i = 1:S
        A_pi(i,:) = I(i,:) - alpha * P(i,:,mu_pi(i));
        b_pi(i) = g_k(Xs(i), Us(mu_pi(i)));
    end
    J_pi = A_pi \ b_pi;

    % --- Policy Improvement ---
    stable = true;
    for i = 1:S
        Q = zeros(A,1);
        for u = 1:A
            Q(u) = g_k(Xs(i), Us(u)) + alpha * (P(i,:,u) * J_pi);
        end
        [~, new_u] = min(Q);
        if new_u ~= mu_pi(i)
            stable = false;
        end
        mu_pi(i) = new_u;
    end
end

fprintf('\nPolicy iteration converged in %d iterations.\n', iter_pi);
disp('Optimal value function Vπ*:');
disp(J_pi);
disp('Optimal policy μπ* (1=do nothing, 2=repair, 3=buy new):');
disp(mu_pi');

%% ===== Visualization =====
% figure;
% subplot(2,1,1);
% bar(Xs, J_pi, 'FaceColor', [0.2 0.6 0.8]);
% title('Optimal Value Function J^*');
% xlabel('Machine State');
% ylabel('Value');
% grid on;
% 
% subplot(2,1,2);
% bar(Xs, mu_pi, 'FaceColor', [0.8 0.4 0.2]);
% title('Optimal Policy μ^*');
% xlabel('Machine State');
% ylabel('Action (1=do nothing, 2=repair, 3=buy new)');
% grid on;
