% Problem 2 – Value Iteration and Policy Iteration
clear; clc;

load("ex_2025_12_09_P.mat");   % loads variable P

% Parameters:
states   = 1:9;      % feasible states
controls = 1:10;     % feasible controls
alpha    = 0.99;     % discount factor
tol      = 1e-5;     % convergence criterium

numX = length(states);
numU = length(controls);

% Cost function g(x,u) = 40*x + u^2
g = zeros(numX, numU);
for x = states
    for u = controls
        g(x, u) = 40*x + u^2;
    end
end

%% VALUE ITERATION
J = zeros(numX,1);   % initial guess
diff = inf;

while diff > tol
    Jnew = zeros(numX,1);

    for x = states
        bestVal = inf;

        for u = controls

            % Transition probability row P(x,:,u) reshaped to 9×1
            px = reshape(P(x, :, u), [numX, 1]);

            % Expected future value
            EJ = px' * J;

            % Bellman update
            val = g(x,u) + alpha * EJ;

            bestVal = min(bestVal, val);
        end

        Jnew(x) = bestVal;
    end

    diff = max(abs(Jnew - J));
    J = Jnew;
end

J_value_iteration = J;

%% POLICY ITERATION
% Initialize arbitrary policy
mu = randi([1 numU], numX, 1);
policyStable = false;

while ~policyStable
    
    A = eye(numX);         % Identity matrix I
    b = zeros(numX,1);     % Vector of immediate costs g_μ(x)

    for x = states
        u = mu(x);         % current policy action μ(x)

        % A(x,:) = I(x,:) - α * P(x,:,μ(x))
        A(x,:) = A(x,:) - alpha * reshape(P(x,:,u), [1, numX]);

        % Immediate cost: g(x, μ(x)) = 40*x + μ(x)^2
        b(x) = g(x,u);
    end

    % Produces J*(x) under current policy μ
    Jp = A \ b;

    %  POLICY IMPROVEMENT STEP
    policyStable = true;   % assume stable unless a state changes its action

    for x = states
        bestVal = inf;
        bestU = mu(x);     % current policy action

        for u = controls
            % Transition probability vector P(x,:,u)
            px = reshape(P(x, :, u), [numX, 1]);

            % Expected value of future costs
            EJ = px' * Jp;

            % Q(x,u) = g(x,u) + α * Σ_j P(x,j,u) * J_μ(j)
            val = g(x,u) + alpha * EJ;

            % Select control minimizing the Bellman expression
            if val < bestVal
                bestVal = val;
                bestU = u;
            end
        end

        % If the best action changes, policy is not yet stable
        if bestU ~= mu(x)
            mu(x) = bestU;
            policyStable = false;
        end
    end
end

% After convergence:
J_policy_iteration = Jp;  % μ = optimal stationary policy μ*(x)
mu_policy_iteration = mu; % Jp = optimal cost-to-go J*(x)

fprintf("Optimal cost-to-go J*(x) from VALUE ITERATION:\n");
disp(J_value_iteration')

fprintf("Optimal cost-to-go J*(x) from POLICY ITERATION:\n");
disp(J_policy_iteration')

fprintf("Optimal policy μ*(x):\n");
disp(mu_policy_iteration')

