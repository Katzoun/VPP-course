%% PROBLEM 2 - Value Iteration + Policy Iteration (discounted cost)
clear; clc;

alpha = 0.99;
eps   = 1e-6;

% Load transition probabilities: P(x, xNext, u) = 9x9x10
S = load('ex_2025_12_09_P.mat');   % uprav název souboru podle sebe
P = S.P;

% Dimensions
[S1,S2,A] = size(P);
assert(S1==9 && S2==9 && A==10, 'Expected P to be 9x9x10.');

% Optional sanity checks (doporučuji)
tol = 1e-10;
assert(all(P(:) >= -tol), 'P contains negative values.');
rowSums = squeeze(sum(P,2));                 % 9x10 (sum over xNext)
assert(all(abs(rowSums(:) - 1) < 1e-8), 'Some rows of P do not sum to 1.');

% Cost function: g(x,u) = 40*x + u^2
g = @(x,u) 40*x + u^2;

%% (A) VALUE ITERATION
J  = zeros(9,1);
mu = ones(9,1);

iter = 0;
while true
    iter = iter + 1;
    Jnew = zeros(9,1);

    for x = 1:9
        Q = inf(10,1);
        for u = 1:10
            pRow = squeeze(P(x,:,u));    % 1x9
            Q(u) = g(x,u) + alpha*(pRow*J);
        end
        [Jnew(x), mu(x)] = min(Q);
    end

    if norm(Jnew - J, inf) < eps
        J = Jnew;
        break;
    end
    J = Jnew;
end

fprintf('Value Iteration finished in %d iterations.\n', iter);
disp('J* (VI):');  disp(J);
disp('mu* (VI):'); disp(mu);

%% (B) POLICY ITERATION
mu = ones(9,1);
iter = 0;

while true
    iter = iter + 1;

    % Policy evaluation: solve (I - alpha*P_mu) J = g_mu
    Pmu = zeros(9,9);
    gmu = zeros(9,1);
    for x = 1:9
        Pmu(x,:) = squeeze(P(x,:,mu(x)));
        gmu(x)   = g(x, mu(x));
    end
    J = (eye(9) - alpha*Pmu) \ gmu;

    % Policy improvement
    mu_new = zeros(9,1);
    for x = 1:9
        Q = inf(10,1);
        for u = 1:10
            pRow = squeeze(P(x,:,u));
            Q(u) = g(x,u) + alpha*(pRow*J);
        end
        [~, mu_new(x)] = min(Q);
    end

    if all(mu_new == mu)
        break;
    end
    mu = mu_new;
end

fprintf('Policy Iteration finished in %d iterations.\n', iter);
disp('J* (PI):');  disp(J);
disp('mu* (PI):'); disp(mu);
