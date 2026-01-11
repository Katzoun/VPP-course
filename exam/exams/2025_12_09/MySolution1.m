%% PROBLEM 1 - Dynamic Programming (finite horizon)
clear; clc;

N = 20;
c = 12;
d = 8;

xs = (1:c)';          % states
us = (1:d)';          % actions
ws = (1:6)';          % disturbances
pw = ones(6,1)/6;     % equal probability

% stage cost and terminal cost
g  = @(x,u,w) x + 2*u - 5*min(x,5);
gN = @(x) 3*x + min(x,3);

% system equation
f = @(x,u,w) min(max(x + 2*u - 2*w - 3, 1), c);

J  = zeros(c, N+1);      % J(x,k) stored as J(x,k+1) with x in 1..c
mu = zeros(c, N);        % optimal action for each state & stage

% terminal cost
for x = 1:c
    J(x, N+1) = gN(x);
end

% backward DP
for k = N:-1:1
    for x = 1:c
        Q = inf(d,1);
        for ui = 1:d
            u = us(ui);
            expVal = 0;
            for wi = 1:6
                w = ws(wi);
                xnext = f(x,u,w);
                expVal = expVal + pw(wi) * ( g(x,u,w) + J(xnext, k+1) );
            end
            Q(ui) = expVal;
        end
        [J(x,k), bestIdx] = min(Q);
        mu(x,k) = us(bestIdx);
    end
end

x0 = 2;
fprintf('J*(x0=2) = %.12f\n', J(x0,1));
