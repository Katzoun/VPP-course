clear;clc;
xs = (1:5)'; us = (1:3)';

% W - psti prechodu z x_k do x_k+1 pri u
W = zeros(5,5,3);
% u=1 : nedelat nic
W(1,1,1) = 0.7; W(1,2,1) = 0.2; W(1,3,1) = 0.1;
W(2,2,1) = 0.7; W(2,3,1) = 0.2; W(2,4,1) = 0.1;
W(3,3,1) = 0.7; W(3,4,1) = 0.2; W(3,5,1) = 0.1;
W(4,4,1) = 0.7; W(4,5,1) = 0.3; 
W(5,5,1) = 1;
% u=2 : oprava
W(2,1,2) = 0.9; W(2,2,2) = 0.1;
W(3,2,2) = 0.9; W(3,3,2) = 0.1;
W(4,3,2) = 0.9; W(4,4,2) = 0.1;
W(5,4,2) = 0.9; W(5,5,2) = 0.1;
% u=3 : koupit novy
W(2:5,1,3) = 1;

g_p = [1;0.9;0.7;0.5;0]; g_a = [0;-0.4;-1];

alpha = 0.9;
eps = 1e-3;

%% value iteration

J = zeros(5,1); mu = zeros(5,1);
iter = 0;
while 1
    iter = iter + 1;
    J_new = zeros(5,1);
    for ii = 1:length(xs)
        x_cur = xs(ii);
        J_temp = zeros(3,1);
        for iii = 1:length(us)
            u_cur = us(iii);
            w_cur = W(x_cur, :, u_cur);
            if sum(w_cur) > 0
                J_temp(iii) = g_p(ii) + g_a(iii) + alpha*w_cur*J;
            else
                J_temp(iii) = -Inf;
            end
        end
        [maxval,maxpos] = max(J_temp);
        J_new(ii) = maxval; mu(ii) = us(maxpos);
    end
    
    if max(abs(J-J_new)) > eps
        J = J_new;
    else
        break
    end
end

J_new, mu 

%% policy iteration

mu = ones(5,1); mu_new = zeros(5,1);
iter = 0;
while 1
    iter = iter + 1;

    P = zeros(5,5);
    g = zeros(5,1);

    for ii = 1: length(xs)
        P(ii,:) = W(ii, :, mu(ii));
        g(ii) = g_p(ii) + g_a(mu(ii));
    end
    A = eye(5)-alpha*P;
    J = A\g; 

    for ii = 1:length(xs)
        x_cur = xs(ii);
        J_temp = zeros(3,1);
        for iii = 1:length(us)
            u_cur = us(iii);
            w_cur = W(x_cur, :, u_cur);
            if sum(w_cur) > 0
                J_temp(iii) = g_p(ii) + g_a(iii) + alpha*w_cur*J;
            else
                J_temp(iii) = -Inf;
            end
        end
        [maxval,maxpos] = max(J_temp);
        mu_new(ii) = us(maxpos);
    end

    if all(mu_new == mu)
        break;
    else
        mu = mu_new;
    end
end

J, mu
