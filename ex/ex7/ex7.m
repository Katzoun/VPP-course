clear; close all; clc

N = 10; 
W = 12;
Xs = 0:W;

w = [2; 3; 7; 1; 5; 4; 3; 5; 6; 2];
c = [5; 6; 11; 3; 7; 7; 6; 7; 8; 5];
p = [0.95; 0.93; 0.7; 0.96; 0.94; 0.96; 0.95; 0.91; 0.9; 0.95];

Us = [0,1];   % decision take or not take

g1_k = @(k,u) -c(k)*u;            
g1_N = @(x) 0;               

g2_k = @(k,u) u*p(k) + (1-u);   
g2_N = @(x) 1 - 0.15*(W - x)/W;

f_k = @(x,u,k) x - w(k)*u;     

JF = cell(length(Xs),N+1);


for i = 1:length(Xs)
    x_cur = Xs(i);
    JF{i,N+1} = [g1_N(x_cur), -log(g2_N(x_cur))];
end

%% --- Dynamic Programming ---
for k = N:-1:1
    for i = 1:length(Xs)
        x_cur = Xs(i);
        new_solutions = [];

        for ii = 1:length(Us)
            u_cur = Us(ii);

            % Feasibility check
            if (u_cur == 1) && (x_cur < w(k))
                continue;
            end

            % next state
            x_next = f_k(x_cur,u_cur,k);
            id_x_next = x_next + 1;

            next_solutions = JF{id_x_next, k+1};

            if isempty(next_solutions)
                continue;
            end

            for j = 1:size(next_solutions,1)
                J1 = g1_k(k,u_cur) + next_solutions(j,1);
                J2 = -log(g2_k(k,u_cur)) + next_solutions(j,2); %logaritmus pro soucet
                new_solutions = [new_solutions; J1, J2];
            end
        end

        % Keep only non-dominated solutions
        JF{i,k} = pareto_prune(new_solutions);
    end
end

%% --- Result ---
PF = JF{length(Xs),1};   % Pareto front for initial state (x0 = W)

profit = -PF(:,1);
prob_catch = exp(-PF(:,2));

% Sort by profit
[profit, idx] = sort(profit,'descend');
prob_catch = prob_catch(idx);

fprintf('Pareto front (non-dominated solutions):\n');
fprintf('  # | Profit | P(not caught)\n');
for i = 1:length(profit)
    fprintf('%3d | %6.2f |   %9.6f\n', i, profit(i), prob_catch(i));
end

figure;
plot(1 - prob_catch, profit, 'bo-','LineWidth',1.2);
xlabel('Probability of being caught');
ylabel('Profit');
title('Pareto Front – Ex7');
grid on;




%% --- Pareto pruning helper ---
function PF = pareto_prune(M)
% Removes dominated solutions (min J1, min J2)
if isempty(M)
    PF = [];
    return;
end

J1 = M(:,1);
J2 = M(:,2);

keep = true(size(M,1),1);
for i = 1:size(M,1)
    if ~keep(i), continue; end
    for j = 1:size(M,1)
        if i ~= j
            if (J1(j) <= J1(i)) && (J2(j) <= J2(i)) ...
                    && ((J1(j) < J1(i)) || (J2(j) < J2(i)))
                keep(i) = false;
                break;
            end
        end
    end
end
PF = M(keep,:);
end
