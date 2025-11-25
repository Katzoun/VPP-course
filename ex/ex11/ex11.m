clear; clc; close all;

n = 20;                 % number of variables
eval_limit = 1e4;       % budget

mu = 20;
lambda = 60;
p_mut = 1/n;

%% hill climbing

evals = 0;
x = randi([0 1],1,n);      % random start
fx = objective(x); 
evals = evals + 1;

best_val_HC = fx;
best_sol_HC = x;

while evals < eval_limit
    % pick one variable to flip
    cand = x;
    idx = randi(n);
    cand(idx) = 1 - cand(idx);

    f_new = objective(cand); 
    evals = evals + 1;

    if f_new >= fx
        x = cand;
        fx = f_new;

        if fx > best_val_HC
            best_val_HC = fx;
            best_sol_HC = x;
        end
    end
end

fprintf('\n---- Hill Climbing ----\n');
fprintf('Best value = %d\n', best_val_HC);
disp(best_sol_HC);


%% (mu,lambda)-ES

% mu = 15;
% lambda = 60;
% p_mut = 1/n;

% initialize parents
Population = randi([0 1], mu, n);
fitP = zeros(mu,1);

evals = 0;
for i = 1:mu
    fitP(i) = objective(Population(i,:));
    evals = evals + 1;
end

[best_val_ES1, idx] = max(fitP);
best_sol_ES1 = Population(idx,:);

while evals < eval_limit
    % Offspring
    Off = zeros(lambda, n);
    fitC = zeros(lambda,1);

    for j = 1:lambda
        parent = Population(randi(mu), :);
        child = parent;

        % bitwise mutation
        mut_mask = rand(1,n) < p_mut;
        child(mut_mask) = 1 - child(mut_mask);

        Off(j,:) = child;
        fitC(j) = objective(child);
        evals = evals + 1;
        if evals >= eval_limit
            break;
        end
    end

    % select best mu offspring
    [~, order] = sort(fitC, 'descend');
    Population = Off(order(1:mu), :);
    fitP = fitC(order(1:mu));

    if fitP(1) > best_val_ES1
        best_val_ES1 = fitP(1);
        best_sol_ES1 = Population(1,:);
    end
end

fprintf('\n---- (mu, lambda)-ES ----\n');
fprintf('Best value = %d\n', best_val_ES1);
disp(best_sol_ES1);


%% (mu+lambda)-ES

% mu = 15;
% lambda = 60;
% p_mut = 1/n;

% init parents
Parents = randi([0 1], mu, n);
fitP = zeros(mu,1);

evals = 0;
for i = 1:mu
    fitP(i) = objective(Parents(i,:));
    evals = evals + 1;
end

[best_val_ES2, idx] = max(fitP);
best_sol_ES2 = Parents(idx,:);

while evals < eval_limit

    % offspring
    Off = zeros(lambda, n);
    fitC = zeros(lambda,1);

    for j = 1:lambda
        par = Parents(randi(mu), :);
        mut = par;

        mut_mask = rand(1,n) < p_mut;
        mut(mut_mask) = 1 - mut(mut_mask);

        Off(j,:) = mut;
        fitC(j) = objective(mut);
        evals = evals + 1;
        if evals >= eval_limit
            break;
        end
    end

    % combine parents + offspring
    All = [Parents; Off];
    fitAll = [fitP; fitC];

    [~, ord] = sort(fitAll, 'descend');

    Parents = All(ord(1:mu), :);
    fitP = fitAll(ord(1:mu));

    if fitP(1) > best_val_ES2
        best_val_ES2 = fitP(1);
        best_sol_ES2 = Parents(1,:);
    end
end

fprintf('\n---- (mu + lambda)-ES ----\n');
fprintf('Best value = %d\n', best_val_ES2);
disp(best_sol_ES2);


%% obj func
function val_sum = objective(x)
    val = zeros(45,1);
    weights = 1:45;
    val(1) = x(1) | not(x(2)) | not(x(20));
    val(2) = x(1) | x(17) | not(x(9));
    val(3) = not(x(1)) | x(2) | x(10);
    val(4) = x(2) | not(x(3)) | not(x(19));
    val(5) = x(3) | x(4) | not(x(5));
    val(6) = x(3) | not(x(7)) | x(11);
    val(7) = x(3) | x(9) | not(x(15));
    val(8) = x(4) | x(7) | not(x(16));
    val(9) = not(x(4)) | not(x(18)) | x(10);
    val(10) = x(4) | not(x(20)) | x(11);
    val(11) = x(5) | not(x(6));
    val(12) = not(x(5)) | not(x(14));
    val(13) = not(x(5)) | not(x(11));
    val(14) = not(x(5)) | x(6);
    val(15) = x(6) | x(7);
    val(16) = x(6) | x(14);
    val(17) = x(6) | not(x(15));
    val(18) = x(7) | x(20);
    val(19) = x(7) | not(x(11));
    val(20) = x(8) | x(13);
    val(21) = not(x(8)) | x(10);
    val(22) = x(9) | not(x(12));
    val(23) = not(x(9)) | not(x(10));
    val(24) = x(10) | x(15);
    val(25) = x(10) | x(16);
    val(26) = not(x(10)) | x(14);
    val(27) = x(11) | x(17);
    val(28) = not(x(11)) | x(20);
    val(29) = not(x(11)) | not(x(19));
    val(30) = x(12) | not(x(18));
    val(31) = x(12) | x(19);
    val(32) = x(12) | x(15);
    val(33) = not(x(13)) | not(x(14));
    val(34) = not(x(13)) | x(16);
    val(35) = x(13) | x(15);
    val(36) = x(14) | x(1);
    val(37) = x(15) | not(x(2));
    val(38) = x(16) | not(x(17)) | x(8);
    val(39) = x(17) | x(18) | not(x(19));
    val(40) = not(x(18)) | not(x(19)) | not(x(5));
    val(41) = not(x(10)) | x(11) | x(12);
    val(42) = not(x(6)) | x(20) | not(x(4));
    val(43) = not(x(2)) | x(19) | not(x(6));
    val(44) = not(x(14))| x(18) | not(x(2));
    val(45) = not(x(3)) | x(11) | x(14);
    val_sum = sum(weights*val);
end
