clear; clc; clf;
load ex12.mat
% dist_matrix  % 30x30 matrix of distances
% positions    % 30x2 city coordinates

nr_cities = 30;

%% ------- baseline: random search (optional) -------
best_perm_rand = [];
best_val_rand  = inf;
for i = 1:1e3
    perm = [1, 1+randperm(nr_cities-1)];
    val  = eval_perm(perm, dist_matrix, nr_cities);
    if val < best_val_rand
        best_val_rand  = val;
        best_perm_rand = perm;
    end
end
fprintf('Random search best: %.2f\n', best_val_rand);

%% ------- GA settings -------
eval_limit_GA = 1e4;
pop_size      = 50;
tour_size     = 4;
pcross        = 0.9;
pmut          = 0.3;   % prob. of applying swap mutation
hc_moves_GA   = 20;    % max 2-opt moves per child

%% ------- GA: initialization -------
evals_GA = 0;
Pop = zeros(pop_size, nr_cities);
fit = zeros(pop_size,1);

for i = 1:pop_size
    Pop(i,:) = [1, 1+randperm(nr_cities-1)];
    fit(i)   = eval_perm(Pop(i,:), dist_matrix, nr_cities);
    evals_GA = evals_GA + 1;
end

[best_val_GA, idx] = min(fit);
best_perm_GA = Pop(idx,:);

%% ------- GA main loop -------
while evals_GA < eval_limit_GA
    NewPop = zeros(size(Pop));
    new_fit = zeros(size(fit));

    for k = 1:2:pop_size
        % tournament selection
        p1 = tournament_select(Pop, fit, tour_size);
        p2 = tournament_select(Pop, fit, tour_size);

        % crossover
        if rand < pcross
            [c1, c2] = pmx_crossover(p1, p2);
        else
            c1 = p1;
            c2 = p2;
        end

        % mutation (swap on positions 2..nr_cities)
        if rand < pmut
            c1 = swap_mutation(c1);
        end
        if rand < pmut
            c2 = swap_mutation(c2);
        end

        % local 2-opt hill-climbing
        [c1, val1] = two_opt_local(c1, dist_matrix, nr_cities, hc_moves_GA);
        [c2, val2] = two_opt_local(c2, dist_matrix, nr_cities, hc_moves_GA);

        % evaluate children (count evaluations)
        val1 = eval_perm(c1, dist_matrix, nr_cities); evals_GA = evals_GA + 1;
        val2 = eval_perm(c2, dist_matrix, nr_cities); evals_GA = evals_GA + 1;

        NewPop(k,:)   = c1;
        new_fit(k)    = val1;
        if k+1 <= pop_size
            NewPop(k+1,:) = c2;
            new_fit(k+1)  = val2;
        end

        % track global best
        if val1 < best_val_GA
            best_val_GA  = val1;
            best_perm_GA = c1;
        end
        if val2 < best_val_GA
            best_val_GA  = val2;
            best_perm_GA = c2;
        end

        if evals_GA >= eval_limit_GA
            break;
        end
    end

    Pop = NewPop;
    fit = new_fit;
end

fprintf('GA best: %.2f (evals: %d)\n', best_val_GA, evals_GA);

%% ------- GRASP settings -------
eval_limit_GRASP = 1e4;
rcl_fraction     = 0.2;   % Restricted Candidate List fraction
hc_moves_GRASP   = 15;   % max 2-opt moves per solution

%% ------- GRASP main loop -------
best_val_GRASP = inf;
best_perm_GRASP = [];
evals_GRASP = 0;

while evals_GRASP < eval_limit_GRASP
    % constructive randomized greedy
    perm = grasp_construct(dist_matrix, nr_cities, rcl_fraction);
    val  = eval_perm(perm, dist_matrix, nr_cities);
    evals_GRASP = evals_GRASP + 1;

    % local improvement by 2-opt
    [perm, val] = two_opt_local(perm, dist_matrix, nr_cities, hc_moves_GRASP);

    if val < best_val_GRASP
        best_val_GRASP  = val;
        best_perm_GRASP = perm;
    end
end

fprintf('GRASP best: %.2f (evals: %d)\n', best_val_GRASP, evals_GRASP);

%% ------- plot results -------
figure;
hold on;
scatter(positions(:,1), positions(:,2), 'filled', 'MarkerFaceColor','red');
text(positions(:,1), positions(:,2), num2str((1:nr_cities)'), 'Color','red');
% GA tour
plot([positions(best_perm_GA,1); positions(best_perm_GA(1),1)], ...
     [positions(best_perm_GA,2); positions(best_perm_GA(1),2)], ...
     'b-', 'LineWidth', 1.5);
% GRASP tour
plot([positions(best_perm_GRASP,1); positions(best_perm_GRASP(1),1)], ...
     [positions(best_perm_GRASP,2); positions(best_perm_GRASP(1),2)], ...
     'g--', 'LineWidth', 1.5);
axis off;
legend({'cities','GA','GRASP'}, 'Location','bestoutside');
title(sprintf('GA best: %.2f   GRASP best: %.2f', best_val_GA, best_val_GRASP));


%% ================= helper functions =================

% evaluation function (given in skeleton)
function val = eval_perm(perm, dist_matrix, nr_cities)
    val = 0;
    for i = 1:nr_cities-1
        cur_city  = perm(i);
        next_city = perm(i+1);
        val = val + dist_matrix(cur_city,next_city);
    end
    val = val + dist_matrix(perm(nr_cities),perm(1));
end

% tournament selection (minimization)
function p = tournament_select(Pop, fit, t_size)
    n = size(Pop,1);
    idx = randi(n, t_size, 1);
    [~, best_idx_local] = min(fit(idx));
    p = Pop(idx(best_idx_local), :);
end

% swap mutation (keeps city 1 in position 1)
function c = swap_mutation(c)
    n = numel(c);
    i = randi([2 n]);
    j = randi([2 n]);
    while j == i
        j = randi([2 n]);
    end
    tmp = c(i); c(i) = c(j); c(j) = tmp;
end

% PMX crossover (Partially Mapped Crossover), cuts in 2..n
function [child1, child2] = pmx_crossover(p1, p2)
    n = numel(p1);
    child1 = zeros(1,n);
    child2 = zeros(1,n);

    % choose cut points avoiding index 1
    c1 = randi([2 n-1]);
    c2 = randi([c1 n]);
    
    % child1 from p1,p2
    child1(c1:c2) = p1(c1:c2);
    for i = c1:c2
        gene = p2(i);
        if ~ismember(gene, child1)
            pos = i;
            while true
                mapped = p1(pos);
                pos = find(p2 == mapped, 1);
                if child1(pos) == 0
                    child1(pos) = gene;
                    break;
                end
            end
        end
    end
    for i = 1:n
        if child1(i) == 0
            child1(i) = p2(i);
        end
    end

    % child2 from p2,p1 (symmetric)
    child2(c1:c2) = p2(c1:c2);
    for i = c1:c2
        gene = p1(i);
        if ~ismember(gene, child2)
            pos = i;
            while true
                mapped = p2(pos);
                pos = find(p1 == mapped, 1);
                if child2(pos) == 0
                    child2(pos) = gene;
                    break;
                end
            end
        end
    end
    for i = 1:n
        if child2(i) == 0
            child2(i) = p1(i);
        end
    end

    % enforce start in city 1
    % rotate so that city 1 is at position 1
    child1 = rotate_to_start(child1, 1);
    child2 = rotate_to_start(child2, 1);
end

function p = rotate_to_start(p, start_city)
    idx = find(p == start_city, 1);
    if idx > 1
        p = [p(idx:end), p(1:idx-1)];
    end
end

% 2-opt hill-climber (greedy, limited moves)
function [perm, val] = two_opt_local(perm, dist_matrix, nr_cities, max_moves)
    % current tour length
    val = eval_perm(perm, dist_matrix, nr_cities);
    moves = 0;
    improved = true;

    while improved && moves < max_moves
        improved = false;
        best_delta = 0;
        best_i = 0; 
        best_j = 0;

        % try all 2-opt moves (except city 1 position)
        for i = 2:nr_cities-2
            for j = i+1:nr_cities-1
                delta = two_opt_delta(perm, dist_matrix, i, j);
                if delta < best_delta
                    best_delta = delta;
                    best_i = i; 
                    best_j = j;
                    improved = true;
                end
            end
        end

        if improved
            % reverse the segment [best_i, best_j]
            perm(best_i:best_j) = perm(best_j:-1:best_i);
            val = val + best_delta;
            moves = moves + 1;
        end
    end
end


% 2-opt delta compute
function delta = two_opt_delta(perm, dist_matrix, i, j)
    n = numel(perm);
    a = perm(i-1);
    b = perm(i);
    c = perm(j);
    if j == n
        d = perm(1);
    else
        d = perm(j+1);
    end
    old_cost = dist_matrix(a,b) + dist_matrix(c,d);
    new_cost = dist_matrix(a,c) + dist_matrix(b,d);
    delta = new_cost - old_cost;
end

% GRASP construction
function perm = grasp_construct(dist_matrix, nr_cities, rcl_fraction)
    perm = zeros(1, nr_cities);
    perm(1) = 1;
    unvisited = 2:nr_cities;
    current = 1;
    for pos = 2:nr_cities
        dists = dist_matrix(current, unvisited);
        [sorted, idx] = sort(dists, 'ascend');
        rcl_size = max(1, ceil(rcl_fraction * numel(unvisited)));
        rcl_idx = idx(1:rcl_size);
        chosen_idx = rcl_idx(randi(numel(rcl_idx)));
        next_city = unvisited(chosen_idx);
        perm(pos) = next_city;
        current = next_city;
        unvisited(chosen_idx) = [];
    end
end
