clear; clc; clf; close all
load cv_06a.mat;
% load cv_06b.mat;

problem.start_node = 1;
end_node = problem.end_node;
node_list = problem.node_list;
node_neighbors = problem.node_neighbors;
neighbors_distance = problem.neighbors_distance;
nr_nodes = length(node_list);
M = problem.M;


figure;
imshow(M, 'InitialMagnification',1200);


%% A* call
fprintf('A* algorithm started...\n');
tic;
[path, path_cost, nodes_expanded] = astar_algorithm(node_list, node_neighbors, neighbors_distance, problem.start_node, end_node);
computation_time = toc;

%% Results
if  not(isempty(path))
    disp('Path found!');
    fprintf('Path cost: %.2f\n', path_cost);
    fprintf('Path length: %d\n', length(path));
    fprintf('No. evaluated nodes: %d\n', nodes_expanded);
    fprintf('Computation time: %.4f s\n', computation_time);

    hold on;
    path_coords = node_list(path, :);
    plot(path_coords(:, 2), path_coords(:, 1), 'r-', 'LineWidth', 2);
else
    disp('Path not found!');
end

hold off;




%% A* Algoritmus Implementation

function [path, path_cost, nodes_expanded] = astar_algorithm(node_list, node_neighbors, neighbors_distance, start_node, end_node)
    nr_nodes = length(node_list);
    
    % Init
    open_set = start_node;
    closed_set = [];
    
    % g_score - distance from start
    g_score = inf(nr_nodes, 1);
    g_score(start_node) = 0;
    
    % f_score - g_score + heuristic
    f_score = inf(nr_nodes, 1);
    f_score(start_node) = heuristic(node_list(start_node, :), node_list(end_node, :));
    
    % came_from - pro rekonstrukci cesty
    came_from = zeros(nr_nodes, 1);

    nodes_expanded = 0;
    
    while ~isempty(open_set)

        % Najdi node s nejmensim f_score v open_set
        [~, min_idx] = min(f_score(open_set));
        current = open_set(min_idx);
        
        % check jestli neni node cilovy
        if current == end_node
            path = reconstruct_path(came_from, current);
            path_cost = g_score(end_node);
            return;
        end
        
        % presun current z open_set do closed_set
        open_set(min_idx) = [];
        closed_set = [closed_set, current];
        nodes_expanded = nodes_expanded + 1;
        
        % Projdi vsechny sousedy
        neighbors = node_neighbors{current};
        distances = neighbors_distance{current};
        
        for i = 1:length(neighbors)
            neighbor = neighbors(i);
            
            % Preskoc pokud je soused uz v closed_set
            if ismember(neighbor, closed_set)
                continue;
            end
            
            temp_g_score = g_score(current) + distances(i);
            
            % Pokud soused neni v open_set, add 
            if ~ismember(neighbor, open_set)
                open_set = [open_set, neighbor];
            
            elseif temp_g_score >= g_score(neighbor)
                continue; % Tahle cesta neni lepsi, skip 
            end
            
            % tato cesta je nejlepsi, save 
            came_from(neighbor) = current;
            g_score(neighbor) = temp_g_score;
            f_score(neighbor) = g_score(neighbor) + heuristic(node_list(neighbor, :), node_list(end_node, :));
        end
    end
    
    % Cesta nebyla nalezena
    path = [];
    path_cost = inf;
end

function h = heuristic(pos1, pos2)
    % euclid
%     h = sqrt(sum((pos1 - pos2).^2));

    %manhattan vykazuje na mrizce lepsi vlastnosti
    h = abs(pos1(1) - pos2(1)) + abs(pos1(2) - pos2(2));

end

function path = reconstruct_path(came_from, current)
    % Rekonstrukce cesty z came_from
    path = current;
    while came_from(current) ~= 0
        current = came_from(current);
        path = [current, path];
    end
end
