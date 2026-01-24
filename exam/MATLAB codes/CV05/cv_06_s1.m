clear;clc;clf;
load cv_06a.mat;
%load cv_06b.mat;
problem.start_node = 1;
end_node = problem.end_node;

node_list = problem.node_list;
node_neighbors = problem.node_neighbors;
neighbors_distance = problem.neighbors_distance;
nr_nodes = length(node_list);
M = problem.M;
imshow(M,'InitialMagnification',1200);











