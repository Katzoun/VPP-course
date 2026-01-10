clear;clc;clf;
load ex12.mat
nr_cities = 30;

%random permutations, starting in 1
best_perm = []; best_val = inf;
for i = 1:1e3
    random_permutation = [1,1+randperm(29)];
    val = eval_perm(random_permutation,dist_matrix,nr_cities);
    if val < best_val
        best_perm = random_permutation;
        best_val = val;
    end
end

hold on;
scatter(positions(:,1),positions(:,2),MarkerFaceColor='red');
text(positions(:,1),positions(:,2),num2str((1:nr_cities)'),color='red');
plot([positions(best_perm,1);positions(1,1)],[positions(best_perm,2);positions(1,2)],'k-');
axis off;
title("Permutation with value " + num2str(best_val))

%evaluation function
function val = eval_perm(perm,dist_matrix,nr_cities)
    val = 0;
    for i=1:nr_cities-1
        cur_city = perm(i);
        next_city = perm(i+1);
        val = val + dist_matrix(cur_city,next_city);
    end
    val = val + dist_matrix(perm(nr_cities),perm(1));
end