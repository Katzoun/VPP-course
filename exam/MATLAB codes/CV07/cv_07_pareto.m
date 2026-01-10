function [idxs] = cv_07_pareto(points)
% points - n x m matice, kde n je pocet bodu, 
%          m je pocet kriterii
% idxs - indexy (radky) kde jsou nedominovane body
[n,~] = size(points);
idxs = true(n,1);
for i=1:n
    for j=1:n
        if all(points(i,:) >= points(j,:)) && any(points(i,:) > points(j,:))
            idxs(i) = false;
            break;
        end
    end
end
end
