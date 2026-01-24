clear;clc;
P = [3,1;1,4;5,1;2,2;4,3];
M = [2,1,4,5,3]; n = 5;
r = length(M);
X = 1:n;
K = setdiff(X,M); P1 = P(K,:);

[cost,times] = p04_schedule_cost(P,M);

[~,sort_ind] = sort(P1(:,1));
S1 = 0;
for k=r+1:n
    for i=M
        S1 = S1 + P(i,1);
    end
    S1 = S1 + (n-k+1)*P1(sort_ind(k-r),1) + ...
        P1(k-r,2);
end
%S1

[~,sort_ind] = sort(P1(:,2));
Cr = times(end,end);
temp = max(Cr, sum(P(M,1)) + min(P(K,1)));
S2 = 0;
for k=r+1:n
    S2 = S2 + temp + (n-k+1)*P1(sort_ind(k-r),2);
end
%S2

cost + max(S1,S2)

%% brute force
all_p = perms(1:5);
all_cost = zeros(factorial(5),1);
for i=1:factorial(5)
    all_cost(i) = p04_schedule_cost(P,all_p(i,:));
end
[minval,minpos] = min(all_cost)
all_p(93,:)

