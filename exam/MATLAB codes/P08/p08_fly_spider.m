clear;clc;
p = 1/4;
u = [0,1]; % 0 - not move, 1 - move
n = 7;
x = 0:n;
P_1 = zeros(n+1);
for i=3:n+1
    P_1(i,i) = p;
    P_1(i,i-1) = 1-2*p;
    P_1(i,i-2) = p;
end
P_1(2,2) = 2*p; P_1(2,1) = 1-2*p;
P_1(1,1) = 1;
P_0 = zeros(n+1);
for i=2:n+1
   P_0(i,min(i+1,n+1)) = p;
   P_0(i,i) = 1-2*p;
   P_0(i,i-1) = p;
end
P_0(end,end) = 1-p;
P_0(1,1) = 1;

J = zeros(n+1,1);
mu = zeros(n+1,1);
tol = 1e-4;
g = @(x) x~=0;
iter = 0;
while 1
    iter = iter + 1;
    J_new = J;
        for i=1:n+1
            x_cur = x(i);
            p_0 = P_0(i,:);
            p_1 = P_1(i,:);
            [minval,minpos] = min([g(x_cur)+p_0*J,g(x_cur)+p_1*J]);
            J_new(i) = minval;
            mu(i) = u(minpos);
        end
    if norm(J_new - J) < tol
        J = J_new;
        break
    end
    J = J_new;
end

[J(2), 1/(1-2*p), 1/p]

