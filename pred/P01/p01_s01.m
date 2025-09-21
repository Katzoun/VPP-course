clear ;clc;

g0 = @(u) 20*sqrt(max(0,u-100));
g1 = @(u) 10*sqrt(u);
g2 = @(x) x/2;
b = 1000;
% x0 = b;
%x_k+1 = x_k - u_k

xs = [0:b]';
J = zeros(length(xs),3);
mu = zeros(length(xs),2);

% J_N(x_N) = g_N(x_N)
for i=1:length(xs)
    state = xs(i);
    J(i,3) = g2(state);
end

for i=1:length(xs)
    state = xs(i);
    us = 0:state;
    J_sub = zeros(length(us),1);
    for ii=1:length(us)
        u = us(ii);
        state_next = state-u;
        state_id = find(state_next == xs);
        J_sub(ii) = g1(u) + J(state_id,3);
    end
    [maxval,maxpos] = max(J_sub);
    J(i,2) = maxval; mu(i,2) = us(maxpos);
end

for i=1:length(xs)
    state = xs(i);
    us = 0:state;
    J_sub = zeros(length(us),1);
    for ii=1:length(us)
        u = us(ii);
        state_next = state-u;
        state_id = find(state_next == xs);
        J_sub(ii) = g0(u) + J(state_id,2);
    end
    [maxval,maxpos] = max(J_sub);
    J(i,1) = maxval; mu(i,1) = us(maxpos);
end

x0 = 1000;
J(x0+1,1)
u0 = mu(x0+1,1)
x1 = x0 - u0;
u1 = (mu(x1+1,2))
x2 = x1 - u1;
[u0,u1,x2]
g0(u0) + g1(u1) + g2(x2)




