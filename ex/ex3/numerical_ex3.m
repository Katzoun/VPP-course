clear;clc;

b = 1;
N = 10; 

x0 = 1
Xs = [1,2,3,4,5];

gp= [1, 0.9, 0.6, 0.4, 0]; %gp(x_k)   
ga = [0, -0.4, -1]; %ga(x_k)  

g_k @(x,u)
g_N = @(x) 0;

f_k = @(x,u,w) x-u;


N = 10; S = 5; A = 3;              % horizon, states, actions



% P(:,:,u+1) -- transition matrix for action u
P = zeros(S,S,A);

% u = 0 (do nothing)
P(:,:,1) = [0.6 0.3 0.1 0   0 ;
            0   0.7 0.2 0.1 0 ;
            0   0   0.75 0.15 0.10;
            0   0   0    0.8 0.2;
            0   0   0    0   1.0];

% u = 1 (repair)
P(:,:,2) = [1   0   0   0   0;
            0.9 0.1 0   0   0;
            0   0.9 0.1 0   0;
            0   0   0.9 0.1 0;
            0   0   0   0.9 0.1];

% u = 2 (buy new)
P(:,:,3) = repmat([1 0 0 0 0], S, 1);

J = zeros(length(Xs),N+1);
mu = zeros(length(Xs),N);

for i=1:length(Xs)
    x_cur = Xs(i);
    J(i,N+1) = g_N(x_cur);
end

for k=N:-1:1
    for i=1:length(Xs)
        x_cur = Xs(i);
        Us = Xs(Xs <= x_cur);
        J_sub = zeros(length(Us),1);
        for ii=1:length(Us)
            u_cur = Us(ii);
            x_next = f_k(x_cur,u_cur,0);

            id_x_next = find(abs(Xs - x_next) < tol);
            J_sub(ii) = (g_k(x_cur,u_cur,0) + alpha*J(id_x_next,k+1));
        end
        [minval,minpos] = min(J_sub);
        J(i,k) = minval; mu(i,k) = Us(minpos);
    end
end

%% simulation

x0 = b;
J(find(abs(Xs - x0) < tol),1) %optimalni cost-to-go z x0
u0 = mu(find(abs(Xs - x0) < tol),1);
x1 = x0 - u0;

u1 = (mu(find(abs(Xs - x1) < tol),2));
x2 = x1 - u1;

u2 = (mu(find(abs(Xs - x2) < tol),3));
x3 = x2 - u2;

% u2 = x2
[u0,u1,u2]







