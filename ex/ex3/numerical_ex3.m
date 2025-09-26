clear;clc;

N = 10; S = 5; A = 3;              % horizon, states, actions

% P(:,:,u+1) -- transition matrix for action u
P = zeros(S,S,A);

% u = 0 (do nothing)
P(:,:,1) = [0.60, 0.30, 0.10, 0.00, 0.00;
                 0.00, 0.70, 0.20, 0.10, 0.00;
                 0.00, 0.00, 0.75, 0.15, 0.10;
                 0.00, 0.00, 0.00, 0.80, 0.20;
                 0.00, 0.00, 0.00, 0.00, 1.00];

% u = 1 (repair)
P(:,:,2) = [1.00, 0.00, 0.00, 0.00, 0.00; 
                 0.90, 0.10, 0.00, 0.00, 0.00;
                 0.00, 0.90, 0.10, 0.00, 0.00;
                 0.00, 0.00, 0.90, 0.10, 0.00;
                 0.00, 0.00, 0.00, 0.90, 0.10];

% u = 2 (buy new)
P(:,:,3) =[1 0 0 0 0;
                1 0 0 0 0;
                1 0 0 0 0;
                1 0 0 0 0;
                1 0 0 0 0;];

Xs = [1,2,3,4,5];
Us = [0,1,2];

gp= [1, 0.9, 0.6, 0.4, 0]; %gp(x_k)   
ga = [0, -0.4, -1]; %ga(u_k)  

g_N = @(x)  gp(x);
g_k = @(x,u)  ( gp(x) + ga(u+1));

f_k = @(x,u,w) w; %x_k+1


J = zeros(length(Xs),N+1);
mu = zeros(length(Xs),N);

for i=1:length(Xs)
    x_cur = Xs(i);
    J(i,N+1) = g_N(x_cur);
end


for k=N:-1:1
    for i=1:length(Xs)
        x_cur = Xs(i);
        J_sub = zeros(length(Us),1);
        for ii=1:length(Us)
            u_cur = Us(ii);
            Ws = P(x_cur,:,u_cur+1); %u_cur+1 kvuli indexum
            J_sub(ii) = g_k(x_cur,u_cur) + dot(Ws,J(:,k+1));

        end
        [maxval,maxpos] = max(J_sub);
        J(i,k) = maxval; mu(i,k) = Us(maxpos);
    end
end

%% simulation
n_sim = 10000; x0 = 1;
J_sim = zeros(n_sim,1); x_sim = zeros(n_sim,N+1);
x_sim(:,1) = x0;
for s=1:n_sim
    for k=1:N
        x_cur_sim = x_sim(s,k);
        u_cur_sim = mu(x_cur_sim,k);
        Ws_sim = P(x_cur_sim,:,u_cur_sim+1); 
        % c =  cumsum(Ws_sim);
        % r = rand < c;
        % ids = find(r);
        %zjednodusene
        ids = find(rand < cumsum(Ws_sim));
        w_sim = Ws_sim(ids(1));

        J_sim(s) = J_sim(s) + g_k(x_cur_sim,u_cur_sim);
        x_sim(s,k+1) = ids(1); 
    end
    J_sim(s) = J_sim(s) + g_N(x_sim(s, N+1));
end

fprintf('Prumerna hodnota zisku ze simulace: %.4d \n' , mean(J_sim)); %zaporne kvuli minimalizaci
fprintf('Optimalni hodnota zisku z DP: %.4d \n' ,  J(1,1)); %zaporne kvuli minimalizaci

