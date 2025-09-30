clear;clc;

M = 4; % inventory capacity
N = 15; % horizon
Xs = 0:M; % states
Ss = 0:M; % yesterday's order (arrives this evening)
Us = 0:M-1; % today's order (arrives tomorrow evening)
Ws = 0:3 ; % demand in a day k
ps = [0.1,0.3,0.5,0.1]; % disturbance probabilities

f_k = @(x,s,w) min(max(x-w,0)+s,M); %x_k+1, s_k = u_k-1
f_s = @(u) u; %s_k+1 = u_k

c1 = @(x) sqrt(x); %cost of inventory
c2 = @(u) 2*u * (u > 0); % cost of order
c3 = @(x,w) -3 * min(x,w); % profit of selling (negative cost)


g_N = @(x)  c1(x);
g_k = @(x,u,w)  c1(x) + c2(u) + c3(x,w);


J = zeros(length(Xs),length(Ss),N+1);
mu = zeros(length(Xs),length(Ss), N);

for ix=1:length(Xs)
    for is=1:length(Ss)
        J(ix,is,N+1) = g_N(Xs(ix));
    end
end


for k=N:-1:1
    for ix=1:length(Xs)
        for is=1:length(Ss)

            x_cur = Xs(ix);
            s_cur = Ss(is);
            J_sub = zeros(length(Us),1);

            for ii=1:length(Us)
                u_cur = Us(ii);
                val = 0.0;
                for iw=1:length(Ws)
                    w_cur = Ws(iw);
                    p_w = ps(iw);
                    x_next = f_k(x_cur,s_cur,w_cur);
                    s_next = f_s(u_cur);
                    ix_next = find(Xs == x_next);
                    is_next = find(Ss == s_next);

                    val = val + p_w * (g_k(x_cur,u_cur,w_cur) + J(ix_next,is_next,k+1));
                end
                J_sub(ii) = val;
            end
            [minval,minpos] = min(J_sub);
            J(ix,is,k) = minval; 
            mu(ix,is,k) = Us(minpos);
        end
    end
end

%% simulation
n_sim = 10000; x0 = 1; s0 = 0;
J_sim = zeros(n_sim,1); 
x_sim = zeros(n_sim,N+1); x_sim(:,1) = x0;
s_sim = zeros(n_sim,N+1); s_sim(:,1) = s0;


cps = cumsum(ps);

for s=1:n_sim
    for k_sim=1:N
        ix_sim = x_sim(s,k_sim) + 1;
        is_sim = s_sim(s,k_sim) + 1;
        u_sim = mu(ix_sim,is_sim,k_sim);


        ids = find(rand <= cps, 1,"first");
        w_sim = Ws(ids);
        J_sim(s) = J_sim(s) + g_k(x_sim(s,k_sim), u_sim, w_sim);
        x_sim(s,k_sim+1) = f_k(x_sim(s,k_sim), s_sim(s,k_sim), w_sim);
        s_sim(s,k_sim+1) = f_s(u_sim);
    end
    J_sim(s) = J_sim(s) + g_N(x_sim(s, N+1)); % Add the terminal cost
end

fprintf('Prumerna hodnota zisku ze simulace: %.4f \n' , mean(J_sim));
fprintf('Optimalni hodnota zisku z DP:       %.4f \n' , J(x0+1, s0+1, 1));

