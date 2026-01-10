clear;close;
N = 3; Xs = [0;1;2];
Ws = [0;1;2]; ps = [0.1;0.7;0.2];
g_k = @(x,u,w) u + (x+u-w)^2;
g_N = @(x) 0;
f_k = @(x,u,w) max(0,x+u-w);

J = zeros(length(Xs),N+1);
mu = zeros(length(Xs),N);

for i=1:length(Xs)
    x_cur = Xs(i);
    J(i,N+1) = g_N(x_cur);
end

for k=N:-1:1
    for i=1:length(Xs)
        x_cur = Xs(i);
        Us = 0:2-x_cur;
        J_sub = zeros(length(Us),1);
        for ii=1:length(Us)
            u_cur = Us(ii);
            for iii=1:length(Ws)
                w_cur = Ws(iii);
                x_next = f_k(x_cur,u_cur,w_cur);
                id_x = find(Xs == x_next);
                J_sub(ii) = J_sub(ii) + ...
                    ps(iii)*(g_k(x_cur,u_cur,w_cur) + J(id_x,k+1));
            end
        end
        [minval,minpos] = min(J_sub);
        J(i,k) = minval; mu(i,k) = Us(minpos);
    end
end

%% simulation
n_sim = 10000; x0 = 2;
J_sim = zeros(n_sim,1); x_sim = zeros(n_sim,N+1);
x_sim(:,1) = x0;
for s=1:n_sim
    for k=1:N
        ids = find(rand < cumsum(ps));
        w_sim = Ws(ids(1));
        id_x_cur = find(Xs == x_sim(s,k));
        u_sim = mu(id_x_cur,k);
        J_sim(s) = J_sim(s) + g_k(x_sim(s,k),u_sim,w_sim);
        x_sim(s,k+1) = f_k(x_sim(s,k),u_sim,w_sim);
    end
end
[mean(J_sim),J(1,1)]





