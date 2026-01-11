clear;clc;
N = 15; M = 4; us = [0:3]'; U = length(us);
ws = [0 1 2 3]';
pw = [0.1 0.3 0.5 0.1]';

xs = [0:M]';

nr_states = (M+1)*(U);
states = zeros(nr_states,2);

iter=0;
for i=1:length(xs)
    for ii=1:length(us)
        iter = iter+1;
        states(iter,:) = [xs(i),us(ii)];
    end
end

c_1 = @(x) x^0.5;
c_2 = @(u) (u ~= 0).*(0.2 + u);
c_3 = @(x,w) -3*min(x,w);

g_k = @(x,u,w) c_1(x) + c_2(u) + c_3(x,w);
g_N = @(x) c_1(x);

J = zeros(nr_states,N+1); mu=zeros(nr_states,N);

for i = 1:nr_states 
    x_cur = states(i,1);
    J(i,N+1) = g_N(x_cur);
end

for i=N:-1:1
    for ii=1:nr_states
        state_cur = states(ii,:);
        J_temp = zeros(length(us),1);
        for iii=1:length(us)
            for iiii = 1:length(ws)
                x_cur = state_cur(1);
                u_cur = us(iii);
                w_cur = ws(iiii);
                s_cur = state_cur(2);

                x_next = min(max(x_cur - w_cur, 0) + s_cur, M);
                s_next = u_cur;

                state_next_row = find(x_next==states(:,1) & s_next==states(:,2));
                
                %average
                J_temp(iii) = J_temp(iii) + pw(iiii)*(g_k(x_cur,u_cur,w_cur)+J(state_next_row,i+1));

                %max-min
                %J_temp(iii) = max(J_temp(iii), g_k(x_cur,u_cur,w_cur)+J(state_next_row,i+1));
              
            end
        end
        [minval,minpos] = min(J_temp);
        J(ii,i) = minval; mu(ii,i)=us(minpos);
    end
end

J



%simulator
start_point = 17;
n = 10000; x0 = states(start_point,:);
J_sim = zeros(n,1); states_sim=cell(n,N+1); states_sim(:,1) = {x0};

for iter = 1:n
    for i = 1:N
        state_cur = states_sim{iter,i};
        state_row = find(state_cur(1)==states(:,1) & state_cur(2)==states(:,2));
    
        x_cur = state_cur(1);
        u_cur = mu(state_row,i);
        s_cur = state_cur(2);

        r = rand;

        idxs = find(r < cumsum(pw));
        w_cur = ws(idxs(1)); 
        
        x_next = min(max(x_cur - w_cur, 0) + s_cur, M);
        s_next = u_cur;

        states_sim(iter,i+1) = {[x_next, s_next]};

        J_sim(iter) = J_sim(iter) + g_k(x_cur,u_cur,w_cur);
   
    end
    J_sim(iter) = J_sim(iter) + g_N(states_sim{iter,end}(1)); 
end

[max(J_sim),mean(J_sim),J(start_point,1)]