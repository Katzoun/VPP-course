clear;clc;
xs = [0,1,2]; %x = [S,H,L]
zs = [1,2,3,4]; % z = [A,C,G,T]
p_S = [0.5,0.5]; % [p_SH, p_SL]
P = [0.5, 0.5;
     0.4, 0.6];
r_H = [0.2,0.3,0.3,0.2];
r_L = [0.3,0.2,0.2,0.3];
Z = [3,3,2,1,2,4,3,1,1];

D = zeros(2,10); D(:,1) = 0;
D(1,2) = -log(p_S(1)*r_H(Z(1)));
D(2,2) = -log(p_S(2)*r_L(Z(1)));

x_pred = zeros(2,10);

for i=3:10
    [D(1,i),x_pred(1,i)] = min([D(1,i-1)-log(P(1,1)*r_H(Z(i-1))),...
                                D(2,i-1)-log(P(2,1)*r_H(Z(i-1)))]);
    [D(2,i),x_pred(2,i)] = min([D(1,i-1)-log(P(1,2)*r_L(Z(i-1))),...
                                D(2,i-1)-log(P(2,2)*r_L(Z(i-1)))]);
end

x_hat = zeros(1,9); [~,x_hat(9)] = min(D(:,end));
for i=9:-1:2
    x_hat(i-1) = x_pred(x_hat(i),i+1);
end
x_hat

x_possible = zeros(2^9,9);
for i=1:2^9 
    x_as_bin = dec2bin(i-1,9);
    for j=1:9
        x_possible(i,j) = str2num(x_as_bin(j));
    end
end
x_possible = x_possible + 1;

prob = ones(2^9,1);
for i=1:2^9
    x_cur = x_possible(i,:);
    if x_cur(1) == 1
        prob(i) = prob(i)*p_S(1)*r_H(Z(1));
    else
        prob(i) = prob(i)*p_S(2)*r_L(Z(1));
    end
    for j=2:9
        if x_cur(j) == 1
            prob(i) = prob(i)*P(x_cur(j-1),x_cur(j))*r_H(Z(j));
        else
            prob(i) = prob(i)*P(x_cur(j-1),x_cur(j))*r_L(Z(j));
        end
    end
end
[maxval,maxpos] = max(prob)
[exp(-D(2,end)), maxval]
x_possible(maxpos,:)



