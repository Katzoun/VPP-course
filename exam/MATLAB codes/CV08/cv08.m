clear;clc;clf;
dt = 0.1; N = 50; gamma = 0.2;
lambda = 0.0001;
sigma = 1*diag([0.005, 0.002, 0.005, 0.001, 0]);
x0 = [2, 0, 2, 0, 9.8]';
A = [1 dt 0 0 0;
     0 1-dt*gamma 0 0 0;
     0 0 1 dt 0;
     0 0 0 1-dt*gamma -dt;
     0 0 0 0 1];
B = [0 0;dt 0;0 0;0 dt;0 0];
%Q_k = zeros(5,5);
Q_k = (1-lambda)*eye(5)/N;
Q_N = (1-lambda)*eye(5);
R_k = lambda*eye(2);

K = zeros(5,5,N+1); K(:,:,N+1) = Q_N;
L = zeros(2,5,N);
for k=N:-1:1
    temp = inv(B'*K(:,:,k+1)*B+R_k);
    L(:,:,k) = -temp*B'*K(:,:,k+1)*A;
    K(:,:,k) = A'*(K(:,:,k+1) - ...
       K(:,:,k+1)*B*temp*B'*K(:,:,k+1))*A+Q_k;
end

J0 = x0'*K(:,:,1)*x0;
for k=1:N
    J0 = J0 + trace(K(:,:,k)*sigma);
end

%% simulace
xs = zeros(5,N+1); xs(:,1) = x0;
us = zeros(2,N);
Cost = 0;
for k=1:N
    us(:,k) = L(:,:,k)*xs(:,k);
    Cost = Cost + xs(:,k)'*Q_k*xs(:,k) + ...
        us(:,k)'*R_k*us(:,k);
    xs(:,k+1) = A*xs(:,k) + B*us(:,k) + ...
        mvnrnd(zeros(5,1),sigma)';
end
Cost = Cost + xs(:,end)'*Q_N*xs(:,end)

subplot(2,1,1)
hold on; plot(xs(1,:),xs(3,:),'b-x');
plot(0,0,'rx')
axis equal; grid on;
subplot(2,1,2)
hold on;
plot(us(1,:)); plot(us(2,:)); grid on;