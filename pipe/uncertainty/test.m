N = 1e6;
% R = rand(N,1);
% save('random_R.mat','R')
load('random_R.mat','R');
x = zeros(N,1);
x(1) = 0.0;
for n = 2:N
    x(n) = 0.9*x(n-1) + 0.1*R(n-1);
end
plot(x)
M = 100;
K = N/M;

x_ave = zeros(K,1);
mu_hat = mean(x,'all');
for k = 1:K
    xk = x((k-1)*M+1:k*M)-mu_hat;
    x_ave(k) = mean(xk,'all');
end

s0 = sum(x_ave.^2,'all');
s1 = sum(x_ave(1:K-1).*x_ave(2:K),'all');

var2 = 1/((K-1)*(K-2))*(s0+2*s1);
disp([s1/s0,N*var2])


