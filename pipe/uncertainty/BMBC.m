tsta = 3800;
tend = 4140;
% data = load('pipe.out');
data = importdata('pipe.out');
t = data(:,1);
for k = 1:size(t,1)
    if t(k) > tsta
        ista = k;
        break;
    end
end
for k = size(t,1):-1:1
    if t(k) < tend
        iend = k;
        break;
    end
end
pgrad = data(ista:iend,5);
ub = data(ista:iend,4);
tauw = -pgrad/2;
ff = 8*tauw./ub.^2;

x = ff;
M = 14000;
N = floor((iend-ista+1)/M)*M;
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
format long
disp([s1/s0,N*var2])
sqrt(var2)/mean(ff,'all')