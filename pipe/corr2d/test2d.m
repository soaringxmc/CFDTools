clear
close all
clc

% FFT2
len1 = 2*pi;
len2 = 4*pi;
m1 = 5;
m2 = 7;
m1m = m1-1;
m2m = m2-1;
m1mh = m1m/2+1;
m2mh = m2m/2+1;
dx = len1/m1m;
dy = len2/m2m;
x1 = (0:m1m-1)*dx;
y1 = (0:m2m-1)*dy;
x = zeros(m1m,m2m);
y = zeros(m1m,m2m);
% meshgrid avoided for clarity
for i1 = 1:m1m
    for i2 = 1:m2m
        x(i1,i2) = dx*(i1-1);
        y(i1,i2) = dy*(i2-1);
    end
end
f = cos(x)+2*sin(y)+cos(x.*y);
fluc = f-mean(f,'all');
figure
contourf(x,y,f)
tmp = fft2(f);
fhat = tmp/(m1m*m2m);
fhat2 = abs(fhat).^2;
fvar = sum(fhat2,'all')-fhat2(1,1);

% symmetric property
kap1 = zeros(m1m,1);
kap2 = zeros(m2m,1);
for k = 1:m1m
    if k <= m1mh
        kap1(k) = (k-1)*2*pi/len1;
    else
        kap1(k) = -(m1m-k+1)*2*pi/len1;
    end
end
for k = 1:m2m
    if k <= m2mh
        kap2(k) = (k-1)*2*pi/len2;
    else
        kap2(k) = -(m2m-k+1)*2*pi/len2;
    end
end
kap1msh = zeros(m1m,m2m);
kap2msh = zeros(m1m,m2m);
r1msh = zeros(m1m,m2mh);
r2msh = zeros(m1m,m2mh);
for k1 = 1:m1m
    for k2 = 1:m2m
        kap1msh(k1,k2) = kap1(k1);
        kap2msh(k1,k2) = kap2(k2);
    end
end
for i1 = 1:m1m
    for i2 = 1:m2mh
        if i1 <= m1mh
            r1msh(i1,i2) = (i1-1)*dx;
        else
            r1msh(i1,i2) = -(m1m-i1+1)*dx;
        end
        r2msh(i1,i2) = (i2-1)*dy;
    end
end
% Correlation via FFT2
cor = zeros(m1m,m2mh);
for i1 = 1:m1m
    for i2 = 1:m2mh
        r1 = r1msh(i1,i2);
        r2 = r2msh(i1,i2);
        cor(i1,i2) = sum(fhat2.*exp(1i*(kap1msh*r1+kap2msh*r2)),'all');
%         cor(i1,i2) = sum(2*fhat2(1:m1m,2:m2mh-1).*cos(kap1msh(1:m1m,2:m2mh-1)*r1+kap2msh(1:m1m,2:m2mh-1)*r2),'all') + ...
%                      sum(2*fhat2(2:m1mh-1,1:m2mh-1:m2mh).*cos(kap1msh(2:m1mh-1,1:m2mh-1:m2mh)*r1+kap2msh(2:m1mh-1,1:m2mh-1:m2mh)*r2),'all') + ...
%                      sum(fhat2(1:m1mh-1:m1mh,1:m2mh-1:m2mh).*cos(kap1msh(1:m1mh-1:m1mh,1:m2mh-1:m2mh)*r1+kap2msh(1:m1mh-1:m1mh,1:m2mh-1:m2mh)*r2),'all');
        cor(i1,i2) = cor(i1,i2) - fhat2(1,1);
        cor(i1,i2) = cor(i1,i2)/fvar;
    end
end
cor1 = real(cor);
figure
contourf(circshift(r1msh,[m1m/2-1,0]),...
         circshift(r2msh,[m1m/2-1,0]),...
         circshift(cor1,[m1m/2-1,0])) % positive/negative real parts

% Correlation
cor = zeros(m1m,m2mh);
for i1 = 1:m1m
    for i2 = 1:m2mh
        r1 = r1msh(i1,i2);
        r2 = r2msh(i1,i2);
        tmp = circshift(fluc,[i1-1 i2-1]);
        cor(i1,i2) = sum(fluc.*tmp,'all')/sum(fluc.*fluc,'all');
    end
end
cor2 = cor;
figure
contourf(circshift(r1msh,[m1m/2-1,0]),...
         circshift(r2msh,[m1m/2-1,0]),...
         circshift(cor2,[m1m/2-1,0])) % positive/negative real parts

% filtering
figure
contourf(circshift(kap1msh,[m1m/2-1,m2m/2-1]),... % circshift for display
     circshift(kap2msh,[m1m/2-1,m2m/2-1]),...
     circshift(fhat2,[m1m/2-1,m2m/2-1]))
% user-specified wavenumber for filtering 
fhatf = zeros(m1m,m2m);
fhatf(2,1) = fhat(2,1);  % first pair
fhatf(m1m,1) = fhat(m1m,1);
fhatf(1,3) = fhat(1,3);  % second pair
fhatf(1,m2m-1) = fhat(1,m2m-1);
fhatf = fhatf*(m1m*m2m);
ff = ifft2(fhatf);
figure
contourf(x,y,ff)

% one-dimensional spectra
fhat2_1d1 = sum(fhat2,2);
fhat2_1d2 = mean(abs(fft(f)/m1m).^2,2);
