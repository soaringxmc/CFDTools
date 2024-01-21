clear
close all
clc

prefix = 'planetz_2_uz';
postfix = '_kap24-36';

fid = fopen('bou.in', 'r');
fgetl(fid);
line = fgetl(fid);
values = sscanf(line, '%d');
m1 = values(1);
m3 = values(3);
m1m=m1-1;
m3m=m3-1;
m1mh=m1m/2+1;
m3mh=m3m/2+1;
for i = 1:3
    fgetl(fid);
end
line = fgetl(fid);
values = sscanf(line, '%d');
laxis = values(1);
ltheta = 2*pi;
dtheta = ltheta/m1m;
dz = laxis/m3m;
fclose(fid);
% wavenumber
kap1 = zeros(m1m,1);
kap3 = zeros(m3m,1);
for k = 1:m1m
    if k <= m1mh
        kap1(k) = (k-1)*2*pi/ltheta;
    else
        kap1(k) = -(m1m-k+1)*2*pi/ltheta;
    end
end
for k = 1:m3m
    if k <= m3mh
        kap3(k) = (k-1)*2*pi/laxis;
    else
        kap3(k) = -(m3m-k+1)*2*pi/laxis;
    end
end
kap1msh = zeros(m1m,m3m);
kap3msh = zeros(m1m,m3m);
r1msh = zeros(m1m,m3mh);
r3msh = zeros(m1m,m3mh);
for k1 = 1:m1m
    for k3 = 1:m3m
        kap1msh(k1,k3) = kap1(k1);
        kap3msh(k1,k3) = kap3(k3);
    end
end
for i1 = 1:m1m
    for i3 = 1:m3mh
        if i1 <= m1mh
            r1msh(i1,i3) = (i1-1)*dtheta;
        else
            r1msh(i1,i3) = -(m1m-i1+1)*dtheta;
        end
        r3msh(i1,i3) = (i3-1)*dz;
    end
end
% read spectra
fid = fopen([prefix,postfix,'_fhat2.bin'], 'rb');
fhat2_ave = fread(fid, [m1m,m3m], 'double');
fclose(fid);
% correlation
cor = ifft2(fhat2_ave*(m1m*m3m));
cor = cor - fhat2_ave(1,1);
cor = cor(:,1:m3mh);

fid = fopen([prefix,postfix,'_corr.bin'], 'wb');
fwrite(fid, cor, 'double');
fclose(fid);

% ridge line
cor_coef = cor/cor(1,1);
i1max = zeros(1,m3mh);
i1max(1) = 1;
for i3 = 2:m3mh
    cormax0 = -1;
    i1l = i1max(i3-1)-5;
    i1r = min(i1max(i3-1)+5,m1m);
    i1range = (i1l:i1r);
    if i1l < 1 
        i1range = [(m1m+i1l:m1m),(1:i1r)];
    end
    for i1 = i1range
        if cor_coef(i1,i3) > cormax0
            cormax0 = cor_coef(i1,i3);
            i1max(i3) = i1;
        end
    end
end
data = [r1msh(i1max(1:m3mh),1),transpose(r3msh(1,1:m3mh))];
fid = fopen([prefix,postfix,'_corrmax.dat'], 'w');
for i = 1:size(data, 1)
    for j = 1:size(data, 2)
        fprintf(fid, '%20.10e', data(i, j));
    end
    fprintf(fid, '\n');
end
fclose(fid);
