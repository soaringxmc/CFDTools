clear
close all
clc

ks = 22960;
ke = 22960;
kinc = 20;
k_tot = floor((ke-ks)/kinc) + 1;
prefix = 'planetz_2_uz';
kapfl = 24;
kapfr = 36;

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
r1msh = zeros(m1m,m3m);
r3msh = zeros(m1m,m3m);
for k1 = 1:m1m
    for k3 = 1:m3m
        kap1msh(k1,k3) = kap1(k1);
        kap3msh(k1,k3) = kap3(k3);
    end
end
for i1 = 1:m1m
    for i3 = 1:m3m
        r1msh(i1,i3) = (i1-1)*dtheta;
        r3msh(i1,i3) = (i3-1)*dz;
    end
end

% read flowfield
for k = ks:kinc:ke
%   fid = fopen([prefix,'_',sprintf('%05d', k),'.q'], 'rb');
    fid = fopen([prefix,'.q'], 'rb');
    tmp = fread(fid, 6, 'int32');
    numBytes = fread(fid, 1, 'int32');
    if numBytes ~= (m1*m3m)*8
        disp('Inconsistency')
        return
    end
    tmp = fread(fid, [m1,m3m], 'double');
    f = tmp(1:m1m,1:m3m);
    fluc = f-mean(f,'all');
    fclose(fid);
    tmp = fft2(f);
    fhat = tmp/(m1m*m3m);
    fhat2 = abs(fhat).^2;
    
    % user-specified wavenumber for filtering
    fhatf = zeros(m1m,m3m);
    ifl1 = kapfl + 1; % one pair
    ifr1 = kapfr + 1;
    ifl2 = m1m-kapfl+1;
    ifr2 = m1m-kapfr+1;
    if kapfr == m1m/2
        ifr2 = m1m-(kapfr-1)+1;
    end
    fhatf(ifl1:ifr1,1:m3m) = fhat(ifl1:ifr1,1:m3m);
    fhatf(ifl2:-1:ifr2,1:m3m) = fhat(ifl2:-1:ifr2,1:m3m);
    fhatf = fhatf*(m1m*m3m);
    ff = ifft2(fhatf);
    % figure
    % contourf(r1msh,r3msh,ff)
    data = zeros(m1,m3m);
    data(1:m1m,:) = ff;
    data(m1,:) = ff(1,:);
    if kapfl == kapfr
        qname = [prefix,'_kap',int2str(kapfl),'.q'];
%         qname = [prefix,'_',sprintf('%05d', k),'_kap',int2str(kapfl),'.q'];
    else
        qname = [prefix,'_kap',int2str(kapfl),'-',int2str(kapfr),'.q'];
%         qname = [prefix,'_',sprintf('%05d', k),'_kap',int2str(kapfl),'-',int2str(kapfr),'.q'];
    end
    fid = fopen(qname, 'wb');
    fwrite(fid, [16,m1,1,m3m,1,16], 'int32');
    fwrite(fid, 8*m1*m3m, 'int32');
    fwrite(fid, data, 'double');
    fwrite(fid, 8*m1*m3m, 'int32');
    fclose(fid);
end
