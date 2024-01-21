clear
close all
clc
%
ks = 22960;
ke = 22960;
kinc = 1;
k_tot = floor((ke-ks)/kinc) + 1;
prefix = 'continua';
kapfl = 3;
kapfr = 3;
%
fid = fopen('bou.in', 'r');
fgetl(fid);
line = fgetl(fid);
values = sscanf(line, '%d');
m1 = values(1);
m2 = values(2);
m3 = values(3);
m1m=m1-1;
m2m=m2-1;
m3m=m3-1;
m1mh=m1m/2+1;
m3mh=m3m/2+1;
% read flowfield
for k = ks:kinc:ke
    fid = fopen([prefix,'.q'], 'rb');
    fread(fid, 6, 'int32');
    numBytes = fread(fid, 1, 'int32');
    if numBytes ~= (m1*m2m*m3m)*8
        disp('Inconsistency')
    end
    tmp = fread(fid, m1*m2m*m3m, 'double');
    tmp = reshape(tmp,[m1,m2m,m3m]);
    f = zeros(m1m,m3m,m2m);
    ff = zeros(m1m,m3m,m2m);
    fhat = zeros(m1m,m3m,m2m);
    for i2 = 1:m2m
        f(1:m1m,1:m3m,i2) = tmp(1:m1m,i2,1:m1m);
    end
    for i2 = 1:m2m
        tmp = f(1:m1m,1:m3m,i2);
        tmp = fft2(tmp);
        tmp = tmp/(m1m*m3m);
        fhat(1:m1m,1:m3m,i2) = tmp;
    end
    % user-specified wavenumber for filtering
    fhatf = zeros(m1m,m3m,m2m);
    ifl1 = kapfl + 1; % one pair
    ifr1 = kapfr + 1;
    ifl2 = m1m-kapfl+1;
    ifr2 = m1m-kapfr+1;
    if kapfr == m1m/2
        ifr2 = m1m-(kapfr-1)+1;
    end
    if kapfl == 0
        ifl2 = m1m;
    end
    fhatf(ifl1:ifr1,1:m3m,1:m2m) = fhat(ifl1:ifr1,1:m3m,1:m2m);
    fhatf(ifl2:-1:ifr2,1:m3m,1:m2m) = fhat(ifl2:-1:ifr2,1:m3m,1:m2m);
    for i2 = 1:m2m
        tmp = fhatf(1:m1m,1:m3m,i2);
        tmp = tmp*(m1m*m3m);
        tmp = ifft2(tmp);
        ff(1:m1m,1:m3m,i2) = tmp;
    end
%
    data = zeros(m1,m2m,m3m);
    for i2 = 1:m2m
        data(1:m1m,i2,1:m3m) = ff(1:m1m,1:m3m,i2);
    end
    data(m1,1:m2m,1:m3m) = data(1,1:m2m,1:m3m);
    if kapfl == kapfr
        qname = [prefix,'_kap',int2str(kapfl),'.q'];
    else
        qname = [prefix,'_kap',int2str(kapfl),'-',int2str(kapfr),'.q'];
    end
    fid = fopen(qname, 'wb');
    fwrite(fid, [16,m1,m2m,m3m,1,16], 'int32');
    fwrite(fid, (m1*m2m*m3m)*8, 'int32');
    fwrite(fid, data, 'double');
    fwrite(fid, (m1*m2m*m3m)*8, 'int32');
    fclose(fid);
end