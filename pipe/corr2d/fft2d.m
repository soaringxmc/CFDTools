clear
close all
clc

lutau2up = 0;
ks = 22960;
ke = 22960;
kinc = 1;
k_tot = floor((ke-ks)/kinc) + 1;
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
%
utau = 1.0;
if lutau2up == 1
    fid = fopen('flowprop.dat', 'r');
    line = fgetl(fid);
    values = sscanf(line, '%f');
    utau = values(4);
    fclose(fid);
end
% compute spectra
for k = ks:kinc:ke
%     fid = fopen([prefix,'_',sprintf('%05d', k),postfix,'.q'], 'rb');
    fid = fopen([prefix,postfix,'.q'], 'rb');
    tmp = fread(fid, 6, 'int32');
    numBytes = fread(fid, 1, 'int32');
    if numBytes ~= (m1*m3m)*8
        disp('Inconsistency')
        return
    end
    tmp = fread(fid, [m1,m3m], 'double');
    f = tmp(1:m1m,1:m3m)/utau;
    fclose(fid);
    tmp = fft2(f);
    fhat = tmp/(m1m*m3m);
    fhat2 = abs(fhat).^2;
%     fid = fopen([prefix,'_',sprintf('%05d', k),postfix,'_fhat2.bin'], 'wb');
    fid = fopen([prefix,postfix,'_fhat2.bin'], 'wb');
    fwrite(fid, fhat2, 'double');
    fclose(fid);
end
