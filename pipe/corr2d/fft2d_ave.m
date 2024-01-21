clear
close all
clc

ks = 22960;
ke = 22960;
kinc = 1;
k_tot = floor((ke-ks)/kinc) + 1;
prefix = 'planetz_2';
postfix = '_kap118-2304';

fid = fopen('bou.in', 'r');
fgetl(fid);
line = fgetl(fid);
values = sscanf(line, '%d');
m1 = values(1);
m3 = values(3);
m1m=m1-1;
m3m=m3-1;
fclose(fid);

% average spectra
fhat2_ave = zeros(m1m,m3m);
for k = ks:kinc:ke
    fid = fopen([prefix,'_',sprintf('%05d', k),postfix,'_fhat2.bin'], 'rb');
    tmp = fread(fid, [m1m,m3m], 'double');
    fhat2_ave = fhat2_ave + tmp;
    fclose(fid);
end
fhat2_ave = fhat2_ave/k_tot;
fid = fopen([prefix,postfix,'_fhat2.bin'], 'wb');
fwrite(fid, fhat2_ave, 'double');
fclose(fid);

