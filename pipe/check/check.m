clear
close all
clc

m1m = 256;
m2m = 55;
m3m = 256;
m1mh = m1m/2 + 1;
m3mh = m3m/2 + 1;
m1 = m1m + 1;
fid = fopen('planetz_03201.bin','rb');
tmp = fread(fid,5*4*m1m*m3m,'double');
fclose(fid);
tmp = reshape(tmp,[5,4,m1m,m3m]);
f(1:m1m,1:m3m) = tmp(3,2,1:m1m,1:m3m);
fhat = fft2(f)/(m1m*m3m);
fhat2_1 = abs(fhat).^2;
fhat2_1(1,1) = 0;

% % test started
% f1 = zeros(m1,m3m);
% f1(1:m1m,1:m3m) = f(1:m1m,1:m3m);
% f1(m1,1:m3m) = f(1,1:m3m);
% fid = fopen('planetz_2_03201.q','wb');
% fwrite(fid, [16,m1,1,m3m,1,16], 'int32');
% fwrite(fid, 8*m1*m3m, 'int32');
% fwrite(fid,f1,'double');
% fwrite(fid, 8*m1*m3m, 'int32');
% % test ended

fid = fopen('spec2d2_3_3201.bin','rb');
fread(fid,1,'int32');
fhat2_2 = fread(fid,[m1mh,m3mh],'double');
fread(fid,1,'int32');
fclose(fid);

figure
levels = (0:0.1:1)*1e-04;
contourf(fhat2_1(1:m1mh,1:m3mh),levels)
xlim([0,20])
ylim([0,20])
figure
contourf(fhat2_2(1:m1mh,1:m3mh),levels)
xlim([0,20])
ylim([0,20])

% 1D spectra
fhat2_1d = sum(fhat2_1,2);
fhat2_1d(2:m1mh-1) = 2.0*fhat2_1d(2:m1mh-1);
fhat2_1d_1 = fhat2_1d(1:m1mh);

fid = fopen('specth_3201.bin','rb');
fread(fid,1,'int32');
tmp = fread(fid,m1mh*m2m*5,'double');
fread(fid,1,'int32');
fclose(fid);
tmp = reshape(tmp,[m1mh,m2m,5]);
fhat2_1d = tmp(1:m1mh,m2m-16,3);
fhat2_1d(2:m1mh-1) = 2.0*fhat2_1d(2:m1mh-1);
fhat2_1d_2 = fhat2_1d;

figure
hold on
plot(fhat2_1d_1(end-20:end),'.')
plot(fhat2_1d_2(end-20:end),'r')