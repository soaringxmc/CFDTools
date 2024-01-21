clear
close all
clc

% user-specified location
ipl = 3;
ivar = 3;

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
for i = 1:3
    fgetl(fid);
end
line = fgetl(fid);
values = sscanf(line, '%d');
laxis = values(1);
ltheta = 2.0*pi;
fclose(fid);
%
fid = fopen('flowprop.dat', 'r');
line = fgetl(fid);
values = sscanf(line, '%f');
retau = values(2);
fclose(fid);
%
data = importdata('radcor.out');
rm = data(:, 3);
%
jplanes(1) = m2m;
jplanes(2) = m2m-16;
for j = 1:m2m
    if rm(j) > 0.8
        break
    end
end
jplanes(3) = j;
for j = 1:m2m
    if rm(j) > 0.5
        break
    end
end
jplanes(4) = j;
for j = 1:m2m
    if rm(j) > 0.2
        break
    end
end
jplanes(5) = j;
jloc = jplanes(ipl);
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

% theta direction
fid = fopen('specth_ave.bin', 'rb');
numBytes = fread(fid, 1, 'int32');
if numBytes ~= (m1mh*m2m*5)*8
    disp('Inconsistency')
    return
end
specth_ave = fread(fid, m1mh*m2m*5, 'double');
fclose(fid);
specth_ave = reshape(specth_ave,[m1mh,m2m,5]);
specth_ave(2:m1mh-1,:,:) = 2.0*specth_ave(2:m1mh-1,:,:);
specth_ave_ener = 2.0*specth_ave;

figure('Name','theta');
loglog(kap1(1:m1mh),specth_ave_ener(1:m1mh,jloc,ivar))
data = [kap1(1:m1mh),specth_ave_ener(1:m1mh,jloc,ivar)];
fid = fopen(['specth_ave_ener_',int2str(ipl),'_',int2str(3),'.dat'], 'w');
for i = 1:size(data, 1)
    for j = 1:size(data, 2)
        fprintf(fid, '%20.10e', data(i, j));
    end
    fprintf(fid, '\n');
end
fclose(fid);

% axial direction
fid = fopen('specze_ave.bin', 'rb');
numBytes = fread(fid, 1, 'int32');
if numBytes ~= (m3mh*m2m*5)*8
    disp('Inconsistency')
    return
end
specz_ave = fread(fid, m3mh*m2m*5, 'double');
fclose(fid);
specz_ave = reshape(specz_ave,[m3mh,m2m,5]);
specz_ave(2:m3mh-1,:,:) = 2.0*specz_ave(2:m3mh-1,:,:);
specz_ave_ener = 2.0*specz_ave;

figure('Name','zeta');
loglog(kap3(1:m3mh),specz_ave_ener(1:m3mh,jloc,ivar))
data = [kap3(1:m3mh),specz_ave_ener(1:m3mh,jloc,ivar)];
fid = fopen(['specze_ave_ener_',int2str(ipl),'_',int2str(3),'.dat'], 'w');
for i = 1:size(data, 1)
    for j = 1:size(data, 2)
        fprintf(fid, '%20.10e', data(i, j));
    end
    fprintf(fid, '\n');
end
fclose(fid);