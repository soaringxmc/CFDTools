clear
close all
clc
% should use ifft  
% user-specified location
ipl = 2;

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
fclose(fid);
%
fid = fopen('flowprop.dat', 'r');
line = fgetl(fid);
values = sscanf(line, '%f');
retau = values(2);
retau = 1.0;
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

fid = fopen('specze_ave.bin', 'rb');
numBytes = fread(fid, 1, 'int32');
if numBytes ~= (m3mh*m2m*5)*8
    disp('Inconsistency')
    return
end
specz_ave = fread(fid, m3mh*m2m*5, 'double');
fclose(fid);
specz_ave = reshape(specz_ave,[m3mh,m2m,5]);
lm3 = zeros(m3m,1);
for i = 1:m3m
    lm3(i) = laxis/(i-1.d0+eps)*retau;
end
specz_ave(2:m3mh-1,:,:) = 2.0*specz_ave(2:m3mh-1,:,:);

% auto-correlation
dz = laxis/m3m;
dist = transpose(dz*(0:m3mh-1));
R = zeros(m3mh,5);
for ivar = 1:5
    for i = 1:m3mh
        s = dz*(i-1);
        for k = 1:m3mh-1
            ome = 2*pi*k/laxis;
            R(i,ivar) = R(i,ivar) + specz_ave(k+1,jloc,ivar)*cos(ome*s);
        end
    end
    R(:,ivar) = R(:,ivar)/R(1,ivar);
end

figure
hold on
plot(dist,R(:,3))

data = [dist,R(:,1:5)];
fid = fopen('correlation.dat', 'w');
for i = 1:size(data, 1)
    for j = 1:size(data, 2)
        fprintf(fid, '%20.10e', data(i, j));
    end
    fprintf(fid, '\n');
end
fclose(fid);

