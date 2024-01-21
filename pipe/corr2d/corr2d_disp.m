clear
close all
clc
%
prefix = 'planetz_2_kap118-2304';
%
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
r1msh = zeros(m1m,m3mh);
r3msh = zeros(m1m,m3mh);
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
%
fid = fopen([prefix,'_corr.bin'], 'rb');
cor = fread(fid, [m1m m3mh], 'double');
fclose(fid);
cor_coef = cor/cor(1,1);
figure
hold on
contourf(circshift(r1msh,[m1m/2-1,0]),... % circshift for display
         circshift(r3msh,[m1m/2-1,0]),...
         circshift(cor_coef,[m1m/2-1,0]))
xlim([-0.3 0.3]);
ylim([0 1.5]);
%
data = importdata([prefix,'_corrmax.dat']);
theta_max = data(:, 1);
zeta_max = data(:, 2);
scatter(theta_max,zeta_max,'.w')
