clear
close all
clc

fid = fopen('flowprop.dat', 'r');
line = fgetl(fid);
values = sscanf(line, '%f');
reb = values(1);
re = 10*floor(reb/10);
ub = 0.5*reb/re;
retau = values(2);
cf0 = values(3)/4;
fclose(fid);

cf_1 = 16./reb;
wvarav = importdata('stats_out.dat');
m2m = size(wvarav,1);
m2 = m2m + 1;
uruz = [wvarav(:,9);0.]/(ub*ub);
rc = [wvarav(:,1);1.];
yc = 1.-rc;
cf_2_cum = -8*cumtrapz(rc(m2:-1:1),rc(m2:-1:1).^2.*uruz(m2:-1:1));
cf_2 = cf_2_cum(m2);
cf = cf_1 + cf_2;
figure(1)
plot(yc,uruz,'b')
hold on
plot(yc,0.5*rc.^2.*uruz,'r')
plot(yc(m2:-1:1),cf_2_cum,'g')
hold off

% write to files
ycp = yc*retau;
file_id = fopen('FIK.dat','w');
for i = 1:m2
    nbytes = fprintf(file_id,'%15.10e %15.10e %15.10e %15.10e %15.10e\n',...
                    [yc(i),ycp(i),uruz(i),0.5*rc(i)*uruz(i),cf_2_cum(m2-i+1)]);
end
fclose(file_id);

fid = fopen('error.dat','w');
fprintf(fid, '%f %f %f %f', [cf0,cf_1,cf_2,(cf_1+cf_2-cf0)/cf0]);
fclose(fid);



