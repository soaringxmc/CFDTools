clear
close all
clc

fid = fopen('flowprop.dat', 'r');
line = fgetl(fid);
values = sscanf(line, '%f');
reb = values(1);
re = 10*round(reb/10);
ub = 0.5*reb/re;
retau = values(2);
cf0 = values(3)/4;
fclose(fid);

fid = fopen('DR.dat','r');
line = fgetl(fid);
values = sscanf(line, '%f');
rot = values(1);
fclose(fid);

cf_1 = 16./reb;
wvarav = importdata('stats_out.dat'); % _interp.dat for RETAU1140/2000_A0
m2m = size(wvarav,1);
m2 = m2m + 1;
uruz = [wvarav(:,9);0.]/(ub*ub);
rc = [wvarav(:,1);1.];
yc = 1.-rc;

vec_sep = [30 68 133 209 290];
vec_re = [5300 17000 44000 82500 133000];
loc = vec_sep(vec_re==re);

cf_2i = -8*trapz(rc(m2:-1:loc),rc(m2:-1:loc).^2.*uruz(m2:-1:loc));
cf_2o = -8*trapz(rc(loc:-1:1),rc(loc:-1:1).^2.*uruz(loc:-1:1));

fid = fopen('FIK_zonal.dat','w');
fprintf(fid, '%15.10e %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e', ...
              [rot,reb,retau,cf_1,cf_2i,cf_2o,cf_2i+cf_2o,cf_1+cf_2i+cf_2o]);
fclose(fid);
