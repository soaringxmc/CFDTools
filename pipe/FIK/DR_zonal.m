clear
close all
clc

fid = fopen('FIK_zonal.dat', 'r');
line = fgetl(fid);
values = sscanf(line, '%f');
fclose(fid);

fid = fopen('FIK_zonal0.dat', 'r');
line = fgetl(fid);
values0 = sscanf(line, '%f');
fclose(fid);

% [n, reb, cf_1,cf_2i,cf_2o,cf_2i+cf_2o,cf_1+cf_2i+cf_2o]
fid = fopen('DR_zonal.dat','w');
fprintf(fid, '%15.10e %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e', ...
    [values(1:3);(values0(4:8)-values(4:8))/values0(8)]);
fclose(fid);

fid = fopen('FIK_zonal_uni.dat','w'); % divided by Cf_A0
fprintf(fid, '%15.10e %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e', ...
    [values(1:3);values(4:8)/values0(8)]);
fclose(fid);

fid = fopen('FIK_zonal_uni_cf.dat','w'); % divided by Cf_A>=0
fprintf(fid, '%15.10e %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e', ...
    [values(1:3);values(4:8)/values(8)]);
fclose(fid);
