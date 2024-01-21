clear 
close all
clc

w0 = importdata('stats_out0.dat'); % A0 case
m2m0 = size(w0,1);
yc0 = w0(:,1);
w = importdata('stats_out.dat');
m2m = size(w,1);
yc = w(:,1);
w(:,3:end) = 0.0;
w(:,9) = interp1(yc0, w0(:,9), yc, 'spline');

% writematrix(w, 'stats_out_interp.dat'); % A0 case
dlmwrite('stats_out_interp.dat', w, 'delimiter', '\t'); % A0 case
