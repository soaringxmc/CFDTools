clear
close all
clc

% fid = fopen('flowprop0.dat', 'r');
% line = fgetl(fid);
% values = sscanf(line, '%f');
% retau0 = values(2);
% fclose(fid);

% fid = fopen('bou.in', 'r');
% for k = 1:9
%     fgetl(fid);
% end
% line = fgetl(fid);
% values = sscanf(line, '%f');
% vt = values(2);
% 
% fid = fopen('flowprop.dat', 'r');
% line = fgetl(fid);
% values = sscanf(line, '%f');
% reb = values(1);
% re = 10*floor(reb/10);
% retau = values(2);
% utau = values(4);
% fclose(fid);

% test case, wei et al.
% influence of grid point number and distribution
reb = 50000;
re = reb;
vt = 0.25;
rotp0 = 120;
lrotp0 = 1;

retau0 = (0.3164/32)^0.5*reb^(7/8); % blasius law

wvarav = importdata('stats_out.dat');
m2m = size(wvarav,1);
m2 = m2m + 1;
uz = wvarav(:,5);
rm = wvarav(:,2);
ym = 1.0-rm;

figure
hold on
% plot(ym,uz,'.')

retau_l = 100;
retau_r = 10000;

for k = 1:1000
    % initialization
    retau = 0.5*(retau_l+retau_r);
    utau = retau/re;
    rm = [0;wvarav(:,2);1.0];
    ym = 1.0-rm;
    ymp = ym*retau;
    % blended law
    kap = 0.4;
    rotp = vt/utau;
    if lrotp0 == 1
        rotp = rotp0;
        vt = rotp*utau;
    end
    ret_rat = retau/retau0;
    y10 = 11;
    y1mp = y10*ret_rat;
    kap1 = kap*ret_rat/((rotp/retau)^2+1)*3.5;
    b = 1/3*sqrt(kap1/kap);
    uzp_pi = 1.0/kap1*log(1+kap1*ymp) + 7.8*ret_rat*(1-exp(-ymp/y1mp)-(ymp/y1mp).*exp(-b*ymp));
    A = sqrt(0.06052*rotp^2+5.4*rotp+25.705);
    B = 1.55+0.338*(1-exp(-0.0553*rotp));
    uclp = reb/(2*retau)+2*A/(2+B); % induce errors
    uclp = 58.284; % @tmp
    uzp_po = uclp-A*rm.^B;
    a1 = 5/ret_rat;
    uzp_p = uzp_pi.*exp(-a1*(1-rm)) + uzp_po.*(1-exp(-a1*(1-rm)));
    % mass rate
    ruzp_p = rm.*uzp_p;
    mass = trapz(rm,ruzp_p);
    x = mass*4*retau;
    diff = (x-reb)/reb;
    disp(diff)
    if abs(diff) < 1e-4
        break
    elseif diff > 0
        retau_r = retau; 
    else
        retau_l = retau;
    end
end

uz_p = uzp_p*utau;
uz_pi = uzp_pi*utau;
uz_po = uzp_po*utau;

uzb_p = uz_p*2*re/reb;
uzb_pi = uz_pi*2*re/reb;
uzb_po = uz_po*2*re/reb;
fname = ['REB',num2str(reb,'%i'),...
         '_Z',num2str(rotp,'%i'),...
         '_RET',num2str(retau,'%.2f'),...
         ];
save([fname,'.mat'],'rm','uzb_p','uzb_pi','uzb_po');

plot(ym,uz_p)
plot(ym,uz_pi)
plot(ym,uz_po)