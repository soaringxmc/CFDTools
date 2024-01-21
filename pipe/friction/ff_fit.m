clear
close all
clc

fid = fopen('drag_moody.dat', 'r');
fgetl(fid);
values = zeros(30,6);
for n = 1:6
    fgetl(fid);
    for k = 1:5
        line = fgetl(fid);
        values(k+(n-1)*5,:) = sscanf(line, '%f');
    end
end
fclose(fid);

ftype = fittype('(ca*rot^2+ce*rot+2.10)*log10(fsre)-(cb*rot^2+cf*rot+1.15)',...
               'independent',{'fsre','rot'},...
               'dependent',{'invfs'},...
               'coefficients',{'ca','cb','ce','cf'});

re = values(:,2);
f = values(:,4);
rot = values(:,1);
% isel = [2:5,7:10,12:15,17:20,22:25,27:30];
isel = 1:30;
re_sel = re(isel);
f_sel = f(isel);
rot_sel = rot(isel);
fsre = f_sel.^0.5.*re_sel;
invfs = 1./f_sel.^0.5;
cfit = fit([fsre,rot_sel],invfs,ftype);
ca = cfit.ca;
cb = cfit.cb;
% cc = cfit.cc;
% cd = cfit.cd;
ce = cfit.ce;
cf = cfit.cf;
figure
lncolor = ['k','r','g','b','c','m'];
for n = 1:6
    rek = re(1+5*(n-1):5*n);
    fk = f(1+5*(n-1):5*n);
    rot0 = rot(1+5*(n-1));
    loglog(rek,fk,['--',lncolor(n)])
    hold on
    func = @(re,f) (ca*rot0^2+ce*rot0+2.10).*log10(re.*f.^0.5)-(cb*rot0^2+cf*rot0+1.15)-1./f.^0.5;
    fimplicit(func,[5e03,1.4e05,0.004,0.04],['-',lncolor(n)])
end
hold off

disp([ca ce])
disp([cb cf])