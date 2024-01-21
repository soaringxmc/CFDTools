close all
clear all
clc

% read data
fid = fopen('drag.dat', 'r');
fgetl(fid);
values = zeros(6,6,5);
for n = 1:5
    fgetl(fid);
    for k = 1:6 % six rows
        line = fgetl(fid);
        values(k,:,n) = sscanf(line, '%f');
    end
end
fclose(fid);

% obtain c(n,reb) based on ubp, ret
cn = zeros(6,5); % 6 n, 5 reb
kap = 0.387;
for n = 1:5
    for k = 1:6
        this_reb = values(k,2,n);
        this_ret = values(k,3,n);
        this_cf  = values(k,4,n)/4;
        this_rot = values(k,1,n);
        this_ubp = sqrt(2/this_cf);
        this_cn  = this_ubp - 1/kap*log(this_ret);
        cn(k,n) = this_cn;
    end
end

% cn vs rot
rot = values(:,1,:);
reb = values(:,2,:);
ret = values(:,3,:);
cf = values(:,4,:)/4;
dr = values(:,5,:);
figure(1)
for n = 1:5
    plot(rot(:,n),cn(:,n)-cn(1,n),'s-')
    hold on
end
% c(0)=1.229 unset
f = fit(rot(1:5,5),cn(1:5,5)-cn(1,5),'poly1');
plot(rot(:,5),f.p1*rot(:,5)+f.p2,'--')
alph = f.p1;
b = f.p2;

% cn vs rotp
figure(2)
rotp = rot/2.*reb./ret;
for n = 1:5
    semilogx(rotp(:,n),cn(:,n)-cn(1,n),'s-')
    hold on
end
f = fit(log(rotp(2:3,5)),cn(2:3,5)-cn(1,5),'poly1');
plot(rotp(:,5),f.p1*log(rotp(:,5))+f.p2,'--')
alphp = f.p1;
bp = f.p2;

% dr vs rotp
figure(3)
for n = 1:5
    semilogx(rotp(:,n),dr(:,n),'s-')
    hold on
end
dr_fit = zeros(6,1,5);
for k = 2:6
    dr_fit(k,5) = 2/(sqrt(2/cf(1,5))+1/kap)*(alphp*log(rotp(k,5))+bp);
end
% semilogx(rotp(2:6,5),dr_fit(2:6,5),'-s')
for k = 2:6
    dr_fit(k,5) = 2/(sqrt(2/cf(1,5))+1/kap)*(cn(k,5)-cn(1,5));
end
% semilogx(rotp(2:6,5),dr_fit(2:6,5),'-g')
for k = 1:6
    dr_fit(k,5) = 2/(sqrt(2/cf(1,5))+1/kap)*(alph*rot(k,5)+0);
end
% semilogx(rotp(1:6,5),dr_fit(1:6,5),'-b')

figure(4)
for n = 1:5
    plot(rot(:,n),dr(:,n),'s-')
    hold on
end
% dr vs rot
for k = 1:6
    dr_fit(k,5) = 2/(sqrt(2/cf(1,5))+1/kap)*(alph*rot(k,5)+b);
end
% plot(rot(1:6,5),dr_fit(1:6,5),'-b')

% cf prediction based on a linear change of c(n)-c(0) with
alph_mean = zeros(5);
cf_pred = zeros(6,1,5);
rotp_pred = zeros(6,1,5);
dr_pred = zeros(6,1,5);
for n = 1:5
    alph_mean(n) = mean((cn(2:5,n)-cn(1,n))./rot(2:5,n));
end
syms y
for n = 1:5
    for k = 1:6
        this_cf0 = cf(1,n);
        this_alpha = alph_mean(n);
        this_rot = rot(k,n);
        eqn = sqrt(2/y) - sqrt(2/this_cf0) - 1/kap*log(sqrt(y/this_cf0)) - this_alpha*this_rot == 0;
        cf_pred(k,n) = vpasolve(eqn,y);
    end
end
for n = 1:5
    for k = 1:6
        dr_pred(k,n) = 1 - cf_pred(k,n)/cf_pred(1,n);
    end
end
figure(4)
plot(rot(1:6,5),dr_pred(1:6,5),'--k','LineWidth',2)
xlabel('N')
ylabel('DR')

figure(3)
rotp_pred(1:6,5) = rot(1:6,5).*sqrt(2./cf_pred(1:6,5));
semilogx(rotp_pred(1:6,5),dr_pred(1:6,5),'--k','LineWidth',2)
xlabel('N+')
ylabel('DR')

% does not make sense, dr is far from small
% the discrepancy is large even using true cn without fitting
% only N=0.25 makes sense

cn(2:5,5)-cn(1,5)
