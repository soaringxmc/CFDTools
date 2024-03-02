clear
close all
clc

syms y

fid = fopen('drag_moody.dat', 'r');
fgetl(fid);
values = zeros(5,6,6);
for n = 1:6
    fgetl(fid);
    for k = 1:5
        line = fgetl(fid);
        values(k,:,n) = sscanf(line, '%f');
    end
end
fclose(fid);

re = values(:,2,:);
rot = values(:,1,:);
f = values(:,4,:);
re1 = logspace(3,6,20);
f1 = zeros(length(re1),6);
for n = 1:6
    for k = 1:length(re1)
        re0 = re1(k);
        rot0 = rot(1,n);
        eqn = (-0.27*rot0^2+2.55*rot0+2.10)*log10(re0*y^0.5)-(-0.70*rot0^2+7.18*rot0+1.15)-1/y^0.5 == 0;
        f1(k,n) = vpasolve(eqn,y);
    end
end

dr = f;
ps = f;
dr1 = f1;
ps1 = f1;
for n = 1:6
    dr(:,n) = 1 - f(:,n)./f(:,1);
    ps(:,n) = dr(:,n) - 16*rot(:,n).^2./(re(:,n).*f(:,1));
    dr1(:,n) = 1 - f1(:,n)./f1(:,1);
    ps1(:,n) = dr1(:,n) - 16*rot(1,n)^2./(re1(:).*f1(:,1));
end

figure
lncolor = ['k','r','g','b','c','m'];
for n = 1:6
    plot(re(:,n),f(:,n),['--',lncolor(n)])
    hold on
    plot(re1(:),f1(:,n),['-',lncolor(n)])
end

xlim([5e03 1.4e05])
ylim([0.006 0.04])
xlabel('Reb')
ylabel('f')
hold off

figure
for n = 1:6
    plot(re(:,n),dr(:,n),['--',lncolor(n)])
    hold on
    plot(re1(:),dr1(:,n),['-',lncolor(n)])
end
xlim([5e03 1.4e05])
xlabel('Reb')
ylabel('DR')
hold off

figure
for n = 1:6
    plot(re(:,n),ps(:,n),['--',lncolor(n)])
    hold on
    plot(re1(:),ps1(:,n),['-',lncolor(n)])
end
xlim([5e03 1.4e05])
xlabel('Reb')
ylabel('PS')
hold off