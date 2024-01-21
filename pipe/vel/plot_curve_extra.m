figure(2)
hold on
load('REB10000_N2.0_RET242.58_KMAX51_NUT.mat','r','nut')
plot(r,nut,'k')
load('REB10000_N2.0_RET226.56_KMAX51_NUT.mat','r','nut')
plot(r,nut,'b')
load('REB10000_N2.0_RET212.11_KMAX51_NUT.mat','r','nut')
plot(r,nut,'g')
legend('REB10000 N2.0 RET242.58 KMAX51 2 NUT',...
       'REB10000 N2.0 RET226.56 KMAX51 1.75 NUT',...
       'REB10000 N2.0 RET212.11 KMAX51 1.5 NUT')
hold off

figure(3)
hold on
load('REB10000_N1.0_RET283.20_KMAX51_RI.mat','r','ri')
plot(r,ri,'c')
load('REB10000_N2.0_RET242.58_KMAX51_RI.mat','r','ri')
plot(r,ri,'k')
load('REB10000_N3.0_RET208.98_KMAX51_RI.mat','r','ri')
plot(r,ri,'r')
legend('REB10000 N1.0 RET283.20 KMAX51 RI',...
       'REB10000 N2.0 RET242.58 KMAX51 RI',...
       'REB10000 N3.0 RET208.98 KMAX51 RI')
hold off
return

figure(4)
hold on
load('REB10000_N2.0_RET242.58_KMAX51_STRAIN.mat','r','strain')
plot(r,strain,'k')
load('REB10000_N2.0_RET226.56_KMAX51_STRAIN.mat','r','strain')
plot(r,strain,'b')
load('REB10000_N2.0_RET212.11_KMAX51_STRAIN.mat','r','strain')
plot(r,strain,'g')
legend('REB10000 N2.0 RET242.58 KMAX51 2 STRAIN',...
       'REB10000 N2.0 RET226.56 KMAX51 1.75 STRAIN',...
       'REB10000 N2.0 RET212.11 KMAX51 1.5 STRAIN')
hold off

figure(5)
hold on
load('REB10000_N2.0_RET242.58_KMAX51_UTH.mat','r','uth')
plot(r,uth,'k')
load('REB10000_N2.0_RET226.56_KMAX51_UTH.mat','r','uth')
plot(r,uth,'b')
load('REB10000_N2.0_RET212.11_KMAX51_UTH.mat','r','uth')
plot(r,uth,'g')
legend('REB10000 N2.0 RET242.58 KMAX51 2 UTH',...
       'REB10000 N2.0 RET226.56 KMAX51 1.75 UTH',...
       'REB10000 N2.0 RET212.11 KMAX51 1.5 UTH')
hold off