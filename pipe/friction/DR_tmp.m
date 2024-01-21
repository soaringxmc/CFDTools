clear all
clc
rot = [0.24999878547154220 0.49999753726463692 0.99999590294757779 1.9999903519028717 3.9999810883361993];
reb = [133000.64613227852 133000.65509083320 133000.54491020468 133000.64160155415 133000.62881579436];
ret = [2820.4250909185075 2586.2950800420476 2275.2027517940624 1882.8919620283230 1699.3876327041698];
dr = [0.13435664910289047 0.27210967447355949 0.43668576237252776 0.61420136130795622 0.68573582329880201];
rot_plus = rot/2.*reb./ret;
logrot = log(rot_plus);

%
kap = 0.387;
ubp = reb./(2*ret);
term1 = 1/kap*log(ret);
cn = ubp - term1;
f=fit([0,rot(1:4)]',[0,cn(1:4)-1.229]','poly1')
% plot([0,rot(1:4)],[0,cn(1:4)-1.229],'*');  % alpha = 8 3020
semilogx(rot_plus(1:4),cn(1:4)-1.229,'*');


rot = [0.24999878547154220 0.49999753726463692 0.99999590294757779 1.9999903519028717 3.9999810883361993];
reb = [82500 82500 82500 82500 82500];
ret = [1861.11 1723.57 1498.29 1316.85 1221.82];
rot_plus = rot/2.*reb./ret;
logrot = log(rot_plus);


%
kap = 0.387;
ubp = reb./(2*ret);
term1 = 1/kap*log(ret);
cn = ubp - term1;
figure(1)
hold on
f=fit([0,rot(1:4)]',[0,cn(1:4)-1.229]','poly1')
% plot([0,rot(1:4)],[0,cn(1:4)-1.229],'*');  % alpha = 7.40 2000
semilogx(rot_plus(1:4),cn(1:4)-1.229,'*');

rot = [0.24999878547154220 0.49999753726463692 0.99999590294757779 1.9999903519028717 3.9999810883361993];
reb = [44000 44000 44000 44000 44000];
ret = [1087.91 1013.63 892.252 808.564 755.777];
rot_plus = rot/2.*reb./ret;
logrot = log(rot_plus);
%
kap = 0.387;
ubp = reb./(2*ret);
term1 = 1/kap*log(ret);
cn = ubp - term1;
figure(1)
hold on
f=fit([0,rot(1:4)]',[0,cn(1:4)-1.229]','poly1')
% plot([0,rot(1:4)],[0,cn(1:4)-1.229],'*');  % alpha = 5.87 1000
semilogx(rot_plus(1:4),cn(1:4)-1.229,'*');

% semilogx([1,10],[1,20],'-')

% % f=fit(logrot',dr','poly1');
% % DR = (0.175ln(N+) - 0.167)*100%
% f=fit(logrot([1,2,3,4])',dr([1,2,3,4])','poly1');
% DR_fit = (0.1941*logrot - 0.2159);
% plot(f,logrot([1,2,3,4])',dr([1,2,3,4])');
% 
% ret0 = 3020.16;
% rot_star = rot/2.*reb/ret0;
% rot_star = rot_plus.*(1-DR_fit).^0.5; %
% 
% p1 = 0.1941;
% p2 = -0.2159;
% term_l = DR_fit + 0.5*p1*log(1-DR_fit);
% term_r = p1*log(rot_star)+p2;
% 
% syms DR0
% k = 1;
% DR_fit = [];
% rot_star = linspace(2,99,40);
% rot_star = linspace(0.1,4.0,40)*133000/(2*3020.16);
% for rot_star0 = rot_star
%     eqn = DR0 + 0.5*p1*log(1-DR0)-p1*log(rot_star0)-p2 == 0;
%     DR_fit0 = vpasolve(eqn,DR0);
%     DR_fit = [DR_fit, DR_fit0]; 
% end
% figure(1)
% semilogx(rot_star,DR_fit*100)
% fprintf('%.1f\n', (linspace(0.1,4.0,40))')
% fprintf('%.4f\n', DR_fit')
% hold on
% rot_star = rot/2.*reb/ret0;
% semilogx(rot_star,dr*100)
% % rot_star = 1:5:205;
% % DR_fit = 1.0/(1-0.5*p1)*(p1*log(rot_star) + p2);
% % semilogx(rot_star,DR_fit*100)


