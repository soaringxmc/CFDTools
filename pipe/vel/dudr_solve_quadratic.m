clear
clc

syms y

n = 2.0;
reb = 10000;
beta = 1./6.;
alph = 2.0;
kap = 0.4;
% kmax depends on Reb
% large kmax needed for high Reb
kmax = 51;
retlim = [300 300]; % 200 400
diff_lim = 1e-03;

r = (0:kmax-1)/(kmax-1);
yr = 1.0-r;
dupdr = zeros(1,kmax); % up->u+
while 1
  ret = 0.5*(retlim(1)+retlim(2));
  zz = n*reb/(2.0*ret);
  for k = 1:kmax
    x = r(k);
    ri = 6*zz^2*x^2/(y^2+zz^2*x^2);
    da = (1-exp(-(1-x)*ret/26));
    l0 = kap*(0.35-0.2*x^2-0.15*x^4)*da;
    l = (1-beta*ri)^alph*l0;
    eqn = l^2*sqrt(y^2+zz^2*x^2)*y + x + 1/ret*y == 0;
    dupdr(k) = vpasolve(eqn,y);
  end
  % solve for up(r), integrate from r=1 to r=0
  % and calculate difference in Reb
  up(kmax:-1:1) = cumtrapz(r(kmax:-1:1),dupdr(kmax:-1:1));
  rup = up.*r;
  fr(1:kmax) = cumtrapz(r(1:kmax),rup(1:kmax));
  x = fr(kmax)*4*ret;
  diff = (x-reb)/reb;
  disp([ret,diff])
  if abs(diff) < diff_lim
    break
  else
    if diff < 0
      retlim = [ret, retlim(2)];
    else
      retlim = [retlim(1), ret];
    end
  end
end
u = up*2*ret/reb;
fname = ['REB',num2str(reb,'%i'),...        
         '_N',num2str(n,'%.1f'),...
         '_RET',num2str(ret,'%.2f'),...
         '.mat'];
save(fname,'r','u');

