close all
clear
clc

syms y

n = 2.0;
reb = 50000;
beta = 1/6;
alph = 2;
kap = 0.4;
lda = 1;
lvis = 0;
locvis = 11.6;
lblend = 0;
locpblend = 10;
leqn = 1;
lusez = 1;
zz = 120;
% larg kmax desired at high Reb
kmax = 201; % (reb,kmax)=(10^4,51),(10^5,201)
retlim = [500 2000];
diff_lim = 1e-04;

% r = (0:kmax-1)/(kmax-1); % linear spacing
ret0 = (0.3164/32)^0.5*reb^(7/8); % blasius law
r = [1.0-logspace(0,log10(1/ret0),kmax-1),1.0];
yr = 1.0-r;
uth1 = n*r;
uth2 = n*r.^2.0;
uth = zeros(1,kmax);
druthdr = zeros(1,kmax);
duthrdr = zeros(1,kmax);
dupdr = zeros(1,kmax);
nut = zeros(1,kmax);
ri = zeros(1,kmax);
strain = zeros(1,kmax);
% u and uth, outer-scaled
% up and uthp, inner-scaled
% r and yr, outer-scaled
while 1
  ret = 0.5*(retlim(1)+retlim(2));
  if lblend == 1
    for k = 1:kmax
      x = yr(k)*ret/locpblend;
      uth(k) = 1/(1+x^2)*(uth1(k) + uth2(k)*x^2);
    end
  else
    uth = uth2;
  end
  % first/second-order differentiation
  % d(r*uth)/d(r)
  druthdr(1) = 0.0;
  druthdr(2:kmax-1) = (r(3:kmax).*uth(3:kmax)-r(1:kmax-2).*uth(1:kmax-2))...
                      ./(r(3:kmax)-r(1:kmax-2));
  druthdr(kmax) = (r(kmax)*uth(kmax)-r(kmax-1)*uth(kmax-1))...
                  ./(r(kmax)-r(kmax-1));
  % d(uth/r)/d(r)
  duthrdr(1) = (uth(2)/r(2)-0.0)...
                  ./(r(2)-r(1));
  duthrdr(2) = (uth(3)/r(3)-0.0)...
                      ./(r(3)-r(1));
  duthrdr(3:kmax-1) = (uth(4:kmax)./r(4:kmax)-uth(2:kmax-2)./r(2:kmax-2))...
                      ./(r(4:kmax)-r(2:kmax-2));
  duthrdr(kmax) = (uth(kmax)/r(kmax)-uth(kmax-1)/r(kmax-1))...
                  ./(r(kmax)-r(kmax-1));
  % d(r*uthp)/d(r)
  % d(uthp/r)/d(r)
  uthp = uth*reb/(2.0*ret);
  druthpdr = druthdr*reb/(2.0*ret);
  duthprdr = duthrdr*reb/(2.0*ret);
  for k = 1:kmax
    x = r(k);
    t1 = uthp(k);
    t2 = druthpdr(k);
    t3 = duthprdr(k);
    t4 = 2.0*t1/x^2*t2/(y^2+(x*t3)^2);
    if lusez == 1
        t3 = zz;
        t4 = 6*(x*zz)^2/(y^2+(x*zz)^2);
    end
    if x == 0
      t4 = 0;
    end
    da = (1-exp(-(1-x)*ret/26));
    if lda == 0
      da = 1.0;
    end
    l0 = kap*(0.35-0.2*x^2-0.15*x^4)*da;
    l = (1-beta*t4)^alph*l0;
    if lvis == 1 && (1.0-x)*ret <= locvis
      l = 0;
    end
    if leqn == 1
      % Wei et al.
      eqn = l^2*sqrt(y^2+(x*t3)^2)*y == -x - 1/ret*y;
    elseif leqn == 2
      % Kiku et al., Reich and Beer
      eqn = l^2*sqrt(1+(x*t3/y)^2)*y^2 == x + 1/ret*y;
    end
    % one negative value
    dupdr(k) = min(vpasolve(eqn,y,-ret));
    % eddy viscosity (nondimensionless)
    nut(k) = ret*subs(l^2*sqrt(y^2+(x*t3)^2),y,dupdr(k));
    % richardson number
    ri(k) = subs(t4,y,dupdr(k));
    % strain rate (nondimensionless)
    strain(k) = 2*ret/reb*subs(sqrt(y^2+(x*t3)^2),y,dupdr(k));
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
% fname = ['REB',num2str(reb,'%i'),... 
%          '_N',num2str(n,'%.1f'),...
%          '_RET',num2str(ret,'%.2f'),...
%          '_KMAX',num2str(kmax,'%i')];
% save([fname,'.mat'],'r','u');
% save([fname,'_NUT','.mat'],'r','nut')
% save([fname,'_RI','.mat'],'r','ri')
% save([fname,'_STRAIN','.mat'],'r','strain')
% save([fname,'_UTH','.mat'],'r','uth')

fname = ['REB',num2str(reb,'%i'),...
         '_Z',num2str(zz,'%i'),...
         '_RET',num2str(ret,'%.2f'),...
         '_KMAX',num2str(kmax,'%i')];
save([fname,'.mat'],'r','u');

fid = fopen(['uclp_REB',num2str(reb),'_Z',num2str(zz,'%03.f'),'.dat'],'w');
fprintf(fid,'%15.10e %15.10e %15.10e %15.10e\n',...
       [reb,zz,ret,up(1)]);
