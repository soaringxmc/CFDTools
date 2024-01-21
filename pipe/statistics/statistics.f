      program avewsca
      parameter(nvar=30)
      integer,parameter :: dp = 8
      dimension y(10000),yc(10000),w(10000),wav(10000)
      dimension wvar(10000,nvar),wvarav(10000,nvar),
     .          prod(10000),diss(10000),prodt(10000)
      dimension cfl(10000,3),cflav(10000,3)
      dimension urms(10000)
      character*3 nc2
      character*4 nc3
      print *,'start,end'
      read(*,*) ibeg,iend
      read(*,*) dt
      read(*,*) f0
      open(11,file='bou.in',form='formatted')
      do l=1,9
       read(11,*)
      enddo
      read(11,*) re, ros, pr
      read(11,*)
      read(11,*) inslwr, icorio
      do l=1,11
       read(11,*)
      enddo
      read(11,*) uwall
      close(11)
!      
      if (icorio.eq.0) ros = 0. 
!     as long as ros.ne.0, rotating frame used
      pe = re*pr ! Peclet number
!     dt = 10.
      tbeg=ibeg*dt
      tend=iend*dt
      open(11,file='pipe.out',form='formatted')
      open(111,file='ff.out',form='formatted')
      iav   = 0 
      tauwav = 0.
      tgradav= 0.
      ub     = 0.
      tb     = 0. ! mean temperature
      open(12,file='cfav.out',form='formatted')
      do
       read(11,*,end=200) telaps,dtt,div,rmfr,pgrad,rmtr,tgrad
       if (telaps.ge.tbeg.and.telaps.le.tend) then
        tauw=-pgrad/2.
        tauwav = (tauwav*iav+tauw)/(iav+1.)
        tgradav = (tgradav*iav+tgrad)/(iav+1.)
        ub = (ub*iav+rmfr)/(iav+1.)  !rmfr is ub
        tb = (tb*iav+rmtr)/(iav+1.)
        cfav = 2.*tauwav/ub**2
        retau= re*sqrt(tauwav)
        write(12,*) telaps,cfav,ub,retau
        iav   = iav  +1   
        ff = 8.*tauw/rmfr**2
        reb = 2*ub*re
        peb = reb*pr
        rnu = -2*peb*tgrad/(4*ub*tb) ! Nusselt number (based on Tb)
        write(111,*) telaps,ff,rnu
       endif
      enddo
 200  close(12)
      open(21,file='cf.dat',form='formatted')
      print *,'Average Cf',cfav      !friction coef Cf = tau/(0.5*ub**2)
      print *,'Average f',4*cfav     !friction factor f = 4*Cf
      write(21,*) 'start,end', ibeg,iend
      write(21,*) 'Average Cf',cfav
      write(21,*) 'Average f',4*cfav
!     ub=0.5
      utau=sqrt(tauwav)
      delv=1./(re*utau)
      ttau=-tgradav/2./utau ! friction temperature
      print *,'utau ',utau
      print *,'Reb  ',2*re*ub
      reb = 2*re*ub
      print *,'Retau',1./delv
      retau = 1./delv
      write(21,*) 'utau  ',utau
      write(21,*) 'Ttau ',ttau
      write(21,*) 'Reb  ',2*re*ub
      write(21,*) 'Retau',1./delv
      print *, 'Turnover times',(iend-ibeg+1)*utau 
      write(21,*) 'Turnover times',(iend-ibeg+1)*utau
      print *,'Ttau',ttau
      close(11)
      write(111,*)
      close(111)
!
!     Start averaging flow samples
!
      wvarav= 0.
      open(22,file='conv_peak.dat',form='formatted') ! check convergence of running averages
      open(23,file='uncertainty.dat',form='formatted') ! write statistics of flow samples for uncertainty analysis
!
      do l=ibeg,iend
       write(nc2,1003) l
       write(nc3,1004) l
!      print *,nc2
!      if (l.le.999) then
!       open(11,file='wmean_'//nc2//'.dat',form='formatted')
!      else
        open(11,file='wmean_'//nc3//'.dat',form='formatted')
!      endif 
!      open(12,file='cfl_'//nc2//'.dat',form='formatted')
       j=0
       do
!       print *,j
        j=j+1
!       read(11,*,end=300) y(j),t(j)
        read(11,*,end=300) y(j),yc(j),(wvar(j,n),n=1,nvar)
!       read(12,*,end=300) y(j),yc(j),(cfl(j,n),n=1,3)
       enddo
 300   continue
!      print *,'lines',j
       nli=j-1  ! nli=m2m
       close(11)
       close(12)
!
       y(1)=0.
       y(nli+1)=1.
!
!      Update averages
!
       wvarav = wvarav+wvar
       cflav = cflav+cfl
!  
!      Running average of u'^2
!
       urms(:)=wvarav(:,6)/(l-ibeg+1)-(wvarav(:,3)/(l-ibeg+1))**2
!
       umax = 0.
       do j=3,nli
        if (urms(j).gt.umax) then
         jj = j
         umax = urms(j) 
        endif
       enddo 
!      jj = j-1 ! peak location
       dym = yc(jj)-yc(jj-1)
       dyp = yc(jj+1)-yc(jj)
       fp  = urms(jj+1)
       fm  = urms(jj-1)
       ff  = urms(jj)
       bb  = ((fp-ff)*dym**2-(fm-ff)*dyp**2)/(dym*dyp*(dyp+dym)) 
       cc  = ((fm-ff)*dyp+(fp-ff)*dym)/(dym*dyp*(dyp+dym)) 
       u2peak = ff-bb**2/(4*cc)
       ypk = (1.-(yc(jj)-bb/(2*cc)))/delv
!      print *,l,urms(jj),maxval(urms),u2peak,jj,-bb/(2*cc)
!      uprmspk=maxval(urms)/utau**2
       uprmspk=u2peak/utau**2
!
!      Running average of ucl
!
       y1=yc(1)
       y2=yc(2)
       u1=wvarav(1,3)/(l-ibeg+1)
       u2=wvarav(2,3)/(l-ibeg+1)
       ucl = (u1*y2**2-u2*y1**2)/(y2**2-y1**2)
       uclp = ucl/utau
!
!      Write running averages
       write(22,*) l,uprmspk,uclp
!
!      Write space averaged statistics for uncertainty estimation
!
       urms(:)=wvar(:,6)-wvar(:,3)**2 ! spatially averaged velocity variance
!
!      Peak velocity variance
!
       umax = 0.
       do j=3,nli
        if (urms(j).gt.umax) then
         jj = j
         umax = urms(j) 
        endif
       enddo 
!      jj = j-1 ! peak location
       dym = yc(jj)-yc(jj-1)
       dyp = yc(jj+1)-yc(jj)
       fp  = urms(jj+1)
       fm  = urms(jj-1)
       ff  = urms(jj)
       bb  = ((fp-ff)*dym**2-(fm-ff)*dyp**2)/(dym*dyp*(dyp+dym)) 
       cc  = ((fm-ff)*dyp+(fp-ff)*dym)/(dym*dyp*(dyp+dym)) 
       u2peak = ff-bb**2/(4*cc)
       ypk = (1.-(yc(jj)-bb/(2*cc)))/delv
       uprmspk=u2peak/utau**2
!
!      Wall dissipation
!
       f1=0.5*urms(nli)/utau**2
       f2=0.5*urms(nli-1)/utau**2
       y1=(1.-yc(nli))/delv
       y2=(1.-yc(nli-1))/delv
       d2u2dy2=2*(-f1*(y2-y1)+(f2-f1)*y1)/((y2-y1)*y1*y2)  
!      print *,'wall dissipation of u',d2u2dy2
!
!      Peak temperature variance
!
       urms(:)=wvar(:,9)-wvar(:,8)**2 ! spatially averaged temperature variance
!
       umax = 0.
       do j=3,nli
        if (urms(j).gt.umax) then
         jj = j
         umax = urms(j) 
        endif
       enddo 
!      jj = j-1 ! peak location
       dym = yc(jj)-yc(jj-1)
       dyp = yc(jj+1)-yc(jj)
       fp  = urms(jj+1)
       fm  = urms(jj-1)
       ff  = urms(jj)
       bb  = ((fp-ff)*dym**2-(fm-ff)*dyp**2)/(dym*dyp*(dyp+dym)) 
       cc  = ((fm-ff)*dyp+(fp-ff)*dym)/(dym*dyp*(dyp+dym)) 
       t2peak = ff-bb**2/(4*cc)
       ytpk = (1.-(yc(jj)-bb/(2*cc)))/delv
       tprmspk=t2peak/ttau**2
!
!      Centerline velocity
!
       y1=yc(1)
       y2=yc(2)
       u1=wvar(1,3)
       u2=wvar(2,3)
       ucl = (u1*y2**2-u2*y1**2)/(y2**2-y1**2)
       uclp = ucl/utau
!
!      Centerline temperature
!
       y1=yc(1)
       y2=yc(2)
       t1=wvar(1,8)
       t2=wvar(2,8)
       tcl = (t1*y2**2-t2*y1**2)/(y2**2-y1**2)
       tclp = tcl/ttau
!
       write(23,*) l,uprmspk, ypk,uclp,d2u2dy2,
     .               tprmspk,ytpk,tclp
!
      enddo  ! end average
!
      close(22)
      write(23,*) 
      close(23)
!
      wvarav=wvarav/(iend-ibeg+1)
      cflav=cflav/(iend-ibeg+1)
!     average of uth and uth2 on the wall
      y1=1.d0-yc(nli)
      y2=1.d0-yc(nli-1)
      y3=1.d0-yc(nli-2)
      f1=wvarav(nli,1)
      f2=wvarav(nli-1,1)
      f3=wvarav(nli-2,1)
      uthw = (f1*y2*y3*(y2-y3)+f2*y1*y3*(y3-y1)+f3*y1*y2*(y1-y2))/
     .        ((y2-y1)*(y2-y3)*(y3-y1))    
      print *,'wall azimuthal velocity',uthw
      duthw = (wvarav(nli,1)-uthw)/y1
      print *,'wall derivative of azimuthal velocity',duthw
!
      f1=wvarav(nli,4)
      f2=wvarav(nli-1,4)
      f3=wvarav(nli-2,4)

!     print *,'f1,f2,f3',f1,f2,f3
!     print *,'y1,y2,y3',y1,y2,y3
      uth2w = (f1*y2*y3*(y2-y3)+f2*y1*y3*(y3-y1)+f3*y1*y2*(y1-y2))/
     .     ((y2-y1)*(y2-y3)*(y3-y1))
      print *,'wall azimuthal velocity variance',uth2w
      duth2w = (wvarav(nli,4)-uth2w)/y1
      print *,'wall derivative of azimuthal velocity variance',duth2w
!
!     Redefine statistical moments
!
      wvarav(:,14)=wvarav(:,10) ! this is now P
      wvarav(:,15)=sqrt(wvarav(:,11)-wvarav(:,14)**2) ! this is now p'
      wvarav(:,4)=sqrt(wvarav(:,4)-wvarav(:,1)**2)
      wvarav(:,5)=sqrt(wvarav(:,5)-wvarav(:,2)**2)
      wvarav(:,6)=sqrt(wvarav(:,6)-wvarav(:,3)**2)
!     wvarav(:,7)=     wvarav(:,7)-wvarav(:,2)*wvarav(:,3)
      wvarav(:,10)=sqrt(wvarav(:,9)-wvarav(:,8)**2) ! T' in now #10
      wvarav(:,9)= wvarav(:,8) ! T is now #9
      do j=2,nli
       uavn=0.5*(wvarav(j,3)+wvarav(j-1,3))
       dtdyp=(wvarav(j+1,9)-wvarav(j,9))/(yc(j+1)-yc(j)) 
       dtdym=(wvarav(j,9)-wvarav(j-1,9))/(yc(j)-yc(j-1))
       dtdy =0.5*(dtdyp+dtdym) ! derivative at cell center
       wvarav(j,13)=dtdy/pe
!      uavn=wp*wvarav(j,3)+wm*wvarav(j-1,3) 
       wvarav(j,7)=wvarav(j,7)-uavn*wvarav(j,2) !defined at grid points
       dudyn=(wvarav(j,3)-wvarav(j-1,3))/(yc(j)-yc(j-1))
       wvarav(j,8)=dudyn/re     !8:tau_z
      enddo
!
!     Outer-scaled statistics
!
      open(11,file='stats_out.dat',form='formatted')
      do j=1,nli
       write(11,100) y(j),yc(j),(wvarav(j,n),n=1,nvar)
      enddo
      close(11)
!     
      open(11,file='friction_out.dat',form='formatted')
!     if (ros.ne.0.) then
!      dudyn  = (wvarav(nli,1)-0.)/(yc(nli)-1.) + ros  !theta-direction
!     elseif (uwall.ne.0.) then
!      dudyn  = (wvarav(nli,1)-uwall)/(yc(nli)-1.)  !theta-direction
!     endif
      dudyn  = (wvarav(nli,1)-uwall)/(yc(nli)-1.) + ros
!     inertial frame, ros=0, uwall>0
!     rotating frame, ros>0, uwall=0
!     rotating frame+travelling waves, ros>0, uwall>0, but
!     the result may be useless if opposite dudyn signs appear      
      tauw_u = dudyn/re
      utau_u = sqrt(tauw_u)
      retau_u= re*utau_u
      cftav = 2.*tauw_u/ub**2
      write(21,*) 'Average Cf_theta',cftav
!
      dudyn  = (wvarav(nli,3)-0.)/(yc(nli)-1.)  !z-direction
      tauw_w = dudyn/re
      utau_w = sqrt(-tauw_w)
      retau_w= re*utau_w
      
      write(11,100) tauw_u,utau_u,retau_u
      write(11,100) tauw_w,utau_w,retau_w
      
      close(11)
!
!     Momentum balance
!
!     Fix issue with first two points ov vT stress
!
!     vtrat= wvarav(nli-3,12)/wvarav(nli-3,7)
      vtrat= wvarav(nli-3,12)/(1.d0-yc(nli-3))**3
      print *,'nli,vtrat',nli,vtrat
      do j=nli-2,nli
!      wvarav(j,12)=vtrat*wvarav(j,7)
       wvarav(j,12)=vtrat*(1.d0-yc(j))**3
      enddo
      uvrat= wvarav(nli-1,7)/(1.d0-y(nli-1))**3
      wvarav(nli,7)=uvrat*(1.d0-y(nli))**3
!
      open(11,file='mombal.dat',form='formatted')
      do j=2,nli
       write(11,100) y(j),wvarav(j,7)/utau**2,wvarav(j,8)/utau**2,
     .    yc(j),wvarav(j,12)/utau/ttau,wvarav(j,13)/utau/ttau
      enddo
      close(11)
!
!     Turbulent Prandtl number
!
      open(11,file='prt.dat',form='formatted')
      do j=2,nli
       dudyc = 0.5d0*(wvarav(j+1,8)+wvarav(j,8)) ! dudy at cell center
!      uvc   = 0.5d0*(wvarav(j+1,7)+wvarav(j,7)) ! shear stress at cell center
       uvc   = 
     .  0.5d0*(wvarav(j+1,7)/(1.d0-y(j+1))**3+
     .         wvarav(  j,7)/(1.d0-y(  j))**3)*(1.d0-yc(j))**3 ! shear stress at cell center
       dtdyc = wvarav(j,13) ! dTdy at cell center
       tvc   = wvarav(j,12) ! heat flux at cell center
       prt = (uvc/dudyc)/(tvc/dtdyc)
       uvcbal = yc(j)+re*dudyc/utau*delv ! uv from momentum balance
       tvcbal = yc(j)+re*dtdyc/ttau*delv ! vT from heat balance
       prt1=(uvcbal/(re*dudyc/utau*delv))/(tvcbal/(re*dtdyc/ttau*delv))
       write(11,100) 1.d0-yc(j),(1.d0-yc(j))/delv,
     .    prt,prt1,
     .    uvc/utau**2,-re*dudyc/utau*delv,tvc/utau/ttau,
     .    -re*dtdyc/ttau*delv
      enddo
      close(11) 
!
!     TKE balance
!
      epsnorm = utau**3/delv ! wall units for dissipation
      print *,'epsnorm',epsnorm
      open(11,file='tkebal.dat',form='formatted')
      do j=2,nli ! interpolate to nodes
       e11 = 0.5*(wvarav(j,19)+wvarav(j-1,19))
       e22 =      wvarav(j,20)
       e33 = 0.5*(wvarav(j,21)-wvarav(j,22)
     .           +wvarav(j-1,21)-wvarav(j-1,22))
       eps = e11+e22+e33
       prod(j)=-wvarav(j,7)*wvarav(j,8)/utau**4
       diss(j)=eps/epsnorm
!      wvarav(j,17)=eps/epsnorm   ! total dissipation
!      wvarav(j,18)=-wvarav(j,7)*wvarav(j,8)/epsnorm  ! production
       ut2 = 0.5*(wvarav(j,4)**2+wvarav(j-1,4)**2)
       ur2 =      wvarav(j,5)**2
       uz2 = 0.5*(wvarav(j,6)**2+wvarav(j-1,6)**2)
       qq = 0.5*(ut2+ur2+uz2) ! kinetic energy
       write(11,100) (1.-y(j))/delv,-wvarav(j,7)*wvarav(j,8)/utau**4,
     .             eps/epsnorm,-qq/eps/delv*utau
      enddo
      close(11)
!
!     Thermal energy balance
!
      epsnorm = utau*ttau**2/delv ! wall units for dissipation
      print *,'epsnorm',epsnorm
      open(11,file='enebal.dat',form='formatted')
      do j=2,nli ! interpolate to nodes
       eps = (wvarav(j,30)-wvarav(j,29)**2)/epsnorm/re
       prodt(j)=-wvarav(j,29)*wvarav(j,12)/epsnorm
       dt2dy =(wvarav(j+1,10)**2-wvarav(j  ,10)**2)/(yc(j+1)-yc(j  ))
       dt2dym=(wvarav(j  ,10)**2-wvarav(j-1,10)**2)/(yc(j  )-yc(j-1))
       d2dt2=(y(j+1)*dt2dy-y(j)*dt2dym)/(y(j+1)-y(j))
       vdiff=-0.5*d2dt2/epsnorm/re/yc(j)
       write(11,100) (1.-yc(j))/delv,prodt(j),eps,vdiff,
     .     wvarav(j,30),wvarav(j,29)**2,wvarav(j,29),wvarav(j,12) 
      enddo
      close(11)
!
!     r.m.s. vorticity
!
      open(11,file='vort.dat',form='formatted')
      do j=2,nli
       omt = sqrt(wvarav(j,26)-wvarav(j,23)**2)/utau*delv
       omr = sqrt(wvarav(j,27)-wvarav(j,24)**2)/utau*delv
       omz = sqrt(wvarav(j,28)-wvarav(j,25)**2)/utau*delv
       write(11,100) (1.-yc(j))/delv,(1.-yc(j))/delv,
     .    omt,omr,omz
      enddo
      close(11)  
!
!     TKE production peak
!
      umax = 0.
      do j=3,nli
       if (prod(j).gt.umax) then
        jj = j
        umax = prod(j) 
       endif
      enddo 
!     jj = j-1 ! peak location
      dym = yc(jj)-yc(jj-1)
      dyp = yc(jj+1)-yc(jj)
      fp  = prod(jj+1)
      fm  = prod(jj-1)
      ff  = prod(jj)
      bb  = ((fp-ff)*dym**2-(fm-ff)*dyp**2)/(dym*dyp*(dyp+dym)) 
      cc  = ((fm-ff)*dyp+(fp-ff)*dym)/(dym*dyp*(dyp+dym)) 
      prodpeak = ff-bb**2/(4*cc)
      ypkprod = (1.-(yc(jj)-bb/(2*cc)))/delv
      print *,'Production peak',prodpeak,ypkprod
!
!     Dissipation peak
!
      umax = 0.
      do j=3,nli
       if ((-diss(j)).gt.umax) then
        jj = j
        umax = -diss(j) 
       endif
      enddo 
!     jj = j-1 ! peak location
      dym = yc(jj)-yc(jj-1)
      dyp = yc(jj+1)-yc(jj)
      fp  = -diss(jj+1)
      fm  = -diss(jj-1)
      ff  = -diss(jj)
      bb  = ((fp-ff)*dym**2-(fm-ff)*dyp**2)/(dym*dyp*(dyp+dym)) 
      cc  = ((fm-ff)*dyp+(fp-ff)*dym)/(dym*dyp*(dyp+dym)) 
      disspeak = ff-bb**2/(4*cc)
      ypkdiss = (1.-(yc(jj)-bb/(2*cc)))/delv
      print *,'Dissipation peak',disspeak,ypkdiss
!
!     Diagnostic function 
!
      open(11,file='diagnostic.dat',form='formatted')
      do j=3,nli-1
       wp = (y(j)-yc(j-1))/(yc(j)-yc(j-1))
       wm = 1.-wp
!      uavn=0.5*(wvarav(j,3)+wvarav(j-1,3))
       uavn=wp*wvarav(j,3)+wm*wvarav(j-1,3) 
       dudyn=-(wvarav(j,3)-wvarav(j-1,3))/(yc(j)-yc(j-1))
       yndudyn=(1.-y(j))*dudyn
       ypdudyp=yndudyn/utau
       ypdudypbyup=ypdudyp/uavn*utau
       tavn=wp*wvarav(j,9)+wm*wvarav(j-1,9) 
       dtdyn=-(wvarav(j,9)-wvarav(j-1,9))/(yc(j)-yc(j-1))
       yndtdyn=(1.-y(j))*dtdyn
       ypdtdyp=yndtdyn/ttau
       ypdtdypbytp=ypdtdyp/tavn*ttau
       write(11,100) (1.-y(j)),(1.-y(j))/delv,ypdudyp,ypdudypbyup,
     .             dudyn/utau*delv,uavn/utau,
     .             ypdtdyp,ypdtdypbytp,
     .             dtdyn/ttau*delv,tavn/ttau
      enddo
      close(11)
!
!     Defect profile and centerline velocity
!
      open(11,file='defect.dat',form='formatted')
      y1=yc(1)
      y2=yc(2)
      u1=wvarav(1,3)
      u2=wvarav(2,3)
      ucl = (u1*y2**2-u2*y1**2)/(y2**2-y1**2)
      uclp = ucl/utau
      t1=wvarav(1,9)
      t2=wvarav(2,9)
      tcl = (t1*y2**2-t2*y1**2)/(y2**2-y1**2)
      tclp = tcl/ttau
!     print *,'y1,y2',y1,y2
!     print *,'u1,u2',u1,u2
!
      do j=1,nli
       write(11,*) (1.-yc(j)),(1.-yc(j))/delv,
     .             (ucl-wvarav(j,3))/ub,
     .             (ucl-wvarav(j,3))/utau,
     .             (tcl-wvarav(j,9))/tb,
     .             (tcl-wvarav(j,9))/ttau
      enddo
      close(11)
!
!     Compute bulk temperature
!
      tm = 0.
      ubb= 0.
      tbb= 0.
      do j=1,nli
       rl = yc(j)
       dy = -y(j)+y(j+1)
       tm = tm + rl*wvarav(j,3)*wvarav(j,9)*dy
       ubb= ubb+ rl*wvarav(j,3)*dy
       tbb= tbb+ rl*wvarav(j,9)*dy
!      print *, rl,dy,wvarav(j,3),wvarav(j,9),tm
      enddo
      tm = 2*tm/ub ! bulk temperature
      ubb= 2*ubb   ! bulk velocity
      tbb= 2*tbb   ! mean temperature
      print *,'Bulk temp',tm
      print *,'Bulk velocity',ubb
      print *,'Mean temperature',tbb
!     stop
!
!     Inner-scaled statistics
!
!     Compute production and dissipation
!
      do j=2,nli ! interpolate to nodes
       e11 =      wvarav(j,19)
       e22 = 0.5*(wvarav(j,20)+wvarav(j,20))
       e33 =      wvarav(j,21)-wvarav(j,22)
       eps = e11+e22+e33
       wvarav(j,17)=eps/epsnorm   ! total dissipation
       wvarav(j,18)=prod(j)  ! production
      enddo
      close(11)
!
!     Wall pressure
!
      f1=wvarav(nli,14)
      f2=wvarav(nli-1,14)
      f3=wvarav(nli-2,14)
      y1=1.-yc(nli)
      y2=1.-yc(nli-1)
      y3=1.-yc(nli-2)
      pw = (f1*y2*y3*(y2-y3)+f2*y1*y3*(y3-y1)+f3*y1*y2*(y1-y2))/
     .     ((y2-y1)*(y2-y3)*(y3-y1))
      print *,'f1,f2,f3',f1,f2,f3
      print *,'wall pressure',pw
      wvarav(:,14)=wvarav(:,14)-pw ! subtract wall pressure
!
!     Wall pressure variance
!
      f1=wvarav(nli,15)**2
      f2=wvarav(nli-1,15)**2
      f3=wvarav(nli-2,15)**2
      y1=1.-yc(nli)
      y2=1.-yc(nli-1)
      y3=1.-yc(nli-2)
      prms2w = (f1*y2*y3*(y2-y3)+f2*y1*y3*(y3-y1)+f3*y1*y2*(y1-y2))/
     .     ((y2-y1)*(y2-y3)*(y3-y1))
      print *,'f1,f2,f3',f1,f2,f3
      print *,'wall pressure variance',prms2w
      prms2w=prms2w/utau**4
!
      open(11,file='stats_inn.dat',form='formatted')
      wvarav(:,1:6)=wvarav(:,1:6)/utau
      wvarav(:,7:8)=wvarav(:,7:8)/utau**2  !7:uruz, 8:tauz
      wvarav(:,9:10)=wvarav(:,9:10)/ttau
      wvarav(:,12)=wvarav(:,12)/utau/ttau
      wvarav(:,13)=wvarav(:,13)/utau/ttau
      wvarav(:,14)=wvarav(:,14)/utau**2    !
      wvarav(:,15)=wvarav(:,15)/utau**2    !
      do j=nli,1,-1
       write(11,100)
     . (1.-y(j))/delv,(1.-yc(j))/delv,(wvarav(j,n),n=1,nvar)
      end do
      close(11)
      
     !! open(11,file='rotating_frame.dat')
     !! do j=nli,1,-1
     !!  write(11,100)
     !!. (1.-yc(j))/delv,wvarav(j,1)-uwall/utau*yc(j)
     !! end do
     !! close(11)
!
      open(11,file='Pipe_profile.dat',form='formatted')
      wvarav(1,5)=2*wvarav(2,5)-wvarav(3,5) !  extrapolate v'
      do j=nli,1,-1
       v2c  = (0.5*(wvarav(j+1,5)+wvarav(j,5)))**2
       uvc  = 0.5*(wvarav(j+1,7)+wvarav(j,7))
       write(11,100)
!    . (1.-y(j))/delv,(1.-yc(j))/delv,
     . (1.-yc(j))/delv,
     .  wvarav(j,3),  ! U
     .  wvarav(j,14), ! P
     .  wvarav(j,6)**2, ! u_z**2
     .  v2c,            ! u_r**2
     .  wvarav(j,4)**2, ! u_t**2
     .  uvc,            ! u_z u_r 
     .  wvarav(j,15)**2  ! p**2
      enddo
      close(11)
!
!     Reynolds stresses at cell centers
!
      open(11,file='reystress_center.dat',form='formatted')
      do j=nli,1,-1
       uc   = wvarav(j,3)
       u2c  = wvarav(j,6)**2
       v2c  = (0.5*(wvarav(j+1,5)+wvarav(j,5)))**2
       w2c  = wvarav(j,4)**2
       uvc  = 0.5*(wvarav(j+1,7)+wvarav(j,7))
       tc   = wvarav(j,9)
       t2c  = wvarav(j,10)**2
       tvc  = wvarav(j,12)
       write(11,100) (1.-yc(j))/delv,
     .          uc,u2c,v2c,w2c,uvc,
     .          tc,t2c,tvc
!     uc : <Uz>
!     u2c: <uzuz>
!     v2c: <urur>
!     w2c: <utut>
!     uvc: <uzur>   0.5*(wvarav(j+1,7)+wvarav(j,7))
!     tc : <T>
!     t2c: <tt>
!     tvc: <tur>
!     uwc: <uzut>   waiting for adding
!     vwc: <urut>   waiting for adding
      enddo
      close(11)
!
!     Time step analysis
!     
      open(11,file='cfl_inn.dat',form='formatted')
      cflav(:,:)=1./cflav(:,:)*utau/delv
      do j=nli,1,-1
       write(11,100) (1.-y(j))/delv,(1.-yc(j))/delv,(cflav(j,n),n=1,3)
      enddo
      close(11)
!
      close(21)
!
!     Wall dissipation
! 
      f1=0.5*(wvarav(nli,4)**2+wvarav(nli,5)**2+wvarav(nli,6)**2)
      f2=0.5*(wvarav(nli-1,4)**2+wvarav(nli-1,5)**2+wvarav(nli-1,6)**2)
      y1=(1.-yc(nli))/delv
      y2=(1.-yc(nli-1))/delv
      d2kdy2=2*(-f1*(y2-y1)+(f2-f1)*y1)/((y2-y1)*y1*y2)  
      etapw = (1./d2kdy2)**(0.25)
!
      f1=0.5*(wvarav(nli,6)**2)
      f2=0.5*(wvarav(nli-1,6)**2)
      y1=(1.-yc(nli))/delv
      y2=(1.-yc(nli-1))/delv
      d2u2dy2=2*(-f1*(y2-y1)+(f2-f1)*y1)/((y2-y1)*y1*y2)  
!
!     Peak velocity variance
!
!     Streamwise velocity
!
      urms(:)=wvarav(:,6)**2
      umax = 0.
      do j=3,nli
       if (urms(j).gt.umax) then
        jj = j
        umax = urms(j) 
       endif
      enddo 
!     jj = j-1 ! peak location
      dym = yc(jj)-yc(jj-1)
      dyp = yc(jj+1)-yc(jj)
      fp  = urms(jj+1)
      fm  = urms(jj-1)
      ff  = urms(jj)
      bb  = ((fp-ff)*dym**2-(fm-ff)*dyp**2)/(dym*dyp*(dyp+dym)) 
      cc  = ((fm-ff)*dyp+(fp-ff)*dym)/(dym*dyp*(dyp+dym)) 
      u2peak = ff-bb**2/(4*cc)
      ypk = (1.-(yc(jj)-bb/(2*cc)))/delv
      uprmspk=u2peak
!
!     Wall-normal velocity
!
      urms(:)=wvarav(:,5)**2
      umax = 0.
      do j=3,nli
       if (urms(j).gt.umax) then
        jj = j
        umax = urms(j) 
       endif
      enddo 
!     jj = j-1 ! peak location
      dym = y(jj)-y(jj-1)
      dyp = y(jj+1)-y(jj)
      fp  = urms(jj+1)
      fm  = urms(jj-1)
      ff  = urms(jj)
      bb  = ((fp-ff)*dym**2-(fm-ff)*dyp**2)/(dym*dyp*(dyp+dym)) 
      cc  = ((fm-ff)*dyp+(fp-ff)*dym)/(dym*dyp*(dyp+dym)) 
      v2peak = ff-bb**2/(4*cc)
      ypkv = (1.-(y(jj)-bb/(2*cc)))/delv
      vprmspk=v2peak
!
!     Azimuthal velocity
!
      urms(:)=wvarav(:,4)**2
      umax = 0.
      do j=3,nli
       if (urms(j).gt.umax) then
        jj = j
        umax = urms(j) 
       endif
      enddo 
!     jj = j-1 ! peak location
      dym = yc(jj)-yc(jj-1)
      dyp = yc(jj+1)-yc(jj)
      fp  = urms(jj+1)
      fm  = urms(jj-1)
      ff  = urms(jj)
      bb  = ((fp-ff)*dym**2-(fm-ff)*dyp**2)/(dym*dyp*(dyp+dym)) 
      cc  = ((fm-ff)*dyp+(fp-ff)*dym)/(dym*dyp*(dyp+dym)) 
      w2peak = ff-bb**2/(4*cc)
      ypkw = (1.-(yc(jj)-bb/(2*cc)))/delv
      wprmspk=w2peak
!
!     Shear stress
!
      urms(:)= wvarav(:,7)
      umax = 0.
      do j=3,nli
       if (urms(j).gt.umax) then
        jj = j
        umax = urms(j) 
       endif
      enddo 
!     jj = j-1 ! peak location
      dym = y(jj)-y(jj-1)
      dyp = y(jj+1)-y(jj)
      fp  = urms(jj+1)
      fm  = urms(jj-1)
      ff  = urms(jj)
      bb  = ((fp-ff)*dym**2-(fm-ff)*dyp**2)/(dym*dyp*(dyp+dym)) 
      cc  = ((fm-ff)*dyp+(fp-ff)*dym)/(dym*dyp*(dyp+dym)) 
      uvpeak = ff-bb**2/(4*cc)
      ypkuv = (1.-(y(jj)-bb/(2*cc)))/delv
      uvppk=uvpeak
!
!     Heat flux
!
      urms(:)= wvarav(:,12)
      umax = 0.
      rewind(88) 
      do j=3,nli
       if (urms(j).gt.umax) then
        jj = j
        umax = urms(j) 
       endif
       write(88,*) yc(j),urms(j)
      enddo 
!     jj = j-1 ! peak location
      dym = yc(jj)-yc(jj-1)
      dyp = yc(jj+1)-yc(jj)
      fp  = urms(jj+1)
      fm  = urms(jj-1)
      ff  = urms(jj)
      bb  = ((fp-ff)*dym**2-(fm-ff)*dyp**2)/(dym*dyp*(dyp+dym)) 
      cc  = ((fm-ff)*dyp+(fp-ff)*dym)/(dym*dyp*(dyp+dym)) 
      vtpk = ff-bb**2/(4*cc)
      yvtpk = (1.-(yc(jj)-bb/(2*cc)))/delv
      print *,'peak of turbulent heat flux',vtpk,yvtpk
!
!     Peak temperature variance
!
      urms(:)=wvarav(:,10)**2
      umax = 0.
      do j=3,nli
       if (urms(j).gt.umax) then
        jj = j
        umax = urms(j) 
       endif
      enddo 
!     jj = j-1 ! peak location
      dym = yc(jj)-yc(jj-1)
      dyp = yc(jj+1)-yc(jj)
      fp  = urms(jj+1)
      fm  = urms(jj-1)
      ff  = urms(jj)
      bb  = ((fp-ff)*dym**2-(fm-ff)*dyp**2)/(dym*dyp*(dyp+dym)) 
      cc  = ((fm-ff)*dyp+(fp-ff)*dym)/(dym*dyp*(dyp+dym)) 
      t2peak = ff-bb**2/(4*cc)
      ytpk = (1.-(yc(jj)-bb/(2*cc)))/delv
      tprmspk=t2peak
!
!     Peak temperature variance
!
      urms(:)=wvarav(:,15)**2
      umax = 0.
      do j=3,nli
       if (urms(j).gt.umax) then
        jj = j
        umax = urms(j) 
       endif
      enddo 
!     jj = j-1 ! peak location
      dym = yc(jj)-yc(jj-1)
      dyp = yc(jj+1)-yc(jj)
      fp  = urms(jj+1)
      fm  = urms(jj-1)
      ff  = urms(jj)
      bb  = ((fp-ff)*dym**2-(fm-ff)*dyp**2)/(dym*dyp*(dyp+dym)) 
      cc  = ((fm-ff)*dyp+(fp-ff)*dym)/(dym*dyp*(dyp+dym)) 
      prmspk = ff-bb**2/(4*cc)
      yprmspk = (1.-(yc(jj)-bb/(2*cc)))/delv
!
!     Heat production peak
!
      umax = 0.
      do j=3,nli
       if (prodt(j).gt.umax) then
        jj = j
        umax = prodt(j) 
       endif
      enddo 
!     jj = j-1 ! peak location
      dym = yc(jj)-yc(jj-1)
      dyp = yc(jj+1)-yc(jj)
      fp  = prodt(jj+1)
      fm  = prodt(jj-1)
      ff  = prodt(jj)
      bb  = ((fp-ff)*dym**2-(fm-ff)*dyp**2)/(dym*dyp*(dyp+dym)) 
      cc  = ((fm-ff)*dyp+(fp-ff)*dym)/(dym*dyp*(dyp+dym)) 
      prodtpeak = ff-bb**2/(4*cc)
      ypkprodt = (1.-(yc(jj)-bb/(2*cc)))/delv
      print *,'Temperature production peak',prodtpeak,ypkprodt
!
      print *,'ucl',ucl
      print *,'ucl+',ucl/utau
      print *,'peak velocity variance',uprmspk
      print *,'peak velocity variance position',ypk
      print *,'wall dissipation',d2kdy2
      print *,'etaplus at wall',(1./d2kdy2)**(0.25) 
      print *,'wall dissipation of u',d2u2dy2
!
      write(21,*) 'ucl',ucl
      write(21,*) 'ucl+',ucl/utau
      write(21,*) 'peak velocity variance',uprmspk
      write(21,*) 'peak velocity variance position',ypk
!
      rnu = -2*peb*tgradav/(4*ub*tm) ! Nusselt number
!
      print *,'Mean temperature',tb
      print *,'Bulk temperature',tm
      print *,'Centerline temperature',tclp
      print *,'Peak temperature variance',tprmspk
      print *,'Nusselt number',rnu
!
      print *,'Wall pressure variance',prms2w
      print *,'Peak pressure variance',prmspk
      print *,'Peak pressure variance position',yprmspk
!
      open(11,file='flowprop.dat',form='formatted')
      write(11,*) reb, retau, 4*cfav, utau,
     .            uclp, 
     .            uprmspk, ypk,
     .            prodpeak,ypkprod, 
     .            disspeak,ypkdiss, 
     .            d2u2dy2, d2kdy2, etapw,
     .            tb,tm,ttau,rnu,tclp,tprmspk,ytpk,
     .            vprmspk,ypkv,
     .            wprmspk,ypkw,
     .            uvpeak,ypkuv,
     .            vtpk,yvtpk,
     .            prms2w,prmspk,yprmspk,pw,
     .            pr,prodtpeak,ypkprodt
      close(11)
!
      open(11,file='DR.dat',form='formatted')
      if (ros.ne.0.) then
!      rotating frame (+travelling waves)
       vrot = ros 
      elseif (uwall.ne.0.) then
!      inertial frame
       vrot = uwall
      endif
      dr=1.-4*cfav/f0      !f=4*cfav, friction factor
      print *,'Drag reduction',100*dr,'%'
      ps=dr-16*(vrot/ub)**2/f0/reb
      write(11,*) vrot/ub, reb, retau, 4*cfav, dr, ps
      print *,'Net power saving (rotation only,linear)',100*ps,'%'
      ps=dr-(vrot/ub)*(4*cftav)/f0
      write(11,*) vrot/ub, reb, retau, 4*cfav, dr, ps
      print *,'Net power saving (rotation only)',100*ps,'%'
!     interpolation of uthw, uth2w
      if (ros.ne.0.) then
!      rotating frame (+travelling waves)
       vrot = ros 
      else
!      inertial frame
       vrot = 0.
      endif      
      cft2_1 = 4.*(vrot/ub)**2/reb
      cft2_2 =-2./reb/ub**2*duth2w
      cft2_3 = 4.*(vrot/ub)/reb/ub*(uthw-duthw)
      ps = dr-4.*cft2_1/f0-4.*cft2_2/f0-4.*cft2_3/f0
      write(11,*) vrot/ub, reb, retau, 4*cfav, dr, ps
      write(11,*) 4.*cft2_1,4.*cft2_2,4.*cft2_3
      print *,'Net power saving (rotation+travelling waves)',100*ps,'%'
      print *,'Rotation term cft2_1', 4.*cft2_1
      print *,'Travelling wave term cft2_2', 4.*cft2_2
      print *,'Cross-influence term cft2_3', 4.*cft2_3
      close(11)
!
  100 format(50E20.10)
 1003 format(I3.3)
 1004 format(I4.4)
 1000 format(30E20.10)
!
      stop
      end
