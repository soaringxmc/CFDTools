      program extract_slices
!
      integer, parameter :: dp = 8
!
      real(dp), parameter :: pi = 4.*atan(1.d0)
!
      integer, parameter :: nplanes=4
      integer, parameter :: nvars=5
!
      real(dp), allocatable, dimension(:,:,:,:) :: planetz
      real(dp), allocatable, dimension(:,:) :: planetzl
      real(dp), allocatable, dimension(:,:,:) :: planerz
      real(dp), allocatable, dimension(:,:,:) :: planetr,planetrl
      real(dp), allocatable, dimension(:) :: th,rc,rm,z,sth,cth,rr
      real(dp), allocatable, dimension(:,:) :: wvarav
      integer :: jplanes(nplanes)
      character(1) :: naplane
      character(100) :: napltr,napltz,naplrz
      logical :: lfluc,luz,luth,lur
      integer :: lnd

!     instantaneous fields
      lfluc = .false.
!     output of uz
      luz = .true.
!     output of uth
      luth =.false.       
!     output of ur
      lur = .false.
!     nondimensionalized by up(0),ub(1),utau(2)
      lnd = 1
!     it is impossible to average to cell centers
!     e.g. ur on planetz 
      if (lfluc) then
       napltr = 'planetr_fluc' 
       napltz = 'planetz_fluc_'
       naplrz = 'planerz_fluc'
      else
       napltr = 'planetr' 
       napltz = 'planetz_'
       naplrz = 'planerz'
      endif
!
!     extract planes for flow visualization
!
      open(11,file='bou.in',form='formatted')
      read(11,*)
      read(11,*) m1, m2, m3
      read(11,*)
      read(11,*)
      read(11,*)
      read(11,*) alx3d
      read(11,*)
      read(11,*)
      read(11,*)
      read(11,*) reb
      close(11)
!
      m1m=m1-1
      m2m=m2-1
      m3m=m3-1
!
      allocate (th(m1),rc(m2m),rm(m2m),z(m3m))
      allocate (sth(m1),cth(m1))
      allocate (rr(2*m2m))
!
      do i = 1,m1
       th(i) = 2*pi*(i-1.d0)/m1m
       sth(i) = sin(th(i))
       cth(i) = cos(th(i))
      enddo
!
      open(unit=98,file='radcor.out',status='unknown')
      do j=1,m2m
       read(98,*) jj,rc(j),rm(j)
      enddo
      close(98)
!
!     define radial positions for theta-z planes
!
      jplanes(1) = m2m    ! first off-wall cell
      jplanes(2) = m2m-16 ! y^+=15
      do j=1,m2m
       if (rm(j).gt.0.8) exit ! y=0.2
      enddo
      jplanes(3) = j
      do j=1,m2m
       if (rm(j).gt.0.5) exit ! y=0.5
      enddo
      jplanes(4) = j
!
      do j=1,2*m2m
       jj=j-m2m
       rr(j)=rm(jj)
       if (j.le.m2m) then
        jj=1+m2m-j
        rr(j)=-rm(jj)
       endif
      enddo
!
      do k=1,m3m
       z(k)=alx3d*(k-1.d0)/m3m
      enddo
!
      if (lnd.eq.3) then
       open(11,file='flowprop.dat',form='formatted')
       read(11,*) aa,aa,aa,utau,aa,aa,aa,aa,aa,
     .            aa,aa,aa,aa,aa,aa,aa,ttau,
     .            aa,aa,aa,aa,aa,aa,aa,
     .            aa,aa,aa,aa,aa,aa,aa,
     .            aa,pw
       print *,'utau,ttau,pw',utau,ttau,pw
       close(11)
      endif
!     read stats along r
      nv = 30
      allocate(wvarav(m2m,nv))
       if (lfluc) then
       open(11,file='stats_out.dat',action='read')
       do j = 1,m2m
        read(11,*) tmp,tmp,(wvarav(j,n),n=1,nv)
       enddo
       close(11)
      endif
!     set wvarav = 0 for instantaneous fields
      if (.not.lfluc) then
       wvarav = 0. 
      endif
! 
      if (lnd.eq.0) then
       utau = 1.
      elseif (lnd.eq.1) then
       utau = .5
      endif
!
      allocate (planetz(nvars,nplanes,m1m,m3m))
      allocate (planetzl(m1,m3m))
!
      open(30,file='planetz.bin',access='stream',form='unformatted')
      read(30) planetz
      close(30)
!     instantaneous to fluctuation
!     pr and dens not considered
      do lplane = 1,nplanes
       jp = jplanes(lplane) 
       planetz(1,lplane,1:m1m,1:m3m) = planetz(1,lplane,1:m1m,1:m3m) - 
     .                                 wvarav(jp,1)
       planetz(2,lplane,1:m1m,1:m3m) = planetz(2,lplane,1:m1m,1:m3m) -
     .                                 wvarav(jp,2)
       planetz(3,lplane,1:m1m,1:m3m) = planetz(3,lplane,1:m1m,1:m3m) -
     .                                 wvarav(jp,3)
      enddo
!
      do lplane=2,3
!
       write(naplane,1001) lplane
 1001  format(I1.1)
!
       if (luz) then
        planetzl(1:m1m,:)=planetz(3,lplane,1:m1m,:)
        planetzl(m1,:)=planetzl(1,:)
        planetzl=planetzl/utau
        open(11,file=trim(napltz)//naplane//'.q',form='unformatted')
        write(11) m1,1,m3m,1
        write(11) (((planetzl(i,k),i=1,m1),j=1,1),k=1,m3m)
        close(11)
       endif
!     
       if (luth) then 
        planetzl(1:m1m,:)=planetz(1,lplane,1:m1m,:)
        planetzl(m1,:)=planetzl(1,:)
        planetzl=planetzl/utau
        open(11,file=trim(napltz)//naplane//'_uth.q',form='unformatted')
        write(11) m1,1,m3m,1
        write(11) (((planetzl(i,k),i=1,m1),j=1,1),k=1,m3m)
        close(11)
       endif
!
       if (lur) then 
        planetzl(1:m1m,:)=planetz(2,lplane,1:m1m,:)
        planetzl(m1,:)=planetzl(1,:)
        planetzl=planetzl/utau
        open(11,file=trim(napltz)//naplane//'_ur.q',form='unformatted')
        write(11) m1,1,m3m,1
        write(11) (((planetzl(i,k),i=1,m1),j=1,1),k=1,m3m)
        close(11)
       endif
!
       open(11,file='planetz_'//naplane//'.x',form='unformatted')
       write(11) m1,1,m3m
       jj = jplanes(lplane)
       write(11) (((rm(j)*cth(i),i=1,m1),j=jj,jj),k=1,m3m),
     .           (((rm(j)*sth(i),i=1,m1),j=jj,jj),k=1,m3m),
     .           (((z(k),i=1,m1),j=jj,jj),k=1,m3m) 
       close(11)
!
       open(11,file='planetz_unrolled.x',form='unformatted')
       write(11) m1,1,m3m
       write(11) (((th(i),i=1,m1),j=m2m,m2m),k=1,m3m),
     .           (((0.d0,i=1,m1),j=m2m,m2m),k=1,m3m),
     .           (((z(k),i=1,m1),j=m2m,m2m),k=1,m3m)
       close(11)

      enddo
!
      deallocate(planetz)
!
      allocate (planerz(nvars,2*m2m,m3m))
!
      open(30,file='planerz.bin',access='stream',form='unformatted')
      read(30) planerz
      close(30)
!     instantaneous to fluctuation
!     pr and dens not considered
      do j = 1,2*m2m
       jj = j-m2m
       if (j.le.m2m) then
        jj = 1+m2m-j
       endif
       if (j.le.m2m) then
        planerz(1,j,1:m3m) = planerz(1,j,1:m3m) - wvarav(jj,1)
        planerz(2,j,1:m3m) = planerz(2,j,1:m3m) - wvarav(jj,2)
       else
!      inst uth, ur are -q1,-q2 on planerz
!      discontinuous at the axis if the sign of uth/ur unchanged
!      i.e. discontinuous at the axis if using -planerz(1,j,1:m3m)-wvarav(jj,1)
        planerz(1,j,1:m3m) = planerz(1,j,1:m3m) + wvarav(jj,1)
        planerz(2,j,1:m3m) = planerz(2,j,1:m3m) + wvarav(jj,2)
       endif
       planerz(3,j,1:m3m) = planerz(3,j,1:m3m) - wvarav(jj,3)
      enddo
!
      planerz=planerz/utau
!
      if (luz) then
       open(11,file=trim(naplrz)//'.q',form='unformatted')
       write(11) 1,2*m2m,m3m,1
       write(11) (((planerz(3,j,k),i=1,1),j=1,2*m2m),k=1,m3m)
       close(11)
      endif
!
      if (luth) then 
       open(11,file=trim(naplrz)//'_uth.q',form='unformatted')
       write(11) 1,2*m2m,m3m,1
       write(11) (((planerz(1,j,k),i=1,1),j=1,2*m2m),k=1,m3m)
       close(11)
      endif
!
      if (lur) then 
       open(11,file=trim(naplrz)//'_ur.q',form='unformatted')
       write(11) 1,2*m2m,m3m,1
       write(11) (((planerz(2,j,k),i=1,1),j=1,2*m2m),k=1,m3m)
       close(11)
      endif

      open(11,file='planerz.x',form='unformatted')
      write(11) 1,2*m2m,m3m
      write(11) (((rr(j)*cth(i),i=1,1),j=1,2*m2m),k=1,m3m),
     .            (((rr(j)*sth(i),i=1,1),j=1,2*m2m),k=1,m3m),
     .            (((z(k),i=1,1),j=1,2*m2m),k=1,m3m)
      close(11)
!
      deallocate(planerz)
!
      allocate (planetr(m1m,m2m,nvars))
      allocate (planetrl(m1,m2m,nvars))
!
!     open(30,file='planetr.bin',access='stream',form='unformatted')
      open(30,file='planetr.bin',form='unformatted')
      read(30) planetr
      close(30)
!     instantaneous to fluctuation
!     pr and dens not considered
      do j = 1,m2m      
       planetr(1:m1m,j,1) = planetr(1:m1m,j,1) - wvarav(j,1)
       planetr(1:m1m,j,2) = planetr(1:m1m,j,2) - wvarav(j,2)
       planetr(1:m1m,j,3) = planetr(1:m1m,j,3) - wvarav(j,3)
      enddo
!
      planetrl(1:m1m,:,:)=planetr(1:m1m,:,:)
      planetrl(m1,:,:)=planetrl(1,:,:)
      planetrl=planetrl/utau
      if (luz) then
       open(11,file=trim(napltr)//'.q',form='unformatted')
       write(11) m1,m2m,1,1
       write(11) (((planetrl(i,j,3),i=1,m1),j=1,m2m),k=1,1)
       close(11)
      endif
!
      if (luth) then
       open(11,file=trim(napltr)//'_uth.q',form='unformatted')
       write(11) m1,m2m,1,1
       write(11) (((planetrl(i,j,1),i=1,m1),j=1,m2m),k=1,1)
       close(11)
      endif
!
      if (lur) then
       open(11,file=trim(napltr)//'_ur.q',form='unformatted')
       write(11) m1,m2m,1,1
       write(11) (((planetrl(i,j,2),i=1,m1),j=1,m2m),k=1,1)
       close(11)
      endif
! 
      open(11,file='planetr.x',form='unformatted')
      write(11) m1,m2m,1
      write(11) (((rm(j)*cth(i),i=1,m1),j=1,m2m),k=1,1),
     .            (((rm(j)*sth(i),i=1,m1),j=1,m2m),k=1,1),
     .            (((z(k),i=1,m1),j=1,m2m),k=1,1)
      close(11)
!
!    
      deallocate(planetr,planetrl)
! 
      stop
      end
