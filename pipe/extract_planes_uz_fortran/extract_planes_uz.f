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
      integer :: jplanes(nplanes)
      character(1) :: naplane
!
!     Extract planes for flow visualization
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
!     print *,'m1m=',m1m
!     print *,'m2m=',m2m
!     print *,'m3m=',m3m
!
      allocate (th(m1),rc(m2m),rm(m2m),z(m3m))
      allocate (sth(m1),cth(m1))
      allocate (rr(2*m2m))
!
      do i=1,m1
       th(i)=2*pi*(i-1.d0)/m1m
       sth(i)=sin(th(i))
       cth(i)=cos(th(i))
      enddo
!
!     print *,sth
!     pause
!     print *,cth
!     pause
!
      open(unit=98,file='radcor.out',status='unknown')
      do j=1,m2m
!      print *,'j=',j
       read(98,*) jj,rc(j),rm(j)
      enddo
      close(98)
!
!     Define radial positions for theta-z planes
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
       if (j.le.m2m) then
        jj=1+m2m-j
        rr(j)=-rm(jj)
       else
        jj=j-m2m
        rr(j)=rm(jj)
       endif
      enddo
!
      do k=1,m3m
       z(k)=alx3d*(k-1.d0)/m3m
      enddo
!     print *,z
!     pause
!
!     print *,'Input utau'
!     read(*,*) utau
!     utau = 2.1130322655568006E-002
!
      open(11,file='flowprop.dat',form='formatted')
      read(11,*) aa,aa,aa,utau,aa,aa,aa,aa,aa,
     .           aa,aa,aa,aa,aa,aa,aa,ttau,
     .           aa,aa,aa,aa,aa,aa,aa,
     .           aa,aa,aa,aa,aa,aa,aa,
     .           aa,pw
      print *,'utau,ttau,pw',utau,ttau,pw
      close(11)
! 
      utau = 0.5
      retau=reb*utau
!     print *,'Retau=',retau
!
      allocate (planetz(nvars,nplanes,m1m,m3m))
      allocate (planetzl(m1,m3m))
!
      open(30,file='planetz.bin',access='stream',form='unformatted')
      read(30) planetz
      close(30)
!
!     lplane = 3 ! index of shell to be extracted
!
      do lplane=1,3
!
       write(naplane,1001) lplane
 1001  format(I1.1)
!
       planetzl(1:m1m,:)=planetz(3,lplane,1:m1m,:)
       planetzl(m1,:)=planetzl(1,:)
!
       planetzl=planetzl/utau
!
       open(11,file='planetz_'//naplane//'.q',form='unformatted')
       rewind (11)
       write  (11) m1,1,m3m,1
       write  (11) (((planetzl(i,k),i=1,m1),j=1,1),k=1,m3m)
       close(11)
! 
       open(11,file='planetz_'//naplane//'.x',form='unformatted')
       rewind (11)
       write  (11) m1,1,m3m
!      write  (11) (((rm(j)*cth(i),i=1,m1),j=m2m,m2m),k=1,m3m),
!    .             (((rm(j)*sth(i),i=1,m1),j=m2m,m2m),k=1,m3m),
!    .             (((z(k),i=1,m1),j=m2m,m2m),k=1,m3m) 
       jj = jplanes(lplane)
       write  (11) (((rm(j)*cth(i),i=1,m1),j=jj,jj),k=1,m3m),
     .             (((rm(j)*sth(i),i=1,m1),j=jj,jj),k=1,m3m),
     .             (((z(k),i=1,m1),j=jj,jj),k=1,m3m)
       close(11)
!
       open(11,file='planetz_unrolled.x',form='unformatted')
       rewind (11)
       write  (11) m1,1,m3m
       write  (11) (((th(i),i=1,m1),j=m2m,m2m),k=1,m3m),
     .             (((0.,i=1,m1),j=m2m,m2m),k=1,m3m),
     .            (((z(k),i=1,m1),j=m2m,m2m),k=1,m3m) 
       close(11)
!
      enddo
!
      deallocate(planetz)
!
      allocate (planerz(nvars,2*m2m,m3m))
!
      open(30,file='planerz.bin',access='stream',form='unformatted')
      read(30) planerz
      write  (100,*) (((planerz(1,j,k),i=1,1),j=1,2*m2m),k=1,m3m)
      close(30)
!
      planerz=planerz/utau
!
      open(11,file='planerz.q',form='unformatted')
      rewind (11)
      write  (11) 1,2*m2m,m3m,1
      write  (11) (((planerz(3,j,k),i=1,1),j=1,2*m2m),k=1,m3m)
      close(11)
! 
      open(11,file='planerz.x',form='unformatted')
      rewind (11)
      write  (11) 1,2*m2m,m3m
      write  (11) (((rr(j)*cth(i),i=1,1),j=1,2*m2m),k=1,m3m),
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
!
!     print *,minval(planetr(:,:,3))
!     print *,maxval(planetr(:,:,3))
!
      planetrl(1:m1m,:,:)=planetr(1:m1m,:,:)
      planetrl(m1,:,:)=planetrl(1,:,:)
!
      planetrl=planetrl/utau
!
      open(11,file='planetr.q',form='unformatted')
      rewind (11)
      write  (11) m1,m2m,1,1
      write  (11) (((planetrl(i,j,3),i=1,m1),j=1,m2m),k=1,1)
      close(11)
! 
      open(11,file='planetr.x',form='unformatted')
      rewind (11)
      write  (11) m1,m2m,1
      write  (11) (((rm(j)*cth(i),i=1,m1),j=1,m2m),k=1,1),
     .            (((rm(j)*sth(i),i=1,m1),j=1,m2m),k=1,1),
     .            (((z(k),i=1,m1),j=1,m2m),k=1,1)
      close(11)
!
!    
      deallocate(planetr,planetrl)
! 
      end
