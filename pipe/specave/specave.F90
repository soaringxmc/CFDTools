    program main
    implicit none
    !---------- basic variables
    integer:: m1m,m2m,m3m,m1,m2,m3,m1mh,m3mh,i,j,k,m,n
    real(kind=8),allocatable,dimension(:,:,:):: spect,specz,spectav,speczav
    real(kind=8),allocatable,dimension(:,:):: spec2d,spec2dav
    real(kind=8),allocatable,dimension(:):: lm1,lm2,lm3
    real(kind=8):: pi,laxis,maxis,tw,utau,ttau,retau,eps=1.0e-30
    !---------- intermediate variables
    integer:: kstart,kend,kstep,kkn
    logical:: alive,lave
    real(kind=8):: tmp1,tmp2,tmp3
    character(4) :: naitav
    character(1) :: navar,naplan
    !--------- end of declaration -------
    
    pi=4.*atan(1.d0)
    
    open(10,file='bou.in',action='read')
    read(10,*)
    read(10,*) m1,m2,m3
    do i=1,3
     read(10,*)
    enddo
    read(10,*) maxis
    close(10)
    
    laxis = 2*pi
    
    print *, 'need time average?'
    read(*,*) lave
    if (lave) then
     print *, 'kstart,kend'
     read(*,*) kstart,kend
    endif
    
    print*,m1,m2,m3
    print*,kstart,kend
    
    m1m=m1-1; m2m=m2-1; m3m=m3-1; m1mh=m1m/2+1; m3mh=m3m/2+1
    allocate(spect(m1mh,m2m,5),specz(m3mh,m2m,5),spectav(m1mh,m2m,5),speczav(m3mh,m2m,5))
    allocate(spec2d(m1mh,m3mh),spec2dav(m1mh,m3mh))
    allocate(lm1(m1m),lm2(m2m),lm3(m3m))
    
    open(13,file='flowprop.dat',action='read')
    read(13,*) tmp1, retau, tmp1, utau, tmp1, tmp1, tmp1, tmp1, tmp1, tmp1, tmp1, tmp1, tmp1, tmp1, tmp1, tmp1, ttau
    close(13)
    
    print *,'utau,ttau',utau,ttau
    if (ttau.eq.0.) ttau=1.
    
    do i=1,m1m
     lm1(i)=laxis/(i-1.d0+eps)*retau !lt+
    enddo
    do i=1,m3m
     lm3(i)=maxis/(i-1.d0+eps)*retau !lz+
    enddo
    
    open(12,file='radcor.out',action='read')
    do i=1,m2m
     read(12,*) tmp1,tmp1,lm2(i) !lm2=rm
    enddo
    close(12)
    
    do i=1,m2m
     lm2(i)=(1.d0-lm2(i))*retau !y+
    enddo
    
    !----- 1D spectra -----------
    
    spectav=0.d0; speczav=0.d0;kkn=0
    
    if (lave) then
     do kstep=kstart,kend
      write(naitav,'(I4.4)') kstep
      inquire(file='specth_'//naitav//'.bin',exist=alive)
      if(alive) then
       print*, kstep
       kkn=kkn+1
       
       !open(21,file='../spec2d'//naplan//'_'//navar//'_'//naitav//'.bin',form='unformatted',action='read')
       !read(21) spec2d !(:,:,k)
       !close(21)
       
       open(20,file='specth_'//naitav//'.bin',form='unformatted',action='read')
       read(20) spect
       close(20)
       
       open(21,file='specze_'//naitav//'.bin',form='unformatted',action='read')
       read(21) specz
       close(21)
       
       speczav=speczav+specz
       spectav=spectav+spect
      endif
     enddo
     tmp1=1.d0/kkn
     speczav=speczav*tmp1; spectav=spectav*tmp1
    else
     open(20, file='specth_ave.bin', form='unformatted')
     read(20) spectav
     close(20)
     open(21, file='specze_ave.bin', form='unformatted')
     read(21) speczav
     close(21)
    endif
    
    !open(22,file='specth-av-tecplot.dat')
    !write(22,*) 'variables = log10(lx+),log10(y+),uth,ur,uz,th,p'
    !write(22,*) 'zone i=',m1mh-2,'j=',m2m-1
    !
    !do j=2,m2m
    !do i=2,m1mh-1
    !write(22,'(20g20.9)') dlog10(lm1(i)),dlog10(lm2(j)),(i-1.d0)*spectav(i,j,1:3)/utau**2,&
    !(i-1.d0)*spectav(i,j,4)/ttau**2,(i-1.d0)*spectav(i,j,5)/utau**4
    !enddo
    !enddo
    !close(22)
    !
    !open(22,file='specz-av-tecplot.dat')
    !write(22,*) 'variables = log10(lz+),log10(y+),uth,ur,uz,th,p'
    !write(22,*) 'zone i=',m3mh-2,'j=',m2m-1
    !
    !do j=2,m2m
    !do i=2,m3mh-1
    !write(22,'(20g20.9)') dlog10(lm3(i)),dlog10(lm2(j)),(i-1.d0)*speczav(i,j,1:3)/utau**2,&
    !(i-1.d0)*speczav(i,j,4)/ttau**2,(i-1.d0)*speczav(i,j,5)/utau**4
    !enddo
    !enddo
    !close(22)
    
    open(23,file='specth-av-plot3d.xyz',form='unformatted')
    write(23) m1mh-2,m2m-1
    write(23) ((lm1(i),i=2,m1mh-1),j=2,m2m),((lm2(j),i=2,m1mh-1),j=2,m2m)
    close(23)
    
    open(23,file='specth-av-plot3d.q',form='unformatted')
    write(23) m1mh-2,m2m-1,5
    
    write(23) (((i-1.d0)*spectav(i,j,1)/utau**2,i=2,m1mh-1),j=2,m2m),&
            (((i-1.d0)*spectav(i,j,2)/utau**2,i=2,m1mh-1),j=2,m2m),&
            (((i-1.d0)*spectav(i,j,3)/utau**2,i=2,m1mh-1),j=2,m2m),&
            (((i-1.d0)*spectav(i,j,4)/ttau**2,i=2,m1mh-1),j=2,m2m),&
            (((i-1.d0)*spectav(i,j,5)/utau**4,i=2,m1mh-1),j=2,m2m)
    
    close(23)
    
    open(23,file='specz-av-plot3d.xyz',form='unformatted')
    write(23) m3mh-2,m2m-1
    write(23) ((lm3(i),i=2,m3mh-1),j=2,m2m),((lm2(j),i=2,m3mh-1),j=2,m2m) 
    close(23)
    
    open(23,file='specz-av-plot3d.q',form='unformatted')
    write(23) m3mh-2,m2m-1,5
    
    write(23) (((i-1.d0)*speczav(i,j,1)/utau**2,i=2,m3mh-1),j=2,m2m), &
            (((i-1.d0)*speczav(i,j,2)/utau**2,i=2,m3mh-1),j=2,m2m), &
            (((i-1.d0)*speczav(i,j,3)/utau**2,i=2,m3mh-1),j=2,m2m), &
            (((i-1.d0)*speczav(i,j,4)/ttau**2,i=2,m3mh-1),j=2,m2m), &
            (((i-1.d0)*speczav(i,j,5)/utau**4,i=2,m3mh-1),j=2,m2m)
    
    close(23)
    
    open(23,file='spec-av-plot3d.nam')
    write(23,*) 'uth'
    write(23,*) 'ur'
    write(23,*) 'uz'
    write(23,*) 'th'
    write(23,*) 'p'
    close(23)
    
    
    !----------- 2D spectra
    
    end
