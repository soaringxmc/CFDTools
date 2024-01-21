    program main
    implicit none
    !---------- basic variables
    integer:: m1m,m2m,m3m,m1,m2,m3,m1mh,m3mh,i,j,k,m,n
    real(kind=8),allocatable,dimension(:,:,:):: spect,specz,spectav,speczav
    real(kind=8),allocatable,dimension(:):: lm1,lm2,lm3
    real(kind=8):: pi,laxis,maxis,tw,utau,ttau,retau,eps=1.0e-30
    !---------- intermediate variables
    integer:: kstart,kend,kstep,kkn
    logical:: alive
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
    read(10,*) maxis    !axial length, maxis = 15 in general
    close(10)
    
    laxis=2*pi
    
    open(11,file='spec.in',action='read')
    read(11,*) kstart,kend
    close(11)
    
    print*,m1,m2,m3
    print*,kstart,kend
    
    m1m=m1-1; m2m=m2-1; m3m=m3-1; m1mh=m1m/2+1; m3mh=m3m/2+1
    allocate(spect(m1mh,m2m,5),specz(m3mh,m2m,5),spectav(m1mh,m2m,5),speczav(m3mh,m2m,5))
    allocate(lm1(m1m),lm2(m2m),lm3(m3m))
    
    open(13,file='flowprop.dat',action='read')
    read(13,*) tmp1, retau, tmp1, utau, tmp1, tmp1, tmp1, tmp1, tmp1, tmp1, tmp1, tmp1, tmp1, tmp1, tmp1, tmp1, ttau
    close(13)
    
    print *,'utau,ttau',utau,ttau
    if (ttau.eq.0.) ttau=1.
    
    do i=1,m3m
     lm3(i)=maxis/(i-1.d0+eps)*retau  !axial length, maxis = 15 in general
    enddo
    do i=1,m1m
     lm1(i)=laxis/(i-1.d0+eps)*retau  !laxis=2*pi, so lm1 assumes r=R
    enddo
    
    open(12,file='radcor.out',action='read')
    do i=1,m2m
     read(12,*) tmp1,tmp1,lm2(i)
    enddo
    close(12)
    
    do i=1,m2m
     lm2(i)=(1.d0-lm2(i))*retau
    enddo
    
    !----- 1D spectra -----------
    
    spectav=0.d0; speczav=0.d0;kkn=0
    
    open(23,file='specth-av-plot3d.xyz',form='unformatted',action='read')
    read(23) i,i
    read(23) ((lm1(i),i=2,m1mh-1),j=2,m2m),((lm2(j),i=2,m1mh-1),j=2,m2m)
    close(23)
    
    open(23,file='specth-av-plot3d.q',form='unformatted',action='read')
    read(23) i,i,i
    read(23) (((spectav(i,j,k),i=2,m1mh-1),j=2,m2m),k=1,5)
    close(23)
    
    open(23,file='specz-av-plot3d.xyz',form='unformatted',action='read')
    read(23) i,i
    read(23) ((lm3(i),i=2,m3mh-1),j=2,m2m),((lm2(j),i=2,m3mh-1),j=2,m2m) 
    close(23)
    
    open(23,file='specz-av-plot3d.q',form='unformatted',action='read')
    read(23) i,i,i
    read(23) (((speczav(i,j,k),i=2,m3mh-1),j=2,m2m),k=1,5)
    close(23)
    
    open(22,file='specz.dat')
    write(22,*) 'variables = lz+,y+,w,v,u,T,p'
    write(22,*) 'zone i=',m1mh-2,'j=',m2m-1
    
    do j=2,m2m
    do i=2,m1mh-1
     write(22,'(20g20.9)') lm1(i),lm2(j),spectav(i,j,1:5)
    enddo
    enddo
    close(22)
    
    !rewind(77)
    !j=m2m-15
    !do i=2,m1mh-1
    ! write(77,*) lm1(i),spectav(i,j,1:5)
    !enddo 
    !close(77)
    
    open(22,file='specx.dat')
    write(22,*) 'variables = lx+,y+,w,v,u,T,p'
    write(22,*) 'zone i=',m3mh-2,'j=',m2m-1
    
    do j=2,m2m
    do i=2,m3mh-1
     write(22,'(20g20.9)') lm3(i),lm2(j),speczav(i,j,1:5)
    enddo
    enddo
    close(22)
    
    do k=1,5
    do j=2,m2m
     tmp1=0.d0
     do i=2,m1mh-1
      spect(i,j,k)=spectav(i,j,k)/(i-1.d0)
      tmp1=tmp1+spect(i,j,k)
     enddo
!    tmp1 is variance, spect denotes theta direction
!    spectav is premultiplied azimulthal spectral densities
!    spect is premultiplied normalized azimulthal spectral densities
     do i=2,m1mh-1
      spect(i,j,k)=spect(i,j,k)/tmp1
      spect(i,j,k)=(i-1.d0)*spect(i,j,k)
     enddo
    enddo
    enddo
    
    do k=1,5
    do j=2,m2m
     tmp1=0.d0
     do i=2,m3mh-1
      specz(i,j,k)=speczav(i,j,k)/(i-1.d0)
      tmp1=tmp1+specz(i,j,k)
     enddo
!    tmp1 is variance, specz denotes axial direction
!    speczav is premultiplied axial spectral densities
!    specz is premultiplied normalized axial spectral densities
     do i=2,m3mh-1
      specz(i,j,k)=specz(i,j,k)/tmp1
      specz(i,j,k)=(i-1.d0)*specz(i,j,k)
     enddo
    enddo
    enddo
    
    !azimulthal direction, specz_norm.dat, spect
    !w,v,u are azimulthal, radial and axial components, respectively
    open(22,file='specz_norm.dat') 
    write(22,*) 'variables = lz+,y+,w,v,u,T,p'
    write(22,*) 'zone i=',m1mh-2,'j=',m2m-1
    
    do j=2,m2m
    do i=2,m1mh-1
     write(22,'(20g20.9)') lm1(i),lm2(j),spect(i,j,1:5)
    enddo
    enddo
    close(22)
    
    !axial direction specx_norm.dat, specz
    open(22,file='specx_norm.dat') 
    write(22,*) 'variables = lx+,y+,w,v,u,T,p'
    write(22,*) 'zone i=',m3mh-2,'j=',m2m-1
    
    do j=2,m2m
    do i=2,m3mh-1
     write(22,'(20g20.9)') lm3(i),lm2(j),specz(i,j,1:5)
    enddo
    enddo
    close(22)
    
    
    !----------- 2D spectra
    
    end
