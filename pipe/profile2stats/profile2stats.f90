!  profile2stats.f90 
!
!  FUNCTIONS:
!  profile2stats - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: profile2stats
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program profile2stats
    character ctmp
    real(8),allocatable :: wmean1(:,:),wmean2(:,:)
!
    print *, 'reb,retau'
    read(*,*) reb,retau
!
    open(11, file='profile.dat')
    do 
     read(11,*) ctmp
     if(ctmp.ne.'%') exit
    enddo
    backspace(11)
    nv1 = 8
    nv2 = 32
    allocate(wmean1(10000,nv1))
    allocate(wmean2(10000,nv2))
    wmean1 = 0.
    i = 0
    do while (.not.eof(11))
     i = i + 1
     read(11,*) wmean1(i,1:8)
    enddo
    nli = i
    print *, nli
    close(11)
!   wmean1(i,1), y^+
!   wmean1(i,2), U_z^+
!   wmean1(i,3), P^+
!   wmean1(i,4), u_z^2
!   wmean1(i,5), u_r^2
!   wmean1(i,6), u_t^2
!   wmean1(i,7), u_r u_z
!   wmean1(i,8), p^2
!
!   outer-scaled
!   all variables located at cell centers
!
    wmean2 = 0.
    do i = 1, nli
     wmean2(i,2) = 1.- wmean1(i,1)/retau !rm/R
     wmean2(i,5) = wmean1(i,2)*retau/reb !Uz/Up
     wmean2(i,6) = sqrt(wmean1(i,6))*retau/reb !(utut)^(1/2)/Up
     wmean2(i,7) = sqrt(wmean1(i,5))*retau/reb !(urur)^(1/2)/Up
     wmean2(i,8) = sqrt(wmean1(i,4))*retau/reb !(uzuz)^(1/2)/Up
     wmean2(i,9) = wmean1(i,7)*(retau/reb)**2  !(uruz)/(UpUp)
     if (wmean2(i,9).ne.wmean2(i,9)) then
      wmean2(i,9) = 0.
     endif
    enddo
    open(11, file='stats_out.dat', form='formatted')
    do i = 1, nli
     write(11,'(50E20.10)') wmean2(i,1:nv2)
    enddo
    close(11)
!
!   inner-scaled
!   all variables located at cell centers
!
    wmean2 = 0.
    do i = 1, nli
     wmean2(i,2) = wmean1(i,1) !ym+
     wmean2(i,5) = wmean1(i,2) !Uz+
     wmean2(i,6) = sqrt(wmean1(i,6)) !(utut)^(1/2)+
     wmean2(i,7) = sqrt(wmean1(i,5)) !(urur)^(1/2)+
     wmean2(i,8) = sqrt(wmean1(i,4)) !(uzuz)^(1/2)+
     wmean2(i,9) = wmean1(i,7)  !(uruz)+
     if (wmean2(i,9).ne.wmean2(i,9)) then
      wmean2(i,9) = 0.
     endif
    enddo
    open(11, file='stats_inn.dat', form='formatted')
    do i = 1, nli
     write(11,'(50E20.10)') wmean2(i,1:nv2)
    enddo
    close(11)
!
    end program profile2stats

