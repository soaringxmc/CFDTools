!  plane_utau2ub.f90 
!
!  FUNCTIONS:
!  plane_utau2ub - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: plane_utau2ub
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program plane_utau2ub

    real(8),allocatable :: uz(:,:,:)
    character(1) :: nap
    
    utau = 2.28E-002 !/diskg/PIPE/RETAU3000/RUN/flowprop.dat
    uclp = 26.417089442982789 !/diskg/PIPE/RETAU3000/RUN/flowprop.dat
    ucl = uclp*utau

    open(11, file='planetr.q', form='unformatted')
    read(11) m1,m2m,itmp,itmp
    allocate(uz(m1,m2m,1))
    read(11) (((uz(i,j,k), i=1,m1), j=1,m2m), k=1,1)
    close(11)
    uz = uz*utau*2
    open(11, file='planetr_1.q', form='unformatted')
    write(11) m1,m2m,1,1
    write(11) (((uz(i,j,k), i=1,m1), j=1,m2m), k=1,1)
    deallocate(uz)
    close(11)
!
    open(11, file='planerz.q', form='unformatted')
    read(11) itmp,m2m,m3m,itmp
    m2m = m2m/2 !real m2m
    allocate(uz(1,2*m2m,m3m))
    read(11) (((uz(i,j,k), i=1,1), j=1,2*m2m), k=1,m3m)
    close(11)
!   uz = uz*utau*2
    uz = uz*2.
    open(11, file='planerz_1.q', form='unformatted')
    write(11) 1,2*m2m,m3m,1
    write(11) (((uz(i,j,k), i=1,1), j=1,2*m2m), k=1,m3m)
    deallocate(uz)
    close(11)
!
    do ip = 2, 3
     write(nap,'(I1)') ip
     open(11, file='planetz_'//nap//'.q', form='unformatted')
     read(11) m1,itmp,m3m,itmp
     allocate(uz(m1,1,m3m))
     read(11) (((uz(i,j,k), i=1,m1), j=1,1), k=1,m3m)
     close(11)
     uz = uz*utau*2
     open(11, file='planetz_'//nap//'_1.q', form='unformatted')
     write(11) m1,1,m3m,1
     write(11) (((uz(i,j,k), i=1,m1), j=1,1), k=1,m3m)
     deallocate(uz)
     close(11)
    enddo
    stop
    end program plane_utau2ub

