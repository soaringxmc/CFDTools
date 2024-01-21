!  Correlation.f90 
!
!  FUNCTIONS:
!  Correlation - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Correlation
!
!  PURPOSE: Spanwise Correlation 
!           data output by charles
!
!****************************************************************************

    program Correlation

    implicit none
    
    ! Variables
    character(len=1000) :: filename_in,filename_out
    integer(kind=4)     :: lsta,lend
    integer(kind=4)     :: np, nts, its, np_half
    integer(kind=4)     :: i, j, jf
    real(kind=8)        :: real_tmp, rtmp1, rtmp2, ave_sp
    real(kind=8),allocatable    :: sp(:,:)
    real(kind=8),allocatable    :: coef(:)
    integer(kind=4),allocatable :: index_p(:)
    

    
    
    ! Body of Correlation
    ! open(11,file = 'input.txt',action = 'read', status = 'old')
    ! open(12,file = 'output.txt', action = 'read', status = 'old')
    ! read(12,*)
    
    open(11,file = 'input.dat',action = 'read', status = 'old')
    read(11,*)
    read(11,*) filename_in
    read(11,*)
    read(11,*) lsta
    read(11,*) lend
    read(11,*)
    read(11,*) filename_out
    close(11)
    
    ! number of time steps
    nts = lend - lsta + 1
    
    open(12, file = filename_in, action = 'read', status = 'old')
    
    ! number of points
    read(12,*) real_tmp, real_tmp, np
    backspace(12)
    
    ! read sample from file
    allocate(sp(np,nts))

    do its = 1, nts
        read(12,*) real_tmp, real_tmp, real_tmp, sp(1:np,its)
    end do
    
    close(12)
    
    ! average in time direction and spanwise dirction
    ave_sp  = sum(sum(sp,1))/(nts*np)
    
    ! instanteneous -> fluctuation
    sp(:,:) = sp(:,:) - ave_sp
    
    ! maximum spanwise length need to be computed
    ! when points number is odd (such as 5), interval number is odd (5), 2 intervals is the maximum spanwise interval needed to be computed
    ! when points number is even (such as 6), interval number is even (6), 3 intervals is the maximum spanwise interval needed to be conputed
    
    np_half = np/2
    
    ! array for find a certain point index
    allocate(index_p(-(np-1):2*np))
    
    do i = 1, np
        index_p(i) = i
    end do
    do i = np+1, np*2
        index_p(i) = i - np
    end do
    do i = -(np-1), 0
        index_p(i) = i + np
    end do
    
    
    !define coef
    allocate(coef(0:np_half))             ! +1 because 0 interval 
    coef    = 0.0
    coef(0) = 1.0
    
    
    ! calculate correlation for interval (0~np_half)
    ! assign correlation = 1 for interval 0
    print *, 'interval', 0,'/', np, 'is done'
    do i = 1, np_half
        print *, 'interval', i, '/', np, 'is done'
        do j = 1, np
            
            jf = index_p(j + i)
            
            real_tmp    =   dot_product(sp(jf,1:nts),sp(j,1:nts))
            rtmp1       =   dot_product(sp(j,1:nts),sp(j,1:nts))
            rtmp2       =   dot_product(sp(jf,1:nts),sp(jf,1:nts))
            !coef(i)     =   coef(i) + real_tmp/dot_product(sp(j,1:nts),sp(j,1:nts))
            !coef(i)     =   coef(i) + real_tmp/dot_product(sp(jf,1:nts),sp(jf,1:nts))
            coef(i)     =   coef(i) + real_tmp/sqrt(rtmp1*rtmp2)
            
        enddo
    
    enddo
    
    !coef(1:np_half)    = coef(1:np_half)/(2*np)
    coef(1:np_half)    = coef(1:np_half)/(np)
    
    
    open(13, file = filename_out, action = 'write', status = 'replace')
    
    do i = 0, np_half
        write(13,'(2e20.10)') dble(i)/np, coef(i)
    end do
    
    close(13)
    
    
    end program Correlation

