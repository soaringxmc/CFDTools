!****************************************************************************
!
!  PROGRAM: Interp
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program Interp
    implicit none

    integer                 ::  nblk_2d,nblk_sp1,nblk_sp2,nblk1,nblk2
    integer                 ::  iblk1,iblk2,nc,n,rt,itmp,idt,ispan1,ispan2,sta1,end1,sta2,end2,ilocal,m
    integer                 ::  itstep,neq,ndt
    real                    ::  tm,pf,rtmp
    character(len=100)      ::  fname
    integer ,   allocatable ::  ni1(:),nj1(:),nk1(:),ni2(:),nj2(:),nk2(:),nsta1(:),nsta2(:),cp(:,:)
    real    ,   allocatable ::  u1(:,:,:),u2(:,:,:)
    
    open ( 11, file = 'input.txt', action = 'read')
    read ( 11, * )  nblk_2d
    read ( 11, * )  nblk_sp1
    read ( 11, * )  nblk_sp2
    
    
    rt      =  nblk_sp2/nblk_sp1
    nblk1   =  nblk_sp1*nblk_2d
    nblk2   =  nblk_sp2*nblk_2d
    
    allocate ( ni1   (nblk1    )  );  ni1    =   0
    allocate ( nj1   (nblk1    )  );  nj1    =   0
    allocate ( nk1   (nblk1    )  );  nk1    =   0
    allocate ( ni2   (nblk2    )  );  ni2    =   0
    allocate ( nj2   (nblk2    )  );  nj2    =   0
    allocate ( nk2   (nblk2    )  );  nk2    =   0
    allocate ( nsta1 (nblk1+1  )  );  nsta1  =   0
    allocate ( nsta2 (nblk2+1  )  );  nsta2  =   0
    allocate ( cp    (nblk2  ,2)  );  cp     =   0
    
    nc      = 0
    nsta1(1)= 1 
    do iblk1    = 1, nblk1
        read( 11, * )   ni1(iblk1),nj1(iblk1),nk1(iblk1)
        itmp            =   (ni1(iblk1)-1)*(nj1(iblk1)-1)*(nk1(iblk1)-1)
        nc              =   nc              + itmp
        nsta1(iblk1+1)  =   nsta1(iblk1)    + itmp
    enddo

    
    allocate ( u1(8,nc,3) );    u1  = 0.0
    allocate ( u2(8,nc,3) );    u2  = 0.0
    
    do  iblk1   = 1, nblk1
        fname   = ''
        call getfilename( iblk1, 22, 'resu_original/d3ns_blk', 4, '.dat', fname )
        open ( 12,  file = fname, action = 'read', form = 'binary' )
        read ( 12 ) itstep, tm                                      !same for all blocks
        read ( 12 ) neq, ndt, pf
        read ( 12 ) ( u1(1:8,n,1    ),   n = nsta1(iblk1),nsta1(iblk1+1)-1              )
        read ( 12 ) ((u1(1:7,n,idt+1),   n = nsta1(iblk1),nsta1(iblk1+1)-1), idt = 1,2  )
        close( 12 )
        write( *,* ) 'iblk1 = ', iblk1
    enddo
    
    do  ispan2              =   1, nblk_sp2
        sta2                =   (ispan2-1)*nblk_2d + 1
        end2                =    ispan2   *nblk_2d
        ispan1              =   (ispan2+rt-1)/rt                !key relation 1
        ilocal              =   mod(ispan2+rt-1,rt) + 1         !key relation 2
        cp(sta2:end2,1)     =   [1:nblk_2d]+nblk_2d*(ispan1-1)  !cp: iblk2 <- iblk1
        cp(sta2:end2,2)     =   ilocal                          !ilocal=1~rt            
    enddo
    
    
    nsta2(1)            =   1
    do  iblk2           =   1, nblk2
        iblk1           =   cp(iblk2,1)
        ilocal          =   cp(iblk2,2)
        
        ni2(iblk2)      =   ni1(iblk1)
        nj2(iblk2)      =   nj1(iblk1)
        nk2(iblk2)      =   (nk1(iblk1)-1)/rt + 1
        
        itmp            =   (ni2(iblk1)-1)*(nj2(iblk1)-1)*(nk2(iblk1)-1)
        sta1            =   nsta1(iblk1)    + itmp*(ilocal-1)
        end1            =   sta1            + (itmp - 1)
        
        itmp            =   (ni2(iblk2)-1)*(nj2(iblk2)-1)*(nk2(iblk2)-1)
        nsta2(iblk2+1)  =   nsta2(iblk2)    + itmp
        
        sta2            =   nsta2(iblk2  )
        end2            =   nsta2(iblk2+1)  - 1
        
        u2(:,sta2:end2,:)   =   u1(:,sta1:end1,:)
        
    enddo
    
    do  iblk2   = 1, nblk2
        fname   = ''
        call getfilename( iblk2, 20, 'resu_interp/d3ns_blk', 4, '.dat', fname )
        open ( 13,  file = fname, action = 'write', form = 'binary' )
        write( 13 ) itstep, tm                                      !same for all blocks
        write( 13 ) neq, ndt, pf
        write( 13 ) ( u2(1:8,n,1  ),     n = nsta2(iblk2),nsta2(iblk2+1)-1              )
        write( 13 ) ((u2(1:7,n,idt+1),   n = nsta2(iblk2),nsta2(iblk2+1)-1), idt = 1,2  )
        close( 13 )
    enddo
    
    do  iblk2   = 1, nblk2
        fname   = ''
        call getfilename( iblk2, 20, 'resu_interp/d3ns_blk', 4, '.dat', fname )
        open ( 13, file = fname, action = 'read', form = 'binary' )
        read( 13 ) itstep, tm                                       !same for all blocks
        read( 13 ) neq, ndt, pf
        read( 13 ) ( (u2(m,n,1  ),     m=1,8), n = nsta2(iblk2),nsta2(iblk2+1)-1              )
        read( 13 ) (((u2(m,n,idt+1),   m=1,7), n = nsta2(iblk2),nsta2(iblk2+1)-1), idt = 1,2  )
        close( 13 )
    enddo
    
    
    end program Interp
    
    


