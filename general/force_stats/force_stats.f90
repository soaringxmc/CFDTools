!  ForceMeanRms.f90 
!
!  FUNCTIONS:
!  ForceMeanRms - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: ForceMeanRms
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program ForceMeanRms

    use Var
    implicit none
    
    integer                 ::  i
    real                    ::  sinA, cosA, deno, a1, a2
    real    ,dimension(3)   ::  fp, fv, ft
    
    open ( 11, file = 'input.txt', action = 'read' )
    read ( 11, * )  fileIn
    read ( 11, * )  stepSta
    read ( 11, * )  iSolver
    read ( 11, * )  fileUnst
    read ( 11, * )  filestat
    read ( 11, * ) 
    read ( 11, * )  rhoRef
    read ( 11, * )  vRef
    read ( 11, * )  AoA
    read ( 11, * )  arRef
    close( 11    )  

    sinA    =   sin( AoA / 180.0 * 3.141592653 )
    cosA    =   cos( AoA / 180.0 * 3.141592653 )
    deno    =   0.5*vRef**2*arRef
    
    call InputData 
    
    if( allocated( lift    ) ) deallocate( lift              )
    if( allocated( drag    ) ) deallocate( drag              )
    if( allocated( mome    ) ) deallocate( mome              )
    if( allocated( liftMov ) ) deallocate( liftMov           )
    if( allocated( dragMov ) ) deallocate( dragMov           )
    if( allocated( momeMov ) ) deallocate( momeMov           )
                                allocate  ( lift   ( 3, nStep ) );  lift    =   0.0
                                allocate  ( drag   ( 3, nStep ) );  drag    =   0.0
                                allocate  ( mome   ( 3, nStep ) );  mome    =   0.0
                                allocate  ( liftMov( 3, nStep ) );  liftMov =   0.0
                                allocate  ( dragMov( 3, nStep ) );  dragMov =   0.0
                                allocate  ( momeMov( 3, nStep ) );  momeMov =   0.0
                                                                    liftRMS =   0.0
                                                                    dragRMS =   0.0
                                                                    momeRMS =   0.0
    
    do  i   =   1, nstep
        
        if( iSolver == 1 ) then
            
            fp  ( 1:3    )  =   force( 1:3, i ) / deno
            fv  ( 1:3    )  =   force( 4:6, i ) / deno
            ft  ( 1:3    )  =   force( 7:9, i ) / deno
            lift(   1, i )  =   fp( 1 )*( -sinA ) + fp( 2 )*( cosA )
            lift(   2, i )  =   fv( 1 )*( -sinA ) + fv( 2 )*( cosA )
            lift(   3, i )  =   ft( 1 )*( -sinA ) + ft( 2 )*( cosA )
            drag(   1, i )  =   fp( 1 )*(  cosA ) + fp( 2 )*( sinA )
            drag(   2, i )  =   fv( 1 )*(  cosA ) + fv( 2 )*( sinA )
            drag(   3, i )  =   ft( 1 )*(  cosA ) + ft( 2 )*( sinA )
            
        elseif( iSolver == 2 ) then
            
            lift(   1, i )  =   0.0
            lift(   2, i )  =   0.0
            lift(   3, i )  =   force( 1, i )
            drag(   1, i )  =   force( 3, i )
            drag(   2, i )  =   force( 4, i )
            drag(   3, i )  =   force( 2, i )
            mome(   1, i )  =   force( 6, i )
            mome(   2, i )  =   force( 7, i )
            mome(   3, i )  =   force( 8, i )
            
        endif
        
        
        if( i >= nsta ) then
            
            a1  =   max( real( i-nSta ), 0.0 )
            a2  =   1.0 / ( a1 + 1.0 )
            
            liftMov( 1:3, i )   =   ( liftMov( 1:3, i-1 ) * a1 + lift( 1:3, i ) ) * a2
            dragMov( 1:3, i )   =   ( dragMov( 1:3, i-1 ) * a1 + drag( 1:3, i ) ) * a2
            momeMov( 1:3, i )   =   ( momeMov( 1:3, i-1 ) * a1 + mome( 1:3, i ) ) * a2
            
        endif

    enddo
    
        
    call OutputUnst
        
    do  i   =   nSta, nStep
        
        liftRMS( 1:3 )  =   liftRMS( 1:3 ) + ( lift( 1:3, i ) - liftMov( 1:3, nStep ) )**2
        dragRMS( 1:3 )  =   dragRMS( 1:3 ) + ( drag( 1:3, i ) - dragMov( 1:3, nStep ) )**2
        momeRMS( 1:3 )  =   momeRMS( 1:3 ) + ( mome( 1:3, i ) - momeMov( 1:3, nStep ) )**2
    enddo
    
        liftRMS( 1:3 )  =   sqrt( liftRMS(1:3) / ( nstep - nsta + 1 ) )
        dragRMS( 1:3 )  =   sqrt( dragRMS(1:3) / ( nstep - nsta + 1 ) )
        momeRMS( 1:3 )  =   sqrt( momeRMS(1:3) / ( nstep - nsta + 1 ) )
        
    call OutputStats
    

    end program ForceMeanRms
    
    
    
    subroutine InputData
    
    use Var
    implicit none
    
    integer                 ::  i,itmp
    
    open( 12, file=filein, status='old', action='read' )
    
    read( 12, * )
    if( iSolver == 1 ) then
        read( 12, * )
        read( 12, * )
    endif
    
    if( allocated( force  ) )   deallocate( force                 )
    if( allocated( stepNo ) )   deallocate( stepNo                )
                                allocate  ( force ( 9, 100000000 ) )
                                allocate  ( stepNo(    100000000 ) )
    
    i   =   0
    do while( .not.eof(12) )  
        
        i   =   i + 1
        if( iSolver == 1 ) then
            read( 12, * )   stepNo(i), itmp, itmp, force( 1:9, i )
        elseif( iSolver == 2 ) then
            read( 12, * )   stepNo(i), itmp, force( 1:9, i )
        endif
        
        if( i > 1 ) then
            if( stepNo( i ) <= stepNo( i-1 ) )  then
                i   =   i - 1
                cycle
            endif
        endif
        if( stepNo( i ) == stepSta )    nSta = i
        

    enddo
    
    nStep=  i
    
    end subroutine InputData
    
    
    
    subroutine OutputUnst
    
    use Var
    implicit none
    
    integer ::  i
    
    open(  13, file = fileUnst, action = 'write' )
    write( 13, '(a)' )  'Variables = "step", "CLp", "CLv", "CL", "CDp", "CDV", "CD", "Mx", "My", "Mz", "CLmov", "CDmov", "MxMov", "MyMov", "MzMov", "CDpMov" ' 
    write( 13, '(i10,15e20.10)' ) ( stepNo(i), lift   ( 1:3, i ), drag   ( 1:3, i ), mome    ( 1:3, i ), &
                                              liftMov( 3  , i ), dragMov( 3  , i ), momeMov ( 1:3, i ), dragMov( 1  , i ), i = 1, nStep )
    close( 13 )
    
    end subroutine OutputUnst

    
    
    subroutine OutputStats
    
    use Var
    implicit none
    
    open ( 14, file = fileStat, action = 'write' )
    write( 14, '(a10,i20)'     )  'stepSta  ', stepNo(     nSta  )
    write( 14, '(a10,i20)'     )  'stepEnd  ', stepNo(     nStep )
    write( 14, '(a10,i20)'     )  'stepStats', stepNo(     nStep ) - stepNo( nSta ) + ( stepNo( nStep ) - stepNo( nStep-1 ) )
    write( 14, '(a10,2e20.10)' )  'CLp      ', liftMov( 1, nStep ), liftRMS( 1 )
    write( 14, '(a10,2e20.10)' )  'CLv      ', liftMov( 2, nStep ), liftRMS( 2 )
    write( 14, '(a10,2e20.10)' )  'CL       ', liftMov( 3, nStep ), liftRMS( 3 )
    write( 14, '(a10,2e20.10)' )  'CDp      ', dragMov( 1, nStep ), dragRMS( 1 )
    write( 14, '(a10,2e20.10)' )  'CDv      ', dragMov( 2, nStep ), dragRMS( 2 )
    write( 14, '(a10,2e20.10)' )  'CD       ', dragMov( 3, nStep ), dragRMS( 3 )
    write( 14, '(a10,2e20.10)' )  'Mx       ', momeMov( 1, nStep ), momeRMS( 1 )
    write( 14, '(a10,2e20.10)' )  'My       ', momeMov( 2, nStep ), momeRMS( 2 )
    write( 14, '(a10,2e20.10)' )  'Mz       ', momeMov( 3, nStep ), momeRMS( 3 )
    close( 14                )
    
    end subroutine OutputStats
    

    
    


