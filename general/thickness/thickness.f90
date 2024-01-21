!  Vorticity_Thickness.f90 
!
!  FUNCTIONS:
!  Vorticity_Thickness - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Vorticity_Thickness
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program Vorticity_Thickness

    implicit none
    
    ! Variables
    integer             i,nLoca,j,idu_max,iu_max,iu_min
    real                u_min,u_max,du_max,v(100,10000),u1,u0,d1,d0,d_omg,u00,loca(1000),theta,spacing,u_dif
    character(len=4)    buf
    character(len=1000) file_out

    ! Body of Vorticity_Thickness
    open( 11, file = 'vort_thickness.inp' )
    read( 11, * ) file_out
    read( 11, * ) nLoca
    do j = 1, nLoca
        read( 11, * ) loca(j)    
    enddo
    close(11)
    
    
    open( 31, file = file_out )
    write(31, '(a)') 'variables = "x", "u_min", "u_max", "du_max", "idu_max", "d_omg", "theta", "ratio", "iu_min", "iu_max"'

    do j = 1, nLoca
        
        write(buf,'(f4.2)') loca(j)
        !open( 21, file = 'E:\04_IceWing\Straight\02_SLA\01_line_'//buf//'_perp.dat' )
        open( 21, file = '.\01_line_'//buf//'.dat' )
        !open( 21, file = 'E:\04_IceWing\Straight\02_SLA\01_line_0.10.dat' )
        
        do i = 1, 55
            read( 21, * )
        enddo
    
        u_min   =  100.0
        u_max   = -100.0
        du_max  = -100.0
        idu_max = -100
        do i = 1, 10000000
            
            if( eof(21) ) exit
            
            read( 21, * ) v( 1:53, i )
            
            if( i > 1 ) then   !estimated via streamwise component
                u1 = v(5, i  )
                d1 = v(53,i  )
                u0 = v(5, i-1)
                d0 = v(53,i-1)
                
                if( u1 < u_min ) then
                    u_min  = u1
                    iu_min = i
                elseif( u1 > u_max ) then
                    u_max  = u1
                    iu_max = i
                endif
                
                if( i > 2 ) then
                    u00 = v(5, i-2 )
                    if( (u1-u0) < 2.0*(u0-u00) ) then           !rule out patches of blocks
                        
                        if( (u1-u0)/(d1-d0) > du_max ) then
                            du_max  = (u1-u0)/(d1-d0)
                            idu_max = i
                        endif
                        
                    endif
                endif
                !du_max  =   max( du_max, (u1-u0)/(d1-d0) )
            endif
            
            !no difference between vertical lines (u component) and perpendicular lines (uv components)
            
            !if( i > 1 ) then    !estimated via 2D components
            !    u1 = sign( sqrt( v(5,i  )**2 + v(6,i  )**2 ), v(5,i  ) )
            !    d1 = v(53,i  )
            !    u0 = sign( sqrt( v(5,i-1)**2 + v(6,i-1)**2 ), v(5,i-1) )
            !    d0 = v(53,i-1)
            !    u_min   =   min( u_min,  u1 )
            !    u_max   =   max( u_max,  u1 )
            !    if( i > 2 ) then
            !        u00 = sign( sqrt( v(5,i-2)**2 + v(6,i-2)**2 ), v(5,i-2) )
            !        if( (u1-u0) < 2.0*(u0-u00) ) then           !rule out patches of blocks
            !            du_max  =   max( du_max, (u1-u0)/(d1-d0) )
            !        endif
            !    endif
            !    !du_max  =   max( du_max, (u1-u0)/(d1-d0) )
            !endif
            
        enddo
        
        d_omg   = (u_max-u_min)/du_max
        
        
        spacing = abs(d1-d0)        !momentum thickness
        theta   = 0.0
        do i = iu_min, iu_max
            
            u1   = v(5, i  )
            d1   = v(53,i  )
            u_dif= (u1-u_min)/(u_max-u_min)
            theta= theta + u_dif*(1.0-u_dif)*spacing
            
        enddo
        
        
        
        close(21)
        
        write(*,'(a,f4.2)') 'Location =', loca(j) 
        !write(*,*) 'umin =', u_min
        !write(*,*) 'umax =', u_max
        !write(*,*) 'dudy_max =', du_max
        !write(*,*) 'd_omg', d_omg
        write(31,'(4f,i,3f,2i)') loca(j),u_min,u_max,du_max,idu_max,d_omg,theta,d_omg/theta,iu_min,iu_max
        
    enddo
    
    close(31)


    end program Vorticity_Thickness
    
    
    !mcr for extracting more x locations
    !what if using lines perpendicular to shear layer, u^2+v^2
    !
    

