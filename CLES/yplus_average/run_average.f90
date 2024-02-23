program run_average
    !parameter(np=100000,ne=100000)
    real,allocatable::xyz(:,:),xyz_avg(:,:)
    integer,allocatable::ele(:,:),ncount(:)
    character*80 filename, file,char
    character*80 varname(100)
    
    open(5,file='run_average.txt',status='old',action='read')
    read(5,*)
    read(5,*)np,ne,nep
    read(5,*)
    read(5,*)nvar
    allocate(xyz(np,nvar),xyz_avg(np,nvar),ncount(np),ele(ne,nep))
    xyz_avg=-1.0e30
    ele=0
    ncount=0
    read(5,*)
    read(5,*)file
    read(5,*)
    read(5,*)nstart,nend,ndelta
    close(5)
    
    do n=nstart,nend,ndelta
        xyz=0.0
        write(char,'(i6.6)')n
        filename=trim(file)//"."//trim(char)//".dat"
        write(*,*)filename
        open(5,file=filename,status='old',action='read')
        read(5,*)
        read(5,*)char,char,varname(1:nvar)
        do i=1,5
            read(5,*)
        end do
        do i=1,np
            read(5,*)xyz(i,1:nvar)
        end do
        if(n==nend)then
            do i=1,ne
                read(5,*)ele(i,1:nep)
            end do
        end if
        close(5)
        
        do i=1,np
            logic=0
            do j=1,i
                dis = sqrt( (xyz_avg(j,1)-xyz(i,1))**2 + (xyz_avg(j,2)-xyz(i,2))**2 )
                if(dis<1.e-5)then
                    xyz_avg(j,4:nvar)=xyz_avg(j,4:nvar)*float(ncount(j))/float(ncount(j)+1) + xyz(i,4:nvar)/float(ncount(j)+1)
                    ncount(j)=ncount(j)+1
                    logic=1
                    exit
                endif
            enddo
            
            if(logic==0) then
                do j=1,i
                    if(ncount(j)==0)then
                        xyz_avg(j,1:nvar)=xyz(i,1:nvar)
                        ncount(j)=ncount(j)+1
                        exit
                    endif
                end do
            end if
        end do
    end do
    
    
    !xyz_avg 多个时刻的表面数据进行展向平均
    
        nx0=0
        ncount_0=ncount(1)
        do i=1,np
            if(ncount(i) .ne. 0) then
                if(ncount(i) .ne. ncount_0) then
                    write(*,'(a,3i6,2f20.8)')"Warning:", i,ncount(i),ncount_0,xyz_avg(i,1:2)
                endif
                nx0=nx0+1
            endif
        enddo
        
        open(51,file='test.dat')
        write(51,'(a11,100a10)')'variables = ',varname(1:nvar)
        do i=1,nx0
            write(51,'(100e20.10)')xyz_avg(i,1:nvar)
        end do
        close(51)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=1,np
            do j=1,nx0
                dis = sqrt( (xyz_avg(j,1)-xyz(i,1))**2 + (xyz_avg(j,2)-xyz(i,2))**2 )
                if(dis<1e-6)then
                    xyz(i,4:nvar)=xyz_avg(j,4:nvar)
                    exit
                endif
            end do
        end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        filename=trim(file)//'.avg.dat'
        open(51,file=filename)
        write(51,'(a12,100a10)')'variables = ',varname(1:nvar)
        write(51,'(a,i8,a,i8,a)')'zone nodes= ',np,' Elements= ',ne, ' ZONETYPE=FEQuadrilateral DATAPACKING=POINT'
        do i=1,np
            write(51,'(100e25.15)')xyz(i,1:nvar)
        end do
        do i=1,ne
            write(51,'(100i8)')ele(i,1:nep)
        end do
        close(51)
    
end program