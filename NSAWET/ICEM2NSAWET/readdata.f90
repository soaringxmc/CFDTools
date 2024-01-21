subroutine readdata(xyzpatch,xyzpatch1,ninter,nblk,nquasi,filegrd,ix,iblock,ireal,iblk,  &
    xyzall_x,xyzall_y,xyzall_z,xyzall1_x,xyzall1_y,xyzall1_z,iblk_index,npointsall)
    character*100 filegrd
    real*4 xyzpatch(ninter,2,2,3),temp
    real*8 xyzpatch1(ninter,2,2,3),temp1
    integer nquasi(2,ninter,10),iblock(iblk,4),ireal
    real,allocatable :: xyz(:,:,:,:)
	real*8,allocatable :: xyz1(:,:,:,:)
    real xyzall_x(npointsall),xyzall_y(npointsall),xyzall_z(npointsall)
    real*8 xyzall1_x(npointsall),xyzall1_y(npointsall),xyzall1_z(npointsall)
    integer iblk_index(iblk)
    
    nqua=nquasi(1,ix,1)
    nblk=nqua/100
    iface1=mod(nqua,100)/10
    iface2=mod(nqua,10)
    idxie1=nquasi(1,ix,2)
    idxie2=nquasi(1,ix,3)
    ideta1=nquasi(1,ix,4)
    ideta2=nquasi(1,ix,5)
    
    if(iface1==2 .or. iface1==3)then
        idxie1=nquasi(1,ix,4)
        idxie2=nquasi(1,ix,5)
        ideta1=nquasi(1,ix,2)
        ideta2=nquasi(1,ix,3)
    endif
!    write(*,*)'nblk',ix,nblk,idxie1,ideta1,idxie2,ideta2
    idim=iblock(nblk,1)
    jdim=iblock(nblk,2)
    kdim=iblock(nblk,3)
    
!    write(*,*)nqua,nblk,idim,jdim,kdim
!    write(*,*)idxie1,idxie2,ideta1,ideta2
!    write(*,*)ix
    
	do i=1,nblk
		ni=iblock(i,1)
		nj=iblock(i,2)
		nk=iblock(i,3)
        istart=iblk_index(i)
		if(ireal==4)then
			allocate(xyz(ni,nj,nk,3))
		    do k1=1,nk
			    do j1=1,nj
				    do i1=1,ni
                        id=istart+(k1-1)*nj*ni+(j1-1)*ni+i1
					    xyz(i1,j1,k1,1)=xyzall_x(id)
                        xyz(i1,j1,k1,2)=xyzall_y(id)
                        xyz(i1,j1,k1,3)=xyzall_z(id)
				    end do
			    end do
		    end do
		endif
		if(ireal==8)then
			!write(*,*)"Block",i," Double precision!"
			allocate(xyz1(ni,nj,nk,3))
		    do k1=1,nk
			    do j1=1,nj
				    do i1=1,ni
                        id=istart+(k1-1)*nj*ni+(j1-1)*ni+i1
					    xyz1(i1,j1,k1,1)=xyzall1_x(id)
                        xyz1(i1,j1,k1,2)=xyzall1_y(id)
                        xyz1(i1,j1,k1,3)=xyzall1_z(id)
				    end do
			    end do
		    end do
        endif
		if(ireal==4 .and. i<nblk)deallocate(xyz)
		if(ireal==8 .and. i<nblk)deallocate(xyz1)
	end do
	close(3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(ireal==4)then
        if(iface1==1)then
            if(iface2==1)idx=1
            if(iface2==2)idx=idim
            xyzpatch(ix,1,1,1:3)=xyz(idx,idxie1,ideta1,1:3)
            xyzpatch(ix,2,1,1:3)=xyz(idx,idxie2,ideta1,1:3)
            xyzpatch(ix,1,2,1:3)=xyz(idx,idxie1,ideta2,1:3)
            xyzpatch(ix,2,2,1:3)=xyz(idx,idxie2,ideta2,1:3)
        elseif(iface1==2)then
            if(iface2==1)idx=1
            if(iface2==2)idx=jdim
            xyzpatch(ix,1,1,1:3)=xyz(idxie1,idx,ideta1,1:3)
            xyzpatch(ix,2,1,1:3)=xyz(idxie2,idx,ideta1,1:3)
            xyzpatch(ix,1,2,1:3)=xyz(idxie1,idx,ideta2,1:3)
            xyzpatch(ix,2,2,1:3)=xyz(idxie2,idx,ideta2,1:3)
        elseif(iface1==3)then
            if(iface2==1)idx=1
            if(iface2==2)idx=kdim
            xyzpatch(ix,1,1,1:3)=xyz(idxie1,ideta1,idx,1:3)
            xyzpatch(ix,2,1,1:3)=xyz(idxie2,ideta1,idx,1:3)
            xyzpatch(ix,1,2,1:3)=xyz(idxie1,ideta2,idx,1:3)
            xyzpatch(ix,2,2,1:3)=xyz(idxie2,ideta2,idx,1:3)            
        endif
    elseif(ireal==8)then
        !write(*,*)iface1,iface2
        if(iface1==1)then
            if(iface2==1)idx=1
            if(iface2==2)idx=idim
        !write(*,*)idx,idxie1,ideta1,idxie2,ideta2
        
            xyzpatch1(ix,1,1,1:3)=xyz1(idx,idxie1,ideta1,1:3)
            xyzpatch1(ix,2,1,1:3)=xyz1(idx,idxie2,ideta1,1:3)
            xyzpatch1(ix,1,2,1:3)=xyz1(idx,idxie1,ideta2,1:3)
            xyzpatch1(ix,2,2,1:3)=xyz1(idx,idxie2,ideta2,1:3)
        elseif(iface1==2)then
            if(iface2==1)idx=1
            if(iface2==2)idx=jdim
            write(*,*)idxie1,idxie2,ideta1,ideta2,idx
            xyzpatch1(ix,1,1,1:3)=xyz1(idxie1,idx,ideta1,1:3)
            xyzpatch1(ix,2,1,1:3)=xyz1(idxie2,idx,ideta1,1:3)
            xyzpatch1(ix,1,2,1:3)=xyz1(idxie1,idx,ideta2,1:3)
            xyzpatch1(ix,2,2,1:3)=xyz1(idxie2,idx,ideta2,1:3)
        elseif(iface1==3)then
            if(iface2==1)idx=1
            if(iface2==2)idx=kdim
            xyzpatch1(ix,1,1,1:3)=xyz1(idxie1,ideta1,idx,1:3)
            xyzpatch1(ix,2,1,1:3)=xyz1(idxie2,ideta1,idx,1:3)
            xyzpatch1(ix,1,2,1:3)=xyz1(idxie1,ideta2,idx,1:3)
            xyzpatch1(ix,2,2,1:3)=xyz1(idxie2,ideta2,idx,1:3)
          !  write(*,*)xyzpatch1(ix,1,1,1:3)
          !  write(*,*)idx,idxie1,ideta1
        endif
    endif
		if(ireal==4 )deallocate(xyz)
		if(ireal==8 )deallocate(xyz1)
            

    end subroutine
    
    
    
subroutine readpoints(xyz_x,xyz_y,xyz_z,xyz1_x,xyz1_y,xyz1_z,ireal,iblk,iblock,npointsall,iblk_index,filegrd)
    character*100 filegrd
    real xyz_x(npointsall),xyz_y(npointsall),xyz_z(npointsall),temp
    real*8 xyz1_x(npointsall),xyz1_y(npointsall),xyz1_z(npointsall),temp1
    integer iblock(iblk,4),ireal,iblk_index(iblk)
    
  !  write(*,*)iblk_index
  !  write(*,*)npointsall
    open(3,file=filegrd,form='binary',status='old')
	read(3)itemp
	read(3) itemp
	!call readword(iblk)
	read(3)itemp
	read(3)itemp
	do j=1,iblk
	do i=1,3
	read(3)itemp
    call readword(itemp)
	!write(*,*)itemp
	end do
	end do
	read(3)itemp
	do i=1,iblk
        Write(*,*)'Read block data: ',i
        istart=iblk_index(i)
		ni=iblock(i,1)
		nj=iblock(i,2)
		nk=iblock(i,3)
		read(3)itemp
		call readword(itemp)
!		write(*,*)itemp
		ireal=itemp/ni/nj/nk/3
        !write(*,*)ireal
       !write(*,*)istart
		if(ireal==4)then
		do k1=1,nk
			do j1=1,nj
				do i1=1,ni
                    id=istart+(k1-1)*nj*ni+(j1-1)*ni+i1
					read(3)temp
					call readword(temp)
                !   if(i>=41) write(*,*)istart,i1,j1,k1,id,temp
					xyz_x(id)=temp
                 !  if(i>=41) write(*,*)ni,nj,nk,xyz_x(id)
				end do
			end do
		end do
		do k1=1,nk
			do j1=1,nj
				do i1=1,ni
                    id=istart+(k1-1)*nj*ni+(j1-1)*ni+i1
					read(3)temp
					call readword(temp)
					xyz_y(id)=temp
				end do
			end do
		end do
		do k1=1,nk
			do j1=1,nj
				do i1=1,ni
                    id=istart+(k1-1)*nj*ni+(j1-1)*ni+i1
					read(3)temp
					call readword(temp)
					xyz_z(id)=temp
				end do
			end do
		end do
		endif
		if(ireal==8)then
			!write(*,*)"Block",i," Double precision!"
		do k1=1,nk
			do j1=1,nj
				do i1=1,ni
                    id=istart+(k1-1)*nj*ni+(j1-1)*ni+i1
					read(3)temp1
					call readword1(temp1)
					xyz1_x(id)=temp1
				end do
			end do
		end do
		do k1=1,nk
			do j1=1,nj
				do i1=1,ni
                    id=istart+(k1-1)*nj*ni+(j1-1)*ni+i1
					read(3)temp1
					call readword1(temp1)
					xyz1_y(id)=temp1
				end do
			end do
		end do
		do k1=1,nk
			do j1=1,nj
				do i1=1,ni
                    id=istart+(k1-1)*nj*ni+(j1-1)*ni+i1
                    !if(i>=39)write(*,*)id
					read(3)temp1
					call readword1(temp1)
					xyz1_z(id)=temp1
				end do
			end do
		end do
		endif

		read(3)itemp
		call readword(itemp)
	end do
	close(3)            
!write(*,*)ireal
!pause
end subroutine