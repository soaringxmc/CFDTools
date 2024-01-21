program readgrid
!!!!!!!!output NSAWET gridfile!!!!!!!!!!!!!!!!!!!!!!!!
	integer,allocatable::  iblock(:,:)
	character*80 filegrd,filegrd1
	real,allocatable:: xyz(:,:,:,:)
	real(8),allocatable:: xyz1(:,:,:,:)
	real(8) temp1
	filegrd='cfl3d.xyz'
	filegrd1='cfl3d.xyz2'
	open(3,file=filegrd,form='binary')
	open(4,file=filegrd1,form='binary')
	read(3)itemp
	read(3)iblk
	call readword(iblk)
	allocate(iblock(iblk,3))
	read(3)itemp
	read(3)itemp
	do j=1,iblk
		do i=1,3
			read(3)iblock(j,i)
			call readword(iblock(j,i))
		end do
	end do
	read(3)itemp

	write(4)iblk
	do j=1,iblk
		do i=1,3
			write(4)iblock(j,i)
		end do
	end do

	do i=1,iblk
		ni=iblock(i,1)
		nj=iblock(i,2)
		nk=iblock(i,3)
		read(3)itemp
		call readword(itemp)
!		write(*,*)itemp
		ireal=itemp/ni/nj/nk/3
		if(ireal==4)then
			write(*,'(A5,i5,a4,i5,a4,i5,a4,i5,a18)')"Block",i,' ni=',ni,' nj=',nj,' nk=',nk," Single precision!"
			allocate(xyz(ni,nj,nk,3))
		do k1=1,nk
			do j1=1,nj
				do i1=1,ni
					read(3)temp
					call readword(temp)
					xyz(i1,j1,k1,1)=temp
				end do
			end do
		end do
		do k1=1,nk
			do j1=1,nj
				do i1=1,ni
					read(3)temp
					call readword(temp)
					xyz(i1,j1,k1,2)=temp
				end do
			end do
		end do
		do k1=1,nk
			do j1=1,nj
				do i1=1,ni
					read(3)temp
					call readword(temp)
					xyz(i1,j1,k1,3)=temp
				end do
			end do
		end do
		write(4)(((xyz(i1,j1,k1,1),i1=1,ni),j1=1,nj),k1=1,nk)
		write(4)(((xyz(i1,j1,k1,2),i1=1,ni),j1=1,nj),k1=1,nk)
		write(4)(((xyz(i1,j1,k1,3),i1=1,ni),j1=1,nj),k1=1,nk)
		deallocate(xyz)

		else if(ireal==8)then
			write(*,'(A5,i5,a4,i5,a4,i5,a4,i5,a18)')"Block",i,' ni=',ni,' nj=',nj,' nk=',nk," Double precision!"
			allocate(xyz1(ni,nj,nk,3))
		do k1=1,nk
			do j1=1,nj
				do i1=1,ni
					read(3)temp1
					call readword1(temp1)
					xyz1(i1,j1,k1,1)=temp1
				end do
			end do
		end do
		do k1=1,nk
			do j1=1,nj
				do i1=1,ni
					read(3)temp1
					call readword1(temp1)
					xyz1(i1,j1,k1,2)=temp1
				end do
			end do
		end do
		do k1=1,nk
			do j1=1,nj
				do i1=1,ni
					read(3)temp1
					call readword1(temp1)
					xyz1(i1,j1,k1,3)=temp1
				end do
			end do
		end do
		write(4)(((xyz1(i1,j1,k1,1),i1=1,ni),j1=1,nj),k1=1,nk)
		write(4)(((xyz1(i1,j1,k1,2),i1=1,ni),j1=1,nj),k1=1,nk)
		write(4)(((xyz1(i1,j1,k1,3),i1=1,ni),j1=1,nj),k1=1,nk)
		deallocate(xyz1)
		endif

		read(3)itemp
		call readword(itemp)
	
	end do
	close(3);close(4)
end program

subroutine readword(i)
  integer i,j
  j=i
  call JMVBITS (j, 0, 8, i, 24)
  call JMVBITS (j, 8, 8, i, 16)
  call JMVBITS (j, 16, 8, i, 8)
  call JMVBITS (j, 24, 8, i, 0)
  return
end subroutine

subroutine readword1(i)
  integer*8 i,j
  j=i
  call kMVBITS (j, 0, 8, i, 56)
  call kMVBITS (j, 8, 8, i, 48)
  call kMVBITS (j, 16, 8, i,40)
  call kMVBITS (j, 24, 8, i,32)
  call kMVBITS (j, 32, 8, i,24)
  call kMVBITS (j, 40, 8, i,16)
  call kMVBITS (j, 48, 8, i,8)
  call kMVBITS (j, 56, 8, i,0)
  return
end subroutine