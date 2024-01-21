program NSAWETinput
	character*100 filegrd,fileinp,outfile
	integer iblk,it(8),iblkcpu(100),iblkcpu1(100),iblkcpu2(100)
	integer,allocatable :: iblock(:,:),ibound(:,:,:),iNSAWET(:,:,:),idim1(:,:),idim2(:,:),inode(:,:)
	real*8,allocatable :: xyz(:,:,:,:)
	
	open(5,file='input.txt')
	read(5,*)fileinp
	read(5,*)filegrd
	read(5,*)ncpu
	close(5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	iblock(iblk,1) : nimax
!		  (iblk,2) : njmax
!		  (iblk,3) : nkmax
!		  (iblk,4) : bounds
!!!!!!!!!GridGen boundary file!!!!!!!!!!!!!!!!!!!!!!!
!	ibound(iblk,100,1) : ia
!		  (iblk,100,2) : ib
!		  (iblk,100,3) : ja
!		  (iblk,100,4) : jb
!		  (iblk,100,5) : ka
!		  (iblk,100,6) : kb
!		  (iblk,100,7) : boundary types ( <0 for patch)
!								  0  none
!								  1  interblock connection
! 								  2  solid surface
!  								  3  symmetry
!								  4  farfield
!								  5  inflow
!								  6  outflow
!								  71 i pole
!								  72 j pole
!								  73 k pole
!								  8  generic #1	 quasi p2p
!								  9  generic #2  surface p2p
!								  10 generic #3
!!!!!!!!!NSAWET boundary file!!!!!!!!!!!!!!!!!!!!!!!
!	iNSAWET(iblk,100, 1) : blkbd_id
!		   (iblk,100, 2) : idxa0
!		   (iblk,100, 3) : idxa1
!		   (iblk,100, 4) : idxb0
!		   (iblk,100, 5) : idxb1
!		   (iblk,100, 6) : isuf
!		   (iblk,100, 7) : kind (1:wall, 2:symmetry, 3:farfield, 4:patch)
!		   (iblk,100, 8) : kdsub
!		   (iblk,100, 9) : lsfpro
!		   (iblk,100,10) : iblkex
!		   (iblk,100,11) : iblkbd_ex
!		   (iblk,100,12) : lcross
!		   (iblk,100,13) : lreva
!		   (iblk,100,14) : lrevb
!		   (iblk,100,15) : npra
!		   (iblk,100,16) : nprb
!		   (iblk,100,17) : v#1
!		   (iblk,100,18) : v#2
!		   (iblk,100,19) : v#3
!		   (iblk,100,20) : v#4
!		   (iblk,100,21) : lhole
!!!!!!!!!!!read input file!!!!!!!!!!!!!!!!!!!!!!!!!!!
	open(5,file=fileinp)
	read(5,*)
	read(5,*)iblk
	allocate(iblock(iblk,4),ibound(iblk,10000,10),iNSAWET(iblk,10000,30))
	iblock=0
	ibound=0
	iNSAWET=0

	do i=1,iblk
		read(5,*)iblock(i,1),iblock(i,2),iblock(i,3)
		read(5,*)
		read(5,*)iblock(i,4)
		k=0
		do j=1,iblock(i,4)
			k=k+1
			read(5,*) (ibound(i,k,ii),ii=1,7)
!			write(*,'(7i5)')(ibound(i,k,ii),ii=1,7)
			if(ibound(i,k,7)<0)then
				k=k+1
				read(5,*) (ibound(i,k,ii),ii=1,7)
!				write(*,'(7i5)')(ibound(i,k,ii),ii=1,7)
			endif
		end do
!		pause
	end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	Sum the total points
	allocate(idim1(iblk,3),idim2(iblk,3),inode(ncpu,2))
	idim1=0
	idim2=0
	npall=0
	do i=1,iblk
		idim1(i,1)=i
		idim1(i,2)=(iblock(i,1)-1)*(iblock(i,2)-1)*(iblock(i,3)-1)
		npall=npall+idim1(i,2)
	end do
	idim2=idim1

	do i=1,iblk
		do j=i,iblk
			if(idim2(i,2) < idim2(j,2))then
				itmp1 = idim2(j,1)
				itmp2 = idim2(j,2)

				idim2(j,1)=idim2(i,1)
				idim2(j,2)=idim2(i,2)

				idim2(i,1)=itmp1
				idim2(i,2)=itmp2
			endif
		end do
	end do

	if(ncpu==1)then
		idim2(:,3)=1
		inode(1,1)=npall
		inode(1,2)=iblk
	else
		inode=0
		do i=1,iblk
			itmp=idim2(i,2)
			imin=1
			do j=1,ncpu
				if(inode(j,1)<inode(imin,1))imin=j
			end do
			inode(imin,1)=inode(imin,1)+itmp
			inode(imin,2)=inode(imin,2)+1
			idim2(i,3)=imin
		end do
	endif

	open(59,file='blockname.txt')
	open(60,file='dimension.txt')
		do icpu=1,ncpu
			write(60,*)inode(icpu,1)
			write(59,*)icpu,inode(icpu,2)
			
			do i=1,iblk
				ixnode=idim2(i,3)
				if(ixnode==icpu)then
					ib=idim2(i,1)

					if(ib<10)then
						write(outfile,'(i1)')ib
						outfile=trim('00')//trim(outfile)
					elseif(ib<100)then
						write(outfile,'(i2)')ib
						outfile=trim('0')//trim(outfile)
					elseif(ib<1000)then
						write(outfile,'(i3)')ib
					elseif(ib<10000)then
						write(outfile,'(i4)')ib
					endif
					outfile=trim('block')//trim(outfile)//trim('.dat')
					ni=iblock(ib,1)
					nj=iblock(ib,2)
					nk=iblock(ib,3)
					write(59,'(a20)')outfile
					write(60,'(a5,i5,a5,i5,a5,i5,a5,i5,a6,i12)')' blk=',ib,' ni=',ni,' nj=',nj,' nk=',nk,' Dim=',(ni-1)*(nj-1)*(nk-1)
				endif
			end do

		end do

	close(59)
	close(60)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	pause
!!!!!!!!!!NSAWET input file!!!!!!!!!!!!!!!!!!!!!!!!!!
	do nblk=1,iblk		! selectface
		k=0
		do nbound=1,iblock(nblk,4)
			k=k+1
			iNSAWET(nblk,k,1)=nbound
			call selectface(iblock,ibound,iNSAWET,iblk,nblk,k)
			if(ibound(nblk,k,7)<0)k=k+1
		end do
	end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do nblk=1,iblk		! selectibjb and kind
		k=0
		do nbound=1,iblock(nblk,4)
			k=k+1
			call selectibjb(iblock,ibound,iNSAWET,iblk,nblk,k)
			call selectkind(iblock,ibound,iNSAWET,iblk,nblk,k)
			if(ibound(nblk,k,7)<0) k=k+1
		end do
	end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do nblk=1,iblk		! select PATCH
		k=0
		do nbound=1,iblock(nblk,4)
			k=k+1
			if(ibound(nblk,k,7)<0)then
				call selectpatch(iblock,ibound,iNSAWET,iblk,nblk,k)
				k=k+1
!				write(*,*)iNSAWET(nblk,nbound,1),iNSAWET(nblk,nbound,10),iNSAWET(nblk,nbound,11)
			endif
		end do
	end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do nblk=1,iblk		! select lcross lreva lrevb
		k=0
		do nbound=1,iblock(nblk,4)
			k=k+1
			if(ibound(nblk,k,7)<0)then
				call selectlcross(iblock,ibound,iNSAWET,iblk,nblk,k)
!				write(*,*) iNSAWET(nblk,nbound,12),'lcross'
			endif
			if(ibound(nblk,k,7)<0) k=k+1
		end do
	end do
!!!!!!!!output NSAWET gridfile!!!!!!!!!!!!!!!!!!!!!!!!
	open(3,file=filegrd)
	read(3,*)
	read(3,*)((itemp,i=1,3),j=1,iblk)
	do i=1,iblk
		ni=iblock(i,1)
		nj=iblock(i,2)
		nk=iblock(i,3)
		allocate(xyz(ni,nj,nk,3))
		read(3,*)(((xyz(i1,j1,k1,1),i1=1,ni),j1=1,nj),k1=1,nk)
		read(3,*)(((xyz(i1,j1,k1,2),i1=1,ni),j1=1,nj),k1=1,nk)
		read(3,*)(((xyz(i1,j1,k1,3),i1=1,ni),j1=1,nj),k1=1,nk)
		if(i<10)then
			write(outfile,'(i1)')i
			outfile=trim('00')//trim(outfile)
		elseif(i<100)then
			write(outfile,'(i2)')i
			outfile=trim('0')//trim(outfile)
		elseif(i<1000)then
			write(outfile,'(i3)')i
!			outfile=trim('0')//trim(outfile)
		endif
		outfile=trim('block')//trim(outfile)//trim('.dat')
		open(4,file=outfile)
		write(4,'(a45)')'#--------------------------------------------'
		write(4,'(a19)')'#   NI     NJ    Nk'
		write(4,'(a1,3i8)')'#',ni,nj,nk
		write(4,'(a45)')'#--------------------------------------------'
		write(4,'(a49)')'#  blockID boundnum holenum gridlayer lengthscale'
		itemp=0
		scale=1.0
		ibnum=iblock(i,4)
		write(4,'(a1,4i8,f11.2)')'#',i,ibnum,itemp,itemp,scale
		write(4,'(a170)')'#   blkbd_id  idxa0   idxa1   idxb0   idxb1   isuf    kind  kdsub   lsfpro  iblkex iblkbd_ex lcross  lreva   lrevb    npra    nprb     v#1     v#2     v#3     v#4   lhole'
		k=0
		do j=1,ibnum
			k=k+1
			write(4,'(a1,21i8)')'#',(iNSAWET(i,k,ii),ii=1,21)
			if(ibound(i,k,7)<0) k=k+1
		end do
		write(4,'(a50)')'#---------name for surface patch file-------------'
		write(4,'(a50)')'#---------name for quasi patch file---------------'
		write(4,'(a50)')'#-------------------------------------------------'
		write(4,*)'  zone i= ',ni,'  j= ',nj,'  k= ',nk
		write(4,'(3e20.10)')((((xyz(i1,j1,k1,ii),ii=1,3),i1=1,ni),j1=1,nj),k1=1,nk)
		close(4)
		deallocate(xyz)
	end do
	close(3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	deallocate(iblock,ibound,iNSAWET)

end program

subroutine selectlcross(iblock,ibound,iNSAWET,iblk,nblk,nbound)
	integer iblock(iblk,4),ibound(iblk,10000,10),iNSAWET(iblk,10000,30)
	integer ipnew(4),ipold(4),lcross

	
	if(iNSAWET(nblk,nbound,6)==1.or.iNSAWET(nblk,nbound,6)==2)then
		do i=3,6
			ipold(i-2)=ibound(nblk,nbound,i)
		end do
	elseif(iNSAWET(nblk,nbound,6)==3.or.iNSAWET(nblk,nbound,6)==4)then
		do i=1,2
			ipold(i)=ibound(nblk,nbound,i)
		end do
		do i=5,6
			ipold(i-2)=ibound(nblk,nbound,i)
		end do
	elseif(iNSAWET(nblk,nbound,6)==5.or.iNSAWET(nblk,nbound,6)==6)then
		do i=1,4
			ipold(i)=ibound(nblk,nbound,i)
		end do
	endif

	if(ibound(nblk,nbound+1,1)==ibound(nblk,nbound+1,2))then
		do i=3,6
			ipnew(i-2)=ibound(nblk,nbound+1,i)
		end do
	elseif(ibound(nblk,nbound+1,3)==ibound(nblk,nbound+1,4))then
		do i=1,2
			ipnew(i)=ibound(nblk,nbound+1,i)
		end do
		do i=5,6
			ipnew(i-2)=ibound(nblk,nbound+1,i)
		end do
	elseif(ibound(nblk,nbound+1,5)==ibound(nblk,nbound+1,6))then
		do i=1,4
			ipnew(i)=ibound(nblk,nbound+1,i)
		end do
	endif
!	write(*,*)ipnew,'new'
!	write(*,*)ipold,'old'

	if(ipold(1)*ipnew(1)>0.and.ipold(3)*ipnew(3)>0)then
		lcross=0
	else
		lcross=1
	endif	
	iNSAWET(nblk,nbound,12)=lcross
!	write(*,*)nblk,nbound
!	write(*,*)ipold
!	write(*,*)ipnew
!pause

	if(lcross==0)then
		lreva=(ipold(2)-ipold(1))/(ipnew(2)-ipnew(1))
		lrevb=(ipold(4)-ipold(3))/(ipnew(4)-ipnew(3))
	endif
	if(lcross==1)then
		lreva=(ipold(2)-ipold(1))/(ipnew(4)-ipnew(3))
		lrevb=(ipold(4)-ipold(3))/(ipnew(2)-ipnew(1))
	endif
	lreva=iabs((lreva-1)/2)
	lrevb=iabs((lrevb-1)/2)
	iNSAWET(nblk,nbound,13)=lreva	
	iNSAWET(nblk,nbound,14)=lrevb

!	write(*,*)lreva,lrevb
	
end subroutine

subroutine selectpatch(iblock,ibound,iNSAWET,iblk,nblk,nbound)
	integer iblock(iblk,4),ibound(iblk,10000,10),iNSAWET(iblk,10000,30)
	integer ipnew(6)
	
	iNSAWET(nblk,nbound,10)=ibound(nblk,nbound+1,7)
	iblknew=iNSAWET(nblk,nbound,10)
	do i=1,6
		ipnew(i)=ibound(nblk,nbound+1,i)
	end do
!	write(*,'(9i5)')nblk,nbound,ipnew,iNSAWET(nblk,nbound,10)

	!write(*,*)iblock(iblknew,4)
	k=0
	do i=1,iblock(iblknew,4)
		!write(*,*)nblk ,i !,iblock(iblknew,4)
		k=k+1
		do j=1,6
			if(ibound(iblknew,k,j).ne.ipnew(j))goto 100
		end do
		iNSAWET(nblk,nbound,11)=i
		write(*,*)' Old',nblk,nbound,' New',iblknew,i

100	ijk=0	
		if(ibound(iblknew,k,7)<0)k=k+1
	end do

end subroutine

subroutine selectkind(iblock,ibound,iNSAWET,iblk,nblk,nbound)
	integer iblock(iblk,4),ibound(iblk,10000,10),iNSAWET(iblk,10000,30)
	
	if(ibound(nblk,nbound,7)==2) iNSAWET(nblk,nbound,7)=1
	if(ibound(nblk,nbound,7)==3) iNSAWET(nblk,nbound,7)=2
	if(ibound(nblk,nbound,7)==4) iNSAWET(nblk,nbound,7)=3
	if(ibound(nblk,nbound,7)< 0) iNSAWET(nblk,nbound,7)=4

	if(iNSAWET(nblk,nbound,7)<4) then
		do ii=8,21
			iNSAWET(nblk,nbound,ii)=0
		end do
	endif

	if(ibound(nblk,nbound,7)==5) then
		iNSAWET(nblk,nbound,7)=3
		iNSAWET(nblk,nbound,8)=1
	endif

	if(ibound(nblk,nbound,7)==6) then
		iNSAWET(nblk,nbound,7)=3
		iNSAWET(nblk,nbound,8)=2
	endif

	if(ibound(nblk,nbound,7)==7) then
		iNSAWET(nblk,nbound,7)=3
		iNSAWET(nblk,nbound,8)=3
	endif


	if(ibound(nblk,nbound,7)==8) then
		iNSAWET(nblk,nbound,7)=4
		iNSAWET(nblk,nbound,8)=2
	endif
	if(ibound(nblk,nbound,7)==9) then
		iNSAWET(nblk,nbound,7)=4
		iNSAWET(nblk,nbound,8)=1
	endif


end subroutine

subroutine selectibjb(iblock,ibound,iNSAWET,iblk,nblk,nbound)
	integer iblock(iblk,4),ibound(iblk,10000,10),iNSAWET(iblk,10000,30)
		
	if(iNSAWET(nblk,nbound,6)==1.or.iNSAWET(nblk,nbound,6)==2)then
		iNSAWET(nblk,nbound,2)=min(iabs(ibound(nblk,nbound,3)),iabs(ibound(nblk,nbound,4)))
		if(iNSAWET(nblk,nbound,2)==1)iNSAWET(nblk,nbound,2)=0
		
		iNSAWET(nblk,nbound,3)=max(iabs(ibound(nblk,nbound,3)),iabs(ibound(nblk,nbound,4)))-1
		if(iNSAWET(nblk,nbound,3)==iblock(nblk,2)-1)iNSAWET(nblk,nbound,3)=0
		
		iNSAWET(nblk,nbound,4)=min(iabs(ibound(nblk,nbound,5)),iabs(ibound(nblk,nbound,6)))
		if(iNSAWET(nblk,nbound,4)==1)iNSAWET(nblk,nbound,4)=0
		
		iNSAWET(nblk,nbound,5)=max(iabs(ibound(nblk,nbound,5)),iabs(ibound(nblk,nbound,6)))-1
		if(iNSAWET(nblk,nbound,5)==iblock(nblk,3)-1)iNSAWET(nblk,nbound,5)=0

	elseif(iNSAWET(nblk,nbound,6)==3.or.iNSAWET(nblk,nbound,6)==4)then
		iNSAWET(nblk,nbound,2)=min(iabs(ibound(nblk,nbound,1)),iabs(ibound(nblk,nbound,2)))
		if(iNSAWET(nblk,nbound,2)==1)iNSAWET(nblk,nbound,2)=0

		iNSAWET(nblk,nbound,3)=max(iabs(ibound(nblk,nbound,1)),iabs(ibound(nblk,nbound,2)))-1
		if(iNSAWET(nblk,nbound,3)==iblock(nblk,1)-1)iNSAWET(nblk,nbound,3)=0

		iNSAWET(nblk,nbound,4)=min(iabs(ibound(nblk,nbound,5)),iabs(ibound(nblk,nbound,6)))
		if(iNSAWET(nblk,nbound,4)==1)iNSAWET(nblk,nbound,4)=0

		iNSAWET(nblk,nbound,5)=max(iabs(ibound(nblk,nbound,5)),iabs(ibound(nblk,nbound,6)))-1
		if(iNSAWET(nblk,nbound,5)==iblock(nblk,3)-1)iNSAWET(nblk,nbound,5)=0

	elseif(iNSAWET(nblk,nbound,6)==5.or.iNSAWET(nblk,nbound,6)==6)then
		iNSAWET(nblk,nbound,2)=min(iabs(ibound(nblk,nbound,1)),iabs(ibound(nblk,nbound,2)))
		if(iNSAWET(nblk,nbound,2)==1)iNSAWET(nblk,nbound,2)=0

		iNSAWET(nblk,nbound,3)=max(iabs(ibound(nblk,nbound,1)),iabs(ibound(nblk,nbound,2)))-1
		if(iNSAWET(nblk,nbound,3)==iblock(nblk,1)-1)iNSAWET(nblk,nbound,3)=0

		iNSAWET(nblk,nbound,4)=min(iabs(ibound(nblk,nbound,3)),iabs(ibound(nblk,nbound,4)))
		if(iNSAWET(nblk,nbound,4)==1)iNSAWET(nblk,nbound,4)=0

		iNSAWET(nblk,nbound,5)=max(iabs(ibound(nblk,nbound,3)),iabs(ibound(nblk,nbound,4)))-1
		if(iNSAWET(nblk,nbound,5)==iblock(nblk,2)-1)iNSAWET(nblk,nbound,5)=0

	endif
endsubroutine


subroutine selectface(iblock,ibound,iNSAWET,iblk,nblk,nbound)
	integer iblock(iblk,4),ibound(iblk,10000,10),iNSAWET(iblk,10000,30)

	if(ibound(nblk,nbound,1)==ibound(nblk,nbound,2))then
		if(ibound(nblk,nbound,1)==1)then
			iNSAWET(nblk,nbound,6)=1
		elseif(ibound(nblk,nbound,1)==iblock(nblk,1))then
			iNSAWET(nblk,nbound,6)=2
		endif
	elseif(ibound(nblk,nbound,3)==ibound(nblk,nbound,4))then
		if(ibound(nblk,nbound,3)==1)then
			iNSAWET(nblk,nbound,6)=3
		elseif(ibound(nblk,nbound,3)==iblock(nblk,2))then
			iNSAWET(nblk,nbound,6)=4
		endif
	elseif(ibound(nblk,nbound,5)==ibound(nblk,nbound,6))then
		if(ibound(nblk,nbound,5)==1)then
			iNSAWET(nblk,nbound,6)=5
		elseif(ibound(nblk,nbound,5)==iblock(nblk,3))then
			iNSAWET(nblk,nbound,6)=6
		endif
	endif
end subroutine