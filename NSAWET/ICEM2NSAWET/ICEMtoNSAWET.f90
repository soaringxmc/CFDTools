program ICEMtoNSAWET
	character*100 filegrd,fileinp,outfile,comm
    real*8,external::dis
	integer iblk,it(8),iblkcpu(100),iblkcpu1(100),iblkcpu2(100)
	real*8 temp1
	integer,allocatable :: iblock(:,:),ibound(:,:,:),iNSAWET(:,:,:),ibc(:,:),ioto(:,:,:),idim1(:,:),idim2(:,:),inode(:,:),idim3(:,:)
    integer,allocatable :: iblk_index(:)
	integer,allocatable :: nload(:),npoints(:),nmpi_np(:,:),npatch(:,:)
    integer,allocatable :: nquasi(:,:,:)
	real,allocatable :: xyz(:,:,:,:)
	real*8,allocatable :: xyz1(:,:,:,:)
    real,allocatable :: xyzall_x(:),xyzall_y(:),xyzall_z(:)
    real*8,allocatable :: xyzall1_x(:),xyzall1_y(:),xyzall1_z(:)
    real,allocatable :: xyzpatch(:,:,:,:)
    real*8,allocatable :: xyzpatch1(:,:,:,:)
    real*8 dmax1,dmax2
	
	open(5,file='ICEMtoNSAWET.txt')
	read(5,*)fileinp
	read(5,*)filegrd
	read(5,*)ncpu
!	read(5,*)tore
	close(5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	iblock(iblk,1) : nimax
!		  (iblk,2) : njmax
!		  (iblk,3) : nkmax
!		  (iblk,4) : bounds
!!!!!!!!!CFL3D boundary file!!!!!!!!!!!!!!!!!!!!!!!
!	ibound(iblk,100,1) : ia
!CFL3D Boundary type.
!
!1000 free stream
!1001 general symmetry plane
!1002 extrapolation
!1003 inflow/outflow
!1004 (no longer available, use 2004 instead)
!1005 inviscid surface
!1008 constant enthalpy and entropy inflow
!1011 singular axis - half-plane symmetry
!1012 singular - full plane
!1013 singular axis - partial plane
!2002 specified pressure ratio
!2003 inflow with specified total conditions
!2004 no-slip wall
!2005 periodic in space
!2006 set pressure to satisfy the radial equilibrium equation
!2007 set all primitive variables
!2102 pressure ratio specified as a sinusoidal function of time
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
	do i=1,21
	read(5,*)
	end do
	read(5,*)
	read(5,*)iblk
	iblk=abs(iblk)
	iblkall=iblk
	allocate(iblock(iblk,4),iblk_index(iblk),ibound(iblk,1000,10),iNSAWET(iblk,1000,30),ibc(iblk,10))
	allocate(npoints(ncpu),nmpi_np(ncpu,ncpu),npatch(ncpu,ncpu),nload(iblk))
	if(iblk>999)then
		write(*,*)'Block larger than 1000, not allowed in NSAWET'
		stop
	endif
	iblock=0
	ibound=0
	iNSAWET=0
	ibc=0
	
	read(5,*)
	do i=1,iblk	!NCG       IEM  IADVANCE    IFORCE  IVISC(I)  IVISC(J)  IVISC(K)
		read(5,*)
	end do
	read(5,*)
    iblk_index(1)=1
	do i=1,iblk	!IDIM      JDIM      KDIM
		read(5,*)iblock(i,1),iblock(i,2),iblock(i,3)
        if(i>=2)iblk_index(i)=iblk_index(i-1)+iblock(i-1,1)*iblock(i-1,2)*iblock(i-1,3)
    end do

    npointsall_index=iblk_index(iblk)+iblock(iblk,1)*iblock(iblk,2)*iblock(iblk,3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	Sum the total points
	allocate(idim1(iblk,3),idim2(iblk,3),inode(ncpu,2),idim3(iblk,2))
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

!	write(*,'(3i10)')(idim1(i,1:3),i=1,iblk)
!	write(*,*)
!	write(*,'(3i10)')(idim2(i,1:3),i=1,iblk)
!	write(*,*)
!	write(*,*)inode


	open(59,file='blockname.txt')
	open(60,file='dimension.txt')
		do icpu=1,ncpu
			write(60,*)icpu,inode(icpu,1)
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
					idim3(ib,1)=icpu
				endif
			end do

		end do

	close(59)
!	close(60)
!	open(61,file='aaa.txt')
!	do i=1,iblk
!		write(61,*)idim3(i,1)-1
!	end do
!	close(61)
!	pause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	pause
!

	read(5,*)	!ILAMLO    ILAMHI    JLAMLO    JLAMHI    KLAMLO    ILAMHI
	do i=1,iblk
		read(5,*)
	end do
	read(5,*)	!INEWG    IGRIDC        IS        JS        KS        IE        JE        KE
	do i=1,iblk
		read(5,*)
	end do
	read(5,*)	!IDIAG(I)  IDIAG(J)  IDIAG(K)  IFLIM(I)  IFLIM(J)  IFLIM(K)
	do i=1,iblk
		read(5,*)
	end do	
	read(5,*)	!IFDS(I)   IFDS(J)   IFDS(K)  RKAP0(I)  RKAP0(J)  RKAP0(K)
	do i=1,iblk
		read(5,*)
	end do	
	read(5,*)	!GRID     NBCIO   NBCIDIM     NBCJO   NBCJDIM     NBCKO   NBCKDIM    IOVRLP
	do i=1,iblk
		read(5,*)itemp,ibc(i,1:6)
		iblock(i,4)=0
		do ix=1,6
			iblock(i,4)=iblock(i,4)+ibc(i,ix)
		end do
	end do
	read(5,*)	!I0:   GRID   SEGMENT    BCTYPE      JSTA      JEND      KSTA      KEND     NDATA
	do i=1,iblk
		iblock(i,4)=0
		do io=1,ibc(i,1)
			read(5,*)it(1:8)
		if(it(3)>0)then
			iblock(i,4)=iblock(i,4)+1
			if(i.ne.it(1))write(*,*)'Error blk id!'
			ibound(i,iblock(i,4),1)=1
			ibound(i,iblock(i,4),2)=1
			ibound(i,iblock(i,4),3)=it(4)
			ibound(i,iblock(i,4),4)=it(5)
			ibound(i,iblock(i,4),5)=it(6)
			ibound(i,iblock(i,4),6)=it(7)
			ibound(i,iblock(i,4),7)=it(3)
		endif
			if(it(8)>0)then
				read(5,*)
				read(5,*)
			end if
		end do
	end do
	read(5,*)	!IDIM: GRID   SEGMENT    BCTYPE      JSTA      JEND      KSTA      KEND     NDATA
	do i=1,iblk
		do io=1,ibc(i,2)
			read(5,*)it(1:8)
		if(it(3)>0)then
			iblock(i,4)=iblock(i,4)+1
			if(i.ne.it(1))write(*,*)'Error blk id!'
			ibound(i,iblock(i,4),1)=iblock(i,1)
			ibound(i,iblock(i,4),2)=iblock(i,1)
			ibound(i,iblock(i,4),3)=it(4)
			ibound(i,iblock(i,4),4)=it(5)
			ibound(i,iblock(i,4),5)=it(6)
			ibound(i,iblock(i,4),6)=it(7)
			ibound(i,iblock(i,4),7)=it(3)
		endif
			if(it(8)>0)then
				read(5,*)
				read(5,*)
			end if
		end do
	end do
	read(5,*)	!J0:   GRID   SEGMENT    BCTYPE      ISTA      IEND      KSTA      KEND     NDATA
	do i=1,iblk
		do io=1,ibc(i,3)
			read(5,*)it(1:8)
		if(it(3)>0)then
			iblock(i,4)=iblock(i,4)+1
			if(i.ne.it(1))write(*,*)'Error blk id!'
			ibound(i,iblock(i,4),1)=it(4)
			ibound(i,iblock(i,4),2)=it(5)
			ibound(i,iblock(i,4),3)=1
			ibound(i,iblock(i,4),4)=1
			ibound(i,iblock(i,4),5)=it(6)
			ibound(i,iblock(i,4),6)=it(7)
			ibound(i,iblock(i,4),7)=it(3)
		endif
			if(it(8)>0)then
				read(5,*)
				read(5,*)
			end if
		end do
	end do
	read(5,*)	!JDIM: GRID   SEGMENT    BCTYPE      ISTA      IEND      KSTA      KEND     NDATA
	do i=1,iblk
		do io=1,ibc(i,4)
			read(5,*)it(1:8)
		if(it(3)>0)then
			iblock(i,4)=iblock(i,4)+1
			if(i.ne.it(1))write(*,*)'Error blk id!'
			ibound(i,iblock(i,4),1)=it(4)
			ibound(i,iblock(i,4),2)=it(5)
			ibound(i,iblock(i,4),3)=iblock(i,2)
			ibound(i,iblock(i,4),4)=iblock(i,2)
			ibound(i,iblock(i,4),5)=it(6)
			ibound(i,iblock(i,4),6)=it(7)
			ibound(i,iblock(i,4),7)=it(3)
		endif
			if(it(8)>0)then
				read(5,*)
				read(5,*)
			end if
		end do
	end do
	read(5,*)	!K0:   GRID   SEGMENT    BCTYPE      ISTA      IEND      JSTA      JEND     NDATA
	do i=1,iblk
		do io=1,ibc(i,5)
			read(5,*)it(1:8)
		if(it(3)>0)then
			iblock(i,4)=iblock(i,4)+1
			if(i.ne.it(1))write(*,*)'Error blk id!'
			ibound(i,iblock(i,4),1)=it(4)
			ibound(i,iblock(i,4),2)=it(5)
			ibound(i,iblock(i,4),3)=it(6)
			ibound(i,iblock(i,4),4)=it(7)
			ibound(i,iblock(i,4),5)=1
			ibound(i,iblock(i,4),6)=1
			ibound(i,iblock(i,4),7)=it(3)
		endif
			if(it(8)>0)then
				read(5,*)
				read(5,*)
			end if
		end do
	end do
	read(5,*)	!KDIM: GRID   SEGMENT    BCTYPE      ISTA      IEND      JSTA      JEND     NDATA
	do i=1,iblk
		do io=1,ibc(i,6)
			read(5,*)it(1:8)
		if(it(3)>0)then
			iblock(i,4)=iblock(i,4)+1
			if(i.ne.it(1))write(*,*)'Error blk id!'
			ibound(i,iblock(i,4),1)=it(4)
			ibound(i,iblock(i,4),2)=it(5)
			ibound(i,iblock(i,4),3)=it(6)
			ibound(i,iblock(i,4),4)=it(7)
			ibound(i,iblock(i,4),5)=iblock(i,3)
			ibound(i,iblock(i,4),6)=iblock(i,3)
			ibound(i,iblock(i,4),7)=it(3)
		endif
			if(it(8)>0)then
				read(5,*)
				read(5,*)
			end if
		end do
	end do

	read(5,*)
	read(5,*)nseq
	do i=1,3
		read(5,*)
	end do

	do i=1,nseq
		read(5,*)
	end do
	read(5,*)
	do i=1,nseq
		read(5,*)
	end do
	
	read(5,*) !1-1 BLOCKING DATA

	read(5,*)!NBLI
	read(5,*)noto
	allocate(ioto(2,noto,20))
	read(5,*)!NUMBER   GRID     :    ISTA   JSTA   KSTA   IEND   JEND   KEND  ISVA1  ISVA2
	do i=1,noto
		read(5,*)ioto(1,i,1:10)
	end do
	read(5,*)!NUMBER   GRID     :    ISTA   JSTA   KSTA   IEND   JEND   KEND  ISVA1  ISVA2
	do i=1,noto
		read(5,*)ioto(2,i,1:10)
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    read(5,*)!PATCH SURFACE DATA:
    read(5,*)!NINTER
    read(5,*)nronnie
    if(nronnie==-1)then
        open(51,file='ronnie.inp',status='old')
        do i=1,8
            read(51,*)
        end do
        read(51,*)!iblk
        read(51,*)!NCG       IEM      IDIM      JDIM      KDIM
        do i=1,iblk
            read(51,*)
        end do
        read(51,*)!NINTER
        read(51,*)ninter        !有多少个搭接面。
        if(mod(ninter,2).ne.0)then
            write(*,*)'Patch Surface Wrong in Ronnie.inp!'
            stop
        endif
        allocate(nquasi(2,ninter,15))
        read(51,*) !INT     IIFIT    LLIMIT    IITMAX   MMCXIE  MMCETA    IIC0  IIORPH
        do i=1,ninter
            read(51,*)
        end do
        read(51,*)!INT    TO   XIE1  XIE2   ETA1   ETA2   NFB (one per int)
        read(51,*)!FROM   XIE1  XIE2   ETA1   ETA2       (repeat nfb times for each int)
        do i=1,ninter
            read(51,*)temp,nquasi(1,i,1:6)
            read(51,*)     nquasi(2,i,1:5)
        end do
        
        do i=1,ninter,2
            do ix=1,5
                i1=nquasi(1,i,ix)
                i2=nquasi(2,i+1,ix)
                if(i1 .ne. i2)then
                    write(*,*)'Patch surface index wrong in ronnie! 1#'
                    stop
                endif
            end do
        end do
        do i=2,ninter,2
            do ix=1,5
                i1=nquasi(1,i,ix)
                i2=nquasi(2,i-1,ix)
                if(i1 .ne. i2)then
                    write(*,*)'Patch surface index wrong in ronnie! 2#'
                    stop
                endif
            end do
        end do
        
        close(51)
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
	close(5)

	ibc=0
	do i=1,noto
		iblk=ioto(1,i,2)

		ibc(iblk,1)=ibc(iblk,1)+1
		iblock(iblk,4)=iblock(iblk,4)+1

		ib1=iblock(iblk,4)+ibc(iblk,1)-1
		ib2=iblock(iblk,4)+ibc(iblk,1)
		
		iflag=ioto(1,i,9)
		ix1=1
		ix2=1
		ix3=1
		if(iflag==1)ix1=-1
		if(iflag==2)ix2=-1
		if(iflag==3)ix3=-1

		ibound(iblk,ib1,1)=ioto(1,i,3)*ix1
		ibound(iblk,ib1,2)=ioto(1,i,6)*ix1
		ibound(iblk,ib1,3)=ioto(1,i,4)*ix2
		ibound(iblk,ib1,4)=ioto(1,i,7)*ix2
		ibound(iblk,ib1,5)=ioto(1,i,5)*ix3
		ibound(iblk,ib1,6)=ioto(1,i,8)*ix3
		ibound(iblk,ib1,7)=-1

		iflag=ioto(2,i,9)
		ix1=1
		ix2=1
		ix3=1
		if(iflag==1)ix1=-1
		if(iflag==2)ix2=-1
		if(iflag==3)ix3=-1

		ibound(iblk,ib2,1)=ioto(2,i,3)*ix1
		ibound(iblk,ib2,2)=ioto(2,i,6)*ix1
		ibound(iblk,ib2,3)=ioto(2,i,4)*ix2
		ibound(iblk,ib2,4)=ioto(2,i,7)*ix2
		ibound(iblk,ib2,5)=ioto(2,i,5)*ix3
		ibound(iblk,ib2,6)=ioto(2,i,8)*ix3
		ibound(iblk,ib2,7)=ioto(2,i,2)
	end do

	do i=1,noto
		iblk=ioto(2,i,2)

		ibc(iblk,1)=ibc(iblk,1)+1
		iblock(iblk,4)=iblock(iblk,4)+1

		ib1=iblock(iblk,4)+ibc(iblk,1)-1
		ib2=iblock(iblk,4)+ibc(iblk,1)
		
		iflag=ioto(2,i,9)
		ix1=1
		ix2=1
		ix3=1
		if(iflag==1)ix1=-1
		if(iflag==2)ix2=-1
		if(iflag==3)ix3=-1

		ibound(iblk,ib1,1)=ioto(2,i,3)*ix1
		ibound(iblk,ib1,2)=ioto(2,i,6)*ix1
		ibound(iblk,ib1,3)=ioto(2,i,4)*ix2
		ibound(iblk,ib1,4)=ioto(2,i,7)*ix2
		ibound(iblk,ib1,5)=ioto(2,i,5)*ix3
		ibound(iblk,ib1,6)=ioto(2,i,8)*ix3
		ibound(iblk,ib1,7)=-1

		iflag=ioto(1,i,9)
		ix1=1
		ix2=1
		ix3=1
		if(iflag==1)ix1=-1
		if(iflag==2)ix2=-1
		if(iflag==3)ix3=-1

		ibound(iblk,ib2,1)=ioto(1,i,3)*ix1
		ibound(iblk,ib2,2)=ioto(1,i,6)*ix1
		ibound(iblk,ib2,3)=ioto(1,i,4)*ix2
		ibound(iblk,ib2,4)=ioto(1,i,7)*ix2
		ibound(iblk,ib2,5)=ioto(1,i,5)*ix3
		ibound(iblk,ib2,6)=ioto(1,i,8)*ix3
		ibound(iblk,ib2,7)=ioto(1,i,2)
    end do

	write(*,*)iblkall
	do i=1,iblkall
		write(*,*)iblock(i,4),ibc(i,1)
		do j=1,iblock(i,4)+ibc(i,1)
			write(*,'(9i5)')i,j,ibound(i,j,1:7)
		end do
	end do
!	pause

	iblk=iblkall
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i=1,iblk
		nload(i)=idim3(i,1)
	end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	统计每一个进程的网格点数
	npointsall=0
	npoints=0
	do i=1,iblk
		icpu=nload(i)
		np=(iblock(i,1)-1)*(iblock(i,2)-1)*(iblock(i,3)-1)
		npoints(icpu)=npoints(icpu)+np
		npointsall=npointsall+np
	end do
!	write(*,*)npoints,npointsall
!	找出最大的和最小的，算出比值
	temp1=-1.0e20
	temp2= 1.0e20
	do i=1,ncpu
		if(npoints(i)>temp1)temp1=npoints(i)
		if(npoints(i)<temp2)temp2=npoints(i)
	!	write(*,*)i,npoints(i)
	end do
	point_max=temp1			!所有节点中网格点数最多的值
	point_min=temp2			!所有节点中网格点数最少的值
	averge=float(npointsall)/float(ncpu)
!	write(*,*)averge,point_max,point_min
	point_ratio_max=point_max/(averge+1e-10)	!网格最多与平均的比值
	point_ratio_min=point_min/(averge+1e-10)	!网格最少与平均的比值

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	nmpi_np=0				!统计所有进程中有多少个点需要通过mpi才传递
	npatch=0
	do i=1,noto
		iblk_old=ioto(1,i,2)	!自己的块序号
		iblk_new=ioto(2,i,2)	!对方的块序号
		np_old=(iabs(ioto(1,i,3)-ioto(1,i,6))+1)*(iabs(ioto(1,i,4)-ioto(1,i,7))+1)*(iabs(ioto(1,i,5)-ioto(1,i,8))+1)
!		np_new=(iabs(ioto(2,i,3)-ioto(2,i,6))+1)*(iabs(ioto(2,i,4)-ioto(2,i,7))+1)*(iabs(ioto(2,i,5)-ioto(2,i,8))+1)
		node_old=nload(iblk_old)!自己所在的进程号
		node_new=nload(iblk_new)!对方所在的进程号
		if(node_old .ne. node_new)then
			nmpi_np(node_old,node_new)=nmpi_np(node_old,node_new)+np_old
			npatch(node_old,node_new)=1
		endif
	end do

	do i=1,noto
		iblk_old=ioto(2,i,2)	!自己的块序号
		iblk_new=ioto(1,i,2)	!对方的块序号
		np_old=(iabs(ioto(2,i,3)-ioto(2,i,6))+1)*(iabs(ioto(2,i,4)-ioto(2,i,7))+1)*(iabs(ioto(2,i,5)-ioto(2,i,8))+1)
!		np_new=(iabs(ioto(2,i,3)-ioto(2,i,6))+1)*(iabs(ioto(2,i,4)-ioto(2,i,7))+1)*(iabs(ioto(2,i,5)-ioto(2,i,8))+1)
		node_old=nload(iblk_old)!自己所在的进程号
		node_new=nload(iblk_new)!对方所在的进程号
		if(node_old .ne. node_new)then
			nmpi_np(node_old,node_new)=nmpi_np(node_old,node_new)+np_old
			npatch(node_old,node_new)=1
		endif
	end do
    max_distance=0
	max_round=0			!记录一个进程所需要传递的最大次数
	nround_all=0		!记录所有进程总共需要传递的次数
	do i=1,ncpu
		nround=0
		idistance=0
		idistance_max=0
		do j=1,ncpu
			if(npatch(i,j).ne.0)then
				nround=nround+1
				nround_all=nround_all+1
				idistance=abs(i-j)
				if(idistance>ncpu/2)then
				    idistance=abs(idistance-ncpu)
				endif
				if(idistance>idistance_max)idistance_max=idistance
!				write(*,*)i,j,idistance,idistance_max,max_distance
!				pause
			endif
		end do
		if(nround>max_round)max_round=nround
		if(idistance_max>max_distance)max_distance=idistance_max
	end do
	nround_all=nround_all/2

!	open(51,file='res.dat')
	write(60,*)	
	write(60,*)'distance ', max_distance
	write(60,*)'roundall ', nround_all
	write(60,*)'Ratio_max ', point_ratio_max
	write(60,*)'Ratio_min ', point_ratio_min

	write(60,*)
	write(60,*)'Mpi Transfer Matrix:'
	do i=1,ncpu
		write(60,'(<ncpu>i2)')npatch(i,:)			!输出迭代交换的面对面矩阵
	end do
	write(60,*)
	nmpi_np_all=0
	do i=1,ncpu
		nodecpu1=0
		nodecpu2=0
		write(60,*)
		do j=1,ncpu
			if(nmpi_np(i,j) .ne. 0)then
				nmpi_np_all=nmpi_np_all+nmpi_np(i,j)
				nodecpu1=nodecpu1+1
				nodecpu2=nodecpu2+nmpi_np(i,j)
				write(60,'(a4,i5,a8,i5,a9,i10,a7)')'Node',i,' to Node',j,' Transfer',nmpi_np(i,j),' Points'
			endif
		end do
		write(60,'(a4,i5,a16,i5,a20,i10,a7)')'Node',i,' connects with',nodecpu1,' Nodes  Transfers' ,nodecpu2,' Points'
	end do
	write(60,*)
	write(60,'(a17,i15,a7)')'Total Transfered ',nmpi_np_all/2,' Points'

	close(60)


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
open(100,file='wall_bc.txt')
nbcwall=0
	do nblk=1,iblk		! selectibjb and kind
		k=0
		do nbound=1,iblock(nblk,4)
			k=k+1
			call selectibjb(iblock,ibound,iNSAWET,iblk,nblk,k)
			call selectkind(iblock,ibound,iNSAWET,iblk,nblk,k)
			if(ibound(nblk,k,7)<0) k=k+1
if(iNSAWET(nblk,nbound,7).eq.1)then
write(100,*) nblk,nbound,1
nbcwall=nbcwall+1
endif
		end do
    end do
write(100,*) nbcwall
close(100)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   处理准点对点搭接的数据
    do i=1,ninter
        nqua=nquasi(1,i,1)
        nblk=nqua/100
        iface1=mod(nqua,100)/10
        iface2=mod(nqua,10)
        !write(*,*)nblk,iface1,iface2
        iblock(nblk,4)=iblock(nblk,4)+1
       ! nbound=0
       ! do ib=1,iblock(nblk,4)
       !     nbound=nbound+1
       !     if(ibound(nblk,ib,7)<0)nbound=nbound+1
       ! end do
       nbound=ibc(nblk,1)+iblock(nblk,4)
!        write(*,*)nblk,iblock(nblk,4),nbound,'testxxx'
!        pause
        call getindex(iblock,iNSAWET,nquasi,ninter,i,iblk,nblk,iface1,iface2,nbound)
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!搜索坐标了！
    if(ninter > 0) then
        allocate(xyzall_x(npointsall_index));   xyzall_x=0.0
        allocate(xyzall1_x(npointsall_index));  xyzall1_x=0.0
        allocate(xyzall_y(npointsall_index));   xyzall_y=0.0
        allocate(xyzall1_y(npointsall_index));  xyzall1_y=0.0
        allocate(xyzall_z(npointsall_index));   xyzall_z=0.0
        allocate(xyzall1_z(npointsall_index));  xyzall1_z=0.0
        
        call readpoints(xyzall_x,xyzall_y,xyzall_z,xyzall1_x,xyzall1_y,xyzall1_z,ireal,iblk,iblock,npointsall_index,iblk_index,filegrd)
    endif
    
allocate(xyzpatch(ninter,2,2,3),xyzpatch1(ninter,2,2,3))
do i=1,ninter
    nblk=nquasi(1,i,7)     !块序号
    nbound=nquasi(1,i,9)    !当前边界在iNSAWET中的位置
    call readdata(xyzpatch,xyzpatch1,ninter,nblk,nquasi,filegrd,i,iblock,ireal,iblk, &
    xyzall_x,xyzall_y,xyzall_z,xyzall1_x,xyzall1_y,xyzall1_z,iblk_index,npointsall_index)
end do
if(ireal==4)xyzpatch1=xyzpatch
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(5,file='patchsurface.dat')
do i=1,ninter
    write(5,*)'zone i= 2 j= 2'
    do j1=1,2
        do i1=1,2
            write(5,'(3e20.10)')xyzpatch1(i,i1,j1,1:3)
        end do
    end do
end do
close(5)
!pause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!分情况讨论应该是lcross/lreva/lrevb
do i=1,ninter
    nblk=nquasi(1,i,7)     !块序号
    nbound=nquasi(1,i,9)    !当前边界在iNSAWET中的位置 
    
    iNSAWET(nblk,nbound,12)=-100
    iNSAWET(nblk,nbound,13)=-100
    iNSAWET(nblk,nbound,14)=-100
    
    if(mod(i,2)==1)i1=i+1
    if(mod(i,2)==0)i1=i-1
    dmax1=dis(xyzpatch1(i,1,1,1:3),xyzpatch1(i,2,2,1:3))
    dmax2=dis(xyzpatch1(i,1,2,1:3),xyzpatch1(i,2,1,1:3))
    dt=max(dmax1,dmax2)*1.0e-5
!    write(*,*)dt

    
    if(dis(xyzpatch1(i,1,1,1:3),xyzpatch1(i1,1,1,1:3))<dt)then
         if    (dis(xyzpatch1(i,1,2,1:3),xyzpatch1(i1,1,2,1:3))<dt)then
               lcross=0; lreva=0; lrevb=0
         elseif(dis(xyzpatch1(i,1,2,1:3),xyzpatch1(i1,2,1,1:3))<dt)then
               lcross=1; lreva=0; lrevb=0
         endif
    elseif(dis(xyzpatch1(i,1,1,1:3),xyzpatch1(i1,2,1,1:3))<dt)then
         if    (dis(xyzpatch1(i,2,1,1:3),xyzpatch1(i1,1,1,1:3))<dt)then
               lcross=0; lreva=1; lrevb=0
         elseif(dis(xyzpatch1(i,2,1,1:3),xyzpatch1(i1,2,2,1:3))<dt)then
               lcross=1; lreva=0; lrevb=1
         endif
    elseif(dis(xyzpatch1(i,1,1,1:3),xyzpatch1(i1,1,2,1:3))<dt)then
         if    (dis(xyzpatch1(i,2,1,1:3),xyzpatch1(i1,1,1,1:3))<dt)then
               lcross=1; lreva=1; lrevb=0
         elseif(dis(xyzpatch1(i,2,1,1:3),xyzpatch1(i1,2,2,1:3))<dt)then
               lcross=0; lreva=0; lrevb=1
         endif
    elseif(dis(xyzpatch1(i,1,1,1:3),xyzpatch1(i1,2,2,1:3))<dt)then
         if    (dis(xyzpatch1(i,2,1,1:3),xyzpatch1(i1,2,1,1:3))<dt)then
               lcross=1; lreva=1; lrevb=1
         elseif(dis(xyzpatch1(i,2,1,1:3),xyzpatch1(i1,2,2,1:3))<dt)then
               lcross=0; lreva=1; lrevb=1
         endif
    else
        write(*,*)'Patch Cross Wrong! stop!'
    endif
    iNSAWET(nblk,nbound,12)=lcross
    iNSAWET(nblk,nbound,13)=lreva
    iNSAWET(nblk,nbound,14)=lrevb
    nquasi(1,i,10)=lcross
end do



    
    
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    
    do i=1,ninter,2
        nblk1=nquasi(1,i,7)
        xb1=nquasi(1,i,8)
        nbound1=nquasi(1,i,9)
        lcross1=nquasi(1,i,10)

        npra1= iNSAWET(nblk1,nbound1,3)-iNSAWET(nblk1,nbound1,2)+1
        nprb1= iNSAWET(nblk1,nbound1,5)-iNSAWET(nblk1,nbound1,4)+1
        
 
        nblk2=nquasi(1,i+1,7)
        xb2=nquasi(1,i+1,8)
        nbound2=nquasi(1,i+1,9)
        lcross2=nquasi(1,i,10)
        npra2= iNSAWET(nblk2,nbound2,3)-iNSAWET(nblk2,nbound2,2)+1
        nprb2= iNSAWET(nblk2,nbound2,5)-iNSAWET(nblk2,nbound2,4)+1
        
        iNSAWET(nblk1,nbound1,10)=nblk2
        iNSAWET(nblk1,nbound1,11)=xb2
        iNSAWET(nblk1,nbound1,15)=npra2
        iNSAWET(nblk1,nbound1,16)=nprb2
        if(lcross1==1)then
            iNSAWET(nblk1,nbound1,15)=nprb2
            iNSAWET(nblk1,nbound1,16)=npra2            
        endif
        
        
        iNSAWET(nblk2,nbound2,10)=nblk1
        iNSAWET(nblk2,nbound2,11)=xb1
        iNSAWET(nblk2,nbound2,15)=npra1
        iNSAWET(nblk2,nbound2,16)=nprb1 
        if(lcross2==1)then
            iNSAWET(nblk2,nbound2,15)=nprb1
            iNSAWET(nblk2,nbound2,16)=npra1            
        endif
        
        if(npra1*nprb1 > npra2*nprb2)then
            iNSAWET(nblk2,nbound2,8)=3
        else
            iNSAWET(nblk1,nbound1,8)=3
        endif
    end do
    
    
    
!!!!!!!!output NSAWET gridfile!!!!!!!!!!!!!!!!!!!!!!!!
!pause
	write(*,*)
	open(3,file=filegrd,form='binary') !,Access = 'Direct',recl=4)
	read(3)itemp
	read(3)iblk
	call readword(iblk)
	read(3)itemp
	read(3)itemp
	do j=1,iblk
	do i=1,3
	read(3)itemp
	!write(*,*)itemp
	end do
	end do
	read(3)itemp



	do i=1,iblk
		ni=iblock(i,1)
		nj=iblock(i,2)
		nk=iblock(i,3)
		nprocess=idim3(i,1)
		read(3)itemp
		call readword(itemp)
!		write(*,*)itemp
		ireal=itemp/ni/nj/nk/3
		if(ireal==4)then
			write(*,*)"Block",i," Single precision!"
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
		endif
		if(ireal==8)then
			write(*,*)"Block",i," Double precision!"
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
		endif

		read(3)itemp
		call readword(itemp)

!		pause
!		read(3)(((xyz(i1,j1,k1,1),i1=1,ni),j1=1,nj),k1=1,nk)
!		read(3)(((xyz(i1,j1,k1,2),i1=1,ni),j1=1,nj),k1=1,nk)
!		read(3)(((xyz(i1,j1,k1,3),i1=1,ni),j1=1,nj),k1=1,nk)
		if(i<10)then
			write(outfile,'(i1)')i
			outfile=trim('00')//trim(outfile)
		elseif(i<100)then
			write(outfile,'(i2)')i
			outfile=trim('0')//trim(outfile)
		elseif(i<1000)then
			write(outfile,'(i3)')i
		elseif(i<10000)then
			write(outfile,'(i4)')i
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
		if(ibnum<6)then
			write(*,*)"Error: Block boundary less than 6!"
			!stop
		endif
		write(4,'(a1,4i8,f11.2)')'#',i,ibnum,itemp,itemp,scale
		write(4,'(a170)')'#   blkbd_id  idxa0   idxa1   idxb0   idxb1   isuf    kind  kdsub   lsfpro  iblkex iblkbd_ex lcross  lreva   lrevb    npra    nprb     v#1     v#2     v#3     v#4   lhole'
		k=0
        
!        write(*,*)'test!',iNSAWET(3,6,1:7)
		do j=1,ibnum
			k=k+1
			if(iNSAWET(i,k,7)<=0)then
				write(*,*)'Block',i,' Boundary Condition Error!'
				!stop
            endif
            !write(*,*)i,k,iNSAWET(i,k,1)
			write(4,'(a1,21i8)')'#',(iNSAWET(i,k,ii),ii=1,21)
            !write(*,'(a1,21i8)')'#',(iNSAWET(i,k,ii),ii=1,21)
			if(ibound(i,k,7)<0) k=k+1
		end do
		write(4,'(a50)')'#---------name for surface patch file-------------'
		write(4,'(a50)')'#---------name for quasi patch file---------------'
		write(4,'(a50)')'#-------------------------------------------------'
!		write(4,*)'  zone i= ',ni,'  j= ',nj,'  k= ',nk
		if(nprocess<10)then
			write(comm,'(i1)')nprocess
		elseif(nprocess<100)then
			write(comm,'(i2)')nprocess
		elseif(nprocess<1000)then
			write(comm,'(i3)')nprocess
		elseif(nprocess<10000)then
			write(comm,'(i4)')nprocess
		endif
		comm=trim('zone T="')//trim(comm)//trim('"')
		write(4,'(a20,a4,i8,a4,i8,a4,i8)')comm,' i= ',ni,' j= ',nj,' k= ',nk

		if(ireal==4)write(4,'(3e20.10)')((((xyz(i1,j1,k1,ii),ii=1,3),i1=1,ni),j1=1,nj),k1=1,nk)
		if(ireal==8)write(4,'(3e20.10)')((((xyz1(i1,j1,k1,ii),ii=1,3),i1=1,ni),j1=1,nj),k1=1,nk)
		close(4)
		if(ireal==4)deallocate(xyz)
		if(ireal==8)deallocate(xyz1)
	end do
	close(3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	deallocate(iblock,ibound,iNSAWET)

end program

subroutine readword(i)
  integer*4 i,j
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

subroutine selectlcross(iblock,ibound,iNSAWET,iblk,nblk,nbound)
	integer iblock(iblk,4),ibound(iblk,1000,10),iNSAWET(iblk,1000,30)
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
	integer iblock(iblk,4),ibound(iblk,1000,10),iNSAWET(iblk,1000,30)
	integer ipnew(6)
	
	iNSAWET(nblk,nbound,10)=ibound(nblk,nbound+1,7)
	iblknew=iNSAWET(nblk,nbound,10)
	do i=1,6
		ipnew(i)=ibound(nblk,nbound+1,i)
	end do
!	write(*,*)
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
	integer iblock(iblk,4),ibound(iblk,1000,10),iNSAWET(iblk,1000,30)
    
	if(iNSAWET(nblk,nbound,7)<4) then
		do ii=8,21
			iNSAWET(nblk,nbound,ii)=0
		end do
    endif	
    
	if(ibound(nblk,nbound,7)==2004) iNSAWET(nblk,nbound,7)=1		!wall
	if(ibound(nblk,nbound,7)==1001) iNSAWET(nblk,nbound,7)=2		!symmetry
	if(ibound(nblk,nbound,7)==1000) iNSAWET(nblk,nbound,7)=3		!farfield
	if(ibound(nblk,nbound,7)==1003) iNSAWET(nblk,nbound,7)=3		!farfield
!    if(ibound(nblk,nbound,7)==1005) then
!            write(*,*)'Invisicd wall 1005 found!'
!            write(*,*)'It is not supported by NSAWET'
!            stop
!    endif
    
!------------------------------------------------------------------------------
	if(ibound(nblk,nbound,7)==1006)then   !inviscid wall 
        iNSAWET(nblk,nbound,7)=1	
        iNSAWET(nblk,nbound,8)=4
        !iNSAWET(nblk,nbound,18)=1
    endif
    
    if(ibound(nblk,nbound,7)==2002)then  !static pressure
        iNSAWET(nblk,nbound,7)=3		
        iNSAWET(nblk,nbound,8)=1		
    endif
    if(ibound(nblk,nbound,7)==1008)then  !tunnel inflow
        iNSAWET(nblk,nbound,7)=3		
        iNSAWET(nblk,nbound,8)=4		
    endif
    if(ibound(nblk,nbound,7)==1002)then  !extrapolation
        iNSAWET(nblk,nbound,7)=3		
        iNSAWET(nblk,nbound,8)=3		
    endif
    if(ibound(nblk,nbound,7)==2003)then  !total parameters
        iNSAWET(nblk,nbound,7)=3		
        iNSAWET(nblk,nbound,8)=5		
    endif
    
    if(ibound(nblk,nbound,7)==1011)then  !singular - half plane
        iNSAWET(nblk,nbound,7)=5		
        iNSAWET(nblk,nbound,8)=1		
    endif
    if(ibound(nblk,nbound,7)==1012)then  !singular - full plane
        iNSAWET(nblk,nbound,7)=5		
        iNSAWET(nblk,nbound,8)=2		
    endif
    if(ibound(nblk,nbound,7)==1013)then  !singular - partial plane
        iNSAWET(nblk,nbound,7)=5				
    endif

    if(ibound(nblk,nbound,7)==1005)then   !span periodic
        iNSAWET(nblk,nbound,7)=4		
        iNSAWET(nblk,nbound,10)=nblk		
		do ii=1,iblock(nblk,4)
        if(ii.eq.nbound) then
            cycle
        elseif(ibound(nblk,ii,7).eq.1005)then
            iNSAWET(nblk,nbound,11)=ii
            exit
        endif
        enddo
    endif
    
     if(ibound(nblk,nbound,7)==2102)then   !moving wall
        iNSAWET(nblk,nbound,7)=1		
        iNSAWET(nblk,nbound,8)=6		
    endif   
!------------------------------------------------------------------------------
	if(ibound(nblk,nbound,7)< 0) iNSAWET(nblk,nbound,7)=4		!1to1patch
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
	integer iblock(iblk,4),ibound(iblk,1000,10),iNSAWET(iblk,1000,30)
		
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
	integer iblock(iblk,4),ibound(iblk,1000,10),iNSAWET(iblk,1000,30)

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
    
    
subroutine getindex(iblock,iNSAWET,nquasi,ninter,inter,iblk,nblk,iface1,iface2,nbound)
    integer iblock(iblk,4),iNSAWET(iblk,1000,30),nquasi(2,ninter,10)
    
!    nbound= 0

    
    idim  = iblock(nblk,1)
    jdim  = iblock(nblk,2)
    kdim  = iblock(nblk,3)
    iNSAWET(nblk,nbound,1)=iblock(nblk,4)

    !write(*,*)nblk,nbound,idim,jdim,kdim
    !write(*,*)iface1,iface2
    !pause 'zyf00'
    
    if(iface1==1)then  !
        if(nquasi(1,inter,2)==0)nquasi(1,inter,2)=1
        if(nquasi(1,inter,3)==0)nquasi(1,inter,3)=jdim
        if(nquasi(1,inter,4)==0)nquasi(1,inter,4)=1
        if(nquasi(1,inter,5)==0)nquasi(1,inter,5)=kdim
        iNSAWET(nblk,nbound,2)=nquasi(1,inter,2)
        iNSAWET(nblk,nbound,3)=nquasi(1,inter,3)-1
        iNSAWET(nblk,nbound,4)=nquasi(1,inter,4)
        iNSAWET(nblk,nbound,5)=nquasi(1,inter,5)-1
        if(iface2==1)iNSAWET(nblk,nbound,6)=1
        if(iface2==2)iNSAWET(nblk,nbound,6)=2
    elseif(iface1==2)then
        if(nquasi(1,inter,2)==0)nquasi(1,inter,2)=1
        if(nquasi(1,inter,3)==0)nquasi(1,inter,3)=kdim
        if(nquasi(1,inter,4)==0)nquasi(1,inter,4)=1
        if(nquasi(1,inter,5)==0)nquasi(1,inter,5)=idim     
        iNSAWET(nblk,nbound,4)=nquasi(1,inter,2)
        iNSAWET(nblk,nbound,5)=nquasi(1,inter,3)-1
        iNSAWET(nblk,nbound,2)=nquasi(1,inter,4)
        iNSAWET(nblk,nbound,3)=nquasi(1,inter,5)-1    
        if(iface2==1)iNSAWET(nblk,nbound,6)=3
        if(iface2==2)iNSAWET(nblk,nbound,6)=4
       ! write(*,*)nblk,idim,jdim,kdim
       ! write(*,*)nblk,nquasi(1,inter,2:5)
       ! pause 'zyf'
    elseif(iface1==3)then
        if(nquasi(1,inter,2)==0)nquasi(1,inter,2)=1
        if(nquasi(1,inter,3)==0)nquasi(1,inter,3)=jdim
        if(nquasi(1,inter,4)==0)nquasi(1,inter,4)=1
        if(nquasi(1,inter,5)==0)nquasi(1,inter,5)=idim     
        iNSAWET(nblk,nbound,4)=nquasi(1,inter,2)
        iNSAWET(nblk,nbound,5)=nquasi(1,inter,3)-1
        iNSAWET(nblk,nbound,2)=nquasi(1,inter,4)
        iNSAWET(nblk,nbound,3)=nquasi(1,inter,5)-1     
        if(iface2==1)iNSAWET(nblk,nbound,6)=5
        if(iface2==2)iNSAWET(nblk,nbound,6)=6
    else
        stop 'iface1 wrong!'
    endif
    
    iNSAWET(nblk,nbound,7)=4
    iNSAWET(nblk,nbound,8)=2
    
    nquasi(1,inter,7) = nblk            !当前的块序号
    nquasi(1,inter,8) = iblock(nblk,4)  !当前的边界编号
    nquasi(1,inter,9) = nbound          !当前的存储编号
!    write(*,*)nblk,'nblk',iNSAWET(nblk,nbound,1:6)
    end subroutine
    
    function dis(a,b)
    real*8 a(3),b(3)
    real*8   dis
    dis=sqrt((a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2)
    return
    end function