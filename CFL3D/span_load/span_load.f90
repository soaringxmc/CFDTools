program Spanload
	parameter(mi=10000)
	real,allocatable:: xyz(:,:,:),SD(:,:,:,:),ppa(:,:,:,:)
	real,allocatable:: xx(:,:,:),yy(:,:,:),zz(:,:,:),Cp(:,:,:),grad(:,:,:)
	real,allocatable:: xc(:,:,:),yc(:,:,:),zc(:,:,:),Cpc(:,:,:),span(:,:,:)
	real,allocatable:: xtemp(:)
	real,allocatable:: xlocation(:),clall(:,:),clall2(:,:),chord(:,:)
	real,allocatable:: x(:),y(:),z(:)
	integer,allocatable:: ixyz(:,:,:),iseg(:,:),ilistzone(:)
	integer locaX,locaY,locaZ,locaCp
	character*100 char,filename
	logical lc
    
    pi=4.*atan(1.)

!!!!!!!!!!!!!!!!!!!!!!!!!读入控制文件，给定变量列数位置!!!!!!!!!!!!!!!!!
	open(5,file='SpanLoad.txt',status='old',action='read')
	read(5,*)
	read(5,*)filename
	read(5,*)
	read(5,*)
	read(5,*)locaX,locaY,locaZ,locaI,locaJ,locaK,locaCp
	mvari=max(locaX,locaY,locaZ,locaI,locaJ,locaK,locaCp)      !不用implicit none，，，有毒，，，
	read(5,*)
	read(5,*)ijkformat
	read(5,*)
	read(5,*)AoA
    aoa=aoa*pi/180.
    read(5,*)
	read(5,*)idirection
	read(5,*)
	read(5,*)ilist
	if(ilist .ne. 0)then
		read(5,*)ilistall
		allocate(ilistzone(ilistall))
		read(5,*)(ilistzone(1:ilistall))
	endif
	read(5,*)
	read(5,*)ilocation
	allocate(xtemp(mvari),xlocation(ilocation),clall(ilocation,3),clall2(ilocation,3),xyz(ilocation,mi,4),ixyz(ilocation,mi,4),iseg(mi,2))
    allocate(chord(ilocation,3))
	do i=1,ilocation
		read(5,*)xlocation(i)
	end do
	close(5)
	clall=0.0
	clall2=0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	open(6,file=filename,status='old',action='read')
	read(6,*)
!	open(7,file='test.dat')
	open(8,file='SpanLoad_slice.dat')
	write(8,*)'variables = x y z cp'
!	open(9,file='clall.dat')
	open(10,file='SpanLoad.dat')
        chord(:,1)=-1e20
        chord(:,2)=1e20
do izone=1,99999900


	dl=1.
	if(ijkformat==2)then
		read(6,*,end=1010)char,char,ni,char,nj,char,nk
		if(ni==1)then
			ni=nj
			nj=nk
			iface=1
			dl=1.
		elseif(nj==1)then
			nj=nk
			dl=1.
			iface=2
		elseif(nk==1)then
			dl=-1.
			iface=3
		else
			write(*,*)'Not a surface. May be a volume.'
			stop
		endif
	elseif(ijkformat==1)then                                !NSAWET format   No iface    dl=1.0
		read(6,*)char,char,ni,char,nj
	else
		write(*,*)'Stop. Zone format wrong!'
		stop
	endif
	nk=1
	ni1=ni-1            !cell number
	nj1=nj-1
	nk1=nk-1
	!write(7,*)'zone i= ',ni1,' j= ',nj1,' k= ',1
	allocate(xx(ni,nj,nk),yy(ni,nj,nk),zz(ni,nj,nk),cp(ni,nj,nk),sd(3,ni,nj,nk),grad(ni,nj,nk))
	allocate(xc(ni,nj,nk),yc(ni,nj,nk),zc(ni,nj,nk),cpc(ni,nj,nk),ppa(3,ni,nj,nk))
    
	!读入坐标点
	do k=1,nk
		do j=1,nj
			do i=1,ni
				read(6,*)xtemp(1:mvari)
				xx(i,j,k)=xtemp(locaX)
				if(idirection==2)then                   
					yy(i,j,k)=xtemp(locaY)
					zz(i,j,k)=xtemp(locaZ)
				elseif(idirection==3)then               !保证展向是y
					yy(i,j,k)=xtemp(locaZ)
					zz(i,j,k)=-xtemp(locaY)
				else
					write(*,*)'Wrong direction!'
					stop
				endif
!				if(i.eq.1 .and. j.eq.1 .and. k.eq.1)write(*,*)izone,iface
				if(iface==1)then			
					if(xtemp(locaI)<1.1)then
						dl=-1.0
					endif
				elseif(iface==2)then
					if(xtemp(locaJ)<1.1)then
						dl=-1.0
					endif
				elseif(iface==3)then
					if(xtemp(locaK)<1.1)then
						dl=1.0
					endif
				endif					

				cp(i,j,k)=xtemp(locaCp)                   !节点Cp
			end do
		end do
	enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(ilist .ne. 0)then
		do il=1,ilistall
			i=ilistzone(il)
			if(i==izone)goto 105
		end do
		goto 106
	endif
				
105 continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	计算面积
	call surface(xx,yy,zz,sd,grad,cp,xc,yc,zc,cpc,ppa,ni,nj,nk,dl)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	

!	do k=1,1
!		do j=1,nj1
!			do i=1,ni1
!				write(7,'(7e20.10)')xc(i,j,k),yc(i,j,k),zc(i,j,k),cpc(i,j,k),ppa(1:3,i,j,k)
!			end do
!		end do
!	end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	求交点
	allocate(x(2),y(2),z(2))
	k=1
	ixyz=-100
	xyz=0.0
	do ilocal=1,ilocation

		itemp=0
		do j=1,nj
			do i=1,ni-1
				y(1)=yy(i  ,j,k)
				y(2)=yy(i+1,j,k)
				yin=xlocation(ilocal)
				if( (y(2)>=yin .and. y(1)<yin) .or. &
					(y(1)>=yin .and. y(2)<yin)) then
					itemp=itemp+1
					x(1)=xx(i  ,j,k)
					x(2)=xx(i+1,j,k)
					z(1)=zz(i  ,j,k)
					z(2)=zz(i+1,j,k)
					c1=cp(i,j,k)
					c2=cp(i+1,j,k)
					c3=(c2-c1)/(y(2)-y(1))*(yin-y(1))+c1
					xyz(ilocal,itemp,4)=c3              !Cp
					call cz_point_z(x,z,y,yin,xout,zout)
					xyz(ilocal,itemp,1)=xout
					xyz(ilocal,itemp,2)=yin
					xyz(ilocal,itemp,3)=zout
					ixyz(ilocal,itemp,1)=i
					ixyz(ilocal,itemp,2)=j-1
					ixyz(ilocal,itemp,3)=i
					ixyz(ilocal,itemp,4)=j
				endif
			end do
        end do


		do i=1,ni
			do j=1,nj-1
				y(1)=yy(i,j  ,k)
				y(2)=yy(i,j+1,k)
				yin=xlocation(ilocal)
				if( (y(2)>=yin .and. y(1)<yin) .or. &
					(y(1)>=yin .and. y(2)<yin)) then
					itemp=itemp+1
					x(1)=xx(i,j  ,k)
					x(2)=xx(i,j+1,k)
					z(1)=zz(i,j  ,k)
					z(2)=zz(i,j+1,k)
					c1=cp(i,j,k)
					c2=cp(i,j+1,k)
					c3=(c2-c1)/(y(2)-y(1))*(yin-y(1))+c1
					xyz(ilocal,itemp,4)=c3
					call cz_point_z(x,z,y,yin,xout,zout)
					xyz(ilocal,itemp,1)=xout
					xyz(ilocal,itemp,2)=yin
					xyz(ilocal,itemp,3)=zout
					ixyz(ilocal,itemp,1)=i-1
					ixyz(ilocal,itemp,2)=j
					ixyz(ilocal,itemp,3)=i
					ixyz(ilocal,itemp,4)=j
				endif
			end do
		end do
		!write(*,*)izone,itemp
		if(itemp>0)then
			isegall=0
			do i=1,itemp
				i1=ixyz(ilocal,i,1)
				j1=ixyz(ilocal,i,2)
				i2=ixyz(ilocal,i,3)
				j2=ixyz(ilocal,i,4)
				do j=i+1,itemp
					lc=.false.
					i3=ixyz(ilocal,j,1)
					j3=ixyz(ilocal,j,2)
					i4=ixyz(ilocal,j,3)
					j4=ixyz(ilocal,j,4)
					if((i1==i3 .and. j1==j3) .or. (i1==i4 .and. j1==j4))then
						ic=i1
						jc=j1
						lc=.true.
					elseif((i2==i3 .and. j2==j3) .or. (i2==i4 .and. j2==j4))then
						ic=i2
						jc=j2
						lc=.true.
					endif
					if(lc)then
						isegall=isegall+1
						iseg(isegall,1)=i
						iseg(isegall,2)=j
						x1=xyz(ilocal,i,1)
						y1=xyz(ilocal,i,2)
						z1=xyz(ilocal,i,3)
						p1=xyz(ilocal,i,4)
						x2=xyz(ilocal,j,1)
						y2=xyz(ilocal,j,2)
						z2=xyz(ilocal,j,3)
						p2=xyz(ilocal,j,4)
					!	write(*,*)i,j
					!	write(*,'(6f12.5)')x1,y1,z1,x2,y2,z2
                        if(x1>chord(ilocal,1))chord(ilocal,1)=x1
                        if(x2>chord(ilocal,1))chord(ilocal,1)=x2
                        if(x1<chord(ilocal,2))chord(ilocal,2)=x1
                        if(x2<chord(ilocal,2))chord(ilocal,2)=x2
                        
						xlength=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
					!	write(*,*)xlength

						 R1X=x1-x2	!XX(I1,J,K)-XX(I,J1,K)
						 R1Y=1.0		!YY(I1,J,K)-YY(I,J1,K)
						 R1Z=z1-z2		!ZZ(I1,J,K)-ZZ(I,J1,K)
						 R2X=x2-x1			!XX(I1,J1,K)-XX(I,J,K)
						 R2Y=1.0				!YY(I1,J1,K)-YY(I,J,K)
						 R2Z=z2-z1					!ZZ(I1,J1,K)-ZZ(I,J,K)


						 xd=(R1Y*R2Z-R1Z*R2Y)*0.5
						 yd=(R1Z*R2X-R1X*R2Z)*0.5
						 zd=(R1X*R2Y-R1Y*R2X)*0.5

						!xd=abs(x1-x2) !/(xlength+1.0e-20)
						!zd=abs(z1-z2) !/(xlength+1.0e-20)
						!yd=abs(y1-y2) !/(xlength+1.0e-20)
					!	write(*,'(6f12.5)')xd,yd,zd
						xd=sign(xd,sd(1,ic,jc,k))
						yd=sign(yd,sd(2,ic,jc,k))
						zd=sign(zd,sd(3,ic,jc,k))
					!	write(*,'(6f12.5)')xd,yd,zd,sd(1,ic,jc,k),sd(2,ic,jc,k),sd(3,ic,jc,k)
					!	pause
						clall(ilocal,1)=clall(ilocal,1)+ppa(1,ic,jc,k)*xlength
						clall(ilocal,2)=clall(ilocal,2)+ppa(2,ic,jc,k)*xlength
						clall(ilocal,3)=clall(ilocal,3)+ppa(3,ic,jc,k)*xlength

						clall2(ilocal,1)=clall2(ilocal,1)+(p1+p2)/2.0*xd !*xlength
						clall2(ilocal,2)=clall2(ilocal,2)+(p1+p2)/2.0*yd !*xlength
						clall2(ilocal,3)=clall2(ilocal,3)+(p1+p2)/2.0*zd !*xlength

					endif
				end do
            end do

				write(8,'(a12,i7,a10,i7,a40)')'zone Nodes=',itemp,' Elements=',isegall, &
						' ZONETYPE=FELineSeg DATAPACKING=POINT'
				if(idirection==2)then
					do i=1,itemp
						write(8,'(4e20.10)')xyz(ilocal,i,1:4) !,ixyz(ilocal,i,1:4)
					end do
				elseif(idirection==3)then
					do i=1,itemp
						write(8,'(4e20.10)')xyz(ilocal,i,1),-xyz(ilocal,i,3),xyz(ilocal,i,2),xyz(ilocal,i,4)
					end do
                endif
				do i=1,isegall
					write(8,'(2i8)')iseg(i,1),iseg(i,2)
				end do
        endif
   
        
	end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	deallocate(x,y,z)

106 continue	

	deallocate(xx,yy,zz,cp,sd,grad,xc,yc,zc,cpc,ppa)
	
end do

1010 continue

     
         do ilocal=1,ilocation
                chord(ilocal,3)=chord(ilocal,1)-chord(ilocal,2)
               ! write(*,*)chord(ilocal,1),chord(ilocal,2),chord(ilocal,3)
    end do
!	write(9,*)'variables = Location Load'
!	do i=1,ilocation
!		write(9,'(f15.5,1f20.8)')xlocation(i),clall(i,3:3)
!	end do
!	close(9)

	write(10,*)'variables = Location Load CY CX CL CD Chord'
	do i=1,ilocation
		if(idirection==2)then
            cy=clall2(i,3)/chord(i,3)
            cx=clall2(i,1)/chord(i,3)
            cl=cy*cos(aoa)-cx*sin(aoa)
            cd=cy*sin(aoa)+cx*cos(aoa)
			write(10,'(f15.5,10f13.8)')abs(xlocation(i)),clall2(i,3:3),cy,cx,cl,cd,chord(i,3)
        elseif(idirection==3)then
            cy=-clall2(i,3)/chord(i,3)
            cx= clall2(i,1)/chord(i,3)
            cl=cy*cos(aoa)-cx*sin(aoa)
            cd=cy*sin(aoa)+cx*cos(aoa)
			write(10,'(f15.5,10f13.8)')abs(xlocation(i)),-clall2(i,3:3),cy,cx,cl,cd,chord(i,3)
		endif
	end do
	close(10)

end program


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cz_point_z(x,y,z,zin,xout,yout)   !cz_point_z(x,z,y,yin,xout,zout)
	real x(2),y(2),z(2),zin,xout,yout
	xout=x(1)+(x(2)-x(1))*(zin-z(1))/(z(2)-z(1))
	yout=y(1)+(y(2)-y(1))*(zin-z(1))/(z(2)-z(1))
	return
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE surface(xx,yy,zz,sd,grad,cp,xc,yc,zc,cpc,ppa,ni,nj,nk,dl)   !网格单元几何信息
	real xx(ni,nj,nk),yy(ni,nj,nk),zz(ni,nj,nk)
	real sd(3,ni,nj,nk),grad(ni,nj,nk),cp(ni,nj,nk),ppa(3,ni,nj,nk)
	real xc(ni,nj,nk),yc(ni,nj,nk),zc(ni,nj,nk),cpc(ni,nj,nk),span(ni,nj,nk)

	ni1=ni-1
	nj1=nj-1
	nk1=nk-1

 
		k=1
        DO I=1,NI1
           I1=I+1
           DO J=1,NJ1
              J1=J+1
			  xc(I,J,K)=(xx(I1,J1,K)+xx(I1,J,K)+xx(I,J,K)+xx(I,J1,K))*0.25
			  yc(I,J,K)=(yy(I1,J1,K)+yy(I1,J,K)+yy(I,J,K)+yy(I,J1,K))*0.25   
			  zc(I,J,K)=(zz(I1,J1,K)+zz(I1,J,K)+zz(I,J,K)+zz(I,J1,K))*0.25
			  cpc(I,J,K)=(cp(I1,J1,K)+cp(I1,J,K)+cp(I,J,K)+cp(I,J1,K))*0.25
	        ENDDO
	     ENDDO

     ! dl=1.
	  xv=xc(2,1,1)-xc(1,1,1)
	  yv=yc(2,1,1)-yc(1,1,1)
	  zv=zc(2,1,1)-zc(1,1,1)

        DO J=1,NJ1
           J1=J+1
           DO I=1,NI1
              I1=I+1
              DO K=1,1
                 R1X=XX(I1,J,K)-XX(I,J1,K)
                 R1Y=YY(I1,J,K)-YY(I,J1,K)
                 R1Z=ZZ(I1,J,K)-ZZ(I,J1,K)
                 R2X=XX(I1,J1,K)-XX(I,J,K)
                 R2Y=YY(I1,J1,K)-YY(I,J,K)
                 R2Z=ZZ(I1,J1,K)-ZZ(I,J,K)
	          ! if(i.eq.1.and.j.eq.1.and.k.eq.1) then
               !    xsv=(R1Y*R2Z-R1Z*R2Y)*0.5
                !   ysv=(R1Z*R2X-R1X*R2Z)*0.5
                 !  zsv=(R1X*R2Y-R1Y*R2X)*0.5
                  ! if(xv*xsv+yv*ysv+zv*zsv.lt.0) dl=-1.
	               
				  ! write(*,*)xv*xsv+yv*ysv+zv*zsv
				  ! pause
				  ! if(dl < -0.5) print*,'Left hand'
	           !endif
                 SD(1,I,J,K)=(R1Y*R2Z-R1Z*R2Y)*0.5*dl
                 SD(2,I,J,K)=(R1Z*R2X-R1X*R2Z)*0.5*dl
                 SD(3,I,J,K)=(R1X*R2Y-R1Y*R2X)*0.5*dl
	        enddo
	     enddo
	  enddo

        do i=1,ni1
	     do j=1,nj1
	        do k=1,1
	              grad(i,j,k)=SQRT(SD(1,I,J,K)**2+SD(2,I,J,K)**2+SD(3,I,J,K)**2)            !面积
				  
				  sd(1,i,j,k)=sd(1,i,j,k)/grad(i,j,k)       !方向
				  sd(2,i,j,k)=sd(2,i,j,k)/grad(i,j,k)
				  sd(3,i,j,k)=sd(3,i,j,k)/grad(i,j,k)

				  ppa(1,i,j,k)=cpc(i,j,k)*sd(1,i,j,k)           !P*方向
				  ppa(2,i,j,k)=cpc(i,j,k)*sd(2,i,j,k)
				  ppa(3,i,j,k)=cpc(i,j,k)*sd(3,i,j,k)

	           enddo
	        enddo
	     enddo

end subroutine
