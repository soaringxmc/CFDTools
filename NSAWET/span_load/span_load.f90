program Spanload
	parameter(mi=10000)
	real,allocatable:: xyz(:,:,:),SD(:,:,:)
	real,allocatable:: xx(:,:),yy(:,:),zz(:,:),Cp(:,:),grad(:,:)
	real,allocatable:: xc(:,:),yc(:,:),zc(:,:),Cpc(:,:),span(:,:),cpc_unused(:,:)
	real,allocatable:: xtemp(:)
	real,allocatable:: zlocation(:),clall2(:,:),chord(:,:)
	real,allocatable:: x(:),y(:),z(:)
	integer,allocatable:: ixyz(:,:,:),iseg(:,:),ilistzone(:)
	integer locaX,locaY,locaZ,locaCp
	character*100 char,filename
	logical lc,alive,stat
    
    pi=4.*atan(1.)

!!!!!!!!!!!!!!!!!!!!!!!!!读入控制文件，给定变量列数位置!!!!!!!!!!!!!!!!!
	open(5,file='SpanLoad.txt',status='old',action='read')
	read(5,*)
	read(5,*)AoA
    aoa=aoa*pi/180.
    read(5,*)
	read(5,*)ilocation
	allocate(xtemp(mvari),zlocation(ilocation),clall2(ilocation,4),xyz(ilocation,mi,4),ixyz(ilocation,mi,4),iseg(mi,2))
    allocate(chord(ilocation,3))
	do i=1,ilocation
		read(5,*)zlocation(i)
	end do
	close(5)
	clall2=0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
	open(8,file='SpanLoad_slice.dat')
	write(8,*)'variables = x y z cp'

    chord(:,1)=-1e20
    chord(:,2)=1e20
        
do ifile=0,1000
    
    call getfilename(ifile,7,'d3ns_Id',4,'.cpa',filename)
    Inquire(file=filename,exist=alive)
    if(alive.eq.0) cycle
    open(6, file=filename, action='read')
    read(6, *)
    
    do izone = 1,100
	!读入坐标点
        
        if(eof(6)) then
            close(6)
            exit
        endif
        
        write(*,*) trim(adjustl(filename)),izone
        
        read(6, *) char,char,char,char,ni,char,nj
        read(6, *)
        read(6, *)
    
        allocate(xx(ni,nj),yy(ni,nj),zz(ni,nj),cp(ni,nj),sd(3,ni,nj),grad(ni,nj))
	    allocate(xc(ni,nj),yc(ni,nj),zc(ni,nj),cpc(ni,nj),cpc_unused(ni,nj))
    
        read(6,*) xx(1:ni, 1:nj  )                  !spanwise direction is z
        read(6,*) yy(1:ni, 1:nj  )
        read(6,*) zz(1:ni, 1:nj  )
        do mvari = 4, 10
            read(6,*) cpc(1:ni-1, 1:nj-1)
        enddo
        do mvari = 11, 14
            read(6,*) cpc_unused(1:ni-1, 1:nj-1)
        enddo
        
        do j = 2, nj-1
            do i = 2, ni-1
                cp(i,j) = 0.25*(cpc(i-1,j-1)+cpc(i-1,j)+cpc(i,j-1)+cpc(i,j))
            enddo
        enddo
        
        
        cp(1:ni-1, 1:nj-1) = cpc(1:ni-1, 1:nj-1)    !node cp is approxiamted by cp at cell centroid
        cp(ni,     1:nj-1) = cpc(ni-1,   1:nj-1)
        cp(1:ni-1, nj    ) = cpc(1:ni-1, nj-1  )
        cp(ni,     nj    ) = cpc(ni-1,   nj-1  )
        
        !j = 1
        !do i = 2, ni-1
        !    cp(i,j) = cpc(i-1,j) + cpc(i,j) - cp(i,j+1)
        !enddo
        !j = nj
        !do i = 2, ni-1
        !    cp(i,j) = cpc(i-1,j-1) + cpc(i,j-1) - cp(i,j-1)
        !enddo
        !i = 1
        !do j = 2, nj-1
        !    cp(i,j) = cpc(i,j-1) + cpc(i,j) - cp(i+1,j)
        !enddo
        !i = ni
        !do j = 2, nj-1
        !    cp(i,j) = cpc(i-1,j-1) + cpc(i-1,j) - cp(i-1,j)
        !enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	计算面积
	call surface(xx,yy,zz,sd,grad,cp,xc,yc,zc,ni,nj,nk,1.0)    !1.0 will change
	
!	求交点
	allocate(x(2),y(2),z(2))

	ixyz=-100
	xyz=0.0
	do ilocal=1,ilocation               

		itemp=0
		do j=1,nj
			do i=1,ni-1                   !loop of edges (i,j)-(i+1,j)
				z(1)=zz(i  ,j)            !spanwise direction: z
				z(2)=zz(i+1,j)
				zin =zlocation(ilocal)
				if( (z(2)>=zin .and. z(1)<zin) .or. &
					(z(1)>=zin .and. z(2)<zin)) then
					itemp=itemp+1         !itemp'th point on the curve 
					x(1)=xx(i  ,j)        !x(1),y(1),z(1)
					x(2)=xx(i+1,j)        !x(2),y(2),z(2)
					y(1)=yy(i  ,j)
					y(2)=yy(i+1,j)
					c1=cp(i,j)          
					c2=cp(i+1,j)
					c3=(c2-c1)/(z(2)-z(1))*(zin-z(1))+c1
					xyz(ilocal,itemp,4)=c3              !Cp for the itemp'th point on the curve 
					call cz_point_z(x,y,z,zin,xout,yout)
					xyz(ilocal,itemp,1)=xout
					xyz(ilocal,itemp,2)=yout
					xyz(ilocal,itemp,3)=zin
					ixyz(ilocal,itemp,1)=i              !what's this?  Do I need to revise it if spanwise direction is z? NO
					ixyz(ilocal,itemp,2)=j-1
					ixyz(ilocal,itemp,3)=i
					ixyz(ilocal,itemp,4)=j
				endif
			enddo
        enddo


		do i=1,ni
			do j=1,nj-1                   !loop of edges (i,j)-(i,j+1)
				z(1)=zz(i,j  )
				z(2)=zz(i,j+1)
				zin=zlocation(ilocal)
				if( (z(2)>zin .and. z(1)<zin) .or. &
					(z(1)>=zin .and. z(2)<zin)) then
					itemp=itemp+1
					x(1)=xx(i,j  )
					x(2)=xx(i,j+1)
					y(1)=yy(i,j  )
					y(2)=yy(i,j+1)
					c1=cp(i,j)
					c2=cp(i,j+1)
					c3=(c2-c1)/(z(2)-z(1))*(zin-z(1))+c1
					xyz(ilocal,itemp,4)=c3
					call cz_point_z(x,y,z,zin,xout,yout)
					xyz(ilocal,itemp,1)=xout
					xyz(ilocal,itemp,2)=yout
					xyz(ilocal,itemp,3)=zin
					ixyz(ilocal,itemp,1)=i-1    !what's this?
					ixyz(ilocal,itemp,2)=j
					ixyz(ilocal,itemp,3)=i
					ixyz(ilocal,itemp,4)=j
				endif
			enddo
        enddo
        
        
		if(itemp>0)then
			isegall=0
			do i=1,itemp                        !all points on one cross section
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
					if(lc)then                      !if point i and and point j are connected
						isegall=isegall+1
						iseg(isegall,1)=i           !start and end of a segment
						iseg(isegall,2)=j
						x1=xyz(ilocal,i,1)
						y1=xyz(ilocal,i,2)
						z1=xyz(ilocal,i,3)
						p1=xyz(ilocal,i,4)
						x2=xyz(ilocal,j,1)
						y2=xyz(ilocal,j,2)
						z2=xyz(ilocal,j,3)
						p2=xyz(ilocal,j,4)

                        
                        if(x1>chord(ilocal,1))chord(ilocal,1)=x1        !so that chord length can be determined chord(ilocal,3)
                        if(x2>chord(ilocal,1))chord(ilocal,1)=x2
                        if(x1<chord(ilocal,2))chord(ilocal,2)=x1
                        if(x2<chord(ilocal,2))chord(ilocal,2)=x2
                        
						xlength=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)


						R1X=x1-x2	!XX(I1,J,K)-XX(I,J1,K)
						R1Y=y1-y2		!YY(I1,J,K)-YY(I,J1,K)
						R1Z=1.0   !ZZ(I1,J,K)-ZZ(I,J1,K)
						R2X=x2-x1   !XX(I1,J1,K)-XX(I,J,K)
						R2Y=y2-y1     !YY(I1,J1,K)-YY(I,J,K)
						R2Z=1.0   !ZZ(I1,J1,K)-ZZ(I,J,K)


						xd=(R1Y*R2Z-R1Z*R2Y)*0.5    !y1-y2     !segment vector and its direction is the same as sd next
						yd=(R1Z*R2X-R1X*R2Z)*0.5    !x2-x1
						zd=(R1X*R2Y-R1Y*R2X)*0.5    !0.0 
                        

						xd=sign(xd,sd(1,ic,jc))
						yd=sign(yd,sd(2,ic,jc))
						zd=sign(zd,sd(3,ic,jc))
                        
                        
						clall2(ilocal,1)=clall2(ilocal,1)+(p1+p2)/2.0*xd         !force element,x
						clall2(ilocal,2)=clall2(ilocal,2)+(p1+p2)/2.0*yd         !force element,y
						clall2(ilocal,3)=clall2(ilocal,3)+(p1+p2)/2.0*zd         !force element,z
                        
                        clall2(ilocal,4)=clall2(ilocal,4)+(p1+p2)/2.0*yd*(0.25-0.5*(x1+x2))   !force moment element
                        

					endif
				enddo
            enddo

				write(8,'(a12,i7,a10,i7,a40)')'zone Nodes=',itemp,' Elements=',isegall, &
						' ZONETYPE=FELineSeg DATAPACKING=POINT'
					do i=1,itemp
						write(8,'(4e20.10)') xyz(ilocal,i,1:4) 
                    end do
                    
				do i=1,isegall
					write(8,'(2i8)')iseg(i,1),iseg(i,2)
				end do
        endif
   
        
	enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	deallocate(x,y,z)
	deallocate(xx,yy,zz,cp,sd,grad,xc,yc,zc,cpc,cpc_unused)
    
    
    enddo
	
enddo

     
    do ilocal=1,ilocation
        !chord(ilocal,3)=chord(ilocal,1)-chord(ilocal,2)             !length of chord
        chord(ilocal,3)=1.0         !we have an ice accretion on the leading edge
    enddo
    
    
    
    open(10,file='SpanLoad.dat')
    write(10,*)'variables = Location Load CX CY CL CD CM Chord'
	do i=1,ilocation
        cx=clall2(i,1)/chord(i,3)
        cy=clall2(i,2)/chord(i,3)
        cl=cy*cos(aoa)-cx*sin(aoa)
        cd=cy*sin(aoa)+cx*cos(aoa)
        cm=clall2(i,4)/chord(i,3)
	    write(10,'(f15.5,10f13.8)') abs(zlocation(i)),clall2(i,3:3),cx,cy,cl,cd,cm,chord(i,3)
    enddo
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

SUBROUTINE surface(xx,yy,zz,sd,grad,cp,xc,yc,zc,ni,nj,nk,dl)   !网格单元几何信息
	real xx(ni,nj),yy(ni,nj),zz(ni,nj)
	real sd(3,ni,nj),grad(ni,nj),cp(ni,nj)
	real xc(ni,nj),yc(ni,nj),zc(ni,nj),cpc(ni,nj),span(ni,nj)

	ni1=ni-1
	nj1=nj-1
	nk1=nk-1

    DO I=1,NI1
        I1=I+1
        DO J=1,NJ1
            J1=J+1
		    xc(I,J)=(xx(I1,J1)+xx(I1,J)+xx(I,J)+xx(I,J1))*0.25
		    yc(I,J)=(yy(I1,J1)+yy(I1,J)+yy(I,J)+yy(I,J1))*0.25   
		    zc(I,J)=(zz(I1,J1)+zz(I1,J)+zz(I,J)+zz(I,J1))*0.25
	    ENDDO
    ENDDO


    DO J=1,NJ1
        J1=J+1
        DO I=1,NI1
            I1=I+1
              
                R1X=XX(I1,J)-XX(I,J1)
                R1Y=YY(I1,J)-YY(I,J1)
                R1Z=ZZ(I1,J)-ZZ(I,J1)
                R2X=XX(I1,J1)-XX(I,J)
                R2Y=YY(I1,J1)-YY(I,J)
                R2Z=ZZ(I1,J1)-ZZ(I,J)
                
                SD(1,I,J)=(R1Y*R2Z-R1Z*R2Y)*0.5*dl      !R1*R2   !d1 will be determined
                SD(2,I,J)=(R1Z*R2X-R1X*R2Z)*0.5*dl
                SD(3,I,J)=(R1X*R2Y-R1Y*R2X)*0.5*dl
                
	    enddo
    enddo
    

    do i=1,ni1
	    do j=1,nj1
	    do k=1,1
	            grad(i,j)=SQRT(SD(1,I,J)**2+SD(2,I,J)**2+SD(3,I,J)**2)  !面积
				  
				sd(1,i,j)=sd(1,i,j)/grad(i,j)           !方向
				sd(2,i,j)=sd(2,i,j)/grad(i,j)
				sd(3,i,j)=sd(3,i,j)/grad(i,j)

	        enddo
	    enddo
    enddo

end subroutine


    subroutine getfilename(iblkid,npre,prename,npost,postname,filename)
    
    implicit none
    integer:: iblkid,npre,npost
    character(len=npre):: prename
    character(len=npost)::postname
    character(len=4)::buf
    character(len=npre+npost+4)::filename
    
	if(iblkid.lt.10)then
		write(buf,'(i1)') iblkid
		filename = trim(prename)//trim('000')//trim(buf)//trim(postname)
	elseif(iblkid.lt.100)then
		write(buf,'(i2)') iblkid
		filename = trim(prename)//trim('00')//trim(buf)//trim(postname)
	elseif(iblkid.lt.1000)then
		write(buf,'(i3)') iblkid
		filename = trim(prename)//trim('0')//trim(buf)//trim(postname)
	elseif(iblkid.lt.10000)then
		write(buf,'(i4)') iblkid
		filename = trim(prename)//trim(buf)//trim(postname)
    else
        write(*,*)'Error: maxblock is more than 10000 @ 2901'
        stop
	endif
	
	end subroutine   
