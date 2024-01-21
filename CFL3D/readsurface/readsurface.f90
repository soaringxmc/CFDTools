program readsurface
	parameter(mi=10000000)
	character*80 filename,filename1,filename2
	character*200 char
	integer ijk(mi,3)
	real    xyz(mi,12)


	filename='cfl3d.prt'
	filename1='surface.dat'
	filename2='surface2.dat'
	open(5,file=filename)
	open(6,file=filename1)
	open(7,file=filename2)
	write(6,*)'variables = X Y Z I J K U V W P T M Cp ut'
	write(7,*)'variables = X Y Z I J K dn P T Cf Ch yplus'

	do !isurface=1,1000000000
			isurface=isurface+1
			read(5,'(a63)',end=1012)char
			if(char=='   I   J   K       X            Y            Z          U/Uinf')then
            write(*,*)'Found Surf!'
			imax=0;imin=1000000;jmax=0;jmin=1000000;kmax=0;kmin=1000000
			i0=0
			do i=1,mi
				read(5,*,err=1011,end=1011)(ijk(i,j),j=1,3),(xyz(i,j),j=1,11)
				i0=i
				imax=max(imax,ijk(i,1))
				imin=min(imin,ijk(i,1))
				jmax=max(jmax,ijk(i,2))
				jmin=min(jmin,ijk(i,2))
				kmax=max(kmax,ijk(i,3))
				kmin=min(kmin,ijk(i,3))
			end do
	1011	continue
			temp=(imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)
			if(i0.eq.temp)then
				if(imax==imin)then
					write(6,*)'zone i= ',imax-imin+1,' j= ',jmax-jmin+1,' k= ',kmax-kmin+1
				elseif(jmax==jmin)then
					write(6,*)'zone i= ',kmax-kmin+1,' j= ',jmax-jmin+1,' k= ',imax-imin+1
				else
					write(6,*)'zone i= ',jmax-jmin+1,' j= ',imax-imin+1,' k= ',kmax-kmin+1
				endif
				do i=1,i0
					write(6,'(3e18.10,3i4,8e18.10)')(xyz(i,j),j=1,3),(ijk(i,j),j=1,3),(xyz(i,j),j=4,11)
				end do
			endif
			write(*,*)'Surface=',isurface
		endif

	enddo 
1012 continue
	rewind(5)

	do !isurface=1,1000000000
			isurface=isurface+1
			read(5,'(a60)',end=1100)char
			if(char=='   I   J   K      X            Y            Z            dn')then
                write(*,*)'Found Turb!'

			imax=0;imin=1000000;jmax=0;jmin=1000000;kmax=0;kmin=1000000
			i0=0
			do i=1,mi
				read(5,*,err=1013,end=1013)(ijk(i,j),j=1,3),(xyz(i,j),j=1,9)
				i0=i0+1
				imax=max(imax,ijk(i,1))
				imin=min(imin,ijk(i,1))
				jmax=max(jmax,ijk(i,2))
				jmin=min(jmin,ijk(i,2))
				kmax=max(kmax,ijk(i,3))
				kmin=min(kmin,ijk(i,3))
			end do
	1013	continue
			temp=(imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)
			if(i0.eq.temp)then
				if(imax==imin)then
					write(7,*)'zone i= ',imax-imin+1,' j= ',jmax-jmin+1,' k= ',kmax-kmin+1
				elseif(jmax==jmin)then
					write(7,*)'zone i= ',kmax-kmin+1,' j= ',jmax-jmin+1,' k= ',imax-imin+1
				else
					write(7,*)'zone i= ',jmax-jmin+1,' j= ',imax-imin+1,' k= ',kmax-kmin+1
				endif

				do i=1,i0
					write(7,'(3e18.10,3i4,6e18.10)')(xyz(i,j),j=1,3),(ijk(i,j),j=1,3),(xyz(i,j),j=4,9)
				end do
			endif
			write(*,*)'Turb.Surface=',isurface
			endif
		enddo 
	1100	continue

end program	

