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