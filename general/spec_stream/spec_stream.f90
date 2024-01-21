!  Spectra_Stream.f90 
!
!  FUNCTIONS:
!  Spectra_Stream - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Spectra_Stream
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program Spectra_Stream
        
    implicit none
    
    integer             ::  nOrder,nPos,nFreq,i,j
    real,allocatable    ::  freq(:),psd(:,:),pos(:)
    character(len=1000) ::  fname(100),ctmp
    
    
    open( 11, file = 'Input_Spectra_Stream.txt', action = 'read' )
    read( 11, * ) ctmp,nOrder
    read( 11, * ) ctmp,nPos
    allocate(pos (nPos))
    read( 11, * ) (fname(i), pos(i), i=1,nPos)
    close(11)
    
    nFreq = (2**nOrder/2)+1
    allocate(freq(nFreq))
    allocate(psd (nFreq,nPos))
    
    do i = 1, nPos
        open( 12, file = fname(i), action = 'read' ) 
        do j = 1, nFreq
            read( 12, * )  freq(j), psd(j,i)
        enddo
        close(12)
    enddo
    
    open(13, file = 'Output_Spectra_Stream.dat', action = 'write')
    write(13, *) 'VARIABLES="x","freq","psd"'
    write(13, *) 'zone i=   ', nPos, 'j=', nFreq-1
    write(13, *) 'DATAPACKING=BLOCK'
    write(13, *) (pos(:), i = 2,nFreq)
    write(13, *) ((freq(i), j = 1,nPos), i = 2,nFreq)
    write(13, *) ((psd(j,i),i = 1,nPos), j = 2,nFreq)
    close(13)
    
    deallocate(freq)
    deallocate(psd)
    deallocate(pos)

    end program Spectra_Stream

