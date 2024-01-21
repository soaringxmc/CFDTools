program Icem2CFL3D
    implicit none
    character*1000      ::  filegrd,fileinp,outfile
    character*1000      ::  temp,temp2
    integer             ::  ipre,iline,ir(2),nblk,i,iprd
    integer             ::  nbli
    integer,allocatable ::  iblock(:,:)
    character*10        ::  chari

!------------- read input file ---------------------------
	fileinp = 'cfl3d.inp'
    iprd = 0
   !iprd = 1 if periodic bc's desired
!------------- calculate nblk from fileinp ----------------------------
    open(6,file=fileinp,action='read',status='old')
    iline = 0
    do while(.not. eof(6))
        iline = iline + 1
        read(6,'(a1000)') temp
        if(temp(1:10) == '       NCG') then
            ir(1) = iline
        else if(temp(1:10) == '      IDIM') then
            ir(2) = iline
            exit
        end if
    end do
    nblk = ir(2) - ir(1) - 1
    allocate(iblock(nblk,3))
    do i = 1,nblk
        read(6,*) iblock(i,:)
    end do
    
!------------ header file -----------------------------------
    open(51,file='cfl3d.inp2')
	write(51,'( A6)')'FILES:'
	write(51,'(a)'  )'cfl3d.x'
	write(51,'(a)'  )'field.g'
	write(51,'(a)'  )'field.q'
	write(51,'( A9)')'cfl3d.out'
	write(51,'( A9)')'resid.out'
	write(51,'(A12)')'cfl3d.turres'
	write(51,'(A12)')'cfl3d.blomax'
	write(51,'(A10)')'cfl3d.2out'
	write(51,'( A9)')'cfl3d.prt'
	write(51,'(A11)')'cfl3d.press'
	write(51,'( A9)')'ovrlp.bin'
	write(51,'( A9)')'patch.bin'
	write(51,'(A13)')'cfl3d.restart'
    write(51,'(a)'  ) '>'
    write(51,'(a)'  ) 'cprec 0'
    write(51,'(a)'  ) 'ides 0'
    write(51,'(a)'  ) 'iblend 0'
    write(51,'(a)'  ) 'iteravg 0'
    write(51,'(a)'  ) 'icoarsemovie 0'
    write(51,'(a)'  ) 'i2dmovie 0'
    write(51,'(a)'  ) 'isample 0'
    write(51,'(a)'  ) '<'
	write(51,'(A42)')'CFL3D V6 INPUT FILE GENERATED WITH GRIDGEN'
	write(51,'(A70)')'    XMACH     ALPHA      BETA  REUE,MIL   TINF,DR     IALPH    IHSTRY'
	write(51,'(A70)')'    0.750     1.000    0.0000      3.00    460.00         1         0'
	write(51,'(A60)')'     SREF      CREF      BREF       XMC       YMC       ZMC'
	write(51,'(A60)')'     1.00      1.00      1.00      0.25      0.00      0.00'
	write(51,'(A60)')'       DT     IREST   IFLAGTS      FMAX     IUNST    CFLTAU'
	write(51,'(A60)')'    -1.50         0       000      1.00         0    7.5000'
	write(51,'(A80)')'    NGRID   NPLOT3D    NPRINT    NWREST      ICHK       I2D    NTSTEP       ITA'
	write(51,'(8i10)')    -nblk,    nblk,       -1,      500,        0,        0,      1,          -2
	write(51,'(A70)')'       NCG       IEM  IADVANCE    IFORCE  IVISC(I)  IVISC(J)  IVISC(K)'
	do i=1,nblk
		write(51,'(7i10)')        2,        0,        0,      333,        7,        7,        7
    end do
    write(51,'(A30)')'      IDIM      JDIM      KDIM'
	do i=1,nblk
		write(51,'(3i10)') iblock(i,1),iblock(i,2),iblock(i,3)
    end do
    	write(51,'(A60)')'    ILAMLO    ILAMHI    JLAMLO    JLAMHI    KLAMLO    ILAMHI'
	do i=1,nblk
		write(51,'(6i10)')   0,        0,        0,        0,        0,        0
	end do
    	write(51,'(A80)')'     INEWG    IGRIDC        IS        JS        KS        IE        JE        KE'
	do i=1,nblk
		write(51,'(8i10)')     0,        0,        0,        0,        0,        0,        0,        0
    end do	
    	write(51,'(A60)')'  IDIAG(I)  IDIAG(J)  IDIAG(K)  IFLIM(I)  IFLIM(J)  IFLIM(K)'
	do i=1,nblk
		write(51,'(6i10)')        1,        1,        1,        4,        4,        4
    end do
    write(51,'(A60)')'   IFDS(I)   IFDS(J)   IFDS(K)  RKAP0(I)  RKAP0(J)  RKAP0(K)'
	do i=1,nblk
		write(51,'(3i10,3f10.4)')         1,        1,        1,   0.3333,   0.3333,   0.3333
    end do
    
    do while(.not. eof(6))
        read(6,'(a1000)') temp
        if(temp(1:10) == '      GRID') exit
    end do
    
    write(51,'(A80)')'      GRID     NBCIO   NBCIDIM     NBCJO   NBCJDIM     NBCKO   NBCKDIM    IOVRLP'
    do i = 1,nblk
        read(6,'(a1000)') temp
        write(51,'(a)') trim(temp)
    end do
!------------------ boundary condition -----------------------------
    !------------------ i0 ---------------------------
    do while(.not. eof(6))
        read(6,'(a1000)') temp
        if(temp(1:10) == '      MSEQ') then
            backspace(6)
            exit
        end if
        if(temp(27:30) == '1000') then
            temp(27:30) = '1003'
            write(51,'(a)') trim(temp)
        else if(temp(27:30) == '1001' .and. iprd == 1) then
            temp(27:30) = '2005'
            temp(71:80) = '         4'
            write(51,'(a)') trim(temp)
            write(51,'(a)') '              hgridp      dthx      dthy      dthz'
            write(51,'(a)') '          '//temp(1:10)//'       0.0       0.0       0.0'
        else if(temp(27:30) == '2004') then
            write(51,'(a)') trim(temp)
            read(6,'(a1000)') temp
            write(51,'(a)') trim(temp)
            read(6,'(a1000)') temp
            write(51,'(a)') '                 0.0       0.0'
        else if(temp(27:30) == '2003') then
            write(51,'(a)') trim(temp)
            read(6,'(a1000)') temp
            write(51,'(a)') trim(temp)
            read(6,'(a1000)') temp
            write(51,'(a)') '                 XXX       XXX       XXX       0.0       0.0'
        else if(temp(27:30) == '2002') then
            write(51,'(a)') trim(temp)
            read(6,'(a1000)') temp
            write(51,'(a)') trim(temp)
            read(6,'(a1000)') temp
            write(51,'(a)') '                 1.0'
        else if(temp(27:30) == '2005') then
            write(51,'(a)') trim(temp)
            temp2 = temp
            read(6,'(a1000)') temp
            write(51,'(a)') trim(temp)
            read(6,'(a1000)') temp
            write(51,'(a)') '          '//temp2(1:10)//'       0.0       0.0       0.0'
        else
            write(51,'(a)') trim(temp)
        end if
    end do
!--------------------- multigrid --------------------------
    write(51,'(a50)')'      MSEQ    MGFLAG    ICONSF       MTT      NGAM'
    write(51,'(a50)')'         3         1         0         0         2'
    write(51,'(a80)')'      ISSC EPSSSC(1) EPSSSC(2) EPSSSC(3)      ISSR EPSSSR(1) EPSSSR(2) EPSSSR(3)'
    write(51,'(a80)')'         0       0.3       0.3       0.3         0       0.3       0.3       0.3'
    write(51,'(a40)')'      NCYC    MGLEVG     NEMGL     NITFO'
    write(51,'(a40)')'       500         1         0         0'
    write(51,'(a40)')'       500         2         0         0'
    write(51,'(a40)')'      1000         3         0         0'
    write(51,'(a80)')'      MIT1      MIT2      MIT3      MIT4      MIT5      MIT6      MIT7      MIT8'
    write(51,'(a80)')'         1         1         1         1         1         1         1         1'
    write(51,'(a80)')'         1         1         1         1         1         1         1         1'
    write(51,'(a80)')'         1         1         1         1         1         1         1         1'
!------------------- point to point patch -------------------
    do while(.not. eof(6))
        read(6,'(a1000)') temp
        if(temp(1:20) == '   1-1 BLOCKING DATA') then
            write(51,'(a)') trim(temp)
            read(6,*)
            read(6,*) nbli
            exit
        end if
    end do
    write(51,'(a10)') '      NBLI'
    write(51,'(i10)') nbli
    
    do i = 1,(nbli+1)*2
        read(6,'(a1000)') temp
        write(51,'(a)') trim(temp)
    end do
!--------------------------------------
      read(6,*)
      write(51,'(a)')'  PATCH SURFACE DATA:'
      read(6,*) 
      write(51,'(a)')'    NINTER'
      read(6,*) i
      write(51,'(i10)') i

      write(51,'(a)')'PLOT3D OUTPUT:'
      write(51,'(a)')'  GRID  IPTYPE  ISTART  IEND  IINC  JSTART  JEND  JINC  KSTART  KEND  KINC'  
      do i = 1,nblk
        write(51,'(i6,a)') i,'      -2       0     0     0       0     0     0       0     0     0'
      end do
      write(51,'(a)')' UPDATE'
      write(51,'(a)')'   500'
      write(51,'(a)')'  PRINT OUT:'
      write(51,'(a)')'  GRID  IPTYPE  ISTART  IEND  IINC  JSTART  JEND  JINC  KSTART  KEND  KINC'
      write(51,'(a)')'     1       0       0     0     0       0     0     0       0     0     0'
      write(51,'(a)') ' SAMPLING POINTS:'
      write(51,'(a)') ' NSMP'
      write(51,'(a)') '    0'
      write(51,'(a)') '  GRID  IPTYPE  ISTART  IEND  IINC  JSTART  JEND  JINC  KSTART  KEND  KINC'
      write(51,'(a)') ' CONTROL SURFACE:'
      write(51,'(a)') ' NCS'
      write(51,'(a)') '   0'
      write(51,'(a)') '  GRID  ISTART    IEND  JSTART    JEND  KSTART    KEND   IWALL   INROM'
      close(6,status='delete')
      close(51)
      call system('rename cfl3d.inp2 cfl3d.inp')
      
    write(*,*)'Input file generated!'

end program Icem2CFL3D