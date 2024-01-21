       program pdfave
!
       integer :: m1,m2,m3
       integer :: m1m,m2m,m3m
       integer :: m1mh,m2mh,m3mh
       integer, parameter :: dp = 8
!
       integer, parameter :: nvar=5
       integer, parameter :: npdf=200
       integer, parameter :: nplane=4
!
       real(dp), allocatable :: rc(:),rm(:)
       real(dp), allocatable, dimension(:,:,:) :: pdf
       real(dp), allocatable, dimension(:,:,:,:,:) :: jpdf
       real(dp) :: pdfmin(nvar),pdfmax(nvar)
       real(dp) :: uu(npdf,nvar)
       real(dp) :: dupdf(nvar)
       integer :: jplane(nplane)
!
       character(1) :: navar,naplan,navar1,navar2
       character(4) :: naitav
!
       open(unit=15,file='bou.in',form='formatted')
       read(15,*)
       read(15,*) m1,m2,m3 
       close(15)
!
       m1m = m1-1; m2m = m2-1; m3m = m3-1
!
       allocate(rc(1:m2))
       allocate(rm(1:m2m))
       allocate(pdf(npdf,m2m,nvar))
       allocate(jpdf(npdf,npdf,nplane,nvar,nvar))
!
       open(unit=98,file='radcor.out',status='unknown')
       do j=1,m2m
        read(98,*) jj,rc(j),rm(j)
       enddo
       close(98)
!
       jplane(1) = m2m    ! first off-wall cell
       jplane(2) = m2m-16 ! y^+=15
       do j=1,m2m
        if (rm(j).gt.0.8) exit ! y=0.2
       enddo
       jplane(3) = j
       do j=1,m2m
        if (rm(j).gt.0.5) exit ! y=0.5
       enddo
       jplane(4) = j
!
       pdfmin(1)=-0.25
       pdfmin(2)=-0.25
       pdfmin(3)=-0.1
       pdfmin(4)=-0.1
       pdfmin(5)=-0.03
       pdfmax(1)= 0.25
       pdfmax(2)= 0.25
       pdfmax(3)= 0.8
       pdfmax(4)= 0.8
       pdfmax(5)= 0.03
       dupdf=(pdfmax-pdfmin)/npdf ! bin width
!
       do iv = 1,nvar
        do ip = 1,npdf
         uu(ip,iv) = pdfmin(iv) + (ip-0.5)*dupdf(iv)
        enddo
       enddo
!
!      pdf average if average is not executed in advance
!
       !print *, 'ibeg,iend'
       !read(*,*) ibeg,iend
       !pdf=0.
       !do itav=ibeg,iend
       ! print *,'pdf sample',itav
       ! write(naitav,'(I4.4)') itav
       ! open(10,file='pdf_'//naitav//'.bin',form='unformatted')
       ! read(10) pdfl
       ! close(10)
       ! pdf=pdf+pdfl
       !enddo
       !pdf = pdf/(iend-ibeg+1)
!
!      average is not executed in advance
!
       open(10,file='pdf_ave.bin',form='unformatted')
       read(10) pdf
       close(10)
!
       do iv = 1,nvar
        do j  = 1,nplane
         write(navar ,'(I1)') iv
         write(naplan,'(I1)') j
         open(20,file='pdf_'//naplan//'_'//navar//'.dat',action='write')
         do ip = 1,npdf
          !probability density function
          write(20,*) uu(ip,iv),pdf(ip,jplane(j),iv)/dupdf(iv)
         enddo
         close(20)
        enddo
       enddo
!
!      Read JPDF
!
       open(10,file='jpdf_ave.bin',form='unformatted')
       read(10) jpdf
       close(10)
!
       do iiv= 1,nvar
       do iv = 1,iiv
        do j = 1,nplane
         write(navar1,'(I1)') iv  !x coordinate
         write(navar2,'(I1)') iiv !y coordinate
         write(naplan,'(I1)') j
         open(20,file='jpdf_'//naplan//'_'//navar1//'_'//navar2//'.dat')
         write(20,*) 'zone i=',npdf,'j=',npdf
         do iip = 1,npdf
         do ip = 1,npdf
          !probability density function
          write(20,*) uu(ip,iv),uu(iip,iiv),jpdf(ip,iip,j,iv,iiv)/(dupdf(iv)*dupdf(iiv))
         enddo
         enddo
         close(20)
        enddo
       enddo
       enddo
!
       end program pdfave
