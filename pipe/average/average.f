       program average 
!
       integer :: m1,m2,m3
       integer :: m1m,m2m,m3m
       integer :: m1mh,m2mh,m3mh
       integer, parameter :: dp = 8
!
       integer, parameter :: nplanes=4
       integer, parameter :: nvars=5
       integer, parameter :: npdf=200
       integer, parameter :: nvwav=30
       integer, parameter :: nbudget=6
!
       real(dp), allocatable :: spec2d(:,:),spec2dl(:,:)
       real(dp), allocatable :: spect(:,:,:),spectl(:,:,:)
       real(dp), allocatable :: spect0(:,:,:),spect0l(:,:,:)
       real(dp), allocatable :: specz(:,:,:),speczl(:,:,:)
       real(dp), allocatable :: wav(:,:),wavl(:,:)
       real(dp), allocatable :: flu(:,:,:),flul(:,:,:)
       real(dp), allocatable :: rc(:),rm(:)
       real(dp), allocatable, dimension(:,:,:) :: pdfl,pdf
       real(dp), allocatable, dimension(:,:,:,:,:) :: jpdfl,jpdf
       real(dp), allocatable :: wvar(:,:),wvarav(:,:)
       real(dp), allocatable :: cfl(:,:),cflav(:,:)
!
       character(4) :: naitav
       character(1) :: navar,naplan,nal,nam
!
       open(unit=15,file='../bou.in',form='formatted')
       read(15,*)
       read(15,*) m1,m2,m3 
       close(15)
!
       m1m = m1-1; m2m = m2-1; m3m = m3-1
       m1mh = m1m/2 + 1
       m2mh = m2m/2 + 1
       m3mh = m3m/2 + 1
!
       allocate(rc(1:m2))
       allocate(rm(1:m2m))
       allocate(spec2d(m1mh,m3mh),spec2dl(m1mh,m3mh))
       allocate(spect(m1mh,m2m,nvars),spectl(m1mh,m2m,nvars))
       allocate(spect0(m1mh,m2m,nvars),spect0l(m1mh,m2m,nvars))
       allocate(specz(m3mh,m2m,nvars),speczl(m3mh,m2m,nvars))
       allocate(wav(m2,nbudget),wavl(m2,nbudget))
       allocate(flu(m2,4,3),flul(m2,4,3))
       allocate (pdfl(npdf,m2m,nvars),pdf(npdf,m2m,nvars))
       allocate(jpdfl(npdf,npdf,nplanes,nvars,nvars),
     .          jpdf(npdf,npdf,nplanes,nvars,nvars))
       allocate(wvar(m2m,30),wvarav(m2m,30))
       allocate(cfl(m2m,3),cflav(m2m,3))
!
       read(*,*) ibeg,iend
!
!      Read wmean
!
       wvarav=0.
       do l=ibeg,iend
        write(naitav,'(I4.4)') l
        open(10,file='wmean_'//naitav//'.dat',form='formatted') 
        do j=1,m2m
         read(10,*) rc(j),rm(j),(wvar(j,n),n=1,30)
        enddo
        close(10)
        wvarav = wvarav+wvar
       enddo
       wvarav = wvarav/(iend-ibeg+1)
       open(10,file='wmean_ave.dat')
       do j=1,m2m
        write(10,'(32E20.10)') rc(j),rm(j),(wvarav(j,n),n=1,30)
       enddo
       close(10)
!
!      Read cfl
!
       cflav=0.
       do l=ibeg,iend
        write(naitav,'(I4.4)') l
        open(10,file='cfl_'//naitav//'.dat',form='formatted')
        do j=1,m2m
         read(10,*) rc(j),rm(j),(cfl(j,n),n=1,3)
        enddo
        close(10)
        cflav = cflav+cfl
       enddo
       cflav = cflav/(iend-ibeg+1)
       open(10,file='cfl_ave.dat')
       do j=1,m2m
        write(10,'(5E20.10)') rc(j),rm(j),(cflav(j,n),n=1,3)
       enddo
       close(10)
!
!      Read PDF
!
       pdf=0.
       do itav=ibeg,iend
        print *,'PDF sample',itav
        write(naitav,1004) itav 
        open(10,file='pdf_'//naitav//'.bin',form='unformatted')
        read(10) pdfl
        close(10)
        pdf=pdf+pdfl
       enddo
       pdf = pdf/(iend-ibeg+1)
       open(10,file='pdf_ave.bin',form='unformatted')
       write(10) pdf
       close(10)
!      
       jpdf=0.
       do itav=ibeg,iend
        print *,'JPDF sample',itav
        write(naitav,1004) itav 
        open(10,file='jpdf_'//naitav//'.bin',form='unformatted')
        read(10) jpdfl
        close(10)
        jpdf=jpdf+jpdfl
       enddo
       jpdf = jpdf/(iend-ibeg+1)
       open(10,file='jpdf_ave.bin',form='unformatted')
       write(10) jpdf
       close(10)
! 
!      Read budgets
!
       do l=1,3
        write(nal,1001) l
        do m=l,3
         write(nam,1001) m
         wav=0.
         do itav=ibeg,iend
          print *,'Budget',l,m,itav
          write(naitav,1004) itav 
          open(10,file='budget_'//nal//nam//'_'//naitav//'.dat',
     .          form='formatted') 
          do j=1,m2m
           read(10,100) rc(j),rm(j),wavl(j,1:nbudget)
          enddo
          close(10) 
          wav=wav+wavl 
         enddo
         wav = wav/(iend-ibeg+1)
         open(10,file='budget_ave_'//nal//nam//'.dat',
     .         form='formatted') 
         do j=1,m2m
          write(10,100) rc(j),rm(j),wav(j,1:nbudget)
         enddo
         close(10) 
        enddo
       enddo
!
       l=4
       m=4
       write(nal,1001) l
       write(nam,1001) m
       wav=0.
       do itav=ibeg,iend
        write(naitav,1004) itav
        open(10,file='budget_'//nal//nam//'_'//naitav//'.dat',
     .        form='formatted')
        do j=1,m2m
         read(10,100) rc(j),rm(j),wavl(j,1:nbudget)
        enddo
        close(10)
        wav=wav+wavl
       enddo
       wav = wav/(iend-ibeg+1)
       open(10,file='budget_ave_'//nal//nam//'.dat',
     .       form='formatted')
       do j=1,m2m
        write(10,100) rc(j),rm(j),wav(j,1:nbudget)
       enddo
       close(10)
!
       do l=1,4
        write(nal,1001) l
        flu=0.
        do itav=ibeg,iend
         print *,'Flux',l,itav
         write(naitav,1004) itav 
         open(10,file='flux_'//nal//'_'//naitav//'.dat',
     .          form='formatted') 
         do j=1,m2m
          read(10,100) rc(j),rm(j),flul(j,l,1:3)
         enddo
         close(10) 
         flu=flu+flul 
        enddo
        flu = flu/(iend-ibeg+1)
        open(10,file='flux_ave_'//nal//'.dat',
     .        form='formatted') 
        do j=1,m2m
         write(10,100) rc(j),rm(j),flu(j,l,1:3)
        enddo
        close(10) 
       enddo
!
!      Read 1d spectra
!
       spect=0.
       do itav=ibeg,iend
        print *,'Theta spectra, sample',itav
        write(naitav,1004) itav 
        open(10,file='specth_'//naitav//'.bin',form='unformatted')
        read(10) spectl
        close(10)
        spect=spect+spectl
       enddo
       spect = spect/(iend-ibeg+1)
       open(10,file='specth_ave.bin',form='unformatted')
       write(10) spect
       close(10)
! 
       specz=0.
       do itav=ibeg,iend
        print *,'Zeta spectra, sample',itav
        write(naitav,1004) itav 
        open(10,file='specze_'//naitav//'.bin',form='unformatted')
        read(10) speczl
        close(10)
        specz=specz+speczl
       enddo
       specz = specz/(iend-ibeg+1)
       open(10,file='specze_ave.bin',form='unformatted')
       write(10) specz
       close(10)
! 
       spect0=0.
       do itav=ibeg,iend
        print *,'Theta0 spectra, sample',itav
        write(naitav,1004) itav 
        open(10,file='specth0_'//naitav//'.bin',form='unformatted')
        read(10) spect0l
        close(10)
        spect0=spect0+spect0l
       enddo
       spect0 = spect0/(iend-ibeg+1)
       open(10,file='specth0_ave.bin',form='unformatted')
       write(10) spect0
       close(10)
! 
!      Read 2d spectra
!
       do m=1,nvars 
        write(navar,1001) m
!
        do l=1,nplanes
         write(naplan,1001) l
!
         spec2d = 0.d0
!
         do itav=ibeg,iend
          write(naitav,1004) itav 
!
          print *,'2D spectra: plane, variable, sample',
     .     l,m,itav
!
          open(10,file='spec2d'//naplan//'_'//navar//'_'//naitav//
     .         '.bin' ,form='unformatted')
          read(10) spec2dl
          print *,minval(spec2dl(:,:)),maxval(spec2dl(:,:))
          close(10)
          spec2d=spec2d+spec2dl 
         enddo
!
         spec2d = spec2d/(iend-ibeg+1)
!
         open(10,file='spec2d_ave_'//naplan//'_'//navar//'.bin',
     .            form='unformatted')
         print *,minval(spec2d(:,:)),maxval(spec2d(:,:))
         write(10) spec2d
         close(10)
!
        enddo
!
       enddo
!
  100 format(40E20.10)
 1001 format(I1.1)
 1004 format(I4.4)    
!
       stop
       end
