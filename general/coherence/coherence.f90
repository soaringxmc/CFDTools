program Correlation              !after Read_interpolation.f90.
  !      Compile the code with:  ifort -r8 UVWP_Corr.f90 -L/pscratch/jhao/fftw/lib/ -lfftw3 -lm 
  implicit none
        integer,parameter                       :: n3m=256,n2m=202
        integer,parameter                       :: skip=4
        integer,parameter                       :: n3mk=int((n3m-1)/skip)+1 !skip points
        integer,parameter                       :: nzband=128
        integer,parameter                       :: ntband=128
        integer,parameter                       :: nf=2*ntband
        
        real(4),allocatable                     :: pnew(:,:)
        real(8)                                 :: dtnew,dz=4./n3m*skip*8.

        real(4)                                 :: pcorr(-ntband:ntband,-nzband:nzband),pcsnd(-ntband:ntband)
        integer                                 :: k,j,it,ntsum
        integer                                 :: nt0,nt0bgn,nt0end
        integer                                 :: nt,nz,nzt,n

        real(8)                                 :: pwfft(nf+2,-nzband:nzband)
        real(8)                                 :: pwspt(0:nf/2,-nzband:nzband) 
        real(8)                                 :: w(nf)
        real(8)                                 :: in(nf)
        complex(8)                              :: out(nf/2+1)
        real(8)                                 :: Repw(nf/2+1),Impw(nf/2+1)
        integer(8)                              :: plan,isym
        integer                                 :: FFTW_ESTIMATE
        real(8)                                 :: pi,scale,pmsf_sp
        
        open(200,file='fluct_pressure.dat',form='unformatted')
        read(200)  ntsum,dtnew
        allocate(pnew(ntsum,n3mk))
        nt0bgn = 1+ntband
        nt0end = ntsum-ntband

        do it=1,ntsum
           read(200)   (pnew(it,k),k=1,n3mk)
        end do
        close(200)

        pcorr = 0.0d0
        pcsnd = 0.0d0 

        do nz=-nzband,nzband
           do nt=-ntband,ntband
              do k=1,n3mk
                 nzt = nz + k
                 if (nzt .gt.n3mk)  nzt = nzt - n3mk
                 if (nzt .lt.   1)  nzt = nzt + n3mk
                 do it = nt0bgn,nt0end
                    pcorr(nt,nz)=pcorr(nt,nz)+pnew(it,k)*pnew(it+nt,nzt)/DBLE(n3mk)/(nt0end-nt0bgn+1)   ! Rpp(dt,dz;z)  
                 end do
              end do
           end do
        end do
        
        do nt=-ntband,ntband
           do k=1,n3mk
              nzt = nz + k
              if (nzt .gt.n3mk)  nzt = nzt - n3mk
              if (nzt .lt.   1)  nzt = nzt + n3mk
              do it = nt0bgn,nt0end
                 pcsnd(nt)=pcsnd(nt)+pnew(it,nzt)*pnew(it+nt,nzt)/DBLE(n3mk)/(nt0end-nt0bgn+1)         ! Rpp(dt,0;z+dz)  
              end do
           end do
        end do

        print*,pcorr(0,0),pcsnd(0)

! ..... symmetrize the correlation data w.r.t. zero separation .....

        if (isym.eq.1) then
           do nz=-nzband,nzband
              do nt=1,ntband
                 pcorr( nt,nz) = 0.5*(pcorr(-nt,nz)+pcorr(nt,nz))
                 pcorr(-nt,nz) = pcorr(nt,nz)
            end do
         end do
      end if

! ..... calculate frequency spectra by Fourier transform of temporal
!       correlation functions .....

      call dfftw_plan_dft_r2c_1d(plan,nf,in,out,FFTW_ESTIMATE)
      print*,plan
! ..... select the window function .....

      do  n=1,nf
         w(n) = 0.5*(1.-cos(2.*pi*(n-1)/nf))
      end do
      do nz=-nzband,nzband
         do n=-nf/2,nf/2-1
            pwfft(n+1+nf/2,nz) = pcorr(n,nz)
         end do
      end do
      
! ..... multiply by the window function .....

      do nz=-nzband,nzband
         do n=1,nf
            pwfft(n,nz) = pwfft(n,nz)*w(n)
         end do
      end do

      scale = 1./nf

      do nz=-nzband,nzband
         do n=1,nf
            in(n)=pwfft(n,nz)
         end do
         call dfftw_execute_dft_r2c(plan,in,out)
         do n=1,nf/2+1
            Repw(n)=(out(n)+conjg(out(n)))/2.0d0
            Impw(n)=(out(n)-conjg(out(n)))/2.0d0
            pwspt(n-1,nz)=sqrt(Repw(n)**2.0d0+Impw(n)**2.0d0)*scale
         end do
      end do
      
      call dfftw_destroy_plan(plan)
      
! ..... check to see if the integral of the spectrum is equal to 
!       the m.s. fluctuations .....

      print*, ' mean-square and rms of p from frequency spectra:'
      do nz=-nzband,nzband
         pmsf_sp = pwspt(0,nz)
         do n=1,nf/2
            pmsf_sp = pmsf_sp + 2.0*pwspt(n,nz)
         end do
         print*, 'nz=',nz,'   pms=', pmsf_sp, '   prms=', sqrt(pmsf_sp)
      end do

      print*, ' mean-square and rms of p from correlation functions:'
      do nz=-nzband,nzband
         print*, 'nz=',nz, '    pms =  ', pcorr(0,nz), '   prms =  ', sqrt(pcorr(0,nz))
      end do

      do nz=-nzband,nzband
         do n=1,nf/2
            pwspt(n,nz) = pwspt(n,nz)*pcorr(0,nz)/pmsf_sp
         end do
      end do


!.....WRITE OUTPUT
      open(84,file='coherence.dat',form='formatted',action='write')
      WRITE(84,*)  'variables = "w","dz","Coherence"'
      WRITE(84,*)'zone t="2D" ,i=',nf/2,', j=',nzband*2+1
      DO nz=-nzband,nzband
         DO n=1,nf/2
            WRITE(84,85) n/(nf*dtnew),dz*nz,pwspt(n,nz)**2./(pwspt(n,0)**2.)
         END DO
      END DO
85    FORMAT(3(2X,E15.8))
      CLOSE(84)
      
      stop
    end program Correlation
