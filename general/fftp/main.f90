program main
    use Variables
    implicit none
    real temp,temp2
    character*100 ctemp
    integer iwindow
!    KofFFT=11
!    CFDpointNum=48
    pi=3.14159265359
    
    open(5,file='fftp.txt')
    read(5,*)ctemp,KofFFT
    read(5,*)ctemp,n_window
    read(5,*)ctemp,n_startline
    read(5,*)ctemp,CFDpointNum
    read(5,*)ctemp,window
    read(5,*)ctemp,filename
    read(5,*)ctemp,fName2
    read(5,*)ctemp,kVar
    read(5,*)ctemp,kCode
    read(5,*)ctemp,kPSD
    read(5,*)ctemp,sound
    read(5,*)ctemp,mach
    read(5,*)ctemp,density
    read(5,*)ctemp,timestep
    read(5,*)ctemp,length
!!!!!!!!!!!!!!量纲
    if ( kCode==1 ) then
        dt  =   timestep*length/sound           !Charles
        vref=   sound
        pref=   density*vref**2
    endif
    
    if ( kCode==2 ) then
        dt  =   timestep*length/(sound*mach)    !NSAWET
        vref=   sound*mach
        pref=   density*vref**2
    endif
!!!!!!!!!!!!!!!
    
    close(5)
    nt=2**KofFFT
    
    allocate(p(nt,CFDpointNum))
    allocate(avg(CFDpointNum))
    allocate(variance(CFDpointNum))
    allocate(pByFreq(nt,CFDpointNum))
    allocate(stepno(nt))
    allocate(timelist(nt))
    allocate(power_avg(nt,CFDpointNum))
    allocate(power_var_avg(nt,CFDpointNum))
    allocate(power_avg2(nt))
    allocate(power_var_avg2(nt))
    power_avg=0.0
    power_var_avg=0.0
    power_avg2=0.0
    power_var_avg2=0.0

    do iwindow=1,n_window*2-1
    
        p=0.0
        avg=0.0
        pbyfreq=0.0
        stepno=0.0
        timelist=0.0
    
        if(iwindow>1)n_startline=n_startline+nt/2
        n_endline = n_startline + nt
    
        write(*,*)'Window',iwindow, n_startline,n_endline
    
        call InputData
        call ApplyWindowFunction
        call FFTsource
    
        do i=1,nt/2+1
            do j=1,CFDpointNum
                power_avg(i,j)=power_avg(i,j)+pByFreq(i,j)/float(n_window*2-1)
                power_var_avg(i,j)=power_var_avg(i,j)+pByFreq(i,j)/float(n_window*2-1)/variance(j)
            enddo
        enddo
    enddo

    !OPEN(12,FILE='test_fft_window.dat',STATUS='UNKNOWN')
    OPEN(12,FILE=fName2,STATUS='UNKNOWN')
    !WRITE(12,*)'VARIABLES =  "freq" "fft_result"'
    
    if( kpsd == 1 ) then
        
        do i=1,nt/2+1
            do j=1,CFDpointNum
                power_avg2(i)=power_avg2(i)+abs(power_avg(i,j))
            enddo
            power_avg2(i)=power_avg2(i)/float(CFDpointNum)
            temp = float(i-1)/float(nt)/dt
            WRITE(12,'(1000e20.10)') temp, power_avg2(i), power_avg(i,1:CFDpointNum)
        enddo

    elseif( kpsd == 2 ) then
        
        do i=1,nt/2+1
            do j=1,CFDpointNum
                power_var_avg2(i)=power_var_avg2(i)+abs(power_var_avg(i,j))
            enddo
            power_var_avg2(i)=power_var_avg2(i)/float(CFDpointNum)
            temp = float(i-1)/float(nt)/dt
            WRITE(12,'(1000e20.10)') temp, temp*power_var_avg2(i),temp*power_var_avg(i,1:CFDpointNum)
        enddo
        
        !!!test if the integration is 1.0
        !!dlof(f) = log(f2)-log(f1)
        !!point 1
        !temp2=0.0
        !do  i=2,nt/2+1
        !    temp = float(i-1)/float(nt)/dt
        !    temp2 = temp2 + power_var_avg2(i)*temp*( log(temp+0.5*1.0/nt/dt)-log(temp-0.5*1.0/nt/dt) )
        !enddo
        !write(*,*) temp2
        !
        !temp2=0.0
        !do  i=1,nt/2+1
        !    temp = float(i-1)/float(nt)/dt
        !    temp2 = temp2 + power_avg(i,1)/variance(1)*(1.0/nt/dt)
        !enddo
        !write(*,*) temp2
        
    endif
    
    close(12)
    
    
    
    
    
end program main
    
subroutine FFTsource
    use Variables
    implicit none
    
    integer tmp

    real,dimension(nt):: FFTinputImag
    real,dimension(nt):: FFToutputReal
    real,dimension(nt):: FFToutputImag
    !print *,'!******* START: Subroutine: FFTsource ********!'
    
    do i=1,CFDpointNum
        tmp=i
        
    FFTinputImag=0
    FFToutputReal=0
    FFToutputImag=0
    !KofFFT:  nt=2^KofFFT
    call KKFFT(p(1:nt,i),FFTinputImag,nt,KofFFT,FFToutputReal,FFToutputImag,0,1)
    !FFToutputReal(1:nt)=FFToutputReal(1:nt)/float(nt)
    !FFToutputReal(1:nt)=FFToutputReal(1:nt)/float(nt)
    
    !pByFreq(1:nt,i)=(FFToutputReal(1:nt)+(0.,1.)*FFToutputImag(1:nt))*2./float(nt)*pi*2.
    !pByFreq(1:nt,i)=p(1:nt,i)**2 *(float(nt)*dt)
    !open(59,file='real.dat')
   ! write(59,*)'zone'
    do j=2,nt/2
        pByFreq(j,i)=(FFToutputReal(j)**2 + FFToutputImag(j)**2)*2.0 * (float(nt)*dt)  /float(nt)/float(nt)
        !pByFreq(j,i)=(FFToutputReal(j)**2 + FFToutputImag(j)**2)*2.0 /float(nt)/float(nt)
        !write(59,'(i5,5f20.10)') j,p(j,i), FFToutputReal(j),FFToutputImag(j)
    end do
    
    pByFreq(nt/2+1,i)=(FFToutputReal(nt/2+1)**2 + FFToutputImag(nt/2+1)**2) * (float(nt)*dt)  /float(nt)/float(nt)
    pByFreq(1,i)=(FFToutputReal(1)**2 + FFToutputImag(1)**2) * (float(nt)*dt)  /float(nt)/float(nt)
    
    !已经核对，没有错误
    !参考matlab计算功率谱密度的方法  https://www.mathworks.com/help/signal/ug/power-spectral-density-estimates-using-fft.html
    
    
    
    !pause
    enddo
    
!    deallocate(p)
        
    !print *,'!******* END: Subroutine: FFTsource ********!'
    RETURN
end subroutine

    
    
subroutine ApplyWindowFunction
    use Variables
    implicit none

    real,dimension(nt):: windowT

    !print *,'!******* START: Subroutine: ApplyWindowFunction ********!'
    !FWH_core中计算p'时没有计算0频率，因此时域源项是否减掉平均值都不影响
    !注意使用窗函数后有energy preserving，需要修正回来
    !!!注意，为了FFT，时域点个数为2^K,但书上的窗函数公式是针对有奇数个时域点的
    if (Window==0) then
        print *,'Window==0, No window function is used.'
    elseif (Window==1) then
        print *,'Window==1, Bartlett window function is used.'
        do j=1,nt/2
            windowT(j)=(2.*float(j))/(float(nt))
            windowT(j+nt/2)=2.-(2.*(float(j)-1.+float(nt)/2.))/(float(nt))
        enddo
        
        do i=1,CFDpointNum
        p(:,i)=p(:,i)*windowT*1.732
        enddo
        
    elseif (Window==2) then
        print *,'Window==2, Hanning window function is used.'
        do j=1,nt
            windowT(j)=1./2.*(1.-cos(2.*pi*(float(j)-1.)/(float(nt)-1.)))
        enddo
        
        do i=1,CFDpointNum
        p(:,i)=p(:,i)*windowT*1.633
        enddo        
        
    elseif (Window==3) then
        print *,'Window==3, Hamming window function is used.'
        do j=1,nt
            windowT(j)=0.54-0.46*cos(2.*pi*(float(j)-1.)/(float(nt)-1.))
        enddo
        
        do i=1,CFDpointNum
        p(:,i)=p(:,i)*windowT*1.586
        enddo
        
    elseif (Window==4) then
        print *,'Window==4, Blackman window function is used.'
        do j=1,nt
            windowT(j)=0.42-0.5*cos(2.*pi*(float(j)-1.)/(float(nt)-1.))+0.08*cos(4.*pi*(float(j)-1.)/(float(nt)-1.))
        enddo
        
        do i=1,CFDpointNum
        p(:,i)=p(:,i)*windowT*1.812
        enddo        
        
    else
        print *,'Input Window is invalid: ',Window
    endif
        
    
    !print *,'!******* END: Subroutine: ApplyWindowFunction ********!'
    end subroutine 

subroutine InputData
    use Variables
    implicit none
    integer tmp48
   ! print *,'Inputing...'
    open(11,FILE=filename,status='old',action='read')
    avg=0.0
    variance=0.0
    
    do i=1,n_startline-1
        read(11,*)
    enddo
    
    i=0
    do 
        i=i+1
        if( kCode == 1 )  read( 11,*,end=100,err=101 )  stepno(i),timelist(i),tmp48,p(i,1:CFDpointNum)  !Charles
        if( kCode == 2 )  read( 11,*,end=100,err=101 )  stepno(i),                  p(i,1:CFDpointNum)  !NSAWET
        
        if(i>1) then
        if(stepno(i)<=stepno(i-1))then
            i=i-1
            write(*,*)i,stepno(i)
           ! pause
            cycle
        endif
        endif
        
        if(kVar==1) then
            do j=1,CFDpointNum
                p(i,j) = p(i,j) * pref          !变成有量纲的Pa
                avg(j) = avg(j)+p(i,j)          !计算平均值
            enddo
        elseif(kVar==2) then
            do j=1,CFDpointNum
                p(i,j) = p(i,j) * vref          !变成有量纲的m/s
                avg(j) = avg(j)+p(i,j)          !计算平均值
            enddo
        endif
        
        if(i==nt)exit
    enddo
    
    do j=1,CFDpointNum
         avg(j)=avg(j)/float(nt)            !计算平均值
    enddo
    
    do i=1,nt
        do j=1,CFDpointNum
            p(i,j)=p(i,j)-avg(j)                !平均值扣除
            variance(j) = variance(j) + p(i,j)*p(i,j)  !方差
        enddo
    enddo
    
    do j=1,CFDpointNum
         variance(j)=variance(j)/float(nt)            !方差
    enddo
    
    
    print *,nt,stepno(nt),timelist(nt),avg(1:CFDpointNum),p(nt,:)
    write(*,*)
    close(11)
  !  open(55,file='pp.dat',access='append')
  !  write(55,*)'zone'
  !  do i=1,nt
  !      write(55,*)float(i)*dt,p(i,1)
  !  end do
  !  close(55)
    
    return
100 continue
101 continue
    stop 'Error read file!'
    
end subroutine 
    
 
    
   SUBROUTINE KKFFT(PR,PI,N,K,FR,FI,L,IL)
	real,dimension(N) :: PR,PI,FR,FI
	real P,Q,S,VR,VI,PODDR,PODDI
	DO 20 IT=0,N-1
	  M=IT
	  IS=0
	  DO 10 I=0,K-1
	    J=M/2
	    IS=2*IS+(M-2*J)
	    M=J
10	  CONTINUE
	  FR(IT+1)=PR(IS+1)
	  FI(IT+1)=PI(IS+1)
20	CONTINUE
	PR(1)=1.0
	PI(1)=0.0
	PR(2)=COS(6.283185306/N)
	PI(2)=-SIN(6.283185306/N)
	IF (L.NE.0) PI(2)=-PI(2)
	DO 30 I=3,N
	  P=PR(I-1)*PR(2)
	  Q=PI(I-1)*PI(2)
	  S=(PR(I-1)+PI(I-1))*(PR(2)+PI(2))
	  PR(I)=P-Q
	  PI(I)=S-P-Q
30	CONTINUE
	DO 40 IT=0,N-2,2
	  VR=FR(IT+1)
	  VI=FI(IT+1)
	  FR(IT+1)=VR+FR(IT+2)
	  FI(IT+1)=VI+FI(IT+2)
	  FR(IT+2)=VR-FR(IT+2)
	  FI(IT+2)=VI-FI(IT+2)
40	CONTINUE
	M=N/2
	NV=2
	DO 70 L0=K-2,0,-1
	  M=M/2
	  NV=2*NV
	  DO 60 IT=0,(M-1)*NV,NV
	  DO 60 J=0,(NV/2)-1
	    P=PR(M*J+1)*FR(IT+J+1+NV/2)
	    Q=PI(M*J+1)*FI(IT+J+1+NV/2)
	    S=PR(M*J+1)+PI(M*J+1)
	    S=S*(FR(IT+J+1+NV/2)+FI(IT+J+1+NV/2))
	    PODDR=P-Q
	    PODDI=S-P-Q
	    FR(IT+J+1+NV/2)=FR(IT+J+1)-PODDR
	    FI(IT+J+1+NV/2)=FI(IT+J+1)-PODDI
	    FR(IT+J+1)=FR(IT+J+1)+PODDR
	    FI(IT+J+1)=FI(IT+J+1)+PODDI
60	  CONTINUE
70	CONTINUE


	IF (L.NE.0) THEN
	  DO 80 I=1,N
	    FR(I)=FR(I)/N
	    FI(I)=FI(I)/N
80	  CONTINUE
	END IF
	IF (IL.NE.0) THEN
	  DO 90 I=1,N
	    PR(I)=SQRT(FR(I)*FR(I)+FI(I)*FI(I))
	    PI(I)=ATAN(FI(I)/FR(I))*360.0/6.283185306
90	  CONTINUE
	END IF
	RETURN
    END
    
