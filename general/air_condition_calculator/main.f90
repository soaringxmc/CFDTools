Program Air_Condition_Calculator
implicit none
!==================================================
!Programmed by Raccoon
!�������ܣ�
!��1���������ں�ƽ��߶ȣ����ձ�׼����ģ�ͼ��㵱�ؿ���״̬���ܶȣ�ѹ�����¶ȣ����٣�����ѧճ�ȣ��˶�ѧճ�ȣ�
!��2�������������ȣ��ٶȣ�����������������ٸ�������������˶�ѧճ�ȣ���ʹ�����������˶�ѧճ�ȣ���yplus������������ŵ����yplus����Ӧ�ı����һ���������
!==================================================

!==================================================
!�����������
real(kind=8)::                      h                   !�߶�
real(kind=8)::                      p                   !ѹǿ
real(kind=8)::                      rho                 !�ܶ�
real(kind=8)::                      t                   !�¶�
real(kind=8)::                      a                   !����
real(kind=8)::                      nu                  !�˶�ѧճ��
real(kind=8)::                      mu                  !����ѧճ��
real(kind=8)::                      l                   !��������
real(kind=8)::                      mach                !���������
real(kind=8)::                      v                   !�����ٶ�
real(kind=8)::                      re                  !��ŵ��
real(kind=8)::                      lturb               !����/������ʶ����=0��ʾ������=1��ʾ����
real(kind=8)::                      yplus               !y+
real(kind=8)::                      d                   !���������
!������ʹ�ñ�������
integer(kind=4)::                   lcase               !caseѡ���־��lcase=1��ʾ����״̬��⣬lcase=2��ʾ��ŵ������������
real(kind=8)::                      temp                !��������ʱ����ʱ����

!==================================================

!==================================================
!���Ƚ�lcase��1�����ٽ���1�ο���״̬����
lcase = 1
do while (lcase>0)
    if (lcase==1) then
        !����߶�
        write(*,*)'Please input the height [m] :'
        read(*,*)h
        !���ݱ�׼����ģ�ͼ���ø߶ȴ�״̬
        call isa_model(h, p, rho, t, a, nu, mu)
        !���
        write(*,*)
        write(*,*)'Pressure [Pa]               :', p
        write(*,*)'Density [kg/m^3]            :', rho
        write(*,*)'Temperature [K]             :', t
        write(*,*)'Air Speed [m/s]             :', a
        write(*,*)'Dynamic Viscosity [Pa*s]    :', mu
        write(*,*)'Kinematic Viscosity [m^2/s] :', nu
        write(*,*)
        !�ٴζ���lcase
        write(*,*)'Input 1 for air status calculation.'
        write(*,*)'Input 2 to calculate reynolds number and estimate wall distance.'
        write(*,*)'Input 0 to quit.'
        read(*,*)lcase
        write(*,*)
    else if (lcase==2) then
        !������ŵ������������
        write(*,*)'Input the characteristic length [m] :'
        read(*,*)l
        !���������/�ٶ�
        write(*,*)'Input the mach number based on the air speed above (',a,' [m/s]) :'
        write(*,*)'(If want to directly input the velocity, input -1 for this value)'
        read(*,*)temp
        if (temp<0) then
            write(*,*)'Input the freestream velocity [m/s] :'
            read(*,*)v
            mach = v/a
        else
            mach = temp
            v = mach * a
        endif
        !�����µ��˶�ѧճ�ȣ������Ҫ��
        write(*,*)'Input kinematic viscocity [m^2/s] :'
        write(*,*)'(If using the nu above (', nu, '[m2/s]), input -1)'
        read(*,*)temp
        if (temp>0) then
            nu = temp
        endif
        !������Ҫ��y+
        write(*,*)'Input desired y+ [1] :'
        read(*,*)yplus
        !��������������ı�ʶ��
        write(*,*)'Input lturb (0 for laminar, 1 for turbulence) :'
        read(*,*)lturb
        !������ŵ��
        re = v * l / nu
        !���ݲ�����������������yplus��Ӧ�ı������
        call wall_distance(nu, v, re, lturb, yplus, d)
        !���
        write(*,*)
        write(*,*)'mach number                 :', mach
        write(*,*)'Velocity [m/s]              :', v
        write(*,*)'Length [m]                  :', l
        write(*,*)'Kinematic viscocity [m^2/s] :', nu
        write(*,*)'Reynolds number             :', re
        write(*,*)'wall distance [m] for y+ = ',yplus
        write(*,*)'                            :', d
        write(*,*)
        !�ٴζ���lcase
        write(*,*)'Input 1 for air status calculation.'
        write(*,*)'Input 2 to calculate reynolds number and estimate wall distance.'
        write(*,*)'Input 0 to quit.'
        read(*,*)lcase
        write(*,*)
    end if
enddo

end program
        
        
