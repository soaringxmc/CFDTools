subroutine isa_model(h, p, rho, t, a, nu, mu)
implicit none
!==================================================
!ʹ�ñ�׼����ģ�ͣ�International Standard Gas������߶�h���Ĵ���״̬
!�ο���https://max.book118.com/html/2019/0115/5342133024002002.shtm
!==================================================

!==================================================
!���������������
real(kind=8)::                      h                   !�߶�
real(kind=8)::                      p                   !ѹǿ
real(kind=8)::                      rho                 !�ܶ�
real(kind=8)::                      t                   !�¶�
real(kind=8)::                      a                   !����
real(kind=8)::                      nu                  !�˶�ѧճ��
real(kind=8)::                      mu                  !����ѧճ��
!��׼����ģ��ʹ�ò���
real(kind=8),parameter::            t0=288.15           !��ƽ�洦���¶�[K]
real(kind=8),parameter::            t1=216.65           !����11000�״����¶�[K]
real(kind=8),parameter::            p0=101325.0         !��ƽ�洦��ѹ��[Pa]
real(kind=8),parameter::            p1=22632.05         !����11000�״���ѹ��[Pa]
real(kind=8),parameter::            alpha=0.0065        !ѹ��ģ��ϵ����alpha [K/m]
real(kind=8),parameter::            g0=9.80665          !ѹ��ģ��ϵ����g0 [m^2/s]
real(kind=8),parameter::            r=287.053           !ѹ��ģ��ϵ����r [J/kg/K]
real(kind=8),parameter::            gamma=1.4           !ѹ��ģ��ϵ����gamma [1]
real(kind=8),parameter::            mu0=1.7894e-5       !��������ʽϵ����mu0 [Pa*s]
real(kind=8),parameter::            b=110.4             !��������ʽϵ����b [K]
!==================================================
if (h<0) then
    h = 0.0
endif
!�����¶�
if (h<11000) then
    t = t0 + (t1-t0)/11000*h
else
    t = t1
endif

!��ѹ
if (h<11000) then
    p = p0 * (1.0 - alpha/t0*h)**(g0/alpha/r)
else
    p = p1 * exp(-g0*(h-11000)/r/t1)
endif

!�ܶ�
rho = p/r/t

!����
a = sqrt(gamma*r*t)

!����ѧճ�ȣ���������ʽ��
mu = mu0 * (t/t0)**1.5 * (t0+b)/(t+b)

!�˶�ѧճ��
nu = mu/rho

end subroutine

