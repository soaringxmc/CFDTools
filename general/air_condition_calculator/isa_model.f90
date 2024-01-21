subroutine isa_model(h, p, rho, t, a, nu, mu)
implicit none
!==================================================
!使用标准大气模型（International Standard Gas）计算高度h处的大气状态
!参考：https://max.book118.com/html/2019/0115/5342133024002002.shtm
!==================================================

!==================================================
!输入输出变量声明
real(kind=8)::                      h                   !高度
real(kind=8)::                      p                   !压强
real(kind=8)::                      rho                 !密度
real(kind=8)::                      t                   !温度
real(kind=8)::                      a                   !声速
real(kind=8)::                      nu                  !运动学粘度
real(kind=8)::                      mu                  !动力学粘度
!标准大气模型使用参数
real(kind=8),parameter::            t0=288.15           !海平面处的温度[K]
real(kind=8),parameter::            t1=216.65           !海拔11000米处的温度[K]
real(kind=8),parameter::            p0=101325.0         !海平面处的压力[Pa]
real(kind=8),parameter::            p1=22632.05         !海拔11000米处的压力[Pa]
real(kind=8),parameter::            alpha=0.0065        !压力模型系数：alpha [K/m]
real(kind=8),parameter::            g0=9.80665          !压力模型系数：g0 [m^2/s]
real(kind=8),parameter::            r=287.053           !压力模型系数：r [J/kg/K]
real(kind=8),parameter::            gamma=1.4           !压力模型系数：gamma [1]
real(kind=8),parameter::            mu0=1.7894e-5       !萨特兰公式系数：mu0 [Pa*s]
real(kind=8),parameter::            b=110.4             !萨特兰公式系数：b [K]
!==================================================
if (h<0) then
    h = 0.0
endif
!大气温度
if (h<11000) then
    t = t0 + (t1-t0)/11000*h
else
    t = t1
endif

!气压
if (h<11000) then
    p = p0 * (1.0 - alpha/t0*h)**(g0/alpha/r)
else
    p = p1 * exp(-g0*(h-11000)/r/t1)
endif

!密度
rho = p/r/t

!声速
a = sqrt(gamma*r*t)

!动力学粘度（萨特兰公式）
mu = mu0 * (t/t0)**1.5 * (t0+b)/(t+b)

!运动学粘度
nu = mu/rho

end subroutine

