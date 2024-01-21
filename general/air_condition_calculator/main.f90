Program Air_Condition_Calculator
implicit none
!==================================================
!Programmed by Raccoon
!本程序功能：
!（1）给定基于海平面高度：按照标准大气模型计算当地空气状态（密度，压力，温度，声速，动力学粘度，运动学粘度）
!（2）给定特征长度，速度（或基于上面计算的声速给出马赫数），运动学粘度（或使用上面计算的运动学粘度），yplus：计算流动雷诺数，yplus数对应的壁面第一层网格距离
!==================================================

!==================================================
!物理变量声明
real(kind=8)::                      h                   !高度
real(kind=8)::                      p                   !压强
real(kind=8)::                      rho                 !密度
real(kind=8)::                      t                   !温度
real(kind=8)::                      a                   !声速
real(kind=8)::                      nu                  !运动学粘度
real(kind=8)::                      mu                  !动力学粘度
real(kind=8)::                      l                   !特征长度
real(kind=8)::                      mach                !来流马赫数
real(kind=8)::                      v                   !来流速度
real(kind=8)::                      re                  !雷诺数
real(kind=8)::                      lturb               !层流/湍流标识符，=0表示层流，=1表示湍流
real(kind=8)::                      yplus               !y+
real(kind=8)::                      d                   !到壁面距离
!程序中使用变量声明
integer(kind=4)::                   lcase               !case选择标志，lcase=1表示空气状态求解，lcase=2表示雷诺数壁面距离求解
real(kind=8)::                      temp                !读入数据时的临时变量

!==================================================

!==================================================
!首先将lcase置1，至少进行1次空气状态计算
lcase = 1
do while (lcase>0)
    if (lcase==1) then
        !读入高度
        write(*,*)'Please input the height [m] :'
        read(*,*)h
        !根据标准大气模型计算该高度处状态
        call isa_model(h, p, rho, t, a, nu, mu)
        !输出
        write(*,*)
        write(*,*)'Pressure [Pa]               :', p
        write(*,*)'Density [kg/m^3]            :', rho
        write(*,*)'Temperature [K]             :', t
        write(*,*)'Air Speed [m/s]             :', a
        write(*,*)'Dynamic Viscosity [Pa*s]    :', mu
        write(*,*)'Kinematic Viscosity [m^2/s] :', nu
        write(*,*)
        !再次读入lcase
        write(*,*)'Input 1 for air status calculation.'
        write(*,*)'Input 2 to calculate reynolds number and estimate wall distance.'
        write(*,*)'Input 0 to quit.'
        read(*,*)lcase
        write(*,*)
    else if (lcase==2) then
        !读入雷诺数的特征长度
        write(*,*)'Input the characteristic length [m] :'
        read(*,*)l
        !读入马赫数/速度
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
        !读入新的运动学粘度（如果需要）
        write(*,*)'Input kinematic viscocity [m^2/s] :'
        write(*,*)'(If using the nu above (', nu, '[m2/s]), input -1)'
        read(*,*)temp
        if (temp>0) then
            nu = temp
        endif
        !读入需要的y+
        write(*,*)'Input desired y+ [1] :'
        read(*,*)yplus
        !读入层流或湍流的标识量
        write(*,*)'Input lturb (0 for laminar, 1 for turbulence) :'
        read(*,*)lturb
        !计算雷诺数
        re = v * l / nu
        !根据层流或湍流计算上面yplus对应的壁面距离
        call wall_distance(nu, v, re, lturb, yplus, d)
        !输出
        write(*,*)
        write(*,*)'mach number                 :', mach
        write(*,*)'Velocity [m/s]              :', v
        write(*,*)'Length [m]                  :', l
        write(*,*)'Kinematic viscocity [m^2/s] :', nu
        write(*,*)'Reynolds number             :', re
        write(*,*)'wall distance [m] for y+ = ',yplus
        write(*,*)'                            :', d
        write(*,*)
        !再次读入lcase
        write(*,*)'Input 1 for air status calculation.'
        write(*,*)'Input 2 to calculate reynolds number and estimate wall distance.'
        write(*,*)'Input 0 to quit.'
        read(*,*)lcase
        write(*,*)
    end if
enddo

end program
        
        
