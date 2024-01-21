subroutine wall_distance(nu, v, re, lturb, yplus, d)
implicit none
!==================================================
!使用布拉修斯解（层流）或平板壁面七分之一律（湍流）计算给定yplus的壁面距离d
!参考：布拉修斯解：吴子牛《空气动力学》P91
!     七分之一律:Frank M. White's Fluid Mechanics 5th edition, P467
!==================================================

!==================================================
!输入输出变量声明
real(kind=8)::                      nu                  !运动学粘度
real(kind=8)::                      v                   !来流速度
real(kind=8)::                      re                  !雷诺数
real(kind=8)::                      lturb               !层流/湍流标识符，=0表示层流，=1表示湍流
real(kind=8)::                      yplus               !y+
real(kind=8)::                      d                   !到壁面距离
!程序中使用临时变量
real(kind=8)::                      cf                  !壁面摩擦系数
real(kind=8)::                      tau                 !壁面切应力
real(kind=8)::                      utau                !壁面磨擦速度
!==================================================

if (lturb==0) then
    !层流
    cf = 0.664/re**0.5
else if (lturb==1) then
    !湍流
    cf = 0.026/re**(1.0/7.0)
endif
tau = cf*v**2/2.0
utau = sqrt(tau)
d = yplus*nu/utau

end subroutine 