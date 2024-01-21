subroutine wall_distance(nu, v, re, lturb, yplus, d)
implicit none
!==================================================
!ʹ�ò�����˹�⣨��������ƽ������߷�֮һ�ɣ��������������yplus�ı������d
!�ο���������˹�⣺����ţ����������ѧ��P91
!     �߷�֮һ��:Frank M. White's Fluid Mechanics 5th edition, P467
!==================================================

!==================================================
!���������������
real(kind=8)::                      nu                  !�˶�ѧճ��
real(kind=8)::                      v                   !�����ٶ�
real(kind=8)::                      re                  !��ŵ��
real(kind=8)::                      lturb               !����/������ʶ����=0��ʾ������=1��ʾ����
real(kind=8)::                      yplus               !y+
real(kind=8)::                      d                   !���������
!������ʹ����ʱ����
real(kind=8)::                      cf                  !����Ħ��ϵ��
real(kind=8)::                      tau                 !������Ӧ��
real(kind=8)::                      utau                !����ĥ���ٶ�
!==================================================

if (lturb==0) then
    !����
    cf = 0.664/re**0.5
else if (lturb==1) then
    !����
    cf = 0.026/re**(1.0/7.0)
endif
tau = cf*v**2/2.0
utau = sqrt(tau)
d = yplus*nu/utau

end subroutine 