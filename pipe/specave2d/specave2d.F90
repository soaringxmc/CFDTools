program main
implicit none
!---------- basic variables
integer:: m1m,m2m,m3m,m1,m2,m3,m1mh,m3mh,i,j,k,m,n
!real(kind=8),allocatable,dimension(:,:,:):: spect,specz,spectav,speczav
real(kind=8),allocatable,dimension(:,:):: spec2d,spec2dav
real(kind=8),allocatable,dimension(:):: lm1,lm2,lm3
real(kind=8):: pi,laxis,tw,utau,qw,retau,eps=1.0e-30
!---------- intermediate variables
integer:: kstart,kend,kstep,kkn,kk1,kk2
logical:: alive
real(kind=8):: tmp1,tmp2,tmp3
character(4) :: naitav
character(1) :: navar,naplan
!--------- end of declaration -------

pi=4.*atan(1.d0)

open(10,file='../bou.in',action='read')
read(10,*)
read(10,*) m1,m2,m3
do i=1,3
read(10,*)
enddo
read(10,*) laxis
close(10)

open(11,file='stats.in',action='read')
read(11,*) kstart,kend
close(11)


m1m=m1-1; m2m=m2-1; m3m=m3-1; m1mh=m1m/2+1; m3mh=m3m/2+1
!allocate(spect(m1mh,m2m,5),specz(m3mh,m2m,5),spectav(m1mh,m2m,5),speczav(m3mh,m2m,5))
allocate(spec2d(m1mh,m3mh),spec2dav(m1mh,m3mh))
allocate(lm1(m1m),lm2(m2m),lm3(m3m))

open(13,file='../flowprop.dat',action='read')
read(13,*) tmp1, retau, tmp1, utau, tmp1, tmp1, tmp1, tmp1, &
tmp1, tmp1, tmp1, tmp1, tmp1, tmp1, tmp1, tmp1, qw
close(13)

!print *,'utau,ttau',utau,ttau
!if (ttau.eq.0.) ttau=1.

do i=1,m3m
lm3(i)=laxis/(i-1.d0 + eps)*retau
enddo
tmp2=datan(1.d0)*8.d0
do i=1,m1m
lm1(i)=tmp2/(i-1.d0+eps)*retau
enddo

open(12,file='../radcor.out',action='read')
do i=1,m2m
read(12,*) tmp1,tmp1,lm2(i)
enddo
close(12)

do i=1,m2m
lm2(i)=(1.d0-lm2(i))*retau
enddo

!----- 2D spectra -----------

spec2dav=0.d0 ;kkn=0

read(*,*) kk1,kk2
!do kk2=1,4
!do kk1=1,5
write(navar,'(I1.1)') kk1
write(naplan,'(I1.1)') kk2
do kstep=kstart,kend
write(naitav,'(I4.4)') kstep
inquire(file='../spec2d'//naplan//'_'//navar//'_'//naitav//'.bin',exist=alive)
if(alive) then
print*, kstep,'var=',kk1,'plane=',kk2
kkn=kkn+1

open(21,file='../spec2d'//naplan//'_'//navar//'_'//naitav//'.bin',form='unformatted',action='read')
read(21) spec2d !(:,:,k)
close(21)

spec2dav=spec2dav+spec2d
endif

enddo

tmp1=1.d0/kkn

spec2dav=spec2dav*tmp1; 

if(kk1.eq.1.or.kk1.eq.2.or.kk1.eq.3) tmp2=utau**2
if(kk1.eq.4) tmp2=qw**4
if(kk1.eq.5) tmp2=utau**4

!open(22,file='spec2d-av'//naplan//'_'//navar//'_tecplot.dat')
!write(22,*) 'variables = log10(lx+),log10(lz+),var'
!write(22,*) 'zone i=',m1mh-2,'j=',m3mh-2
!
!do j=2,m3mh-1
!do i=2,m1mh-1
!write(22,'(20g20.9)') dlog10(lm1(i)),dlog10(lm3(j)),(i-1.d0)*(j-1.d0)*spec2dav(i,j)/tmp2
!enddo
!enddo
!close(22)

open(23,file='spec2d-av-plot3d.xyz',form='unformatted')
write(23) (m1mh-2),(m3mh-2)
write(23) ((lm1(i),i=2,m1mh-1),j=2,m3mh-1),((lm3(j),i=2,m1mh-1),j=2,m3mh-1)
close(23)

open(23,file='spec2d-av'//naplan//'_'//navar//'.q',form='unformatted')
write(23) (m1mh-2),(m3mh-2),1
!write(23) 0.0,0.0,0.0,0.0
write(23) (((i-1.d0)*(j-1.d0)*spec2dav(i,j)/tmp2,i=2,m1mh-1),j=2,m3mh-1)
close(23)

open(23,file='spec-av-plot3d.nam')
write(23,*) 'uth'
write(23,*) 'ur'
write(23,*) 'uz'
write(23,*) 'th'
write(23,*) 'p'
close(23)

!enddo
!enddo
end
