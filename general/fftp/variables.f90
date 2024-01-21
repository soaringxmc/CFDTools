module Variables
    implicit none
    real,allocatable,dimension(:,:):: p,power_avg,power_var_avg
    real,allocatable,dimension(:):: avg,variance,power_avg2,power_var_avg2
    integer,allocatable,dimension(:):: stepno
    real,allocatable,dimension(:):: timelist
    complex*16,allocatable,dimension(:,:):: pByFreq
    
    integer Window
    integer KofFFT, nt, CFDpointNum
    integer i, j, k, ii, jj, kk
    real pi
    integer n_startline,n_points,n_window,n_endline
    real dt
    real sound, mach, density,timestep,length
    real pref,vref,gamma
    character*100 filename
    character*100 fName2  
    integer kVar,kCode,kPSD
    
    
    end module Variables   