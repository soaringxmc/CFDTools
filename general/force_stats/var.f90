    module Var
    
    implicit none
    
    real    ,allocatable, dimension( :, : ) ::  force, lift, drag, mome, liftMov, dragMov, momeMov
    real                , dimension( 3    ) ::  liftRMS, dragRMS, momeRMS
    real                                    ::  rhoRef, vRef, AoA, arRef
    integer ,allocatable, dimension( :    ) ::  stepNo
    integer                                 ::  stepSta, nSta, nStep, iSolver
    character(len=1000)                     ::  fileIn, fileUnst, fileStat
    
    end module Var