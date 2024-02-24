!  plot3d.f90 
!
!  FUNCTIONS:
!  plot3d - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: plot3d
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program plot3d
!
    use hdf5
    include 'mpif.h'   !defined variables
!
    integer, parameter :: dp = 8
    real(dp), allocatable, dimension(:,:,:) :: q1,q1m,u1m,uxm
    real(dp), allocatable, dimension(:,:,:) :: q2,q2m,u2m,uym
    real(dp), allocatable, dimension(:,:,:) :: q3,q3m,u3m
    real(dp), allocatable, dimension(:,:,:) :: dens,densm
    real(dp), allocatable, dimension(:,:,:) :: pr,prm
    real(dp), allocatable, dimension(:,:,:) :: xm,ym,zm
    real(dp), allocatable, dimension(:)     :: rc,rm,zz,zzm,thetac,thetam
    real(dp), allocatable, dimension(:,:)   :: planetzl(:,:),planetz(:,:,:)
    real(dp), allocatable, dimension(:,:)   :: planetrl(:,:),planetr(:,:)
!    
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
!
    open(11,file='bou.in')
    read(11,*)
    read(11,*) m1,m2,m3
    read(11,*)
    read(11,*)
    read(11,*)
    read(11,*) alx3d
    close(11)
!
    m1m = m1 - 1
    m2m = m2 - 1
    m3m = m3 - 1
    kstart = 1 + (m3m/nproc)*nrank 
    kend   = (m3m/nproc)*(nrank+1)
    jplane = m2m-20
    l3d = 0
!
    if (mod(m3m,nproc).ne.0) then
      print *, 'error: mod(m3m,nproc).ne.0'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
    endif
!    
    allocate(planetz(m1m,1,1-1:m3m+1))
    allocate(planetzl(m1,m3m))
    allocate(planetrl(m1,m2m))
    allocate(planetr(m1m,m2m))
!
    call h5open_f(ihdf_error)
!
!   grid
!
    allocate(rc (1:m2m))
    allocate(rm (1:m2m))
    allocate(zz (kstart:kend))
    allocate(zzm(kstart:kend))
    allocate(thetac(1:m1))
    allocate(thetam(1:m1))
!
    do k=kstart,kend
      zz(k) = alx3d*real(k-1,dp)/real(m3m,dp)
      zzm(k)= alx3d*(real(k-1,dp)+0.5d0)/real(m3m,dp) 
    enddo
!
    open(20,file='radcor.out')
    do j=1,m2m
      read(20,*) itmp,rc(j),rm(j)
    enddo
    close(20)
!
    pi  = 3.141592653589
    dx1 = 1.0d0/(2.d0*pi/real(m1m, dp))
    do i=1,m1m
      thetac(i) = real(i-1,dp)/dx1
      thetam(i) = (real(i-1,dp)+0.5d0)/dx1
    enddo
    thetac(m1) = real(m1-1,dp)/dx1
!
    call mpi_read_continua_surf(m1m,m2m,m3,1,m3m,jplane,3,planetz)
!
    if (nrank.eq.0) then
      planetzl(1:m1m,1:m3m) = planetz(1:m1m,1,1:m3m)
      planetzl(m1,1:m3m) = planetz(1,1,1:m3m)
      ! utau = 1.9445827669489081E-002
      utau = 0.5
      planetzl = planetzl/utau
      open(11,file='planetz.q',form='unformatted')
      write(11) m1,1,m3m,1
      write(11) (((planetzl(i,k),i=1,m1),j=1,1),k=1,m3m)
      close(11)
      open(11,file='planetz.x',form='unformatted')
      write(11) m1,1,m3m
      jj = jplane
      ! thetac used here for consistency with extract_plane_uz.f
      write(11) (((rm(j)*cos(thetac(i)),i=1,m1),j=jj,jj),k=1,m3m), &
                (((rm(j)*sin(thetac(i)),i=1,m1),j=jj,jj),k=1,m3m), &
                (((zz(k),i=1,m1),j=jj,jj),k=1,m3m)
      close(11)
      open(11,file='planetz_unrolled.x',form='unformatted')
      rewind (11)
      write  (11) m1,1,m3m
      write  (11) (((thetac(i),i=1,m1),j=m2m,m2m),k=1,m3m), &
                  (((0.,i=1,m1),j=m2m,m2m),k=1,m3m), &
                  (((zz(k),i=1,m1),j=m2m,m2m),k=1,m3m)
      close(11)
  !
      call mpi_read_continua_surf_tr(m1m,m2m,m3,1,1,1,3,planetr)
      planetrl(1:m1m,1:m2m) = planetr(1:m1m,1:m2m)
      planetrl(m1,1:m2m) = planetr(1,1:m2m)
      planetrl = planetrl/utau
      open(11,file='planetr.x',form='unformatted')
      write(11) m1,m2m,1
      write(11) (((rm(j)*cos(thetac(i)),i=1,m1),j=1,m2m),k=1,1), &
                (((rm(j)*sin(thetac(i)),i=1,m1),j=1,m2m),k=1,1), &
                (((zz(k),i=1,m1),j=1,m2m),k=1,1)
      close(11)
      open(11,file='planetr.q',form='unformatted')
      write(11) m1,m2m,1,1
      write(11) ((planetrl(i,j),i=1,m1),j=1,m2m)
      close(11)
    end if
!
    if (l3d.eq.0) then
      call MPI_Barrier(MPI_COMM_WORLD)
      stop
    end if
!
!   continue if 3D field output
!
!   the dimensions are consistent with pipeflow, 
!   note: each process reads all the q3/Uz data for 
!   its continuity between consecutive ranks
    allocate(q1  (1:m1m,1:m2m,kstart-1:kend+1))
    allocate(q2  (1:m1 ,1:m2 ,kstart-1:kend+1))
    allocate(q3  (1:m1m,1:m2m,     1-1: m3m+1))
    allocate(dens(1:m1m,1:m2m,kstart-1:kend+1))
    allocate(pr  (1:m1m,1:m2m,kstart-1:kend))
!  
    allocate(q1m  (1:m1,1:m2m,kstart-1:kend+1))
    allocate(q2m  (1:m1,1:m2m,kstart-1:kend+1))
    allocate(q3m  (1:m1,1:m2m,kstart-1:kend+1))
    allocate(u1m  (1:m1,1:m2m,kstart-1:kend+1))
    allocate(u2m  (1:m1,1:m2m,kstart-1:kend+1))
    allocate(u3m  (1:m1,1:m2m,kstart-1:kend+1))   
    allocate(densm(1:m1,1:m2m,kstart-1:kend+1))
    allocate(prm  (1:m1,1:m2m,kstart-1:kend+1))
    allocate(uxm  (1:m1,1:m2m,kstart-1:kend+1))
    allocate(uym  (1:m1,1:m2m,kstart-1:kend+1))
    allocate(xm   (1:m1,1:m2m,kstart-1:kend+1))
    allocate(ym   (1:m1,1:m2m,kstart-1:kend+1))
    allocate(zm   (1:m1,1:m2m,kstart-1:kend+1))
!
    do k=kstart,kend
    do j=1,m2m
    do i=1,m1m
      zm(i,j,k) = zzm(k)
      xm(i,j,k) = rm(j)*cos(thetam(i))
      ym(i,j,k) = rm(j)*sin(thetam(i))
    enddo
    enddo
    enddo
!   call mpi_read_continua(m1m,m2m,m3,kstart,kend,4,dens)
    call mpi_read_continua(m1m,m2m,m3,kstart,kend,1,q1)
    call mpi_read_continua(m1 ,m2 ,m3,kstart,kend,2,q2)
    call mpi_read_continua(m1m,m2m,m3,     1, m3m,3,q3)
    call mpi_read_continua(m1m,m2m,m3,kstart,kend,5,pr)
!
    do k = kstart,kend
    do j = 1,m2m
    do i = 1,m1m-1
      q1m (i,j,k) = 0.5*(q1 (i,j,k) + q1 (i+1,j,k))
    enddo
      q1m (m1m,j,k) = 0.5*(q1 (m1m,j,k) + q1 (1,j,k))
    enddo
    enddo
!
    do k = kstart,kend
    do j = 1,m2m
    do i = 1,m1m
      q2m (i,j,k) = 0.5*(q2 (i,j,k) + q2 (i,j+1,k))
    enddo
    enddo
    enddo
!
    do k = kstart,kend-1
    do j = 1,m2m
    do i = 1,m1m
      q3m (i,j,k) = 0.5*(q3 (i,j,k) + q3 (i,j,k+1))
    enddo
    enddo
    enddo
    if(nrank.eq.nproc-1) then
      do k = kend,kend
      do j = 1,m2m
      do i = 1,m1m
        q3m (i,j,k) = 0.5*(q3 (i,j,k) + q3 (i,j,1))
      enddo
      enddo
      enddo
    else
      do k = kend,kend
      do j = 1,m2m
      do i = 1,m1m
        q3m (i,j,k) = 0.5*(q3 (i,j,k) + q3 (i,j,k+1))
      enddo
      enddo
      enddo
    endif
!
    do k = kstart,kend
    do j = 1,m2m
    do i = 1,m1m
      u1m (i,j,k) = q1m (i,j,k) 
      u2m (i,j,k) = q2m (i,j,k)/rm(j)
      u3m (i,j,k) = q3m (i,j,k)    
      densm (i,j,k) = dens (i,j,k)
      prm (i,j,k) = pr (i,j,k)    
      uxm (i,j,k) = -sin(thetam(i))*q1m (i,j,k) &
                    +cos(thetam(i))*q2m (i,j,k)/rm(j) 
      uym (i,j,k) =  cos(thetam(i))*q1m (i,j,k) & 
                    +sin(thetam(i))*q2m (i,j,k)/rm(j)
    enddo
    enddo
    enddo
    do k = kstart,kend
    do j = 1,m2m
    do i = m1,m1
      u1m (i,j,k) = u1m (1,j,k)
      u2m (i,j,k) = u2m (1,j,k)
      u3m (i,j,k) = u3m (1,j,k)
      densm (i,j,k) = densm (1,j,k)
      prm (i,j,k) = prm (1,j,k)
      uxm (i,j,k) = uxm (1,j,k)
      uym (i,j,k) = uym (1,j,k)
      xm (i,j,k) = xm (1,j,k)
      ym (i,j,k) = ym (1,j,k)
      zm (i,j,k) = zm (1,j,k)
    enddo
    enddo
    enddo
!
!   hdf5 format
    call mpi_write_continua(m1,m2m,m3,kstart,kend,6,xm,1)
    call mpi_write_continua(m1,m2m,m3,kstart,kend,7,ym,1)
    call mpi_write_continua(m1,m2m,m3,kstart,kend,8,zm,1)
!   call mpi_write_continua(m1,m2m,m3,kstart,kend,4,densm,1)
    call mpi_write_continua(m1,m2m,m3,kstart,kend,1,u1m,1)
    call mpi_write_continua(m1,m2m,m3,kstart,kend,2,u2m,1)
    call mpi_write_continua(m1,m2m,m3,kstart,kend,3,u3m,1)
    call mpi_write_continua(m1,m2m,m3,kstart,kend,5,prm,1)
    call mpi_write_continua(m1,m2m,m3,kstart,kend,9,uxm,1)
    call mpi_write_continua(m1,m2m,m3,kstart,kend,10,uym,1)
!
!   plot3d format
    open(20,file='continua.g',form='unformatted')
    write(20) m1,m2m,m3m
    write(20) ((((xm(i,j,k),i=1,m1),j=1,m2m),k=1,m3m), &
               (((ym(i,j,k),i=1,m1),j=1,m2m),k=1,m3m), &
               (((zm(i,j,k),i=1,m1),j=1,m2m),k=1,m3m))
    close(20)
!
    open(30,file='continua.q',form='unformatted')
    write(30) m1,m2m,m3m,int(1)
    write(30) (((u3m(i,j,k),i=1,m1),j=1,m2m),k=1,m3m)
    close(30)

    print *, nrank,'/',nproc
    call MPI_FINALIZE(ierr)
    end program
! 
! 
    subroutine mpi_read_continua(n1o,n2o,n3o,ks,ke,intvar,qua)
    !use constants
    !use mpih, only: MPI_COMM_WORLD,MPI_INFO_NULL
    use hdf5  !this module contains all necessary HDF5 modules
    implicit none
    include 'mpif.h'
    integer,intent(in) :: ks,ke,n2o,n1o,n3o,intvar
    integer,parameter :: dp = 8
    real(dp),dimension(1:n1o,1:n2o,ks-1:ke+1)::qua
    integer(HID_T) :: file_id,slabspace,memspace,dset_qua,plist_id
    integer(HSIZE_T) :: dims(3)
    integer(HSIZE_T), dimension(3) :: data_count
    integer(HSSIZE_T), dimension(3) :: data_offset 
    integer :: ndims,hdf_error
    character*70 :: filnam1
    character*10 :: dsetname

!   Select file and dataset based on intvar
    select case (intvar)
      case (1)
        dsetname = trim('Vth')
        filnam1 = trim('continua_q1.h5')
      case (2)
        dsetname = trim('Vr')
        filnam1 = trim('continua_q2.h5')
      case (3)
        dsetname = trim('Vz')
        filnam1 = trim('continua_q3.h5')
      case (4)
        dsetname = trim('dens')
        filnam1 = trim('continua_dens.h5')
      case (5)
        dsetname = trim('Pr')
        filnam1 = trim('continua_pr.h5')
    end select

!   Set offsets and element counts
    ndims = 3

    dims(1)=n1o
    dims(2)=n2o
    dims(3)=n3o-1

    data_count(1) = n1o
    data_count(2) = n2o
    data_count(3) = ke-ks+1

    data_offset(1) = 0
    data_offset(2) = 0
    data_offset(3) = ks-1

    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,hdf_error)
    call h5fopen_f(filnam1,H5F_ACC_RDONLY_F,file_id,hdf_error,access_prp=plist_id)
    call h5dopen_f(file_id,dsetname,dset_qua,hdf_error)

    call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

    call h5dget_space_f(dset_qua,slabspace,hdf_error)
    call h5sselect_hyperslab_f(slabspace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
    call h5dread_f(dset_qua,H5T_NATIVE_DOUBLE,&
         qua(1:n1o,1:n2o,ks:ke),data_count,&
         hdf_error,file_space_id=slabspace,mem_space_id=memspace,xfer_prp=plist_id)

    call h5pclose_f(plist_id, hdf_error)
    call h5dclose_f(dset_qua, hdf_error)
    call h5sclose_f(memspace, hdf_error)
    call h5fclose_f(file_id, hdf_error)
    
    end subroutine mpi_read_continua
! 
! 
      subroutine mpi_write_continua(n1o,n2o,n3o,ks,ke,intvar,qua,version)
      !use constants
      use hdf5
      implicit none
      include 'mpif.h'
      integer,parameter :: dp = 8
      real(dp) time
      integer itav,itmov
      integer, intent(in) :: ks,ke,n2o,n1o,n3o,intvar,version
      real(dp), dimension(1:n1o,1:n2o,ks-1:ke+1)::qua
      integer(HID_T) :: file_id,filespace,slabspace,memspace,dset_qua,plist_id
      integer(HSIZE_T) :: dims(3)
      integer(HSIZE_T), dimension(3) :: data_count
      integer(HSSIZE_T), dimension(3) :: data_offset
      integer :: ndims,hdf_error,itime
      character*70 :: filnam1
      character*10 :: dsetname
      character*6 :: ipfi
      
!     Form the name of the file
      itime=nint(time)
      write(ipfi,82)itime
      !print *,'ipfi',ipfi
   82 format(i6.6)

!     Select file and dataset based on intvar
      select case (intvar)
        case (1)
          dsetname = trim('Vth')
          if (version.eq.1) filnam1 = 'continua_u1m.h5'
          if (version.eq.2) filnam1 =trim('continua_u1m'//ipfi//'.h5') 
        case (2)
          dsetname = trim('Vr')
          if (version.eq.1) filnam1 = 'continua_u2m.h5'
          if (version.eq.2) filnam1 =trim('continua_u2m'//ipfi//'.h5')
        case (3)
          dsetname = trim('Vz')
          if (version.eq.1) filnam1 = 'continua_u3m.h5'
          if (version.eq.2) filnam1 =trim('continua_u3m'//ipfi//'.h5')
        case (4)
          dsetname = trim('dens')
          if (version.eq.1) filnam1 = 'continua_densm.h5'
          if (version.eq.2) filnam1 =trim('continua_densm'//ipfi//'.h5')
        case (5)
          dsetname = trim('Pr')
          if (version.eq.1) filnam1 = 'continua_prm.h5'
          if (version.eq.2) filnam1 =trim('continua_prm_'//ipfi//'.h5')
        case (6)
          dsetname = trim('X')
          if (version.eq.1) filnam1 = 'continua_xm.h5'
          if (version.eq.2) filnam1 =trim('continua_xm_'//ipfi//'.h5')
        case (7)
          dsetname = trim('Y')
          if (version.eq.1) filnam1 = 'continua_ym.h5'
          if (version.eq.2) filnam1 =trim('continua_ym_'//ipfi//'.h5')
        case (8)
          dsetname = trim('Z')
          if (version.eq.1) filnam1 = 'continua_zm.h5'
          if (version.eq.2) filnam1 =trim('continua_zm_'//ipfi//'.h5')
        case (9)
          dsetname = trim('Vx')
          if (version.eq.1) filnam1 = 'continua_uxm.h5'
          if (version.eq.2) filnam1 =trim('continua_uxm_'//ipfi//'.h5')
        case (10)
          dsetname = trim('Vy')
          if (version.eq.1) filnam1 = 'continua_uym.h5'
          if (version.eq.2) filnam1 =trim('continua_uym_'//ipfi//'.h5')
      end select

!     Set offsets and element counts
      ndims = 3

      dims(1)=n1o  
      dims(2)=n2o  
      dims(3)=n3o-1

      data_count(1) = n1o
      data_count(2) = n2o
      data_count(3) = ke-ks+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = ks-1 

      call h5screate_simple_f(ndims,dims,filespace,hdf_error)

      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,hdf_error)
      call h5fcreate_f(filnam1,H5F_ACC_TRUNC_F,file_id,hdf_error,access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_qua,hdf_error)
      call h5screate_simple_f(ndims,data_count,memspace,hdf_error) 
      call h5dget_space_f(dset_qua,slabspace,hdf_error)
      call h5sselect_hyperslab_f(slabspace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error)
      call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_qua, H5T_NATIVE_DOUBLE,&
           qua(1:n1o,1:n2o,ks:ke),data_count,& 
           hdf_error,file_space_id=slabspace,mem_space_id=memspace,xfer_prp=plist_id)

      call h5pclose_f(plist_id ,hdf_error)
      call h5dclose_f(dset_qua ,hdf_error)
      call h5sclose_f(memspace ,hdf_error)
      call h5fclose_f(file_id  ,hdf_error)
      
      end subroutine mpi_write_continua


    subroutine mpi_read_continua_surf(n1o,n2o,n3o,ks,ke,jplane,intvar,qua)
    !use constants
    !use mpih, only: MPI_COMM_WORLD,MPI_INFO_NULL
    use hdf5  !this module contains all necessary HDF5 modules
    implicit none
    include 'mpif.h'
    integer,intent(in) :: ks,ke,jplane,n2o,n1o,n3o,intvar
    integer,parameter :: dp = 8
    real(dp),dimension(1:n1o,1,ks-1:ke+1)::qua
    integer(HID_T) :: file_id,slabspace,memspace,dset_qua,plist_id
    integer(HSIZE_T) :: dims(3)
    integer(HSIZE_T), dimension(3) :: data_count
    integer(HSSIZE_T), dimension(3) :: data_offset 
    integer :: ndims,hdf_error
    character*70 :: filnam1
    character*10 :: dsetname

!   Select file and dataset based on intvar
    select case (intvar)
      case (1)
        dsetname = trim('Vth')
        filnam1 = trim('continua_q1.h5')
      case (2)
        dsetname = trim('Vr')
        filnam1 = trim('continua_q2.h5')
      case (3)
        dsetname = trim('Vz')
        filnam1 = trim('continua_q3.h5')
      case (4)
        dsetname = trim('dens')
        filnam1 = trim('continua_dens.h5')
      case (5)
        dsetname = trim('Pr')
        filnam1 = trim('continua_pr.h5')
    end select

!   Set offsets and element counts
    ndims = 3

    dims(1)=n1o
    dims(2)=1
    dims(3)=n3o-1

    data_count(1) = n1o  
    data_count(2) = 1  
    data_count(3) = ke-ks+1

    data_offset(1) = 0 
    data_offset(2) = jplane-1
    data_offset(3) = ks-1

    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,hdf_error)
    call h5fopen_f(filnam1,H5F_ACC_RDONLY_F,file_id,hdf_error,access_prp=plist_id)
    call h5dopen_f(file_id,dsetname,dset_qua,hdf_error)

    call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

    call h5dget_space_f(dset_qua,slabspace,hdf_error)
    call h5sselect_hyperslab_f(slabspace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
    call h5dread_f(dset_qua,H5T_NATIVE_DOUBLE,&
         qua(1:n1o,1,ks:ke),data_count,&
         hdf_error,file_space_id=slabspace,mem_space_id=memspace,xfer_prp=plist_id)

    call h5pclose_f(plist_id, hdf_error)
    call h5dclose_f(dset_qua, hdf_error)
    call h5sclose_f(memspace, hdf_error)
    call h5fclose_f(file_id, hdf_error)
    
    end subroutine mpi_read_continua_surf


    subroutine mpi_read_continua_surf_tr(n1o,n2o,n3o,ks,ke,jplane,intvar,qua)
    !use constants
    !use mpih, only: MPI_COMM_WORLD,MPI_INFO_NULL
    use hdf5  !this module contains all necessary HDF5 modules
    implicit none
    include 'mpif.h'
    integer,intent(in) :: ks,ke,jplane,n2o,n1o,n3o,intvar
    integer,parameter :: dp = 8
    real(dp),dimension(1:n1o,1:n2o)::qua
    integer(HID_T) :: file_id,slabspace,memspace,dset_qua,plist_id
    integer(HSIZE_T) :: dims(3)
    integer(HSIZE_T), dimension(3) :: data_count
    integer(HSSIZE_T), dimension(3) :: data_offset 
    integer :: ndims,hdf_error
    character*70 :: filnam1
    character*10 :: dsetname

!   Select file and dataset based on intvar
    select case (intvar)
      case (1)
        dsetname = trim('Vth')
        filnam1 = trim('continua_q1.h5')
      case (2)
        dsetname = trim('Vr')
        filnam1 = trim('continua_q2.h5')
      case (3)
        dsetname = trim('Vz')
        filnam1 = trim('continua_q3.h5')
      case (4)
        dsetname = trim('dens')
        filnam1 = trim('continua_dens.h5')
      case (5)
        dsetname = trim('Pr')
        filnam1 = trim('continua_pr.h5')
    end select

!   Set offsets and element counts
    ndims = 3

    dims(1)=n1o
    dims(2)=n2o
    dims(3)=1

    data_count(1) = n1o  
    data_count(2) = n2o 
    data_count(3) = 1

    data_offset(1) = 0 
    data_offset(2) = 0
    data_offset(3) = 0

    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,hdf_error)
    call h5fopen_f(filnam1,H5F_ACC_RDONLY_F,file_id,hdf_error,access_prp=plist_id)
    call h5dopen_f(file_id,dsetname,dset_qua,hdf_error)

    call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

    call h5dget_space_f(dset_qua,slabspace,hdf_error)
    call h5sselect_hyperslab_f(slabspace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
    call h5dread_f(dset_qua,H5T_NATIVE_DOUBLE,&
         qua(1:n1o,1:n2o),data_count,&
         hdf_error,file_space_id=slabspace,mem_space_id=memspace,xfer_prp=plist_id)

    call h5pclose_f(plist_id, hdf_error)
    call h5dclose_f(dset_qua, hdf_error)
    call h5sclose_f(memspace, hdf_error)
    call h5fclose_f(file_id, hdf_error)
    
    end subroutine mpi_read_continua_surf_tr
