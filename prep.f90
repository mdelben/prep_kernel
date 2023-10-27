! # module load oneapi/eng-compiler/2023.10.15.002  e4s/23.05/2023.05.15.006-2 hdf5/1.14.2-oneapi-mpich mpich/52.2
! mpif90 -fc=ifx -free  -O2 -g -traceback -check shape -fp-model precise -no-ipo -align array64byte -fiopenmp -fopenmp-targets=spir64 -qmkl -lmkl_sycl -lsycl -lOpenCL  aux_mod.f90  prep_mod.f90  prep.f90 -o prep.x
 
    program prep
     use aux_mod
     use prep_mod
     implicit none
     type (scalapack) :: scal
     integer              :: nrk
     integer, allocatable :: nst(:) !< (nrk)
     double complex, allocatable :: gmetempr(:,:),gmetempc(:,:)
     double complex, allocatable :: chilocal(:,:)
     type (polarizability) :: pol
     integer, dimension(:,:,:), allocatable :: indt
     double complex, dimension(:,:,:), allocatable :: pht
     integer :: ipe
     integer :: ispin

     call init_mpi()

     if(peinfo%inode == 0) write(*,*) 'Read in data'
     call alloc_data(scal,nrk,nst,gmetempr,gmetempc,chilocal,pol,indt,pht,ipe,ispin)
     if(peinfo%inode == 0) write(*,*) 'Done!'
     if(peinfo%inode == 0) write(*,*) 

     if (peinfo%algo == OMP_TARGET_ALGO ) then
       if(peinfo%inode == 0) write(*,*) 'Copy data to GPU'
       !$omp target enter data map (to: chilocal, gmetempr, gmetempc)
       if ( peinfo%full_offload ) then
#ifdef ONE_API
         call accel_enter_data_map_to_r6(pol%gme)
#endif
#ifdef NVHPC
         !$omp target enter data map (to: pol)
         !$omp target enter data map (to: pol%gme)
         !$omp target enter data map (to: scal)
#endif
         !$omp target enter data map (to: scal%nprd, scal%npcd, scal%imyrowd, scal%imycold)
         !$omp target enter data map (to: indt)
         !$omp target enter data map (to: pht)
         !$omp target enter data map (to: peinfo)
       end if
       if(peinfo%inode == 0) write(*,*) 'Done!'
       if(peinfo%inode == 0) write(*,*) 
     end if

     if(peinfo%inode == 0) write(*,*) 'Prep data'
     call chi_sum_prep(scal,nrk,nst,gmetempr,gmetempc,chilocal,pol%gme,pol,indt,pht,ipe,ispin)
     if(peinfo%inode == 0) write(*,*) 'Done!'
     if(peinfo%inode == 0) write(*,*) 

     if(peinfo%inode == 0) write(*,*) 'Mult data'
     call chi_sum_mult(gmetempr,gmetempc,chilocal,pol)
     if(peinfo%inode == 0) then
       write(*,*) 'Done!'
       write(*,*) chilocal(1,1)
       write(*,*) chilocal(peinfo%s1,peinfo%s2)
       write(*,*) 
     end if

     if(peinfo%inode == 0) write(*,*) 'Cleanup'
     if (peinfo%algo == OMP_TARGET_ALGO ) then
       !$omp target exit data map(delete: chilocal, gmetempr, gmetempc)
       if ( peinfo%full_offload ) then
#ifdef ONE_API
         call accel_exit_data_map_delete_r6(pol%gme)
#endif
         !$omp target exit data map(delete: scal%nprd, scal%npcd, scal%imyrowd, scal%imycold)
#ifdef NVHPC
         !$omp target exit data map(delete: scal)
         !$omp target exit data map(delete: pol%gme)
         !$omp target exit data map(delete: pol)
#endif
         !$omp target exit data map(delete: indt)
         !$omp target exit data map(delete: pht)
         !$omp target exit data map(delete: peinfo)
       end if
     end if
     call dealloc_data(scal,nst,gmetempr,gmetempc,chilocal,pol,indt,pht)
     if(peinfo%inode == 0) write(*,*) 'Done!'
     if(peinfo%inode == 0) write(*,*) 

     call clean_mpi()

    end program prep

