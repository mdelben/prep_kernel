#ifdef ONE_API
include "mkl_omp_offload.f90"
#endif
        module aux_mod
         use mpi
         use omp_lib
#ifdef NVHPC
         use cublas
         use cudafor
#elif defined(HIP_API)
         use ISO_C_BINDING
         use hipfort
         use hipfort_hipblas
#elif defined(ONE_API)
         use onemkl_blas_omp_offload_lp64
#endif

         implicit none

         public :: peinfo, init_mpi, clean_mpi, scalapack, polarizability
         public :: alloc_data, dealloc_data
         public :: accel_enter_data_map_to_r6, accel_exit_data_map_delete_r6

         integer, parameter, public :: OMP_TARGET_ALGO = 1
         integer, parameter, public :: CPU_ALGO = 0
         double complex, parameter, public :: ZERO = (0.0D+00,0.0D+00)

         type scalapack
           integer :: nprow  !< the number of processors in a row of your processor grid
           integer :: npcol  !< the number of processors in a column of your processor grid
           integer :: nbl    !< the linear dimension of a block of a distributed matrix
           integer :: myprow !< the processor`s row coordinate in your processor grid
           integer :: mypcol !< the processor`s column coordinate in your processor grid
           integer :: npr    !< the number of rows of the matrix the processor owns
           integer :: npc    !< the number of columns of the matrix the processor owns
           integer :: icntxt !< BLACS context; see BLACS documentation 
           integer, pointer :: npcd(:)       !< global list of the number of cols of the matrix owned by all processors
           integer, pointer :: nprd(:)       !< global list of the number of rows of the matrix owned by all processors
           integer, pointer :: isrtxrow(:)   !< isrtxrow/isrtxcol give the sorted index of the gvector in a given block  
           integer, pointer :: isrtxcol(:)   !! owned by a processor in terms of the whole list of gvectors 
           integer, pointer :: imycol(:)     !< imyrow/imycol give the row/column index of a g-vector owned by a given
           integer, pointer :: imyrow(:)     !! processor in the whole matrix
           integer, pointer :: imycolinv(:)  !< inverse of imycol
           integer, pointer :: imyrowinv(:)  !! inverse of imyrow
           integer, pointer :: imycold(:,:)  !< imycold/imyrowd are global lists of the row/column index of g-vectors 
           integer, pointer :: imyrowd(:,:)  !! owned by all the processors in the whole matrix
         end type scalapack

         type polarizability
           double complex, pointer :: gme(:,:,:,:,:,:)
           double complex :: negfact
           integer :: nfreq_group
         end type polarizability

         type :: peinf
           integer :: inode, npe
           integer :: comm_global
           integer :: err
           integer :: num_threads
           integer :: num_devices, my_active_device
           integer :: nv_block, ncownactual
           integer :: ntot, nrow_max, ncol_max
           integer :: s1, s2
           integer :: algo
           logical :: full_offload
         end type peinf

         type(peinf) :: peinfo

        contains

#if defined(HIP_API)
         function hipblas_t2op(tran) result(ret)
             use hipfort_hipblas_enums
             implicit none
             integer(kind(HIPBLAS_OP_N)) :: ret
             character                   :: tran
             !
             select case (TRIM(tran))
               case ('N','n')
                 ret = HIPBLAS_OP_N
               case ('T','t')
                 ret = HIPBLAS_OP_T
               case ('C','c')
                 ret = HIPBLAS_OP_C
               case default
                 STOP "HIP BLAS operator not mapper "
             end select
         end function
#endif

         subroutine init_mpi()
           implicit none
           character(len=32) :: arg
           integer :: ii, algo, full_offload

           call MPI_INIT(peinfo%err)
           peinfo%comm_global = MPI_COMM_WORLD
           call MPI_COMM_RANK(peinfo%comm_global, peinfo%inode, peinfo%err)
           call MPI_COMM_SIZE(peinfo%comm_global, peinfo%npe, peinfo%err)

           ! threads
           peinfo%num_threads = 1
           !$OMP parallel
           peinfo%num_threads = OMP_GET_NUM_THREADS()
           !$OMP end parallel

           ! initialize GPUs (OpenMP target)
           peinfo%num_devices      = omp_get_num_devices()
           peinfo%my_active_device = MOD(peinfo%inode, peinfo%num_devices)
           call omp_set_default_device( peinfo%my_active_device )

           if (peinfo%inode == 0) then
             write(*,"(A,I4)") "MPI Tasks:         ", peinfo%npe
             write(*,"(A,I4)") "OMP threads:       ", peinfo%num_threads
             write(*,"(A,I4)") "Num active Devices:", peinfo%num_devices
             write(*,*)
           end if
           do ii = 0, peinfo%npe-1
             if ( ii == peinfo%inode ) then
               write(*,"(2(A,I4))") "My rank:", peinfo%inode,"   My active device: ", peinfo%my_active_device
             end if
             call MPI_BARRIER(peinfo%comm_global, peinfo%err)
           end do
       
           if (peinfo%inode == 0) write(*,*)

           peinfo%algo = 0
           peinfo%full_offload = .false.

           call getarg(1, arg)
           read(arg,*) algo
           call getarg(2, arg)
           read(arg,*) full_offload

           if(algo == 1)         peinfo%algo = OMP_TARGET_ALGO
           if(full_offload == 1) peinfo%full_offload = .true.

           if ( algo /= 0 .and. algo /= 1 ) then
             stop 'First argumet shoul be 0 for CPU or 1 for GPU run'
           end if
           if ( full_offload /= 0 .and. full_offload /= 1 ) then
             stop 'Second argumet shoul be 0 for submat or 1 for fulloffload GPU run '
           end if
           if ( peinfo%full_offload .and. peinfo%algo/=OMP_TARGET_ALGO ) then
             stop 'Full Offload only works with GPU algo.'
           end if

           if (peinfo%inode == 0) then
             if ( peinfo%algo == OMP_TARGET_ALGO ) then
               write(*,*) "GPU Run"
             else
               write(*,*) "CPU Run"
             end if
             if ( peinfo%full_offload ) then
               write(*,*) "Full Offload Algorithm"
             else
               write(*,*) "Submat Algorithm"
             end if
             write(*,*)
           end if

         end subroutine init_mpi

         subroutine clean_mpi()
           implicit none

           call MPI_BARRIER(peinfo%comm_global, peinfo%err)
           call MPI_Finalize(peinfo%err)

         end subroutine clean_mpi

         subroutine alloc_data(scal,nrk,nst,gmetempr,gmetempc,chilocal,pol,indt,pht,ipe,ispin)
           implicit none
           type (scalapack), intent(out) :: scal
           integer, intent(out) :: nrk
           integer, allocatable, intent(out) :: nst(:) !< (nrk)
           double complex, allocatable, intent(out) :: gmetempr(:,:),gmetempc(:,:)
           double complex, allocatable, intent(out) :: chilocal(:,:)
           type (polarizability), intent(out) :: pol
           integer, allocatable, dimension(:,:,:), intent(out) :: indt
           double complex, allocatable, dimension(:,:,:), intent(out) :: pht
           integer, intent(out) :: ipe
           integer, intent(out) :: ispin

           integer :: d1, d2, d3, d4, d5, d6

           open(6000,file="prep_kern.dat",form='unformatted')
           
           read(6000) peinfo%ntot, peinfo%nrow_max, peinfo%ncol_max, pol%negfact
           read(6000) d1, d2
           allocate(scal%nprd(d1))
           allocate(scal%npcd(d2))
           read(6000) scal%nprd
           read(6000) scal%npcd
           read(6000) nrk,d1,ipe,ispin
           allocate(nst(d1))
           read(6000) nst
           read(6000) d1, d2
           allocate(scal%imyrowd(d1,d2))
           read(6000) scal%imyrowd
           read(6000) d1, d2
           allocate(scal%imycold(d1, d2))
           read(6000) scal%imycold
           read(6000) d1, d2, d3
           allocate(indt(d1, d2, d3))
           read(6000) indt
           read(6000) d1, d2, d3
           allocate(pht(d1, d2, d3))
           read(6000) pht
           read(6000) peinfo%nv_block, peinfo%ncownactual
           read(6000) d1, d2
           allocate(gmetempr(d1,d2))
           read(6000) d1, d2
           allocate(gmetempc(d1,d2))
           read(6000) d1, d2
           allocate(chilocal(d1,d2))
           peinfo%s1 = d1
           peinfo%s2 = d2
           read(6000) d1, d2, d3, d4, d5, d6
           allocate(pol%gme(d1, d2, d3, d4, d5, d6))
           read(6000) pol%gme

           close(6000)

           pol%nfreq_group = d6
           gmetempr = ZERO
           gmetempc = ZERO
           chilocal = ZERO

         end subroutine alloc_data

         subroutine dealloc_data(scal,nst,gmetempr,gmetempc,chilocal,pol,indt,pht)
           implicit none
           type (scalapack) :: scal
           integer, allocatable :: nst(:) !< (nrk)
           double complex, allocatable :: gmetempr(:,:),gmetempc(:,:)
           double complex, allocatable :: chilocal(:,:)
           type (polarizability) :: pol
           integer, allocatable, dimension(:,:,:)        :: indt
           double complex, allocatable, dimension(:,:,:) :: pht

           deallocate(scal%nprd)
           deallocate(scal%npcd)
           deallocate(nst)
           deallocate(scal%imyrowd)
           deallocate(scal%imycold)
           deallocate(indt)
           deallocate(pht)
           deallocate(gmetempr)
           deallocate(gmetempc)
           deallocate(chilocal)
           deallocate(pol%gme)

         end subroutine dealloc_data

         subroutine accel_zgemm(transa, transb, &
                            m, n, k, &
                            alpha, &
                            a, lda, &
                            b, ldb, &
                            beta, &
                            c, ldc, &
                            algo)
           implicit none
      
           character :: transa, transb
           integer :: m, n, k, lda, ldb, ldc
           double complex :: alpha, beta
           double complex :: a(:,:), b(:,:), c(:,:)
           integer :: algo, ierr
#ifdef NVHPC
           type(cublasHandle) :: accel_blas_handle
#elif defined(HIP_API)
           type(c_ptr)        :: accel_blas_handle = c_null_ptr
#endif

           if (algo == OMP_TARGET_ALGO) then
             if (peinfo%inode == 0) write(*,*) 'GPU ZGEMM'
#ifdef NVHPC
               ierr = cublasCreate(accel_blas_handle)
               !$omp target data use_device_ptr(a, b, c)
               call cublasZgemm(transa, transb, &
                                m, n, k, &
                                alpha, a, lda, b, ldb, &
                                beta, c, ldc)
               !$omp end target data
               ierr = cublasDestroy(accel_blas_handle)
#elif defined(HIP_API)
               ierr = hipblasCreate(accel_blas_handle)
               !$omp target data use_device_ptr(a, b, c)
               ierr = hipblasZgemm(accel_blas_handle, &
                                   hipblas_t2op(transa), hipblas_t2op(transb), &
                                   m, n, k, &
                                   alpha, c_loc(a), lda, c_loc(b), ldb, &
                                   beta, c_loc(c), ldc)
               !$omp end target data
               ierr = hipblasDestroy(accel_blas_handle)
#elif ONE_API
             !$omp target data use_device_ptr(a, b, c)
             !$omp target variant dispatch
             call zgemm(transa, transb, &
                        m, n, k, &
                        alpha, &
                        a, lda, &
                        b, ldb, &
                        beta, &
                        c, ldc)
             !$omp end target variant dispatch
             !$omp end target data
#endif
           else if (algo == CPU_ALGO) then
             if (peinfo%inode == 0) write(*,*) 'CPU ZGEMM'
             call zgemm(transa, transb, &
                        m, n, k, &
                        alpha, &
                        a, lda, &
                        b, ldb, &
                        beta, &
                        c, ldc)
           else
             STOP "Invalid algo for zgemm specified"
           end if
         end subroutine accel_zgemm

         subroutine accel_enter_data_map_to_r6(array)
           implicit none
           double complex, intent(inout) :: array(:, :, :, :, :, :)
           !$omp target enter data map(to: array)
         end subroutine accel_enter_data_map_to_r6

         subroutine accel_exit_data_map_delete_r6(array)
           implicit none
           double complex, intent(inout) :: array(:, :, :, :, :, :)
           !$omp target exit data map(delete: array)
         end subroutine accel_exit_data_map_delete_r6

        end module aux_mod
