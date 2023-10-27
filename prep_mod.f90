  module prep_mod
  use aux_mod
  implicit none

  public :: chi_sum_prep, chi_sum_mult

  contains
 
  subroutine chi_sum_mult(gmetempr,gmetempc,chilocal,pol)
    double complex, dimension(:,:), intent(inout) :: gmetempr(:,:),gmetempc(:,:)
    double complex, dimension(:,:), intent(inout) :: chilocal(:,:)
    type (polarizability), intent(inout) :: pol

    call accel_zgemm('n', 't', &
                     peinfo%nrow_max, peinfo%ncol_max, peinfo%ntot, &
                     pol%negfact, &
                     gmetempr, peinfo%nrow_max, &
                     gmetempc, peinfo%ncol_max, &
                     ZERO, &
                     chilocal, peinfo%nrow_max, &
                     peinfo%algo)

    if ( peinfo%algo == OMP_TARGET_ALGO ) then
      !$omp target update from(chilocal)
    end if

  end subroutine chi_sum_mult

  subroutine chi_sum_prep(scal,nrk,nst,gmetempr,gmetempc,chilocal,polgme,pol,indt,pht,ipe,ispin)
    type (scalapack), intent(in) :: scal
    integer, intent(in) :: nrk
    integer, intent(in) :: nst(:) !< (nrk)
    double complex, dimension(:,:), intent(inout) :: gmetempr(:,:),gmetempc(:,:)
    double complex, dimension(:,:), intent(inout) :: chilocal(:,:)
    double complex, intent(in) :: polgme(:,:,:,:,:,:)
    type (polarizability), intent(inout) :: pol
    integer, dimension(:,:,:), intent(in) :: indt
    double complex, dimension(:,:,:), intent(in) :: pht
    integer, intent(in) :: ipe
    integer, intent(in) :: ispin

    integer, allocatable :: tmprowindex(:),tmpcolindex(:)
    double complex, allocatable :: tmprowph(:),tmpcolph(:)
    integer :: irk, iv, j, it, icurr, itot, mytot
    logical :: do_write

    ALLOCATE(tmprowindex (scal%nprd(ipe)))
    ALLOCATE(tmpcolindex (scal%npcd(ipe)))
    ALLOCATE(tmprowph (scal%nprd(ipe)))
    ALLOCATE(tmpcolph (scal%npcd(ipe)))

    if ( peinfo%algo == OMP_TARGET_ALGO .and. peinfo%full_offload ) then
      !$omp target enter data map(alloc:tmprowindex,tmpcolindex,tmprowph,tmpcolph)
    end if

    do_write = .true.
    itot = 0
    do irk = 1, nrk
      do it = 1, nst(irk)

        if ( peinfo%algo == OMP_TARGET_ALGO .and. peinfo%full_offload ) then

          if ( peinfo%inode == 0 .and. do_write ) then
            write(*,*) 'Full OffLoad Prep'
            do_write = .false.
          end if

          !$omp target teams loop default(shared) private(icurr)
          do icurr=1,scal%nprd(ipe)
            tmprowindex(icurr) = indt(scal%imyrowd(icurr,ipe),it,irk)
            tmprowph(icurr) = pht(scal%imyrowd(icurr,ipe),it,irk)
          enddo
          !$omp end target teams loop

          !$omp target teams loop default(shared) private(icurr)
          do icurr=1,scal%npcd(ipe)
            tmpcolindex(icurr) = indt(scal%imycold(icurr,ipe),it,irk)
            tmpcolph(icurr) = pht(scal%imycold(icurr,ipe),it,irk)
          enddo
          !$omp end target teams loop

          !$omp target teams loop default(shared) &
          !$omp& private(iv, j, mytot, icurr) &
          !$omp& collapse(2)
          do iv = 1,peinfo%nv_block
            do j = 1, peinfo%ncownactual

              mytot = itot + (iv-1)*peinfo%ncownactual + j

              gmetempr(:,mytot) = ZERO
              gmetempc(:,mytot) = ZERO

              do icurr=1,scal%nprd(ipe)
                gmetempr(icurr,mytot)=polgme( &
                  tmprowindex(icurr),j,iv, &
                  ispin,irk,pol%nfreq_group)* &
                  tmprowph(icurr)
              enddo

              do icurr=1,scal%npcd(ipe)
                gmetempc(icurr,mytot)= &
                  CONJG(polgme(tmpcolindex(icurr),j,iv,ispin,irk,pol%nfreq_group)*tmpcolph(icurr))
              enddo

            enddo ! j
          enddo ! iv
          !$omp end target teams loop

        else

          if ( peinfo%inode == 0 .and. do_write ) then
            write(*,*) 'Submat (CPU) Prep'
            do_write = .false.
          end if

          do icurr=1,scal%nprd(ipe)
            tmprowindex(icurr) = indt(scal%imyrowd(icurr,ipe),it,irk)
            tmprowph(icurr) = pht(scal%imyrowd(icurr,ipe),it,irk)
          enddo
          do icurr=1,scal%npcd(ipe)
            tmpcolindex(icurr) = indt(scal%imycold(icurr,ipe),it,irk)
            tmpcolph(icurr) = pht(scal%imycold(icurr,ipe),it,irk)
          enddo

          !$OMP PARALLEL DO collapse(2) private (mytot,iv,j,icurr)
          do iv = 1,peinfo%nv_block
            do j = 1, peinfo%ncownactual
              mytot = itot + (iv-1)*peinfo%ncownactual + j

              gmetempr(:,mytot) = ZERO
              gmetempc(:,mytot) = ZERO

              do icurr=1,scal%nprd(ipe)
                gmetempr(icurr,mytot)=polgme( &
                  tmprowindex(icurr),j,iv, &
                  ispin,irk,pol%nfreq_group)* &
                  tmprowph(icurr)
              enddo

              do icurr=1,scal%npcd(ipe)
                gmetempc(icurr,mytot)= &
                  CONJG(polgme(tmpcolindex(icurr),j,iv,ispin,irk,pol%nfreq_group)*tmpcolph(icurr))
              enddo ! icurr
            enddo ! j
          enddo ! iv
          !$OMP END PARALLEL DO

        end if

      end do
    end do

    if ( peinfo%algo == OMP_TARGET_ALGO ) then
      if ( peinfo%full_offload ) then
        !$omp target exit data map(delete:tmprowindex,tmpcolindex,tmprowph,tmpcolph)
      else
        !$omp target update to(gmetempr, gmetempc, chilocal)
      end if
    end if

    DEALLOCATE(tmprowindex)
    DEALLOCATE(tmpcolindex)
    DEALLOCATE(tmprowph)
    DEALLOCATE(tmpcolph)

  end subroutine chi_sum_prep

  end module prep_mod
