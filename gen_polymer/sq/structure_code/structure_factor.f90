module structure_factor

contains

  !*******************************
  ! if optional argument traj_file is present, then
  ! structure factor is averaged over the
  ! trajectory
  !*******************************
  subroutine generate_structure_factor( n_atom, xyz, box, dfti_desc,dfti_desc_inv, traj_file )
    use global_variables
    use routines
    use pme_routines
    use MKL_DFTI
    integer, intent(in) :: n_atom
    real*8,dimension(:,:),intent(inout) :: xyz
    real*8,dimension(3,3),intent(inout) :: box
    TYPE(DFTI_DESCRIPTOR), pointer,intent(in):: dfti_desc,dfti_desc_inv
    character(*), intent(in) , optional :: traj_file
    complex*16,dimension(:,:,:),allocatable::FQ
    real*8, dimension(:,:,:),allocatable :: SQ2
    real*8, dimension(:,:,:,:,:),allocatable :: SQ2_store
    complex*16, dimension(:,:,:,:),allocatable :: SQ_store
    real*8,dimension(:), allocatable::q_1r
    complex*16,dimension(:), allocatable::q_1d
    real*8,dimension(3,3) :: kk, kk_avg
    integer :: n, K
    real*8,dimension(:,:),allocatable :: xyz_scale
    integer :: n_atom_kind, i_type, status, n_traj, i_step, ifile=99, i_atom

    n=spline_order
    K=pme_grid
    kk_avg=0d0


    allocate( FQ(pme_grid,pme_grid,pme_grid), SQ2(pme_grid,pme_grid,pme_grid), q_1r(pme_grid**3), q_1d(pme_grid**3) )
    SQ2=0d0

    Select Case(partial_structure_factor)
    Case("no")
       allocate(  SQ_store(n_atom_type,pme_grid,pme_grid,pme_grid), SQ2_store(n_atom_type,n_atom_type,pme_grid,pme_grid,pme_grid) )
       SQ_store=0d0
       SQ2_store=0d0
    End Select


    !************* if we are averaging over a trajectory ******************
    if ( present(traj_file) ) then
       ! gen ntraj
       call get_n_trajectories( n_traj, n_atom, traj_file )
       ! open trajectory file, so as we read it in on the fly
       open( ifile, file=traj_file, status='old' )
    else
       n_traj=1
    end if


    call construct_reciprocal_lattice_vector(kk, box)
    ! create scaled coordinates
    allocate(xyz_scale(n_atom,3))

    ! now loop over trajectories

    do i_step =1, n_traj
       if ( present(traj_file) ) then
          ! read coordinates from trajectory file
          ! in case this is npt simulation, we average kvectors (box lengths) over the different snapshots
          call read_trajectory_snapshot( ifile , xyz , n_atom, box )
          call construct_reciprocal_lattice_vector(kk, box)
          kk_avg = kk_avg + kk
       endif

       Select Case(partial_structure_factor)
       Case("no")
          !****************** computing full structure factor
          ! loop over  atom types for structure factor
          do i_type = 1 , n_atom_type
             ! note this only creates coordinates for atomtype "i_type"
             call create_scaled_direct_coordinates(i_type, xyz_scale, xyz, n_atom, n_atom_kind, kk, K)
             call grid_Q(Q_grid,xyz_scale,n_atom_kind,K,n)
             q_1r=RESHAPE(Q_grid, (/K**3/) )
             q_1d=cmplx(q_1r,0.,16)
             status=DftiComputeForward(dfti_desc, q_1d)
             FQ=RESHAPE(q_1d, (/K,K,K/) )
             ! structure factor = B * FQ
             SQ = FQ*B
             SQ_store(i_type,:,:,:) = SQ(:,:,:)
          enddo

          ! now create all the cross SQ2 structure factors for this snapshot
          call combine_partial_structure_factors( SQ2_store , SQ_store )

       Case("yes")
          !********************* computing partial structure factor
          i_type = partial_structure_factor_index
          call create_scaled_direct_coordinates(i_type, xyz_scale, xyz, n_atom, n_atom_kind, kk, K)
          call grid_Q(Q_grid,xyz_scale,n_atom_kind,K,n)
          q_1r=RESHAPE(Q_grid, (/K**3/) )
          q_1d=cmplx(q_1r,0.,16)
          status=DftiComputeForward(dfti_desc, q_1d)
          FQ=RESHAPE(q_1d, (/K,K,K/) )
          ! structure factor = B * FQ
          SQ = FQ*B
          SQ2 = SQ2 + dble(SQ)**2+aimag(SQ)**2
       End Select

    enddo


    ! now average
    Select Case(partial_structure_factor)
    Case("no")
       SQ2_store = SQ2_store / n_traj
       kk_avg = kk_avg / n_traj
       ! now add partial structure factors
       call add_partial_structure_factors( SQ2 , SQ2_store, kk_avg )
    Case("yes")
       SQ2 = SQ2 / n_traj
       kk_avg = kk_avg / n_traj
    end Select

    ! now print, note because we took the FT in reduced coordinates, we need to convert to 
    ! the physical wavevectors
    write(*,*) " k vec ,  SQ^2 "
    call print_SQ( SQ2 , kk_avg, pme_max_print )


    deallocate( FQ, SQ2, q_1r, q_1d )

    Select Case(partial_structure_factor)
    Case("no")
       deallocate( SQ_store, SQ2_store)
    End Select

    if ( present(traj_file) ) then
       close( ifile )
    endif

  end subroutine generate_structure_factor



  subroutine add_partial_structure_factors( SQ2 , SQ2_store,kk )
    use global_variables
    real*8,dimension(:,:,:), intent(out) :: SQ2
    real*8,dimension(:,:,:,:,:), intent(in) :: SQ2_store
    real*8,dimension(3,3) , intent(in) :: kk

    integer :: i_type, j_type, i_k, j_k, l_k
    real*8  :: fi, fj, k_vec(3), k_mag

    SQ2=0d0

    ! here we add only those structure factors that will be printed, to save time
    do l_k = 1, pme_max_print
       do j_k = 1, pme_max_print
          do i_k = 1, pme_max_print

             ! convert wavevector, note the reciprocal lattice vectors kk don't have the 2*pi factor
             k_vec(:) = 2 * pi * ( dble(i_k-1) * kk(1,:) +  dble(j_k-1) * kk(2,:) +  dble(l_k-1) * kk(3,:)  )
             k_mag = sqrt(dot_product(k_vec,k_vec))

             do i_type=1,n_atom_type
                do j_type = i_type, n_atom_type

                   fi = get_form_fac( i_type, k_mag )
                   fj = get_form_fac( j_type, k_mag )

                   SQ2(i_k,j_k,l_k) = SQ2(i_k,j_k,l_k) + fi * fj * SQ2_store(i_type,j_type,i_k,j_k,l_k)
		enddo
             enddo
          enddo
       enddo
    enddo

  end subroutine add_partial_structure_factors




  subroutine  combine_partial_structure_factors( SQ2 , SQ_store )
    use global_variables
    real*8,dimension(:,:,:,:,:),intent(inout) :: SQ2
    complex*16,dimension(:,:,:,:),intent(in) :: SQ_store

    integer :: i_type, j_type

    do i_type=1,n_atom_type
       do j_type = i_type, n_atom_type

          if ( i_type == j_type ) then
             SQ2(i_type,j_type,:,:,:) = SQ2(i_type,j_type,:,:,:) + dble(SQ_store(i_type,:,:,:))**2+aimag(SQ_store(i_type,:,:,:))**2
          else
             ! here we have unlike types, so we want SQi(c.c.) * SQj + SQi * SQj(c.c.) , where
             ! (c.c.) is complex conjugate, which equals 2 * real ( SQi * SQj(c.c.) )
             SQ2(i_type,j_type,:,:,:) = SQ2(i_type,j_type,:,:,:) + 2d0 * dble( SQ_store(i_type,:,:,:) * conjg((SQ_store(j_type,:,:,:))) )
          end if
       enddo
    enddo

  end subroutine combine_partial_structure_factors



  subroutine print_SQ( SQ2 , kk, K )
    real*8,dimension(:,:,:),intent(in) :: SQ2
    real*8,dimension(3,3),intent(in) :: kk
    integer, intent(in) :: K

    integer :: i, j , l , n
    real*8, dimension(3) :: k_vec
    real*8, parameter :: pi = 3.141592654

    do i=1, K
       do j=1,K
          do l=1,K
             ! convert wavevector, note the reciprocal lattice vectors kk don't have the 2*pi factor
             k_vec(:) = 2 * pi * ( dble(i-1) * kk(1,:) +  dble(j-1) * kk(2,:) +  dble(l-1) * kk(3,:)  )
             write(*,'(3F14.6, F20.6)') k_vec, SQ2(i,j,l)
          enddo
       enddo
    enddo

  end subroutine print_SQ


  ! *************************
  ! if q_mag is greater than the longest qvec for which we
  ! have the form factor, then return 0
  !*************************
  real*8 function get_form_fac( i_index, q_mag )
    use global_variables
    integer, intent(in) :: i_index
    real*8, intent(in) :: q_mag

    integer :: i, index, flag

    !    flag=0
    !    do i=1, max_q_form-1
    !       if ( ( q_grid_form(i) < q_mag ) .and. ( q_mag < q_grid_form(i+1)  ) ) then
    !          index = i
    !          flag=1
    !          exit
    !       endif
    !    enddo

    !    if ( flag == 0 ) then
    !       get_form_fac=0d0
!!$       write(*,*) "couldn't find atomic form factor for q value ", q_mag
!!$       stop
    !   else
    index = ceiling( q_mag / dq_form )
    if ( index < max_q_form ) then
       get_form_fac = atomtype_form( i_index, index )
    else
       get_form_fac = 0d0
    endif
    !   endif

  end function get_form_fac




end module structure_factor
