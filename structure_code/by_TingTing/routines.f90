module routines

contains


  !*******************************************
  ! this subroutine decides which command line arguments
  ! correspond to which input files, depending on parameter settings
  !******************************************
  subroutine sort_input_files( ifile_conf, traj_file )
    use global_variables
    character(*), intent(out) :: ifile_conf, traj_file
    character(5) :: temp
    integer :: n_files

    n_files = IARGC ()

    ! first file is always .gro file
    call getarg( 1, ifile_conf )

    Select Case(traj_avg)
    Case("yes")
       ! here trajectory file should be present
       call getarg( 2, traj_file )
       ! if an additional argument, compute partial structure factor for requested atom index
       if ( n_files == 3 ) then
          partial_structure_factor="yes"
          call getarg(3, temp )
          read(temp,'(I5)') partial_structure_factor_index
       else
          partial_structure_factor="no"
       end if

    Case default
       ! if an additional argument, compute partial structure factor for requested atom index
       if ( n_files == 2 ) then
          partial_structure_factor="yes"
          call getarg(2, temp )
          read(temp,'(I5)') partial_structure_factor_index
       else
          partial_structure_factor="no"
       end if

    End Select


  end subroutine sort_input_files




  !***********************************************************************
  ! this subroutine reads a .gro file and stores the atomic names and coordinates
  !***********************************************************************
  subroutine read_gro( ifn, n_atom, xyz, alist, box )
    use global_variables
    implicit none
    character(MAX_FN), intent(in) :: ifn
    character(MAX_ANAME), intent(out), dimension(MAX_N_ATOM) :: alist
    integer, intent(out) :: n_atom
    real*8, intent(out), dimension(MAX_N_ATOM,3) :: xyz
    real*8, intent(out), dimension(3,3) :: box
    integer :: ifile, i_atom, i_mole, index
    character(MAX_ANAME) :: aname, mname
    character(50)::line
    real*8 :: vel(3), tmp1, tmp2, tmp3
    integer :: ios

    ifile=99
    open( ifile, file=ifn, status='old' )
    read( ifile, * ) line
    read( ifile, '(I)' ) n_atom
    do i_atom = 1 , n_atom
       if ( i_atom > MAX_N_ATOM ) then
          stop " please increase setting of MAX_N_ATOM"
       endif
       Select Case(traj_format)
       Case("gro")
          ! here we assume trjconv was used, which makes a gro file with higher precision F9.4 coordinates
          read( ifile, '(I5,2A5,I5,3F9.4)' ), i_mole, mname, aname, index, xyz(i_atom,1), xyz(i_atom,2), xyz(i_atom,3)
       case default
          read( ifile, '(I5,2A5,I5,3F8.3,3F8.4)' ), i_mole, mname, aname, index, xyz(i_atom,1), xyz(i_atom,2), xyz(i_atom,3), vel(1), vel(2), vel(3)
       end Select
       ! convert nm to angstoms
       xyz(i_atom,:) = xyz(i_atom,:) * 10d0
       call trim_end( aname )
       alist( i_atom ) = aname
    end do

    box=0d0

    read(ifile,*) tmp1, tmp2, tmp3
    read( ifile, *, iostat=ios ) box(2,1), box(2,2), box(2,3)
    if ( ios < 0 ) then ! end of the file
       box = 0.0d0
       box(1,1) = tmp1
       box(2,2) = tmp2
       box(3,3) = tmp3
    else
       box(1,1) = tmp1
       box(1,2) = tmp2
       box(1,3) = tmp3
       read( ifile, * ) box(3,1), box(3,2), box(3,3)
    endif

    ! convert to angstroms
    box(:,:) = box(:,:) * 10d0


    close( ifile )
  end subroutine read_gro


  !*******************************************
  ! this subroutine reads the number of trajectories from a 
  ! xtc file that has been dumped, i.e. gmxdump -f file.xtc > traj
  !*******************************************
  subroutine get_n_trajectories( n_traj, n_atom, traj_file )
    use global_variables
    integer, intent(out) :: n_traj
    character(*), intent(in) :: traj_file
    integer, intent(in) :: n_atom

    character(50) :: line
    integer :: ifile=99, inputstatus, i_line, i_atom

    write(*,*) " figuring out number of snapshots from trajectory ...."
    open( ifile, file=traj_file, status='old' )

    n_traj = 0

    Select Case(traj_format)
    Case("gro")
       do
          Read(ifile,'(A)',Iostat=inputstatus) line
          If(inputstatus < 0) Exit
          n_traj = n_traj + 1
          Read(ifile,'(A)',Iostat=inputstatus) line
          do i_atom=1,n_atom
             Read(ifile,'(A)',Iostat=inputstatus) line
             If(inputstatus < 0) Exit
          enddo
          Read(ifile,'(A)',Iostat=inputstatus) line
       enddo
    Case default
       do
          Read(ifile,'(A)',Iostat=inputstatus) line
          If(inputstatus < 0) Exit
          n_traj = n_traj + 1
          do i_line=1,6
             Read(ifile,'(A)',Iostat=inputstatus) line
             If(inputstatus < 0) Exit
          enddo
          do i_atom=1,n_atom
             Read(ifile,'(A)',Iostat=inputstatus) line
             If(inputstatus < 0) Exit
          enddo
       enddo
    End select

    write(*,*) " number of snapshots : ", n_traj
    close( ifile )

  end subroutine get_n_trajectories



  !*******************************************
  ! this subroutine reads xyz coordinates for a snapshot of a trajectory  from a 
  ! xtc file that has been dumped, i.e. gmxdump -f file.xtc > traj
  !*******************************************
  subroutine read_trajectory_snapshot( ifile, xyz, n_atom, box )
    use global_variables
    integer,intent(in) :: ifile
    integer, intent(in) :: n_atom
    real*8,dimension(:,:),intent(out) :: xyz
    real*8,dimension(:,:),intent(inout) :: box
    character(20) :: junk1, junk2, junk3, junk4

    character(50) :: line
    integer :: inputstatus, i_line, i_atom, i_mole, index
    character(5) :: mname, aname

    Select Case(traj_format)
    Case("gro")
       Read(ifile,'(A)',Iostat=inputstatus) line
       if ( inputstatus < 0 ) stop "error reading trajectory file"
       Read(ifile,'(A)',Iostat=inputstatus) line
       do i_atom = 1 , n_atom
          ! here we assume trjconv was used, which makes a gro file with higher precision F9.4 coordinates
          read( ifile, '(I5,2A5,I5,3F9.4)' ), i_mole, mname, aname, index, xyz(i_atom,1), xyz(i_atom,2), xyz(i_atom,3)
          ! convert nm to angstoms
          xyz(i_atom,:) = xyz(i_atom,:) * 10d0
       enddo
       Read(ifile,'(A)',Iostat=inputstatus) line

    Case default
       Read(ifile,'(A)',Iostat=inputstatus) line
       if ( inputstatus < 0 ) stop "error reading trajectory file"
       Read(ifile,'(A)',Iostat=inputstatus) line
       Read(ifile,'(A)',Iostat=inputstatus) line
       ! get box
       do i=1,3
          read(ifile,'(A19,E11.5,A3,E11.5,A3,E11.5,A3)') junk1, box(i,1) , junk2, box(i,2), junk3, box(i,3), junk4
       enddo
       ! convert to angstrom
       box = box * 10d0
       Read(ifile,'(A)',Iostat=inputstatus) line
       do i_atom=1,n_atom
          read(ifile,'(A16,E12.5,A2,E12.5,A2,E12.5,A3)') junk1, xyz(i_atom,1) , junk2, xyz(i_atom,2), junk3, xyz(i_atom,3), junk4
          ! convert nm to angstoms
          xyz(i_atom,:) = xyz(i_atom,:) * 10d0
          if ( inputstatus < 0 ) stop "error reading trajectory file"
       enddo
    End Select

  end subroutine read_trajectory_snapshot




  subroutine create_atom_index( n_atom, alist )
    use global_variables
    character(MAX_ANAME), intent(in), dimension(MAX_N_ATOM) :: alist
    integer, intent(in) :: n_atom   

    integer :: i_atom, i_type, flag_new

    ! get first atom type
    n_atom_type = 1

    atom_index(1) = 1
    atype_name(1) = alist(1)

    do i_atom = 2, n_atom
       ! see if this a new atom type
       flag_new=1
       do i_type = 1, n_atom_type
          if ( alist( i_atom ) .eq. atype_name(i_type ) ) then
             flag_new=0
             atom_index(i_atom ) = i_type
             exit
          end if
       enddo

       ! if new atom type
       if ( flag_new == 1 ) then
          n_atom_type = n_atom_type + 1
          atype_name(n_atom_type) = alist( i_atom )
          atom_index( i_atom ) = n_atom_type
       end if
    enddo

  end subroutine create_atom_index



  !*********************************
  ! this subroutine gets the atomic form factors
  ! for every atomtype in the system
  !*********************************

  subroutine get_atomic_form_factor
    use global_variables

    character(MAX_FN) :: form_file
    character(4) :: prefix, suffix
    character(MAX_ANAME) :: aname

    integer :: i, i_type, ifile=99

    prefix="AFF_"
    suffix=".out"

    ! loop over atom types
    do i_type = 1, n_atom_type
       aname = atype_name(i_type)
       call trim_head( aname )
       form_file = prefix//trim(aname)//suffix

       write(*,*) "getting form factor for atom type ", aname , "  ....."

       open( ifile, file=form_file, status='old' )

       do i=1, max_q_form
          read(ifile,*) q_grid_form(i) , atomtype_form(i_type, i )

          ! convert from nm^-1, to angstrom^-1
          !q_grid_form(i) = q_grid_form(i) / 10d0

       enddo
       dq_form = q_grid_form(2) - q_grid_form(1)

       close(ifile)
    enddo


  end subroutine get_atomic_form_factor



  !********************************************
  ! this subroutine creates direct coordinates, scaled
  ! by input integer "K" (pme_grid), using the
  ! reciprocal lattice vectors
  ! 
  ! note, only coordinates are returned for atoms of type "i_type"
  !********************************************
  subroutine create_scaled_direct_coordinates(i_type,xyz_scale, xyz, n_atom, n_atom_kind, kk, K)
    use global_variables
    integer  , intent(in)  :: i_type
    real*8,dimension(:,:),intent(out) :: xyz_scale
    real*8,dimension(:,:),intent(in) :: xyz
    integer, intent(in) :: n_atom
    integer, intent(out) :: n_atom_kind
    real*8,dimension(:,:),intent(in) :: kk
    integer, intent(in) :: K

    integer :: i_atom,j,l, index
    real*8,parameter :: small=1D-6

    n_atom_kind=0
    do j=1,n_atom
       index = atom_index(j)
       ! if desired atom type
       if ( index == i_type ) then
          n_atom_kind = n_atom_kind + 1
          i_atom = n_atom_kind
          do l=1,3
             xyz_scale(i_atom,l)=dble(K)*dot_product(kk(l,:),xyz(j,:))
             ! if atoms are out of grid, shift them back in
             if (xyz_scale(i_atom,l)<0.) then
                xyz_scale(i_atom,l)=xyz_scale(i_atom,l)+dble(K)
             else if(xyz_scale(i_atom,l)>= dble(K)) then
                xyz_scale(i_atom,l)=xyz_scale(i_atom,l)-dble(K)
             endif
             ! make sure scaled coordinates are not numerically equal to zero, otherwise this will screw up Q grid routine
             if ( abs(xyz_scale(i_atom,l)) < small ) then
                xyz_scale(i_atom,l) = small
             end if
          enddo
       end if
    enddo

  end subroutine create_scaled_direct_coordinates



  !******************************************
  ! reciprocal lattice vector.  This is essentially the same
  ! as initialize_non_orth_transform subroutine, but we keep both
  ! in for compatibility with older code
  !******************************************
  subroutine construct_reciprocal_lattice_vector(kk,box)
    real*8,dimension(:,:),intent(out) :: kk
    real*8,dimension(:,:),intent(in) :: box

    real*8 :: a(3), b(3), c(3), ka(3), kb(3), kc(3), vol

    a(:) = box(1,:)
    b(:) = box(2,:)
    c(:) = box(3,:)

    ! calculate the volume and the reciprocal vectors (notice no 2pi)
    vol = volume( a, b, c )
    call crossproduct( a, b, kc ); kc = kc /vol 
    call crossproduct( b, c, ka ); ka = ka /vol
    call crossproduct( c, a, kb ); kb = kb /vol
    kk(1,:)=ka(:);kk(2,:)=kb(:);kk(3,:)=kc(:)

  end subroutine construct_reciprocal_lattice_vector



  subroutine crossproduct( a,b,ans )
    implicit none
    real*8, intent(in) :: a(3),b(3)
    real*8, intent(out) :: ans(3)
    ans(1) = a(2)*b(3)-a(3)*b(2)
    ans(2) = -a(1)*b(3)+a(3)*b(1)
    ans(3) = a(1)*b(2)-a(2)*b(1)
  end subroutine crossproduct

  real function volume( a, b, c )
    implicit none
    real*8, intent(in) :: a(3), b(3), c(3)
    volume = a(1) * (b(2)*c(3)-b(3)*c(2)) - a(2) * (b(1)*c(3)-b(3)*c(1)) + a(3) * (b(1)*c(2)-b(2)*c(1))
    volume = abs(volume)
  end function volume


  !*************************************************************************
  ! This subroutine moves all spaces at the beginning of a string to the end
  !*************************************************************************
  subroutine trim_end( aname )
    implicit none
    character(*), intent(inout) :: aname
    integer :: i, n, n_space, flag
    n = len( aname )
    n_space=0
    flag = 0
    do i = n, 1, -1
       if ( flag == 0 .and. aname(i:i) == ' ' ) then
          n_space = n_space + 1
       else if ( flag == 0 ) then
          flag = 1
          aname(i+n_space:i+n_space) = aname(i:i)
       else
          aname(i+n_space:i+n_space) = aname(i:i)
       end if
    end do
    do i = 1, n_space
       aname(i:i) = ' '
    end do
  end subroutine trim_end


  !*************************************************************************
  ! This subroutine moves all spaces at the beginning of a string to the end
  !*************************************************************************
  subroutine trim_head( aname )
    implicit none
    character(*), intent(inout) :: aname
    integer :: i, n, n_space, flag
    n = len( aname )
    n_space = 0
    flag = 0
    do i = 1, n
       if ( flag == 0 .and. aname(i:i) == ' ' ) then
          n_space = n_space + 1
       else if ( flag == 0 ) then
          flag = 1
          aname(i-n_space:i-n_space) = aname(i:i)
       else 
          aname(i-n_space:i-n_space) = aname(i:i)
       end if
    end do
    do i = n-n_space+1, n
       aname(i:i) = ' ' 
    end do
  end subroutine trim_head


end module routines
