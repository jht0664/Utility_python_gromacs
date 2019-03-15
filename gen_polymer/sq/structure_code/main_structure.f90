program compute_structure_factor
  use routines
  use global_variables
  use structure_factor
  use MKL_DFTI
  use pme_routines

  character(MAX_FN) :: ifile_conf, traj_file
  character(MAX_ANAME), dimension(MAX_N_ATOM) :: alist
  integer      :: n_atom
  real*8, dimension(MAX_N_ATOM,3) :: xyz
  real*8, dimension(3,3) :: box, kk
  TYPE(DFTI_DESCRIPTOR), POINTER :: dfti_desc,dfti_desc_inv

  integer :: i_atom, i_type

  call sort_input_files( ifile_conf, traj_file )
  call read_gro( ifile_conf, n_atom, xyz, alist, box )
  call create_atom_index( n_atom, alist )
  call get_atomic_form_factor
  !call construct_reciprocal_lattice_vector(kk,box)

  ! initialize bspline interpolation and FFT
  call initialize_spline_FFT(dfti_desc,dfti_desc_inv)

  if ( traj_avg == "yes" ) then
  call generate_structure_factor( n_atom, xyz, box, dfti_desc,dfti_desc_inv, traj_file )
  else
  call generate_structure_factor( n_atom, xyz, box, dfti_desc,dfti_desc_inv )
  endif


end program compute_structure_factor
