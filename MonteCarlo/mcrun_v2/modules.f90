MODULE pos
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: x, y, z
  CHARACTER(LEN=4), DIMENSION(:), ALLOCATABLE, SAVE :: typs
  double precision, save :: pos_scaling, pos_scaling2, pos_scaling3
  DOUBLE PRECISION, DIMENSION(3), save :: box
  double precision, dimension(3), save :: coord_try
  CHARACTER(LEN=4), ALLOCATABLE, SAVE :: typs_try
END MODULE pos

MODULE time
  double precision :: time1, time2, time3, time4, time11, time22, time33
end MODULE time

MODULE inp
  INTEGER(kind=8), save :: nncon, ncon, nskip
  integer, save :: openmp_thread_num
  CHARACTER(LEN=3), save :: iconfig, ipres, igr, iex
  logical, save :: check_overlap_log
END MODULE inp

MODULE traj
  INTEGER, save :: traj_nstep ! # step skip to save trajectory
END MODULE traj

MODULE trans
  DOUBLE PRECISION, SAVE :: trans_prob, dlr_a, dlr_b
END MODULE trans

MODULE calc_pres
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, save :: print_press, scale_length2
  INTEGER, DIMENSION(:), ALLOCATABLE, save :: pres_noverlap
  DOUBLE PRECISION, SAVE :: ratio_dv_v
  INTEGER, SAVE :: n_dv
END MODULE calc_pres

MODULE calc_ex
  INTEGER, SAVE :: ex_ntry
  integer, save :: ex_noverlap
  character(LEN=8), save :: ex_solvent, ex_solute
END MODULE calc_ex

MODULE ints
  INTEGER, save :: nptot, nmol_a, nmon_a, nmol_b, nmon_b
  double precision, save :: a_frac
END MODULE ints

MODULE cellmap
  INTEGER, DIMENSION(:), ALLOCATABLE, save :: lead, map, list, jnear
  INTEGER, save :: mcx, mcy, mcz, niter
  DOUBLE PRECISION, save :: cellx, celly, cellz
END MODULE cellmap

MODULE movetype
  logical :: ensemble_pres, ensemble_temp, ensemble_exch
  integer, save :: nmovetypes
  CHARACTER(LEN=5), DIMENSION(:), ALLOCATABLE, save :: movetype_name
  INTEGER(kind=8), dimension(:), allocatable, save :: movetype_ntry, movetype_nsuccess
  INTEGER, dimension(:), allocatable, save :: movetype_i_try, movetype_i_success
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, save :: movetype_prob
END MODULE

MODULE coupling_exch
  character(len=4), dimension(:), allocatable, save :: exch_tcomp
  integer, save :: exch_ncomp
  double precision, save :: xi_val, log_xi
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: exch_typs
  logical, save :: xi_val_opt
END MODULE

MODULE coupling_pres
  integer, save :: semiiso
  DOUBLE precision, save :: press_val, dvx, tinv
  DOUBLE precision, save :: delta, expd
END MODULE

MODULE sigmas
  DOUBLE PRECISION :: sigma_b, sigma_a, sigma_ab
  SAVE sigma_b, sigma_a, sigma_ab
END MODULE sigmas

MODULE rdf
  DOUBLE PRECISION :: bsz ! bin size
  INTEGER :: rdf_nstep, nbin
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gr12
  SAVE nbin,bsz,gr12
END MODULE rdf

MODULE filenames
  CHARACTER(LEN=100),save :: filename_pres, filename_dens, filename_ex
  CHARACTER(LEN=100),save :: filename_traj_gro, filename_nmol
END MODULE filenames

!MODULE minheap
!  type short_dist_in_cell
!    integer :: index_cell
!    double precision :: distance
!  end type short_dist_in_cell
!  type(short_dist_in_cell), dimension(:), allocatable, save :: minheap_dist
!  integer, save :: minheap_size
!  integer, dimension(:), ALLOCATABLE, save :: minheap_address
!  double precision, save :: dist_criteria
!END MODULE
!