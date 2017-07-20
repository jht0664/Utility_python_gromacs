MODULE pos
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: x, y, z
  CHARACTER(LEN=4), DIMENSION(:), ALLOCATABLE, SAVE :: typs
  DOUBLE PRECISION, DIMENSION(3), save :: box
END MODULE pos

MODULE inp
  INTEGER, save :: nncon, ncon, nskip, openmp_thread_num
  CHARACTER(LEN=3), save :: iconfig, ipres, igr, iex
END MODULE inp

MODULE traj
  INTEGER, save :: traj_nstep ! # step skip to save trajectory
END MODULE traj

MODULE trans
  DOUBLE PRECISION, SAVE :: trans_prob, dlr_a, dlr_b
END MODULE trans

MODULE calc_pres
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, save :: print_press, scale_length
  INTEGER, DIMENSION(:), ALLOCATABLE, save :: pres_noverlap
  DOUBLE PRECISION, SAVE :: ratio_dv_v
  INTEGER, SAVE :: n_dv
END MODULE calc_pres

MODULE calc_ex
  INTEGER, SAVE :: ex_nstep, ex_ntry
  integer, save :: ex_noverlap
  character(LEN=4), save :: ex_solvent, ex_solute
END MODULE calc_ex

MODULE ints
  INTEGER, save :: nptot, nmol_a, nmon_a, nmol_b, nmon_b
END MODULE ints

MODULE cellmap
  INTEGER, DIMENSION(:), ALLOCATABLE :: lead, map, list, jnear
  INTEGER :: mcx, mcy, mcz
  DOUBLE PRECISION, save :: cellx, celly, cellz
  SAVE lead, map, list, mcx, mcy, mcz
END MODULE cellmap

MODULE try
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, save :: xitr, yitr, zitr
  CHARACTER(LEN=4), DIMENSION(:), ALLOCATABLE, save :: titr
  DOUBLE PRECISION, DIMENSION(3), save :: boxitr
END MODULE try

MODULE cellmap_try
  INTEGER, DIMENSION(:), ALLOCATABLE, save :: lead_itr, map_itr, list_itr, jnear_itr
  INTEGER, save :: mcx_itr, mcy_itr, mcz_itr
  DOUBLE PRECISION, save :: cellx_itr, celly_itr, cellz_itr
END MODULE

MODULE movetype
  logical :: ensemble_pres, ensemble_temp, ensemble_exch
  integer, save :: nmovetypes
  CHARACTER(LEN=5), DIMENSION(:), ALLOCATABLE, save :: movetype_name
  INTEGER, dimension(:), allocatable, save :: movetype_ntry, movetype_nsuccess
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, save :: movetype_prob
END MODULE

MODULE coupling_exch
  character(len=1), save :: exch_comp
  character(len=4), dimension(:), allocatable, save :: exch_tcomp
  integer, save :: exch_ncomp
  double precision, save :: xi_val, log_xi2_of_xi1
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

!MODULE ens
!  INCLUDE 'parameter.h'
!  DOUBLE PRECISION :: ener,eold,enew
!  SAVE ener,eold,enew
!END MODULE ens
!
!MODULE ene
!  INCLUDE 'parameter.h'
!  DOUBLE PRECISION :: enernew,enerold,eninew,eniold,eninew2,eniold2,enewtot,eoldtot
!  SAVE enernew,enerold,eninew,eniold,eninew2,eniold2,enewtot,eoldtot
!END MODULE ene
!
!MODULE vcm
!  INCLUDE 'parameter.h'
!  DOUBLE PRECISION :: pres,dvx,delta
!  SAVE pres,dvx,delta
!END MODULE vcm
!
!MODULE pot
!  INCLUDE 'parameter.h'
!  DOUBLE PRECISION :: energy
!  SAVE energy
!END MODULE pot
!
!MODULE initial
!  INCLUDE 'parameter.h'
!  DOUBLE PRECISION, DIMENSION(NBMAX) :: xini,yini,zini,xcmini,ycmini,zcmini
!  SAVE xini,yini,zini,xcmini,ycmini,zcmini
!END MODULE initial

!MODULE sqwell
!  INCLUDE 'parameter.h'
!  DOUBLE PRECISION :: beps, rcut2, rcut, bepspol, rcut3, rcutpol,bepsinter, rcut4, rcutinter
!  SAVE beps, rcut2, rcut, bepspol, rcut3, rcutpol,bepsinter, rcut4, rcutinter
!END MODULE sqwell

