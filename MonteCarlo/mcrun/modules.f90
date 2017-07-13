MODULE pos
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: x, y, z
  CHARACTER(LEN=4), DIMENSION(:), ALLOCATABLE, SAVE :: typs
  DOUBLE PRECISION, DIMENSION(3) :: box
END MODULE pos

MODULE inp
  INTEGER :: ncon, nskip
  CHARACTER(LEN=3) :: iconfig, ipres, igr, iex
  save ncon, nskip, iconfig, ipres, igr, iex
END MODULE inp

MODULE traj
  INTEGER, save :: istep, nconfig ! # step skip to save trajectory
END MODULE traj

MODULE trans
  DOUBLE PRECISION, SAVE :: trans_prob, dlr_a, dlr_b
END MODULE trans

MODULE calc_pres
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, save :: pressure, scale_length
  INTEGER, DIMENSION(:), ALLOCATABLE, save :: noverlap
  DOUBLE PRECISION, SAVE :: ratio_dv_v
  INTEGER, SAVE :: n_dv
END MODULE calc_pres

MODULE calc_ex
  INTEGER, SAVE :: ex_ntry
END MODULE calc_ex

MODULE ints
  INTEGER, save :: nptot, nmol_a, nmon_a, nmol_b, nmon_b
END MODULE ints

MODULE cellmap
  INTEGER, DIMENSION(:), ALLOCATABLE :: lead
  INTEGER, DIMENSION(:), ALLOCATABLE :: map
  INTEGER, DIMENSION(:), ALLOCATABLE :: list
  INTEGER :: mcx, mcy, mcz
  DOUBLE PRECISION, save :: cellx, celly, cellz
  SAVE lead, map, idi, list, mcx, mcy, mcz
END MODULE cellmap

MODULE try
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xitr, yitr, zitr
  CHARACTER(LEN=4), DIMENSION(:), ALLOCATABLE :: titr
  SAVE xitr, yitr, zitr, titr
END MODULE try

MODULE movetype
  CHARACTER(LEN=5), save :: ensemble_name
  INTEGER, dimension(:), allocatable, save :: movetype_ntry, movetype_nsuccess
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, save :: movetype_prob
  integer, save :: nmovetypes
END MODULE

MODULE sigmas
  DOUBLE PRECISION :: sigma_b, sigma_a, sigma_ab
  SAVE sigma_b, sigma_a, sigma_ab
END MODULE sigmas

MODULE rdf
  INCLUDE 'parameter.h'
  DOUBLE PRECISION :: bsz ! bin size
  INTEGER :: nbin
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: fr11,rint,gr11,gr12,gr22
  SAVE nbin,bsz,fr11,rint,gr11,gr12,gr22
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

