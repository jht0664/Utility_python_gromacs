PROGRAM WR_MCRUN
!     -----------------------------------------------------
!     Original Source
!     MONTE CARLO SIMULATION OF YUKAWA CHAINS IN 
!     THE NPT ENSEMBLE
!     WRITTEN BY ARUN YETHIRAJ
!     ALL RIGHTS RESERVED
!
!     This program is for Widom-Rowlinson LJ mixture system in NVT emsemble
!      upgraded by Hyuntae Jung
!     07/05/2017 converted to modern fortran version by Hyun Tae Jung
!     07/14/2017 add calculation of pressure in NVT and pressure coupling
!
!     -----------------------------------------------------

  use omp_lib
  ! 1. Use the xdr interface
  USE xdr, only: xtcfile
  use inp
  use movetype
  use traj, only: traj_nstep
  use rdf, only: rdf_nstep
  use calc_ex, only: ex_nstep

  IMPLICIT NONE
  DOUBLE PRECISION :: rand, guess
  INTEGER :: itype, iicon, icon, naver
  double precision :: start, finish
  CHARACTER(LEN=100) :: filename_pres, filename_dens, filename_ex, filename_traj_gro, filename_nmol
  real :: speed, intermediate, itemp1, itemp2
  ! 2. Declare a variable of type xtcfil
  type(xtcfile) :: xtcf

  start = omp_get_wtime ()
  call omp_set_num_threads(openmp_thread_num)
! load data
  call read_ic('composite.ic')
  call read_inp('mcrun.inp')
  call cellmap_init

! initialize properties and files
  !  energy not necessary for hard sphere
  !  OPEN(UNIT=26,FILE='energy.out',STATUS='UNKNOWN')
  if(igr == 'YES') call gr12_init
  IF(ipres .EQ. 'YES') then
    filename_pres = 'pressure.out'
    call newfile_del_oldfile(filename_pres)
  endif
  if(ensemble_pres) then
    filename_dens = 'density.out'
    call newfile_del_oldfile(filename_dens)
  endif
  if(iex == 'YES') then
    filename_ex = 'hifj.out'
    call newfile_del_oldfile(filename_ex)
  endif
  if(ensemble_exch) then
    filename_traj_gro = 'traj.gro'
    filename_nmol = 'nmol.out'
    call newfile_del_oldfile(filename_traj_gro)
    call newfile_del_oldfile(filename_nmol)
  endif
! configure
  naver = 0
  call save_gro('conf.gro','o',naver) ! overlap if already exist
  IF(iconfig .EQ. 'YES') THEN
    call xtcf % init('traj.xtc','w')
    call save_xtc(xtcf,naver)
  ENDIF
! run simulation during do loop (end statement label is "1000 continue")
  DO IICON = 1, NNCON
  DO ICON = 1, NCON
    guess = rand() ! pick movement algorithm 
    call movetype_select(guess,itype) ! return itype; i-th movement type
    call movetype_run(itype)
    ! run every delta time
    IF(MOD(ICON,NSKIP).EQ.0)THEN
      WRITE(*,*) ">>>>> simulating ", ICON, "th steps"
      naver = naver + 1
      ! call energy_write
      intermediate = real(omp_get_wtime () - start)
      itemp1 = real(icon)/1000.0
      itemp2 = real(ncon)/1000.0
      speed = real((iicon-1)*itemp2+itemp1)*0.0036/intermediate
      write(*,'(1x,a,F10.3,1x,a,F10.3,a)') " >> ", speed, " (10^9 steps/hours) or ", 1.0/speed, "(hours/10^9 steps)"
      write(*,'(1x,a,F10.3,1x,a,F10.3,a)') " >> ", real(nncon*itemp2)/1000000/speed, " hours left"
      IF((ipres .EQ. 'YES')) CALL print_pressure_dV(filename_pres) ! pressure calculation for WR model
      IF(ensemble_pres) CALL print_density(filename_dens) ! density calculation
      IF(ensemble_exch) CALL print_nmol(filename_nmol) ! density calculation
    ENDIF
    IF((iconfig .EQ. 'YES') .and. (mod(icon,traj_nstep) == 0)) then
      if(ensemble_exch) then
        call save_gro(filename_traj_gro,'a',naver) ! trajectory in gro
      else
        call save_xtc(xtcf,naver) ! save trajectory in xtc
      endif
    endif
    IF((igr .EQ. 'YES') .and. (mod(icon,rdf_nstep) == 0)) CALL gr12_cal ! r(g) calculation
    IF((iex .EQ. 'YES') .and. (mod(icon,ex_nstep) == 0)) CALL print_hifj(filename_ex) ! r(g) calculation
  ENDDO
  ENDDO

! end simulation and average properties
  call save_log('mcrun.out')
  call save_ic('composite.fc')
  call save_gro('confout.gro','o',naver)
  IF(IGR .EQ.'YES') call gr12_save('gr.out',ncon*nncon/rdf_nstep)

! exit
  call xtcf % close
  finish = omp_get_wtime () - start
  write(*,'(a,1x,F10.3,1x,a)') "final speed = ", real(NCON)*0.0036/finish, " (10^6 steps/hours)"
! deallocate
  call movetype_exit
  call pos_exit
  call cellmap_exit
!  if (ensemble_pres) call cellmap_itr_exit
  STOP
END PROGRAM WR_MCRUN

