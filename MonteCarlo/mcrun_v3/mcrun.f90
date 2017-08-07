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
!     08/03/2017 add calculation of NPTX couplings
!     -----------------------------------------------------
!  use omp_lib
  ! 1. Use the xdr interface
  USE xdr, only: xtcfile
  use inp
  use ints
  use pos
  use filenames
  use movetype
  use coupling_pres
  use coupling_exch
  use traj, only: traj_nstep
  use rdf, only: rdf_nstep
  use time
  use sigmas

  IMPLICIT NONE
  DOUBLE PRECISION :: rand, guess
  INTEGER :: itype, naver, i
  integer(kind=8) :: icon
  double precision :: start, speed, before_run, inter_time1, inter_time2
  external rand
  ! 2. Declare a variable of type xtcfile
  type(xtcfile) :: xtcf

  call time_marker(start)
  !call omp_set_num_threads(openmp_thread_num)
! load data
  call read_ic('composite.ic')
  call read_inp('mcrun.inp')
! initialize cell list and map. If necessary, find the shortest distance of pairs.
  call cellmap_init
!  initialize properties and files
  call files_init
! save initial configure
  naver = 0
  call save_gro('conf.gro','o',naver) ! overlap if already exist
  IF( (iconfig .eq. 'YES') .and. (.not. ensemble_exch) ) THEN
    call xtcf % init('traj.xtc','w')
    call save_xtc(xtcf,naver)
  ENDIF

  a_frac = dble(nmol_a)/dble(nptot)

! run simulation during do loop
  call time_marker(before_run)
  inter_time1 = before_run
  time11 = 0.0D0
  time22 = 0.0D0
  time33 = 0.0D0
  DO ICON = 1, NCON
    guess = rand() ! pick movement algorithm 
    call movetype_select(guess,itype) ! return itype; i-th movement type
    call movetype_run(itype)
! run every delta time
    IF(MOD(ICON,NSKIP).EQ.0)THEN
      WRITE(*,'(a,i0,a)') ">>>>> simulating ", ICON, "th steps"
      naver = naver + 1
      call save_ic('composite.tmp') ! temporary ic file
!     call energy_write
! total time
      call time_marker(inter_time2)
      inter_time1 = inter_time2 - inter_time1
      if(inter_time1 > 1.0D0 )  then
        speed = dble(nskip*36)/10000.0D0/dble(inter_time1)
        write(*,'(a,F10.5,1x,a,F10.3,a)') " >> total speed  = ", speed, " (10^6 steps/hours) "
      endif
      inter_time1 = inter_time2 ! update time
! properties
      IF(ipres .EQ. 'YES') CALL print_pressure_dV(filename_pres) ! pressure calculation for WR model
      IF(ensemble_pres) then
        call update_dvx() ! modify dvx value success fraction to be around 0.5
        CALL print_density(filename_dens) ! density calculation
      endif
      IF(ensemble_exch) then
        if(xi_val_opt) call update_xi_val() ! modify xi_val to be stable with initial a_frac
        CALL print_nmol(filename_nmol) ! print nmol
      endif
      IF((iex .EQ. 'YES')) CALL print_hifj(filename_ex) ! hi/fj calculation
! update nsuccess
      write(*,'(A)') " ====== movetypes info ====== "
      do i=1,nmovetypes
        write(*,'(A,A,I,A,F10.5)') movetype_name(i), &
        " nsucc => ", movetype_i_success(i), &
        ", frac => ", real(movetype_i_success(i))/real(movetype_i_try(i))
        movetype_nsuccess(i) = movetype_nsuccess(i) + movetype_i_success(i)
        movetype_ntry(i) = movetype_ntry(i) + movetype_i_try(i)
        movetype_i_success(i) = 0
        movetype_i_try(i) = 0
      enddo
    ENDIF
! save trajectory
    IF((iconfig .EQ. 'YES') .and. (mod(icon,traj_nstep) == 0)) then
      if(ensemble_exch) then
        call save_gro(filename_traj_gro,'a',naver) ! trajectory in gro
      else
        call save_xtc(xtcf,naver) ! save trajectory in xtc
      endif
    endif
! r(g) calculation
    IF((igr .EQ. 'YES') .and. (mod(icon,rdf_nstep) == 0)) CALL gr12_cal 
  ENDDO

! end simulation and average properties
  call save_log('mcrun.out')
  call save_ic('composite.fc')
  call save_gro('confout.gro','o',naver)
!  IF(IGR .EQ.'YES') call gr12_save('gr.out',ncon*nncon/rdf_nstep)
  IF(IGR .EQ.'YES') call gr12_save('gr.out',ncon/rdf_nstep)

! exit
  call xtcf % close
! deallocate
  call movetype_exit
  call pos_exit
  IF(ipres .EQ. 'YES') call calc_pres_exit
  call cellmap_exit
  STOP
END PROGRAM WR_MCRUN

subroutine files_init()
  use filenames
  use movetype
  use inp
  implicit none
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
  return
end subroutine files_init

subroutine time_marker(current_time)
  implicit none
  integer, dimension(8) :: time_array_0
  double precision :: current_time
  call date_and_time(values=time_array_0)
  current_time = (time_array_0 (4) * 24.0 + time_array_0 (5)) * 3600.0 &
    + time_array_0 (6) * 60.0 + time_array_0 (7) + 0.001 * dble(time_array_0 (8))
  return
end subroutine time_marker
