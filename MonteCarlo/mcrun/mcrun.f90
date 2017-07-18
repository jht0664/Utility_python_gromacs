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
  INTEGER :: itype, icon, naver
  double precision :: start, finish, intermediate
  CHARACTER(LEN=100) :: filename_pres, filename_dens, filename_ex, filename_traj_gro, filename_nmol
  ! 2. Declare a variable of type xtcfil
  type(xtcfile) :: xtcf

  start = omp_get_wtime ()
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
  DO ICON = 1, NCON
    guess = rand() ! pick movement algorithm 
    call movetype_select(guess,itype) ! return itype; i-th movement type
    call movetype_run(itype)
    ! run every delta time
    IF(MOD(ICON,NSKIP).EQ.0)THEN
      WRITE(*,*) ">>>>> simulating ", ICON, "th steps"
      naver = naver + 1
      ! call energy_write
      IF(ipres .EQ. 'YES') CALL print_pressure_dV(filename_pres) ! pressure calculation for WR model
      IF(ensemble_pres) CALL print_density(filename_dens) ! density calculation
      IF(ensemble_exch) CALL print_nmol(filename_nmol) ! density calculation
      intermediate = omp_get_wtime () - start
      write(*,'(a,F10.3,1x,a)') " >> current speed ", real(icon)*0.0036/intermediate, &
        " (10^6 steps/hours)"
    ENDIF
    IF((iconfig .EQ. 'YES') .and. (mod(icon,traj_nstep) == 0)) then
      if(ensemble_exch) then
        call save_gro(filename_traj_gro,'a',icon) ! trajectory in gro
      else
        call save_xtc(xtcf,icon) ! save trajectory in xtc
      endif
    endif
    IF((igr .EQ. 'YES') .and. (mod(icon,rdf_nstep) == 0)) CALL gr12_cal ! r(g) calculation
    IF((iex .EQ. 'YES') .and. (mod(icon,ex_nstep) == 0)) CALL print_hifj(filename_ex) ! r(g) calculation
  ENDDO

! end simulation and average properties
  call save_log('mcrun.out')
  call save_ic('composite.fc')
  call save_gro('confout.gro','o',icon)
  IF(IGR .EQ.'YES') call gr12_save('gr.out',ncon/rdf_nstep)

! exit
  call xtcf % close
  finish = omp_get_wtime () - start
  write(*,'(a,1x,F10.3,1x,a)') "final speed = ", real(NCON)*0.0036/finish, " (10^6 steps/hours)"
! deallocate
  call movetype_exit
  call pos_exit
  call cellmap_exit
  if (ensemble_pres) call cellmap_itr_exit

113 FORMAT(1X,'PACKING FRACTION ',3X,F6.4 / &
           1X,'        PRESSURE ',4X,F5.2 )
151 FORMAT (1X / 5X, 'ENERGY AVERAGE', D20.10 / &
                 5X, 'MSR E  AVERAGE', D20.10 / &
                 5X, 'E DIV  AVERAGE', D20.10 )
171 FORMAT(1X,A3,2X,3E18.10)
119 FORMAT(1X,F6.3,2X,6(F15.7,1X))
129 FORMAT(F10.5,2X,6(F15.7,1X))
149 FORMAT(I10,I10,6(F15.7,1X))
158 FORMAT('#',3X,'Energy',8X,'Enthalpy',14X,'Rg',14X,'Re',14X,'Rd',8X,'Rd_of_CM')
159 FORMAT(I10,6(F15.7,1X))
169 FORMAT(20(F8.4))
179 FORMAT(20(I5))
189 FORMAT('   PRESSURE : ',3(F15.7,2X))
199 FORMAT(5(F15.7,1X))
    STOP
END PROGRAM WR_MCRUN

