PROGRAM WR_MCRUN
!     -----------------------------------------------------
!     This program is for Widom-Rowlinson LJ mixture system in NVT emsemble
!     MONTE CARLO SIMULATION OF YUKAWA CHAINS IN 
!     THE NPT ENSEMBLE
!     WRITTEN BY ARUN YETHIRAJ
!     ALL RIGHTS RESERVED
!     07/05/2017 converted to modern fortran version by Hyun Tae Jung
!     -----------------------------------------------------

  use omp_lib
  ! 1. Use the xdr interface
  USE xdr, only: xtcfile

  use inp
  use movetype, only: nmovetypes
  use ints, only: nptot

  IMPLICIT NONE
  DOUBLE PRECISION :: rand
  INTEGER :: imol, itype, icon, naver
  LOGICAL :: success
  real :: start, finish
  ! 2. Declare a variable of type xtcfil
  type(xtcfile) :: xtcf

  call cpu_time(start)
! load data
  call read_ic('composite.ic')
  call read_inp('mcrun.inp')
  call movetype_init
! configure
  call save_gro('conf.gro','o') ! overlap if already exist
  IF(iconfig .EQ. 'YES') THEN
    call xtcf % init('traj.xtc', 'w')
    call save_xtc(xtcf)
  ENDIF
  call cellmap_init

  ! initialize properties
  !  energy not necessary for hard sphere
  !  OPEN(UNIT=26,FILE='energy.out',STATUS='UNKNOWN')
  call gr12_init
  IF(ipres .EQ. 'YES') call pressure_dV_init('pressure.log')
    
! run simulation during do loop (end statement label is "1000 continue")
  naver = 0
  DO ICON = 1, NCON
    ! run movement
    imol = INT( NPTOT*RAND() ) + 1 ! pick random particle, (0 <= RAND < 1)
    itype = int(nmovetypes*rand()) + 1 ! pick movement algorithm 
    call try_move(imol,itype)
!    WRITE(*,*) "A"
    ! run every delta time
    IF(MOD(ICON,NSKIP).EQ.0)THEN
!      write(*,*) "L"
      WRITE(*,*) ">>>>> Calculating", ICON, "th configuration"
!      write(*,*) "J"
      naver = naver + 1
!      write(*,*) "G"
      ! call energy_write
      IF(IGR.EQ.'YES') CALL gr12_cal ! r(g) calculation
!      write(*,*) "H"
      IF(ipres .EQ. 'YES') CALL pressure_dV('pressure.log') ! pressure calculation for WR model
!      write(*,*) "I"
      ! save trajectory
      IF(ICONFIG .EQ. 'YES') THEN
        call save_xtc(xtcf)
        !WRITE(*,*) "not supported yet"
        !STOP
      ENDIF   
!      write(*,*) "F"
    ENDIF
  ENDDO

! end simulation and average properties
  call save_log('colloid.out')
  call save_ic('composite.fc')
  call save_gro('confout.gro','o')
  IF(IGR .EQ.'YES') call gr12_save('gr.out',naver)

! exit
  call xtcf % close
  call cpu_time(finish)
  print '("Time = ",f10.3," seconds.")',real(finish-start)
! deallocate
  call movetype_exit
  call pos_exit
  call cellmap_exit

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

