PROGRAM DENSITYPROFILE

!  I follow the instruction of James' code! More detail you can see:  
!  XDR Fortran Interface XTC Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/

! Summary: 2D-Density Profile on XY(R) and Z, setting interface should be on center of box (x,y,z=0)

! Use the xdr interface
  USE xtc_mod, only: xtcfile
  IMPLICIT NONE
! Declare a variable of type xtcfile
  type(xtcfile) :: xtc
! command
  INTEGER :: NUM_OF_ARGS, I_ARG
  CHARACTER(LEN=20) :: FILE_NAME, INPUT_NAME, OUTPUT_NAME1,OUTPUT_NAME2,OUTPUT_NAME3 !MOL_NAME 
  LOGICAL :: LookForInput=.FALSE.,LookForOutput1=.FALSE.,LookForOutput2=.FALSE.,LookForOutput3=.FALSE.
  LOGICAL :: LookForDZ=.FALSE.,LookForBLOCK=.FALSE. !LookForMOL=.FALSE.
  LOGICAL :: FileExist
  INTEGER :: SIZE_BLOCK
  REAL :: D_Z
! Read
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: AR_DATA, TEMP
  INTEGER :: MAX_DATA = 5000, ICR = 2000, N_DATA, N_ATOM, IERR
  REAL, DIMENSION(:,:), ALLOCATABLE :: BOX, TEMP_BOX
! Calculate concentation in slab
  REAL, DIMENSION(:), ALLOCATABLE :: T_SLICE
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CONC
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: COM
  INTEGER :: N_INTRA, N_SLICE, I_COM
  REAL :: T_SLICE_AVG


! Set parameter
  INPUT_NAME = "D-in.xtc"
  OUTPUT_NAME1 = "D-prof-time.dat" ! raw save file
  OUTPUT_NAME2 = "D-prof-pos-int.dat" ! interface fluctuation 
  OUTPUT_NAME3 = "D-prof-pos.dat" ! set x,y,or z = 0 has interface (assume only one interface)
  D_Z = 0.12D0
  SIZE_BLOCK = 5

! =========== START(Update 12/27/2015) =============
! Command line
  NUM_OF_ARGS = command_argument_count()
  IF (NUM_OF_ARGS == 0) THEN
    WRITE(*,*) "USAGE : ./D-prof.x -i [xtc file] -dz [real] -block [integer] -o1 [file] -o2 [file] -o3 [file]" 
    WRITE(*,*) "-i    D-in.xtc         		Input xtc file"
    WRITE(*,*) "-dz   0.12             		Slab thickness" !(usually cut-off electrostatic interaction is 1.4)
    WRITE(*,*) "-block 5               		Size of block (5 = every 5 steps (time) average) for cal. density"
    WRITE(*,*) "-o1   D-prof-time.dat  		Density-position at various time (raw save file)"
    WRITE(*,*) "-o2   D-prof-pos-int.dat 	Position of interface at various time"
    WRITE(*,*) "-o3   D-prof-pos.dat        Average density-position"
    WRITE(*,*) "It calcualtes density profiles."
    WRITE(*,*) "Assume that it has only one interface"
    WRITE(*,*) "xtc file should be created by g_traj -jump option"
    WRITE(*,*) "EXAMPLE: make_ndx_mpi -f conf.gro (select only one kind of particles);"
    WRITE(*,*) "         g_traj_mpi -jump -oxt D-in.xtc; ./D-prof.x -i D-in.xtc"
    STOP
    ELSE IF(NUM_OF_ARGS > 0 ) Then
    !loop across options
      DO I_ARG=1,NUM_OF_ARGS
        CALL get_command_argument(I_ARG,FILE_NAME)
        FILE_NAME = adjustl(FILE_NAME)
        SELECT CASE(FILE_NAME)
          CASE("-i")
            LookForInput=.TRUE. 
            CYCLE
          CASE("-o1")
            LookForOutput1=.TRUE. 
            CYCLE
          CASE("-o2")
            LookForOutput2=.TRUE. 
            CYCLE
          CASE("-o3")
            LookForOutput3=.TRUE. 
            CYCLE
!          CASE("-mol")
!            LookForMOL=.TRUE. 
!            CYCLE
          CASE("-dz")
            LookForDZ=.TRUE.
            CYCLE
          CASE("-block")
            LookForBLOCK=.TRUE.
            CYCLE
          CASE default
            IF(LookForInput) THEN
              inquire(file=TRIM(FILE_NAME),exist=FileExist)
              IF(.NOT.FileExist) THEN
                WRITE(*,*) 'Input File ',FILE_NAME,' is not found'
                STOP
              EndIF
              INPUT_NAME = FILE_NAME
              LookForInput = .FALSE.
              CYCLE
              ELSE IF(LookForOutput1) THEN
                inquire(file=TRIM(FILE_NAME),exist=FileExist)
                IF(FileExist) THEN
                  WRITE(*,*) 'Output File1 ',FILE_NAME,' already exists'
                EndIF
                OUTPUT_NAME1 = FILE_NAME
                LookForOutput1 = .FALSE.
                CYCLE
              ELSE IF(LookForOutput2) THEN
                inquire(file=TRIM(FILE_NAME),exist=FileExist)
                IF(FileExist) THEN
                  WRITE(*,*) 'Output File2 ',FILE_NAME,' already exists'
                EndIF
                OUTPUT_NAME2 = FILE_NAME
                LookForOutput2 = .FALSE.
                CYCLE
              ELSE IF(LookForOutput3) THEN
                inquire(file=TRIM(FILE_NAME),exist=FileExist)
                IF(FileExist) THEN
                  WRITE(*,*) 'Output File3 ',FILE_NAME,' already exists'
                EndIF
                OUTPUT_NAME3 = FILE_NAME
                LookForOutput3 = .FALSE.
                CYCLE
!              ELSE IF(LookForMOL) THEN
!                READ(FILE_NAME,*) MOL_NAME
!                LookForMOL = .FALSE.
!                CYCLE
              ELSE IF(LookForDZ) THEN
                READ(FILE_NAME,*) D_Z
                LookForDZ = .FALSE.
                CYCLE
              ELSE IF(LookForBLOCK) THEN
                READ(FILE_NAME,*) SIZE_BLOCK
                LookForBLOCK = .FALSE.
                CYCLE
              ELSE
                PRINT*, "Option ", FILE_NAME, "is unknown"
                STOP
            ENDIF
         END SELECT
      ENDDO
  ENDIF
  PRINT*, "Your command line: D-prof.x -i ",TRIM(INPUT_NAME)," -o1 ",TRIM(OUTPUT_NAME1)," -o2 ",TRIM(OUTPUT_NAME2),&
          " -o3 ",TRIM(OUTPUT_NAME3)," -dz ",D_Z," -block ",SIZE_BLOCK!" -mol ",TRIM(MOL_NAME),
! =========== END(Update 12/27/2015) =============

! =========== START(Update 12/22/2015) =============
! the reading xtc code is for semi-isotropic box (a,a,c)
! Initially, read Input file
  PRINT*, "========================"
  PRINT*, "Reading Input File, ",TRIM(INPUT_NAME)
  call xtc%init(TRIM(INPUT_NAME)) ! Initialize it with the name of xtc file you want to read in.
  I = 0
  call xtc%read
  IF ( xtc%STAT == 0 ) THEN
    IF( xtc%prec /= 1000.00D0 ) THEN ! check precision of xtc file
      PRINT*, "xtc file precision is not 1000.0 (real)!, but ",xtc%prec
      PRINT*, "this program does only support real precision."
      STOP
    ENDIF
    IF(xtc%box(1,1) /= xtc%box(2,2)) THEN ! only for semiisotropic box and NAPZT ensemble
      PRINT*, "pbc box might not for NPZAT. This program does not support this system."
      STOP
    ENDIF
    !IF((xtc%box(1,1) /= xtc%box(2,2)) .or. &
    !   (xtc%box(1,1) /= xtc%box(3,3)) .or. &
    !   (xtc%box(2,2) /= xtc%box(3,3))) THEN ! only for cubic box and NPT ensemble
    !  PRINT*, "pbc box might not for NPZAT. This program does not support this system."
    !  STOP
    !ENDIF
    ! save original position 
    N_ATOM = xtc%NATOMS
    ALLOCATE(AR_DATA(1:3,1:N_ATOM,0:MAX_DATA),stat=IERR)
    IF (IERR /= 0) STOP "DATA: problem in attempt to allocate memory"
    AR_DATA = 0.0D0
    AR_DATA(:,:,0) = xtc%pos(:,:)
    ALLOCATE(BOX(1:3,0:MAX_DATA),stat=IERR)
    DO J=1,3
      BOX(J,I) = xtc%box(J,J)
    EndDO
    PRINT*, "This initial BOX size (t=0) (",(BOX(J,I),J=1,3),")" 
    ELSE
      PRINT*, TRIM(INPUT_NAME)," has no trajectory. Check again!"
      STOP
  ENDIF
  DO ! save position and pbc cell
    IF(MOD(I+1,500) == 0) PRINT*,"Reading ",I+1,"th time"
    call xtc%read ! read data
    IF ( xtc%STAT /= 0 ) THEN !checking end of the xtc file
      PRINT*, "Reading COMPLETE"
      EXIT
    ENDIF
    IF((xtc%box(1,1) /= BOX(1,0)) .or. (xtc%box(2,2) /= BOX(2,0))) THEN ! check fixed box lengths for NAPZT
      PRINT*, "pbc box might be changed. this program does only support NPZAT system."
      STOP
    ENDIF
    I = I + 1
    IF(I == MAX_DATA) THEN ! check array memory for coordinates (AR_DATA) and box (BOX)
      PRINT*, "Memory increases from size of 3 x ",N_ATOM," x ",SIZE(AR_DATA)/(3*N_ATOM)
      ALLOCATE(TEMP(1:3,1:N_ATOM,0:MAX_DATA+ICR),stat=IERR)
      ALLOCATE(TEMP_BOX(1:3,0:MAX_DATA+ICR),stat=IERR)
      IF (IERR /= 0) STOP "reallocate_temp: problem in attempt to allocate memory"
      TEMP = 0.0D0
      TEMP_BOX = 0.0D0
      DO J=1,3
        DO K=1,N_ATOM
          DO L=0,MAX_DATA
            TEMP(J,K,L) = AR_DATA(J,K,L)
          ENDDO
        ENDDO
      ENDDO
      DO J=1,3
        DO L=0,MAX_DATA
          TEMP_BOX(J,L) = BOX(J,L)
        ENDDO
      ENDDO
      DEALLOCATE(AR_DATA)
      DEALLOCATE(BOX)
      ALLOCATE(AR_DATA(1:3,1:N_ATOM,0:MAX_DATA+ICR),stat=IERR)
      ALLOCATE(BOX(1:3,0:MAX_DATA+ICR),stat=IERR)
      IF (IERR /= 0) STOP "reallocate_array: problem in attempt to allocate memory"
      DO J=1,3
        DO K=1,N_ATOM
          DO L=0,MAX_DATA+ICR
            AR_DATA(J,K,L) = TEMP(J,K,L)
          ENDDO
        ENDDO
      ENDDO
      DO J=1,3
        DO L=0,MAX_DATA+ICR
          BOX(J,L) = TEMP_BOX(J,L)
        ENDDO
      ENDDO
      DEALLOCATE(TEMP)
      DEALLOCATE(TEMP_BOX)
      MAX_DATA = MAX_DATA+ICR
    ENDIF
    AR_DATA(:,:,I) = xtc%pos(:,:) ! save position into newly expanded matrix
    DO J=1,3
      BOX(J,I) = xtc%box(J,J)
    EndDO
  ENDDO
  N_DATA = I
  PRINT*, "Total frame is ",N_DATA
  call xtc % close ! Close the file
! Just an example to show what was read in
! write(*,'(a,f12.6,a,i0)') " Time (ps): ", xtc % time, "  Step: ", xtc % STEP
! write(*,'(a,f12.6,a,i0)') " Precision: ", xtc % prec, "  No. Atoms: ", xtc % NATOMS
! write(*,'(3f9.3)') xtc % pos
! This is the same order as found in the GRO format fyi
! write(*,'(11f9.5)') xtc % box(1,1), xtc % box(2,2), xtc % box(3,3), &
!                     xtc % box(1,2), xtc % box(1,3), & 
!                     xtc % box(2,1), xtc % box(2,3), &
!                     xtc % box(3,1), xtc % box(3,2) 
! =========== END (Update 12/22/2015) =============

! =========== START (Update 1/29/2016) ===========
! Slicing slabs in NAPZT ensemble and get concentration, C(i,t)
  PRINT*, "=================================================="
  PRINT*, "Slicing the simulation box and counting molecules"
  OPEN(UNIT=11,FILE=TRIM(OUTPUT_NAME1),STATUS='unknown')  
  IF(D_Z <= 0.1D0) THEN ! # of Slices from initial box_z length
    PRINT*, "WARNING: Thickness of slices, ",D_Z,"(nm) may be too small, < 1 A"
    PRINT*, "Because of single precision simulation, very small thickness is not recommended."
  ENDIF
  N_SLICE = CEILING(BOX(3,0)/D_Z)
  PRINT*, "set ",N_SLICE,"# of slices"
  ALLOCATE(CONC(1:N_SLICE,0:N_DATA),stat=IERR)
  IF (IERR /= 0) STOP "CONC: problem in attempt to allocate memory"
  CONC = 0
  ALLOCATE(T_SLICE(0:N_DATA),stat=IERR)
  IF (IERR /= 0) STOP "T_SLICE: problem in attempt to allocate memory"
  T_SLICE = 0
  DO I=0,N_DATA
    T_SLICE(I) = BOX(3,I)/REAL(N_SLICE)
    DO J=1,N_ATOM/N_INTRA
      COM(3,J,I)=COM(3,J,I)-BOX(3,I)*FLOOR(COM(3,J,I)/BOX(3,I))
      INDEX = CEILING(COM(3,J,I)/T_SLICE(I))
      IF(COM(3,J,I) == 0.0D0) INDEX = N_SLICE
!      PRINT*, INDEX, COM(3,J,I)/T_SLICE, CEILING(COM(3,J,I)/T_SLICE), COM(3,J,I)/BOX(3,I),FLOOR(COM(3,J,I)/BOX(3,I))
      IF(INDEX <= 0) THEN
        PRINT*, COM(:,J,I),I,J
        STOP "index of slabs is less than 0 or zero"
        ELSE IF(INDEX > N_SLICE) THEN
          STOP "index of slabs is greater than MAX"
      ENDIF
      CONC(INDEX,I) = CONC(INDEX,I) + 1.0D0
!      PRINT*, INDEX,I,CONC(INDEX,I),INT(CONC(INDEX,I))
    ENDDO
    DO K=1,N_SLICE
      CONC(K,I) = CONC(K,I)/(T_SLICE(I)*BOX(2,I)*BOX(2,I))
      WRITE(11,'(I8,1X,I3,1X,I5)') I,K,INT(CONC(K,I))
    ENDDO
    WRITE(11,*)
    WRITe(11,*)
  ENDDO
  CLOSE(11)
  PRINT*, "Writing ",TRIM(OUTPUT_NAME1)," COMPLETE"
  PRINT*, "print CONC: Time, i-th slab, concentration in the i-th slab"
  PRINT*, "The following information is not considering block average"
  T_SLICE_AVG = SUM(T_SLICE)/(N_DATA+1)
  PRINT*, "avg thickness of slice (nm) =",SUM(T_SLICE)/(N_DATA+1)
  ! Block average (time average)
  IF(SIZE_BLOCK /= 1) THEN
    PRINT*, "========================"
    PRINT*, "Block averaging"
    DO K=1,N_SLICE
      DO I=MOD(N_DATA+1,SIZE_BLOCK),N_DATA-SIZE_BLOCK+1,SIZE_BLOCK! skip first time steps of each block
        DO J=I+SIZE_BLOCK-1,I+1,-1! summation from last teim step to first time step within size of block
          CONC(K,J-1)=CONC(K,J-1)+CONC(K,J)
        ENDDO
        CONC(K,I)=CONC(K,I)/DBLE(SIZE_BLOCK)
!        PRINT*,K,I,CONC(K,I)
      ENDDO
    ENDDO
    ELSE
      PRINT*, "Skip Block Averaging"
  ENDIF
! =========== END (Update 1/29/2016) ===========


