PROGRAM NFLUC

!  I follow the instruction of James' code! More detail you can see:  
!  XDR Fortran Interface XTC Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/

! Summary: Concentration (or density) Fluctuation (dela_N) in slabs

! Use the xdr interface
  USE xtc_mod, only: xtcfile

  IMPLICIT NONE
! Declare a variable of type xtcfile
  type(xtcfile) :: xtc

! command
  INTEGER :: NUM_OF_ARGS, I_ARG, N_SKIP
  CHARACTER(LEN=20) :: FILE_NAME,INPUT_NAME,OUTPUT_NAME1,OUTPUT_NAME2,MOL_NAME
  LOGICAL :: LookForInput=.FALSE.,LookForOutput1=.FALSE.,LookForOutput2=.FALSE.
  LOGICAL :: LookForMOL=.FALSE.,LookForDZ=.FALSE.,LookForBLOCK=.FALSE.
  LOGICAL :: FileExist
  
  REAL :: D_Z
  INTEGER :: SIZE_BLOCK

! Read
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: AR_DATA, TEMP
  INTEGER :: MAX_DATA = 5000, ICR = 2000, N_DATA, N_ATOM, IERR
  REAL, DIMENSION(:,:), ALLOCATABLE :: BOX, TEMP_BOX

! Calculate density fluctuation
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: MASS, NPAR_AVG, NPAR_STD2
  REAL, DIMENSION(:), ALLOCATABLE :: T_SLICE
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: NPAR
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: COM
  DOUBLE PRECISION :: PI = 3.141592653589793
  INTEGER :: N_INTRA, N_SLICE, I_COM
  INTEGER :: I, J, K, L, INDEX
  REAL :: T_SLICE_AVG, T_SLICE_STD

! Set default
  INPUT_NAME = "trajout.xtc"
  OUTPUT_NAME1 = "N-slab.dat"
  OUTPUT_NAME2 = "N-fluc.dat"
  D_Z = 0.12D0
  SIZE_BLOCK = 5
  MOL_NAME = "bf4"
    
! Command
  NUM_OF_ARGS = command_argument_count()

  IF (NUM_OF_ARGS == 0) THEN
    WRITE(*,*) "USAGE : ./N-fluc.x -i [xtc file] -mol [bf4/pf6] -o1 [file] -o2 [file] -dz [real] -block [integer]"
    WRITE(*,*) "-i trajout.xtc  xtc input file"
    WRITE(*,*) "-o1 N-slab.dat  # of particles in the slab, N(i,t), output file"
    WRITE(*,*) "-o2 N-fluc.dat  Concentation Fluctuation output file"
    WRITE(*,*) "-mol bf4        Which molecules position you saved in the xtc input file (bf4/pf6)"
    WRITE(*,*) "-dz 0.12        Slab thickness (usually cut-off electrostatic interaction is 1.5)"
    WRITE(*,*) "-block 5        size of block (5 = every 5 steps (time) average)"
    WRITE(*,*) "It calcualtes concentration (or denstiy) average and fluctuation in slabs"
    WRITE(*,*) "xtc file should be created by g_traj -jump option"
    WRITE(*,*) "EXAMPLE: make_ndx_mpi -f conf.gro (select only BF4);"
    WRITE(*,*) "         g_traj_mpi -jump -oxt trajout.xtc; ./N-fluc.x -i trajout.xtc"
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
          CASE("-mol")
            LookForMOL=.TRUE. 
            CYCLE
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
              ELSE IF(LookForMOL) THEN
                READ(FILE_NAME,*) MOL_NAME
                LookForMOL = .FALSE.
                CYCLE
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
  PRINT*, "Your command line: N-fluc.x -i ",TRIM(INPUT_NAME),&
   " -o1 ",TRIM(OUTPUT_NAME1)," -o2 ",TRIM(OUTPUT_NAME2),&
   " -mol ",TRIM(MOL_NAME)," -dz ",D_Z,"-block",SIZE_BLOCK


! Read Input file
  PRINT*, "========================"
  PRINT*, "Reading Input File, ",TRIM(INPUT_NAME)

! Initialize it with the name of xtc file you want to read in.
  call xtc%init(TRIM(INPUT_NAME))
  I = 0
  call xtc%read
  IF ( xtc%STAT == 0 ) THEN
    ! check precision of xtc file
    IF( xtc%prec /= 1000.00D0 ) THEN
      PRINT*, "xtc file precision is not 1000.0 (real)!, but ",xtc%prec
      PRINT*, "this program does only support real precision."
      STOP
    ENDIF
!   check x and y of unit cell is same! (I only consider semiisotropic NPZAT ensemble)
    IF(xtc%box(1,1) /= xtc%box(2,2)) THEN
      PRINT*, "pbc box might not for NPZAT. This program does not support this system."
      STOP
    ENDIF
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
    PRINT*, "semiisotropic initial BOX size (",(BOX(J,I),J=1,3),")" 
    ELSE
      PRINT*, TRIM(INPUT_NAME)," has no trajectory. Check again!"
      STOP
  ENDIF
! save coordinates
  DO
    IF(MOD(I+1,500) == 0) PRINT*,"Reading ",I+1,"th time"
    ! read data
    call xtc%read
    !checking
    IF ( xtc%STAT /= 0 ) THEN
      PRINT*, "Reading COMPLETE"
      EXIT
    ENDIF
    IF((xtc%box(1,1) /= BOX(1,0)) .or. (xtc%box(2,2) /= BOX(2,0))) THEN
      PRINT*, "pbc box might be changed. this program does only support NPZAT system."
      STOP
    ENDIF
    I = I + 1
    ! check array memory
    IF(I == MAX_DATA) THEN
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
    AR_DATA(:,:,I) = xtc%pos(:,:)
    DO J=1,3
      BOX(J,I) = xtc%box(J,J)
    EndDO
  ENDDO
  N_DATA = I
  PRINT*, "Total frame is ",N_DATA

! Close the file
  call xtc % close

! Just an example to show what was read in
! write(*,'(a,f12.6,a,i0)') " Time (ps): ", xtc % time, "  Step: ", xtc % STEP
! write(*,'(a,f12.6,a,i0)') " Precision: ", xtc % prec, "  No. Atoms: ", xtc % NATOMS
! write(*,'(3f9.3)') xtc % pos
! This is the same order as found in the GRO format fyi
! write(*,'(11f9.5)') xtc % box(1,1), xtc % box(2,2), xtc % box(3,3), &
!                     xtc % box(1,2), xtc % box(1,3), & 
!                     xtc % box(2,1), xtc % box(2,3), &
!                     xtc % box(3,1), xtc % box(3,2) 

! Calculate COM of the molecule
  PRINT*, "========================"
  PRINT*, "Calculating COM of molecules"

  SELECT CASE(MOL_NAME)
    CASE("bf4") !the molecules are BF4
      N_INTRA = 5
      IF(MOD(N_ATOM,N_INTRA) /= 0) THEN
        PRINT*, "Molecules are not ",TRIM(MOL_NAME)
        STOP
      ENDIF
      ALLOCATE(MASS(1:N_INTRA+1),stat=IERR)
      MASS=0.0D0
      ! mass information of BF4
      MASS(1)=10.811 ! mass of B
      DO I=2,N_INTRA 
        MASS(I)=18.9984
      ENDDO
      ! sum of mass of a BF4 molecule
      DO I=1,N_INTRA 
        MASS(N_INTRA+1)=MASS(N_INTRA+1)+MASS(I)
      ENDDO
    CASE("pf6") !the molecules are PF6
      N_INTRA = 7
      IF(MOD(N_ATOM,N_INTRA) /= 0) THEN
        PRINT*, "Molecules are not ",TRIM(MOL_NAME)
      ENDIF
      ALLOCATE(MASS(1:N_INTRA+1),stat=IERR)
      MASS=0.0D0
      ! mass information of BF4
      MASS(1)=30.97376 ! mass of P
      DO I=2,N_INTRA 
        MASS(I)=18.9984
      ENDDO
      ! sum of mass of a BF4 molecule
      DO I=1,N_INTRA 
        MASS(N_INTRA+1)=MASS(N_INTRA+1)+MASS(I)
      ENDDO
    CASE default
      PRINT*, "Cannot match the name of molecules you set!"
      STOP
  END SELECT

! COM 
  ALLOCATE(COM(1:3,1:N_ATOM/N_INTRA,0:N_DATA),stat=IERR)
  IF (IERR /= 0) STOP "COM: problem in attempt to allocate memory"
  COM = 0.0D0
  DO I=0,N_DATA
    DO J=1,N_ATOM
      ! find centeral atom
      IF(MOD(J-1,N_INTRA) == 0) THEN
        AR_DATA(:,J,I)=AR_DATA(:,J,I)-BOX(:,I)*FLOOR(AR_DATA(:,J,I)/BOX(:,I))
        ! translate long distance atom 
        ELSE 
          AR_DATA(:,J,I)=AR_DATA(:,J,I)-BOX(:,I)*NINT((AR_DATA(:,J,I)-AR_DATA(:,J-MOD(J-1,N_INTRA),I))/BOX(:,I))
      ENDIF
      I_COM = (J-1)/N_INTRA+1
      COM(:,I_COM,I)=COM(:,I_COM,I)+AR_DATA(:,J,I)*MASS(MOD(J-1,N_INTRA)+1)
    EndDO
  EndDO
  COM=COM/MASS(N_INTRA+1)
  PRINT*, "calculating COMPLETE"

! Slicing slabs 
  PRINT*, "========================"
  PRINT*, "Slicing the simulation box and counting molecules"

  OPEN(UNIT=11,FILE=TRIM(OUTPUT_NAME1),STATUS='unknown')  
! # of Slices from initial box_z length
  IF(D_Z <= 0.1D0) THEN
    PRINT*, "thickness of slices, ",D_Z," is too small, < 1 A"
  ENDIF
  N_SLICE = CEILING(BOX(3,0)/D_Z)
  PRINT*, N_SLICE,"# of slices"
  ALLOCATE(NPAR(1:N_SLICE,0:N_DATA),stat=IERR)
  IF (IERR /= 0) STOP "NPAR: problem in attempt to allocate memory"
  NPAR = 0
  ALLOCATE(T_SLICE(0:N_DATA),stat=IERR)
  IF (IERR /= 0) STOP "T_SLICE: problem in attempt to allocate memory"
  T_SLICE = 0
  DO I=0,N_DATA
    T_SLICE(I) = BOX(3,I)/REAL(N_SLICE)
    DO J=1,N_ATOM/N_INTRA
      COM(3,J,I)=COM(3,J,I)-BOX(3,I)*FLOOR(COM(3,J,I)/BOX(3,I))
      INDEX = CEILING(COM(3,J,I)/T_SLICE(I))
!      PRINT*, INDEX, COM(3,J,I)/T_SLICE, CEILING(COM(3,J,I)/T_SLICE), COM(3,J,I)/BOX(3,I),FLOOR(COM(3,J,I)/BOX(3,I))
      IF(INDEX <= 0) THEN
        PRINT*, COM(:,J,I),I,J
        STOP "index of slabs is less than 0 or zero"
      ENDIF
      IF(INDEX > N_SLICE) STOP "index of slabs is greater than MAX"
      NPAR(INDEX,I) = NPAR(INDEX,I) + 1.0D0
!      PRINT*, INDEX,I,NPAR(INDEX,I),INT(NPAR(INDEX,I))
    ENDDO
    DO K=1,N_SLICE
      WRITE(11,'(I8,1X,I3,1X,I5)') I,K,INT(NPAR(K,I))
    ENDDO
    WRITE(11,*)
    WRITe(11,*)
  ENDDO
  CLOSE(11)
  PRINT*, "Writing ",TRIM(OUTPUT_NAME1)," COMPLETE"
  PRINT*, "print NPAR: Time, i-th slab, N_particles in the i-th slab"

  !Block average (time average)
  PRINT*, "========================"
  PRINT*, "Block averaging"
  DO K=1,N_SLICE
    DO I=MOD(N_DATA+1,SIZE_BLOCK),N_DATA-SIZE_BLOCK+1,SIZE_BLOCK ! neglect frist steps in each block
      DO J=I+SIZE_BLOCK-1,I+1,-1 ! sum from last number to first number within size of block
        NPAR(K,J-1) = NPAR(K,J-1)+NPAR(K,J)
      ENDDO
      NPAR(K,I)=NPAR(K,I)/DBLE(SIZE_BLOCK)
!      PRINT*,K,I,NPAR(K,I)
    ENDDO
  ENDDO

  !Calculating avg and std
  PRINT*, "========================"
  PRINT*, "Calculating average and std"
  ALLOCATE(NPAR_AVG(1:N_SLICE),stat=IERR)
  ALLOCATE(NPAR_STD2(1:N_SLICE),stat=IERR)
  NPAR_AVG = 0.0D0
  NPAR_STD2 = 0.0D0
  DO K=1,N_SLICE
    DO I=MOD(N_DATA+1,SIZE_BLOCK),N_DATA-SIZE_BLOCK+1,SIZE_BLOCK
!      PRINT*, "I",I
      NPAR_AVG(K) = NPAR_AVG(K)+NPAR(K,I)
      NPAR_STD2(K) = NPAR_STD2(K)+NPAR(K,I)**2
    ENDDO
    NPAR_AVG(K) = NPAR_AVG(K)/DBLE((N_DATA+1)/SIZE_BLOCK)
    NPAR_STD2(K) = NPAR_STD2(K)/DBLE((N_DATA+1)/SIZE_BLOCK)-NPAR_AVG(K)**2
  ENDDO

  OPEN(UNIT=12,FILE=TRIM(OUTPUT_NAME2),STATUS='unknown')
  DO K=1,N_SLICE
    WRITE(12,'(I5,3(1X,F10.5))') K, NPAR_AVG(K), NPAR_STD2(K), NPAR_STD2(K)/(NPAR_AVG(K)**2)
  ENDDO
  PRINT*, "Writing ",TRIM(OUTPUT_NAME2)," COMPLETE"
  PRINT*, "print NFLUC: i-th slab, avg, std**2, std**2/avg**2" 
  CLOSE(12)


  DO I=0,N_DATA
    T_SLICE_AVG = T_SLICE_AVG+T_SLICE(I)
    T_SLICE_STD = T_SLICE_STD+T_SLICE(I)**2
  EndDO
  T_SLICE_AVG = T_SLICE_AVG/(N_DATA+1)
  T_SLICE_STD = SQRT(T_SLICE_STD/(N_DATA+1)-T_SLICE_AVG**2)
  
  PRINT*, "The following information is not considering block average"
  PRINT*, "avg thickness of slice (nm) =",SUM(T_SLICE)/(N_DATA+1)
  PRINT*, "STD thickness of slice (nm) =",T_SLICE_STD,"(",T_SLICE_STD*100/(SUM(T_SLICE)/(N_DATA+1)),"%)"
  PRINT*, "EXIT N-FLUC PROGRAM"

END PROGRAM NFLUC
