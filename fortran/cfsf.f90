PROGRAM CFSF

!  I follow the instruction of James' code! More detail you can see:  
!  XDR Fortran Interface XTC Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/

! Summary: Structure factor of concentration fluctuation for BF4/PF6 molecules

! Use the xdr interface
  USE xtc_mod, only: xtcfile

  IMPLICIT NONE
! Declare a variable of type xtcfile
  type(xtcfile) :: xtc

! command
  INTEGER :: NUM_OF_ARGS, I_ARG
  CHARACTER(LEN=20) :: MOL_NAME, FILE_NAME, INPUT_NAME, OUTPUT_NAME1,OUTPUT_NAME2,OUTPUT_NAME3
  LOGICAL :: LookForInput=.FALSE.,LookForOutput1=.FALSE.,LookForOutput2=.FALSE.,LookForOutput3=.FALSE.
  LOGICAL :: LookForMOL=.FALSE.,LookForDZ=.FALSE.,LookForBLOCK=.FALSE.
  LOGICAL :: FileExist
  INTEGER :: SIZE_BLOCK
  REAL :: D_Z

! Read
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: AR_DATA, TEMP
  INTEGER :: MAX_DATA = 5000, ICR = 2000, N_DATA, N_ATOM, IERR
  REAL, DIMENSION(:,:), ALLOCATABLE :: BOX, TEMP_BOX

! Calculate concentation
  REAL, DIMENSION(:), ALLOCATABLE :: MASS, T_SLICE
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CONC
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: COM
  INTEGER :: N_INTRA, N_SLICE, I_COM
  REAL :: T_SLICE_AVG

! Calculate concentration fluctuation
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CONC_AVG, CONC_DEL_AVG, CONC_DEL_STD
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CONC_DEL
  DOUBLE PRECISION :: PI = 3.141592653589793
  INTEGER :: I, J, K, L, INDEX, N_BLOCK, R_LAG
  
! Fourier transform
  INTEGER :: NBIN_R
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AUTOCORR
  DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: CCFSF
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TFKT

!  For FFTW3
  INCLUDE "fftw3.f"
  DOUBLE PRECISION :: PLAN

! Set default
  INPUT_NAME = "cfsf.xtc"
  OUTPUT_NAME1 = "conc.dat"
  OUTPUT_NAME2 = "delta_conc.dat"
  OUTPUT_NAME3 = "cfsf.dat"
  D_Z = 0.12D0
  SIZE_BLOCK = 1
  MOL_NAME = "bf4"

! =========== START(Update 12/27/2015) =============
! Command line
  NUM_OF_ARGS = command_argument_count()
  IF (NUM_OF_ARGS == 0) THEN
    WRITE(*,*) "USAGE : ./cfsf.x -i [xtc file] -mol [bf4/pf6/atom] -dz [real] -block [integer] -o1 [file] -o2 [file] -o3 [file]"
    WRITE(*,*) "-i    cfsf.xtc         Input xtc file"
    WRITE(*,*) "-mol  bf4              Which molecule (bf4/pf6) or atom you saved in xtc file"
    WRITE(*,*) "-dz   0.12             Slab thickness (usually cut-off electrostatic interaction is 1.5)"
    WRITE(*,*) "-block 1               Size of block (5 = every 5 steps (time) average)"
    WRITE(*,*) "-o1   conc.dat         Concentration of the molecules in i-th slabs"
    WRITE(*,*) "-o2   delta_conc.dat   Delta_conc ot those"
    WRITE(*,*) "-o3   cfsf.dat         Structure factor of concentration fluctuation"
    WRITE(*,*) "It calcualtes structure factor of concentration fluctuation (from correlation fn. of conc.)"
    WRITE(*,*) "Assume that all particle is identical and only concern z-direction (slices)"
    WRITE(*,*) "xtc file should be created by g_traj -jump option"
    WRITE(*,*) "EXAMPLE: make_ndx_mpi -f conf.gro (select only BF4);"
    WRITE(*,*) "         g_traj_mpi -jump -oxt trajout.xtc; ./cfsf.x -i trajout.xtc"
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
              ELSE IF(LookForOutput3) THEN
                inquire(file=TRIM(FILE_NAME),exist=FileExist)
                IF(FileExist) THEN
                  WRITE(*,*) 'Output File3 ',FILE_NAME,' already exists'
                EndIF
                OUTPUT_NAME3 = FILE_NAME
                LookForOutput3 = .FALSE.
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
  PRINT*, "Your command line: cfsf.x -i ",TRIM(INPUT_NAME)," -o1 ",TRIM(OUTPUT_NAME1)," -o2 ",TRIM(OUTPUT_NAME2),&
          " -o3 ",TRIM(OUTPUT_NAME3)," -mol ",TRIM(MOL_NAME)," -dz ",D_Z," -block ",SIZE_BLOCK
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

! =========== START (Update 12/22/2015) ===========
! You need to set MOL_NAME to determine what molecule you saved in xtc file
! Then, you will get center of mass of molecules into variable COM
! Calculate mass of the molecule
  PRINT*, "======================================"
  PRINT*, "Calculating COM of identical molecules"
  SELECT CASE(MOL_NAME)
    CASE("bf4") ! the molecule is BF4
      N_INTRA = 5
      IF(MOD(N_ATOM,N_INTRA) /= 0) THEN
        PRINT*, "Molecules are not ",TRIM(MOL_NAME)
      ENDIF
      ALLOCATE(MASS(1:N_INTRA+1),stat=IERR)
      MASS=0.0D0
      MASS(1)=10.811 ! mass of B
      DO I=2,N_INTRA 
        MASS(I)=18.9984 ! mass of F
      ENDDO
      ! sum of mass of a BF4 molecule
      DO I=1,N_INTRA 
        MASS(N_INTRA+1)=MASS(N_INTRA+1)+MASS(I)
      ENDDO
    CASE("pf6") ! the molecule is PF6
      N_INTRA = 7
      IF(MOD(N_ATOM,N_INTRA) /= 0) THEN
        PRINT*, "Molecules are not ",TRIM(MOL_NAME)
      ENDIF
      ALLOCATE(MASS(1:N_INTRA+1),stat=IERR)
      MASS=0.0D0
      MASS(1)=30.97376 ! mass of P
      DO I=2,N_INTRA 
        MASS(I)=18.9984 ! mass of F
      ENDDO
      ! sum of mass of a BF4 molecule
      DO I=1,N_INTRA 
        MASS(N_INTRA+1)=MASS(N_INTRA+1)+MASS(I)
      ENDDO
    CASE("atom") ! select only atom (ex. Oxygen of PEO chains)
      N_INTRA = 1
      ALLOCATE(MASS(1:N_INTRA+1),stat=IERR)
      MASS(1)=1.D0
      DO I=1,N_INTRA 
        MASS(N_INTRA+1)=MASS(N_INTRA+1)+MASS(I)
      ENDDO
    CASE default
      PRINT*, "Cannot match the name of molecules you set!"
      STOP
  END SELECT
! calculate COM
  ALLOCATE(COM(1:3,1:N_ATOM/N_INTRA,0:N_DATA),stat=IERR)
  IF (IERR /= 0) STOP "COM: problem in attempt to allocate memory"
  COM = 0.0D0
  DO I=0,N_DATA
    DO J=1,N_ATOM
      IF(MOD(J-1,N_INTRA) == 0) THEN  ! find central atom (B or P) to avoid beyond box size
        AR_DATA(:,J,I)=AR_DATA(:,J,I)-BOX(:,I)*FLOOR(AR_DATA(:,J,I)/BOX(:,I))
        ELSE ! If it is not central atom (B of P), translate rest of atoms depending on central atom's coordinate         
          AR_DATA(:,J,I)=AR_DATA(:,J,I)-BOX(:,I)*NINT((AR_DATA(:,J,I)-AR_DATA(:,J-MOD(J-1,N_INTRA),I))/BOX(:,I))
      ENDIF
      I_COM = (J-1)/N_INTRA+1
      COM(:,I_COM,I)=COM(:,I_COM,I)+AR_DATA(:,J,I)*MASS(MOD(J-1,N_INTRA)+1)
    EndDO
  EndDO
  COM=COM/MASS(N_INTRA+1)
  PRINT*, "calculating COMPLETE for COM"
! =========== END (Update 12/22/2015) ===========

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

! =========== START (Update 12/23/2015) ===========
  ! Calculating concentration fluctuation (auto correlation)
  PRINT*, "========================"
  PRINT*, "Calculating spatial average of <C(t)> and delta_C(i,t)"
  OPEN(UNIT=21,FILE=TRIM(OUTPUT_NAME2),STATUS='unknown') 
  N_BLOCK = FLOOR(DBLE(N_DATA+1)/DBLE(SIZE_BLOCK)) ! number of time blocks
  ALLOCATE(CONC_AVG(1:N_BLOCK),stat=IERR)
  ALLOCATE(CONC_DEL(1:N_SLICE,1:N_BLOCK),stat=IERR)
  ALLOCATE(CONC_DEL_AVG(1:N_BLOCK),stat=IERR)
  ALLOCATE(CONC_DEL_STD(1:N_BLOCK),stat=IERR)
  CONC_AVG = 0.0D0
  CONC_DEL = 0.0D0
  CONC_DEL_AVG = 0.0D0
  CONC_DEL_STD = 0.0D0
  J=0
  DO I=MOD(N_DATA+1,SIZE_BLOCK),N_DATA-SIZE_BLOCK+1,SIZE_BLOCK
    J=J+1
    DO K=1,N_SLICE ! spatial sum
      CONC_AVG(J) = CONC_AVG(J)+CONC(K,I)
    ENDDO
    CONC_AVG(J) = CONC_AVG(J)/DBLE(N_SLICE)
    DO K=1,N_SLICE ! calculate delta_C(i,t)
      CONC_DEL(K,J) = CONC(K,I)-CONC_AVG(J)
      CONC_DEL_AVG(J) = CONC_DEL_AVG(J)+CONC_DEL(K,J)
      CONC_DEL_STD(J) = CONC_DEL_STD(J)+CONC_DEL(K,J)**2
    ENDDO
    CONC_DEL_AVG(J) = CONC_DEL_AVG(J)/DBLE(N_SLICE)
    CONC_DEL_STD(J) = SQRT(CONC_DEL_STD(J)/DBLE(N_SLICE)-CONC_DEL_AVG(J)**2)
    WRITE(21,'(I5,1X,F10.5,1X,F10.5)') J,CONC_DEL_AVG(J),CONC_DEL_STD(J)
  ENDDO
  IF(J /= N_BLOCK) STOP "ERROR: wrong assignment for memory of CONC_AVG and CONC_DEL. Check your code!"
  PRINT*, "Writing ",TRIM(OUTPUT_NAME2)," COMPLETE"
  PRINT*, "print delta_conc: block index, avg. delta_conc, std. delta_conc" 
  CLOSE(21)
  ! Equal-time (spatial) correlation Function
  ! \left< \Delta c({ r }_{ 1 })\Delta c({ r }_{ 2 }) \right> 
  ! i.e. \Delta c(r)=c(r)-\left< c \right> 
  R_LAG = FLOOR(DBLE(N_SLICE-1)/2.0D0)
  ALLOCATE(AUTOCORR(1:R_LAG),stat=IERR)
  AUTOCORR = 0.0D0
  DO J=1,N_BLOCK ! total time blocks
    DO K=1,R_LAG ! length of slices (=slice lag)
      DO I=1,N_SLICE
        AUTOCORR(K)=AUTOCORR(K)+(CONC_DEL(I,J)-CONC_DEL_AVG(J))*(CONC_DEL(MOD(I+K-1,N_SLICE)+1,J)-CONC_DEL_AVG(J))/(CONC_DEL_STD(J)**2) ! MOD(I+K-1,N_SLICE)+1 due to matrix index 1:N_SLICE, not 0:N_SLICE-1
      EndDO
    EndDO
  EndDO
  DO K=1,R_LAG
    AUTOCORR(K)=AUTOCORR(K)/DBLE(N_SLICE)/DBLE(N_BLOCK)
    PRINT*, K,AUTOCORR(K)
  EndDO
! =========== END (Update 12/23/2015) ===========

! =========== START (Update 12/23/2015) ===========
  ! Fourier Transform of the correlation function (above)
  NBIN_R = FLOOR(DBLE(N_SLICE-1)/2.0D0)
  ALLOCATE(CCFSF(1:NBIN_R/2+1),stat=IERR)
  IF (IERR /= 0) STOP "CCFSF: problem in attempt to allocate memory"
  CCFSF = 0.0D0
  ALLOCATE(TFKT(1:NBIN_R/2+1),stat=IERR)
  IF (IERR /= 0) STOP "TFKT: problem in attempt to allocate memory"
  TFKT = 0.0D0
  CALL dfftw_plan_dft_r2c_1d(PLAN,NBIN_R,AUTOCORR,CCFSF,FFTW_ESTIMATE) 
  CALL dfftw_execute_dft_r2c(PLAN,AUTOCORR,CCFSF)
  CALL dfftw_destroy_plan(PLAN)
  DO K=1,NBIN_R/2+1
    TFKT(K) = CDABS(CCFSF(K))
  ENDDO
  OPEN(UNIT=31,FILE=TRIM(OUTPUT_NAME3),STATUS='unknown')  
  DO K=1,NBIN_R/2+1 ! The value at the longest distance seems to be wrong, so I remove k=1
    WRITE(31,'(F10.5,1X,F10.5)') DBLE(K),REAL(TFKT(K))
  ENDDO
  PRINT*, "Writing ",TRIM(OUTPUT_NAME3)," COMPLETE"
  PRINT*, "print abs(CCFSF): q, cdabs(CCFSF)" 
  CLOSE(31)
! =========== END (Update 12/23/2015) ===========
  PRINT*, "EXIT CFSF PROGRAM"

END PROGRAM CFSF
