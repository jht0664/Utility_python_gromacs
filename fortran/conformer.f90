PROGRAM conformer

! I follow the instruction of James' code! More detail you can see:  
! XDR Fortran Interface XTC Example Program with Wrapper
! 2014 (c) James W. Barnett <jbarnet4@tulane.edu>
! https://github.com/wesbarnett/

! Summary: Dihedral angles of PEOs, e.g. C-O-C-C, O-C-C-O, and C-C-O-C of C-O-C-C-O-C.

! Use the xdr interface
  USE xtc_mod, only: xtcfile

  IMPLICIT NONE
! Declare a variable of type xtcfile
  type(xtcfile) :: xtc

! command
  INTEGER :: NUM_OF_ARGS, I_ARG
  CHARACTER(LEN=20) :: FILE_NAME,INPUT_NAME,OUTPUT_NAME1
  LOGICAL :: LookForInput=.FALSE.,LookForOutput1=.FALSE.
  LOGICAL :: LookForBLOCK=.FALSE.,LookForPOLYMERIZATION=.FALSE.
  LOGICAL :: FileExist
  INTEGER :: N_DOP ! Degree of polymerization, ! SIZE_BLOCK

! Read
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: AR_DATA, TEMP
  INTEGER :: MAX_DATA = 5000, ICR = 2000, N_DATA, N_ATOM, IERR
  REAL, DIMENSION(:,:), ALLOCATABLE :: BOX, TEMP_BOX

! Calculate dihedral angle
  DOUBLE PRECISION, DIMENSION(1:3) :: B1B2, B2B3, IMG_Y
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: DIH
  DOUBLE PRECISION :: PI = 3.141592653589793
  DOUBLE PRECISION :: RE_X
  INTEGER :: I, J, K, L, INDEX_J, NCOMP_MON
  REAL :: NOR_VEC
  character(len=40) :: myfmt

! Analysis conformation of PEO
  INTEGER, DIMENSION(1:3) :: NGNG, PGPG, TRTR, NGPG, NGTR, PGTR
  REAL, DIMENSION(1:3) :: FORM

! Set default
  INPUT_NAME = "dih.xtc"
  OUTPUT_NAME1 = "dih.dat"
  N_DOP = 9
!  SIZE_BLOCK = 5
    
! =========== START(Update 12/27/2015) =============
! Command
  NUM_OF_ARGS = command_argument_count()
  IF (NUM_OF_ARGS == 0) THEN
    WRITE(*,*) "USAGE : ./conformer.x -i [xtc file] -o1 [file] -npol [integer]" !-block [integer]
    WRITE(*,*) "-i      dih.xtc   xtc input file(only contain C1, C2, and O of PEO chains)"
    WRITE(*,*) "-o1     dih.dat   Dihedral output file for PEO"
!    WRITE(*,*) "-block  5         Size of block (5 = every 5 steps (time) average)"
    WRITE(*,*) "-npol   9         Number of degree of polymerization of a (monodisperse) PEO chain (9 = 9-mers)"
    WRITE(*,*) "It calculates dihedral angels of PEOs"
    WRITE(*,*) "Assume that xtc file created by g_traj -nojump option"
    WRITE(*,*) "EXAMPLE: make_ndx_mpi -f conf.gro (for PEO, a select group is C1, C2, and O);"
    WRITE(*,*) "         g_traj_mpi -nojump -n index.ndx -oxt dih.xtc; ./conformer.x -i dih.xtc"
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
          CASE("-block")
            LookForBLOCK=.TRUE. 
            CYCLE
          CASE("-npol")
            LookForPOLYMERIZATION=.TRUE. 
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
!              ELSE IF(LookForBLOCK) THEN
!                READ(FILE_NAME,*) SIZE_BLOCK
!                LookForBLOCK = .FALSE.
!                CYCLE
              ELSE IF(LookForPOLYMERIZATION) THEN
                READ(FILE_NAME,*) N_DOP
                LookForPOLYMERIZATION = .FALSE.
                CYCLE
              ELSE
                PRINT*, "Option ", FILE_NAME, "is unknown"
                STOP
            ENDIF
         END SELECT
      ENDDO
  ENDIF
  PRINT*, "Your command line: conformer.x -i ",TRIM(INPUT_NAME),&
          " -o1 ",TRIM(OUTPUT_NAME1)," -npol ",N_DOP!," -block ",SIZE_BLOCK
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

! =========== START (Update 1/13/2016) ===========
! Calculate Dihedral Angle of PEO
! You need to set N_DOP to determine how long the length of a chain.
  PRINT*, "========================"
  PRINT*, "Calculating Dihedral Angle"
  OPEN(UNIT=12,FILE=TRIM(OUTPUT_NAME1),STATUS='unknown')
! OPEN(UNIT=13,FILE=TRIM(OUTPUT_NAME1)//".2",STATUS='unknown')
  NCOMP_MON = 3 ! number of particles in a monomer (for PEO, 3 because of C1, C2, and O)
  PRINT*, "# of PEO polymers", N_ATOM/(NCOMP_MON*(N_DOP+1))
  PRINT*, "# of degree of polymerization of PEO", N_DOP
! Unit bond vectors (normalize)
  DO I=0,N_DATA ! Time steps
    DO J=1,N_ATOM/(NCOMP_MON*(N_DOP+1)) ! # of polymer chains
      INDEX_J=(J-1)*NCOMP_MON*N_DOP  
      DO K=1,NCOMP_MON*N_DOP-1 ! index of particles within each polymer chain
        AR_DATA(:,INDEX_J+K,I)=AR_DATA(:,INDEX_J+K+1,I)-AR_DATA(:,INDEX_J+K,I)
        AR_DATA(:,INDEX_J+K,I)=AR_DATA(:,INDEX_J+K,I)-BOX(:,I)*NINT(AR_DATA(:,INDEX_J+K,I)/BOX(:,I))
        NOR_VEC=SQRT(AR_DATA(1,INDEX_J+K,I)**2+AR_DATA(2,INDEX_J+K,I)**2+AR_DATA(3,INDEX_J+K,I)**2)
        AR_DATA(:,INDEX_J+K,I)=AR_DATA(:,INDEX_J+K,I)/NOR_VEC
      ENDDO
    ENDDO
  ENDDO
  ALLOCATE(DIH(1:N_ATOM/(NCOMP_MON*(N_DOP+1)),1:NCOMP_MON*N_DOP-3,0:N_DATA),stat=IERR)
  IF (IERR /= 0) STOP "DIH: problem in attempt to allocate memory"
! dihedral angle using atan2 function. dihedral angle ranges [-pi:pi]
  DO I=0,N_DATA
    DO J=1,N_ATOM/(NCOMP_MON*(N_DOP+1))
      INDEX_J=(J-1)*NCOMP_MON*N_DOP  
      DO K=1,NCOMP_MON*N_DOP-3
        ! b1 x b2
        B1B2(1) = AR_DATA(2,INDEX_J+K,I)*AR_DATA(3,INDEX_J+K+1,I)-AR_DATA(3,INDEX_J+K,I)*AR_DATA(2,INDEX_J+K+1,I)
        B1B2(2) = AR_DATA(3,INDEX_J+K,I)*AR_DATA(1,INDEX_J+K+1,I)-AR_DATA(1,INDEX_J+K,I)*AR_DATA(3,INDEX_J+K+1,I)
        B1B2(3) = AR_DATA(1,INDEX_J+K,I)*AR_DATA(2,INDEX_J+K+1,I)-AR_DATA(2,INDEX_J+K,I)*AR_DATA(1,INDEX_J+K+1,I)
        ! b2 x b3
        B2B3(1) = AR_DATA(2,INDEX_J+K+1,I)*AR_DATA(3,INDEX_J+K+2,I)-AR_DATA(3,INDEX_J+K+1,I)*AR_DATA(2,INDEX_J+K+2,I)
        B2B3(2) = AR_DATA(3,INDEX_J+K+1,I)*AR_DATA(1,INDEX_J+K+2,I)-AR_DATA(1,INDEX_J+K+1,I)*AR_DATA(3,INDEX_J+K+2,I)
        B2B3(3) = AR_DATA(1,INDEX_J+K+1,I)*AR_DATA(2,INDEX_J+K+2,I)-AR_DATA(2,INDEX_J+K+1,I)*AR_DATA(1,INDEX_J+K+2,I)
        ! for y of atan2(y,x) function
        IMG_Y(1) = (B1B2(2)*B2B3(3)-B1B2(3)*B2B3(2))*AR_DATA(1,INDEX_J+K+1,I)
        IMG_Y(2) = (B1B2(3)*B2B3(1)-B1B2(1)*B2B3(3))*AR_DATA(2,INDEX_J+K+1,I)
        IMG_Y(3) = (B1B2(1)*B2B3(2)-B1B2(2)*B2B3(1))*AR_DATA(3,INDEX_J+K+1,I)
        ! for x of atan2(y,x) function
        RE_X = B1B2(1)*B2B3(1)+B1B2(2)*B2B3(2)+B1B2(3)*B2B3(3)
        DIH(J,K,I) = atan2(IMG_Y(1)+IMG_Y(2)+IMG_Y(3),RE_X)
      ENDDO
    ENDDO
  ENDDO
  write(myfmt, '("(I5,1X,I4,",I0,"(1X,F10.5))")') NCOMP_MON*N_DOP-3 !unfixed write format
  PRINT*, "format",TRIM(myfmt)
  DO I=0,N_DATA
    DO J=1,N_ATOM/(NCOMP_MON*(N_DOP+1))
      WRITE(12,myfmt) I,J,(DIH(J,K,I)*180.0D0/PI,K=1,NCOMP_MON*N_DOP-3) ! check evernote!
    ENDDO
  EndDO
  PRINT*, "print DIH: Time, Polymer Index, Dihedral sets (bond index [1:",NCOMP_MON*N_DOP-3,"])"
  PRINT*, "Writing ",TRIM(OUTPUT_NAME1)," COMPLETE"
!  DO I=0,N_DATA
!    DO K=1,N_ATOM-3
!      WRITE(13,'(2(I5,1X),F10.5)') I,K,DIH(K,I)*180.0D0/PI
!    ENDDO
!    WRITE(13,*)
!    WRITE(13,*)
!    IF(MOD(I,1000) == 0) PRINT*,  "Writing on ",I, "th line..."
!  EndDO
!  PRINT*, "print DIH: Time, Bond index, Dihedral angles"
!  PRINT*, "Writing ",TRIM(OUTPUT_NAME1),".2 COMPLETE"
  CLOSE(12)
  CLOSE(13)
! =========== END (Update 1/13/2016) =============

! =========== START (Update 1/14/2016) ===========
! Correlation between trans and gaushe
  ! translate the range of dihedral angles from [-180:180] to [0:360]
  DO I=0,N_DATA
    DO J=1,N_ATOM/(NCOMP_MON*(N_DOP+1))
      DO K=1,NCOMP_MON*N_DOP-3
        IF(DIH(J,K,I) < 0.0D0) THEN
          DIH(J,K,I) = DIH(J,K,I) + PI
        EndIF
      ENDDO
    ENDDO
  ENDDO
  ! logical dividing 3 connection
  DO I=0,N_DATA
    DO J=1,N_ATOM/(NCOMP_MON*(N_DOP+1))
      DO K=1,NCOMP_MON*N_DOP-3,3
        FORM = 0.0D0
        IF(DIH(J,K+1,I)*180.0D0/PI < 120.0D0) THEN !second dihedral angle in each monomer, 120 comes from anlges at lowest peak heights
          FORM(2) = 1 ! guashe -
          ELSE IF(DIH(J,K+1,I)*180.0D0/PI < 300.0D0) THEN
            FORM(2) = 3 ! trans
          ELSE
            FORM(2) = 2 ! gaushe +
        ENDIF
        IF(DIH(J,K,I)*180.0D0/PI < 125.0D0) THEN !first dihedral angle in each monomer, 125 comes from anlges at lowest peak heights
          FORM(1) = -1 ! guashe -
          ELSE IF(DIH(J,K,I)*180.0D0/PI < 305.0D0) THEN
            FORM(1) = 0.5 ! trans
          ELSE 
            FORM(1) = 1 ! gaushe +
        ENDIF
        IF(DIH(J,K+2,I)*180.0D0/PI < 125.0D0) THEN !first dihedral angle in each monomer, 125 comes from anlges at lowest peak heights
          FORM(3) = -1 ! guashe -
          ELSE IF(DIH(J,K+2,I)*180.0D0/PI < 305.0D0) THEN
            FORM(3) = 0.5 ! trans
          ELSE 
            FORM(3) = 1 ! gaushe +
        ENDIF
        IF(FORM(2) == 0) STOP "ERROR: Something missing in calculating dihedral angles"
        IF(FORM(1)+FORM(3) == -2) THEN
          NGNG(INT(FORM(2))) = NGNG(INT(FORM(2)))+1 ! NG is negative gaushe
          ELSE IF (FORM(1)+FORM(3) == 2) THEN
            PGPG(INT(FORM(2))) = PGPG(INT(FORM(2)))+1 ! PG is posive gaushe
          ELSE IF (FORM(1)+FORM(3) == 1) THEN
            TRTR(INT(FORM(2))) = TRTR(INT(FORM(2)))+1 ! TR is trans
          ELSE IF (FORM(1)+FORM(3) == 0) THEN
            NGPG(INT(FORM(2))) = NGPG(INT(FORM(2)))+1
          ELSE IF (FORM(1)+FORM(3) == -0.5) THEN
            NGTR(INT(FORM(2))) = NGTR(INT(FORM(2)))+1
          ELSE IF (FORM(1)+FORM(3) == 1.5) THEN
            PGTR(INT(FORM(2))) = PGTR(INT(FORM(2)))+1
          ELSE
            PRINT*, I, J, K, FORM(1)+FORM(3) 
            STOP "ERROR: Something missing in calculating dihedral angles2"
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  ! normalize
  PRINT*, "======================================"
  PRINT*, "Probability of conformations of PEO"
  PRINT*, "(g-)(g-)(g-) = ",REAL(NGNG(1))/REAL((N_DATA+1)*(N_ATOM/(NCOMP_MON*(N_DOP+1)))*(NCOMP_MON*N_DOP-3)/3)
  PRINT*, "(g-)(g+)(g-) = ",REAL(NGNG(2))/REAL((N_DATA+1)*(N_ATOM/(NCOMP_MON*(N_DOP+1)))*(NCOMP_MON*N_DOP-3)/3)
  PRINT*, "(g-)(t )(g-) = ",REAL(NGNG(3))/REAL((N_DATA+1)*(N_ATOM/(NCOMP_MON*(N_DOP+1)))*(NCOMP_MON*N_DOP-3)/3)
  PRINT*, ""
  PRINT*, "(g+)(g-)(g+) = ",REAL(PGPG(1))/REAL((N_DATA+1)*(N_ATOM/(NCOMP_MON*(N_DOP+1)))*(NCOMP_MON*N_DOP-3)/3)
  PRINT*, "(g+)(g+)(g+) = ",REAL(PGPG(2))/REAL((N_DATA+1)*(N_ATOM/(NCOMP_MON*(N_DOP+1)))*(NCOMP_MON*N_DOP-3)/3)
  PRINT*, "(g+)(t )(g+) = ",REAL(PGPG(3))/REAL((N_DATA+1)*(N_ATOM/(NCOMP_MON*(N_DOP+1)))*(NCOMP_MON*N_DOP-3)/3)
  PRINT*, ""
  PRINT*, "(t )(g-)(t ) = ",REAL(TRTR(1))/REAL((N_DATA+1)*(N_ATOM/(NCOMP_MON*(N_DOP+1)))*(NCOMP_MON*N_DOP-3)/3)
  PRINT*, "(t )(g+)(t ) = ",REAL(TRTR(2))/REAL((N_DATA+1)*(N_ATOM/(NCOMP_MON*(N_DOP+1)))*(NCOMP_MON*N_DOP-3)/3)
  PRINT*, "(t )(t )(t ) = ",REAL(TRTR(3))/REAL((N_DATA+1)*(N_ATOM/(NCOMP_MON*(N_DOP+1)))*(NCOMP_MON*N_DOP-3)/3)
  PRINT*, ""
  PRINT*, "(g-)(g-)(g+) = ",REAL(NGPG(1))/REAL((N_DATA+1)*(N_ATOM/(NCOMP_MON*(N_DOP+1)))*(NCOMP_MON*N_DOP-3)/3)
  PRINT*, "(g-)(g+)(g+) = ",REAL(NGPG(2))/REAL((N_DATA+1)*(N_ATOM/(NCOMP_MON*(N_DOP+1)))*(NCOMP_MON*N_DOP-3)/3)
  PRINT*, "(g-)(t )(g+) = ",REAL(NGPG(3))/REAL((N_DATA+1)*(N_ATOM/(NCOMP_MON*(N_DOP+1)))*(NCOMP_MON*N_DOP-3)/3)
  PRINT*, ""
  PRINT*, "(g-)(g-)(t ) = ",REAL(NGTR(1))/REAL((N_DATA+1)*(N_ATOM/(NCOMP_MON*(N_DOP+1)))*(NCOMP_MON*N_DOP-3)/3)
  PRINT*, "(g-)(g+)(t ) = ",REAL(NGTR(2))/REAL((N_DATA+1)*(N_ATOM/(NCOMP_MON*(N_DOP+1)))*(NCOMP_MON*N_DOP-3)/3)
  PRINT*, "(g-)(t )(t ) = ",REAL(NGTR(3))/REAL((N_DATA+1)*(N_ATOM/(NCOMP_MON*(N_DOP+1)))*(NCOMP_MON*N_DOP-3)/3)
  PRINT*, ""
  PRINT*, "(g+)(g-)(t ) = ",REAL(PGTR(1))/REAL((N_DATA+1)*(N_ATOM/(NCOMP_MON*(N_DOP+1)))*(NCOMP_MON*N_DOP-3)/3)
  PRINT*, "(g+)(g+)(t ) = ",REAL(PGTR(2))/REAL((N_DATA+1)*(N_ATOM/(NCOMP_MON*(N_DOP+1)))*(NCOMP_MON*N_DOP-3)/3)
  PRINT*, "(g+)(t )(t ) = ",REAL(PGTR(3))/REAL((N_DATA+1)*(N_ATOM/(NCOMP_MON*(N_DOP+1)))*(NCOMP_MON*N_DOP-3)/3)
  PRINT*, "END PROGRAM"


! Construct Profile of conformers in PEO
!  ALLOCATE(CONFORM(1:N_ATOM-3,0:N_DATA),stat=IERR)
!  IF (IERR /= 0) STOP "CONFORM: problem in attempt to allocate memory"
!  DO I=0,N_DATA
!    DO K=1,N_ATOM-3
!      IF(DIH(K,I)*180.0D0/PI >= 150.0D0) THEN
!        CONFORM(K,I) = 1 ! trans
!      ENDIF
!      IF(DIH(K,I)*180.0D0/PI <= -150.0D0) THEN
!          CONFORM(K,I) = 1 ! trans
!      EndIF
!      IF((DIH(K,I)*180.0D0/PI >= 30.0D0) .and. (DIH(K,I)*180.0D0/PI <= 90.0D0)) THEN
!        CONFORM(K,I) = 1/2 ! gaushe
!      ENDIF
!      IF((DIH(K,I)*180.0D0/PI >= -90.0D0) .and. (DIH(K,I)*180.0D0/PI <= -30.0D0)) THEN
!        CONFORM(K,I) = -1/2 ! -gaushe
!      ENDIF
!      IF(CONFORM(K,I) == 0) THEN
!        PRINT*, "STRANGE Site",I,"th line and ",K,"th bond = ",DIH(K,I)*180.0D0/PI
!      ENDIF
!    ENDDO      
!  ENDDO

END PROGRAM conformer
