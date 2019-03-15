PROGRAM SQ

  IMPLICIT NONE
! command
  INTEGER :: NUM_OF_ARGS, I_ARG
  CHARACTER(LEN=20) :: FILE_NAME, INPUT_NAME, OUTPUT_NAME
  LOGICAL :: LookForInput=.FALSE.,LookForOutput=.FALSE.
  LOGICAL :: LookForStart=.FALSE.,LookForEnd=.FALSE.,LookForDq=.FALSE.
  LOGICAL :: FileExist
  DOUBLE PRECISION :: START_Q, END_Q, D_Q

! Read .gro file
!  INTEGER :: RES_ID, INDEX, MAX_DATA
!  CHARACTER(LEN=5) :: RES_NAME, ATOM_TYPE
!  REAL :: X, Y, Z, VX, VY, VZ
!  REAL, DIMENSION(1:3) :: BOX

! Read
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AR_DATA,TEMP
  INTEGER :: MAX_DATA = 5000, ICR = 2000
  REAL, DIMENSION(1:3) :: BOX
  REAL :: X, Y, Z

! calculate #1
  INTEGER :: N_DATA, I, J, K, IERR
  INTEGER, DIMENSION(1:3) :: MAX_Q
  DOUBLE PRECISION :: PI = 3.141592653589793
  DOUBLE PRECISION, DIMENSION(1:3) :: DIST
  COMPLEX, DIMENSION(1:3) :: IKR
  COMPLEX :: XI = (0.0D0,1.0D0)
  COMPLEX, DIMENSION(:), ALLOCATABLE :: SQ1

! Set default
  INPUT_NAME = "sq.xyz"
  OUTPUT_NAME = "sq.dat"
  START_Q = 0
  END_Q = 30
  D_Q = 10

  
! Command
  NUM_OF_ARGS = command_argument_count()

  IF (NUM_OF_ARGS == 0) THEN
    WRITE(*,*) "USAGE : ./sq.x -i [xtc file]] -o [sq.xvg] -b [real] -e [real] -dq [integer]"
    WRITE(*,*) "-i sq.xyz      input file (only contains (x,y,z) on coloumn 1, 2, and 3)"
    WRITE(*,*) "-o sq.dat      output file"
    WRITE(*,*) "-b 0           starting q (1/nm)"
    WRITE(*,*) "-e 30          ending q (1/nm)"
    WRITE(*,*) "-dq 10         number of q in [0:min_q]"
    WRITE(*,*) "It calcualtes structure factor of Rouse Spring Model."
    WRITE(*,*) "Assume that all particles are identical and box is cubic."
    WRITE(*,*) "EXAMPLE: tail -1 confout.gro > sq.xyz for get box size"
    WRITE(*,*) "EXAMPLE: awk '{ if ( $2 == ""C1"" || $2 == ""C2"" || $2 == ""O"" ) print $4,"" "",$5,"" "",$6 }' confout.gro >> sq.xyz"
    WRITE(*,*) "EXAMPLE: ./sq.x -i sq.xyz -dq 20"
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
          CASE("-o")
            LookForOutput=.TRUE. 
            CYCLE
          CASE("-b")
            LookForStart=.TRUE. 
            CYCLE
          CASE("-e")
            LookForEnd=.TRUE. 
            CYCLE
          CASE("-dq")
            LookForDq=.TRUE.
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
              ELSE IF(LookForOutput) THEN
                inquire(file=TRIM(FILE_NAME),exist=FileExist)
                IF(FileExist) THEN
                  WRITE(*,*) 'Output File ',FILE_NAME,' already exists'
                EndIF
                OUTPUT_NAME = FILE_NAME
                LookForOutput = .FALSE.
                CYCLE
              ELSE IF(LookForStart) THEN
                READ(FILE_NAME,*) START_Q
                LookForStart = .FALSE.
                CYCLE
              ELSE IF(LookForEnd) THEN
                READ(FILE_NAME,*) END_Q
                LookForEnd = .FALSE.
                CYCLE
              ELSE IF(LookForDq) THEN
                READ(FILE_NAME,*) D_Q
                LookForDq = .FALSE.
                CYCLE
              ELSE
                PRINT*, "Option ", FILE_NAME, "is unknown"
                STOP
            ENDIF
         END SELECT
         PRINT*, "Your command line: sq.x -i ",TRIM(INPUT_NAME),&
          " -o ",TRIM(OUTPUT_NAME)," -b ",START_Q," -e ",END_Q
      ENDDO
  ENDIF

! Read Input file
  PRINT*, "========================"
  PRINT*, "Reading Input File, ",TRIM(INPUT_NAME)
  OPEN(UNIT=11,FILE=TRIM(INPUT_NAME),STATUS='OLD')

  ALLOCATE(AR_DATA(3,MAX_DATA),stat=IERR) 
  IF (IERR /= 0) THEN
    PRINT*, "DATA: problem in attempt to allocate memory"
    STOP
  ENDIF
  READ(11,*) (BOX(I),I=1,3)
  I = 0
  DO 
    READ(11,*,IOSTAT=IERR) X,Y,Z
    IF(IERR < 0) THEN
      PRINT*, "End of line is ",I+1
      EXIT
      ELSE IF(IERR > 0) THEN
        PRINT*, "# line is not correct. Should Check!"
        STOP
      ELSE ! if the array size has to increase
        I = I + 1
        IF(I > MAX_DATA) THEN
          ALLOCATE(TEMP(3,SIZE(AR_DATA)/3+ICR),stat=IERR)
          IF (IERR /= 0) THEN
            PRINT*, "reallocate_temp: problem in attempt to allocate memory"
            STOP
            ELSE
              TEMP = 0.0D0
              TEMP=AR_DATA
              DEALLOCATE(AR_DATA)
          ENDIF
          ALLOCATE(AR_DATA(3,SIZE(TEMP)/3+ICR),stat=IERR)
          IF (IERR /= 0) THEN
            PRINT*, "reallocate_array: problem in attempt to allocate memory"
            STOP
            ELSE
              AR_DATA=TEMP
              DEALLOCATE(TEMP)
          ENDIF
          MAX_DATA = MAX_DATA+ICR
        ENDIF
        AR_DATA(1,I) = X
        AR_DATA(2,I) = Y
        AR_DATA(3,I) = Z
    ENDIF
    IF(MOD(I+1,1000) == 0) PRINT*,"Reading ",I+1,"th line"
  ENDDO

! for reading gro file
!  READ(11,*)
!  READ(11,*) MAX_DATA
!  ALLOCATE(AR_DATA(3,MAX_DATA),stat=IERR) 
!  AR_DATA = 0.0D0
!  IF (IERR /= 0) THEN
!    PRINT*, "AR_DATA: problem in attempt to allocate memory"
!    STOP
!  EndIF
!  DO 
!    READ(11,"(i5,2a5,i5,3f8.3,3f8.4)",IOSTAT=IERR) RES_ID,RES_NAME,ATOM_TYPE,INDEX,X,Y,Z,VX,VY,VZ
!    IF(IERR < 0) THEN
!      PRINT*, "Check! number of atoms in line 2 and real ddata are not same, so end ahead"
!      EXIT
!      ELSE IF(IERR > 0) THEN
!        PRINT*, "# line is not correct. Should Check!"
!        STOP
!      ELSE ! if the array size has to increase
!        AR_DATA(1,I) = X
!        AR_DATA(2,I) = Y
!        AR_DATA(3,I) = Z
!    ENDIF
!    IF(MOD(I,4000) == 0) PRINT*,"Reading ",I,"th line"
!  ENDDO
!  PRINT*, "Loaded particle's number is ",I-1
!  READ(11,*) (BOX(I),I=1,3)
!  CLOSE(11)

! Calculate structure factor 
! ==============================
!  S(k) = (1/N) * (sum_over)(exp(ik(r_j-r_i))
!  S(k)=\frac { 1 }{ N } \sum _{ i,j=1 }^{ N }{ { e }^{ ik({ r }_{ i }-{ r }_{ j }) } } 
! ==============================

  PRINT*, "========================"
  PRINT*, "Calculating distances of pair atoms"
  IF((BOX(1) /= BOX(2)) .or. (BOX(1) /= BOX(3)) .or. (BOX(2) /= BOX(3))) THEN
    PRINT*, "ERROR: this is not isotopic system which is not supported by current program"
    STOP
  ENDIF
  ! for cubic system
  N_DATA = I
  MAX_Q(1) = INT(END_Q*(BOX(1)*D_Q)/PI/2.0D0)
  PRINT*, "Number of k vector = ",MAX_Q(1)
  ALLOCATE(SQ1(1:MAX_Q(1)),stat=IERR) 
  IF (IERR /= 0) THEN
    PRINT*, "SQ: problem in attempt to allocate memory"
    STOP
  EndIF
  SQ1 = (0.0D0,0.0D0)
  DO I=1,N_DATA
    IF(MOD(I,100) == 0) PRINT*,"Calculating neighbor distances from ",I,"th particle" 
    DO J=1,N_DATA
      DIST(:) = DABS(AR_DATA(:,I) - AR_DATA(:,J))
      DIST(:) = DIST(:)-BOX(:)*DINT(DIST(:)/BOX(:))
      DO K=1,MAX_Q(1)      
        IKR(:) = CDEXP(XI*DIST(:)*2*PI*K/(BOX(:)*D_Q))
        SQ1(K) = SQ1(K) + IKR(1)*IKR(2)*IKR(3)
      ENDDO
    ENDDO
  ENDDO
  
! write structure factor
  PRINT*, "========================"
  PRINT*, "writing structure factor"

  OPEN(UNIT=12,FILE=TRIM(OUTPUT_NAME),STATUS='unknown')
  WRITE(12,'(A19,1X,F13.6)') "# refer resolution ",2*PI/BOX(1)
  DO K=1,MAX_Q(1)      
!    WRITE(12,'(3(F13.6,1X),F13.6)') 2*PI*K/(BOX(:)*D_Q), CABS(SQ1(K))/REAL(N_DATA)
     WRITE(12,'(F13.6,1X,F13.6)') 2*PI*K/(BOX(1)*D_Q), CABS(SQ1(K))/REAL(N_DATA)
  ENDDO
  CLOSE(12)

  PRINT*, "Writing ",TRIM(OUTPUT_NAME)," COMPLETE"

END PROGRAM SQ
