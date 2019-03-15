PROGRAM SQ

  IMPLICIT NONE
! command
  INTEGER :: NUM_OF_ARGS, I_ARG
  CHARACTER(LEN=20) :: FILE_NAME, INPUT_NAME, OUTPUT_NAME
  LOGICAL :: LookForInput=.FALSE.,LookForOutput=.FALSE.
  LOGICAL :: LookForStart=.FALSE.,LookForEnd=.FALSE.
  LOGICAL :: FileExist
  DOUBLE PRECISION :: START_Q, END_Q

! Read .gro file
  INTEGER :: RES_ID, INDEX, MAX_DATA
  CHARACTER(LEN=5) :: RES_NAME, ATOM_TYPE
  REAL :: X, Y, Z, VX, VY, VZ
  REAL, DIMENSION(1:3) :: BOX

! practice
  DOUBLE PRECISION, DIMENSION(1:1818) :: DATA1, DATA2
  DOUBLE PRECISION :: TEMP1, TEMP2, TEMP3
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: SQ2

! Read
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AR_DATA 
  INTEGER :: N_DATA, I, J, K, IERR
  INTEGER, DIMENSION(1:3) :: MAX_Q
  DOUBLE PRECISION :: PI = 3.141592653589793
  DOUBLE PRECISION, DIMENSION(1:3) :: DIST
  COMPLEX, DIMENSION(1:3) :: IKR
  COMPLEX :: XI = (0.0D0,1.0D0)
  COMPLEX, DIMENSION(:), ALLOCATABLE :: SQ1


! Set default
  INPUT_NAME = "confout.gro"
  OUTPUT_NAME = "sq.xvg"
  START_Q = 0
  END_Q = 60
  
! Command
  NUM_OF_ARGS = command_argument_count()

  IF (NUM_OF_ARGS == 0) THEN
    WRITE(*,*) "USAGE : ./sq.x -i [xtc file]] -o [sq.xvg] -b [real] -e [real]"
    WRITE(*,*) "-i confout.gro input file"
    WRITE(*,*) "-o sq.xvg      output file"
    WRITE(*,*) "-b 0           Starting q (1/nm)"
    WRITE(*,*) "-e 60          Ending q (1/nm)"
    WRITE(*,*) "It calcualtes general structure factor"
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

  READ(11,*)
  READ(11,*) MAX_DATA
  ALLOCATE(AR_DATA(3,MAX_DATA),stat=IERR) 
  AR_DATA = 0.0D0
  IF (IERR /= 0) THEN
    PRINT*, "AR_DATA: problem in attempt to allocate memory"
    STOP
  EndIF
  DO I=1,MAX_DATA
    READ(11,"(i5,2a5,i5,3f8.3,3f8.4)",IOSTAT=IERR) RES_ID,RES_NAME,ATOM_TYPE,INDEX,X,Y,Z,VX,VY,VZ
    IF(IERR < 0) THEN
      PRINT*, "Check! number of atoms is not correct, so end ahead"
      EXIT
      ELSE IF(IERR > 0) THEN
        PRINT*, "# line is not correct. Should Check!"
        STOP
      ELSE ! if the array size has to increase
        AR_DATA(1,I) = X
        AR_DATA(2,I) = Y
        AR_DATA(3,I) = Z
    ENDIF
    IF(MOD(I,4000) == 0) PRINT*,"Reading ",I,"th line"
  ENDDO
  IF(I /= INDEX) PRINT*, "Warning: index is not the same as array index"
  READ(11,*) (BOX(I),I=1,3)
  CLOSE(11)

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
  MAX_Q(:) = INT(END_Q*BOX(:)/PI/2.0D0)
  ! for isotopic system
  ALLOCATE(SQ1(0:MAX_Q(1)),stat=IERR) 
  IF (IERR /= 0) THEN
    PRINT*, "SQ: problem in attempt to allocate memory"
    STOP
  EndIF
  SQ1 = (0.0D0,0.0D0)
  ALLOCATE(SQ2(0:MAX_Q(1)),stat=IERR) 
  OPEN(20,file='rdf.xvg') 
  DO I=1,1818
    READ(20,*) DATA1(I), DATA2(I)
    DATA1(I) = DATA1(I) + 0.001D0
  ENDDO
  DO K=1,MAX_Q(1)
    TEMP1 = 0.0D0
    DO I=1,1818
      TEMP2=(DATA2(I)-1.0D0)*dsin(DATA1(I)*PI*DBLE(2*K)/BOX(1))*(DATA1(I)**2)*(DATA1(2)-DATA1(1))/(DATA1(I)*PI*DBLE(2*K)/BOX(1))
      PRINT*,"2",TEMP2
      TEMP1=TEMP1+TEMP2
      PRINT*,"1",TEMP1
    ENDDO
!    PRINT*, TEMP1
    SQ2(K) = TEMP1*4.0D0*PI*10/BOX(1)/BOX(2)/BOX(3)+1.0D0
  ENDDO
  close(20)
!  DO I=1,MAX_DATA
!    IF(MOD(I,100) == 0) PRINT*,"Calculating neighbor distance from ",I,"th atom" 
!    DO J=I+1,MAX_DATA
!      DIST(:) = DABS(AR_DATA(:,I) - AR_DATA(:,J))
!      DIST(:) = DIST(:)-BOX(:)*DINT(DIST(:)/BOX(:))
!      DO K=0,MAX_Q(1)      
!        IKR(:) = CDEXP(XI*DIST(:)*PI*DBLE(2*K)/BOX(:))
!        SQ1(K) = SQ1(K) + IKR(1)*IKR(2)*IKR(3)
!      ENDDO
!    ENDDO
!  ENDDO
  
!  for anisotropic system
!  DO I=0,MAX_Q(1)
!    DO J=0,MAX_Q(2)
!      DO K=0,MAX_Q(3)
!        BASIS_Q(1,L)=2*PI*I/BOX(1)
!        BASIS_Q(2,L)=2*PI*J/BOX(2)
!        BASIS_Q(3,L)=2*PI*K/BOX(3)
!        IF((I*2*PI/BOX(1))**2+(J*2*PI/BOX(2))**2+(K*2*PI/BOX(3))**2 >= END_Q) L = L + 1
!      ENDDO
!    ENDDOㅌ
!  ENDDO
!  IF(L == 1) "ERROR: too small END_Q"
!  SQ = (0.0D0,0.0D0)
!  DO I=1,MAX_DATA
!    DO J=I+1,MAX_DATA
!      DIST(:) = DABS(AR_DATA(:,I) - AR_DATA(:,J))
!      DIST(:) = DIST(:)-DINT(DIST(:)/BOX(:))*BOX(:)
!      DO K=1,L-1      
!        IKR(:) = CDEXP(XI*DIST(:)*BASIS_Q(:,K))
!        SQ(K) = SQ(K) + IKR(1)*IKR(2)*IKR(3)
!      ENDDO
!    ENDDO
!  ENDDO
!  DO K=1,L-1
!    SQ(K) = SQ(K)/DCOMPLEX(MAX_DATA,0.0D0)
!  ENDDO


! write structure factor
  PRINT*, "========================"
  PRINT*, "writing structure factor"

  OPEN(UNIT=12,FILE=TRIM(OUTPUT_NAME),STATUS='unknown')
!  WRITE(12,'(2F13.6)') 0.0D0, 0.0D0
  DO K=0,MAX_Q(1)      
    WRITE(12,'(2F13.6)') 2*PI*K/BOX(1), ABS(SQ2(K))**2/REAL(MAX_DATA)
!    WRITE(12,'(2F13.6)') 2*PI*K/BOX(1), CABS(SQ1(K))/REAL(MAX_DATA)
!    WRITE(*,'(2F13.6)') REAL(SQ1(K)), CABS(SQ1(K))
  ENDDO
  CLOSE(12)

  PRINT*, "Writing ",TRIM(OUTPUT_NAME)," COMPLETE"

END PROGRAM SQ

SUBROUTINE D1_ARRAY_INC(ICR,IN_ARRAY)
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: IN_ARRAY
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TEMP
  INTEGER :: SIZE_IN, IERR
  INTEGER, INTENT(IN) :: ICR

  SIZE_IN = SIZE(IN_ARRAY)
  ALLOCATE(TEMP(SIZE_IN+ICR),stat=IERR)
  IF (IERR /= 0) THEN
    PRINT*, "reallocate_temp: problem in attempt to allocate memory"
    STOP
    ELSE
      TEMP = 0.0D0
      TEMP(1:SIZE_IN)=IN_ARRAY(1:SIZE_IN)
      DEALLOCATE(IN_ARRAY)
  ENDIF
  ! copy temp to output
  ALLOCATE(IN_ARRAY(SIZE_IN+ICR),stat=IERR)
  IF (IERR /= 0) THEN
    PRINT*, "reallocate_array: problem in attempt to allocate memory"
    STOP
    ELSE
      IN_ARRAY=TEMP
      DEALLOCATE(TEMP)
  ENDIF
END SUBROUTINE

