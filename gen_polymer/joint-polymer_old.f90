PROGRAM JOINT_PEO

  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793D0
  INTEGER, PARAMETER :: NBMAX = 10000
  CHARACTER(LEN=5) :: POL_N1, POL_N2
  DOUBLE PRECISION :: ROTA, RAND
  INTEGER :: I, J

  CHARACTER(LEN=5) :: C1, C2
  CHARACTER(LEN=5),DIMENSION(NBMAX) :: CC1, CC2, CC1_NEW, CC2_NEW
  CHARACTER(LEN=7) :: F1, F2

  REAL :: UNIT1, UNIT2, UNIT3, UNITV, TOL, SELE_ROT
  REAL :: D1, D2, D3, MIN_BOX1, MIN_BOX2, MIN_BOX3, MAX_BOX1, MAX_BOX2, MAX_BOX3
  REAL :: ORIGINX, ORIGINY, ORIGINZ, ROT1, ROT2, ROT3
  REAL, DIMENSION(3,3) :: ROTM
  REAL, DIMENSION(NBMAX) :: DD1, DD2, DD3, DD1_NEW, DD2_NEW, DD3_NEW

  INTEGER :: A5, B5, STATUS, NATOM1, NATOM2, CHECK_D, NSEED, TRY_ROT
  INTEGER, DIMENSION(NBMAX) :: A5_NEW, B5_NEW

  OPEN(Unit=10,File='conf1.gro',STATUS='OLD')
  OPEN(Unit=11,file='conf2.gro',status='old')
  OPEN(Unit=12,File='initial.gro',Status='NEW')

  PRINT*, "This program makes a chain of Polyethyl oxide."
  PRINT*, "2 chians would be jointed."
  PRINT*, "conf1.gro & conf2.gro - input file"
  PRINT*, 'Must be removed solvents in input files'

  PRINT*, "Minimum distance? (nm), prefered value = 0.15 nm"
  READ(*,*) TOL

  READ(10,*) F1
  READ(10,*) NATOM1
  READ(11,*) F1
  READ(11,*) NATOM2

  WRITE(12,"(A5)") F1
  WRITE(12,"(I5)") NATOM1+NATOM2-9

      NSEED=2*INT(SECNDS(0.0)) + 2937
      CALL SRAND(NSEED)

  DO
    ! A5 (Group Index) C1 (Group Name) C2 (Atom Name) B5 (Atom Index) D1-D3 (X,Y,Z [nm])
    READ(10,905,IOSTAT=STATUS) A5, C1, C2, B5, D1, D2, D3

    IF(STATUS < 0) THEN
      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
      PRINT*, "Abnormal Finish! Check your code and files!"
      STOP
      ELSE IF(STATUS > 0) THEN
        PRINT*, "WHY ERROR?"
        PRINT*, "the total numbers of atoms in conf1.gro and conf2.gro would be not proper."
        PRINT*, "Or your file format is not from trjconv_d."
        PRINT*, "This program is only for one step from trjconv_d"
        PRINT*, "CHECK!!!!"
        CYCLE
    ENDIF

    DD1(B5) = D1
    DD2(B5) = D2
    DD3(B5) = D3
    CC1(B5) = C1
    CC2(B5) = C2

    ! write 1st chain
    IF(B5 < NATOM1-7) THEN
      IF(C1 == "PEOE ") THEN
        C1 = "PEO  "
      ENDIF
      WRITE(12,905) A5, C1, C2, B5, D1, D2, D3
      CYCLE
    ENDIF

    ! for connecting 2st chain
    IF(B5 == NATOM1-7) THEN
      ! the position of all atoms of 2nd chain is represented by vector where origin is C2 of head gorup of 2nd chain
      DO J=1, NATOM2
        READ(11,905,IOSTAT=STATUS) A5_NEW(J), CC1_NEW(J), CC2_NEW(J), B5_NEW(J), DD1_NEW(J), DD2_NEW(J), DD3_NEW(J)
        IF(J==1) THEN
          ORIGINX = DD1_NEW(1)
          ORIGINY = DD2_NEW(1)
          ORIGINZ = DD3_NEW(1)
        ENDIF
        DD1_NEW(J) = DD1_NEW(J) - ORIGINX
        DD2_NEW(J) = DD2_NEW(J) - ORIGINY
        DD3_NEW(J) = DD3_NEW(J) - ORIGINZ
      ENDDO

      ! All atoms of 2nd chain are tranformed into the new coordinate system where orgin is C2 of head group of 2nd chain
      DO J=1, NATOM2
        DD1_NEW(J) = DD1(NATOM1-7) + DD1_NEW(J)
        DD2_NEW(J) = DD2(NATOM1-7) + DD2_NEW(J)
        DD3_NEW(J) = DD3(NATOM1-7) + DD3_NEW(J)
      ENDDO

      TRY_ROT = 0



      DO
        ! check the distance
        CALL DISTANCE(TOL,NATOM1,NATOM2,DD1,DD2,DD3,DD1_NEW,DD2_NEW,DD3_NEW,CHECK_D)

        IF(CHECK_D == 0) THEN
          ! If okay, we write the new coordinate.
          PRINT*, "Trial number of rotation :", TRY_ROT
          DO J=1, NATOM2
            IF(CC1_NEW(J) == "PEOB ") THEN
              CC1_NEW(J) = "PEO  "
            ENDIF
            IF(J < 4) THEN
              WRITE(12,905) A5+A5_NEW(J)-1, CC1_NEW(J), CC2_NEW(J), NATOM1-8+J, DD1_NEW(J), DD2_NEW(J), DD3_NEW(J)
              ELSE IF(J > 4) THEN
                WRITE(12,905) A5+A5_NEW(J)-1, CC1_NEW(J), CC2_NEW(J), NATOM1-9+J, DD1_NEW(J), DD2_NEW(J), DD3_NEW(J)
            ENDIF
          ENDDO

          !Calculate a new box size
          MAX_BOX1 = DD1(1)
          MIN_BOX1 = DD1(1)
          MAX_BOX2 = DD2(1)
          MIN_BOX2 = DD2(1)
          MAX_BOX3 = DD3(1)
          MIN_BOX3 = DD3(1)

          DO I=1, NATOM1-7
            IF(DD1(I) < MIN_BOX1) THEN
              MIN_BOX1 = DD1(I)
            ENDIF
            IF(DD1(I) > MAX_BOX1) THEN
              MAX_BOX1 = DD1(I)
            ENDIF
            IF(DD2(I) < MIN_BOX2) THEN
              MIN_BOX2 = DD2(I)
            ENDIF
            IF(DD2(I) > MAX_BOX2) THEN
              MAX_BOX2 = DD2(I)
            ENDIF
            IF(DD3(I) < MIN_BOX3) THEN
              MIN_BOX3 = DD3(I)
            ENDIF
            IF(DD3(I) > MAX_BOX3) THEN
              MAX_BOX3 = DD3(I)
            ENDIF
          ENDDO

          DO J=1, NATOM2
            IF(DD1_NEW(J) < MIN_BOX1) THEN
              MIN_BOX1 = DD1_NEW(J)
            ENDIF
            IF(DD1_NEW(J) > MAX_BOX1) THEN
              MAX_BOX1 = DD1_NEW(J)
            ENDIF
            IF(DD2_NEW(J) < MIN_BOX2) THEN
              MIN_BOX2 = DD2_NEW(J)
            ENDIF
            IF(DD2_NEW(J) > MAX_BOX2) THEN
              MAX_BOX2 = DD2_NEW(J)
            ENDIF
            IF(DD3_NEW(J) < MIN_BOX3) THEN
              MIN_BOX3 = DD3_NEW(J)
            ENDIF
            IF(DD3_NEW(J) > MAX_BOX3) THEN
              MAX_BOX3 = DD3_NEW(J)
            ENDIF
          ENDDO

          WRITE(12,904) MAX_BOX1-MIN_BOX1+0.1D0, MAX_BOX2-MIN_BOX2+0.1D0, MAX_BOX3-MIN_BOX3+0.1D0
          CLOSE(10)
          CLOSE(11)
          CLOSE(12)
          PRINT*, "Finish! Output file is initial.gro."
          STOP
        ENDIF

        IF(CHECK_D == 1) THEN
          TRY_ROT = TRY_ROT + 1

          ! GENERATE RANDOM NUMBER SEED
          ROTA = RAND()*2*PI
          PRINT*, "Try angle (deg.) = ", INT(ROTA*180/PI)

          ! axis unit vector
          SELE_ROT = RAND()
          ! axis unit vector is C1 of A polymer and C2 of B polymer (actually H13 of A polymer)
          IF(RAND() >= 0.5D0) THEN
            UNIT1 = DD1(NATOM1-10) - DD1(NATOM1-7)
            UNIT2 = DD2(NATOM1-10) - DD2(NATOM1-7)
            UNIT3 = DD3(NATOM1-10) - DD3(NATOM1-7)
            ! axis unit vector is C2-O of B polymer
            ELSE
              UNIT1 = DD1_NEW(5) - DD1(NATOM1-7)
              UNIT2 = DD2_NEW(5) - DD2(NATOM1-7)
              UNIT3 = DD3_NEW(5) - DD3(NATOM1-7)
          ENDIF
          UNITV = SQRT(UNIT1**2+UNIT2**2+UNIT3**2)
          UNIT1 = UNIT1/UNITV
          UNIT2 = UNIT2/UNITV
          UNIT3 = UNIT3/UNITV

          ! rotation matrix by an angle of ROT_ANGLE about an axis in the direction of axis unit vector
          ROTM(1,1) = REAL(COS(ROTA)+UNIT1**2*(1-COS(ROTA)))
          ROTM(1,2) = REAL(UNIT1*UNIT2*(1-COS(ROTA))-UNIT3*SIN(ROTA))
          ROTM(1,3) = REAL(UNIT1*UNIT3*(1-COS(ROTA))+UNIT2*SIN(ROTA))
          ROTM(2,1) = REAL(UNIT2*UNIT1*(1-COS(ROTA))+UNIT3*SIN(ROTA))
          ROTM(2,2) = REAL(COS(ROTA)+(1-COS(ROTA))*UNIT2**2)
          ROTM(2,3) = REAL(UNIT2*UNIT3*(1-COS(ROTA))-UNIT1*SIN(ROTA))
          ROTM(3,1) = REAL(UNIT3*UNIT1*(1-COS(ROTA))-UNIT2*SIN(ROTA))
          ROTM(3,2) = REAL(UNIT3*UNIT2*(1-COS(ROTA))+UNIT1*SIN(ROTA))
          ROTM(3,3) = REAL(COS(ROTA)+(1-COS(ROTA))*UNIT3**2)

          ! translation and rotation
          DO J=1,NATOM2
            DD1_NEW(J) = DD1_NEW(J) - DD1(NATOM1-7)
            DD2_NEW(J) = DD2_NEW(J) - DD2(NATOM1-7)
            DD3_NEW(J) = DD3_NEW(J) - DD3(NATOM1-7)
            ROT1 = ROTM(1,1)*DD1_NEW(J)+ROTM(1,2)*DD2_NEW(J)+ROTM(1,3)*DD3_NEW(J)
            ROT2 = ROTM(2,1)*DD1_NEW(J)+ROTM(2,2)*DD2_NEW(J)+ROTM(2,3)*DD3_NEW(J)
            ROT3 = ROTM(3,1)*DD1_NEW(J)+ROTM(3,2)*DD2_NEW(J)+ROTM(3,3)*DD3_NEW(J)
            DD1_NEW(J) = ROT1 + DD1(NATOM1-7)
            DD2_NEW(J) = ROT2 + DD2(NATOM1-7)
            DD3_NEW(J) = ROT3 + DD3(NATOM1-7)
          ENDDO
        ENDIF
      ENDDO
    ENDIF

    PRINT*, "ERROR!!1"

  ENDDO

! READ FORMAT
901 FORMAT(i5,2a5,i5,3f8.3,3f8.4)
! WRITE FORMAT for PEOB
902 FORMAT(i5,2a5,i5,a24)
! WRITE FORMAT for PEO
903 FORMAT(i5,2a5,i5,f8.3,a16)
904 FORMAT(3f8.3)
! READ from trajconv_d
905 FORMAT(i5,2a5,i5,3f9.4)

END PROGRAM

SUBROUTINE DISTANCE(TOL,NATOM1,NATOM2,DD1,DD2,DD3,DD1_NEW,DD2_NEW,DD3_NEW,CHECK_D)

  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793D0
  INTEGER, PARAMETER :: NBMAX = 10000
  REAL :: ID0, ID_SAVE, TOL
  INTEGER :: I, J
  REAL, DIMENSION(NBMAX) :: DD1, DD2, DD3, DD1_NEW, DD2_NEW, DD3_NEW
  INTEGER :: NATOM1, NATOM2, CHECK_D
  INTEGER, DIMENSION(NBMAX) :: A5_NEW, B5_NEW

  ID0 = 0.0D0
  ID_SAVE = 100
  CHECK_D = 0

  ! check the distance between 1st chain and 2nd chain
  DO I=1, NATOM1-8
    DO J=2, NATOM2
      ! except H23 of head group of 2nd chain
      IF(J == 4) THEN
        CYCLE
      ENDIF

      ! calculate distance of atoms
      ID0 = SQRT((DD1_NEW(J)-DD1(I))**2+(DD2_NEW(J)-DD2(I))**2+(DD3_NEW(J)-DD3(I))**2)

      ! save minimum distance
      IF(ID_SAVE > ID0) THEN
        ID_SAVE = ID0
      ENDIF

      IF(ID0 <= TOL) THEN
        PRINT*, "WARNING! Due to close distance, we have to rotate 2nd chain!"
        PRINT*, "Distance between ", I, "th atom in one chain and", J, "th atoms in other chain"
        PRINT*, "a distance(nm) = ", ID0, "tolerence = ", TOL
        CHECK_D = 1
        EXIT
        ELSE
          CYCLE
      ENDIF
      PRINT*, "Error!2"
      STOP
    ENDDO
    IF(CHECK_D == 1) THEN
      EXIT
    ENDIF
  ENDDO
  PRINT*, CHECK_D, ";If 0, Successful! Or not, Fail!!"
  PRINT*, "Minimum distance(nm) = ", ID_SAVE
  RETURN
END