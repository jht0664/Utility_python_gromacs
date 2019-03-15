PROGRAM Delete_IL

	IMPLICIT NONE
	CHARACTER(LEN=5) :: C1, C2
	CHARACTER(LEN=7) :: F1, F2
	REAL :: D1, D2, D3, E1, E2, E3
	INTEGER :: A5, B5, STATUS, A51, B51

	OPEN(Unit=10,File='conf.gro',STATUS='OLD')
	OPEN(Unit=11,FIle='conf-out.gro',STATUS='UNKNOWN')
	READ(10,"(a7)") F1
	READ(10,"(a7)") F2
	WRITE(11,"(a7)") F1
	WRITE(11,"(a7)") F2

	A51 = 0
	B51 = 0
	adshfjahjdklfhkalhsdfkjlhajklsdhfhkadskjfkljadshfjkhklja 

	DO
		READ(10,901,IOSTAT=STATUS) A5, C1, C2, B5, D1, D2, D3, E1, E2, E3
	
		IF(STATUS < 0) THEN
			CLOSE(10)
			CLOSE(11)
			PRINT*, "Finish!"
			EXIT
			ELSE IF(STATUS > 0) THEN
				PRINT*, "WHY?"
				CYCLE
		ENDIF

		IF(C1 /= 'BF4') THEN
			WRITE(11,901) A5, C1, C2, B5, D1, D2, D3, E1, E2, E3
			A51 = A5
			B51 = B5
			CYCLE
		END IF

		IF(C1 == 'BF4') THEN
		    IF (C2 == '    B') THEN
				A51 = A51 + 1
			ENDIF
			B51 = B51 + 1
			WRITE(11,901) A51, C1, C2, B51, D1, D2, D3, E1, E2, E3 
			CYCLE
		ENDIF
	ENDDO

901	FORMAT(i5,2a5,i5,3f8.3,3f8.4)

END PROGRAM
