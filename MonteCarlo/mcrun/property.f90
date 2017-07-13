
SUBROUTINE pressure_dV_init(outputfile)
  CHARACTER(LEN=100) :: outputfile
  LOGICAL :: file_exist
  INQUIRE(file=outputfile, exist=file_exist)
  IF(file_exist) THEN
    OPEN(UNIT=1,FILE=outputfile,status='old')
    close(1,status='delete')
  ENDIF
END SUBROUTINE pressure_dV_init

! pressure calculation using volume change in NVT ensemble
SUBROUTINE pressure_dV(outputfile)
  use calc_pres
  use ints, only: nptot
  use pos
 
  IMPLICIT NONE
  CHARACTER(LEN=100) :: outputfile
  INTEGER :: i, ierr
  LOGICAL :: file_exist, overlap_q
  double precision :: vol

!  write(*,*) "pressure_dV:"
  ALLOCATE(pressure(n_dv),STAT=ierr) 
  ALLOCATE(scale_length(n_dv),STAT=ierr)
  ALLOCATE(noverlap(n_dv),STAT=ierr)
  pressure = 0.0D0
  scale_length = 0.0D0
  noverlap = 0
  ! compute scale_length, delta_l
  DO i = 1, n_dv
    if(1+ratio_dv_v*i < 0) THEN
      write(*,*) "wrong input ratio_dv_v (big number)"
      stop
    endif
    scale_length(i) = (1+ratio_dv_v*i)**(1.0D0/3.0D0)
!    write(*,*) "scale_length ", scale_length(i)
  ENDDO
!  write(*,*) "1-1"
  ! count number of overlap paris
  DO i = 1, nptot
    call overlap(i,x(i),y(i),z(i),typs(i),'c',overlap_q) !count overlap pairs
  ENDDO
!  write(*,*) "2"
  ! calculate pressure of WR model
  !  beta*Pressure = limit(del_V -> 0) <number density - (number of ovelapping A and B pair) / del_V >_ensemble average
  !   to remove duplicates by dividing by 2
  vol = box(1)*box(2)*box(3)
  DO i = 1, n_dv
    !WRITE(*,*) DBLE(noverlap(i))/(ratio_dv_v*i)/2.0D0
    pressure(i) = DBLE(nptot)/(box(1)*box(2)*box(3)) - DBLE(noverlap(i))/(ratio_dv_v*i*vol)/2.0D0
  ENDDO
!  write(*,*) "3"
  ! save data
  INQUIRE(file=outputfile, exist=file_exist)
  IF(file_exist) THEN
    OPEN(UNIT=1,FILE=outputfile,status='old', position='append', action='write')
  ELSE
    OPEN(UNIT=1,FILE=outputfile,status='new', action='write')  
    WRITE(1,*) "# first line is delta_V's"
    WRITE(1,'(*(E16.8,2x))') (ratio_dv_v*i, i=1, n_dv)
  ENDIF
  WRITE(1,'(*(E16.8,2x))') ( pressure(i), i=1,n_dv )
  WRITE(*,'(*(E16.8,2x))') ( pressure(i), i=1,n_dv )
  CLOSE(1)
!  write(*,*) "4"
  DEALLOCATE(pressure,stat=ierr)
  DEALLOCATE(noverlap,STAT=ierr)
  DEALLOCATE(scale_length,STAT=ierr)
  RETURN
END SUBROUTINE pressure_dV

SUBROUTINE save_log(filename)
  use pos
  use ints
  use sigmas
  use inp
  use trans
  use rdf
  use movetype
  IMPLICIT NONE
  CHARACTER(LEN=100) :: filename
  DOUBLE PRECISION :: VOL, PFC, RGAV, R2AV

    rgav = 0
    r2av = 0
  !    IF(NAVER.GT.0)THEN
!      RGAV = RGAV/DBLE(NAVER*nmol_a)
!      R2AV = R2AV/DBLE(NAVER*nmol_a)
!    ENDIF

    OPEN(UNIT=1,FILE=filename,STATUS='UNKNOWN')
  ! PACKING FRACTION
    VOL = box(1)*box(2)*box(3)
    PFC = (nmol_a*nmon_a*sigma_a**3+nmol_b*nmon_b*sigma_b**3)*PI/(6.*VOL)
    WRITE(1,*)"PFC,RGAV,R2AV"
    WRITE(1,117) PFC,RGAV,R2AV
117 FORMAT(1X,F6.4,2X,5(F9.6,2X),5(F12.6,2X) )

!   result
    WRITE(1,*)'RESULTS OF FREE ROTATED CHAIN POLYMER ', ensemble_name, ' RUN'
    WRITE(1,*)
    WRITE(1,111) nmol_a,nmol_b,sigma_b,sigma_a,box(1),box(2),box(3)
    WRITE (1,114) NCON, NSKIP
    WRITE (1,115) dlr_a, dlr_b,BSZ,NBIN
    WRITE (1,116) movetype_ntry(1),0, &
                  movetype_nsuccess(1),0, &
                  DBLE(movetype_ntry(1)/movetype_nsuccess(1)), &
                  0
    close(1)

111 FORMAT(2X, I10, 6X, I10, 6X, 'NUMBER OF A'  / &
           2X, I10, 6X, I10, 6X, 'NUMBER OF B'  / & 
           2X,F6.4, 6X,'DIAMETER OF A' / &
           2X,F6.4, 6X,'DAIMETER OF B' / &       
           4X, E25.18, 2X, &
           4X, E25.18, 2X, &
           4X, E25.18, 2X)
114 FORMAT (1X / 1X, 'NO. TOTAL MOVES           ', 2X, I15 / &
                 1X, 'NO. MOVES PER ACCUMULATION', 2X, I15)
115 FORMAT (1X / 1X, 'DELTA CHAIN   ', 6X, F5.3 / &    
                 1X, 'DELTA MEDIA   ', 4X, F7.4 / &
                 1X, 'BIN SIZES     ', 6X, F5.3 / &
                 1X, 'NUMBER OF BINS', 7X, I4)
116 FORMAT (1X / 13X, 'MOVE STATISTICS' / &
     1X,       14X,'TRANS',6X,'GDI' /&
     1X,' ATTEMPTED',1X,I8,1X,I8 /&
     1X,'SUCCESSFUL',1X,I8,1X,I8 /&
     1X,'  FRACTION',3X,F6.4,3X,F6.4)
END SUBROUTINE save_log