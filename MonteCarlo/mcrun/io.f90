! read initial configuration (save 'pos' block)
! input: filename
SUBROUTINE read_ic(filename)
  USE pos
  USE ints
  USE sigmas
  
  IMPLICIT NONE
  CHARACTER(LEN=100) :: filename
  INTEGER :: i, id, jd, ierr

  WRITE(*,*) "read_ic:"
! read coordinates info
  OPEN(UNIT=2,FILE=filename,STATUS='OLD')
  READ(2,*) nmol_a, nmon_a ! NUMBER OF A PARTICLE
  READ(2,*) nmol_b, nmon_b ! NUMBER OF B PARTICLE
  READ(2,*) sigma_a ! SIZE OF THE A
  READ(2,*) sigma_b ! SIZE OF THE B
  READ(2,*) box(1),box(2),box(3) ! BOX LENGTH

  sigma_ab = (sigma_a+sigma_b)/2.0D0    
  nptot = nmol_a * nmon_a + nmol_b * nmon_b
  
  ALLOCATE(x(nptot),STAT=ierr)
  IF (ierr /= 0) STOP
  ALLOCATE(y(nptot),STAT=ierr)
  IF (ierr /= 0) STOP
  ALLOCATE(z(nptot),STAT=ierr)
  IF (ierr /= 0) STOP
  ALLOCATE(typs(nptot),STAT=ierr)
  IF (ierr /= 0) STOP

  DO i = 1, nptot
    READ(2,*) ID,JD,typs(i),X(I),Y(I),Z(I)
  ENDDO
  CLOSE(2) 
  RETURN
END SUBROUTINE read_ic

SUBROUTINE save_ic(filename)
    use pos
    use ints
    use sigmas
    IMPLICIT NONE
    integer :: i, ij
    DOUBLE Precision :: xcnt, ycnt, zcnt
    CHARACTER(LEN=100) :: filename

    OPEN(UNIT=2,FILE=filename,STATUS='UNKNOWN')
    WRITE(2,111)nmol_a,nmon_a,nmol_b,nmon_b,sigma_a,sigma_b,box(1),box(2),box(3)
    i = 1
    DO ij = 1, nptot
      !WRITE(2,112)I,J,X(IJ),Y(IJ),Z(IJ)
      XCNT = X(IJ) - box(1)*DNINT(X(IJ)/box(1))
      YCNT = Y(IJ) - box(2)*DNINT(Y(IJ)/box(2))
      ZCNT = Z(IJ) - box(3)*DNINT(Z(IJ)/box(3))
      WRITE(2,112)MOD(IJ,100000),I,TYPS(IJ),XCNT,YCNT,ZCNT
    ENDDO
    CLOSE(UNIT=2,STATUS='KEEP')

111 FORMAT(2X, I10, 6X, I10, 6X, 'NUMBER OF A'  / &
           2X, I10, 6X, I10, 6X, 'NUMBER OF B'  / & 
           2X,F6.4, 6X,'DIAMETER OF A' / &
           2X,F6.4, 6X,'DAIMETER OF B' / &       
           4X, E25.18, 2X, &
           4X, E25.18, 2X, &
           4X, E25.18, 2X)
112 FORMAT(1X,I5,1X,I5,1X,A4,2X,3E18.10)
END SUBROUTINE save_ic

! read simulation settings
! 
SUBROUTINE read_inp(filename)
  use omp_lib
  use pos, only: box
  use inp
  use traj
  use trans
  use rdf 
  use calc_pres
  use calc_ex
  use movetype
  
  IMPLICIT NONE
  CHARACTER(LEN=100) :: filename
  INTEGER :: openmp_thread_num, nseed, movetype_size
  INTEGER :: ierr, i
  DOUBLE PRECISION :: boxmin

  write(*,*) "read_inp:"
  ! read simulation settings
  OPEN(UNIT=1,FILE=filename,STATUS='OLD')
  
  ! openmp setting
  READ(1,*) openmp_thread_num ! openmp thread number
  write(*,*) ' #processors available = ', omp_get_num_procs ( )
  write(*,*) ' #threads = ', omp_get_max_threads ( )
  write(*,*) ' #threads you set = ', openmp_thread_num
  IF( openmp_thread_num .GT. omp_get_num_procs ( )) THEN
    openmp_thread_num = omp_get_num_procs ( )
    write(*,*) ' thread number is set to ', openmp_thread_num
  ENDIF
  call omp_set_num_threads ( openmp_thread_num )
  
  ! simulation setting  
  READ(1,*) ensemble_name ! ensemble (NVT, NPT, NPAT, NPTX)
  IF( ensemble_name == 'NVT') THEN
    nmovetypes = 1
    call movetype_init ! call movetype_init(number of movement types)
  ELSE 
    WRITE(*,*) " Not supported yet, ", ensemble_name, " ensemble."
    STOP
  ENDIF
 
  ! INITIALIZE THE SEED
  READ(1,*) nseed 
  !nseed = 2*INT(SECNDS(0.0)) + 2937
  CALL SRAND(nseed)
  
  ! trajectory
  READ(1,*) ncon, nskip ! number of trial moves, number of skip for print
  ! translation
  READ(1,*) dlr_a, dlr_b ! translation increment for components
  trans_prob = 1.0

  ! probability of movetype
  READ(1,*) (movetype_prob(i), i=1,nmovetypes) 
  
  ! calculate g(r)?
  READ(1,*) igr, bsz ! bin size
  boxmin = MIN(box(1),box(2),box(3))/2.0D0
  nbin = INT(boxmin/bsz)

  ! printing configuration?
  READ(1,*) iconfig, nconfig ! number of configuration_save skip
  
  ! printing pressure?, delta_V/V = ratio_dv_v * nwindow_pressure 
  !  (for example, ratio_dv_v = -0.0025 and nwindow_pressure = 5) 
  READ(1,*) ipres, ratio_dv_v, n_dv
!  write(*,*) "input ", ratio_dv_v, n_dv
  
  ! calculate the ratio of Henry constant H_i to fugacity f_j by identity exchange
  READ(1,*) iex, ex_ntry ! number of trial

  istep = 0
  CLOSE(1)
  RETURN
END SUBROUTINE read_inp

! save coordinate in gromac format
SUBROUTINE save_gro(filename,action)
  use pos
  use traj
  use ints
  IMPLICIT NONE
  CHARACTER(LEN=1) :: action
  CHARACTER(LEN=100) :: filename

  DOUBLE PRECISION :: xi, yi, zi
  INTEGER :: IJ

  write(*,*) "save_gro:"
  IF(action == 'o') THEN
    OPEN(UNIT=3,FILE=filename,STATUS='UNKNOWN',action='write') ! overwrite if already exist
    WRITE(3,*) '# initial coordinate'
  ELSE IF(action == 'a') THEN
    OPEN(UNIT=3,FILE=filename,STATUS='UNKNOWN',position='append',action='write') ! append if already exist
    WRITE(3,*) '# ', istep, ' frame'
  ELSE
    write(*,*) ' wrong action'
    STOP
  ENDIF
  
  WRITE(3,*) nptot ! total #particles
  DO IJ = 1, nptot
    !WRITE(2,112)I,J,X(IJ),Y(IJ),Z(IJ)
    xi = X(IJ) - box(1)*FLOOR(X(IJ)/box(1))
    yi = Y(IJ) - box(2)*FLOOR(Y(IJ)/box(2))
    zi = Z(IJ) - box(3)*FLOOR(Z(IJ)/box(3))
    IF(typs(IJ) .EQ. 'A') THEN
      WRITE(3,211)MOD(IJ,100000),'LJA',typs(IJ),MOD(IJ,100000),xi,yi,zi
    ELSE
      WRITE(3,211)MOD(IJ,100000),'LJB',typs(IJ),MOD(IJ,100000),xi,yi,zi
    ENDIF
  ENDDO
  ! BOX
  WRITE(3,212) box(1),box(2),box(3)
  CLOSE(3)
211 FORMAT(i5,2a5,i5,3f8.3,3f8.4)
212 FORMAT(3f10.5)
END SUBROUTINE save_gro

! save xtc trajectory file
SUBROUTINE save_xtc(xtcf)
  USE xdr, only: xtcfile
  use pos
  use ints, only: nptot
  use traj, only: istep
  IMPLICIT NONE
  type(xtcfile) :: xtcf
  REAL, DIMENSION(3,nptot) :: coord
  REAL, DIMENSION(3,3) :: boxl
  integer :: i

  DO I = 1, nptot
    COORD(1,I) = X(I) - box(1)*FlOOR(X(I)/box(1))
    COORD(2,I) = Y(I) - box(2)*FLOOR(Y(I)/box(2))
    COORD(3,I) = Z(I) - box(3)*FLOOR(Z(I)/box(3))
  ENDDO
  
  ! This is the same order as found in the GRO format
  !write(*,'(11f9.5)') box(1), box(2), box(3), 0, 0, 0, 0, 0, 0 
  boxl(1,1) = box(1)
  boxl(2,2) = box(2)
  boxl(3,3) = box(3)
  boxl(1,2) = 0.0
  boxl(1,3) = 0.0
  boxl(2,1) = 0.0
  boxl(2,3) = 0.0
  boxl(3,1) = 0.0
  boxl(3,2) = 0.0

  !WRITE(*,*) " ENTER ICONFIG"
  ! Just an example to show what was read in
  !write(*,'(a,f12.6,a,i0)') " Time (ns): ", 
  ! istep/1000000.0, "  Step: ", istep
  !write(*,'(a,f12.6,a,i0)') " Precision: ", 
  ! 1000.0, "  No. Atoms: ", NBTOT+ NMEDIA
  !write(*,'(3f9.3)') xtcf % pos
  call xtcf % write(Nptot, istep, istep/1000000.0, boxl, COORD, 1000.0)
END SUBROUTINE save_xtc

SUBROUTINE pos_exit
  use pos
  integer :: ierr
  deallocate(x,stat=ierr)
  deallocate(y,stat=ierr)
  deallocate(z,stat=ierr)
  deallocate(typs,stat=ierr)
END SUBROUTINE pos_exit
