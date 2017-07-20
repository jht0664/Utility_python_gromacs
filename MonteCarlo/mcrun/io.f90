SUBROUTINE newfile_del_oldfile(outputfile)
  CHARACTER(LEN=100) :: outputfile
  LOGICAL :: file_exist
  INQUIRE(file=outputfile, exist=file_exist)
  if (file_exist) THEN
    OPEN(UNIT=1,FILE=outputfile,status='old')
    close(1,status='delete')
  endif
END SUBROUTINE newfile_del_oldfile

! read initial configuration (save 'pos' block)
! input: filename
SUBROUTINE read_ic(filename)
  USE pos
  USE ints
  USE sigmas
  
  IMPLICIT NONE
  CHARACTER(LEN=*) :: filename
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
  ALLOCATE(x(nptot),y(nptot),z(nptot),typs(nptot),STAT=ierr)
  call print_over_memory(ierr)
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
    CHARACTER(LEN=*) :: filename

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
  use coupling_pres
  
  IMPLICIT NONE
  CHARACTER(LEN=*) :: filename
  CHARACTER(LEN=5) :: name
  INTEGER :: nseed
  DOUBLE PRECISION :: boxmin, SECNDS

  write(*,*) "read_inp:"
  OPEN(UNIT=1,FILE=filename,STATUS='OLD')
  ! openmp setting
  READ(1,*) openmp_thread_num ! openmp thread number
  write(*,*) ' #processors available = ', omp_get_num_procs()
  write(*,*) ' #threads = ', omp_get_max_threads()
  write(*,*) ' #threads you set = ', openmp_thread_num
  IF( openmp_thread_num .GT. omp_get_num_procs()) THEN
    openmp_thread_num = omp_get_num_procs()
    write(*,*) ' thread number is set to ', openmp_thread_num
  ENDIF
  call omp_set_num_threads ( openmp_thread_num )
  ! INITIALIZE THE SEED
  READ(1,*) nseed 
  write(*,*) " input seed = ", nseed
  if (nseed == 0) then
    nseed = 2*INT(SECNDS(0.0)) + 2937
    write(*,*) " random seed generated with ", nseed
  endif
  CALL SRAND(nseed)
  ! time
  READ(1,*) nncon, ncon, nskip ! number of trial moves, number of skip for print
  write(*,*) " total simulation steps = ", nncon, " x ", ncon, ", skip time = ", nskip
  READ(1,*)
  
  ! simulation setting  
  READ(1,*) name ! ensemble name (NVT, NPT, NPTX2)
  write(*,*) " ensemble = ", name
  call movetype_set_ensemble(name) ! turn on couplings depending on name
  call movetype_init ! call movetype_init(number of movement types)
  call movetype_read(1) !read detail setting of movement types
  READ(1,*)
  ! calculate g(r)?
  READ(1,*) igr, rdf_nstep, bsz ! bin size
  if(igr == 'YES') then
    boxmin = MIN(box(1),box(2),box(3))/2.0D0
    nbin = INT(boxmin/bsz)
    write(*,*) " activate g(r) calculation, rdf_nstep = ", rdf_nstep, &
               ", bin size = ", real(bsz), ", nbin = ", nbin
  endif
  ! printing configuration?
  READ(1,*) iconfig, traj_nstep ! number of configuration_save skip
  ! printing pressure?, delta_V/V = ratio_dv_v * nwindow_pressure 
  !  (for example, ratio_dv_v = -0.0025 and nwindow_pressure = 5) 
  READ(1,*) ipres, ratio_dv_v, n_dv
  if (ipres == 'YES') THEN
    if (ensemble_pres == .true.) then
      write(*,*) "print pressure only works in NVT"
      ipres = 'NO'
    else
      write(*,*) " activate pressure calculation, ratio_dv_v = ", real(ratio_dv_v), &
                 ", n_dv = ", n_dv
    endif
  endif
  ! calculate the ratio of Henry constant H_i to fugacity f_j by identity exchange
  ! argument: Yes/no, atomname of solvent, solute, skip time, number of trial exchanges
  READ(1,*) iex, ex_solvent, ex_solute, ex_nstep, ex_ntry 
  if ((ensemble_exch == .true.) .and. iex == 'YES') THEN
    write(*,*) "Hi/fj calculation does not work in NPTX$"
    iex = 'NO'
  else if ((ensemble_pres == .true.) .and. iex == 'YES') THEN
    write(*,*) "Hi/fj calculation does not support NPT yet"
    iex = 'NO'
  endif
  if (iex == 'YES') THEN
    write(*,*) " activated Hi/fj calculation with ex_solvent ", ex_solvent, &
               ", ex_solute ", ex_solute, ", ex_nstep", ex_nstep, ", ex_ntry ", ex_ntry
    write(*,*) " Keep in mind that it only works for the system which follows Rault's law"
  endif
  CLOSE(1)
  RETURN
END SUBROUTINE read_inp

! save coordinate in gromac format
SUBROUTINE save_gro(filename,action,istep)
  use pos
  use ints
  IMPLICIT NONE
  CHARACTER(LEN=*) :: filename, action

  DOUBLE PRECISION :: xi, yi, zi
  INTEGER :: IJ, istep

!  write(*,*) "save_gro:"
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
SUBROUTINE save_xtc(xtcf,istep)
  USE xdr, only: xtcfile
  use pos
  use ints, only: nptot
  IMPLICIT NONE
  type(xtcfile) :: xtcf
  REAL, DIMENSION(3,nptot) :: coord
  REAL, DIMENSION(3,3) :: boxl
  integer :: i, istep

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
