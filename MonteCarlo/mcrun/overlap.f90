! check overlap
! input: itry is index of molecule at original position
!        xo, yo, zo, to is position and typs of trial molecule
! output: overlap_q is output of overlap; if ovelaps, return true. Otherwise, false
SUBROUTINE overlap(ITRY,xo,yo,zo,to,action,overlap_q)
  use pos
  use cellmap
  use try
  use ints, only: nptot
  use sigmas
  use calc_ex
  IMPLICIT NONE
  integer :: itry, icell, id, niter
  integer :: ierr, nc, jcell, im, nn
  DOUBLE PRECISION :: xo, yo, zo
  DOUBLE PRECISION :: x1, y1, z1, xt, yt, zt, rt, dcomp
  CHARACTER(LEN=4) :: to, tt
  CHARACTER(LEN=1) :: action
  logical :: overlap_q
  INTEGER, DIMENSION(:), ALLOCATABLE :: JNEAR ! save the index of all neighbor particles

  if ( action /= 'x' .and. action /= 'c') THEN
    write(*,*) "wrong action argument"
    stop
  endif
!  write(*,*) "overlap:"
  ALLOCATE(jnear(nptot),stat=ierr)
!  IF(ierr /= 0) THEN
!    write(*,*) "oversize"
!  ENDIF
  jnear = 0
!write(*,*) "oversize2"
  ! check index of cell of trials
  ! find particles in neighbor cells
  overlap_q = .false.
!  write(*,*) "oversize1-1"
  X1 = XO - box(1)*DNINT(XO/box(1)-0.50D0) 
  Y1 = YO - box(2)*DNINT(YO/box(2)-0.50D0) 
  Z1 = ZO - box(3)*DNINT(ZO/box(3)-0.50D0) 
!  write(*,*) "oversize1-2"
  ICELL = 1 + INT( X1*CELLX ) &
            + INT( Y1*CELLY ) *MCX &
            + INT( Z1*CELLZ ) *MCY*MCX
!  write(*,*) "oversize1-2"
  ID = (ICELL-1)*27
  NITER = 0
!$omp parallel private ( nc ) shared ( niter, jnear )
  !$omp do
  DO NC = 1, 27
    JCELL = MAP(ID+NC)
    IF(JCELL.EQ.0) CONTINUE
    IM = LEAD(JCELL)
    DO WHILE (IM.NE.0)
      ! skip the identical polymer when saving neighbor particles, but count duplicates
      IF(IM .NE. ITRY) THEN
        NITER = NITER + 1
        JNEAR(NITER) = IM
      ENDIF
      IM = LIST(IM)
    ENDDO
  ENDDO
  !$omp end do
!$omp end parallel
  DO NN = 1, NITER ! # possible overlap particles calculated in previous loop
    IM = JNEAR(NN) 
    XT = X(IM) - X1
    YT = Y(IM) - Y1
    ZT = Z(IM) - Z1
    TT = typs(IM)
    XT = XT - box(1)*DNINT(XT/box(1))
    YT = YT - box(2)*DNINT(YT/box(2))
    ZT = ZT - box(3)*DNINT(ZT/box(3))
    RT = XT*XT+YT*YT+ZT*ZT
    IF(TT == TO) THEN ! for same atom type
      cycle
    ELSE
      DCOMP = sigma_ab**2
    ENDIF
    if (action == 'x' .and. rt < dcomp) then ! to avoid overlapping
      overlap_q = .true.
      exit
    else if (action == 'c') then ! to count overlapping for pressure
      call count_overlap_pres(rt, dcomp)
    endif
  ENDDO
  DEALLOCATE(jnear)
  RETURN
END SUBROUTINE overlap

! count number of overlapping pairs (finally, double counting)
SUBROUTINE count_overlap_pres(rt, dcomp)
  use calc_pres
  IMPLICIT NONE
  double precision :: rt, dcomp
  integer :: k
  DO k = 1, n_dv
    IF(RT*scale_length(k)**2 .LT. DCOMP) THEN
      pres_noverlap(k) = pres_noverlap(k) + 1 
    ENDIF
  ENDDO
  return
END SUBROUTINE count_overlap_pres

! check overlap
! input: itry is index of molecule at original position
!        jtry is index of molecule at trial position
! output: overlap_q is output of overlap; if ovelaps, return true. Otherwise, false
SUBROUTINE overlap_pres(itry,overlap_q)
  use try
  use cellmap_try
  use ints, only: nptot
  use sigmas
  IMPLICIT NONE
  
  integer :: itry, icell, id, niter
  integer :: ierr, nc, jcell, im, nn
  DOUBLE PRECISION :: xo, yo, zo, x1, y1, z1, xt, yt, zt, rt, dcomp
  CHARACTER(LEN=4) :: to, tt
  logical :: overlap_q
  INTEGER, DIMENSION(:), ALLOCATABLE :: JNEAR ! save the index of all neighbor particles

!  write(*,*) "overlap_pres:"
  ALLOCATE(jnear(nptot),stat=ierr)
  jnear = 0
! check index of cell of trials
! find particles in neighbor cells
  overlap_q = .false.
  xo = xitr(itry)
  yo = yitr(itry)
  zo = zitr(itry)
  to = titr(itry)
  X1 = XO - boxitr(1)*DNINT(XO/boxitr(1)-0.50D0) 
  Y1 = YO - boxitr(2)*DNINT(YO/boxitr(2)-0.50D0) 
  Z1 = ZO - boxitr(3)*DNINT(ZO/boxitr(3)-0.50D0) 
  ICELL = 1 + INT( X1*CELLX_itr ) &
            + INT( Y1*CELLY_itr ) *MCX_itr &
            + INT( Z1*CELLZ_itr ) *MCY_itr*MCX_itr
  ID = (ICELL-1)*27
  NITER = 0
!$omp parallel private ( nc ) shared ( niter, jnear )
  !$omp do
  DO NC = 1, 27
    JCELL = MAP_itr(ID+NC)
    IF(JCELL.EQ.0) CONTINUE
    IM = LEAD_itr(JCELL)
    DO WHILE (IM.NE.0)
      ! skip the identical polymer when saving neighbor particles, but count duplicates
      IF(IM .NE. ITRY) THEN
        NITER = NITER + 1
        JNEAR(NITER) = IM
      ENDIF
      IM = LIST_itr(IM)
    ENDDO
  ENDDO
  !$omp end do
!$omp end parallel
  DO NN = 1, NITER ! # possible overlap particles calculated in previous loop
    IM = JNEAR(NN) 
    XT = Xitr(IM) - X1
    YT = Yitr(IM) - Y1
    ZT = Zitr(IM) - Z1
    TT = titr(IM)
    XT = XT - boxitr(1)*DNINT(XT/boxitr(1))
    YT = YT - boxitr(2)*DNINT(YT/boxitr(2))
    ZT = ZT - boxitr(3)*DNINT(ZT/boxitr(3))
    RT = XT*XT+YT*YT+ZT*ZT
    IF(TT == TO) THEN ! for same atom type
      cycle
    ELSE
      DCOMP = sigma_ab**2
    ENDIF
    if (rt < dcomp) then ! to avoid overlapping
      overlap_q = .true.
      exit
    endif
  ENDDO
  DEALLOCATE(jnear)
  RETURN
END SUBROUTINE overlap_pres
