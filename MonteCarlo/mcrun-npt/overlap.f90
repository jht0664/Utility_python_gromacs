! check overlap
! input: itry is index of molecule at original position
!        coord_new, to is position and typs of trial molecule
! output: overlap_q is output of overlap; if ovelaps, return true. Otherwise, false
SUBROUTINE overlap(ITRY,coord_new,to,action,overlap_q)
  use pos
  use cellmap
  use sigmas
  use calc_ex
  use calc_pres
  use coupling_pres
  use movetype
  IMPLICIT NONE
  integer :: itry, id, icell, jcell
  integer :: im, nn, k, nc
  DOUBLE PRECISION :: rt, dcomp, x1, y1, z1, xt, yt, zt
  double precision, dimension(3) :: coord_new
  CHARACTER(LEN=4) :: to
  CHARACTER(LEN=1) :: action
  logical :: overlap_q

!  write(*,*) "overlap:"
  if ( action /= 'x' .and. action /= 'c' .and. action /= 'p') THEN
    write(*,*) "wrong action argument"
    stop
  endif
! initialize (not necessary, it slows down to half)
  overlap_q = .false.
  X1 = coord_new(1) - box(1)*DNINT(coord_new(1)/box(1)-0.50D0)
  Y1 = coord_new(2) - box(2)*DNINT(coord_new(2)/box(2)-0.50D0) 
  Z1 = coord_new(3) - box(3)*DNINT(coord_new(3)/box(3)-0.50D0) 
  ICELL = 1 + INT( X1*CELLX ) &
            + INT( Y1*CELLY ) *MCX &
            + INT( Z1*CELLZ ) *MCY*MCX
  id = (icell-1)*27
! find particles in neighbor cells
  NITER = 0
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
  DO NN = 1, NITER ! # possible overlap particles calculated in previous loop
    IM = JNEAR(NN) 
    if(typs(IM) == TO) then
      cycle ! for same atom type
    else
      DCOMP = sigma_ab**2 ! for different atomtype
    endif
    XT = X(IM) - X1
    YT = Y(IM) - Y1
    ZT = Z(IM) - Z1
    XT = XT - box(1)*DNINT(XT/box(1))
    YT = YT - box(2)*DNINT(YT/box(2))
    ZT = ZT - box(3)*DNINT(ZT/box(3))
    RT = (XT*XT+YT*YT+ZT*ZT)*pos_scaling2
    !RT = XT*XT+YT*YT+ZT*ZT
    if (action == 'x' .and. rt < dcomp) then ! to avoid overlapping for try_conf (translation)
      overlap_q = .true.
      exit
    else if (action == 'c') then ! to count overlapping for pressure in NVT
      DO k = 1, n_dv  ! count number of overlapping pairs (finally, double counting)
        IF(RT*scale_length(k)**2 < DCOMP) THEN
          pres_noverlap(k) = pres_noverlap(k) + 1 
        ENDIF
      ENDDO
    else if (action == 'p' .and. (rt*(expd**2)) < dcomp) then ! to count overlapping for volume change in NPT
      overlap_q = .true.
      exit
    endif
  ENDDO
!  write(*,*) "overlap_q = ", overlap_q
  RETURN
END SUBROUTINE overlap
