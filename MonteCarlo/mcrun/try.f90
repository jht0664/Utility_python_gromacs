! try movement
SUBROUTINE try_conf(imol,itype)
  use movetype
  use ints
  use coupling_pres
  IMPLICIT NONE
  integer :: imol, itype
  logical :: success
!  write(*,*) "try_move:"
  movetype_i_try(itype) = movetype_i_try(itype) + 1
  call conf_tran(imol,success)
  if(success) then
    movetype_i_success(itype) = movetype_i_success(itype) + 1
    call cellmap_update(imol) ! xyz and cell list update
  endif
  RETURN
END SUBROUTINE try_conf


! translation algorithm
SUBROUTINE conf_tran(imol,success)
  use pos
  use trans
  use ints
  IMPLICIT NONE
  integer :: imol
  logical :: success, overlap_q
  double precision :: dlr, rand
  if (typs(imol) == "A") THEN
    dlr = dlr_a
  else
    dlr = dlr_b
  endif
  coord_try(1) = x(imol) + dlr*(rand()-0.50d0)
  coord_try(2) = y(imol) + dlr*(rand()-0.50d0)
  coord_try(3) = z(imol) + dlr*(rand()-0.50d0)
  typs_try = typs(imol)
  call overlap(imol,coord_try,typs_try,'x',overlap_q) ! if overlaps, stop the subroutine
  if ( .not. overlap_q ) then
    success = .true.
  else
    success = .false.
  endif
  return
end subroutine conf_tran

! try movement
SUBROUTINE try_pres(itype)
  use omp_lib
  use coupling_pres
  use pos
  use ints, only: nptot
  use inp
  use cellmap
  use movetype
  use sigmas
  IMPLICIT NONE
  double precision :: oldvol, rand, xo, yo, zo, expd2
  double precision :: delh, boltz, xt, yt, zt, rt
  integer :: itype, i, j
  logical :: success
  integer :: id, start_k, end_k, kmax, k ! for openmp
  logical :: overlap_q        ! for openmp
  external rand
!  write(*,*) "try_pres:"
  movetype_i_try(itype) = movetype_i_try(itype) + 1
  expd2 = expd**2
! check overlaps
  if (expd < 1.0) then
!!    single node job
!    do i=1,nptot-1 
!    do j=i+1, nptot
!    if(typs(i) == typs(j)) cycle
!    xt = x(i) - x(j)
!    yt = y(i) - y(j)
!    zt = z(i) - z(j)
!    XT = XT - box(1)*DNINT(XT/box(1))
!    YT = YT - box(2)*DNINT(YT/box(2))
!    ZT = ZT - box(3)*DNINT(ZT/box(3))
!    rt = (XT*XT+YT*YT+ZT*ZT)*expd2*pos_scaling2
!    if( rt < sigma_ab**2 ) return
!    enddo
!    enddo
!  endif
!   openmp version
    overlap_q = .false.
    start_k = 0
    kmax = 16 ! might be better performance if the optimized value used.
    do k=1,kmax
    end_k = k*(nptot-1)/kmax
!$omp parallel private (xt, yt, zt, rt, id, i, j) shared (overlap_q)
    id = omp_get_thread_num() + 1
    do i=start_k+id,end_k,openmp_thread_num
      do j=i+1, nptot
        if(typs(i) == typs(j)) cycle
        xt = x(i) - x(j)
        yt = y(i) - y(j)
        zt = z(i) - z(j)
        XT = XT - box(1)*DNINT(XT/box(1))
        YT = YT - box(2)*DNINT(YT/box(2))
        ZT = ZT - box(3)*DNINT(ZT/box(3))
        RT = (XT*XT+YT*YT+ZT*ZT)*expd2*pos_scaling2
        if( (rt < sigma_ab**2) .or. overlap_q) then
          overlap_q = .true.
          exit
        endif
      enddo
      if(overlap_q) exit
    enddo
!$omp end parallel
    if(overlap_q) return
    start_k = end_k
    enddo
  endif
! accepted or denied?
  success = .false.
  ! if volume change is on 3-axis
  !oldvol = box(1)*box(2)*box(3)*(pos_scaling**3)
  oldvol = box(1)*box(2)*box(3)*pos_scaling3
  delh = tinv*press_val*oldvol*(expd**3-1.0D0) - dble(nptot+1)*DLOG(expd)*3.0D0
  if(delh <= 0.0D0) then
    !write(*,*) "negative", delh
    success = .true.
  else 
    boltz = dexp(-delh)

    if(rand() < boltz) then
      success = .true.
    endif
  endif  
! accepted
  if(success) then
    movetype_i_success(itype) = movetype_i_success(itype) + 1
    !call energy_update
! position and box update
    pos_scaling  = pos_scaling*expd
    pos_scaling2 = pos_scaling**2
    pos_scaling3 = pos_scaling**3
  endif
  return
END SUBROUTINE try_pres

subroutine try_exch(imol,itype)
  use pos
  use movetype
  use ints
  use coupling_exch
  IMPLICIT NONE
  integer :: imol, itype, try_type
  double precision :: m_val, delxi, rand, boltz
  logical :: success, overlap_q
  external rand
  
!  write(*,*) "try_exch:"
  movetype_i_try(itype) = movetype_i_try(itype) + 1
  !write(*,*) "movetype_i_try exch", movetype_i_try(itype)
  DO while (.true.)
    try_type = INT(exch_ncomp*RAND() ) + 1
    if(typs(imol) == exch_tcomp(try_type)) then
      cycle
    else
      exit ! continue
    endif
  enddo
  write(*,*) "exchange =>", typs(imol), exch_tcomp(try_type)
  coord_try(1) = x(imol)
  coord_try(2) = y(imol)
  coord_try(3) = z(imol)
  CALL overlap(imol,coord_try,exch_tcomp(try_type),'x',overlap_q) ! if overlaps, stop the subroutine
  if(overlap_q) return
! determine m value
  if ( (typs(imol) == exch_tcomp(1)) .and. (exch_tcomp(try_type) == exch_tcomp(2)) ) then
    m_val = -1.0D0

  else
    m_val = 1.0D0
  endif
  !write(*,*) "m val = ", m_val, "due to", typs(imol), "=",exch_tcomp(1),"and",exch_tcomp(try_type),"=",exch_tcomp(2)
! accepted or denied?
  success = .false.
  delxi = m_val*log_xi
  if(delxi >= 0.0D0) then
    success = .true.
  else 
    boltz = dexp(delxi)
    !write(*,*) "boltz = ", boltz
    if(rand() < boltz) then
      success = .true.
    endif
  endif
  ! update new particle
  if(success) then
    !write(*,*) "exchange result: ",imol," th particle ",typs(imol), " => ",exch_tcomp(try_type)
    movetype_i_success(itype) = movetype_i_success(itype) + 1
    !write(*,*) movetype_i_success(itype)
    if (typs(imol) == 'A') then
      nmol_a = nmol_a - 1
    else if (typs(imol) == 'B') then
      nmol_b = nmol_b - 1
    else 
      write(*,*) "no info component"
      stop
    endif
    typs(imol) = exch_tcomp(try_type)
    if (typs(imol) == 'A') then
      nmol_a = nmol_a + 1
    else if (typs(imol) == 'B') then
      nmol_b = nmol_b + 1
    else 
      write(*,*) "no info component"
      stop
    endif
  endif
  RETURN
end subroutine try_exch