SUBROUTINE try_conf_init()
  use try
  use ints, only: nmon_a, nmon_b
  implicit none
  integer :: array_size, ierr
  array_size = max(nmon_a,nmon_b)
  ALLOCATE(xitr(array_size), stat=ierr)
  ALLOCATE(yitr(array_size), stat=ierr)
  ALLOCATE(zitr(array_size), stat=ierr)
  ALLOCATE(titr(array_size), stat=ierr)
  return
end subroutine try_conf_init

! try movement
SUBROUTINE try_conf(imol,itype)
  use movetype
  use try
  use ints
  IMPLICIT NONE
  integer :: imol, itype
  logical :: success

!  write(*,*) "try_move:"
  if( itype > nmovetypes) THEN
    write(*,*) "wrong itype argument."
    STOP
  endif

  movetype_ntry(itype) = movetype_ntry(itype) + 1
  call tran(imol,success)
!  write(*,*) "B"
  if(success) then
 !   write(*,*) "success"
    movetype_nsuccess(itype) = movetype_nsuccess(itype) + 1
    !call energy_update
    call cellmap_update(imol)
!    write(*,*) "success update"
  endif
  RETURN
END SUBROUTINE try_conf

! translation algorithm
SUBROUTINE tran(imol,success)
  use pos
  use trans
  use try
  use ints
  IMPLICIT NONE
  integer :: imol
  logical :: success, overlap_q
  double precision :: dlr, dx, dy, dz, rand
  success = .false.
  if (typs(imol) == "A") THEN
    dlr = dlr_a
  else
    dlr = dlr_b
  endif
  DX = dlr*(RAND()-0.50D0)
  DY = dlr*(RAND()-0.50D0)
  DZ = dlr*(RAND()-0.50D0)
  XITR(1) = X(imol) + DX
  YITR(1) = Y(imol) + DY
  ZITR(1) = Z(imol) + DZ
  TITR(1) = TYPS(imol)
  CALL overlap(imol,xitr(1),yitr(1),zitr(1),titr(1),'x',overlap_q) ! if overlaps, stop the subroutine
  if ( .not. overlap_q ) success = .true.
  return
END SUBROUTINE tran

subroutine try_conf_exit()
  use try
  implicit none
  DEALLOCATE(xitr)
  DEALLOCATE(yitr)
  DEALLOCATE(zitr)
  DEALLOCATE(titr)
  return   
end subroutine try_conf_exit
  

! try movement
SUBROUTINE try_pres(itype)
  use coupling_pres
  use pos
  use ints, only: nptot
  use try
  use cellmap_try
  use movetype
  use omp_lib
  use inp, only: openmp_thread_num
  IMPLICIT NONE
  double precision :: rand
  double precision :: delh, boltz
  integer :: itype, i
  logical :: overlap_q, success

!  write(*,*) "try_pres:"
! rescale box size
!  do i=1,3
!    boxitr(i) = box(i)*expd
!  enddo
!  if (semiiso /= 0) then
!    boxitr(semiiso) = box(semiiso)*expd
!  endif
! try volume change
  movetype_ntry(itype) = movetype_ntry(itype) + 1
  ! for semi pressure system
  !call try_pres_init(semiiso,expd) ! rescale positions
  !call cellmap_itr_init ! initialize new cellmap_itr
! check overlap
 ! call omp_set_num_threads(openmp_thread_num)
  if (expd < 1.0) then
!!$omp parallel do reduction(.OR.:overlap_q)   
    do i=1, nptot
      !call overlap_pres(i,overlap_q)
!      write(*,*) "pres", OMP_get_thread_num(), " I ", i
      CALL overlap(i,x(i),y(i),z(i),typs(i),'p',overlap_q)
      if(overlap_q) exit
    enddo
!!$omp end parallel do
!    if(.not. overlap_q) write(*,*) "overlap_q = ", overlap_q
  else
    overlap_q = .false.
  endif
! accepted or denied?
  success = .false.
  if(.not. overlap_q) then
    if (semiiso /= 0) then ! if volume change is on 1-axis
      delh = tinv*press_val*box(1)*box(2)*box(3)*(expd-1.0D0) - dble(nptot+1)*delta
    else ! if volume change is on 3-axis
      delh = tinv*press_val*box(1)*box(2)*box(3)*(expd**3-1.0D0) - dble(nptot+1)*delta
    endif
    if(delh <= 0.0D0) then
      success = .true.
    else 
      boltz = dexp(-delh)
      if(rand() < boltz) then
        success = .true.
      endif
    endif
  endif
  if(success) then
 !   write(*,*) "success"
    movetype_nsuccess(itype) = movetype_nsuccess(itype) + 1
    !call energy_update
    call pos_itr_update
    !call cellmap_itr_update
    call cellmap_init
!    write(*,*) "success update"
  endif
  !call try_pres_exit
  !call cellmap_itr_exit
  return
END SUBROUTINE try_pres

subroutine try_pres_init(ibox,expd)
  use try
  use pos
  use ints, only: nptot
  implicit none
  integer :: ierr, ibox
  double precision :: expd, xo, yo, zo, i
  allocate(xitr(nptot),yitr(nptot),zitr(nptot),titr(nptot),stat=ierr)
  ! rescale positions
  do i=1, nptot
    if (ibox == 1) then
      xo = x(i)
      xo = xo - box(ibox)*DNINT(xo/box(ibox)-0.50D0)
      xitr(i) = xo*expd
      yitr(i) = y(i)
      zitr(i) = z(i)
      titr(i) = typs(i)
    else if (ibox == 2) then
      xitr(i) = x(i)
      yo = y(i)
      yo = yo - box(ibox)*DNINT(yo/box(ibox)-0.50D0)
      yitr(i) = yo*expd
      zitr(i) = z(i)
      titr(i) = typs(i)
    else if (ibox == 3) then
      xitr(i) = x(i)
      yitr(i) = y(i)
      zo = z(i)
      zo = zo - box(ibox)*DNINT(zo/box(ibox)-0.50D0)
      zitr(i) = zo*expd
      titr(i) = typs(i)
    else if (ibox == 0) then
      xo = x(i)
      xo = xo - box(1)*DNINT(xo/box(1)-0.50D0)
      xitr(i) = xo*expd
      yo = y(i)
      yo = yo - box(2)*DNINT(yo/box(2)-0.50D0)
      yitr(i) = yo*expd
      zo = z(i)
      zo = zo - box(3)*DNINT(zo/box(3)-0.50D0)
      zitr(i) = zo*expd
      titr(i) = typs(i)
    endif
  enddo
  return
end subroutine try_pres_init

subroutine try_pres_exit()
  use try
  implicit none
  DEALLOCATE(xitr)
  DEALLOCATE(yitr)
  DEALLOCATE(zitr)
  DEALLOCATE(titr)
end subroutine try_pres_exit

subroutine try_exch(imol,itype)
  use pos
  use movetype
  use ints
  use coupling_exch
  IMPLICIT NONE
  integer :: imol, itype, try_type
  double precision :: m_val, delxi, rand, boltz
  logical :: success, overlap_q

!  write(*,*) "try_exch:"
  if( itype > nmovetypes) THEN
    write(*,*) "wrong itype argument."
    STOP
  endif
  movetype_ntry(itype) = movetype_ntry(itype) + 1
  DO while (.true.)
    try_type = INT(exch_ncomp*RAND() ) + 1
    if(typs(imol) == exch_tcomp(try_type)) then
      cycle
    else
      CALL overlap(imol,x(imol),y(imol),z(imol),exch_tcomp(try_type),'x',overlap_q) ! if overlaps, stop the subroutine
      exit
    endif
  enddo
! determine m value
  if (typs(imol) == exch_tcomp(1) .and. exch_tcomp(try_type) == exch_tcomp(2)) then
    m_val = -1.0
  else
    m_val = 1.0
  endif
! accepted or denied?
  success = .false.
  if(.not. overlap_q) then
    delxi = m_val*log_xi2_of_xi1
    if(delxi <= 0.0D0) then
      success = .true.
    else 
      boltz = dexp(-delxi)
      if(rand() < boltz) then
        success = .true.
      endif
    endif
  endif
  ! update new particle
  if(success) then
    movetype_nsuccess(itype) = movetype_nsuccess(itype) + 1
    if (typs(imol) == 'A') then
      nmol_a = nmol_a - 1
    else
      nmol_b = nmol_b - 1
    endif
    typs(imol) = exch_tcomp(try_type)
    if (typs(imol) == 'A') then
      nmol_a = nmol_a + 1
    else
      nmol_b = nmol_b + 1
    endif
    if( nptot /= nmol_a*nmon_a+nmol_b*nmon_b) then
      write(*,*) "wrong total particle number during exchange"
      stop
    endif
  endif
  RETURN
end subroutine try_exch