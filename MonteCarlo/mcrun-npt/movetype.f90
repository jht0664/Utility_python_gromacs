! initilize movement types
SUBROUTINE movetype_init()
  use movetype
  IMPLICIT NONE
  integer :: ierr
  write(*,*) "movetype_init:"
  ALLOCATE(movetype_name(nmovetypes),stat=ierr)
  ALLOCATE(movetype_prob(nmovetypes),stat=ierr)
  ALLOCATE(movetype_ntry(nmovetypes),stat=ierr)
  ALLOCATE(movetype_nsuccess(nmovetypes),stat=ierr)
  ALLOCATE(movetype_i_success(nmovetypes),stat=ierr)
  movetype_ntry = 0
  movetype_nsuccess = 0
  movetype_i_success = 0
END SUBROUTINE

! initilize movement types
SUBROUTINE movetype_set_ensemble(name)
  use movetype
  use coupling_pres
  use coupling_exch
  IMPLICIT NONE
  CHARACTER(LEN=5) :: name
  integer :: ind_str

  write(*,*) "movetype_set_ensemble:"
  ensemble_temp = .False.
  ensemble_pres = .False.
  ensemble_exch = .False.
  name = adjustl(name)
  
  nmovetypes = 0
  IF (scan(name,"T") /= 0) THEN ! if ensemble_name contain "T"
    ensemble_temp = .True.
    write(*,*) " Turn on temperature coupling."
    nmovetypes = nmovetypes + 1
  ENDIF
  If (scan(name,"P") /= 0) THEN
    ensemble_pres = .True.
    nmovetypes = nmovetypes + 1
    if (scan(name,"Z") /= 0) THEN
      semiiso = 3
      write(*,*) " Turn on semiisotropic pressure coupling on z-direction."
      write(*,*) " But, not supported yet"
      stop
    else
      semiiso = 0
      write(*,*) " Turn on isotropic pressure coupling."
    endif
  ENDIF
  IF (scan(name,"X") /= 0) THEN
    ensemble_exch = .True.
    nmovetypes = nmovetypes + 1
    ind_str = scan(name,"X")+1
    exch_comp = name(ind_str:ind_str)
    write(*,*) " Turn on exchange identity coupling based on the component, ", exch_comp
  ENDIF
  if (nmovetypes == 0) THEN
    WRITE(*,*) " no coupling? because of your ", name, " ensemble."
    STOP
  ENDIF
  write(*,*) " total #couplings = ", nmovetypes
  return
END SUBROUTINE movetype_set_ensemble

SUBROUTINE movetype_read(fileunit)
  use movetype
  use trans
  use coupling_pres
  use coupling_exch

  IMPLICIT NONE
  CHARACTER(len=5) :: coupling
  integer :: fileunit, i, j, ierr
  double precision :: temp, tprob

  write(*,*) "movetype_read:"
  do i=1, nmovetypes
    read(fileunit,*) coupling
    write(*,*) coupling
    ! translation
    if (trim(coupling) == 'trans') THEN
      movetype_name(i) = 'trans'
      read(fileunit,*) movetype_prob(i), dlr_a, dlr_b ! probability, translation increment for components
      write(*,'(1x,A,1X,3(F10.5),1x,A)') " trans = ",movetype_prob(i), dlr_a, dlr_b, "(prob, dlr_a, dlr_b)"
      cycle
    ! pressure coupling
    else if (trim(coupling) == 'press') THEN
      movetype_name(i) = 'press'
      read(fileunit,*) movetype_prob(i), press_val, dvx, temp ! probability, pressure, volume change, temperature
      write(*,'(1x,A,1X,4(F10.5),1x,A)') " press = ", movetype_prob(i), press_val, dvx, temp, "(prob, press., volume change, temp.)"
      tinv = 1.0/temp
      cycle
    ! particle exchange coupling
    else if (trim(coupling) == 'exch') THEN
      movetype_name(i) = 'exch'
      read(fileunit,*) movetype_prob(i), xi_val, exch_ncomp ! probability, xi value, number of component
      write(*,'(1x,A,1X,2(F10.5),I5,1x,A)') " exch = ", movetype_prob(i), xi_val, exch_ncomp, "(prob, xi_val, exchange component)"
      if(exch_ncomp <= 1) then
        write(*,*) " exch_ncomp should be greater than one component"
        stop
      else if (xi_val >= 1.0) then
        write(*,*) " xi_val, fugacity ratio, should be less than one"
        stop
      end if
      log_xi2_of_xi1 = dlog(xi_val/(1.0D0-xi_val))
      ALLOCATE(exch_tcomp(exch_ncomp),stat=ierr)
      read(fileunit,*) (exch_tcomp(j), j=1, exch_ncomp)
      write(*,*) " component => ", (exch_tcomp(j), j=1, exch_ncomp)
      write(*,*) " From your setting, component change xi_(", exch_comp, &
                "), attempts from species ", exch_tcomp(1), " to species ", exch_tcomp(2)
      write(*,*) "which gives m = -1. Otherwise, m = 1"
      write(*,*) " (in other words, ", exch_comp, " is the same as ", exch_tcomp(1)
      cycle
    else
      write(*,*) " wrong some lines"
      stop
    endif
  enddo
  tprob = 0.0D0
  do i=1, nmovetypes
    tprob = tprob + movetype_prob(i)
  enddo
  if (tprob < 1.0D0 .or. tprob > 1.0D0) then
    write(*,*) " total probability is not 1, ", tprob, ". check your input file."
    stop
  endif
  return
END SUBROUTINE movetype_read

! exit movetype
SUBROUTINE movetype_exit
  use movetype
  use coupling_exch
  IMPLICIT NONE
  integer :: ierr
  write(*,*) "movetype_exit:"
  DEALLOCATE(movetype_name,stat=ierr)
  DEALLOCATE(movetype_prob,stat=ierr)
  DEALLOCATE(movetype_ntry,stat=ierr)
  DEALLOCATE(movetype_nsuccess,stat=ierr)
  DEALLOCATE(movetype_i_success,stat=ierr)
  if (ensemble_exch) DEALLOCATE(exch_tcomp,stat=ierr)
END SUBROUTINE movetype_exit

subroutine movetype_select(guess,select_type)
    use movetype, only: nmovetypes, movetype_prob
    implicit none
    logical :: check_movetype
    double precision :: guess, tprob
    integer :: itype, select_type
    check_movetype = .false.
    tprob = 0.0D0
    do itype = 1, nmovetypes
      tprob = tprob + movetype_prob(itype)
      if (guess < tprob) then 
        check_movetype = .true.
        exit
      endif
    enddo
    if (.not. check_movetype) THEN
      write(*,*) " wrong check_movetype"
      STOP
    endif
    select_type = itype
    return
end subroutine movetype_select

subroutine movetype_run(itype)
  use movetype, only: movetype_name
  use pos, only: box
  use ints, only: nptot
  use coupling_pres, only: dvx, delta, expd
  implicit none
  integer :: itype, imol
  double precision :: rand
  external rand
  if (movetype_name(itype) == 'trans') then
    imol = INT( NPTOT*RAND() ) + 1 ! pick random particle, (0 <= RAND < 1)
    call try_conf(imol,itype)
  else if (movetype_name(itype) == 'press') then
    delta = (2.0D0*RAND()-1.0D0)*dvx   ! pick volume change direction or total volume
    expd  = (box(3)+delta)/box(3)      ! get expansion ratio on z-axis (use this factor for other axes)
    call try_pres(itype)
  else if (movetype_name(itype) == 'exch') then
    imol = INT( NPTOT*RAND() ) + 1 ! pick random particle, (0 <= RAND < 1)
    call try_exch(imol,itype)
  else
    write(*,*) " wrong argument itype"
    stop
  endif
  return
end subroutine movetype_run
