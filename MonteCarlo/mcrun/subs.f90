! initilize movement types
SUBROUTINE movetype_init
  use movetype
  IMPLICIT NONE
  integer :: ierr
  write(*,*) "movetype_init:"
  ALLOCATE(movetype_prob(nmovetypes), stat=ierr)
  ALLOCATE(movetype_ntry(nmovetypes), stat=ierr)
  ALLOCATE(movetype_nsuccess(nmovetypes), stat=ierr)
  movetype_prob = 0.0D0
  movetype_ntry = 0
  movetype_nsuccess = 0
END SUBROUTINE

! exit movetype
SUBROUTINE movetype_exit
  use movetype
  IMPLICIT NONE
  integer :: ierr
  write(*,*) "movetype_exit:"
  DEALLOCATE(movetype_prob,stat=ierr)
  DEALLOCATE(movetype_ntry,stat=ierr)
  DEALLOCATE(movetype_nsuccess,stat=ierr)
END SUBROUTINE movetype_exit

! try movement
SUBROUTINE try_move(imol,itype)
  use pos, only: typs
  use movetype
  use try
  use ints
  IMPLICIT NONE
  integer :: imol, itype, array_size, ierr
  logical :: success

!  write(*,*) "try_move:"
  if( itype > nmovetypes) THEN
    write(*,*) "wrong itype argument."
    STOP
  endif

  if (typs(imol) == "A") THEN
    array_size = nmon_a
  else
    array_size = nmon_b
  endif
  ALLOCATE(xitr(array_size), stat=ierr)
  ALLOCATE(yitr(array_size), stat=ierr)
  ALLOCATE(zitr(array_size), stat=ierr)
  ALLOCATE(titr(array_size), stat=ierr)
    
  movetype_ntry(itype) = movetype_ntry(itype) + 1
  if (itype == 1) THEN
    call tran(imol,success)
  endif
!  write(*,*) "B"
  if(success) then
 !   write(*,*) "success"
    movetype_nsuccess(itype) = movetype_nsuccess(itype) + 1
    !call energy_update
    call cellmap_update(imol)
!    write(*,*) "success update"
  endif
  DEALLOCATE(xitr)
  DEALLOCATE(yitr)
  DEALLOCATE(zitr)
  DEALLOCATE(titr)

END SUBROUTINE try_move

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
  integer :: ierr
  
  success = .false.
  
  if (typs(imol) == "A") THEN
      dlr = dlr_a
    else
      dlr = dlr_b
    endif
    
  DX = dlr*(RAND()-0.50D0)
  DY = dlr*(RAND()-0.50D0)
  DZ = dlr*(RAND()-0.50D0)
!  write(*,*) "E"
  
  XITR(1) = X(imol) + DX
  YITR(1) = Y(imol) + DY
  ZITR(1) = Z(imol) + DZ
  TITR(1) = TYPS(imol)
!  write(*,*) "D"
  CALL overlap(imol,xitr(1),yitr(1),zitr(1),titr(1),'x',overlap_q) ! if overlaps, stop the subroutine

  if ( .not. overlap_q ) success = .true.

!  write(*,*) "C"
  return
END SUBROUTINE tran

! gibbs-duhem integration (exchange), semigrand ensemble
!SUBROUTINE GDI_EX(I,ISUC)
!  USE pos, only: x, y, z, typs, box
!  USE try, only: xitr, yitr, zitr, titr, titr
!  USE ints, only: nmol_a, nmon_a, nmol_b, nmon_b
!  USE trans, only: trans_prob, dlr_a, dlr_b
!  USE ens, only: ener,eold,enew
!  USE ene, only: enernew,enerold,eninew,eniold,eninew2,eniold2,enewtot,eoldtot
!  USE sigmas, only: sigma_b,sigma_a,sigma_ab
!  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!  include 'parameter.h'
!  
!  ISUC = 1
!  nptot = nmol_a + nmol_b
!  
!  ENEW = 0.0D0
!  EOLD = 0.0D0
!  XITR(1) = X(I)
!  YITR(1) = Y(I)
!  ZITR(1) = Z(I)
!  IF (TYPS(I) == 'A') THEN
!    TITR(1) = 'B'
!  ELSE
!    TITR(1) = 'A'
!  ENDIF
!  CALL overlap(I,1)
!  IF(NVL.EQ.1) RETURN
!  ISUC = 0
!  RETURN
!END SUBROUTINE GDI_EX

! save gro format trajectory
!SUBROUTINE save_gro_traj

