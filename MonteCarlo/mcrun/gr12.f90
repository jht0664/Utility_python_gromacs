SUBROUTINE gr12_init
  use rdf, only: nbin, gr12
  IMPLICIT NONE
  integer :: ierr
  write(*,*) "gr12_init: "
  allocate(gr12(nbin),stat=ierr)
  gr12 = 0.0D0

END SUBROUTINE gr12_init

SUBROUTINE gr12_cal
  USE pos
  USE ints
  USE rdf
  use omp_lib
  IMPLICIT NONE
  INTEGER, DIMENSION(nbin) :: ndum
  integer :: i, k, kbin
  double precision :: dxo, dyo, dzo, dx, dy, dz, dist
  double precision :: rhot, const, x1, x2, rlower, rupper, rideal, vol12
  DOUBLE PRECISION, PARAMETER :: pi=3.141592653589793D0

!  write(*,*) "gr12_cal: "
  ndum = 0

! count all atoms in a shell
!$omp parallel private ( i, k, kbin, dxo, dyo, dzo, dx, dy, dz, dist ) shared ( ndum )  
!$omp do
  DO I = 1, NPTOT-1
!    write(*,*) "gr12", OMP_get_thread_num(), " I ", I
    DO K = I+1, NPTOT
      IF (TYPS(I) .NE. TYPS(K)) THEN
        DXO = X(K)-X(I)
        DYO = Y(K)-Y(I)
        DZO = Z(K)-Z(I)
        DX = DXO - box(1)*DNINT(DXO/box(1))
        DY = DYO - box(2)*DNINT(DYO/box(2))
        DZ = DZO - box(3)*DNINT(DZO/box(3))
        DIST = DSQRT(DX*DX+DY*DY+DZ*DZ)
        KBIN = 1 + INT(DIST/BSZ)
        IF(KBIN.LE.NBIN) THEN
          ndum(KBIN) = ndum(KBIN) + 2
        ENDIF
      ENDIF
    ENDDO
  ENDDO
!$omp end do
!$omp end parallel

  RHOT = DBLE(NPTOT)/(box(1)*box(2)*box(3))
  CONST = 4.0*PI*RHOT/3.0
  X1 = DBLE(nmol_a)/DBLE(NPTOT)
  X2 = 1.0D0 - X1
  DO I= 1, NBIN
    RLOWER = DBLE(I-1)*BSZ
    RUPPER = RLOWER + BSZ
    RIDEAL = CONST*(RUPPER**3 - RLOWER**3)
    VOL12 = DBLE(NPTOT) * RIDEAL *2.0D0* X1*X2
    GR12(I) = GR12(I) + NDUM(I)/VOL12
  ENDDO
  RETURN
END SUBROUTINE gr12_cal

SUBROUTINE gr12_save(filename, naver)
  use rdf
  IMPLICIT NONE
  CHARACTER(LEN=*) :: filename
  integer :: i, naver
  double precision :: distr

  write(*,*) "gr12_save: "
  ! save g(r)_12.out
  OPEN(UNIT=12,FILE=filename,STATUS='UNKNOWN')
  DO I= 1, NBIN                              
    GR12(I) = GR12(I)/DBLE(NAVER)
    DISTR = (DBLE(I)*2.0 - 1.0)*BSZ/2.0
    !GDIFF = GR11(I)+GR22(I)-2.0D0*GR12(I)
    !WRITE (12,119) DISTR, GR11(I), GR12(I), GR22(I), GDIFF 
    WRITE (12,120) DISTR, GR12(I)
  ENDDO
  CLOSE(12)
119 FORMAT(1X,F6.3,2X,6(F15.7,1X))
120 FORMAT(1X,F6.3,2X,F15.7)  
END SUBROUTINE gr12_save