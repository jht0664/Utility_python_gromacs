SUBROUTINE cellmap_init()
  use pos, only: box
  use sigmas
  use cellmap
  use trans
  use ints, only: nptot
  IMPLICIT NONE
  DOUBLE PRECISION :: DCMAX
  integer :: ierr

  write(*,*) "cellmap_init:"
! initialize cell_list
  DCMAX = MAX(sigma_a,sigma_b,dlr_a,dlr_b)
  
  mcx = INT(box(1)/DCMAX)
  mcy = INT(box(2)/DCMAX)
  mcz = INT(box(3)/DCMAX)

  CELLX = DBLE(MCX)/box(1)
  CELLY = DBLE(MCY)/box(2)
  CELLZ = DBLE(MCZ)/box(3)
  CALL MAPS
  CALL LINKS
  Allocate(jnear(nptot),stat=ierr)
END SUBROUTINE cellmap_init

! initialize map array
! MAP has the indice of 27 neighbor cells of a given cell
! e.g. MAP(1) has the index of the first neighbor cell of 1st cell
! MAP(2) has the index of the 2nd neighbor cell of 1st cell
! MAP(3) has the index of the 3rd neighbor cell of 1st cell
! MAP(28) has the index of the first neighbor cell of 2nd cell
! MAP(55) has the index of the first neighbor cell of 3rd cell
SUBROUTINE MAPS
  use cellmap
  IMPLICIT NONE
  INTEGER :: ncell, mapsiz
  INTEGER :: ierr, icell, icx, icy, icz, id, kcell
  INTEGER :: ix, iy, iz, ixr, iyr, izr, ict, icn

  NCELL = MCX*MCY*MCZ ! total #cells of box
  write(*,*) " NCELL = ", ncell
  MAPSIZ = 27*NCELL ! mapping memory array size due to 27 neighbor cells of a given cell
  ! initialize map array
  Allocate(map(mapsiz),stat=ierr)
  MAP = 0

  DO ICELL = 1, NCELL ! index of cell
    ICX = MOD(ICELL,MCX)
    IF(ICX.EQ.0)ICX=MCX
    ICY = MOD( (ICELL-ICX)/MCX , MCY ) + 1 
    ICZ = (ICELL-ICX-MCX*(ICY-1))/MCY/MCX+1
    ID = (ICELL-1)*27 
    KCELL = 1
    ! do loop for neighbor cells of i-th cell
    DO IX = ICX-1,ICX+1
      IXR=IX
      IF(IXR.EQ.0)IXR=MCX
      IF(IXR.EQ.MCX+1)IXR=1
      DO IY = ICY-1,ICY+1
        IYR = IY
        IF(IYR.EQ.0)IYR=MCY
        IF(IYR.EQ.MCY+1)IYR=1
        DO IZ = ICZ-1,ICZ+1
          IZR = IZ
          IF(IZR.EQ.0)IZR=MCZ
          IF(IZR.EQ.MCZ+1)IZR=1
          ICT = MCY*MCX*(IZR-1)+MCX*(IYR-1)+IXR
          ICN = ID + KCELL
          MAP(ICN) = ICT
          KCELL = KCELL + 1
        ENDDO
      ENDDO
   ENDDO    
  ENDDO
  RETURN
END SUBROUTINE MAPS

! initialized the linked list and head element for checking overlap
! output: LIST has the next particle occupied in the same cell
! e.g. If the linked list is [43 | 9 | 2 | 1 | 0],
!       LIST(43) = 9 (= 24th particle also occupied in the same cell)
!       LIST(9) = 2, 
!       LIST(1) = 0 means 1nd particle is the last particle occupied in the same cell
!        In other words, it is the end element
!         LEAD has the index of head element of the linked list, LIST, in a given cell
! e.g. LEAD(4) returns index of head element of the 4th cell
!        IF 4th cell has 5th and 3th particle, LEAD(4) == 5, not 3.
SUBROUTINE LINKS
  USE pos
  USE ints
  use cellmap
  IMPLICIT NONE
  INTEGER :: ncell 
  INTEGER :: i, ierr, ic
  double precision :: x1, y1, z1
! allocate
  NCELL = MCX*MCY*MCZ
  Allocate(lead(ncell),stat=ierr)
  Allocate(list(nptot),stat=ierr)
  LEAD = 0
! MAKE A NEW LIST
  DO I = 1, NPTOT
    X1 = X(I)
    Y1 = Y(I)
    Z1 = Z(I)
    X1 = X1 - box(1)*DNINT(X1/box(1)-0.50D0)
    Y1 = Y1 - box(2)*DNINT(Y1/box(2)-0.50D0)
    Z1 = Z1 - box(3)*DNINT(Z1/box(3)-0.50D0)
    IC = 1 + INT( X1*CELLX ) &
           + INT( Y1*CELLY ) *MCX &
           + INT( Z1*CELLZ ) *MCY*MCX
    LIST(I) = LEAD(IC)
    LEAD(IC) = I
  ENDDO
  RETURN
END SUBROUTINE LINKS

SUBROUTINE cellmap_update(I)
    use pos
    use cellmap
    IMPLICIT NONE
    DOUBLE PRECISION :: xo, yo, zo, xn, yn, zn
    INTEGER :: ico, icn, jmol
    INTEGER :: i, j

    ! remove old X(i) FROM THE LIST
    XO = X(i)
    YO = Y(i)
    ZO = Z(i)
    XO = XO - box(1)*DNINT(XO/box(1)-0.50D0)
    YO = YO - box(2)*DNINT(YO/box(2)-0.50D0)
    ZO = ZO - box(3)*DNINT(ZO/box(3)-0.50D0)
    ICO = 1 + INT( XO*CELLX ) + INT( YO*CELLY ) *MCX + INT( ZO*CELLZ ) *MCY*MCX
    JMOL = LEAD(ICO)
    IF(JMOL.EQ.0) STOP
    IF(JMOL.NE.i)THEN
      DO WHILE (LIST(JMOL).NE.i)
        JMOL = LIST(JMOL)
      ENDDO
      LIST(JMOL) = LIST(i) 
    ELSE
      LEAD(ICO) = LIST(JMOL)
    ENDIF

    !   ADD XITR(J) TO THE LIST
    J = 1
    XN = coord_try(1) 
    YN = coord_try(2) 
    ZN = coord_try(3) 
    XN = XN - box(1)*DNINT(XN/box(1)-0.50D0)
    YN = YN - box(2)*DNINT(YN/box(2)-0.50D0)
    ZN = ZN - box(3)*DNINT(ZN/box(3)-0.50D0)
    ICN = 1 + INT( XN*CELLX ) + INT( YN*CELLY ) *MCX + INT( ZN*CELLZ ) *MCY*MCX
    JMOL = LEAD(ICN)
    IF(JMOL.NE.0)THEN
      DO WHILE (LIST(JMOL).NE.0)
        JMOL = LIST(JMOL)
      ENDDO
      LIST(JMOL) = i
    ELSE
      LEAD(ICN) = i
    ENDIF
    LIST(i) = 0
    X(i) = coord_try(1) 
    Y(i) = coord_try(2) 
    Z(i) = coord_try(3) 
END SUBROUTINE cellmap_update

SUBROUTINE cellmap_exit
  use cellmap
  INTEGER :: ierr
  DEALLOCATE(lead,stat=ierr)
  DEALLOCATE(map,stat=ierr)
  DEALLOCATE(list,stat=ierr)
  DEALLOCATE(jnear,stat=ierr)
END SUBROUTINE cellmap_exit

