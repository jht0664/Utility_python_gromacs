! pressure calculation using volume change in NVT ensemble
SUBROUTINE print_pressure_dV(outputfile)
  use omp_lib
  use calc_pres
  use ints, only: nptot
  use inp, only : openmp_thread_num
  use pos
  use sigmas
  IMPLICIT NONE
  CHARACTER(LEN=*) :: outputfile
  LOGICAL :: file_exist
  double precision :: xt, yt, zt, rt, vol, dcomp
  integer :: id, i, j, k
  integer, dimension(openmp_thread_num,n_dv) :: ompv

!  write(*,*) "pressure_dV:"
  pres_noverlap = 0
  ompv = 0
! count number of overlap paris under virtual volume
  DCOMP = sigma_ab**2
!$omp parallel private (xt, yt, zt, rt, id, i, j, k) shared (ompv)
  id = omp_get_thread_num() + 1
  do i=id, nptot-1, openmp_thread_num
    do j=i+1, nptot
      if(typs(i) == typs(j)) cycle
      xt = x(i) - x(j)
      yt = y(i) - y(j)
      zt = z(i) - z(j)
      XT = XT - box(1)*DNINT(XT/box(1))
      YT = YT - box(2)*DNINT(YT/box(2))
      ZT = ZT - box(3)*DNINT(ZT/box(3))
      RT = (XT*XT+YT*YT+ZT*ZT) !*pos_scaling2
      DO k=1, n_dv  ! count number of overlapping pairs (finally, double counting)
        IF(RT*scale_length(k) < DCOMP) THEN
          ompv(id,k) = ompv(id,k) + 1 
          !write(*,*) id,k,ompv(id,k)
        ENDIF
      ENDDO
    enddo
  enddo
!$omp end parallel
  do k=1, n_dv
    do j=1, openmp_thread_num
      pres_noverlap(k) = pres_noverlap(k) + ompv(j,k)  
    enddo
  enddo
  ! calculate pressure of WR model
  !  beta*Pressure = limit(del_V -> 0) <number density - (number of ovelapping A and B pair) / del_V >_ensemble average
  !   to remove duplicates by dividing by 2
  vol = box(1)*box(2)*box(3)*(pos_scaling**3)
  !vol = box(1)*box(2)*box(3)
  DO i = 1, n_dv
    !WRITE(*,*) DBLE(noverlap(i))/(ratio_dv_v*i)/2.0D0
    print_press(i) = DBLE(nptot)/vol - DBLE(pres_noverlap(i))/(ratio_dv_v*dble(i)*vol)
  ENDDO
  ! save data
  INQUIRE(file=outputfile, exist=file_exist)
  IF(file_exist) THEN
    OPEN(UNIT=1,FILE=outputfile,status='old', position='append', action='write')
  ELSE
    OPEN(UNIT=1,FILE=outputfile,status='new', action='write')  
    WRITE(1,'(A,2x,*(E16.8,1x))') "# delta_V's", (ratio_dv_v*i, i=1, n_dv)
  ENDIF
  WRITE(1,'(*(E16.8,1x))') ( print_press(i), i=1,n_dv )
  WRITE(*,*) "========= Instance Pressures with virtual volume change ======="
  WRITE(*,'(*(E16.8,2x))') ( print_press(i), i=1,n_dv )
  CLOSE(1)
  RETURN
END SUBROUTINE print_pressure_dV

! surface tension calculation using volume change in NVT ensemble
! be aware that surface tension distribution code does NOT consider using openmp
SUBROUTINE print_surface_tension_dV(outputfile1, outputfile2)
  use omp_lib
  use calc_pres
  use ints, only: nptot
  use inp, only : openmp_thread_num
  use pos
  use sigmas
  IMPLICIT NONE
  CHARACTER(LEN=*) :: outputfile1, outputfile2
  LOGICAL :: file_exist
  double precision :: xt, yt, zt, rt, vol, dcomp
  integer :: id, i, j, k
  integer, dimension(openmp_thread_num,n_dv) :: ompvxy, ompvz
  integer, parameter :: distr_nbin = 200
  integer :: distr_ibin 
  double precision :: distr_binsize, zcen, zpos
  integer, dimension(n_dv,distr_nbin) :: vxy_z, vz_z
  double precision, dimension(n_dv,distr_nbin) :: distr_xy, distr_z
  CHARACTER(LEN=100) :: file3, file4

!  write(*,*) "pressure_dV:"
  presxy_noverlap = 0
  presz_noverlap = 0
  ompvxy = 0
  ompvz = 0
  distr_binsize = box(3)/DBLE(distr_nbin)
  vxy_z = 0
  vz_z = 0
! count number of overlap pairs under virtual volume in area
  DCOMP = sigma_ab**2
!$omp parallel private (xt, yt, zt, rt, id, i, j, k) shared (ompvxy, ompvz)
  id = omp_get_thread_num() + 1
  do i=id, nptot-1, openmp_thread_num
    do j=i+1, nptot
      if(typs(i) == typs(j)) cycle
      xt = x(i) - x(j)
      yt = y(i) - y(j)
      zt = z(i) - z(j)
      XT = XT - box(1)*DNINT(XT/box(1))
      YT = YT - box(2)*DNINT(YT/box(2))
      ZT = ZT - box(3)*DNINT(ZT/box(3))
      DO k=1, n_dv  ! count number of overlapping pairs (finally, double counting)
        RT = (XT*XT+YT*YT)*scale_length_xy(k)+(ZT*ZT) !*pos_scaling2
        IF(RT < DCOMP) THEN
          ompvxy(id,k) = ompvxy(id,k) + 1 
          !write(*,*) "xy",rt,k,ompvxy(id,k),i,j
          zpos = z(i) - box(3)*DNINT(z(i)/box(3)-0.50D0) 
          zcen = zpos + dabs(ZT)
          zcen = zcen - box(3)*DNINT(zcen/box(3)-0.50D0)
          distr_ibin = int(zcen/distr_binsize) + 1
          vxy_z(k,distr_ibin) = vxy_z(k,distr_ibin) + 1
        ENDIF
        RT = XT*XT+YT*YT+ZT*ZT*scale_length_z(k) !*pos_scaling2
        IF(RT < DCOMP) THEN
          ompvz(id,k) = ompvz(id,k) + 1 
          zpos = z(i) - box(3)*DNINT(z(i)/box(3)-0.50D0) 
          zcen = zpos + dabs(ZT)
          zcen = zcen - box(3)*DNINT(zcen/box(3)-0.50D0)
          distr_ibin = int(zcen/distr_binsize) + 1
          vz_z(k,distr_ibin) = vz_z(k,distr_ibin) + 1
          !write(*,*) "z ",rt,k,ompvz(id,k),i,j
        ENDIF
      ENDDO
    enddo
  enddo
!$omp end parallel
  do k=1, n_dv
    do j=1, openmp_thread_num
      presxy_noverlap(k) = presxy_noverlap(k) + ompvxy(j,k)  
      presz_noverlap(k) = presz_noverlap(k) + ompvz(j,k)  
    enddo
  enddo
  ! calculate pressure of WR model
  !  beta*Pressure = limit(del_V -> 0) <number density - (number of ovelapping A and B pair) / del_V >_ensemble average
  !   to remove duplicates by dividing by 2
  vol = box(1)*box(2)*box(3)*(pos_scaling**3)
  !vol = box(1)*box(2)*box(3)
  DO i = 1, n_dv
    !WRITE(*,*) DBLE(noverlap(i))/(ratio_dv_v*i)/2.0D0
    print_presxy(i) = DBLE(nptot)/vol - DBLE(presxy_noverlap(i))/(ratio_dv_v*dble(i)*vol)
    print_presz(i) = DBLE(nptot)/vol - DBLE(presz_noverlap(i))/(ratio_dz_z*dble(i)*vol)
    DO j = 1, distr_nbin
      distr_xy(i,j) = DBLE(nptot)/vol - DBLE(vxy_z(i,j))/(ratio_dv_v*dble(i)*vol)
      distr_z(i,j) = DBLE(nptot)/vol - DBLE(vz_z(i,j))/(ratio_dv_v*dble(i)*vol)
    ENDDO
  ENDDO
  ! save data1
  INQUIRE(file=outputfile1, exist=file_exist)
  IF(file_exist) THEN
    OPEN(UNIT=1,FILE=outputfile1,status='old', position='append', action='write')
  ELSE
    OPEN(UNIT=1,FILE=outputfile1,status='new', action='write')  
    WRITE(1,'(A,2x,*(E16.8,1x))') "# delta_V's", (ratio_dv_v*i, i=1, n_dv)
  ENDIF
  WRITE(1,'(*(E16.8,1x))') ( print_presxy(i), i=1,n_dv )
  CLOSE(1)
  ! save data2
  INQUIRE(file=outputfile2, exist=file_exist)
  IF(file_exist) THEN
    OPEN(UNIT=2,FILE=outputfile2,status='old', position='append', action='write')
  ELSE
    OPEN(UNIT=2,FILE=outputfile2,status='new', action='write')  
    WRITE(2,'(A,2x,*(E16.8,1x))') "# delta_V's", (ratio_dz_z*i, i=1, n_dv)
  ENDIF
  WRITE(2,'(*(E16.8,1x))') ( print_presz(i), i=1,n_dv )
  CLOSE(2)

  ! save data1
  file3 = 'surf_dist_xy_z.out'
  INQUIRE(file=file3, exist=file_exist)
  IF(file_exist) THEN
    OPEN(UNIT=3,FILE=file3,status='old', position='append', action='write')
  ELSE
    OPEN(UNIT=3,FILE=file3,status='new', action='write')  
    !WRITE(3,'(A,2x,*(E16.8,1x))') "# delta_V's", (ratio_dv_v*i, i=1, n_dv)
  ENDIF
  do i = 1, n_dv
    WRITE(3,'(*(E16.8,1x))') ( distr_xy(i,j), j=1, distr_nbin)
  enddo
  CLOSE(3)
  ! save data2
  file4 = 'surf_dist_z_z.out'
  INQUIRE(file=file4,exist=file_exist)
  IF(file_exist) THEN
    OPEN(UNIT=4,FILE=file4,status='old', position='append', action='write')
  ELSE
    OPEN(UNIT=4,FILE=file4,status='new', action='write')  
    !WRITE(4,'(A,2x,*(E16.8,1x))') "# delta_V's", (ratio_dz_z*i, i=1, n_dv)
  ENDIF
  do i = 1, n_dv
    WRITE(4,'(*(E16.8,1x))') ( distr_z(i,j), j=1, distr_nbin)
  enddo
  CLOSE(4)
  
  !WRITE(*,*) "========= Instance Pressures with virtual volume change (xy, z) ======="
  !WRITE(*,'(*(E16.8,2x))') ( print_presxy(i), i=1,n_dv )
  !WRITE(*,'(*(E16.8,2x))') ( print_presz(i), i=1,n_dv )
  RETURN
END SUBROUTINE print_surface_tension_dV

SUBROUTINE check_overlap()
  use omp_lib
  use coupling_pres
  use pos
  use ints, only: nptot
  use inp
  use cellmap
  use movetype
  use sigmas
  IMPLICIT NONE
  double precision ::  xt, yt, zt, rt
  integer :: i, j
  integer :: id, start_k, end_k, kmax, k ! for openmp
  logical :: overlap_q        ! for openmp

  write(*,*) "check_overlap:"
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
      RT = (XT*XT+YT*YT+ZT*ZT)*pos_scaling2
      if( (rt < sigma_ab**2) .or. overlap_q) then
        overlap_q = .true.
        exit
      endif
    enddo
    if(overlap_q) exit
  enddo
!$omp end parallel
  if(overlap_q) exit
  start_k = end_k
  enddo
  if(overlap_q) then
    write(*,*) " your coordinate has overlapped particles. Maybe wrong initial coord."
    write(*,*) " So, it always check until no overlapping"
    check_overlap_log = .true.
  else
    write(*,*) " Passed."
    check_overlap_log = .false.
  endif
  return
END SUBROUTINE check_overlap

SUBROUTINE print_density(outputfile)
  use ints
  use pos
 
  IMPLICIT NONE
  CHARACTER(LEN=*) :: outputfile
  double precision :: vol
  double precision, dimension(3) :: dens
  LOGICAL :: file_exist
  integer :: i
  
  ! calculate density
  !vol = box(1)*box(2)*box(3)
  vol = box(1)*box(2)*box(3)*pos_scaling3
  dens(1) = dble(nptot)/vol
  dens(2) = dble(nmol_a*nmon_a)/vol
  dens(3) = dble(nmol_b*nmon_b)/vol

  ! save data
  INQUIRE(file=outputfile, exist=file_exist)
  IF(file_exist) THEN
    OPEN(UNIT=1,FILE=outputfile,status='old', position='append', action='write')
  ELSE
    OPEN(UNIT=1,FILE=outputfile,status='new', action='write')  
    WRITE(1,*) "# density files"
  ENDIF
!  WRITE(*,*) "========= Density ======="
!  WRITE(*,'(*(E16.8,2x))') ( dens(i), i=1,3 )
  WRITE(1,'(*(E16.8,2x))') ( dens(i), i=1,3 )
  CLOSE(1)
  RETURN
END SUBROUTINE print_density

subroutine update_dvx()
  use inp, only: nskip
  use movetype
  use coupling_pres
  implicit none
  integer :: i
  double precision :: frac
  do i=1, nmovetypes
    if( movetype_name(i) == 'press') then
      frac = dble(movetype_i_success(i))/dble(nskip)/movetype_prob(i)
      if( frac > 0.55 ) then
        dvx = dvx*1.05
        write(*,*) "increase dvx by factor 0.05 =>", real(dvx), "due to frac ", real(frac)
      else if( frac < 0.45) then
        dvx = dvx*0.95
        write(*,*) "decrease dvx by factor 0.05 => ", real(dvx), "due to frac ", real(frac)
      endif
      exit
    endif
  enddo
  return
end subroutine update_dvx

! pressure calculation using volume change in NVT ensemble
SUBROUTINE print_hifj(outputfile)
  use calc_ex
  use ints, only: nptot
  use pos
 
  IMPLICIT NONE
  CHARACTER(LEN=*) :: outputfile
  double precision :: rand
  INTEGER :: imol, itry
  LOGICAL :: file_exist, overlap_q
  double precision, dimension(3), save :: coord
  ex_noverlap = 0
  itry = 0
  do while (itry < ex_ntry)
    imol = INT( NPTOT*RAND() ) + 1 ! pick a molecule
    if(typs(imol) == ex_solvent) THEN
      itry = itry + 1
      coord(1) = x(imol)
      coord(2) = y(imol)
      coord(3) = z(imol)
      call overlap(imol,coord,ex_solute,'x',overlap_q) !check overlap
      if(overlap_q) THEN
        ex_noverlap = ex_noverlap + 1
      endif
    else
      cycle
    endif
  ENDDO
  ! calculate fj/Hi and save data
  !  using exchange method by converting solvent to solute
  !  See eq (3.1) in "Coexistence diagrams of mixtures by molecular simulation" M. Mehta and D. A. Kofke.
  !   or see eq (10) in "Direct Evaluation of Vapour-Liquid Equilibria of Mixtures by Molecular Dynamics Using Gibbs-Duhem Integration" 
  !                     Martin Lísal & Václav Vacek
  INQUIRE(file=outputfile, exist=file_exist)
  IF(file_exist) THEN
    OPEN(UNIT=1,FILE=outputfile,status='old', position='append', action='write')
  ELSE
    OPEN(UNIT=1,FILE=outputfile,status='new', action='write')  
    WRITE(1,*) "# ntry = ", ex_ntry
  ENDIF
  WRITE(*,*) " >> fj/Hi = ", 1.0D0-DBLE(ex_noverlap)/DBLE(ex_ntry)
  WRITE(1,'(E16.8)') 1.0D0-DBLE(ex_noverlap)/DBLE(ex_ntry)
  CLOSE(1)
  RETURN
END SUBROUTINE print_hifj

SUBROUTINE print_nmol(outputfile)
  use ints
  use coupling_exch
  IMPLICIT NONE
  CHARACTER(LEN=*) :: outputfile
  LOGICAL :: file_exist
  ! save data
  INQUIRE(file=outputfile, exist=file_exist)
  IF(file_exist) THEN
    OPEN(UNIT=1,FILE=outputfile,status='old', position='append', action='write')
  ELSE
    OPEN(UNIT=1,FILE=outputfile,status='new', action='write')  
    WRITE(1,*) "# x_val, nptot, nmol_a, nmol_b"
  ENDIF
  WRITE(*,*) "===== xi_val, nmol, nmol_a, nmol_b ====="
  WRITE(*,'(F10.5,1x,3(I8,1x))') xi_val, nptot, nmol_a, nmol_b
  WRITE(1,*) xi_val, nptot, nmol_a, nmol_b
  CLOSE(1)
  RETURN
END SUBROUTINE print_nmol

subroutine update_xi_val()
  use coupling_exch
  use sigmas
  use ints
  implicit none
  double precision :: a_frac_compare, ratio
  a_frac_compare = dble(nmol_a)/dble(nptot)
  ratio = a_frac_compare/a_frac - 1.0D0
  if( dabs(ratio) > 0.0001 ) then
    xi_val = xi_val*(1.0D0-ratio)
    write(*,*) " adjust xi_val to ", real(xi_val), "due to difference ", real(ratio)
  endif

  if( xi_val >= 1.0D0) then
    write(*,*) "too big xi_val"
    stop
  else if (xi_val <= 0.0D0) then
    write(*,*) "too small xi_val"
    stop
  endif
  log_xi = dlog(xi_val/(1.0D0-xi_val))
  return
end subroutine update_xi_val

SUBROUTINE save_log(filename)
  use pos
  use ints
  use sigmas
  use inp
  use trans
  use rdf
  use movetype
  IMPLICIT NONE
  CHARACTER(LEN=*) :: filename
  integer :: i
  DOUBLE PRECISION :: VOL, PFC, RGAV, R2AV
  DOUBLE PRECISION, PARAMETER :: pi=3.141592653589793D0

  rgav = 0
  r2av = 0
  !    IF(NAVER.GT.0)THEN
  !      RGAV = RGAV/DBLE(NAVER*nmol_a)
  !      R2AV = R2AV/DBLE(NAVER*nmol_a)
  !    ENDIF
  write(*,*) "save_log:"
  OPEN(UNIT=1,FILE=filename,STATUS='UNKNOWN')
!  PACKING FRACTION
  !VOL = box(1)*box(2)*box(3)
  VOL = box(1)*box(2)*box(3)*(pos_scaling**3)
  PFC = (nmol_a*nmon_a*sigma_a**3+nmol_b*nmon_b*sigma_b**3)*PI/(6.*VOL)
  WRITE(1,*)"PFC,RGAV,R2AV"
  WRITE(1,117) PFC,RGAV,R2AV
117 FORMAT(1X,F6.4,2X,F9.6,2X,F12.6,2X)
!   result
  WRITE(1,*)'RESULTS OF WR particles RUN'
  WRITE(1,*)
  WRITE(1,111) nmol_a,nmon_a,nmol_b,nmon_b,sigma_b,sigma_a,box(1),box(2),box(3)
  WRITE (1,114) NCON, NSKIP
  WRITE (1,115) dlr_a, dlr_b,BSZ,NBIN
  WRITE (1,*)
  WRITE(1,*) "NAME =>  ATTEMPTED  SUCCESSFUL   FRACTION"
  do i=1, nmovetypes
    WRITE (1,'(A5,A,2X,I0,2X,I0,2X,F9.5)') movetype_name(i)," => ",movetype_ntry(i),movetype_nsuccess(i),&
          DBLE(movetype_nsuccess(i))/DBLE(movetype_ntry(i))
  enddo
111 FORMAT(2X, I10, 6X, I10, 6X, 'NUMBER OF A'  / &
           2X, I10, 6X, I10, 6X, 'NUMBER OF B'  / & 
           2X,F6.4, 6X,'DIAMETER OF A' / &
           2X,F6.4, 6X,'DAIMETER OF B' / &       
           4X, E25.18, 2X, &
           4X, E25.18, 2X, &
           4X, E25.18, 2X)
114 FORMAT (1X / 1X, 'NO. TOTAL MOVES           ', 2X, I15 / &
                 1X, 'NO. MOVES PER ACCUMULATION', 2X, I15)
115 FORMAT (1X / 1X, 'DELTA TRANS A   ', 6X, F5.3 / &    
                 1X, 'DELTA TRANS B   ', 6X, F5.3 / &
                 1X, 'BIN SIZES       ', 6X, F5.3 / &
                 1X, 'NUMBER OF BINS  ', 7X, I4)
  close(1)
END SUBROUTINE save_log

subroutine check_shortest_distance()
  use pos
  use cellmap
  use coupling_pres
  use ints
  implicit none
  integer :: i, j, k
  double precision, dimension(3) :: coord_temp
  integer, dimension(2) :: pair
  double precision :: rt, dist

  !dist = minval(box)*pos_scaling/2.0D0
  dist = minval(box)/2.0D0
  do i=1,nptot-1
    do j=i+1,nptot
      if (typs(i) == typs(j)) cycle ! only count different kinds of particles
      coord_temp(1) = x(i) - x(j)
      coord_temp(2) = y(i) - y(j)
      coord_temp(3) = z(i) - z(j)
      do k=1,3
        coord_temp(k) = coord_temp(k) - box(k)*DNINT(coord_temp(k)/box(k))
      enddo
      rt = coord_temp(1)*coord_temp(1)+coord_temp(2)*coord_temp(2)+coord_temp(3)*coord_temp(3)
      rt = dsqrt(rt)*pos_scaling
      rt = dsqrt(rt)
      if(rt < dist) then
        dist = rt
        pair(1) = i
        pair(2) = j
      endif
    enddo
  enddo
  write(*,*) "shortest_dist = ", dist, "with a pair of ", (pair(i),i=1,2)
end subroutine check_shortest_distance
