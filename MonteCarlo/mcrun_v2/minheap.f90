! save all shortest distance in cells
subroutine minheap_dist_init()
  use cellmap
  use pos
  use minheap
  use ints, only: nptot
  use sigmas
  use coupling_pres
  implicit none
  integer :: i, j, k, ierr, ncell, icell, jcell
  double precision :: rt
  double precision, dimension(3) :: coordi, coordj, coord_temp
  write(*,*) "minheap_dist_init:"
! initialize
  ncell = mcx*mcy*mcz
  allocate(minheap_dist(ncell),stat=ierr)
  if (ierr /= 0) write(*,*) "insufficient space."
  allocate(minheap_address(ncell),stat=ierr)
  if (ierr /= 0) write(*,*) "insufficient space."
  minheap_address = 0 ! minheap_address(icell) = array position of minheap_dist
  dist_criteria = sigma_ab*dexp(dvx)
! update shortest distances every cell
  do i=1,nptot-1
    coordi(1) = x(i) - box(1)*DNINT(x(i)/box(1)-0.50D0) 
    coordi(2) = y(i) - box(2)*DNINT(y(i)/box(2)-0.50D0) 
    coordi(3) = z(i) - box(3)*DNINT(z(i)/box(3)-0.50D0) 
    icell = 1 + INT( coordi(1)*CELLX ) &
            + INT( coordi(2)*CELLY ) *MCX &
            + INT( coordi(3)*CELLZ ) *MCY*MCX
    do j=i+1,nptot
      if (typs(i) == typs(j)) cycle ! only count different kinds of particles
      coord_temp(1) = x(i) - x(j)
      coord_temp(2) = y(i) - y(j)
      coord_temp(3) = z(i) - z(j)
      do k=1,3
        coord_temp(k) = coord_temp(k) - box(k)*DNINT(coord_temp(k)/box(k))
      enddo
      rt = coord_temp(1)*coord_temp(1)+coord_temp(2)*coord_temp(2)+coord_temp(3)*coord_temp(3)
      rt = dsqrt(rt)
      if(rt > dist_criteria) cycle
      call minheap_manage_data(rt, icell)
      coordj(1) = x(j) - box(1)*DNINT(x(j)/box(1)-0.50D0) 
      coordj(2) = y(j) - box(2)*DNINT(y(j)/box(2)-0.50D0) 
      coordj(3) = z(j) - box(3)*DNINT(z(j)/box(3)-0.50D0) 
      jcell = 1 + INT( coordj(1)*CELLX ) &
            + INT( coordj(2)*CELLY ) *MCX &
            + INT( coordj(3)*CELLZ ) *MCY*MCX
      if( icell /= jcell ) call minheap_manage_data(rt, jcell)
    enddo
  enddo
  write(*,*) " total minheap_dist #elements = ", minheap_size
  return
end subroutine minheap_dist_init

subroutine minheap_manage_data(data, icell)
  use minheap
  integer :: icell, inode
  double precision :: data
!  write(*,*) "minheap_manage_data:"
  if( data > dist_criteria) return ! reduce memory usage
  inode = minheap_address(icell)
! if no exist in minheap_address
  if( inode == 0) then
!    write(*,*) "insert!"
    call minheap_dist_insert(data,icell)
! if node exist in minheap_address and data is smaller, then update
  else if( data < minheap_dist(inode)%distance ) then
    call minheap_dist_update(data,inode)
  endif
end subroutine minheap_manage_data

subroutine minheap_print(text,text2)
  use minheap
  implicit none
  integer :: i, text2
  CHARACTER(LEN=*) :: text
!  write(*,*) "minheap_print:"
  OPEN(UNIT=3,FILE='heap.array',STATUS='UNKNOWN',position='append',action='write') ! append if already exist
  WRITE(3,'(a,1x,a,1x,I0,1x,a)') '#', text, text2, ' frame'
  do i=1,minheap_size
    write(3,*) i, minheap_dist(i)
  enddo
  return
end subroutine minheap_print

subroutine minheap_size_print()
  use minheap
  implicit none
  integer :: i,j
  j = 0
  do i=1,size(minheap_address)
    if(minheap_address(i) /= 0) then
!      write(*,*) i, minheap_address(i)
      j = j + 1
    endif
  enddo
  if(minheap_size /= j) then
    write(*,*) minheap_size, "should be same ", j
    stop
  endif
  return
end subroutine minheap_size_print


subroutine short_dist_update()
  use minheap
  use coupling_pres, only: expd
  implicit none
  integer :: i
  do i=1, minheap_size
    minheap_dist(i)%distance = minheap_dist(i)%distance*expd
  enddo
  return
end subroutine short_dist_update

! update minheap if imol th particle succeed in translation movement
subroutine minheap_update_by_tran(imol)
  use pos
  use minheap
  use cellmap
  implicit none
  integer :: imol, i, ii, j, jj, k, icell, jcell
  double precision :: rt
  double precision, dimension(3) :: coordi, coordj, coord_temp
!  write(*,*) "minheap_update_by_tran:"
  niter = niter + 1
  jnear(niter) = imol
! update shorter distance in neighbor cell
  do i=1,niter-1
    ii = jnear(i)
    coordi(1) = x(ii) - box(1)*DNINT(x(ii)/box(1)-0.50D0) 
    coordi(2) = y(ii) - box(2)*DNINT(y(ii)/box(2)-0.50D0) 
    coordi(3) = z(ii) - box(3)*DNINT(z(ii)/box(3)-0.50D0) 
    icell = 1 + INT( coordi(1)*CELLX ) &
            + INT( coordi(2)*CELLY ) *MCX &
            + INT( coordi(3)*CELLZ ) *MCY*MCX
    do j=i+1,niter
      if(i == j) cycle ! avoid self pair
      jj = jnear(j)
      if (typs(ii) == typs(jj)) cycle ! only count different kinds of particles
      coord_temp(1) = x(ii) - x(jj)
      coord_temp(2) = y(ii) - y(jj)
      coord_temp(3) = z(ii) - z(jj)
      do k=1,3
        coord_temp(k) = coord_temp(k) - box(k)*DNINT(coord_temp(k)/box(k))
      enddo
      rt = coord_temp(1)*coord_temp(1)+coord_temp(2)*coord_temp(2)+coord_temp(3)*coord_temp(3)
      rt = dsqrt(rt)
      if(rt > dist_criteria) cycle
      call minheap_manage_data(rt, icell)
      coordj(1) = x(j) - box(1)*DNINT(x(j)/box(1)-0.50D0) 
      coordj(2) = y(j) - box(2)*DNINT(y(j)/box(2)-0.50D0) 
      coordj(3) = z(j) - box(3)*DNINT(z(j)/box(3)-0.50D0) 
      jcell = 1 + INT( coordj(1)*CELLX ) &
            + INT( coordj(2)*CELLY ) *MCX &
            + INT( coordj(3)*CELLZ ) *MCY*MCX
      if( icell /= jcell ) call minheap_manage_data(rt, jcell)
    enddo
  enddo
  return
end subroutine minheap_update_by_tran

subroutine minheap_dist_insert(key, icell)
  use minheap
  implicit none
  integer :: icell
  double precision :: key ! value will be stored in minheap_dist
! insert a node in the fartest left location of the lowest level
  minheap_size = minheap_size + 1
  minheap_dist(minheap_size)%distance = key
  minheap_dist(minheap_size)%index_cell = icell
  minheap_address(icell) = minheap_size
! Compare the values of the inserted node with its parent node
  call minheap_dist_compare_parent(minheap_size)
  return
end subroutine minheap_dist_insert

recursive subroutine minheap_dist_compare_parent(child_node)
  use minheap
  implicit none
  integer :: parent_node, child_node
  if( child_node == 1 ) then
    return ! root node does not have parent node, so finish comparing
  endif
! if parent node exists,
  parent_node = floor(real(child_node)/2.0)
  if( minheap_dist(child_node)%distance < minheap_dist(parent_node)%distance ) then 
    call minheap_dist_swap(parent_node,child_node) ! swap
    call minheap_dist_compare_parent(parent_node) ! call compare with grandparent
  endif
  return  
end subroutine minheap_dist_compare_parent

subroutine minheap_dist_swap(inode,jnode)
  use minheap
  implicit none
  integer :: inode, jnode, temp_icell, icell, jcell
  double precision :: temp_dist
! swap
  temp_dist = minheap_dist(jnode)%distance
  temp_icell = minheap_dist(jnode)%index_cell
  minheap_dist(jnode)%distance    = minheap_dist(inode)%distance
  minheap_dist(jnode)%index_cell  = minheap_dist(inode)%index_cell
  minheap_dist(inode)%distance    = temp_dist
  minheap_dist(inode)%index_cell  = temp_icell
! update address
  jcell = minheap_dist(jnode)%index_cell
  icell = minheap_dist(inode)%index_cell
  minheap_address(jcell) = jnode
  minheap_address(icell) = inode
  return
end subroutine minheap_dist_swap

recursive subroutine minheap_dist_update(key,inode)
  use minheap
  implicit none
  integer :: inode
  double precision :: key
  minheap_dist(inode)%distance = key
  call minheap_dist_compare_parent(inode)
  call minheap_dist_compare_child(inode)
end subroutine minheap_dist_update

recursive subroutine minheap_dist_compare_child(parent_node)
  use minheap
  implicit none
  integer :: parent_node, child_node1, child_node2
  double precision :: child_dist1, child_dist2, parent_dist
  child_node1 = parent_node*2
  child_node2 = parent_node*2+1
! if it does not have children node
  if( child_node1 > minheap_size) then 
    return 
! if only one children node exists,
  else if (child_node1 == minheap_size) then
    if( minheap_dist(parent_node)%distance > minheap_dist(child_node1)%distance ) then
      call minheap_print("before swap with child1 and parent ",parent_node)
      call minheap_dist_swap(parent_node,child_node1) ! swap
      call minheap_print("after swap with child1 and parent ",parent_node)
    endif
! if two childlen exist
  else 
    child_dist1 = minheap_dist(child_node1)%distance
    child_dist2 = minheap_dist(child_node2)%distance
    parent_dist = minheap_dist(parent_node)%distance
    if( parent_dist < child_dist1 .and. parent_dist < child_dist2 ) then
      return
    else if ( child_dist1 < parent_dist .and. child_dist2 > parent_dist) then
      call minheap_dist_swap(parent_node,child_node1)
      call minheap_dist_compare_child(child_node1)
    else if ( child_dist1 > parent_dist .and. child_dist2 < parent_dist) then
      call minheap_dist_swap(parent_node,child_node2)
      call minheap_dist_compare_child(child_node2)
    else if ( child_dist1 < parent_dist .and. child_dist2 < parent_dist) then
      if (child_dist1 <= child_dist2) then
        call minheap_dist_swap(parent_node,child_node1)
        call minheap_dist_compare_child(child_node1)
      else
        call minheap_dist_swap(parent_node,child_node2)
        call minheap_dist_compare_child(child_node2)
      endif
    endif
  endif
  return  
end subroutine minheap_dist_compare_child

! http://www.mathcs.emory.edu/~cheung/Courses/171/Syllabus/9-BinTree/heap-delete.html
subroutine minheap_dist_delete(inode)
  use minheap, only: minheap_address, minheap_size, minheap_dist
  implicit none
  integer, intent(in) :: inode
  integer :: temp
!  write(*,*) "minheap_dist_delete: ", inode, " th node in minheap"
  temp = minheap_dist(inode)%index_cell
  minheap_address(temp) = 0
  if(inode == minheap_size) then
    minheap_size = minheap_size - 1
    return ! just delete for the last element 
  else
    minheap_dist(inode)%distance = minheap_dist(minheap_size)%distance
    minheap_dist(inode)%index_cell = minheap_dist(minheap_size)%index_cell
    minheap_size = minheap_size - 1
    temp = minheap_dist(inode)%index_cell
    minheap_address(temp) = inode
    call minheap_dist_compare_parent(inode)
    call minheap_dist_compare_child(inode)
  endif
  return
end subroutine minheap_dist_delete

subroutine minheap_dist_exit()
  use minheap
  implicit none
  integer :: ierr
  deallocate(minheap_dist,stat=ierr)
  deallocate(minheap_address,stat=ierr)
end subroutine minheap_dist_exit


subroutine add_jnear_oldcell(imol)
  use cellmap, only: jnear, niter
  use minheap
  implicit none
  integer :: imol, id, icell, i, j, k, ierr, node
  integer, dimension(:), allocatable :: res  
  call get_icell_imol(imol,icell)
  node = minheap_address(icell)
!  write(*,*) "add_jnear_oldcell:",icell, node
  if (node /= 0) then
    call minheap_dist_delete(node)
  endif
  id = (icell-1)*27
  call build_jnear(imol,id,"a")
  allocate(res(niter),stat=ierr)
  k = 1
  res(1) = jnear(1)
! remove duplicates of particles
  outer: do i=2,niter
    do j=1,k
      if( res(i) == jnear(j)) then
        cycle outer
      endif
    enddo
    k = k+1
    res(k) = jnear(i)
  end do outer
  niter = k
! we got a res array(1:k)
  do i=1,k
    jnear(i) = res(i)
  enddo
  deallocate(res,stat=ierr)
  return
end subroutine add_jnear_oldcell

! get map_id for the coord argumnet
subroutine get_icell_xyz(coord, icell)
  use pos, only: box
  use cellmap, only: cellx, celly, cellz, mcx, mcy
  implicit none
  double precision, dimension(3) :: coord
  integer :: i, icell
  do i=1,3
    coord(i) = coord(i) - box(i)*DNINT(coord(i)/box(i)-0.50D0) 
  enddo
  icell = 1 + INT( coord(1)*CELLX ) &
            + INT( coord(2)*CELLY ) *MCX &
            + INT( coord(3)*CELLZ ) *MCY*MCX
  return
end subroutine get_icell_xyz

subroutine get_icell_imol(imol,icell)
  use pos
  use cellmap, only: cellx, celly, cellz, mcx, mcy
  implicit none
  double precision, dimension(3) :: coord
  integer :: imol, icell
  coord(1) = x(imol) - box(1)*DNINT(x(imol)/box(1)-0.50D0) 
  coord(2) = y(imol) - box(2)*DNINT(y(imol)/box(2)-0.50D0) 
  coord(3) = z(imol) - box(3)*DNINT(z(imol)/box(3)-0.50D0) 
  icell = 1 + INT( coord(1)*CELLX ) &
            + INT( coord(2)*CELLY ) *MCX &
            + INT( coord(3)*CELLZ ) *MCY*MCX
  return
end subroutine get_icell_imol

subroutine get_icell_try(icell)
  use pos
  use cellmap, only: cellx, celly, cellz, mcx, mcy
  implicit none
  double precision, dimension(3) :: coord
  integer :: i, icell
  do i=1,3
    coord(i) = coord_try(i) - box(i)*DNINT(coord_try(i)/box(i)-0.50D0) 
  enddo
  icell = 1 + INT( coord(1)*CELLX ) &
            + INT( coord(2)*CELLY ) *MCX &
            + INT( coord(3)*CELLZ ) *MCY*MCX
  return
end subroutine get_icell_try

subroutine build_jnear(imol,map_id,action)
  use cellmap, only: map, lead, jnear, list, niter
  implicit none
  integer, intent(in) :: map_id, imol ! index of a cell, index of particle to avoid couting
  integer :: nc, im, jcell
  CHARACTER(LEN=1) :: action
  if(action == 'o') then
    NITER = 0 ! overwrite
  else if (action == 'a') then
    ! append (no change niter)
  else
    write(*,*) " wrong argument in build_jnear"
    stop
  endif
  DO NC = 1, 27
    JCELL = MAP(map_id+NC)
    IF(JCELL.EQ.0) CONTINUE
    IM = LEAD(JCELL)
    DO WHILE (IM.NE.0)
      ! skip the identical polymer when saving neighbor particles, but count duplicates
      IF(IM .NE. imol) THEN
        NITER = NITER + 1
        JNEAR(NITER) = IM
      ENDIF
      IM = LIST(IM)
    ENDDO
  ENDDO    
  return
end subroutine build_jnear

subroutine cal_dist_xyz(coord1, coord2, rt)
  use pos, only: box
  implicit none
  double precision, dimension(3) :: coord1, coord2, coord_temp
  integer :: i
  double precision :: rt
  do i=1,3
    coord_temp(i) = coord1(i) - coord2(i)
    coord_temp(i) = coord_temp(i) - box(i)*DNINT(coord_temp(i)/box(i))
  enddo
  rt = coord_temp(1)*coord_temp(1)+coord_temp(2)*coord_temp(2)+coord_temp(3)*coord_temp(3)
  return
end subroutine cal_dist_xyz

subroutine cal_dist_imolxyz(imol, coord2, rt)
  use pos, only: x,y,z,box
  implicit none
  double precision, dimension(3) :: coord2, coord_temp
  integer :: i
  integer :: imol
  double precision :: rt
  coord_temp(1) = x(imol) - coord2(1)
  coord_temp(2) = y(imol) - coord2(2)
  coord_temp(3) = z(imol) - coord2(3)
  do i=1,3
    coord_temp(i) = coord_temp(i) - box(i)*DNINT(coord_temp(i)/box(i))
  enddo
  rt = coord_temp(1)*coord_temp(1)+coord_temp(2)*coord_temp(2)+coord_temp(3)*coord_temp(3)
  return
end subroutine cal_dist_imolxyz

subroutine cal_dist_imol(imol, jmol, rt)
  use pos, only: x,y,z,box
  implicit none
  double precision, dimension(3) :: coord_temp
  integer :: i, imol, jmol
  double precision :: rt
  coord_temp(1) = x(imol) - x(jmol)
  coord_temp(2) = y(imol) - y(jmol)
  coord_temp(3) = z(imol) - z(jmol)
  do i=1,3
    coord_temp(i) = coord_temp(i) - box(i)*DNINT(coord_temp(i)/box(i))
  enddo
  rt = coord_temp(1)*coord_temp(1)+coord_temp(2)*coord_temp(2)+coord_temp(3)*coord_temp(3)
  return
end subroutine cal_dist_imol