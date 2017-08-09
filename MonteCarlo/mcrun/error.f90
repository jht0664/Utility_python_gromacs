subroutine print_over_memory(ierr)
  implicit none
  integer :: ierr
  IF(ierr /= 0) then
    write(*,*) " exceed memory size."
    STOP
  endif
  return
end subroutine print_over_memory
  
