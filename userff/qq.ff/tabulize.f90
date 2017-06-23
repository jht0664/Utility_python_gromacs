program gen_table
implicit none
!real,parameter :: delr=0.002,rcut=1.0
real, parameter :: delr=0.0005, rcut=10.0
real :: r
integer :: nbins,j

nbins=int( (rcut+1)/delr) + 1
open(unit=6,file='table_QA_QA.xvg')

do j=0,nbins
    r=delr*j
    ! for table_QA_QA
    write(6,*) r, 1/r, 1/(r*r), (-1.0)/(r**6), (-6.0)/(r**7), 1/(r**12), 12/(r**13)
    
    !! for table_QA_QB
    !IF (r <= 1.12246) THEN
    !    write(6,*) r*1.0D0, 1.0D0/r, 1.0D0/(r*r), (-4.0D0)/(r**6), (-24.0D0)/(r**7), 4.0D0/(r**12)+1.0D0, 48.0D0/(r**13) 
    !ELSE
    !    write(6,*) r*1.0D0, 1.0D0/r, 1.0D0/(r*r), 0.0D0, 0.0D0, 0.0D0, 0.0D0
    !ENDIF
    !!  for table_QA_QB
end do

close(6)

End program
