program gen_table
implicit none
!real,parameter :: delr=0.002,rcut=1.0
real, parameter :: delr=0.0005, rcut=10.0
real :: r
integer :: nbins,j

nbins=int( (rcut+1)/delr) + 1
open(unit=6,file='table_A1_a1.xvg')

do j=0,nbins
    r=delr*j
!    write(6,*) r, 1/r, 1/(r*r),−1/(r**6),−6/(r**7), 1/(r**9), 9/(r**10)
!    write(6,*) r, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0/(r**40), 40.0D0/(r**41)
    write(6,*) r, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0
end do

close(6)

End program
