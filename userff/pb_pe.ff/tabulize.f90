program gen_table
implicit none
!real,parameter :: delr=0.002,rcut=1.0
DOUBLE PRECISION, parameter :: delr=0.0005D0, rcut=3.0D0
DOUBLE PRECISION :: r
integer :: nbins,j

nbins=int( (rcut+1)/delr) + 1
open(unit=6,file='table_A_B.xvg')

do j=0,nbins
    r=delr*j
    
    ! for table_A_A
    IF (r < 0.96D0) THEN
        write(6,*) r*1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 10**5, 10**5
    else if (r <= 1.57D0) then
        !write(6,*) r*1.0D0, 0.0D0, 0.0D0, (-4.0D0)/(r**6), (-24.0D0)/(r**7), &
        !        4.0D0/(r**12)-(4.0D0/(1.57D0**12)-4.0D0/(1.57D0**6)), 48.0D0/(r**13)  ! A-A or B-B
        write(6,*) r*1.0D0, 0.0D0, 0.0D0, (-4.0D0)/((r/1.03D0)**6), (-24.0D0)/((r/1.03D0)**7), &
                4.0D0/((r/1.03D0)**12)-(4.0D0/((1.57D0/1.03D0)**12)-4.0D0/((1.57D0/1.03D0)**6)),&
                48.0D0/((r/1.03D0)**13)   ! A-B
    ELSE
        write(6,*) r*1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0
    ENDIF
    
end do

close(6)

End program
