rm *.o *.mod

for f90 in modules.f90 subs.f90 gr12.f90 ranlux.f xdr.f90 io.f90 cell_map.f90 property.f90 overlap.f90 mcrun.f90
do 
	ifort  -xSSE4.2 -O3 -c $f90
done

ifort -qopenmp -xSSE4.2 -O3 -o mcrun.x mcrun.o subs.o gr12.o ranlux.o xdr.o io.o cell_map.o property.o overlap.o modules.o -L/home/htjung/xdr/lib -lxdrfile

