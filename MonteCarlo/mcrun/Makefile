OBJS = subs.o gr12.o ranlux.o xdr.o io.o cell_map.o property.o overlap.o error.o try.o modules.o mcrun.o
CC = ifort -xSSE4.2 -O2
DEBUG = -g -check all -fpe0 -warn -traeback -debug extended
OPENMP = -qopenmp
XDR = -L/home/htjung/xdr/lib -lxdrfile

all: subs.o gr12.o ranlux.o xdr.o io.o cell_map.o property.o overlap.o error.o try.o modules.o mcrun.o

clean:
	rm *.o *.mod *.x *genmod.f90

debug:
	$(CC) = $(CC)$(DEBUG)

subs.o: subs.f90 modules.o
	$(CC) -c subs.f90
gr12.o: gr12.f90 modules.o
	$(CC) -c gr12.f90
ranlux.o: ranlux.f
	$(CC) -c ranlux.f
xdr.o: xdr.f90
	$(CC) -c xdr.f90
io.o: io.f90 modules.o
	$(CC) -c io.f90
cell_map.o: cell_map.f90 modules.o
	$(CC) -c cell_map.f90
property.o: property.f90 modules.o
	$(CC) -c property.f90
overlap.o: overlap.f90 modules.o
	$(CC) -c overlap.f90
error.o: error.f90 modules.o
	$(CC) -c error.f90
try.o: try.f90 modules.o
	$(CC) -c try.f90
modules.o: modules.f90
	$(CC) -c modules.f90
mcrun.o: mcrun.f90 subs.o gr12.o ranlux.o xdr.o io.o cell_map.o property.o overlap.o error.o try.o modules.o
	$(CC) -c mcrun.f90

install: $(OBJS)
	$(CC) $(OPENMP) $(OBJS) -o mcrun.x $(XDR) 


