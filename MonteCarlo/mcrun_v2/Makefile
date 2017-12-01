## How to compile
## 1. use ifort
## make;make install
##
## 2. use gfortran
## make gfort;make finstall
##
## Do you want to debug?
## make debug
##

OBJS = movetype.o gr12.o ranlux.o xdr.o io.o cell_map.o property.o overlap.o error.o try.o modules.o mcrun.o
CC = ifort -xSSE4.2 -O3
FCC = gfortran 
#DEBUG = -g -check all -fpe0 -warn -traeback -debug extended
# check all option makes problem for minheap.f90
DEBUG = -g -fpe0 -warn all -traceback -debug extended 
FDEBUG = -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan 
OPENMP = -qopenmp
FOPENMP = -fopenmp
XDR = -L/home/htjung/xdr/lib -lxdrfile

all: movetype.o gr12.o ranlux.o xdr.o io.o cell_map.o property.o overlap.o error.o try.o modules.o mcrun.o

clean:
	rm *.o *.mod *.x *genmod.f90

debug: CC += $(DEBUG)
debug: all

gfort: CC = $(FCC)
gfort: CC += $(FDEBUG)
gfort: OPENMP = $(FOPENMP)
gfort: all

movetype.o: movetype.f90 modules.o
	$(CC) -c movetype.f90
gr12.o: gr12.f90 modules.o
	$(CC) -c gr12.f90
ranlux.o: ranlux.f
	$(CC) -c ranlux.f
xdr.o: xdr.f90
	$(CC) -c xdr.f90
io.o: io.f90 modules.o
	$(CC) $(OPENMP) -c io.f90
cell_map.o: cell_map.f90 modules.o
	$(CC) -c cell_map.f90
property.o: property.f90 modules.o
	$(CC) $(OPENMP) -c property.f90
overlap.o: overlap.f90 modules.o
	$(CC) -c overlap.f90
error.o: error.f90 modules.o
	$(CC) -c error.f90
try.o: try.f90 modules.o
	$(CC) $(OPENMP) -c try.f90
modules.o: modules.f90
	$(CC) -c modules.f90
#minheap.o: minheap.f90 modules.o
#!	$(CC) -c minheap.f90

mcrun.o: mcrun.f90 movetype.o gr12.o ranlux.o xdr.o io.o cell_map.o property.o overlap.o error.o try.o modules.o
	$(CC) -c mcrun.f90

install: $(OBJS)
	$(CC) $(OPENMP)  $(OBJS) -o mcrun.x $(XDR)

finstall: $(OBJS)
	$(FCC) $(FOPENMP) $(OBJS) -o mcrun.x $(XDR)
