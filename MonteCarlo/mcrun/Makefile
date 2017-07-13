OBJS = colloid.o subs.o gr11.o gr12.o gr22.o ranlux.o xdr.o modules.o
F77  = ifort -xSSE4.2 -O3  #-fast 

install : $(OBJS) 
	$(F77) -o colloid.x $(OBJS) -L/home/htjung/xdr/lib -lxdrfile

%.o : %.f90 
	$(F77) -c $< 

wrap : 
	$(F77) -c xdr.f90 -L/home/htjung/xdr/lib -lxdrfile

clean :
	rm *.o

