FC = ifort
FFLAGS = -g -check all -traceback -O3
LDFLAGS = 

all: solver

solver: equations.o observer.o solve.o main.o
	$(FC) $(LDFLAGS) equations.o observer.o solve.o main.o -o solver -I${NETCDF}/include -L${NETCDF}/lib -lnetcdf

equations.o: equations.f90
	$(FC) $(FFLAGS) -c equations.f90
	
solve.o: solve.f90
	$(FC) $(FFLAGS) -c solve.f90

observer.o: observer.f90
	$(FC) $(FFLAGS) -c observer.f90 -I$(NETCDF)/include
	
main.o: main.f90 observer.f90 equations.f90
	$(FC) $(FFLAGS) -c main.f90 

clean:
	rm *.o solver *.mod *.nc