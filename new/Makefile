FC = ifort
FFLAGS = -g -check all -traceback -O3
LDFLAGS = 

all: solver_parallel

solver_parallel: equations.o observer.o driver.o
	$(FC) $(LDFLAGS) equations.o observer.o driver.o -o solver_parallel -I${NETCDF}/include -L${NETCDF}/lib -lnetcdf

equations.o: equations.f90
	$(FC) $(FFLAGS) -c equations.f90
	
observer.o: observer.f90
	$(FC) $(FFLAGS) -c observer.f90 -I$(NETCDF)/include
	
driver.o: driver.f90 observer.f90 equations.f90
	$(FC) $(FFLAGS) -c driver.f90 

clean:
	rm *.o solver *.mod *.nc
