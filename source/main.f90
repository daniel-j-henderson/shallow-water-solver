!main

program ode_solver
	
	use equations
	use observer
	use solve
	
	implicit none
	
	real, dimension(:), allocatable :: state
	real, dimension(:), allocatable :: tendency
	integer :: num, system_size, i, j
	real :: step_size, time=0
	
	print *, "How many steps through time?"
	read(*,*) num
	
	print *, "How big is the step size?"
	read(*,*) step_size
	
	system_size = get_system_size()
	print *, system_size
	
	allocate(state(system_size))
	allocate(tendency(system_size))	
	
	call set_initial_state(state, wrock)
	
	call observer_init()
	
	do i=1, num
	
		call observer_write(state)
		!Euler's method does not work for this, so I used Michael's RK setup
		!tendency = f(time, state)
		!state = state + step_size * tendency
		!time = time + step_size
		call advance_state(state, time, step_size)
		
		!call observer_write(tendency)
		
	end do
	
	call observer_finalize()
	
	deallocate(state)
	deallocate(tendency)
	
	write (*,*), "A netCDF file called solver_data.nc has been created for shallow water model over ", time, "time"
	
end program ode_solver
			
			
