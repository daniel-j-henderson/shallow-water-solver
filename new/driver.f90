!driver

program ode_solver
	
	use equations
	use observer
	use solve
	
	implicit none
	
	real, dimension(:), allocatable :: state, tendency, q1, q2, q3, q4
	real, dimension(:,:,:,:), allocatable :: statechunks, tendchunks, qchunks
	real, dimension(:,:,:), allocatable :: s_chunk, tend_chunk
	integer :: num, system_size, i, j, my_id, ierr, num_procs
	real :: dt, time=0
	
	call mpi_init()
	call mpi_comm_rank(MPI_COMM_WORLD, my_id, ierr)
	call mpi_comm_size(MPI_COMM_WORLD, num_procs, ierr)
	
	if (my_id == 0) then
		print *, "How many steps through time?"
		read(*,*) num
	
		print *, "How big is the step size?"
		read(*,*) dt
	
		system_size = get_system_size()
		print *, system_size
	
		allocate(state(system_size))
		allocate(tendency(system_size))	
	
		call set_initial_state(state)
	
		call observer_init()
		
		allocate(state(system_size), tendency(system_size), q1(system_size), q2(system_size), q3(system_size), q4(system_size))
		allocate(statechunks(m+4, n+4, 3, num_procs), qchunks(m+4, n+4, 3, num_procs))
		allocate(tendchunks(m, n, 3, num_procs))
		
	else
	
		allocate(s_chunk(m, n, 3), tend_chunk(m, n, 3))
		
	end if
	
	call mpi_bcast(num, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
	call mpi_bcast(dt, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
	call mpi_bcast(system_size, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
	
	
	
	if (my_id == 0) then
	
	do i=1, num
	
		if (my_id == 0) then
	
			call observer_write(state)
			
			call disassemble_state(state, statechunks)
			do j=1, num_procs - 1
				call mpi_send(statechunks(:,:,:,j+1), (n+4)*(m+4)*3, MPI_REAL, j, 1, MPI_COMM_WORLD)
			end do
			call f_chunk(statechunks(:,:,1,1), statechunks(:,:,2,1), statechunks(:,:,3,1), &
				         tendchunks(:,:,1,1), tendchunks(:,:,2,1), tendchunks(:,:,3,1))
			do j=1, num_procs - 1
				call mpi_recv(tendchunks(:,:,:,j+1), n*m*3, MPI_REAL, j, 1, MPI_COMM_WORLD, ierr)
			end do
			call reassemble_tend(tendchunks, tendency)
			q1 = dt * tendency
			
			call disassemble_state(q1, qchunks)
			qchunks = .5 * qchunks + statechunks !qchunks now holds state + .5*q1 in chunk form
			do j=1, num_procs - 1
				call mpi_send(qchunks(:,:,:,j+1), (n+4)*(m+4)*3, MPI_REAL, j, 2, MPI_COMM_WORLD)
			end do
			call f_chunk(qchunks(:,:,1,1), q(:,:,2,1), q(:,:,3,1), &
				         tendchunks(:,:,1,1), tendchunks(:,:,2,1), tendchunks(:,:,3,1))
			do j=1, num_procs - 1
				call mpi_recv(tendchunks(:,:,:,j+1), n*m*3, MPI_REAL, j, 2, MPI_COMM_WORLD, ierr)
			end do
			call reassemble_tend(tendchunks, tendency)
			q2 = dt * tendency
			
			call disassemble_state(q2, qchunks)
			qchunks = .5 * qchunks + statechunks !qchunks now holds state + .5*q2 in chunk form
			do j=1, num_procs - 1
				call mpi_send(qchunks(:,:,:,j+1), (n+4)*(m+4)*3, MPI_REAL, j, 3, MPI_COMM_WORLD)
			end do
			call f_chunk(qchunks(:,:,1,1), q(:,:,2,1), q(:,:,3,1), &
				         tendchunks(:,:,1,1), tendchunks(:,:,2,1), tendchunks(:,:,3,1))
			do j=1, num_procs - 1
				call mpi_recv(tendchunks(:,:,:,j+1), n*m*3, MPI_REAL, j, 3, MPI_COMM_WORLD, ierr)
			end do
			call reassemble_tend(tendchunks, tendency)
			q3 = dt * tendency
			
			call disassemble_state(q3, qchunks)
			qchunks = qchunks + statechunks !qchunks now holds state + q3 in chunk form
			do j=1, num_procs - 1
				call mpi_send(qchunks(:,:,:,j+1), (n+4)*(m+4)*3, MPI_REAL, j, 4, MPI_COMM_WORLD)
			end do
			call f_chunk(qchunks(:,:,1,1), q(:,:,2,1), q(:,:,3,1), &
				         tendchunks(:,:,1,1), tendchunks(:,:,2,1), tendchunks(:,:,3,1))
			do j=1, num_procs - 1
				call mpi_recv(tendchunks(:,:,:,j+1), n*m*3, MPI_REAL, j, 4, MPI_COMM_WORLD, ierr)
			end do
			call reassemble_tend(tendchunks, tendency)
			q4 = dt * tendency
			
			state = state + (q1 + 2.0*q2 + 2.0*q3 + q4)/6.0
			
			!PROCESS:
			!disassemble state into chunks
			!send out chunk
			!do chunk
			!receive chunk
			!reassemble
			!calculate q1
			!repeat for q2, q3, q4, state
			!calculate new state
		
		else
			
			call mpi_recv(s_chunk(:,:,:), (n+4)*(m+4)*3, MPI_REAL, my_id, 1, MPI_COMM_WORLD, ierr)
			call f_chunk(s_chunk(:,:,1), s_chunk(:,:,2), s_chunk(:,:,3), &
						 tend_chunk(:,:,1), tend_chunk(:,:,2), tend_chunk(:,:,3))
			call mpi_send(tenc_chunk(:,:,:), n*m*3, MPI_REAL, 0, 1, MPI_COMM_WORLD, ierr)
			
			call mpi_recv(s_chunk(:,:,:), (n+4)*(m+4)*3, MPI_REAL, my_id, 2, MPI_COMM_WORLD, ierr)
			call f_chunk(s_chunk(:,:,1), s_chunk(:,:,2), s_chunk(:,:,3), &
						 tend_chunk(:,:,1), tend_chunk(:,:,2), tend_chunk(:,:,3))
			call mpi_send(tenc_chunk(:,:,:), n*m*3, MPI_REAL, 0, 2, MPI_COMM_WORLD, ierr)
			
			call mpi_recv(s_chunk(:,:,:), (n+4)*(m+4)*3, MPI_REAL, my_id, 3, MPI_COMM_WORLD, ierr)
			call f_chunk(s_chunk(:,:,1), s_chunk(:,:,2), s_chunk(:,:,3), &
						 tend_chunk(:,:,1), tend_chunk(:,:,2), tend_chunk(:,:,3))
			call mpi_send(tenc_chunk(:,:,:), n*m*3, MPI_REAL, 0, 3, MPI_COMM_WORLD, ierr)
			
			call mpi_recv(s_chunk(:,:,:), (n+4)*(m+4)*3, MPI_REAL, my_id, 4, MPI_COMM_WORLD, ierr)
			call f_chunk(s_chunk(:,:,1), s_chunk(:,:,2), s_chunk(:,:,3), &
						 tend_chunk(:,:,1), tend_chunk(:,:,2), tend_chunk(:,:,3))
			call mpi_send(tenc_chunk(:,:,:), n*m*3, MPI_REAL, 0, 4, MPI_COMM_WORLD, ierr)
			!receive chunks
			!do chunks
			!send chunks back
		end if
		time = time + dt
	end do
	
	
	if (my_id == 0) then
	
		call observer_finalize()
	
		deallocate(state)
		deallocate(tendency)
	
		write (*,*), "A netCDF file called solver_data.nc has been created for shallow water model over ", time, "time"
	end if
	
	!clean up all threads here
	
end program ode_solver
			
			