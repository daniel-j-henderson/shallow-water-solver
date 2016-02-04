!driver

program ode_solver
    use equations
    use observer
    
    implicit none
    
    include 'mpif.h'

    real, dimension(:), allocatable :: state, tendency, q1, q2, q3, q4
    real, dimension(:,:,:,:), allocatable :: statechunks, tendchunks, qchunks
    real, dimension(:,:,:), allocatable :: s_chunk, tend_chunk
    integer :: num, system_size, i, j, my_id, ierr, num_procs, m, n
    integer, dimension(MPI_STATUS_SIZE) :: status
    real :: dt, time=0
    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, my_id, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, num_procs, ierr)
    if (num_procs == 1) then !Just to get started, I figured I'd just deal with m and n for specific cases and handle general number of processes after I get that working.
        m = matSize
        n = matSize
    else if (num_procs == 4) then
        m = matSize/2
        n = matSize/2
    else if (num_procs == 8) then
        m = matSize/4
        n = matSize/2
    else if (num_procs == 16) then
        m = matSize/4
        n = matSize/4
    else if (num_procs == 32) then
        m = matSize/8
        n = matSize/4
    else
        write (*,*) "Invalid number of threads, stopping all processes"
        stop
    end if
    
    if (my_id == 0) then
        print *, "How many steps through time?"
        read(*,*) num
    
        print *, "How big is the step size?"
        read(*,*) dt
    
        system_size = get_system_size()
        write (*,*) 'The system size is:', system_size
        allocate(state(system_size))   
    
        call set_initial_state(state)
    
        call observer_init()
        call observer_write(state)
        
    end if
        allocate(s_chunk(m+4, n+4, 3), tend_chunk(m, n, 3))
        
    call mpi_bcast(num, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(dt, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(system_size, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(state, system_size, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    
    !extract initial state chunk sans halo
    if (my_id == 0) then
    	allocate(statechunks(m+4, n+4, 3, num_procs))
    	call disassemble_state(state, statechunks)
    	do j=1, num_procs - 1
            call mpi_send(statechunks(:,:,:,j+1), (n+4)*(m+4)*3,  MPI_REAL, j, 1, MPI_COMM_WORLD, ierr)
        end do
        s_chunk = statechunks(:,:,:,1)
    else
    	call mpi_recv(s_chunk(:,:,:), (n+4)*(m+4)*3, MPI_REAL, 0, 1, MPI_COMM_WORLD, status, ierr)
    end if
    
    
    do i=1, num
    
        call calculate_RK4(s_chunk)
        call observer_write_chunk(tend_chunk, my_id)
        time = time + dt
        
    end do
    
    if (my_id == 0) then
    
        call observer_finalize()
    
        deallocate(state)
        deallocate(tendency)
    	!deallocate more things
        write (*,*), 'A netCDF file called', filename, 'has been created for shallow water model over ', time, 'time'
    else
    	!deallocate some stuff
    end if
    call mpi_finalize(ierr)    
    !clean up all threads here
    
end program ode_solver
        
        