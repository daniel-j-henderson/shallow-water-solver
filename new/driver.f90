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
        allocate(state(system_size), tendency(system_size))   
    
        call set_initial_state(state)
    
        call observer_init()
        allocate(q1(system_size), q2(system_size), q3(system_size), q4(system_size))
        allocate(statechunks(m+4, n+4, 3, num_procs), qchunks(m+4, n+4, 3, num_procs))
        allocate(tendchunks(m, n, 3, num_procs))
        
    end if
        allocate(s_chunk(m+4, n+4, 3), tend_chunk(m, n, 3))
        
    call mpi_bcast(num, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(dt, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(system_size, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
    
    
    
    do i=1, num
    
        if (my_id == 0) then
            call observer_write(state)
            
            call disassemble_state(state, statechunks) !turns the state vector into an array of state chins to be sent out to each thread
            do j=1, num_procs - 1
                call mpi_send(statechunks(:,:,:,j+1), (n+4)*(m+4)*3,  MPI_REAL, j, 1, MPI_COMM_WORLD, ierr)
            end do
            call f_chunk(statechunks(:,:,1,1), statechunks(:,:,2,1), statechunks(:,:,3,1), &
                         tendchunks(:,:,1,1), tendchunks(:,:,2,1), tendchunks(:,:,3,1))
            do j=1, num_procs - 1
                call mpi_recv(tendchunks(:,:,:,j+1), n*m*3, MPI_REAL, j, 1, MPI_COMM_WORLD, status, ierr)
            end do
            call reassemble_tend(tendchunks, tendency)
            
            q1 = dt * tendency !RK steps
            
            call disassemble_state(q1, qchunks) 
            qchunks = .5 * qchunks + statechunks !qchunks now holds state + .5*q1 in chunk form
            do j=1, num_procs - 1
                call mpi_send(qchunks(:,:,:,j+1), (n+4)*(m+4)*3, MPI_REAL, j, 2, MPI_COMM_WORLD, ierr)
            end do
            call f_chunk(qchunks(:,:,1,1), qchunks(:,:,2,1), qchunks(:,:,3,1), &
                         tendchunks(:,:,1,1), tendchunks(:,:,2,1), tendchunks(:,:,3,1))
            do j=1, num_procs - 1
                call mpi_recv(tendchunks(:,:,:,j+1), n*m*3, MPI_REAL, j, 2, MPI_COMM_WORLD, status, ierr)
            end do
            call reassemble_tend(tendchunks, tendency)
            
            q2 = dt * tendency
            
            call disassemble_state(q2, qchunks)
            qchunks = .5 * qchunks + statechunks !qchunks now holds state + .5*q2 in chunk form
            do j=1, num_procs - 1
                call mpi_send(qchunks(:,:,:,j+1), (n+4)*(m+4)*3, MPI_REAL, j, 3, MPI_COMM_WORLD, ierr)
            end do
            call f_chunk(qchunks(:,:,1,1), qchunks(:,:,2,1), qchunks(:,:,3,1), &
                         tendchunks(:,:,1,1), tendchunks(:,:,2,1), tendchunks(:,:,3,1))
            do j=1, num_procs - 1
                call mpi_recv(tendchunks(:,:,:,j+1), n*m*3, MPI_REAL, j, 3, MPI_COMM_WORLD, status, ierr)
            end do
            call reassemble_tend(tendchunks, tendency)
            
            q3 = dt * tendency
            
            call disassemble_state(q3, qchunks)
            qchunks = qchunks + statechunks !qchunks now holds state + q3 in chunk form
            do j=1, num_procs - 1
                call mpi_send(qchunks(:,:,:,j+1), (n+4)*(m+4)*3, MPI_REAL, j, 4, MPI_COMM_WORLD, ierr)
            end do
            call f_chunk(qchunks(:,:,1,1), qchunks(:,:,2,1), qchunks(:,:,3,1), &
                         tendchunks(:,:,1,1), tendchunks(:,:,2,1), tendchunks(:,:,3,1))
            do j=1, num_procs - 1
                call mpi_recv(tendchunks(:,:,:,j+1), n*m*3, MPI_REAL, j, 4, MPI_COMM_WORLD, status, ierr)
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
            call mpi_recv(s_chunk(:,:,:), (n+4)*(m+4)*3, MPI_REAL, 0, 1, MPI_COMM_WORLD, status, ierr)
            call f_chunk(s_chunk(:,:,1), s_chunk(:,:,2), s_chunk(:,:,3), &
                         tend_chunk(:,:,1), tend_chunk(:,:,2), tend_chunk(:,:,3))
            call mpi_send(tend_chunk(:,:,:), n*m*3, MPI_REAL, 0, 1, MPI_COMM_WORLD, status, ierr)
            
            call mpi_recv(s_chunk(:,:,:), (n+4)*(m+4)*3, MPI_REAL, 0, 2, MPI_COMM_WORLD, status, ierr)
            call f_chunk(s_chunk(:,:,1), s_chunk(:,:,2), s_chunk(:,:,3), &
                         tend_chunk(:,:,1), tend_chunk(:,:,2), tend_chunk(:,:,3))
            call mpi_send(tend_chunk(:,:,:), n*m*3, MPI_REAL, 0, 2, MPI_COMM_WORLD, ierr)
            
            call mpi_recv(s_chunk(:,:,:), (n+4)*(m+4)*3, MPI_REAL, 0, 3, MPI_COMM_WORLD, status, ierr)
            call f_chunk(s_chunk(:,:,1), s_chunk(:,:,2), s_chunk(:,:,3), &
                         tend_chunk(:,:,1), tend_chunk(:,:,2), tend_chunk(:,:,3))
            call mpi_send(tend_chunk(:,:,:), n*m*3, MPI_REAL, 0, 3, MPI_COMM_WORLD, ierr)
            
            call mpi_recv(s_chunk(:,:,:), (n+4)*(m+4)*3, MPI_REAL, 0, 4, MPI_COMM_WORLD, status, ierr)
            call f_chunk(s_chunk(:,:,1), s_chunk(:,:,2), s_chunk(:,:,3), &
                         tend_chunk(:,:,1), tend_chunk(:,:,2), tend_chunk(:,:,3))
            call mpi_send(tend_chunk(:,:,:), n*m*3, MPI_REAL, 0, 4, MPI_COMM_WORLD, ierr)
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
    
        write (*,*), 'A netCDF file called', filename, 'has been created for shallow water model over ', time, 'time'
    end if
    call mpi_finalize(ierr)    
    !clean up all threads here
    
end program ode_solver
            
            
