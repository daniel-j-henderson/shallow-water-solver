!equations


module equations
    
    implicit none    
    private
    integer :: i, j, im, ip, jm, jp 
    real :: dx = 100 ! dx and dy are required for the numerical differentiation I used to get the divergences/gradients
    real :: dy= 100 ! for arbitrary dx, step_size should be about 1/10-1/5 of dx to yield results that don't blow up. 
    
    public :: get_system_size, set_initial_state, f_chunk, disassemble_state, reassemble_tend, calculate_RK4
    integer, parameter, public :: matSize = 100 !set the dimensions of the grid of fluid. Interestingly, it seg faults for 400+
    integer, public :: arrSize, my_id, num_procs, ierr
    real, public :: dt
    integer :: count = 0
    
    
    contains
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    integer function get_system_size()
    
        arrSize = matSize**2
        get_system_size = 3*arrSize
        return
        
    end function get_system_size
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine calculate_RK4(s_chunk, tend_chunk)
		implicit none
        include 'mpif.h'
		real, dimension(:,:,:), intent(inout) :: s_chunk !(m+4)x(n+4)x3
		real, dimension(:,:,:), intent(inout) :: tend_chunk !mxnx3, always h:u:v
		real, dimension(:,:,:), allocatable :: qchunk1, qchunk2, qchunk3, qchunk4
		integer, dimension(3) :: stuff 
		integer :: m, n, dest, source, width, height
		integer :: request1, request2, request3, request4
		integer :: status
        count = count + 1
        stuff = shape(s_chunk)
        m = stuff(1)
        n = stuff(2)
        
        allocate(qchunk1(m, n, 3), qchunk2(m, n, 3), qchunk3(m, n, 3), qchunk4(m, n, 3))
        m = m-4
        n = n-4
        width = matSize / m
        height = matSize / n
        call f_chunk(s_chunk(:,:,1), s_chunk(:,:,2), s_chunk(:,:,3), tend_chunk(:,:,1), tend_chunk(:,:,2), tend_chunk(:,:,3))
        write (*,*) 'count = ',count
        dest = mod(my_id, width) + mod((my_id/width + 1), height)*width !send top halo
        write(*,*) my_id, 'sending top halo to',dest
        call mpi_isend(tend_chunk(1:m, n-1:n, :), 2*m, MPI_REAL, dest, 1, MPI_COMM_WORLD, request1, ierr)
        
        dest = mod(my_id, width) + mod((my_id/width - 1) + height, height)*width !send bottom halo
        write(*,*) my_id, 'sending bottom halo to',dest
        call mpi_isend(tend_chunk(1:m, 1:2, :), 2*m, MPI_REAL, dest, 2, MPI_COMM_WORLD, request2, ierr)
        
        dest = mod(my_id-1 + width, width) + (my_id/width)*width !send left halo
        write(*,*) my_id, 'sending left halo to',dest
        call mpi_isend(tend_chunk(1:2, 1:n, :), 2*n, MPI_REAL, dest, 3, MPI_COMM_WORLD, request3, ierr)
        
        dest = mod(my_id+1, width) + (my_id/width)*width !send right halo
        write(*,*) my_id, 'sending right halo to',dest
        call mpi_isend(tend_chunk(m-1:m, 1:n, :), 2*n, MPI_REAL, dest, 4, MPI_COMM_WORLD, request4, ierr)
        
        !recv top halo
        source = mod(my_id, width) + mod((my_id/width + 1), height)*width
        write(*,*) my_id, 'receiving top halo from',source
        call mpi_recv(qchunk1(3:m+2, n+3:n+4, :), 2*m, MPI_REAL, source, 1, MPI_COMM_WORLD, status, ierr)
        !recv bottom halo
        source = mod(my_id, width) + mod((my_id/width - 1) + height, height)*width
        write(*,*) my_id, 'receiving bottom halo from',source
        call mpi_recv(qchunk1(3:m+2, 1:2, :), 2*m, MPI_REAL, source, 2, MPI_COMM_WORLD, status, ierr)
        !recv left halo
        source = mod(my_id-1 + width, width) + (my_id/width)*width
        write(*,*) my_id, 'receiving left halo from',source
        call mpi_recv(qchunk1(1:2, 3:n+2, :), 2*n, MPI_REAL, source, 3, MPI_COMM_WORLD, status, ierr)
        !recv right halo
        source = mod(my_id+1, width) + (my_id/width)*width
        write(*,*) my_id, 'receiving right halo from',source
        call mpi_recv(qchunk1(m+3:m+4, 3:n+2, :), 2*n, MPI_REAL, source, 4, MPI_COMM_WORLD, status, ierr)
        print *, 'A'
        !send tend_chunk halo bits
        !receive halo bits into qchunk1
        qchunk1(3:m+2, 3:n+2, :) = tend_chunk(:,:,:)
        qchunk1 = qchunk1*dt
        qchunk2 = .5 * qchunk1 + s_chunk
        !ensure tend_chunks has been received
        call mpi_wait(request1, status, ierr)
        call mpi_wait(request2, status, ierr)
        call mpi_wait(request3, status, ierr)
        call mpi_wait(request4, status, ierr)
        print *, 'B'
        call f_chunk(qchunk2(:,:,1), qchunk2(:,:,2), qchunk2(:,:,3), tend_chunk(:,:,1), tend_chunk(:,:,2), tend_chunk(:,:,3))
        !send tend_chunk halo bits
        !receive halo bits into qchunk2
        dest = mod(my_id, width) + mod((my_id/width + 1), height)*width !send top halo
        call mpi_isend(tend_chunk(1:m, n-1:n, :), 2*m, MPI_REAL, dest, 1, MPI_COMM_WORLD, request1, ierr)
        
        dest = mod(my_id, width) + mod((my_id/width - 1) + height, height)*width !send bottom halo
        call mpi_isend(tend_chunk(1:m, 1:2, :), 2*m, MPI_REAL, dest, 2, MPI_COMM_WORLD, request2, ierr)
        
        dest = mod(my_id-1 + width, width) + (my_id/width)*width !send left halo
        call mpi_isend(tend_chunk(1:2, 1:n, :), 2*n, MPI_REAL, dest, 3, MPI_COMM_WORLD, request3, ierr)
        
        dest = mod(my_id+1, width) + (my_id/width)*width !send right halo
        call mpi_isend(tend_chunk(m-1:m, 1:n, :), 2*n, MPI_REAL, dest, 4, MPI_COMM_WORLD, request4, ierr)
        
        !recv top halo
        source = mod(my_id, width) + mod((my_id/width + 1), height)*width
        call mpi_recv(qchunk2(3:m+2, n+3:n+4, :), 2*m, MPI_REAL, source, 1, MPI_COMM_WORLD, status, ierr)
        !recv bottom halo
        source = mod(my_id, width) + mod((my_id/width - 1) + height, height)*width
        call mpi_recv(qchunk2(3:m+2, 1:2, :), 2*m, MPI_REAL, source, 2, MPI_COMM_WORLD, status, ierr)
        !recv left halo
        source = mod(my_id-1 + width, width) + (my_id/width)*width
        call mpi_recv(qchunk2(1:2, 3:n+2, :), 2*n, MPI_REAL, source, 3, MPI_COMM_WORLD, status, ierr)
        !recv right halo
        source = mod(my_id+1, width) + (my_id/width)*width
        call mpi_recv(qchunk2(m+3:m+4, 3:n+2, :), 2*n, MPI_REAL, source, 4, MPI_COMM_WORLD, status, ierr)
        
        qchunk2(3:m+2, 3:n+2, :) = tend_chunk(:,:,:)
        qchunk2 = qchunk2*dt
        qchunk3 = qchunk2 * .5 + s_chunk
        call mpi_wait(request1, status, ierr)
        call mpi_wait(request2, status, ierr)
        call mpi_wait(request3, status, ierr)
        call mpi_wait(request4, status, ierr)
        
        call f_chunk(qchunk3(:,:,1), qchunk3(:,:,2), qchunk3(:,:,3), tend_chunk(:,:,1), tend_chunk(:,:,2), tend_chunk(:,:,3))
		!send tend_chunk halo bits
        !receive halo bits into qchunk2
        dest = mod(my_id, width) + mod((my_id/width + 1), height)*width !send top halo
        call mpi_isend(tend_chunk(1:m, n-1:n, :), 2*m, MPI_REAL, dest, 1, MPI_COMM_WORLD, request1, ierr)
        
        dest = mod(my_id, width) + mod((my_id/width - 1) + height, height)*width !send bottom halo
        call mpi_isend(tend_chunk(1:m, 1:2, :), 2*m, MPI_REAL, dest, 2, MPI_COMM_WORLD, request2, ierr)
        
        dest = mod(my_id-1 + width, width) + (my_id/width)*width !send left halo
        call mpi_isend(tend_chunk(1:2, 1:n, :), 2*n, MPI_REAL, dest, 3, MPI_COMM_WORLD, request3, ierr)
        
        dest = mod(my_id+1, width) + (my_id/width)*width !send right halo
        call mpi_isend(tend_chunk(m-1:m, 1:n, :), 2*n, MPI_REAL, dest, 4, MPI_COMM_WORLD, request4, ierr)
        
        !recv top halo
        source = mod(my_id, width) + mod((my_id/width + 1), height)*width
        call mpi_recv(qchunk3(3:m+2, n+3:n+4, :), 2*m, MPI_REAL, source, 1, MPI_COMM_WORLD, status, ierr)
        !recv bottom halo
        source = mod(my_id, width) + mod((my_id/width - 1) + height, height)*width
        call mpi_recv(qchunk3(3:m+2, 1:2, :), 2*m, MPI_REAL, source, 2, MPI_COMM_WORLD, status, ierr)
        !recv left halo
        source = mod(my_id-1 + width, width) + (my_id/width)*width
        call mpi_recv(qchunk3(1:2, 3:n+2, :), 2*n, MPI_REAL, source, 3, MPI_COMM_WORLD, status, ierr)
        !recv right halo
        source = mod(my_id+1, width) + (my_id/width)*width
        call mpi_recv(qchunk3(m+3:m+4, 3:n+2, :), 2*n, MPI_REAL, source, 4, MPI_COMM_WORLD, status, ierr)
        
        qchunk3(3:m+2, 3:n+2, :) = tend_chunk(:,:,:)
        qchunk3 = qchunk3*dt
        qchunk4 = qchunk3 + s_chunk
        call mpi_wait(request1, status, ierr)
        call mpi_wait(request2, status, ierr)
        call mpi_wait(request3, status, ierr)
        call mpi_wait(request4, status, ierr)
        
        call f_chunk(qchunk4(:,:,1), qchunk4(:,:,2), qchunk4(:,:,3), tend_chunk(:,:,1), tend_chunk(:,:,2), tend_chunk(:,:,3))
        !send tend_chunk halo bits
        !receive halo bits into qchunk2
        dest = mod(my_id, width) + mod((my_id/width + 1), height)*width!send top halo
        call mpi_isend(tend_chunk(1:m, n-1:n, :), 2*m, MPI_REAL, dest, 1, MPI_COMM_WORLD, request1, ierr)
        
        dest = mod(my_id, width) + mod((my_id/width - 1) + height, height)*width !send bottom halo
        call mpi_isend(tend_chunk(1:m, 1:2, :), 2*m, MPI_REAL, dest, 2, MPI_COMM_WORLD, request2, ierr)
        
        dest = mod(my_id-1 + width, width) + (my_id/width)*width !send left halo
        call mpi_isend(tend_chunk(1:2, 1:n, :), 2*n, MPI_REAL, dest, 3, MPI_COMM_WORLD, request3, ierr)
        
        dest = mod(my_id+1, width) + (my_id/width)*width !send right halo
        call mpi_isend(tend_chunk(m-1:m, 1:n, :), 2*n, MPI_REAL, dest, 4, MPI_COMM_WORLD, request4, ierr)
        
        !recv top halo
        source = mod(my_id, width) + mod((my_id/width + 1), height)*width
        call mpi_recv(qchunk4(3:m+2, n+3:n+4, :), 2*m, MPI_REAL, source, 1, MPI_COMM_WORLD, status, ierr)
        !recv bottom halo
        source = mod(my_id, width) + mod((my_id/width - 1) + height, height)*width
        call mpi_recv(qchunk4(3:m+2, 1:2, :), 2*m, MPI_REAL, source, 2, MPI_COMM_WORLD, status, ierr)
        !recv left halo
        source = mod(my_id-1 + width, width) + (my_id/width)*width
        call mpi_recv(qchunk4(1:2, 3:n+2, :), 2*n, MPI_REAL, source, 3, MPI_COMM_WORLD, status, ierr)
        !recv right halo
        source = mod(my_id+1, width) + (my_id/width)*width
        call mpi_recv(qchunk4(m+3:m+4, 3:n+2, :), 2*n, MPI_REAL, source, 4, MPI_COMM_WORLD, status, ierr)
        qchunk4(3:m+2, 3:n+2, :) = tend_chunk(:,:,:)
        qchunk4 = qchunk4*dt
        print *, 'last'
        call mpi_wait(request1, status, ierr)
        call mpi_wait(request2, status, ierr)
        call mpi_wait(request3, status, ierr)
        call mpi_wait(request4, status, ierr)
        print *, 'D'
		!s_chunk = s_chunk + (qchunk1 + 2.0*qchunk2 + 2.0*qchunk3 + qchunk4)/6.0
        !add and such to get next state(chunk)
                
        
    end subroutine calculate_RK4
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine set_initial_state(s)
    !these initial conditions are borrowed from Michael
        
      real, dimension(:), intent(out) :: s
      
      

      integer :: i, j
      real, dimension(matSize,matSize) :: h, u, v
      real :: d


      h(:,:) = 1.0
      u(:,:) = 0.0
      v(:,:) = 0.0

      do i=1,matSize
      do j=1,matSize
         d = sqrt((real(i)-real(matSize)/2.0)**2 + (real(j)-real(matSize)/2.0)**2) * 3.14159 / (real(matSize)/5.0)
         if (d <= 3.14159) then
            h(i,j) = h(i,j) + (cos(d) + 1.0)/20.0
         end if
      end do
      end do
      
      
     
      s(1:arrSize) = reshape(h,(/arrSize/))
      s(arrSize+1:2*arrSize) = reshape(u,(/arrSize/))
      s(2*arrSize+1:3*arrSize) = reshape(v,(/arrSize/))


    
        
                
    end subroutine set_initial_state
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine f_chunk(h_chunk, u_chunk, v_chunk, dh_chunk, du_chunk, dv_chunk)
    
        real, dimension(:,:), intent(in) :: h_chunk, u_chunk, v_chunk
        real, dimension(:,:), intent(out) :: dh_chunk, du_chunk, dv_chunk
        real :: divV, gradxh, gradyh
        integer :: m, n, i, im, ip, im2, ip2, j, jm, jp, jm2, jp2
        integer, dimension(2) :: stuff 
        
        stuff = shape(dh_chunk)
        m = stuff(1)
        n = stuff(2)  
        
        do i=3,m+2
            im = i-1 !since we do a 4-point central difference, need 4 values of i and j
            ip = i+1
            im2 = i-2
            ip2 = i+2
            
            do j=3,n+2
                jm = j-1
                jp = j+1
                jm2 = j-2
                jp2 = j+2
                
                divV = (-u_chunk(ip2,j)+8*u_chunk(ip,j)-8*u_chunk(im,j)+u_chunk(im2,j))/(12*dx) + (-v_chunk(i,jp2)+8*v_chunk(i,jp)-8*v_chunk(i,jm)+v_chunk(i,jm2))/(12*dy)
                gradxh = (-h_chunk(ip2,j)+8*h_chunk(ip,j)-8*h_chunk(im,j)+h_chunk(im2,j))/(12*dx) !these three quantities are 4-point central finite difference.
                gradyh = (-h_chunk(i,jp2)+8*h_chunk(i,jp)-8*h_chunk(i,jm)+h_chunk(i,jm2))/(12*dy)
                dh_chunk(i-2,j-2) = -h_chunk(i,j)*divV
                du_chunk(i-2,j-2) = -(divV*u_chunk(i,j) + 9.8*gradxh)
                dv_chunk(i-2,j-2) = -(divV*v_chunk(i,j) + 9.8*gradyh)
            end do
        end do
        
                
        return
    
    
    end subroutine f_chunk

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! DISASSEMBLE_STATE
    ! Turns state into <num_procs> different chunks, each chunk dimensioned by m+4, n+4
    ! where m, n are the optimal dimensions for breaking up the matrix into squarish pieces.
    ! m, n must be already known, and passed in by way of the size of statechunks.
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    subroutine disassemble_state(state, statechunks)
        implicit none
        real, dimension(:), intent(in) :: state
        real, dimension(:,:,:,:), intent(out) :: statechunks
        real, dimension(matsize, matsize) :: h, u, v
        integer :: i, j, p, q, l, m, n, x, num_procs, xdivs, ydivs
        integer, dimension(4) :: stuff
        stuff = shape(statechunks)
        m = stuff(1) - 4
        n = stuff(2) - 4
        num_procs = stuff(4)
        xdivs = matsize/m
        ydivs = matsize/n
        h = reshape(state(1:arrSize), (/matSize, matSize/))
                u = reshape(state(arrSize+1:2*arrSize), (/matSize, matSize/))
                v = reshape(state(2*arrSize+1:3*arrSize), (/matSize, matSize/))
        l = 0
        do j=1, ydivs
            do i=1, xdivs
                l = l + 1
            	
            	statechunks(3:m+2, 3:n+2, 1, l) = h((i-1)*m+1:i*m, (j-1)*n+1:j*n)
                statechunks(1:2, 3:n+2, 1, l) = h(mod((i-1)*m-1+xdivs*m, xdivs*m):mod((i-1)*m-1+xdivs*m, xdivs*m)+1, (j-1)*n+1:j*n)
                statechunks(m+3:m+4, 3:n+2, 1, l) = h(mod(i*m+1, xdivs*m):mod(i*m+2, xdivs*m), (j-1)*n+1:j*n)
                statechunks(3:m+2, 1:2, 1, l) = h((i-1)*m+1:i*m, mod((j-1)*n-1+ydivs*n, ydivs*n):mod((j-1)*n-1+ydivs*n, ydivs*n)+1)
                statechunks(3:m+2, n+3:n+4, 1, l) = h((i-1)*m+1:i*m, mod(j*n+1, ydivs*n):mod(j*n+2, ydivs*n))

                statechunks(3:m+2, 3:n+2, 2, l) = u((i-1)*m+1:i*m, (j-1)*n+1:j*n)
                statechunks(1:2, 3:n+2, 2, l) = u(mod((i-1)*m-1+xdivs*m, xdivs*m):mod((i-1)*m-1+xdivs*m, xdivs*m)+1, (j-1)*n+1:j*n)
                statechunks(m+3:m+4, 3:n+2, 2, l) = u(mod(i*m+1, xdivs*m):mod(i*m+2, xdivs*m), (j-1)*n+1:j*n)
                statechunks(3:m+2, 1:2, 2, l) = u((i-1)*m+1:i*m, mod((j-1)*n-1+ydivs*n, ydivs*n):mod((j-1)*n-1+ydivs*n, ydivs*n)+1)
                statechunks(3:m+2, n+3:n+4, 2, l) = u((i-1)*m+1:i*m, mod(j*n+1, ydivs*n):mod(j*n+2, ydivs*n))
                
                statechunks(3:m+2, 3:n+2, 3, l) = v((i-1)*m+1:i*m, (j-1)*n+1:j*n)
                statechunks(1:2, 3:n+2, 3, l) = v(mod((i-1)*m-1+xdivs*m, xdivs*m):mod((i-1)*m-1+xdivs*m, xdivs*m)+1, (j-1)*n+1:j*n)
                statechunks(m+3:m+4, 3:n+2, 3, l) = v(mod(i*m+1, xdivs*m):mod(i*m+2, xdivs*m), (j-1)*n+1:j*n)
                statechunks(3:m+2, 1:2, 3, l) = v((i-1)*m+1:i*m, mod((j-1)*n-1+ydivs*n, ydivs*n):mod((j-1)*n-1+ydivs*n, ydivs*n)+1)
                statechunks(3:m+2, n+3:n+4, 3, l) = v((i-1)*m+1:i*m, mod(j*n+1, ydivs*n):mod(j*n+2, ydivs*n))
            	
            
            end do
        end do    
                !this way we'll need like 27 calculations in total to make each chunk for h, u, v...
    
    end subroutine disassemble_state
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! RESASSEMBLE_TEND
    ! Constructs a new tendency vector out of the chunks in tendchunks, which are dimensioned
    ! by m, n. 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine reassemble_tend(tendchunks, tendency)
        implicit none
        real, dimension(:,:,:,:), intent(in) :: tendchunks
        real, dimension(:), intent(out) :: tendency
        real, dimension(matSize, matSize) :: h, u, v
        integer :: i, j, k, m, n, x, num_procs, xdivs, ydivs
        integer, dimension(4) :: stuff
        stuff = shape(tendchunks)
        m = stuff(1)
        n = stuff(2)
        num_procs = stuff(4)
        xdivs = matsize/m
        ydivs = matsize/n
        i = 0
        do j=1, ydivs
            do k=1, xdivs
                i = i+1
                h((k-1)*m+1:k*m, (j-1)*n+1:j*n) = tendchunks(:,:,1,i)
                u((k-1)*m+1:k*m, (j-1)*n+1:j*n) = tendchunks(:,:,2,i)
                v((k-1)*m+1:k*m, (j-1)*n+1:j*n) = tendchunks(:,:,3,i)
            end do
        end do
        tendency(1:arrSize) = reshape(h, (/arrSize/))
        tendency(arrSize+1:2*arrSize) = reshape(u, (/arrSize/))
        tendency(2*arrSize+1:3*arrSize) = reshape(v, (/arrSize/))
        !put tendchunks(:,:,1,i) into tendency
        !put tendchunks(:,:,2,i) into tendency
        !put tendchunks(:,:,3,i) into tendency
    
    end subroutine reassemble_tend
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
end module equations
