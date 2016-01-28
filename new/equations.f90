!equations


module equations
	
	
	private
	integer :: i, j, im, ip, jm, jp 
	real :: dx = 100 ! dx and dy are required for the numerical differentiation I used to get the divergences/gradients
	real :: dy= 100 ! for arbitrary dx, step_size should be about 1/10-1/5 of dx to yield results that don't blow up. 
	
	public :: get_system_size, set_initial_state, f
	integer, parameter, public :: matSize = 100 !set the dimensions of the grid of fluid. Interestingly, it seg faults for 400+
	integer, public :: arrSize
	
	
	contains
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	integer function get_system_size()
	
		arrSize = matSize**2
		get_system_size = 3*arrSize
		return
		
	end function get_system_size
	
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
		integer :: m, n
		
		(/m, n/) = shape(dh_chunk)
				
		
		
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
				
				divV = (-u_chunk(ip2,j)+8*u_chunk(ip,j)-8*u_chunk(im,j)+u_chunk(im2,j))/(12*dx) + (-v_chunk(i,jp2)+8*v_chunk(i,jp)-8*v(i,jm)+v_chunk(i,jm2))/(12*dy)
				gradxh = (-h_chunk(ip2,j)+8*h_chunk(ip,j)-8*h_chunk(im,j)+h_chunk(im2,j))/(12*dx) !these three quantities are 4-point central finite difference.
				gradyh = (-h_chunk(i,jp2)+8*h_chunk(i,jp)-8*h_chunk(i,jm)+h_chunk(i,jm2))/(12*dy)
				dh(i-2,j-2) = -h_chunk(i,j)*divV
				du(i-2,j-2) = -(divV*u_chunk(i,j) + 9.8*gradxh)
				dv(i-2,j-2) = -(divV*v_chunk(i,j) + 9.8*gradyh)
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
		real, dimension(:,:,:,:), intent(inout) :: statechunks
		real, dimension(matsize, matsize) :: h, u, v
		integer :: i, j, l, m, n, x, num_procs, xdivs, ydivs
		(/m, n, x, num_procs/) = shape(statechunks)
		xdivs = matsize/m
		ydivs = matzsize/n
		
		l = 0
		do j=1, ydivs
			do i=1, xdivs
				l = l + 1
				statechunks((3:m+2, 3:n+2, 1, l) = h((i-1)*m+1:i*m, (j-1)*n+1:j*n)
				statechunks(1:2, 1:n+4, 1, l) = h((
				
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
		real, dimension(:), intent(inout) :: tendency
		integer :: i, m, n, x, num_procs
		(/m, n, x, num_procs/) = shape(tendchunks)
		do i=1, num_procs
			!put tendchunks(:,:,1,i) into tendency
			!put tendchunks(:,:,2,i) into tendency
			!put tendchunks(:,:,3,i) into tendency
		end do
	
	end subroutine reassemble_tend
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
end module equations
