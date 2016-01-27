!equations


module equations
	
	
	private
	integer :: i, j, im, ip, jm, jp 
	real :: dx = 100 ! dx and dy are required for the numerical differentiation I used to get the divergences/gradients
	real :: dy= 100 ! for arbitrary dx, step_size should be about 1/10-1/5 of dx to yield results that don't blow up. 
	real :: divV, gradxh, gradyh
	
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
	
	subroutine set_initial_state(s, setup)
	!these initial conditions are borrowed from Michael
		
	  real, dimension(:), intent(out) :: s
	  integer, intent(in) :: setup
	  
	  

      integer :: i, j
      real, dimension(matSize,matSize) :: h, u, v
      real :: d

	  if (setup == 1) then

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
      
      end if
      
      if (setup == 2) then
      
      h(:,:) = 1.0
      u(:,:) = 1.0
      v(:,:) = 0.0
      
      do i=matsize/2-9, matsize/2+9
      	do j=matsize/2 - modulo(i, matsize/2-10), matsize/2 + modulo(i, matsize/2-9)
      		u(i,j) = 0.0
      	end do
      end do     	
      
      end if
      
      if (setup == 3) then
        
      h(:,:) = 1.0
      u(:,:) = 1.0
      v(:,:) = 0.0

      u(40:60, 40:60) = 0.0

      end if 
     
      s(1:arrSize) = reshape(h,(/arrSize/))
      s(arrSize+1:2*arrSize) = reshape(u,(/arrSize/))
      s(2*arrSize+1:3*arrSize) = reshape(v,(/arrSize/))


	
		
				
	end subroutine set_initial_state
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	function f(t,s,setup)
	
		real, intent(in) :: t
		real, dimension(:), intent(in) :: s	
		integer, intent(in) :: setup
		
		real, dimension(size(s)) :: f
		
		real, dimension(matSize, matSize) :: h,u,v, dh, du, dv
				
		h = reshape(s(1:arrSize), (/matSize, matSize/))
		u = reshape(s(arrSize+1:2*arrSize), (/matSize, matSize/))
		v = reshape(s(2*arrSize+1:3*arrSize), (/matSize, matSize/))
		
		do i=1,matSize
			im = mod(i-1+matSize-1,matSize)+1 !since we do a 4-point central difference, need 4 values of i and j
			ip = mod(i,matSize)+1
			im2 = mod(i-2+matSize-1,matSize)+1
			ip2 = mod(i+1,matSize)+1
			
			do j=1,matSize
				jm = mod(j-1+matSize-1,matSize)+1
				jp = mod(j,matSize)+1
				jm2 = mod(j-2+matSize-1,matSize)+1
				jp2 = mod(j+1,matSize)+1
				
				divV = (-u(ip2,j)+8*u(ip,j)-8*u(im,j)+u(im2,j))/(12*dx) + (-v(i,jp2)+8*v(i,jp)-8*v(i,jm)+v(i,jm2))/(12*dy)
				gradxh = (-h(ip2,j)+8*h(ip,j)-8*h(im,j)+h(im2,j))/(12*dx) !these three quantities are 4-point central finite difference.
				gradyh = (-h(i,jp2)+8*h(i,jp)-8*h(i,jm)+h(i,jm2))/(12*dy)
				dh(i,j) = -h(i,j)*divV
				du(i,j) = -(divV*u(i,j) + 9.8*gradxh)
				dv(i,j) = -(divV*v(i,j) + 9.8*gradyh)
			end do
		end do
		
		if (setup == 2) then
		do i=matsize/2-9, matsize/2+9
      		do j=matsize/2 - modulo(i, matsize/2-10), matsize/2 + modulo(i, matsize/2-9)
      			du(i,j) = 0.0
      			!dh(i,j) = 0.0
      		end do
        end do 
       end if

        if (setup == 3) then
        du(40:60, 40:60) = 0.0
 
        end if
				
		f(1:arrSize) = reshape(dh, (/arrSize/))
		f(arrSize+1:2*arrSize) = reshape(du, (/arrSize/))
		f(2*arrSize+1:3*arrSize) = reshape(dv, (/arrSize/))
				
		return
	
	
	end function f

	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
end module equations
