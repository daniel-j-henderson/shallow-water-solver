function f(t,s)
	
		real, intent(in) :: t
		real, dimension(:), intent(in) :: s	
		
		real, dimension(size(s)) :: f
		
		real, dimension(matSize, matSize) :: h,u,v, dh, du, dv
				
		h = reshape(s(1:10000), (/matSize, matSize/))
		u = reshape(s(10001:20000), (/matSize, matSize/))
		v = reshape(s(20001:30000), (/matSize, matSize/))
		
		do i=2,matSize-1
			!im = mod(i-1+100-1,100)+1
			!ip = mod(i+100,100)+1
			im = i-1
			ip = i+1
			do j=2,matSize-1
				!jm = mod(j-1+100-1,100)+1
				!jp = mod(j+100,100)+1
				jm = j-1
				jp = j+1
				divV = (u(ip, j) - u(im, j)) / (2*dx) + (v(i,jp) - v(i,jm))/(2*dy)
				gradxh = (h(ip,j)-h(im,j))/(2*dx)
				gradyh = (h(i,jp)-h(i,jm))/(2*dy)
				dh(i,j) = -h(i,j)*divV
				du(i,j) = -(divV*u(i,j) + 9.8*gradxh)
				dv(i,j) = -(divV*v(i,j) + 9.8*gradyh)
			end do
		end do
				
				
		f(1:10000) = reshape(dh, (/10000/))
		f(10001:20000) = reshape(du, (/10000/))
		f(20001:30000) = reshape(dv, (/10000/))
				
		return
	
	
	end function f
	
	!!
		
	  real, dimension(:), intent(out) :: s

      integer :: i, j
      real, dimension(100,100) :: h, u, v
      real :: d

      h(:,:) = 1.0
      u(:,:) = 1.0
      v(:,:) = 0.0

      do i=1,100
      do j=1,100
         d = sqrt((real(i)-50.0)**2 + (real(j)-50.0)**2) * 3.14159 / 20.0
         if (d <= 3.14159) then
            h(i,j) = h(i,j) + (cos(d) + 1.0)/20.0
         end if
      end do
      end do

      s(1:10000) = reshape(h,(/10000/))
      s(10001:20000) = reshape(u,(/10000/))
      s(20001:30000) = reshape(v,(/10000/))

	
	!!
	
	real, dimension(:), intent(out) :: s
		real, dimension(matSize, matSize) :: h,u,v
		
		h(:, :) = 1.0
		u(:, :) = 1.0
		v(:, :) = 0.0
		h(50, 50) = 2.0

		
	
		
		s(1:arrSize) = reshape(h, (/arrSize/))
		s(arrSize+1:2*arrSize) = reshape(u, (/arrSize/))
		s(2*arrSize+1:3*arrSize) = reshape(v, (/arrSize/))
	
	
	
	
	
	
	
	
	
	
	function f(t,s)
	
		real, intent(in) :: t
		real, dimension(:), intent(in) :: s	
		
		real, dimension(size(s)) :: f
		
		real, dimension(matSize, matSize) :: h,u,v, dh, du, dv
				
		h = reshape(s(1:10000), (/matSize, matSize/))
		u = reshape(s(10001:20000), (/matSize, matSize/))
		v = reshape(s(20001:30000), (/matSize, matSize/))
		
		do i=2,matSize-1
			!im = mod(i-1+100-1,100)+1
			!ip = mod(i+100,100)+1
			im = i-1
			ip = i+1
			do j=2,matSize-1
				!jm = mod(j-1+100-1,100)+1
				!jp = mod(j+100,100)+1
				jm = j-1
				jp = j+1
				divV = (u(ip, j) - u(im, j)) / (2*dx) + (v(i,jp) - v(i,jm))/(2*dy)
				gradxh = (h(ip,j)-h(im,j))/(2*dx)
				gradyh = (h(i,jp)-h(i,jm))/(2*dy)
				dh(i,j) = -h(i,j)*divV
				du(i,j) = -(divV*u(i,j) + 9.8*gradxh)
				dv(i,j) = -(divV*v(i,j) + 9.8*gradyh)
			end do
		end do
				
				
		f(1:10000) = reshape(dh, (/10000/))
		f(10001:20000) = reshape(du, (/10000/))
		f(20001:30000) = reshape(dv, (/10000/))
				
		return
	
	
	end function f