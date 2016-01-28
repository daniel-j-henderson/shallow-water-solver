module solve

   use equations
	
   contains
	

   subroutine advance_state(s, t, dt)

      real, dimension(:), intent(inout) :: s
      real, intent(inout) :: t
      real, intent(in) :: dt

      real, dimension(size(s)) :: q1, q2, q3, q4

! RK4
      q1(:) = dt * f(t,s(:)            )
      q2(:) = dt * f(t,s(:) + 0.5*q1(:))
      q3(:) = dt * f(t,s(:) + 0.5*q2(:))
      q4(:) = dt * f(t,s(:) +     q3(:))

      s(:) = s(:) + (q1(:) + 2.0*q2(:) + 2.0*q3(:) + q4(:))/6.0

! Forward Euler
!      s(:) = dt * f(t,s(:)) + s(:)

      t = t + dt

   end subroutine advance_state

end module solve
