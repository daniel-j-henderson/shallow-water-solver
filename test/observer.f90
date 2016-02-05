!observer

module observer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Module: observer
!
! This module is used to "observe" the state vector of a 
! system of equations by writing that state in some format to a file
! for later viewing or plotting by the user.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use netcdf
	use equations
	implicit none

    private
    
	character (len=*), parameter, public :: filename = 'shallow_water_data.nc'   ! Name of NetCDF file to be created
                 
    integer, private :: cnt=1      ! Return error code from NetCDF calls, track num of writes
    integer :: ncid       ! Handle to the NetCDF file
    integer :: xdimID, ydimID, tdimID    ! ID of the time and space dimensions
    integer :: hvarID, uvarID, vvarID    ! IDs of the spatial and velocity variables

    public :: observer_init, observer_write, observer_finalize, observer_write_chunk   


    contains


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Name: observer_init
    !
    ! Description: Initializes the observer module by, e.g., opening 
    !   files for later writing. This routine must be called before the 
    !   first call to observer_write().
    !
    ! Input: none
    !
    ! Output: none
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine observer_init()

        implicit none
        integer :: ierr
        ierr = nf90_create(filename, NF90_CLOBBER, ncid)
    	if (ierr /= NF90_NOERR) then
        	write(0,*) '*********************************************************************************'
        	write(0,*) 'Error creating NetCDF file '//filename
        	write(0,*) 'ierr = ', ierr
        	write(0,*) '*********************************************************************************'
        	stop
    	end if
    	
    	ierr = nf90_def_dim(ncid, 'timeDim', NF90_UNLIMITED, tdimID)
    	if (ierr /= NF90_NOERR) then
     	   write(0,*) '*********************************************************************************'
      	   write(0,*) 'Error defining dimension timeDim of unlimited length in file '//filename
       	   write(0,*) 'ierr = ', ierr
      	   write(0,*) '*********************************************************************************'
           stop
  		end if
  		
  		ierr = nf90_def_dim(ncid, 'xDim', matSize, xdimID)
    	if (ierr /= NF90_NOERR) then
     	   write(0,*) '*********************************************************************************'
      	   write(0,*) 'Error defining dimension xDim of unlimited length in file '//filename
       	   write(0,*) 'ierr = ', ierr
      	   write(0,*) '*********************************************************************************'
           stop
  		end if
  		
  		ierr = nf90_def_dim(ncid, 'yDim', matSize, ydimID)
    	if (ierr /= NF90_NOERR) then
     	   write(0,*) '*********************************************************************************'
      	   write(0,*) 'Error defining dimension timeDim of unlimited length in file '//filename
       	   write(0,*) 'ierr = ', ierr
      	   write(0,*) '*********************************************************************************'
           stop
  		end if
  		
  		
        ierr = nf90_def_var(ncid, 'height', NF90_REAL, (/xdimID, ydimID, tdimID/), hvarID)
	    if (ierr /= NF90_NOERR) then
 	       write(0,*) '*********************************************************************************'
 	       write(0,*) 'Error defining height variable in file '//filename
 	       write(0,*) 'ierr = ', ierr
   	       write(0,*) '*********************************************************************************'
           stop
        end if
        
        ierr = nf90_def_var(ncid, 'x component of velocity (u)', NF90_REAL, (/xdimID, ydimID, tdimID/), uvarID)
	    if (ierr /= NF90_NOERR) then
 	       write(0,*) '*********************************************************************************'
 	       write(0,*) 'Error defining uvar variable in file '//filename
 	       write(0,*) 'ierr = ', ierr
   	       write(0,*) '*********************************************************************************'
           stop
        end if
        
        ierr = nf90_def_var(ncid, 'y component of velocity', NF90_REAL, (/xdimID, ydimID, tdimID/), vvarID)
	    if (ierr /= NF90_NOERR) then
 	       write(0,*) '*********************************************************************************'
 	       write(0,*) 'Error defining vvar variable in file '//filename
 	       write(0,*) 'ierr = ', ierr
   	       write(0,*) '*********************************************************************************'
           stop
        end if
        
        ierr = nf90_enddef(ncid)
	    if (ierr /= NF90_NOERR) then
 	       write(0,*) '*********************************************************************************'
 	       write(0,*) 'Error ending def in file '//filename
 	       write(0,*) 'ierr = ', ierr
   	       write(0,*) '*********************************************************************************'
           stop
        end if

    end subroutine observer_init


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Name: observer_write
    !
    ! Description: Formats and writes the contents of the state vector s
    !   to a file.
    !
    ! Input: s -- the state vector
    !
    ! Output: none
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine observer_write(s)
    
    	!use netcdf
		
        implicit none
        real, dimension(:), intent(in) :: s
		real, dimension(matSize,matSize) :: h, u, v
        integer :: ierr
		h = reshape(s(1:arrSize), (/matSize, matSize/))
		u = reshape(s(arrSize+1:2*arrSize), (/matSize, matSize/))
		v = reshape(s(2*arrSize+1:3*arrSize), (/matSize, matSize/))

		
		ierr = nf90_put_var(ncid, hvarID, h, (/1, 1, cnt/), (/matSize, matSize, 1/))
	    if (ierr /= NF90_NOERR) then
 	       write(0,*) '*********************************************************************************'
 	       write(0,*) 'Error putting h in file '//filename
 	       write(0,*) 'ierr = ', ierr
   	       write(0,*) '*********************************************************************************'
           stop
        end if

        ierr = nf90_put_var(ncid, uvarID, u, (/1, 1, cnt/), (/matSize, matSize, 1/))
	    if (ierr /= NF90_NOERR) then
 	       write(0,*) '*********************************************************************************'
 	       write(0,*) 'Error putting u in file '//filename
 	       write(0,*) 'ierr = ', ierr
   	       write(0,*) '*********************************************************************************'
           stop
        end if
        
        ierr = nf90_put_var(ncid, vvarID, v, (/1, 1, cnt/), (/matSize, matSize, 1/))
	    if (ierr /= NF90_NOERR) then
 	       write(0,*) '*********************************************************************************'
 	       write(0,*) 'Error putting v in file '//filename
 	       write(0,*) 'ierr = ', ierr
   	       write(0,*) '*********************************************************************************'
           stop
        end if
        
        cnt=cnt+1 ! j keeps track of the next spot in each variable to write to

    end subroutine observer_write
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Name: observer_write_chunk
    !
    ! Description: Formats and writes the contents of the state chunk s
    !   to a file.
    !
    ! Input: s -- the state chunk
    !
    ! Output: none
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine observer_write_chunk(s, i, j)
    	
    	implicit none
        real, dimension(:,:,:), intent(in) :: s
        integer, intent(in) :: i, j
    	integer :: ierr, m, n
    	integer, dimension(3) :: stuff
    	
    	write (*,*) 'my i = ',i ,'and j = ',j
    	
    	! stuff = shape(s)
!     	m = stuff(1)
!     	n = stuff(2)
! 		
! 		ierr = nf90_put_var(ncid, hvarID, s(:,:,1), (/i*m+1, i*n+1, cnt/), (/m, n, 1/))
! 	    if (ierr /= NF90_NOERR) then
!  	       write(0,*) '*********************************************************************************'
!  	       write(0,*) 'Error putting h-chunk in file '//filename
!  	       write(0,*) 'ierr = ', ierr
!    	       write(0,*) '*********************************************************************************'
!            stop
!         end if
! 
!         ierr = nf90_put_var(ncid, uvarID, s(:,:,2), (/i*m+1, i*n+1, cnt/), (/m, n, 1/))
! 	    if (ierr /= NF90_NOERR) then
!  	       write(0,*) '*********************************************************************************'
!  	       write(0,*) 'Error putting u-chunk in file '//filename
!  	       write(0,*) 'ierr = ', ierr
!    	       write(0,*) '*********************************************************************************'
!            stop
!         end if
!         
!         ierr = nf90_put_var(ncid, vvarID, s(:,:,3), (/i*m+1, i*n+1, cnt/), (/m, n, 1/))
! 	    if (ierr /= NF90_NOERR) then
!  	       write(0,*) '*********************************************************************************'
!  	       write(0,*) 'Error putting v-chunk in file '//filename
!  	       write(0,*) 'ierr = ', ierr
!    	       write(0,*) '*********************************************************************************'
!            stop
!         end if
!         
!         cnt=cnt+1
    
    end subroutine observer_write_chunk


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Name: observer_finalize
    !
    ! Description: Finalizes the observer module by, e.g., closing any
    !   files that were opened by the module. This routine must be called 
    !   only once after all calls to observer_write() have been made.
    !
    ! Input: none
    !
    ! Output: none
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine observer_finalize()

        implicit none
        integer :: ierr
    ierr = nf90_close(ncid)
    if (ierr /= NF90_NOERR) then
        write(0,*) '*********************************************************************************'
        write(0,*) 'Error while closing NetCDF file '//filename
        write(0,*) 'ierr = ', ierr
        write(0,*) '*********************************************************************************'
        stop
    end if

    end subroutine observer_finalize

end module observer