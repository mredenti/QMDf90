module sssh

   use gaussian

   implicit none

contains

   function sssh_hop(param) result(hop)
      ! detects crossing, changes state and applies SA formula
      type(gaussian_param(*,*,*)), intent(in) :: param

      logical :: hop

      if ( (param%gap_old - param%gap_old_old) * (param%gap - param%gap_old)  < 0) then
         hop = .true.
      else
         hop = .false.
      end if

   end function sssh_hop

end module sssh
