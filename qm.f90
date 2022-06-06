module qm

   use gaussian

   implicit none

contains

   subroutine qm_init(param)

      type(gaussian_param(*,*,*)), intent(inout) :: param

   end subroutine

end module qm
