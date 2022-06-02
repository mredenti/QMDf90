module potential_constant_module

   use potential_module
   use setup, only : mykind

   implicit none
   public

   !> @brief derived type in which the LZ parameters are set-up.
   type, extends (potential_type) :: potential_constant_type

      real(mykind) :: delta
      real(mykind) :: alpha 
      real(mykind) :: c

   contains
      ! diabatic potential matrix
      procedure :: v11
      procedure :: v12
      procedure :: v21
      procedure :: v22
      !
      procedure :: v11d
      procedure :: v12d
      procedure :: v21d
      procedure :: v22d
      !
      procedure :: v11dd
      procedure :: v12dd
      procedure :: v21dd
      procedure :: v22dd
      ! theta
      procedure :: theta
      procedure :: thetad
      procedure :: thetadd

   end type potential_constant_type

contains

   elemental function theta(this, x) result(v)

      class(potential_constant_type), intent(in) :: this
      real(mykind), intent(in) :: x

      ! return value
      real(mykind) :: v

      v = this%c / this%alpha * atan(tanh(this%alpha / 2.0 * x))

   end function theta

   subroutine thetad(this, x, grad)

      class(potential_constant_type), intent(in) :: this
      real(mykind), intent(in) :: x
      ! i am not sure about this dimension
      real(mykind), intent(out), dimension(:) :: grad

      grad = this%alpha

   end subroutine thetad

   subroutine thetadd(this, x, grad)

      class(potential_constant_type), intent(in) :: this
      real(mykind), intent(in) :: x
      ! i am not sure about this dimension
      real(mykind), intent(out), dimension(:) :: grad

      grad = this%alpha

   end subroutine thetadd

   elemental function v11(this, x) result(v)

      class(potential_constant_type), intent(in) :: this
      real(mykind), intent(in) :: x

      ! return value
      real(mykind) :: v

      v = this%delta

   end function v11

   elemental function v12(this, x) result(v)

      class(potential_constant_type), intent(in) :: this
      real(mykind), intent(in) :: x

      ! return value
      real(mykind) :: v

      v = 0.0

   end function v12

   elemental function v21(this, x) result(v)

      class(potential_constant_type), intent(in) :: this
      real(mykind), intent(in) :: x

      ! return value
      real(mykind) :: v

      v = 0.0

   end function v21

   elemental function v22(this, x) result(v)

      class(potential_constant_type), intent(in) :: this
      real(mykind), intent(in) :: x

      ! return value
      real(mykind) :: v

      v = - this%delta 

   end function v22

   subroutine v11d(this, x, grad)

      class(potential_constant_type), intent(in) :: this
      real(mykind), intent(in) :: x
      ! i am not sure about this dimension
      real(mykind), intent(out), dimension(:) :: grad

      grad = 0.0

   end subroutine v11d

   subroutine v12d(this, x, grad)

      class(potential_constant_type), intent(in) :: this
      real(mykind), intent(in) :: x
      ! i am not sure about this dimension
      real(mykind), intent(out), dimension(:) :: grad

      grad = 0.0

   end subroutine v12d

   subroutine v21d(this, x, grad)

      class(potential_constant_type), intent(in) :: this
      real(mykind), intent(in) :: x
      ! i am not sure about this dimension
      real(mykind), intent(out), dimension(:) :: grad

      grad = 0.0

   end subroutine v21d

   subroutine v22d(this, x, grad)

      class(potential_constant_type), intent(in) :: this
      real(mykind), intent(in) :: x
      ! i am not sure about this dimension
      real(mykind), intent(out), dimension(:) :: grad

      grad = 0.0

   end subroutine v22d


   subroutine v11dd(this, x, hessian)

      class(potential_constant_type), intent(in) :: this
      real(mykind), intent(in) :: x
      ! i am not sure about this dimension
      real(mykind), intent(out), dimension(:,:) :: hessian

      hessian = 0.0

   end subroutine v11dd

   subroutine v12dd(this, x, hessian)

      class(potential_constant_type), intent(in) :: this
      real(mykind), intent(in) :: x
      ! i am not sure about this dimension
      real(mykind), intent(out), dimension(:,:) :: hessian

      hessian = 0.0

   end subroutine v12dd

   subroutine v21dd(this, x, hessian)

      class(potential_constant_type), intent(in) :: this
      real(mykind), intent(in) :: x
      ! i am not sure about this dimension
      real(mykind), intent(out), dimension(:,:) :: hessian

      hessian = 0.0

   end subroutine v21dd

   subroutine v22dd(this, x, hessian)

      class(potential_constant_type), intent(in) :: this
      real(mykind), intent(in) :: x
      ! i am not sure about this dimension
      real(mykind), intent(out), dimension(:,:) :: hessian

      hessian = 0.0

   end subroutine v22dd

end module potential_constant_module
