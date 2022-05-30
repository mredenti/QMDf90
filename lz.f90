module lz_module 

    use setup, only : mykind 

    implicit none 
    private 

    type, public :: lz_type

        real(mykind) :: delta
        real(mykind) :: alpha 

    contains 
        procedure :: gap => lz_gap
        procedure :: print => lz_print
    
    end type lz_type

contains 

    elemental function lz_gap(this, x) result(v_x)

        class(lz_type), intent(in) :: this 
        real(mykind), intent(in) :: x 

        real(mykind) :: v_x 

        v_x = sqrt(this%alpha**2*x**2 + this%delta**2)

    end function lz_gap

    subroutine lz_print(this)

        class(lz_type), intent(in) :: this 

        print *, 'Landau-Zener potential! delta = ', this%delta

    end subroutine lz_print

end module lz_module