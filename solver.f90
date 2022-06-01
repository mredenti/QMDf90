module solver 

    use setup, only : mykind 
    use gaussian 
    use potential_module
    
    implicit none 
    
    type, public :: time_type 

        real(mykind) :: t 
        real(mykind) :: dt 
        integer :: itr = 1

    end type time_type   

    interface
        ! defer the implementation to the particular solver
        module subroutine do_step(param, time, potential)

            type (gaussian_param(*,*,*)), intent(inout) :: param
            type (time_type), intent(in) :: time
            class (potential_type), intent(in) :: potential

        end subroutine do_step
    end interface

end module solver