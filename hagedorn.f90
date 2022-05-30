module hagedorn 

    use setup, only : mykind 
    use gaussian  
    ! the type might be different for the Hagedorn wavepacket ?
    use potential_module

    implicit none 
    private 

    ! change name - define it somewhere else
    type, public :: time 

        real(mykind) :: dt 
        real(mykind) :: t 
    
    end type time 

    public :: do_step

contains 

    subroutine do_step(param, time_type, potential)

        type (gaussian_param(*,*,*)), intent(inout) :: param
        type (time), intent(in) :: time_type 
        class (potential_type), intent(in) :: potential

        integer :: i 

        print *, potential%v11(2.0_mykind)
        print *, potential%dim
        !do i = 1, param%nq 
            !param%q(:, i) = 2 * time_type%dt
        !end do 

    end subroutine do_step

end module hagedorn 