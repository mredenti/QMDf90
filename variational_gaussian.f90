module variational_gaussian 

    use setup, only : mykind 
    !use gaussian  
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

    subroutine do_step(time_type, potential)
        
        ! Numerical integrator based on variational splitting
        ! norm-preserving, symplectic, time-reversible 
        ! Lubich's blue book : Algorithm 2.1, 
        ! Section 4.4: VARIATIONAL SPLITTING FOR GAUSSIAN WAVE PACKETS

        !type (gaussian_param(*,*,*)), intent(inout) :: param
        type (time), intent(in) :: time_type 
        class (potential_type), intent(in) :: potential

        integer :: i 

        print *, potential%v11(2.0_mykind)
        print *, potential%dim
        !do i = 1, param%nq 
            !param%q(:, i) = 2 * time_type%dt
        !end do 

    end subroutine do_step

end module variational_gaussian 