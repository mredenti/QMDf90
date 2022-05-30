submodule(solver) variational_gaussian 

    implicit none  


contains 

    module subroutine do_step(param, time, potential)
        ! implementation of do_step declared in solver
        type (gaussian_param(*,*,*)), intent(inout) :: param
        type (time_type), intent(in) :: time
        class (potential_type), intent(in) :: potential
        ! Numerical integrator based on variational splitting
        ! norm-preserving, symplectic, time-reversible 
        ! Lubich's blue book : Algorithm 2.1, 
        ! Section 4.4: VARIATIONAL SPLITTING FOR GAUSSIAN WAVE PACKETS

        !integer :: i 
        
        ! potential half-step 
        ! compute average gradient and average hessian
        param%p = param%p - time%dt / 2.0
        param%C = param%C
        param%s = param%s - time%dt / 2.0

        ! kinetic full-step 
        param%q = param%q + time%dt * param%p
        param%C = param%C
        param%s = param%s - time%dt

        ! potential half-step
        ! compute average gradient and average hessian
        param%p = param%p - time%dt / 2.0
        param%C = param%C
        param%s = param%s - time%dt / 2.0
        

    end subroutine do_step

end submodule variational_gaussian 