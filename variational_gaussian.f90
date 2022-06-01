submodule(solver) variational_gaussian 

    use hermite

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

        integer :: i
        integer, parameter :: N = 50

        ! average potential 
        !real(mykind), dimension(potential%dim) :: pot_avg
        real(mykind) :: V_avg
        ! GH weights
        real (mykind),      dimension(N)    :: w ! GH weights
        ! GH nodes
        real (mykind),      dimension(N)    :: y ! GH weights
        real (mykind),      dimension(N)    :: x ! GH weights
        ! gradient of potential 
        !real(mykind), dimension(potential%dim) :: potup_grad
        real(mykind), dimension(potential%dim) :: V_grad, V_grad_avg
        real(mykind), dimension(potential%dim, potential%dim) :: V_hessian, V_hessian_avg 

        call cgqf (N, y, w )

        V_avg = 0
        V_grad_avg = 0
        V_hessian_avg = 0
        do i = 1, N 
            x(i) = y(i) * param%eps / param%C(1,1,1,1)%im + param%q(1,1,1)
            call potential%potup_d(x(i), V_grad)
            call potential%potup_dd(x(i), V_hessian)
            V_avg = V_avg + w(i)*potential%potup(x(i))
            V_grad_avg = V_grad_avg + w(i)*V_grad
            V_hessian_avg = V_hessian_avg + w(i)*V_hessian
        end do 
        
        ! potential half-step 

        ! compute average gradient and average hessian
        param%p(:,1,1) = param%p(:,1,1) - time%dt / 2.0 * V_grad_avg
        param%s = param%s - time%dt / 2.0 * V_avg &
                + time%dt * param%eps / 8.0 / param%C(1,1,1,1)%im * V_hessian_avg(1,1)
        param%C(:,:,1,1) = param%C(:,:,1,1) - time%dt / 2.0 * V_hessian_avg

        ! kinetic full-step 
        param%q = param%q + time%dt * param%p
        param%s = param%s + time%dt / 2.0 * param%p(1,1,1)**2 & 
                + cmplx(0.0, 1.0) / 2.0 * param%eps * log(1 + time%dt * param%C(1,1,1,1))
        param%C = param%C / (1 + time%dt * param%C)
     

        ! potential half-step
        ! compute average gradient and average hessian
        V_avg = 0
        V_grad_avg = 0
        V_hessian_avg = 0
        do i = 1, N 
            x(i) = y(i) * param%eps / param%C(1,1,1,1)%im + param%q(1,1,1)
            call potential%potup_d(x(i), V_grad)
            call potential%potup_dd(x(i), V_hessian)
            V_avg = V_avg + w(i)*potential%potup(x(i))
            V_grad_avg = V_grad_avg + w(i)*V_grad
            V_hessian_avg = V_hessian_avg + w(i)*V_hessian
        end do 
        param%p(:,1,1) = param%p(:,1,1) - time%dt / 2.0 * V_grad_avg
        param%s = param%s - time%dt / 2.0 * V_avg &
                + time%dt * param%eps / 8.0 / param%C(1,1,1,1)%im * V_hessian_avg(1,1)
        param%C(:,:,1,1) = param%C(:,:,1,1) - time%dt / 2.0 * V_hessian_avg


    end subroutine do_step

end submodule variational_gaussian 