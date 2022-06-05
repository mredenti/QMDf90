submodule(solver) variational_gaussian

   use hermite
   use constants, only : PI

   implicit none


contains

   subroutine do_step(param, time, potential)
      ! implementation of do_step declared in solver
      type (gaussian_param(*,*,*)), intent(inout) :: param
      type (time_type), intent(in) :: time
      class (potential_type), intent(in) :: potential
      ! Numerical integrator based on variational splitting
      ! norm-preserving, symplectic, time-reversible
      ! Lubich's blue book : Algorithm 2.1,
      ! Section 4.4: VARIATIONAL SPLITTING FOR GAUSSIAN WAVE PACKETS

      call potential_step(param, time%dt/2.0, potential)
      call kinetic_step(param, time%dt)
      call potential_step(param, time%dt/2.0, potential)

   end subroutine do_step

   subroutine kinetic_step(param, dt)
      ! implementation of do_step declared in solver
      type (gaussian_param(*,*,*)), intent(inout) :: param
      real(mykind), intent(in) :: dt

      ! kinetic full-step
      param%q = param%q + dt * param%p
      param%s = param%s + dt / 2.0 * param%p(1,1,1)**2 &
         + cmplx(0.0, 1.0) / 2.0 * param%eps * log(1 + dt * param%C(1,1,1,1))
      param%C = param%C / (1 + dt * param%C)

   end subroutine kinetic_step

   subroutine potential_step(param, dt, potential)
      ! PROBLEM: you are generating the samples twice
      ! for no reason

      ! implementation of do_step declared in solver
      type (gaussian_param(*,*,*)), intent(inout) :: param
      real(mykind), intent(in) :: dt
      class (potential_type), intent(in) :: potential

      real(mykind) :: V_avg
      real(mykind), dimension(potential%dim) :: V_grad, V_grad_avg
      real(mykind), dimension(potential%dim, potential%dim) :: V_hessian, V_hessian_avg


      integer :: i
      integer, parameter :: N = 50 ! THIS SHOULD BE AN ADDITIONAL PARAMETER - HOW TO CHOOSE IT?

      ! average potential
      !real(mykind), dimension(potential%dim) :: pot_avg
      ! GH weights
      real (mykind),      dimension(N)    :: w ! GH weights
      ! GH nodes
      real (mykind),      dimension(N)    :: y ! GH weights
      real (mykind),      dimension(N)    :: x ! GH weights
      ! gradient of potential
      !real(mykind), dimension(potential%dim) :: potup_grad

      call cgqf (N, y, w )

      V_avg = 0.0
      V_grad_avg = 0.0
      V_hessian_avg = 0.0
      do i = 1, N
         x(i) = y(i) * param%eps / param%C(1,1,1,1)%im + param%q(1,1,1)
         !if (param%state)
         call potential%potup_d(x(i), V_grad)
         call potential%potup_dd(x(i), V_hessian)
         V_avg = V_avg + w(i)*potential%potup(x(i))
         V_grad_avg = V_grad_avg + w(i)*V_grad
         V_hessian_avg = V_hessian_avg + w(i)*V_hessian
      end do

      V_avg = V_avg * (PI * param%eps)**(- 0.5) * param%C(1,1,1,1)%im**(0.5)
      V_grad_avg = V_grad_avg * (PI * param%eps)**(- 0.5) * param%C(1,1,1,1)%im**(0.5)
      V_hessian_avg = V_hessian_avg * (PI * param%eps)**(- 0.5) * param%C(1,1,1,1)%im**(0.5)


      ! potential half-step
      !print *, V_grad_avg
      ! compute average gradient and average hessian
      param%p = param%p - dt / 2.0 * V_grad_avg(1)
      param%s = param%s - dt / 2.0 * V_avg &
         + dt * param%eps / 8.0 / param%C(1,1,1,1)%im * V_hessian_avg(1,1)
      param%C(:,:,1,1) = param%C(:,:,1,1) - dt / 2.0 * V_hessian_avg

   end subroutine potential_step

end submodule variational_gaussian
