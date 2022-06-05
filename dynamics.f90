module dynamics

   use solver
   use gaussian
   use potential_module

   implicit none

contains

   subroutine qm_propa(time, param, pot)

      type(time_type), intent(inout) :: time
      type(gaussian_param(*,*,*)), intent(inout) :: param
      class(potential_type), intent(in) :: pot

      integer :: myunit
      ! open file to store data
      open(newunit = myunit, file = 'parameters.txt', form = 'formatted', &
         action = 'write', status = 'replace')

      write(myunit, *) time%itr * time%dt, param%q(1,1,1), param%p(1,1,1), &
         param%C(1,1,1,1)%re, param%C(1,1,1,1)%im, &
         param%a(1,1)%re, param%a(1,1)%im, &
         param%s(1,1)%re, param%s(1,1)%im

      do while ( time%itr < time%t / time%dt)
         ! later - detect crossing, detect decomposition, re-evolve,
         ! hop, re-evolve
         ! rename it to 
         if  (mod(time%itr, 40) == 0) then
            ! call a write to file function
            write(myunit, *) time%itr * time%dt, param%q(1,1,1), param%p(1,1,1), &
               param%C(1,1,1,1)%re, param%C(1,1,1,1)%im, &
               param%a(1,1)%re, param%a(1,1)%im, &
               param%s(1,1)%re, param%s(1,1)%im
         end if
         !rename it to solver_do_step()
         call do_step(param, time, pot)
         ! detect crossing 
         ! call sssh_hop(param, pot)?
         ! call hopper_hop()
         ! there will be a decomposition function inside sssh_hop


         time%itr = time%itr + 1

         ! something interesting to point out is that the
         ! specification of the solvers using this submodule approach
         ! will be done at compilation time rather than specifying a
         ! set of parameteres from a txt. file. Consider the
         ! comparison with the pontial OOP implementation
      end do

   end subroutine qm_propa

end module dynamics
