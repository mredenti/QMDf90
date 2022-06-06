module dynamics

   use save 
   use solver
   use gaussian
   use sssh
   use potential_module

   implicit none

contains

   subroutine qm_propa(file_name, time, param, pot)

      character(len=*), intent(in)    :: file_name
      type(time_type), intent(inout) :: time
      type(gaussian_param(*,*,*)), intent(inout) :: param
      class(potential_type), intent(in) :: pot
      
      ! define variables
      call save_open_gaussian(file_name, param, time)

      do while ( time%itr < time%t / time%dt)
         
         if  (mod(time%itr, 40) == 0) then
            call save_write_gaussian(file_name, param, time)
         end if
         !rename it to solver_do_step()
         call do_step(param, time, pot)
         ! update gap information 
         param%gap_old_old = param%gap_old
         param%gap_old = param%gap
         param%gap = pot%gap(param%q(1,1,1))
         
         print *, param%gap, param%gap_old
         

         ! later - detect crossing, detect decomposition, re-evolve,
         ! hop, re-evolve
         ! rename it to 
         ! update gap information ?
         ! detect crossing 
         if (sssh_hop(param)) then 
            param%state = .not. param%state
            ! call superadiabatic_transition()
         end if 
         ! call sssh_hop(param, pot)?
         ! call hopper_hop('sssh')
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
