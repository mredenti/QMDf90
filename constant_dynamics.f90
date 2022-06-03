program main

   use potential_constant_module
   use gaussian
   use dynamics
   use setup, only : mykind

   use, intrinsic :: iso_fortran_env, only: stderr => error_unit
   implicit none

   type(gaussian_param(1,1,1)) :: param0
   type(potential_constant_type) :: pot
   real(mykind) :: alpha, delta, c
   type(time_type) :: time

   ! QM_INIT()
   ! CALL Init_Model(QModel,pot_name='Tully',Print_init=.FALSE.)
   ! initial condition parameters
   param0%eps = 0.1
   param0%q = - 10
   param0%p = 4
   param0%C = (0.0, 1.0)
   param0%a = (1.0, 0.0)
   param0%s = (0.0, 0.0)

   ! potential parameters
   alpha = 0.5
   delta = 0.5
   c = 0.5

   ! set up LZ potential parameters
   ! case(potential_name) set up corresponding potential
   pot = potential_constant_type(1, alpha, delta, c)

   ! set up time
   time = time_type(5.0, 0.01)

   ! ######################################################################
   ! qm_read(input_txt)
   ! qm_init()
   call qm_propa(time, param0, pot)
   ! qm_finalise()

end program main
