program main

    use potential_lz_module
    use gaussian
    use solver
    use setup, only : mykind 

    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    implicit none

    type(gaussian_param(1,1,1)) :: param0 
    type(potential_lz_type) :: pot 
    real(mykind) :: alpha
    real(mykind) :: delta 
    type(time_type) :: time 
    integer :: myunit

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

    ! set up LZ potential parameters 
    ! case(potential_name) set up corresponding potential
    pot = potential_lz_type(1, alpha, delta)

    ! open file to store data 
    open(newunit = myunit, file = 'parameters.txt', form = 'formatted', &
      action = 'write', status = 'replace')

    ! set up time 
    time = time_type(5.0, 0.01)

    ! ######################################################################

    ! QM_RUN() ! perhaps passing some parameters to the run subroutine?
    do while ( time%itr < time%t / time%dt) 
        ! later - detect crossing, detect decomposition, re-evolve, 
        ! hop, re-evolve
        ! rename it to solver_do_step()
        call do_step(param0, time, pot) 
        
        time%itr = time%itr + 1

        if  (mod(time%itr, 10) == 0) then
            ! call a write to file function 
            write(myunit, *) param0%q(1,1,1), param0%p(1,1,1), &
                             param0%C(1,1,1,1)%re, param0%C(1,1,1,1)%im, & 
                             param0%a(1,1)%re, param0%a(1,1)%im
        end if  
        
        ! something interesting to point out is that the 
        ! specification of the solvers using this submodule approach 
        ! will be done at compilation time rather than specifying a 
        ! set of parameteres from a txt. file. Consider the 
        ! comparison with the pontial OOP implementation 
    end do 

end program main
