program test_potential

    use setup, only : mykind 
    use potential_lz_module 
    use variational_gaussian 

    implicit none 

    type(potential_lz_type) :: pot

    real(mykind) :: x = 0.5
    real(mykind) :: v_x

    type(time) :: t = time(0.5,5)

    pot = potential_lz_type(1, 0.5, 0.5)
    ! can pass in the values directly

    v_x = pot%v11(x)
    print *, v_x
    v_x = pot%gap(x)
    print *, v_x

    call do_step(t, pot)


end program test_potential
