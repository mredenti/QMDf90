module constants 

    use setup, only : mykind
    
    implicit none 
    public

    real(mykind),  parameter :: PI = 4.D0*DATAN(1.D0)
    ! add chemistry constants

end module constants