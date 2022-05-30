! one approach to implement the actual potentials would be to use 
! submodules
module potential_module

    use setup, only : mykind

    implicit none 
    private 

    type, public, abstract :: potential_type 
        
        integer :: dim 

    contains 
        ! diabatic potential matrix 
        procedure(interface_v11), deferred :: v11
        procedure(interface_v12), deferred :: v12
        procedure(interface_v21), deferred :: v21
        procedure(interface_v22), deferred :: v22
        ! derivatives 
        procedure(interface_v11d), deferred :: v11d
        procedure(interface_v12d), deferred :: v12d
        procedure(interface_v21d), deferred :: v21d
        procedure(interface_v22d), deferred :: v22d
        procedure :: trace
        procedure :: rho
        procedure :: gap
        procedure :: potup
        procedure :: potdown
    
    !contains
        !procedure :: rho => mod_potential_rho
        !procedure :: vd => mod_potential_vd
    
    end type potential_type

    interface 
        elemental function interface_v11(this, x) result(v)
            
            use setup, only : mykind
            import potential_type
            
            class(potential_type), intent(in) :: this
            real(mykind), intent(in) :: x 
            
            ! return value
            real(mykind) :: v 

        end function interface_v11
        elemental function interface_v12(this, x) result(v)
            
            use setup, only : mykind
            import potential_type
            
            class(potential_type), intent(in) :: this
            real(mykind), intent(in) :: x 
            
            ! return value
            real(mykind) :: v
        end function interface_v12
        elemental function interface_v21(this, x) result(v)
            
            use setup, only : mykind
            import potential_type
            
            class(potential_type), intent(in) :: this
            real(mykind), intent(in) :: x 
            
            ! return value
            real(mykind) :: v

        end function interface_v21
        elemental function interface_v22(this, x) result(v)
            
            use setup, only : mykind
            import potential_type
            
            class(potential_type), intent(in) :: this
            real(mykind), intent(in) :: x 
            
            ! return value
            real(mykind) :: v

        end function interface_v22
    
        subroutine interface_v11d(this, x, grad)
                
            use setup, only : mykind
            import potential_type
            
            class(potential_type), intent(in) :: this
            real(mykind), intent(in) :: x 
            ! i am not sure about this dimension
            real(mykind), intent(out), dimension(:) :: grad 

        end subroutine interface_v11d

        subroutine interface_v12d(this, x, grad)
            
            use setup, only : mykind
            import potential_type
            
            class(potential_type), intent(in) :: this
            real(mykind), intent(in) :: x 
            ! i am not sure about this dimension
            real(mykind), intent(out), dimension(:) :: grad 
    
    
        end subroutine interface_v12d

        subroutine interface_v21d(this, x, grad)
            
            use setup, only : mykind
            import potential_type
            
            class(potential_type), intent(in) :: this
            real(mykind), intent(in) :: x 
            ! i am not sure about this dimension
            real(mykind), intent(out), dimension(:) :: grad 
    
        end subroutine interface_v21d

        subroutine interface_v22d(this, x, grad)
            
            use setup, only : mykind
            import potential_type
            
            class(potential_type), intent(in) :: this
            real(mykind), intent(in) :: x 
            ! i am not sure about this dimension
            real(mykind), intent(out), dimension(:) :: grad 
    
        end subroutine interface_v22d

    end interface

contains 

    elemental function trace(this, x) result(u)   

        class(potential_type), intent(in) :: this
        real(mykind), intent(in) :: x 
            
        ! return value
        real(mykind) :: u

        u = ( this%v11(x) + this%v22(x) ) / 2.0 
        
    end function trace

    elemental function rho(this, x) result(u)   

        class(potential_type), intent(in) :: this
        real(mykind), intent(in) :: x 
            
        ! return value
        real(mykind) :: u

        u = sqrt( ((this%v11(x) - this%v22(x))/2.0)**2 - this%v12(x)**2 )
        
    end function rho

    elemental function gap(this, x) result(u)

        class(potential_type), intent(in) :: this
        real(mykind), intent(in) :: x 
            
        ! return value
        real(mykind) :: u

        u = this%potup(x) - this%potdown(x)
        
    end function gap 

    elemental function potup(this, x) result(u)

        class(potential_type), intent(in) :: this
        real(mykind), intent(in) :: x 
            
        ! return value
        real(mykind) :: u

        u = this%trace(x) + this%rho(x)
        
    end function potup

    elemental function potdown(this, x) result(u)

        class(potential_type), intent(in) :: this
        real(mykind), intent(in) :: x 
            
        ! return value
        real(mykind) :: u

        u = this%trace(x) - this%rho(x)
        
    end function potdown

end module potential_module