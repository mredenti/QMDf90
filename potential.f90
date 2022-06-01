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
        ! 1st order derivatives 
        procedure(interface_v11d), deferred :: v11d
        procedure(interface_v12d), deferred :: v12d
        procedure(interface_v21d), deferred :: v21d
        procedure(interface_v22d), deferred :: v22d
        ! 2nd order derivatives 
        procedure(interface_v11dd), deferred :: v11dd
        procedure(interface_v12dd), deferred :: v12dd
        procedure(interface_v21dd), deferred :: v21dd
        procedure(interface_v22dd), deferred :: v22dd
        ! trace of matrix 
        procedure :: trace
        ! eigenvalues 
        procedure :: rho
        procedure :: trace_d
        procedure :: rho_d
        procedure :: trace_dd
        procedure :: rho_dd
        procedure :: gap
        procedure :: potup
        procedure :: potdown
        procedure :: potup_d
        procedure :: potdown_d
        procedure :: potup_dd
        procedure :: potdown_dd
        
    
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

        subroutine interface_v11dd(this, x, hessian)
                
            use setup, only : mykind
            import potential_type
            
            class(potential_type), intent(in) :: this
            real(mykind), intent(in) :: x 
            ! i am not sure about this dimension
            real(mykind), intent(out), dimension(:,:) :: hessian 

        end subroutine interface_v11dd

        subroutine interface_v12dd(this, x, hessian)
            
            use setup, only : mykind
            import potential_type
            
            class(potential_type), intent(in) :: this
            real(mykind), intent(in) :: x 
            ! i am not sure about this dimension
            real(mykind), intent(out), dimension(:,:) :: hessian 
    
    
        end subroutine interface_v12dd

        subroutine interface_v21dd(this, x, hessian)
            
            use setup, only : mykind
            import potential_type
            
            class(potential_type), intent(in) :: this
            real(mykind), intent(in) :: x 
            ! i am not sure about this dimension
            real(mykind), intent(out), dimension(:,:) :: hessian
    
        end subroutine interface_v21dd

        subroutine interface_v22dd(this, x, hessian)
            
            use setup, only : mykind
            import potential_type
            
            class(potential_type), intent(in) :: this
            real(mykind), intent(in) :: x 
            ! i am not sure about this dimension
            real(mykind), intent(out), dimension(:,:) :: hessian 
    
        end subroutine interface_v22dd

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

        u = sqrt( ((this%v11(x) - this%v22(x))/2.0)**2 + this%v12(x)**2 )
        
    end function rho

    subroutine trace_d(this, x, grad)   

        class(potential_type), intent(in) :: this
        real(mykind), intent(in) :: x 
        real(mykind), dimension(:), intent(out) :: grad
            
        real(mykind), dimension(this%dim) :: v11_grad, v22_grad
        
        call this%v11d(x, v11_grad)
        call this%v22d(x, v22_grad)

        grad = ( v11_grad + v22_grad ) / 2.0
        
    end subroutine trace_d

    subroutine rho_d(this, x, grad)   

        class(potential_type), intent(in) :: this
        real(mykind), intent(in) :: x 
        real(mykind), dimension(:), intent(out) :: grad
        
        real(mykind), dimension(this%dim) :: v11_grad, v22_grad, v12_grad
        
        call this%v11d(x, v11_grad)
        call this%v11d(x, v22_grad)
        call this%v12d(x, v12_grad)
  
        grad = ( (this%v11(x) - this%v22(x)) / 2.0 * (v11_grad - v22_grad) / 2.0 &
                    + this%v12(x) * v12_grad ) / this%rho(x)
        
    end subroutine rho_d

    subroutine trace_dd(this, x, hessian)   

        class(potential_type), intent(in) :: this
        real(mykind), intent(in) :: x 
        real(mykind), dimension(:,:), intent(out) :: hessian
            
        real(mykind), dimension(this%dim, this%dim) :: v11_hessian, v22_hessian
        
        call this%v11dd(x, v11_hessian)
        call this%v22dd(x, v22_hessian)

        hessian = ( v11_hessian + v22_hessian ) / 2.0
        
    end subroutine trace_dd

    subroutine rho_dd(this, x, hessian)   

        class(potential_type), intent(in) :: this
        real(mykind), intent(in) :: x 
        real(mykind), dimension(:,:), intent(out) :: hessian
       
        real(mykind) :: Z
        real(mykind), dimension(this%dim) :: v11_grad, v22_grad, v12_grad, Z_grad, C_grad
        real(mykind), dimension(this%dim, this%dim) :: v11_hessian, v22_hessian, v12_hessian, Z_hessian
        
        call this%v11d(x, v11_grad)
        call this%v11d(x, v22_grad)
        call this%v12d(x, v12_grad)

        call this%v11dd(x, v11_hessian)
        call this%v11dd(x, v22_hessian)
        call this%v22dd(x, v12_hessian)

        Z = ( this%v11(x) - this%v22(x) ) / 2.0
        Z_grad = ( v11_grad - v22_grad ) / 2.0
        Z_hessian = (v11_hessian - v12_hessian) / 2.0
        C_grad = Z * Z_grad + this%v12(x) * v12_grad
        
        hessian = ( matmul( reshape(Z_grad , (/ this%dim, 1 /)), reshape(Z_grad, (/ 1, this%dim /)) ) &
                    + Z *  Z_hessian &
                    + matmul( reshape(v12_grad , (/ this%dim, 1 /)), reshape(v12_grad, (/ 1, this%dim /))) &
                    + this%v12(x) *  v12_hessian) &
                    / this%rho(x) &
                    - matmul( reshape(C_grad , (/ this%dim, 1 /)), reshape(C_grad, (/ 1, this%dim /))) & 
                    / this%rho(x)**3
        
    end subroutine rho_dd

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

    subroutine potup_d(this, x, grad)   

        class(potential_type), intent(in) :: this
        real(mykind), intent(in) :: x 
        real(mykind), dimension(:), intent(out) :: grad
        
        real(mykind), dimension(this%dim) :: trace_grad
        real(mykind), dimension(this%dim) :: rho_grad

        call this%trace_d(x, trace_grad)
        call this%rho_d(x, rho_grad)

        grad = trace_grad + rho_grad
        
    end subroutine potup_d

    subroutine potdown_d(this, x, grad)   

        class(potential_type), intent(in) :: this
        real(mykind), intent(in) :: x 
        real(mykind), dimension(:), intent(out) :: grad
            
        real(mykind), dimension(this%dim) :: trace_grad
        real(mykind), dimension(this%dim) :: rho_grad

        call this%trace_d(x, trace_grad)
        call this%rho_d(x, rho_grad)

        grad = trace_grad - rho_grad
        
    end subroutine potdown_d

    subroutine potup_dd(this, x, hessian)   

        class(potential_type), intent(in) :: this
        real(mykind), intent(in) :: x 
        real(mykind), dimension(:,:), intent(out) :: hessian
            
        real(mykind), dimension(this%dim, this%dim) :: trace_hessian
        real(mykind), dimension(this%dim, this%dim) :: rho_hessian

        call this%trace_dd(x, trace_hessian)
        call this%rho_dd(x, rho_hessian)

        hessian = trace_hessian + rho_hessian
        
    end subroutine potup_dd

    subroutine potdown_dd(this, x, hessian)   

        class(potential_type), intent(in) :: this
        real(mykind), intent(in) :: x 
        real(mykind), dimension(:,:), intent(out) :: hessian
            
        real(mykind), dimension(this%dim, this%dim) :: trace_hessian
        real(mykind), dimension(this%dim, this%dim) :: rho_hessian

        call this%trace_dd(x, trace_hessian)
        call this%rho_dd(x, rho_hessian)

        hessian = trace_hessian - rho_hessian
        
    end subroutine potdown_dd

    

end module potential_module