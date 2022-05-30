!> Interface to Gaussian processing library.
!>
!> ...
module gaussian

   use hermite
   use setup, only : mykind
   use constants, only : PI

   implicit none
   public

   ! parameter is like const in C
   !integer, parameter :: mykind = kind(1.d0)
   integer, parameter :: mydim = 1

   ! define a type for the Gaussian basis

   type, public :: gaussian_param(dim, nq, np)

      integer, len :: dim
      integer, len :: nq
      integer, len :: np
      real(mykind)                           :: eps
      complex(mykind), dimension(nq,np)     :: a
      real(mykind), dimension(dim, nq, np)   :: q
      real(mykind), dimension(dim, nq, np)   :: p
      complex(mykind), dimension(dim, dim, 1) :: C

   end type gaussian_param

   ! An example of the time-sliced thawed Gaussian approximation

   ! The motivation behind using Gaussian wave functions or Hagedorn
   ! wavepackets is that the time propagations consists in solving
   ! a system of ODEs and it is exact for quadratic potentials

   ! Find out the easiest one to implement among the following three:

   ! 1) MCTDH
   ! 2) Variational Gaussian wave-packets
   ! 3) Time-sliced thawed Gaussian propagation https://arxiv.org/pdf/2108.12182.pdf (why would you re-initialise the evolved basis in the time-independent approximation space?)
   ! 4) Variational multi-configurational Gaussian wave packets

contains

   function gaussian_evaluate(x, q, p, C, eps) result(y) ! return value of Gaussian on a grid

      real (mykind),      dimension(:),   intent(in) :: x ! grid
      real (mykind),      dimension(:),   intent(in) :: q ! position
      real (mykind),      dimension(:),   intent(in) :: p ! momentum
      real (mykind),                      intent(in) :: eps ! momentum
      complex (mykind),   dimension(:,:), intent(in) :: C ! covariance matrix
      complex (mykind),                   parameter  :: j = (0.0, 1.0)
      complex (mykind),   dimension(size(x))         :: y

      y = (PI*eps)**(- 0.25) * (C(1,1)%im)**(0.25) * exp(j / eps * (0.5*(x - q(1))**2*C(1,1) + p(1)*(x-q(1))))

   end function gaussian_evaluate
   
   !> Decomposes a Gaussian as a linear combination of Gaussians
   subroutine gaussian_decompose_uniform(param0, param1)

      ! the big question here is how to choose the variance

      ! https://arxiv.org/pdf/2010.03478.pdf

      ! param0 is the Gaussian wavepacket to be decomposed
      type (gaussian_param(*,*,*)),      intent(in) :: param0
      ! param1 is the collection of Gaussian wavepackets
      type (gaussian_param(*,*,*)),      intent(inout) :: param1

      complex (mykind),      dimension(param1%nq,param1%np)  :: b ! coefficients

      ! fix a square box around q (might want to consider rectangular boxes in higher dimension)
      real(mykind) :: L_q
      ! fix a square box around p
      real(mykind) :: L_p
      ! spacing
      real(mykind) :: dq, dp

      integer :: j, k

      L_q = 8.0 ! it should depend on the variance of psi_tc - call a subroutine
      L_p = 4.0*PI
      dq = 2.0*L_q / param1%nq
      dp = 2.0*L_p / param1%np
      param1%C = (0.0, 1.0) ! gamma

      do j = 1, param1%np
         ! midpoint rule
         param1%p(1, : , j) = param0%p(1, 1, 1) - L_p + dp / 2.0 * (2.0*j - 1.0)
      end do

      do k = 1, param1%nq
         ! midpoint rule
         param1%q(1, k, :) = param0%q(1,1,1) - L_q + dq / 2.0 * (2.0*k - 1.0)
      end do

      call gaussian_inner_product(param0, param1, b)

      do j = 1, param1%nq ! should be careful how to loop - this is wrong
         do k = 1, param1%np
            ! this product occurs also in GH quadrature rules
            param1%a(j,k) = dp / (2.0*PI*param0%eps) * dq * b(j,k)
         end do
      end do


   end subroutine gaussian_decompose_uniform

   subroutine gaussian_decompose_gauss_hermite(param0, param1)

      ! the big question here is how to choose the variance

      ! https://arxiv.org/pdf/2010.03478.pdf

      ! param0 is the Gaussian wavepacket to be decomposed
      type (gaussian_param(*,*,*)),      intent(in) :: param0
      ! param1 is the collection of Gaussian wavepackets
      type (gaussian_param(*,*,*)),      intent(inout) :: param1
      complex (mykind),      dimension(param1%nq,param1%np)  :: b ! coefficients
      
      real (mykind),      dimension(param1%np)    :: w ! GH weights

      ! fix a square box around q (might want to consider rectangular boxes in higher dimension)
      real(mykind) :: L_q ! it should depend on the variance of psi_tc
      ! spacing
      real(mykind) :: dq

      integer :: j, k

      L_q = 8
      dq = 2.0*L_q / param1%nq

      param1%C = (0.0, 1.0) ! gamma
     
      do k = 1, param1%nq
         ! midpoint rule
         param1%q(1,k) = param0%q(1,1) - L_q + dq / 2.0 * (2.0*k - 1.0)
      end do

      ! call to obtain N nodes and N weights
      call cgqf ( param1%np, param1%p(1,:), w )

      ! transform nodes and weights
      do j = 1, param1%np
         ! midpoint rule
         w(j) = exp(param1%p(1,j)**2)*w(j)*sqrt(2.0 * param0%eps) ! prior to updating pj
         param1%p(1,j) = param0%p(1,1) + param1%p(1,j)*sqrt(2.0 * param0%eps)  ! sampled from ?
      end do

      call gaussian_inner_product(param0, param1, b)

      do k = 1, param1%np ! should be careful how to loop - this is wrong
         do j = 1, param1%nq
            param1%a(j,k) = dq / (2.0*PI*param0%eps) * b(j,k) * w(k)
         end do
      end do

   end subroutine gaussian_decompose_gauss_hermite

   subroutine gaussian_inner_product(param0, param1, b)

      ! subroutine to compute the inner product between two gaussians
      ! param0 is the Gaussian wavepacket to be decomposed
      type (gaussian_param(*,*,*)),      intent(in) :: param0
      ! param1 is the collection of Gaussian wavepackets
      type (gaussian_param(*,*,*)),      intent(inout) :: param1

      complex (mykind),      dimension(param1%nq,param1%np), intent(out):: b ! coefficients
      ! nodes (assume d=1)
      complex(mykind) :: A, D1, D2
      ! coefficients
      complex(mykind) :: i = (0.0, 1.0)
      integer :: j, k

      A = i / (param0%C(1,1,1) - conjg(param1%C(1,1,1)))

      D1 = 1.0 / (1.0 / param0%C(1,1,1) - 1.0 / (conjg(param1%C(1,1,1))))
      D2 = - 1.0 / ( param0%C(1,1,1) - (conjg(param1%C(1,1,1))) )


      do j = 1, param1%nq ! should be careful how to loop - this is wrong
         do k = 1, param1%np
            ! this product occurs also in GH quadrature rules
            b(j,k) = 2**(0.5)*(param1%C(1,1,1)%im * param0%C(1,1,1)%im)**(0.25) / sqrt(1.0 / A) &
               * exp( i / 2.0 / param0%eps * (param1%p(1,k) + param0%p(1,1)) * (param1%q(1,j) - param0%q(1,1)) ) &
               * exp( 1.0 / 2.0 / param0%eps * (param1%p(1,k) - param0%p(1,1)) * A * (param0%C(1,1,1) + conjg(param1%C(1,1,1))) &
               * (param1%q(1,j) - param0%q(1,1)) ) &
               * exp(i / 2.0 / param0%eps * ( (param1%p(1,k) - param0%p(1,1))**2 * D2 + (param1%q(1,j) - param0%q(1,1))**2 * D1 ))
         end do
      end do

   end subroutine gaussian_inner_product

   function gaussian_variance() result(var)

      real(mykind) :: var

      var = 2

   end function gaussian_variance


   subroutine print_kind()

      print *, "mykind = ", mykind

   end subroutine print_kind

end module gaussian


