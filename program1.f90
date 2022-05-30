program program1

   use gaussian
   use setup, only : mykind
   
   implicit none


   integer, parameter :: M = 16 ! grid points in position space
   integer, parameter :: N = 20 ! grid points in momentum space
   integer :: myunit
   integer :: j,k
   !real (mykind),  dimension(N) :: x
   !real (mykind),  dimension(M) :: q
   !real (mykind),  dimension(N) :: p
   !complex (mykind),   dimension(1,1) :: C
   !complex (mykind),  dimension(M,N) :: r

   type (gaussian_param(:,:,:)), allocatable :: param0 
   type (gaussian_param(:,:,:)), allocatable :: param1

   ! allocate memory at run time - init_subroutine?
   allocate(gaussian_param(dim=1,nq=1,np=1) :: param0)
   allocate(gaussian_param(dim=1,nq=M,np=N) :: param1)
   
   ! set wavepacket parameters - d =1      
   param0%eps = 1.0
   param0%a = 1
   param0%q = 0.0
   param0%p = 0.0
   param0%C = cmplx(0.0, 1.0, mykind)

   ! project onto Gaussians using the truncate midpoint riemann summation
   call gaussian_decompose_uniform(param0, param1)

   ! project onto Gaussians using GH rule
   !call gaussian_decompose_gauss_hermite(param0, param1)

   !print *, y

   open(newunit = myunit, file = 'parameters.txt', form = 'formatted', &
      action = 'write', status = 'replace')

   ! save routine - to save data to a filename
   do j=1,M
      do k=1,N
         write(myunit, *) param1%q(1,j), param1%p(1,k), param1%C(1,1,1)%re, param1%C(1,1,1)%im, param1%a(j,k)%re, param1%a(j,k)%im
      end do
   end do

   !close(unit = myunit, status = 'keep')

   ! For the dynamics I would have the following routines 
   ! qm_setup()
      ! initial_condition from a file where you specify parameters
   ! qm_run()
      ! qm_propa()
      ! qm_hop()
      ! qm_observables
   !qm_finalise()
   deallocate(param0)
   deallocate(param1)
end program program1
