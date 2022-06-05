module save

   use netcdf
   use solver
   use gaussian

   implicit none

   character (len = *), parameter  :: TIME_NAME = "time"
   character (len = *), parameter  :: ITR_NAME = "itr"
   character (len = *), parameter  :: DIM_NAME = "dim"
   character (len = *), parameter :: POS_NAME = "position"
   character (len = *), parameter :: MOM_NAME = "momentum"
   character (len = *), parameter :: COVARIANCE_IM_NAME = "covariance.im"
   character (len = *), parameter :: COVARIANCE_RE_NAME = "covariance.re"
   character (len = *), parameter :: A_RE_NAME = "a.re"
   character (len = *), parameter :: A_IM_NAME = "a.im"
   character (len = *), parameter :: EPS_NAME = "eps"
   character (len = *), parameter :: UNITS = "units"
   character (len = *), parameter :: POS_UNITS = "a.u."
   character (len = *), parameter :: AU_UNITS = "a.u."
   character (len = *), parameter :: TIME_STEP = "dt"
   character (len = *), parameter :: TIME_T = "t"
   integer :: REC_COUNT = 1 

contains

   subroutine check(stat)
      integer, intent(in) :: stat

      if (stat /= NF90_NOERR) then
         print '(a)', trim(nf90_strerror(stat))
         stop
      end if
   end subroutine check

   subroutine save_open_gaussian(file_name, param, time)

      character(len=*), intent(in)    :: file_name
      type(gaussian_param(*,*,*)), intent(inout) :: param
      type(time_type), intent(in) :: time

      ! ncid is the file id
      integer                         :: ncid, pos_varid, mom_varid, time_varid, itr_varid
      integer                         :: cov_im_varid, cov_re_varid, are_varid, aim_varid


      !
      integer                         :: dim_dimid, time_dimid

      ! We will create two netCDF variables, one each for temperature and
      ! pressure fields.
      integer                        :: pos_dimid, mom_dimid
      
      ! Create the NetCDF file. Override file, if it already exists.
      call check(nf90_create(file_name, NF90_CLOBBER, ncid))

      ! The record dimension is defined to have
      ! unlimited length - it can grow as needed. 
      call check(nf90_def_dim(ncid, DIM_NAME, param%dim, dim_dimid))
      call check(nf90_def_dim(ncid, TIME_NAME, NF90_UNLIMITED, time_dimid))
      !call check(nf90_def_dim(ncid, 'p', size(array, 2), p_dimid))
      
      
      call check( nf90_def_var(ncid, ITR_NAME, NF90_INT, time_dimid, itr_varid) )
      call check( nf90_def_var(ncid, TIME_NAME, NF90_DOUBLE, time_dimid, time_varid) )

      ! Define the variable type (NF90_INT: 4-byte integer).
      call check(nf90_def_var(ncid, POS_NAME, NF90_DOUBLE, (/ dim_dimid, time_dimid /), pos_varid))
      call check(nf90_def_var(ncid, MOM_NAME, NF90_DOUBLE, (/ dim_dimid, time_dimid /), mom_varid))
      call check(nf90_def_var(ncid, COVARIANCE_RE_NAME, NF90_DOUBLE, (/ dim_dimid, dim_dimid, time_dimid /), cov_re_varid))
      call check(nf90_def_var(ncid, COVARIANCE_IM_NAME, NF90_DOUBLE, (/ dim_dimid, dim_dimid, time_dimid /), cov_im_varid))
      call check(nf90_def_var(ncid, A_RE_NAME, NF90_DOUBLE, (/ time_dimid /), are_varid))
      call check(nf90_def_var(ncid, A_IM_NAME, NF90_DOUBLE, (/ time_dimid /), aim_varid))


      ! ! Assign units attributes to the netCDF variables.
      call check( nf90_put_att(ncid, pos_varid, UNITS, AU_UNITS) )

      ! time attributes
      call check( nf90_put_att(ncid, time_varid, UNITS, AU_UNITS) )
      call check( nf90_put_att(ncid, time_varid, TIME_STEP, time%dt) )
      call check( nf90_put_att(ncid, time_varid, TIME_T, time%t) )

      ! global attributes (includes potential information, eps, ...)
      call check (nf90_put_att(ncid, NF90_GLOBAL, DIM_NAME, param%dim) )
      call check (nf90_put_att(ncid, NF90_GLOBAL, EPS_NAME, param%eps) )

      ! End define mode.
      call check(nf90_enddef(ncid))

   end subroutine save_open_gaussian

   subroutine save_write_gaussian(file_name, param, time)
      
      character(len=*), intent(in)    :: file_name
      type(gaussian_param(*,*,*)), intent(inout) :: param
      type(time_type), intent(in) :: time

      integer           :: ncid, pos_varid, mom_varid, time_varid, itr_varid
      integer           :: covre_varid, covim_varid, are_varid, aim_varid

      ! Open the file. 
      call check( nf90_open(file_name, nf90_write, ncid) )

      ! Get the `varid` of the data variable, based on its name.
      call check( nf90_inq_varid(ncid, ITR_NAME, itr_varid) )
      call check( nf90_inq_varid(ncid, TIME_NAME, time_varid) )
      call check( nf90_inq_varid(ncid, POS_NAME, pos_varid) )
      call check( nf90_inq_varid(ncid, MOM_NAME, mom_varid) )
      call check( nf90_inq_varid(ncid, COVARIANCE_RE_NAME, covre_varid) )
      call check( nf90_inq_varid(ncid, COVARIANCE_IM_NAME, covim_varid) )
      call check( nf90_inq_varid(ncid, A_RE_NAME, are_varid) )
      call check( nf90_inq_varid(ncid, A_IM_NAME, aim_varid) )


      ! Write the data to the file. (do n = 1, max itr)
      ! time 
      call check(nf90_put_var(ncid, itr_varid, time%itr, start=[REC_COUNT] ))
      call check(nf90_put_var(ncid, time_varid, time%itr * time%dt, start=[REC_COUNT] ))
      ! gaussian parameters
      call check(nf90_put_var(ncid, pos_varid, param%q(:,1,1), start=[param%dim, REC_COUNT] ))
      call check(nf90_put_var(ncid, mom_varid, param%p(:,1,1), start=[param%dim, REC_COUNT] ))
      call check(nf90_put_var(ncid, covre_varid, param%C(:,:,1,1)%re, start=[param%dim, param%dim , REC_COUNT] ))
      call check(nf90_put_var(ncid, covim_varid, param%C(:,:,1,1)%im, start=[param%dim, param%dim , REC_COUNT] ))
      call check(nf90_put_var(ncid, are_varid, param%a(1,1)%re, start=[REC_COUNT] ))
      call check(nf90_put_var(ncid, aim_varid, param%a(1,1)%im, start=[REC_COUNT] ))
      

      REC_COUNT = REC_COUNT + 1
      ! Close the file.
      call check(nf90_close(ncid))

   end subroutine save_write_gaussian

end module save
