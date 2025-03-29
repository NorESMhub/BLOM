module mod_output_forcing

   use netcdf,    only: nf90_64bit_offset, nf90_global, nf90_noerr, nf90_nofill, nf90_def_dim,   &
                        nf90_enddef, nf90_close, nf90_create, nf90_strerror, &
                        nf90_double, nf90_def_var, nf90_put_var, nf90_put_att
   use mod_types, only: r8
   use mod_xc

   implicit none
   private

   public :: output_forcing

contains

   subroutine output_forcing(filename, varname, kpie, kpje, vardata)

      ! Arguments
      character(len=*) , intent(in) :: filename
      character(len=*) , intent(in) :: varname
      integer          , intent(in) :: kpie, kpje
      real(r8)         , intent(in) :: vardata(kpie,kpje)

      ! Local variables
      integer :: ncid,ncvarid,ncstat,ncdims(2),nclatid,nclonid
      integer :: i,j
      integer :: start(2),count(2)
      real    :: arr_g(itdm,jtdm)
      real    :: arr_l(1-nbdy:idm+nbdy, 1-nbdy:jdm+nbdy, 1)

      if (mnproc==1) then
         write(lp,'(a)') 'creating netcdf file '//trim(filename)

         ! open file
         ncstat = nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid)
         if ( ncstat  /=  NF90_NOERR ) then
            call xchalt('(cannot open output forcing data file)'); stop
         endif

         ! defined dimensions
         ncstat = nf90_def_dim(ncid, 'lon', itdm, nclonid)
         if ( ncstat  /=  NF90_NOERR ) then
            call xchalt('(cannot create lon dim)'); stop
         endif
         ncstat = nf90_def_dim(ncid, 'lat', jtdm, nclatid)
         if ( ncstat  /=  NF90_NOERR ) then
            call xchalt('(cannot create lat dim)'); stop
         endif
         ncdims(1) = nclonid
         ncdims(2) = nclatid

         ! define variable
         ncstat = NF90_DEF_VAR(ncid, trim(varname), NF90_DOUBLE, ncdims, ncvarid)
         if ( ncstat /=  NF90_NOERR ) then
            call xchalt('(cannot define variable)'); stop
         endif

         ! end definition
         ncstat = nf90_enddef(ncid)
      end if

      arr_l(:,:,:) = 0.0
      start(1) = 1; count(1) = itdm
      start(2) = 1; count(2) = jtdm

      do j=1,jdm
         do i=1,idm
            arr_l(i,j,1) = vardata(i,j)
         enddo
      enddo
      call xcaget(arr_g, arr_l, 1)

      if (mnproc == 1) then
         ! Output variable data
         ncstat = nf90_put_var(ncid, ncvarid, arr_g, start, count)
         if (ncstat /= nf90_noerr) then
            write(lp,'(4a)') 'nf90_put_var: ',trim(varname),': ',nf90_strerror(ncstat)
            call xchalt('(write_netcdf_var)'); stop
         endif

         ! Close file
         ncstat = nf90_close(ncid)
         if ( ncstat  /=  NF90_NOERR ) then
            call xchalt('(cannot close file)'); stop
         endif
      end if

    end subroutine output_forcing

 end module mod_output_forcing
