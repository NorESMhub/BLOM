! Copyright (C) 2020  I. Bethke
!
! This file is part of BLOM/iHAMOCC.
!
! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free 
! Software Foundation, either version 3 of the License, or (at your option) 
! any later version. 
!
! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details. 
!
! You should have received a copy of the GNU Lesser General Public License 
! along with BLOM. If not, see https://www.gnu.org/licenses/.


module ncutils

   use netcdf
   use types

   implicit none

   interface ncgetvar
      module procedure &
         ncgetvar_i2_scalar, ncgetvar_i2_1d, ncgetvar_i2_2d, ncgetvar_i2_3d, &
         ncgetvar_i2_4d, &
         ncgetvar_i4_scalar, ncgetvar_i4_1d, ncgetvar_i4_2d, ncgetvar_i4_3d, &
         ncgetvar_i4_4d, &
         ncgetvar_r4_scalar, ncgetvar_r4_1d, ncgetvar_r4_2d, ncgetvar_r4_3d, &
         ncgetvar_r4_4d, &
         ncgetvar_r8_scalar, ncgetvar_r8_1d, ncgetvar_r8_2d, ncgetvar_r8_3d, &
         ncgetvar_r8_4d
   end interface ncgetvar

contains

   subroutine ncopen(filename, ncid)

      character (len = *), intent(in) :: filename
      integer, intent(out) :: ncid

      integer status

      status = nf90_open(filename, nf90_nowrite, ncid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncopen: nf90_open: ', trim(filename), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
   
   end subroutine ncopen

   subroutine ncclose(ncid)

      integer, intent(in) :: ncid

      integer status

      status = nf90_close(ncid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncclose: nf90_close: ', trim(nf90_strerror(status))
         stop
      endif
   
   end subroutine ncclose

   subroutine ncgetdimlen(ncid, dimname, dimlen)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: dimname
      integer, intent(out) :: dimlen

      integer status, dimid

      status = nf90_inq_dimid(ncid, dimname, dimid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetdimlen: nf90_inq_dimid: ', trim(dimname), &
            ': ', trim(nf90_strerror(status))
         stop
      endif
      status = nf90_inquire_dimension(ncid, dimid, len = dimlen)
      if (status /= nf90_noerr) then
         write (*,'(2a)') ' nf90_inquire_dimensions: x: ', &
            trim(nf90_strerror(status))
         write (*,'(4a)') 'ncgetdimlen: nf90_inquire_dimension: ', &
            trim(dimname), ': ', trim(nf90_strerror(status))
         stop
      endif
      
   end subroutine ncgetdimlen

   subroutine ncgetvar_i2_scalar(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      integer (i2), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_i2_scalar

   subroutine ncgetvar_i2_1d(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      integer (i2), dimension(:), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_i2_1d

   subroutine ncgetvar_i2_2d(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      integer (i2), dimension(:,:), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_i2_2d

   subroutine ncgetvar_i2_3d(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      integer (i2), dimension(:,:,:), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_i2_3d

   subroutine ncgetvar_i2_4d(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      integer (i2), dimension(:,:,:,:), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_i2_4d

   subroutine ncgetvar_i4_scalar(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      integer (i4), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_i4_scalar

   subroutine ncgetvar_i4_1d(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      integer (i4), dimension(:), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_i4_1d

   subroutine ncgetvar_i4_2d(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      integer (i4), dimension(:,:), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_i4_2d

   subroutine ncgetvar_i4_3d(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      integer (i4), dimension(:,:,:), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_i4_3d

   subroutine ncgetvar_i4_4d(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      integer (i4), dimension(:,:,:,:), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_i4_4d

   subroutine ncgetvar_r4_scalar(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      real (r4), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_r4_scalar

   subroutine ncgetvar_r4_1d(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      real (r4), dimension(:), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_r4_1d

   subroutine ncgetvar_r4_2d(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      real (r4), dimension(:,:), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_r4_2d

   subroutine ncgetvar_r4_3d(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      real (r4), dimension(:,:,:), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_r4_3d

   subroutine ncgetvar_r4_4d(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      real (r4), dimension(:,:,:,:), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_r4_4d

   subroutine ncgetvar_r8_scalar(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      real (r8), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_r8_scalar

   subroutine ncgetvar_r8_1d(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      real (r8), dimension(:), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_r8_1d

   subroutine ncgetvar_r8_2d(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      real (r8), dimension(:,:), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_r8_2d

   subroutine ncgetvar_r8_3d(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      real (r8), dimension(:,:,:), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_r8_3d

   subroutine ncgetvar_r8_4d(ncid, varname, var)

      integer, intent(in) :: ncid
      character (len = *), intent(in) :: varname
      real (r8), dimension(:,:,:,:), intent(out) :: var

      integer status, varid

      status = nf90_inq_varid(ncid, varname, varid)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_inq_varid: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif
      status = nf90_get_var(ncid, varid, var)
      if (status /= nf90_noerr) then
         write (*,'(4a)') 'ncgetvar: nf90_get_var: ', trim(varname), ': ', &
            trim(nf90_strerror(status))
         stop
      endif

   end subroutine ncgetvar_r8_4d

end module ncutils
