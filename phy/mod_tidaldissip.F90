! ------------------------------------------------------------------------------
! Copyright (C) 2015-2025 Mats Bentsen, Mehmet Ilicak
!
! This file is part of BLOM.
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
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

module mod_tidaldissip
! ------------------------------------------------------------------------------
! This module contains variables and procedures related to tidal wave energy
! dissipation.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_constants, only: spval
   use mod_xc
   use mod_checksum, only: csdiag, chksum
   use mod_utility, only: fnmlen
   use netcdf

   implicit none

   private

   real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy) :: &
      twedon ! Tidal wave energy dissipation over buoyancy frequency [g s-2].

   character(len = fnmlen) :: &
      tdfile ! Name of file containing tidal wave energy dissipation divided by
             ! by bottom buoyancy frequency.

   public :: twedon, tdfile, inivar_tidaldissip, read_tidaldissip
   
contains

   subroutine inivar_tidaldissip
   ! ---------------------------------------------------------------------------
   ! Initialize variables related to forcing.
   ! ---------------------------------------------------------------------------

      integer :: i, j, l

      twedon(:,:) = spval

   !$omp parallel do private(l, i)
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
            twedon(i, j) = 0._r8
         enddo
         enddo
      enddo
   !$omp end parallel do

   end subroutine inivar_tidaldissip

   subroutine read_tidaldissip
   ! ---------------------------------------------------------------------------
   ! Read tidal wave energy dissipation divided by bottom buoyancy frequency.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(itdm,jtdm) :: tmpg
      integer :: i, j, l, errstat, ncid, dimid, varid

      if (mnproc == 1) then
         write (lp, '(2a)') ' reading tidal dissipation data from ', &
                            trim(tdfile)
         call flush(lp)

         ! Open netcdf file.
         errstat = nf90_open(tdfile, nf90_nowrite, ncid)
         if (errstat /= nf90_noerr) then
          write(lp, '(4a)') ' nf90_open: ', trim(tdfile), ': ', &
                            nf90_strerror(errstat)
          call xchalt('(read_tidaldissip)')
                 stop '(read_tidaldissip)'
         endif

         ! Check dimensions.
         errstat = nf90_inq_dimid(ncid, 'x', dimid)
         if (errstat /= nf90_noerr) then
            write(lp, '(2a)') ' nf90_inq_dimid: x: ', nf90_strerror(errstat)
            call xchalt('(read_tidaldissip)')
                   stop '(read_tidaldissip)'
         endif
         errstat = nf90_inquire_dimension(ncid, dimid, len = i)
         if (errstat /= nf90_noerr) then
            write(lp, '(2a)') ' nf90_inquire_dimension: x: ', &
                              nf90_strerror(errstat)
            call xchalt('(read_tidaldissip)')
                   stop '(read_tidaldissip)'
         endif
         errstat = nf90_inq_dimid(ncid, 'y', dimid)
         if (errstat /= nf90_noerr) then
            write(lp, '(2a)') ' nf90_inq_dimid: y: ', nf90_strerror(errstat)
            call xchalt('(read_tidaldissip)')
                   stop '(read_tidaldissip)'
         endif
         errstat = nf90_inquire_dimension(ncid, dimid, len = j)
         if (errstat /= nf90_noerr) then
            write(lp, '(2a)') ' nf90_inquire_dimension: y: ', &
                              nf90_strerror(errstat)
            call xchalt('(read_tidaldissip)')
                   stop '(read_tidaldissip)'
         endif
         if (i /= itdm .or. j /= jtdm) then
            write (lp, '(2a)') ' wrong dimensions in ', trim(tdfile)
            call xchalt('(read_tidaldissip)')
                   stop '(read_tidaldissip)'
         endif

         ! Read tidal dissipation data.
         errstat = nf90_inq_varid(ncid, 'twedon', varid)
         if (errstat /= nf90_noerr) then
            write(lp, '(2a)') ' nf90_inq_varid: twedon: ', &
                              nf90_strerror(errstat)
            call xchalt('(read_tidaldissip)')
                   stop '(read_tidaldissip)'
         endif
         errstat = nf90_get_var(ncid, varid, tmpg)
         if (errstat /= nf90_noerr) then
            write(lp, '(2a)') ' nf90_get_var: twedon: ', nf90_strerror(errstat)
            call xchalt('(read_tidaldissip)')
                   stop '(read_tidaldissip)'
         endif

         errstat = nf90_close(ncid)
         if (errstat /= nf90_noerr) then
            write(lp, '(4a)') ' nf90_close: ', trim(tdfile), ': ', &
                             nf90_strerror(errstat)
            call xchalt('(read_tidaldissip)')
                   stop '(read_tidaldissip)'
         endif
      endif

      call xcaput(tmpg, twedon, 1)
      call xctilr(twedon,  1, 1, nbdy, nbdy, halo_ps)

      if (csdiag) then
         if (mnproc == 1) then
            write (lp, *) 'read_tidaldissip:'
         endif
         call chksum(twedon, 1, halo_ps, 'twedon')
      endif

   end subroutine read_tidaldissip

end module mod_tidaldissip
