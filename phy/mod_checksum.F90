! ------------------------------------------------------------------------------
! Copyright (C) 2006-2024 Mats Bentsen, Mariana Vertenstein
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

module mod_checksum

   use mod_types, only: r8
   use mod_xc

   implicit none

   private

   ! Options with default values, modifiable by namelist.
   logical :: &
      csdiag = .false. ! Flag that indicates whether checksums are written.

   public :: csdiag, chksum

contains

   ! ---------------------------------------------------------------------------
   ! Public procedures.
   ! ---------------------------------------------------------------------------

   subroutine chksum(a, kcsd, itype, text)
   ! ---------------------------------------------------------------------------
   ! Compute checksum of model field after multiplying with mask.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: kcsd, itype
      real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kcsd), &
         intent(inout) :: a
      character(len = *), intent(in) :: text

      integer :: crc

      select case (itype)
         case (halo_ps, halo_pv)
            call xccrc(crc, a, kcsd, ip, itype)
         case (halo_qs, halo_qv)
            call xccrc(crc, a, kcsd, iq, itype)
         case (halo_us, halo_uv)
            call xccrc(crc, a, kcsd, iu, itype)
         case (halo_vs, halo_vv)
            call xccrc(crc, a, kcsd, iv, itype)
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') ' chksum: itype = ', itype, ' is unsupported!'
            call xcstop('(chksum)')
                   stop '(chksum)'
      end select

      if (mnproc == 1) then
         write (lp,'(3a,z8.8)') ' chksum: ', trim(text), ': 0x', crc
         call flush(lp)
      endif

   end subroutine chksum

end module mod_checksum
