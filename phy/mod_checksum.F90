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

   interface crc32
      module procedure crc32_1d_integer, crc32_2d_r8
   end interface crc32

   public :: csdiag, chksum, chksummsk

contains

   ! ---------------------------------------------------------------------------
   ! Private procedures.
   ! ---------------------------------------------------------------------------

   function crc32_1d_integer(iarr)

      integer, dimension(:), intent(in) :: iarr

      integer :: crc32_1d_integer

      integer :: crcfast
      external :: crcfast

      real(r8), dimension((size(iarr) + 1)/2) :: rarr

      rarr = transfer(iarr, rarr)

      crc32_1d_integer = crcfast(rarr, size(iarr)*4)

   end function crc32_1d_integer

   function crc32_2d_r8(rarr)

      real(r8), dimension(:,:), intent(in) :: rarr

      integer :: crc32_2d_r8

      integer :: crcfast
      external :: crcfast

      crc32_2d_r8 = crcfast(rarr, size(rarr)*8)

   end function crc32_2d_r8

   ! ---------------------------------------------------------------------------
   ! Public procedures.
   ! ---------------------------------------------------------------------------

   subroutine chksum(a, kcsd, text)
   ! ---------------------------------------------------------------------------
   ! Compute checksum of model field.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: kcsd
      real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy, kcsd), &
         intent(in) :: a
      character(len = *), intent(in) :: text

      real(r8), dimension(itdm, jtdm) :: aa
      integer, dimension(kcsd) :: cslist
      integer :: kcs

      do kcs = 1, kcsd
         call xcaget(aa, a(1 - nbdy, 1 - nbdy, kcs), 1)
         cslist(kcs) = crc32(aa)
      enddo

      if (mnproc == 1) then
         if (kcsd == 1) then
            write (lp,'(3a,z8.8)') ' chksum: ', text, ': 0x', cslist(1)
         else
            write (lp,'(3a,z8.8)') ' chksum: ', text, ': 0x', crc32(cslist)
         endif
      endif

   end subroutine chksum

   subroutine chksummsk(a, msk, kcsd, text)
   ! ---------------------------------------------------------------------------
   ! Compute checksum of model field after multiplying with mask.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: kcsd
      real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy, kcsd), &
         intent(in) :: a
      integer, dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy), &
         intent(in) :: msk
      character(len = *), intent(in) :: text

      real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy) :: amsk
      real(r8), dimension(itdm, jtdm) :: aa
      integer, dimension(kcsd) :: cslist
      integer :: ics, jcs, kcs

      do kcs = 1, kcsd
        !$omp parallel do private(ics)
        do jcs = 1, jj
          do ics = 1, ii
            if (msk(ics, jcs) == 0) then
              amsk(ics, jcs) = 0._r8
            else
              amsk(ics, jcs) = a(ics, jcs, kcs)
            endif
          enddo
        enddo
        !$omp end parallel do
        call xcaget(aa, amsk, 1)
        cslist(kcs) = crc32(aa)
      enddo

      if (mnproc == 1) then
         if (kcsd == 1) then
            write (lp,'(3a,z8.8)') ' chksum: ', text, ': 0x', cslist(1)
         else
            write (lp,'(3a,z8.8)') ' chksum: ', text, ': 0x', crc32(cslist)
         endif
      endif

   end subroutine chksummsk

end module mod_checksum
