! ------------------------------------------------------------------------------
! Copyright (C) 2006-2020 Mats Bentsen
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

   ! Constants.
   logical :: &
      csdiag = .false. ! Flag that indicates whether checksums are written.

   integer :: crcfast
   external :: crcfast

   public :: csdiag, chksum, chksummsk

contains

   subroutine chksum(a, kcsd, text)
   ! ---------------------------------------------------------------------------
   ! Compute checksum of model field.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: kcsd
      real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy, kcsd), &
         intent(in) :: a
      character(len = *), intent(in) :: text

      real(r8), dimension(itdm, jtdm, kcsd) :: aa
      integer :: kcs

      do kcs = 1, kcsd
        call xcaget(aa(1, 1, kcs), a(1 - nbdy, 1 - nbdy, kcs), 1)
      enddo

      if (mnproc == 1) then
         write (lp,'(3a,z8.8)') ' chksum: ', text, ': 0x', &
            crcfast(aa, itdm*jtdm*kcsd*8)
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

      real(r8), dimension(itdm, jtdm, kcsd) :: aa
      real(r8), dimension(itdm, jtdm) :: rrmsk
      real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy) :: rmsk
      integer :: ics, jcs, kcs

      do jcs = 1, jj
         do ics = 1, ii
            rmsk(ics, jcs) = msk(ics, jcs)
         enddo
      enddo

      do kcs = 1, kcsd
         call xcaget(aa(1 , 1, kcs), a(1 - nbdy, 1 - nbdy, kcs), 1)
      enddo
      call xcaget(rrmsk, rmsk, 1)

      if (mnproc == 1) then
         do kcs = 1, kcsd
            do jcs = 1, jtdm
              do ics = 1, itdm
                 if (rrmsk(ics, jcs) < .5_r8) then
                    aa(ics, jcs, kcs) = 0._r8
                 endif
               enddo
            enddo
         enddo
         write (lp,'(3a,z8.8)') ' chksum: ', text, ': 0x', &
            crcfast(aa, itdm*jtdm*kcsd*8)
      endif

   end subroutine chksummsk

end module mod_checksum
