! ------------------------------------------------------------------------------
! Copyright (C) 2006-2008 Mats Bentsen
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

      subroutine chksum(a,kcsd,text)
c
c --- compute checksum of model field
c
      use mod_xc
c
      implicit none
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kcsd) :: a
      integer kcsd
      character text*(*)
c
      real, dimension(itdm,jtdm,kcsd) :: aa
      integer kcs
c
      integer crcfast
      external crcfast
c
      do kcs=1,kcsd
        call xcaget(aa(1,1,kcs),a(1-nbdy,1-nbdy,kcs),1)
      enddo
c
      if (mnproc.eq.1) then
        write (lp,'(3a,z8.8)') ' chksum: ',text,': 0x',
     .    crcfast(aa,itdm*jtdm*kcsd*8)
      endif
c
      return
      end
