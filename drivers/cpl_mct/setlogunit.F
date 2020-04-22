! ------------------------------------------------------------------------------
! Copyright (C) 2015 Mats Bentsen
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

      subroutine setlogunit

      ! Uses modules

      use mod_xc
      use shr_file_mod, only: shr_file_getUnit, shr_file_setIO

      implicit none

      logical exists

      ! ----------------------------------------------------------------
      ! Redirect standard out to a log file if requested
      ! ----------------------------------------------------------------
      if (mnproc.eq.1) then
         inquire(file = 'ocn_modelio.nml', exist = exists)
         if (exists) then
            lp = shr_file_getUnit()
            call shr_file_setIO('ocn_modelio.nml', lp)
         endif
      endif

      end subroutine setlogunit
