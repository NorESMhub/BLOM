! ------------------------------------------------------------------------------
! Copyright (C) 2012 Mats Bentsen
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

      subroutine getprecipfact_mct(lsend_precip_fact, precip_fact)

      ! Uses modules
      
      use types, only : r8
      use mod_xc

      implicit none

      ! Input/output arguments

      logical, intent(out)  :: lsend_precip_fact
      real(r8), intent(out) :: precip_fact

      ! Include files

#include "common_forc.h"

      lsend_precip_fact = sprfac
      precip_fact = prfac

      end subroutine getprecipfact_mct
