! ------------------------------------------------------------------------------
! Copyright (C) 2008 Mats Bentsen
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

      subroutine modeltime(ymd, tod)

      ! Uses modules

      use mod_xc
      
      implicit none

      ! Input/output arguments

      integer, intent(out)    :: ymd, tod

      ! Include files

#include "common_blocks.h"
#include "common_clndr.h"

      ! Local variables


      ! External functions


      ymd = nyear*10000 + nmonth*100 + nday
      tod = nint(mod(nstep,nstep_in_day)*baclin)

      end subroutine modeltime
