! ------------------------------------------------------------------------------
! Copyright (C) 2006-2014 Mats Bentsen, Alok Kumar Gupta
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

      real function gettime()
c
c --- ------------------------------------------------------------------
c --- return time in seconds since last call to either initime of
c --- gettime
c --- ------------------------------------------------------------------
c
      use mod_xc, only: mnproc
      implicit none
c
#include "timing.h"
c
      if (mnproc.eq.1) then
        gettime=wtime()-wtimeold
        wtimeold=gettime+wtimeold
      else
        gettime=1.
      endif
c
      return
      end
