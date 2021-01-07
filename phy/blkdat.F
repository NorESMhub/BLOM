! ------------------------------------------------------------------------------
! Copyright (C) 2000 HYCOM Consortium and contributors
! Copyright (C) 2001-2020 Mats Bentsen
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

      block data blkdat
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
c
c --- 'mltmin' = minimum mixed-layer thickness (m)
c --- 'thktop' = maximum thickness of top layer (m)
      data mltmin/5./,thktop/10./
c
c --- slip=+1  for free-slip boundary cond., slip=-1  for non-slip cond.
c --- 'thkbot' = thickness of bottom boundary layer (m)
      data slip/-1./,thkbot/10./
c
c --- weights for time smoothing
ccc      data wuv1,wuv2/.5,.25/
      data wuv1,wuv2/.75,.125/
      data wts1,wts2/.875,.0625/
ccc      data wts1,wts2/.9375,.03125/
ccc      data wts1,wts2/.96875,.015625/
      data wbaro/.125/
ccc      data wbaro/.0625/
      data wpgf/.25/
c
c --- 'itest,jtest' = grid point where detailed diagnostics are desired
      data itest,jtest/60,60/
c
c ---  s w i t c h e s    (if set to .true., then...)
c
      data csdiag/.false./,cnsvdi/.true./
c
      end
