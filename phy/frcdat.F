! ------------------------------------------------------------------------------
! Copyright (C) 2002-2015 Mats Bentsen
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

      block data frcdat
c
      use mod_xc
c
      implicit none
c
#include "common_forc.h"
c
      data
     . albw_d   /.065/    ! --                 albedo over water; diff light
     .,rhowat   /1024./   ! kg / m^3           app. density of surf. sea water
     .,sref     /34.65/   ! g / kg             global ref. surface salinity
     .,nrfets   /7/       ! days               e-folding time scale for
                          !                    detrainment of runoff reservoar
c
      end
