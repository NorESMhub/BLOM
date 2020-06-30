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

      subroutine inifrc_cesm
c
c --- ------------------------------------------------------------------
c --- Initialize climatological fields for surface restoring and
c --- interpolation of CESM forcing fields
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mod_ben02, only: initai, rdcsic, rdctsf, fnlzai
c
      implicit none
c
#include "common_forc.h"
#include "common_cesm.h"
c
c --- If SST restoring is requested, prepare interpolation of surface
c --- fields and read climatological sea-ice concentration and surface
c --- temperature.
      if (trxday.gt.0.) then
        call initai
        call rdcsic
        call rdctsf
      endif
c
c --- If SSS restoring is requested, read climatological sea surface
c --- salinity.
      if (srxday.gt.0.) call rdcsss
c
c --- Initialize diagnosing/application of relaxation fluxes
      call idarlx
c
c --- Deallocate memory used for interpolation of surface fields.
      if (trxday.gt.0.) then
        call fnlzai
      endif
c
c --- Initialize time level indexes
      l1ci=1
      l2ci=1
c
      return
      end
