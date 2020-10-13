! ------------------------------------------------------------------------------
! Copyright (C) 2020 Aleksi Nummelin, Mats Bentsen
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

      subroutine icaux_channel
c
c --- ------------------------------------------------------------------
c --- THIS IS AN EMPTY FUNCTION - COULD BE USED TO SETUP THE SEA ICE
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mod_ben02, only: ntda
c
      implicit none
c
#include "common_blocks.h"
#include "common_forc.h"
c#include "common_seaice.h"
c
      ntda=0
      if (mnproc.eq.1) then
          write(lp,*) 'ntda setup done!'
      endif
c
      return
      end
