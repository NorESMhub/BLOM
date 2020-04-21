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

      subroutine getfrc
c
c --- ------------------------------------------------------------------
c --- Get forcing
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mod_ben02, only: getfrc_ben02clim, getfrc_ben02syn
c
      implicit none
c
#include "common_blocks.h"
c
      if     (expcnf.eq.'cesm') then
        call getfrc_cesm
      elseif (expcnf.eq.'ben02clim') then
        call getfrc_ben02clim
      elseif (expcnf.eq.'ben02syn') then
        call getfrc_ben02syn
      elseif (expcnf.eq.'isomip1') then
c       call getfrc_isomip1
      elseif (expcnf.eq.'isomip2') then
c       call getfrc_isomip2
      else
        if (mnproc.eq.1) then
          write (lp,'(3a)') ' expcnf=',trim(expcnf),' is unsupported!'
        endif
        call xcstop('(getfrc)')
               stop '(getfrc)'
      endif
c
      return
      end
