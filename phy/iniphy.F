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

      subroutine iniphy
c
c --- ------------------------------------------------------------------
c --- Initialize physical parameterizations
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
c
      if     (expcnf.eq.'cesm'.or.
     .        expcnf.eq.'ben02clim'.or.expcnf.eq.'ben02syn') then
        call read_twedon
      elseif (expcnf.eq.'isomip1') then
      elseif (expcnf.eq.'isomip2') then
      else
        if (mnproc.eq.1) then
          write (lp,'(3a)') ' expcnf=',trim(expcnf),' is unsupported!'
        endif
        call xcstop('(iniphy)')
               stop '(iniphy)'
      endif
c
      return
      end
