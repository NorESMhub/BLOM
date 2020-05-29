! ------------------------------------------------------------------------------
! Copyright (C) 2007-2018 Mats Bentsen
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

      subroutine restart_trcrd(rstfnm_ocn)
c
c --- ------------------------------------------------------------------
c --- Read tracer state from restart files
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      character rstfnm_ocn*(*)
c
      call restart_ocntrcrd(rstfnm_ocn)
c
#ifdef HAMOCC
      call hamocc_init(1,rstfnm_ocn)
#endif
c
      return
      end
