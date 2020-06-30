! ------------------------------------------------------------------------------
! Copyright (C) 2007-2015 Mats Bentsen
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

      subroutine updtrc(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- update tracers due to non-passive processes
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "param_trc.h"
#include "common_trc.h"
c
      integer m,n,mm,nn,k1m,k1n
c
#ifdef HAMOCC
      call hamocc_step(m,n,mm,nn,k1m,k1n)
#endif
#ifdef IDLAGE
      call idlage_step(m,n,mm,nn,k1m,k1n)
#endif
c
      return
      end
