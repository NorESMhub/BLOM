! ------------------------------------------------------------------------------
! Copyright (C) 2015-2018 Mats Bentsen, Alf Grini
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

      subroutine icaux_cesm
c
c --- ------------------------------------------------------------------
c --- Define initial conditions specifically when coupled to CESM.
c --- ------------------------------------------------------------------
c
c
      use mod_xc
c
      implicit none
c
#include "common_forc.h"
c
      integer i,j,l
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          frzpot(i,j)=0.
          mltpot(i,j)=0.
          flxco2(i,j)=0.
          flxdms(i,j)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
