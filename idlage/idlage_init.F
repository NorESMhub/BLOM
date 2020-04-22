! ------------------------------------------------------------------------------
! Copyright (C) 2010-2018 Mats Bentsen
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

      subroutine idlage_init
c
c --- ------------------------------------------------------------------
c --- initialization of ideal age tracer
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "param_trc.h"
#include "common_trc.h"
c
      integer i,j,k,l
c
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=1-nbdy,jj+nbdy
        do k=1,kk
          do l=1,isp(j)
          do i=ifp(j,l),ilp(j,l)
            trc(i,j,k   ,itriag)=0.
            trc(i,j,k+kk,itriag)=0.
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
