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

      subroutine idlage_step(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- update ideal age tracer
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
#include "param_trc.h"
#include "common_trc.h"
#include "common_clndr.h"
c
      integer i,j,k,l,kn
      real q
c
c --- set ideal age to zero in the surface layer and increment the age
c --- in the subsurface layers
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1-nbdy,jj+nbdy
        do l=1,isp(j)
        do i=ifp(j,l),ilp(j,l)
          trc(i,j,k1n,itriag)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      q=delt1/(86400.*nday_in_year)
c$OMP PARALLEL DO PRIVATE(k,kn,l,i)
      do j=1-nbdy,jj+nbdy
        do k=2,kk
          kn=k+nn
          do l=1,isp(j)
          do i=ifp(j,l),ilp(j,l)
            trc(i,j,kn,itriag)=trc(i,j,kn,itriag)+q
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
