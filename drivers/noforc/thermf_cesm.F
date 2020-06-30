! ------------------------------------------------------------------------------
! Copyright (C) 2009-2015 Mats Bentsen
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

      subroutine thermf_cesm(m,n,mm,nn,k1m,k1n)
c
c --- NERSC version of thermf. To be used when coupled to CESM
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
#include "common_forc.h"
c
      integer i,j,k,l
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      integer nt
#endif
c
c$OMP PARALLEL DO PRIVATE(l,i
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
          rnfflx(i,j)=0.
          rfiflx(i,j)=0.
          brnflx(i,j)=0.
          frzpot(i,j)=0.
          mltpot(i,j)=0.
          hmltfz(i,j)=0.
          surflx(i,j)=0.
          sswflx(i,j)=0.
          surrlx(i,j)=0.
          salrlx(i,j)=0.
          salflx(i,j)=0.
#ifdef TRC
          do nt=1,ntrc
            trflx(nt,i,j)=0.
          enddo
#endif
          ustar(i,j)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
