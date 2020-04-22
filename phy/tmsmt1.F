! ------------------------------------------------------------------------------
! Copyright (C) 2005-2018 Mats Bentsen
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

      subroutine tmsmt1(m,n,mm,nn,k1m,k1n)
c
c --- save old layer thickness, temperature and salinity for time
c --- smoothing
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      integer nt
#endif
c
      integer i,j,k,l,kn
c
c$OMP PARALLEL DO PRIVATE(k,kn,l,i
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )
      do j=1,jj
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            dpold(i,j,kn)=dp(i,j,kn)
            told(i,j,k)=temp(i,j,kn)
            sold(i,j,k)=saln(i,j,kn)
#ifdef TRC
            do nt=1,ntr
              trcold(i,j,k,nt)=trc(i,j,kn,nt)
            enddo
#endif
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'tmsmt1:'
        endif
        call chksummsk(dpold,ip,2*kk,'dpold')
        call chksummsk(told,ip,kk,'told')
        call chksummsk(sold,ip,kk,'sold')
#ifdef TRC
        do nt=1,ntr
          call chksummsk(trcold(1-nbdy,1-nbdy,1,nt),ip,kk,'trcold')
        enddo
#endif
      endif
c
      return
      end
