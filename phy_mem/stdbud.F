! ------------------------------------------------------------------------------
! Copyright (C) 2007-2018 Mats Bentsen, Mehmet Ilicak
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

      subroutine stdbud(ncall,n,nn)
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "param_bud.h"
#include "common_bud.h"
c
      integer ncall,n,nn
c
      integer i,j,k,l,kn
      real q
c
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
#endif
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=0.
          util2(i,j)=0.
#if defined(TRC) && defined(TKE)
          util3(i,j)=0.
#  ifdef GLS
          util4(i,j)=0.
#  endif
#endif
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(k,kn,l,i,q)
      do j=1,jj
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            q=dp(i,j,kn)*scp2(i,j)
            util1(i,j)=util1(i,j)+saln(i,j,kn)*q
            util2(i,j)=util2(i,j)+temp(i,j,kn)*q
#if defined(TRC) && defined(TKE)
            util3(i,j)=util3(i,j)+trc(i,j,kn,itrtke)*q
#  ifdef GLS
            util4(i,j)=util4(i,j)+trc(i,j,kn,itrgls)*q
#  endif
#endif
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call xcsum(sdp(ncall,n),util1,ips)
      call xcsum(tdp(ncall,n),util2,ips)
#if defined(TRC) && defined(TKE)
      call xcsum(tkedp(ncall,n),util3,ips)
#  ifdef GLS
      call xcsum(glsdp(ncall,n),util4,ips)
#  endif
#endif
c
#ifdef TRC
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(k,kn,l,i,q)
      do j=1,jj
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            q=dp(i,j,kn)*scp2(i,j)
            util1(i,j)=util1(i,j)+trc(i,j,kn,1)*q
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call xcsum(trdp(ncall,n),util1,ips)
#endif
c
      return
      end
