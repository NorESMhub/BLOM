! ------------------------------------------------------------------------------
! Copyright (C) 2005-2019 Mats Bentsen
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

      subroutine tmsmt2(m,n,mm,nn,k1m,k1n)
c
c --- time smoothing of layer thickness, temperature and salinity
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
c
      real, dimension(1-nbdy:idm+nbdy) :: pbfaco,pbfacn
      integer i,j,k,l,kn,km,kp
      real pold,pmid,pnew,q
c
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      integer nt
#endif
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+ l,i,pbfaco,pbfacn,k,kn,km,kp,pold,pmid,pnew,q
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          pbfaco(i)=0.
          pbfacn(i)=0.
        enddo
        enddo
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            pbfaco(i)=pbfaco(i)+dpold(i,j,kn)
            pbfacn(i)=pbfacn(i)+dp(i,j,kn)
          enddo
          enddo
        enddo
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          pbfaco(i)=pb(i,j,m)/pbfaco(i)
          pbfacn(i)=pb(i,j,m)/pbfacn(i)
        enddo
        enddo
        do k=1,kk
          km=k+mm
          kn=k+nn
          kp=min(k+1,kk)
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            pold=max(0.,dpold(i,j,kn)*pbfaco(i))
            pmid=max(0.,dp(i,j,km))
            pnew=max(0.,dp(i,j,kn)*pbfacn(i))
            dp(i,j,km)=wts1*pmid+wts2*(pold+pnew)
            pold=pold+epsil
            pmid=pmid+epsil
            pnew=pnew+epsil
            temp(i,j,km)=(wts1*pmid*temp(i,j,km)
     .                   +wts2*(pold*told(i,j,k)+pnew*temp(i,j,kn)))
     .                   /(dp(i,j,km)+epsil)
            saln(i,j,km)=(wts1*pmid*saln(i,j,km)
     .                   +wts2*(pold*sold(i,j,k)+pnew*saln(i,j,kn)))
     .                   /(dp(i,j,km)+epsil)
#ifdef TRC
            do nt=1,ntr
              trc(i,j,km,nt)=(wts1*pmid*trc(i,j,km,nt)
     .                       +wts2*(pold*trcold(i,j,k,nt)
     .                             +pnew*trc(i,j,kn,nt)))
     .                       /(dp(i,j,km)+epsil)
            enddo
#endif
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call xctilr(dp(1-nbdy,1-nbdy,k1m), 1,kk, 3,3, halo_ps)
c
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=-2,jj+2
        do k=1,kk
          do l=1,isp(j)
          do i=max(-2,ifp(j,l)),min(ii+2,ilp(j,l))
            p(i,j,k+1)=p(i,j,k)+dp(i,j,k+mm)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(k,km,l,i,q)
      do j=-1,jj+2
        do k=1,kk
          km=k+mm
          do l=1,isu(j)
          do i=max(-1,ifu(j,l)),min(ii+2,ilu(j,l))
            q=min(p(i,j,kk+1),p(i-1,j,kk+1))
            dpu(i,j,km)=
     .        .5*((min(q,p(i-1,j,k+1))-min(q,p(i-1,j,k)))
     .           +(min(q,p(i  ,j,k+1))-min(q,p(i  ,j,k))))
          enddo
          enddo
          do l=1,isv(j)
          do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            q=min(p(i,j,kk+1),p(i,j-1,kk+1))
            dpv(i,j,km)=
     .        .5*((min(q,p(i,j-1,k+1))-min(q,p(i,j-1,k)))
     .           +(min(q,p(i,j  ,k+1))-min(q,p(i,j  ,k))))
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'tmsmt2:'
        endif
        call chksummsk(dp,ip,2*kk,'dp')
        call chksummsk(temp,ip,2*kk,'temp')
        call chksummsk(saln,ip,2*kk,'saln')
        call chksummsk(dpu,iu,2*kk,'dpu')
        call chksummsk(dpv,iv,2*kk,'dpv')
#ifdef TRC
        do nt=1,ntr
          call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,2*kk,'trc')
        enddo
#endif
      endif
c
      return
      end
