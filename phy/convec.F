! ------------------------------------------------------------------------------
! Copyright (C) 2009-2020 Mats Bentsen
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

      subroutine convec(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Remove static instabilitites between the mixed layer and interior
c --- layers
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
c
      real, dimension(kdm+1) :: po,pn
      real, dimension(kdm) :: ttem,ssal,delp,dens,densr,uo,un
      real tdps,sdps,dps,ttmp,stmp,dtmp,q,udpn
      integer i,j,k,l,kn,kfpl,kfplo,kmix,ko
      logical done
c
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      real, dimension(ntr,kdm) :: ttrc
      real, dimension(ntr) :: trdps
      integer nt
#endif
c
      real sofsig,sig,rho
      external sofsig,sig,rho
c
      integer niter
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+ l,i,k,kn,ttem,ssal,delp,dens,densr,dps,kfpl,kfplo,tdps,sdps,q,
c$OMP+ ttmp,stmp,dtmp,done,niter,kmix
#ifdef TRC
c$OMP+ ,nt,ttrc,trdps
#endif
c$OMP+ )
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c --- --- Copy variables into 1d arrays
          do k=1,kk
            kn=k+nn
            ttem(k)=temp(i,j,kn)
            ssal(k)=saln(i,j,kn)
            delp(k)=dp(i,j,kn)
            dens(k)=sigma(i,j,kn)
            densr(k)=sigmar(i,j,k)
#ifdef TRC
            do nt=1,ntr
              ttrc(nt,k)=trc(i,j,kn,nt)
            enddo
#endif
          enddo
c
c --- ------------------------------------------------------------------
c --- --- Define first physical interior layer
c --- ------------------------------------------------------------------
c
          k=3
          dps=0.
          do while (delp(k).lt.epsil)
            dps=dps+delp(k)
            delp(k)=0.
            k=k+1
            if (k.gt.kk) exit
          enddo
          if (k.gt.kk) then
            delp(2)=delp(2)+dps
          else
            delp(k)=delp(k)+dps
          endif
          kfpl=k
          kfplo=kfpla(i,j,n)
          if (kfpl.lt.kfplo) then
            if (kfplo.le.kk) then
              tdps=0.
              sdps=0.
              dps=0.
#ifdef TRC
              do nt=1,ntr
                trdps(nt)=0.
              enddo
#endif
              do k=kfpl,kfplo
                tdps=tdps+ttem(k)*delp(k)
                sdps=sdps+ssal(k)*delp(k)
                dps=dps+delp(k)
#ifdef TRC
                do nt=1,ntr
                  trdps(nt)=trdps(nt)+ttrc(nt,k)*delp(k)
                enddo
#endif
              enddo
              q=1./dps
              ttmp=tdps*q
              stmp=sdps*q
              dtmp=sig(ttmp,stmp)
              if (dtmp.gt.densr(kfplo)) then
                do k=kfpl,kfplo-1
                  delp(k)=0.
                enddo
                kfpl=kfplo
                ttem(kfpl)=ttmp
                ssal(kfpl)=stmp
                dens(kfpl)=dtmp
                delp(kfpl)=dps
#ifdef TRC
                do nt=1,ntr
                  ttrc(nt,kfpl)=trdps(nt)*q
                enddo
#endif
              endif
            else
              tdps=0.
              sdps=0.
              dps=0.
#ifdef TRC
              do nt=1,ntr
                trdps(nt)=0.
              enddo
#endif
              do k=kfpl,kk
                tdps=tdps+ttem(k)*delp(k)
                sdps=sdps+ssal(k)*delp(k)
                dps=dps+delp(k)
#ifdef TRC
                do nt=1,ntr
                  trdps(nt)=trdps(nt)+ttrc(nt,k)*delp(k)
                enddo
#endif
                delp(k)=0.
              enddo
              q=1./dps
              ttmp=tdps*q
              stmp=sdps*q
              dtmp=sig(ttmp,stmp)
              kfpl=kk
              do while (dtmp.lt.densr(kfpl))
                if (kfpl.eq.3) exit
                kfpl=kfpl-1
              enddo
              ttem(kfpl)=ttmp
              ssal(kfpl)=stmp
              dens(kfpl)=dtmp
              delp(kfpl)=dps
#ifdef TRC
              do nt=1,ntr
                ttrc(nt,kfpl)=trdps(nt)*q
              enddo
#endif
            endif
          endif
c
          if (kfpl.le.kk) then
c
c --- ------------------------------------------------------------------
c --- ----- Remove static instabilities
c --- ------------------------------------------------------------------
c
            done=.false.
c
            niter=0
            do while (.not.done)
              niter=niter+1
              if (niter.eq.100) then
                write (lp,*) 'blom: convec: no convergence!',i+i0,j+j0
                exit
              endif
c
              done=.true.
c
c --- ------- Remove instabilities between the lower mixed layer and
c --- ------- interior layers by considering the potential density jump
c --- ------- across the mixed layer base with reference pressure at the
c --- ------- interface
              tdps=ttem(2)*delp(2)
              sdps=ssal(2)*delp(2)
              dps=delp(2)
#ifdef TRC
              do nt=1,ntr
                trdps(nt)=ttrc(nt,2)*delp(2)
              enddo
#endif
              ttmp=ttem(2)
              stmp=ssal(2)
              k=kfpl
              do while (rho(dps,ttmp,stmp).gt.
     .                  rho(dps,ttem(k),ssal(k)).or.
     .                  delp(k).lt.epsil)
                tdps=tdps+ttem(k)*delp(k)
                sdps=sdps+ssal(k)*delp(k)
                dps=dps+delp(k)
                q=1./dps
                ttmp=tdps*q
                stmp=sdps*q
#ifdef TRC
                do nt=1,ntr
                  trdps(nt)=trdps(nt)+ttrc(nt,k)*delp(k)
                enddo
#endif
                k=k+1
                if (k.gt.kk) exit
              enddo
              kmix=k-1
              if (kmix.ge.kfpl) then
                ttem(2)=ttmp
                ssal(2)=stmp
                dens(2)=sig(ttem(2),ssal(2))
#ifdef TRC
                do nt=1,ntr
                  ttrc(nt,2)=trdps(nt)*q
                enddo
#endif
                dps=0.
                do k=kfpl,kmix
                  dps=dps+delp(k)
                  delp(k)=0.
                enddo
                k=kmix
                do while (dens(2).lt.densr(k))
                  if (k.eq.3) exit
                  k=k-1
                enddo
                kfpl=k
                ttem(kfpl)=ttem(2)
                ssal(kfpl)=ssal(2)
                dens(kfpl)=dens(2)
                delp(kfpl)=dps
#ifdef TRC
                do nt=1,ntr
                  ttrc(nt,kfpl)=ttrc(nt,2)
                enddo
#endif
                do k=kfpl+1,kmix
                  ttem(k)=ttem(2)
                  dens(k)=densr(k)
                  ssal(k)=sofsig(dens(k),ttem(k))
                enddo
              endif
c
            enddo
c
          endif
c
          kfpla(i,j,n)=kfpl
c
c --- --- Copy 1d arrays to 3d arrays
          do k=1,kk
            kn=k+nn
            temp(i,j,kn)=ttem(k)
            saln(i,j,kn)=ssal(k)
            sigma(i,j,kn)=dens(k)
            dp(i,j,kn)=delp(k)
            p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
#ifdef TRC
            do nt=1,ntr
              trc(i,j,kn,nt)=ttrc(nt,k)
            enddo
#endif
          enddo
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Redistribute momentum
c --- ------------------------------------------------------------------
c
      call xctilr(p, 1,kk+1, 1,1, halo_ps)
c
c$OMP PARALLEL DO PRIVATE(l,i,k,kn,uo,po,pn,ko,un,udpn)
      do j=1,jj
c
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          do k=1,kk
            kn=k+nn
            uo(k)=u(i,j,kn)
          enddo
          po(1)=0.
          pn(1)=0.
          do k=2,kk+1
            po(k)=pu(i,j,k)
            pn(k)=.5*(min(pu(i,j,kk+1),p(i  ,j,k))
     .               +min(pu(i,j,kk+1),p(i-1,j,k)))
          enddo
c
          ko=1
          do kn=1,kk
            if (pn(kn+1)-pn(kn).eq.0.) then
              un(kn)=0.
            else
              udpn=0.
              do while (pn(kn+1).gt.po(ko+1))
                udpn=udpn+uo(ko)*(po(ko+1)-max(po(ko),pn(kn)))
                ko=ko+1
              enddo
              un(kn)=(udpn+uo(ko)*(pn(kn+1)-max(po(ko),pn(kn))))
     .               /(pn(kn+1)-pn(kn))
            endif
          enddo
          do k=1,kk
            kn=k+nn
            u(i,j,kn)=un(k)
          enddo
c
        enddo
        enddo
c
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          do k=1,kk
            kn=k+nn
            uo(k)=v(i,j,kn)
          enddo
          po(1)=0.
          pn(1)=0.
          do k=2,kk+1
            po(k)=pv(i,j,k)
            pn(k)=.5*(min(pv(i,j,kk+1),p(i,j  ,k))
     .               +min(pv(i,j,kk+1),p(i,j-1,k)))
          enddo
c
          ko=1
          do kn=1,kk
            if (pn(kn+1)-pn(kn).eq.0.) then
              un(kn)=0.
            else
              udpn=0.
              do while (pn(kn+1).gt.po(ko+1))
                udpn=udpn+uo(ko)*(po(ko+1)-max(po(ko),pn(kn)))
                ko=ko+1
              enddo
              un(kn)=(udpn+uo(ko)*(pn(kn+1)-max(po(ko),pn(kn))))
     .               /(pn(kn+1)-pn(kn))
            endif
          enddo
          do k=1,kk
            kn=k+nn
            v(i,j,kn)=un(k)
          enddo
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(k,kn,l,i,q)
      do j=1,jj
        do k=1,kk
          kn=k+nn
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            q=min(p(i,j,kk+1),p(i-1,j,kk+1))
            dpu(i,j,kn)=
     .        .5*((min(q,p(i-1,j,k+1))-min(q,p(i-1,j,k)))
     .           +(min(q,p(i  ,j,k+1))-min(q,p(i  ,j,k))))
          enddo
          enddo
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            q=min(p(i,j,kk+1),p(i,j-1,kk+1))
            dpv(i,j,kn)=
     .        .5*((min(q,p(i,j-1,k+1))-min(q,p(i,j-1,k)))
     .           +(min(q,p(i,j  ,k+1))-min(q,p(i,j  ,k))))
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c     do j=1,jj
c       do l=1,isu(j)
c       do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
c         q=0.
c         do k=1,kk
c           kn=k+nn
c           q=q+u(i,j,kn)*dpu(i,j,kn)
c         enddo
c         if (abs(q).gt.1.e-4) then
c           write (lp,*) 'convec: u imbalance:',q,i,j
c         endif
c       enddo
c       enddo
c       do l=1,isv(j)
c       do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
c         q=0.
c         do k=1,kk
c           kn=k+nn
c           q=q+v(i,j,kn)*dpv(i,j,kn)
c         enddo
c         if (abs(q).gt.1.e-4) then
c           write (lp,*) 'convec: v imbalance:',q,i,j
c         endif
c       enddo
c       enddo
c     enddo
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'convec:'
        endif
        call chksummsk(dp,ip,2*kk,'dp')
        call chksummsk(temp,ip,2*kk,'temp')
        call chksummsk(saln,ip,2*kk,'saln')
        call chksummsk(sigma,ip,2*kk,'sigma')
        call chksummsk(u,iu,2*kk,'u')
        call chksummsk(v,iv,2*kk,'v')
#ifdef TRC
        do nt=1,ntr
          call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,2*kk,'trc')
        enddo
#endif
      endif
c
      return
      end
