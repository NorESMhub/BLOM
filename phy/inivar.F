! ------------------------------------------------------------------------------
! Copyright (C) 2015-2018 Mats Bentsen, Jerry Tjiputra
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

      subroutine inivar
c
c --- ------------------------------------------------------------------
c --- Initialize various arrays
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "common_forc.h"
#include "common_seaice.h"
c
      integer i,j,k,l
c
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      integer nt
#endif
c
c --- ------------------------------------------------------------------
c --- initialize some arrays
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(i,k
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          p(i,j,1)=huge
          pu(i,j,1)=huge
          pv(i,j,1)=huge
          phi(i,j,1)=huge
          do k=1,3
            ubflxs(i,j,k)=huge
            vbflxs(i,j,k)=huge
          enddo
          do k=1,2
            pb(i,j,k)=huge
            ubflx(i,j,k)=huge
            vbflx(i,j,k)=huge
            pb_mn(i,j,k)=huge
            ubflx_mn(i,j,k)=huge
            vbflx_mn(i,j,k)=huge
            ub(i,j,k)=huge
            vb(i,j,k)=huge
            ubflxs_p(i,j,k)=huge
            vbflxs_p(i,j,k)=huge
            pbu(i,j,k)=huge
            pbv(i,j,k)=huge
            pvtrop(i,j,k)=huge
          enddo
          pb_p(i,j)=huge
          pbu_p(i,j)=huge
          pbv_p(i,j)=huge
          ubcors_p(i,j)=huge
          vbcors_p(i,j)=huge
          defor1(i,j)=huge
          defor2(i,j)=huge
          pvtrop_o(i,j)=huge
          utotm(i,j)=huge
          vtotm(i,j)=huge
          utotn(i,j)=huge
          vtotn(i,j)=huge
          uflux (i,j)=huge
          vflux (i,j)=huge
          uflux2(i,j)=huge
          vflux2(i,j)=huge
          uflux3(i,j)=huge
          vflux3(i,j)=huge
          uja(i,j)=huge
          ujb(i,j)=huge
          via(i,j)=huge
          vib(i,j)=huge
          pgfxm(i,j,1)=huge
          pgfxm(i,j,2)=huge
          pgfxm_o(i,j)=huge
          xixp(i,j,1)=huge
          xixp(i,j,2)=huge
          xixp_o(i,j)=huge
          xixm(i,j,1)=huge
          xixm(i,j,2)=huge
          xixm_o(i,j)=huge
          pgfym(i,j,1)=huge
          pgfym(i,j,2)=huge
          pgfym_o(i,j)=huge
          xiyp(i,j,1)=huge
          xiyp(i,j,2)=huge
          xiyp_o(i,j)=huge
          xiym(i,j,1)=huge
          xiym(i,j,2)=huge
          xiym_o(i,j)=huge
          difwgt(i,j)=huge
          sealv(i,j)=huge
          twedon(i,j)=huge
          taux(i,j)=huge
          tauy(i,j)=huge
          swa(i,j)=huge
          nsf(i,j)=huge
          hmltfz(i,j)=huge
          hmlt(i,j)=huge
          dfl(i,j)=huge
          lip(i,j)=huge
          sop(i,j)=huge
          eva(i,j)=huge
          rnf(i,j)=huge
          rfi(i,j)=huge
          fmltfz(i,j)=huge
          sfl(i,j)=huge
          ztx(i,j)=huge
          mty(i,j)=huge
          ustarw(i,j)=huge
          slp(i,j)=huge
          abswnd(i,j)=huge
          albw(i,j)=huge
          frzpot(i,j)=huge
          mltpot(i,j)=huge
          atmco2(i,j)=huge
          flxco2(i,j)=huge
          flxdms(i,j)=huge
          alb(i,j)=huge
          rnfres(i,j)=huge
          rnfflx(i,j)=huge
          rfiflx(i,j)=huge
          eiacc(i,j)=huge
          pracc(i,j)=huge
          ficem(i,j)=huge
          hicem(i,j)=huge
          hsnwm(i,j)=huge
          ustari(i,j)=huge
          tauxice(i,j)=huge
          tauyice(i,j)=huge
          uicem(i,j)=huge
          vicem(i,j)=huge
          iagem(i,j)=huge
          do k=1,kk
            u  (i,j,k   )=huge
            u  (i,j,k+kk)=huge
            v  (i,j,k   )=huge
            v  (i,j,k+kk)=huge
            uflx(i,j,k   )=huge
            uflx(i,j,k+kk)=huge
            utflx(i,j,k   )=huge
            utflx(i,j,k+kk)=huge
            usflx(i,j,k   )=huge
            usflx(i,j,k+kk)=huge
            vflx(i,j,k   )=huge
            vflx(i,j,k+kk)=huge
            vtflx(i,j,k   )=huge
            vtflx(i,j,k+kk)=huge
            vsflx(i,j,k   )=huge
            vsflx(i,j,k+kk)=huge
            umfltd(i,j,k   )=huge
            umfltd(i,j,k+kk)=huge
            vmfltd(i,j,k   )=huge
            vmfltd(i,j,k+kk)=huge
            dp (i,j,k   )=huge
            dp (i,j,k+kk)=huge
            dpu(i,j,k   )=huge
            dpu(i,j,k+kk)=huge
            dpuold(i,j,k)=huge
            dpv(i,j,k   )=huge
            dpv(i,j,k+kk)=huge
            dpvold(i,j,k)=huge
            p (i,j,k+1)=huge
            pu(i,j,k+1)=huge
            pv(i,j,k+1)=huge
            phi(i,j,k+1)=huge
c
            temp(i,j,k)=huge
            saln(i,j,k)=huge
            sigma(i,j,k   )=huge
            sigma(i,j,k+kk)=huge
            pgfx(i,j,k   )=huge
            pgfx(i,j,k+kk)=huge
            pgfy(i,j,k   )=huge
            pgfy(i,j,k+kk)=huge
            pgfxo(i,j,k)=huge
            pgfyo(i,j,k)=huge
#ifdef TRC
            do nt=1,ntr
              trc(i,j,k   ,nt)=huge
              trc(i,j,k+kk,nt)=huge
            enddo
#endif
          enddo
          do k=1,4
            uml(i,j,k)=huge
            vml(i,j,k)=huge
          enddo
          do k=1,2
            umlres(i,j,k)=huge
            vmlres(i,j,k)=huge
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,i,k)
      do j=1,jj+1
        do l=1,isq(j)
        do i=max(1,ifq(j,l)),min(ii+1,ilq(j,l))
          do k=1,2
            pb(i  ,j  ,k)=0.
            pb(i-1,j  ,k)=0.
            pb(i  ,j-1,k)=0.
            pb(i-1,j-1,k)=0.
            pb_mn(i  ,j  ,k)=0.
            pb_mn(i-1,j  ,k)=0.
            pb_mn(i  ,j-1,k)=0.
            pb_mn(i-1,j-1,k)=0.
          enddo
          pb_p(i  ,j  )=0.
          pb_p(i-1,j  )=0.
          pb_p(i  ,j-1)=0.
          pb_p(i-1,j-1)=0.
          p(i  ,j  ,1)=0.
          p(i-1,j  ,1)=0.
          p(i  ,j-1,1)=0.
          p(i-1,j-1,1)=0.
          do k=1,kk
            dp(i  ,j  ,k   )=0.
            dp(i  ,j  ,k+kk)=0.
            dp(i-1,j  ,k   )=0.
            dp(i-1,j  ,k+kk)=0.
            dp(i  ,j-1,k   )=0.
            dp(i  ,j-1,k+kk)=0.
            dp(i-1,j-1,k   )=0.
            dp(i-1,j-1,k+kk)=0.
          enddo
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call xctilr(pb,   1,   2, nbdy,nbdy, halo_ps)
      call xctilr(pb_mn,1,   2, nbdy,nbdy, halo_ps)
      call xctilr(pb_p, 1,   1, nbdy,nbdy, halo_ps)
      call xctilr(p,    1,   1, nbdy,nbdy, halo_ps)
      call xctilr(dp,   1,2*kk, nbdy,nbdy, halo_ps)
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1-nbdy,jj+nbdy
        do l=1,isp(j)
        do i=ifp(j,l),ilp(j,l)
          p(i,j,1)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          defor2(i  ,j  )=0.
          defor2(i+1,j  )=0.
          defor2(i  ,j+1)=0.
          defor2(i+1,j+1)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,i,k)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          twedon(i,j)=0.
          ustari(i,j)=0.
          tauxice(i,j)=0.
          tauyice(i,j)=0.
          uicem(i,j)=0.
          vicem(i,j)=0.
        enddo
        enddo
        do k=1,4
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
            uml(i,j,k)=0.
          enddo
          enddo
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            vml(i,j,k)=0.
          enddo
          enddo
        enddo
        do k=1,2
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            umlres(i,j,k)=0.
          enddo
          enddo
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            vmlres(i,j,k)=0.
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- initialize  u,ub,utotm,uflx,uflux,uflux2/3,uja,ujb  at points
c --- located upstream and downstream (in i direction) of p points.
c --- initialize  pbu,dpu,utotn  upstream and downstream of p
c --- points as well as at lateral neighbors of interior u points.
c
c$OMP PARALLEL DO PRIVATE(l,i,k)
      do j=0,jj+1
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          pu(i,j,1)=0.
          utotn(i,j-1)=0.
          utotn(i,j+1)=0.
          do k=1,2
            pbu(i,j-1,k)=0.
            pbu(i,j+1,k)=0.
          enddo
          pbu_p(i,j-1)=0.
          pbu_p(i,j+1)=0.
          do k=1,kk
            dpu(i,j-1,k   )=0.
            dpu(i,j-1,k+kk)=0.
            dpuold(i,j-1,k)=0.
            dpu(i,j+1,k   )=0.
            dpu(i,j+1,k+kk)=0.
            dpuold(i,j+1,k)=0.
          enddo
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,i,k)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l)+1)
          do k=1,3
            ubflxs(i,j,k)=0.
          enddo
          do k=1,2
            ubflx(i,j,k)=0.
            ubflx_mn(i,j,k)=0.
            ub(i,j,k)=0.
            ubflxs_p(i,j,k)=0.
            pbu(i,j,k)=0.
          enddo
          pbu_p(i,j)=0.
          ubcors_p(i,j)=0.
          utotn (i,j)=0.
          utotm(i,j)=0.
          uflux (i,j)=0.
          uflux2(i,j)=0.
          uflux3(i,j)=0.
          uja(i,j)=0.
          ujb(i,j)=0.
          do k=1,kk
            dpu(i,j,k   )=0.
            dpu(i,j,k+kk)=0.
            dpuold(i,j,k)=0.
            uflx(i,j,k   )=0.
            uflx(i,j,k+kk)=0.
            utflx(i,j,k   )=0.
            utflx(i,j,k+kk)=0.
            usflx(i,j,k   )=0.
            usflx(i,j,k+kk)=0.
            umfltd(i,j,k   )=0.
            umfltd(i,j,k+kk)=0.
            u(i,j,k   )=0.
            u(i,j,k+kk)=0.
          enddo
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call xctilr(pu,     1,   1, nbdy,nbdy, halo_us)
      call xctilr(utotn,  1,   1, nbdy,nbdy, halo_us)
      call xctilr(pbu,    1,   2, nbdy,nbdy, halo_us)
      call xctilr(pbu_p,  1,   1, nbdy,nbdy, halo_us)
      call xctilr(dpu,    1,2*kk, nbdy,nbdy, halo_us)
      call xctilr(dpuold, 1,  kk, nbdy,nbdy, halo_us)
      call xctilr(ub,     1,   2, nbdy,nbdy, halo_us)
      call xctilr(ubflx,  1,   2, nbdy,nbdy, halo_us)
      call xctilr(ubflx_mn, 1,   2, nbdy,nbdy, halo_us)
      call xctilr(ubflxs, 1,   3, nbdy,nbdy, halo_us)
      call xctilr(ubflxs_p, 1, 2, nbdy,nbdy, halo_us)
      call xctilr(ubcors_p, 1, 1, nbdy,nbdy, halo_us)
      call xctilr(utotm,  1,   1, nbdy,nbdy, halo_us)
      call xctilr(uflux,  1,   1, nbdy,nbdy, halo_us)
      call xctilr(uflux2, 1,   1, nbdy,nbdy, halo_us)
      call xctilr(uflux3, 1,   1, nbdy,nbdy, halo_us)
      call xctilr(uja,    1,   1, nbdy,nbdy, halo_us)
      call xctilr(ujb,    1,   1, nbdy,nbdy, halo_us)
      call xctilr(uflx,   1,2*kk, nbdy,nbdy, halo_us)
      call xctilr(utflx,  1,2*kk, nbdy,nbdy, halo_us)
      call xctilr(usflx,  1,2*kk, nbdy,nbdy, halo_us)
      call xctilr(umfltd, 1,2*kk, nbdy,nbdy, halo_us)
      call xctilr(u   ,   1,2*kk, nbdy,nbdy, halo_us)
#ifdef TRC
c$OMP PARALLEL DO PRIVATE(i,nt)
      do j=1-nbdy,jdm+nbdy
        do i=1-nbdy,idm+nbdy
          do nt=1,ntr
            uflxtr(nt,i,j)=uflux(i,j)
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
#endif
c
c --- initialize  v,vb,vtotm,vflx,vflux,vflux2/3,via,vib  at points
c --- located upstream and downstream (in j direction) of p points.
c --- initialize  pby,dpv,vtotn  upstream and downstream of p
c --- points as well as at lateral neighbors of interior v points.
c
c$OMP PARALLEL DO PRIVATE(l,j,k)
      do i=0,ii+1
        do l=1,jsv(i)
        do j=max(1,jfv(i,l)),min(jj,jlv(i,l))
          pv(i,j,1)=0.
          vtotn (i-1,j)=0.
          vtotn (i+1,j)=0.
          do k=1,2
            pbv(i-1,j,k)=0.
            pbv(i+1,j,k)=0.
          enddo
          pbv_p(i-1,j)=0.
          pbv_p(i+1,j)=0.
          do k=1,kk
            dpv(i-1,j,k   )=0.
            dpv(i-1,j,k+kk)=0.
            dpvold(i-1,j,k)=0.
            dpv(i+1,j,k   )=0.
            dpv(i+1,j,k+kk)=0.
            dpvold(i+1,j,k)=0.
          enddo
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,j,k)
      do i=1,ii
        do l=1,jsp(i)
        do j=max(1,jfp(i,l)),min(jj,jlp(i,l)+1)
          do k=1,3
            vbflxs(i,j,k)=0.
          enddo
          do k=1,2
            vbflx(i,j,k)=0.
            vbflx_mn(i,j,k)=0.
            vb(i,j,k)=0.
            vbflxs_p(i,j,k)=0.
            pbv(i,j,k)=0.
          enddo
          pbv_p(i,j)=0.
          vbcors_p(i,j)=0.
          vtotn (i,j)=0.
          vtotm(i,j)=0.
          vflux (i,j)=0.
          vflux2(i,j)=0.
          vflux3(i,j)=0.
          via(i,j)=0.
          vib(i,j)=0.
          do k=1,kk
            dpv(i,j,k   )=0.
            dpv(i,j,k+kk)=0.
            dpvold(i,j,k)=0.
            vflx(i,j,k   )=0.
            vflx(i,j,k+kk)=0.
            vtflx(i,j,k   )=0.
            vtflx(i,j,k+kk)=0.
            vsflx(i,j,k   )=0.
            vsflx(i,j,k+kk)=0.
            vmfltd(i,j,k   )=0.
            vmfltd(i,j,k+kk)=0.
            v(i,j,k   )=0.
            v(i,j,k+kk)=0.
          enddo
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call xctilr(pv,     1,   1, nbdy,nbdy, halo_vs)
      call xctilr(vtotn,  1,   1, nbdy,nbdy, halo_vs)
      call xctilr(pbv,    1,   2, nbdy,nbdy, halo_vs)
      call xctilr(pbv_p,  1,   1, nbdy,nbdy, halo_vs)
      call xctilr(dpv,    1,2*kk, nbdy,nbdy, halo_vs)
      call xctilr(dpvold, 1,  kk, nbdy,nbdy, halo_vs)
      call xctilr(vb,     1,   2, nbdy,nbdy, halo_vs)
      call xctilr(vbflx,  1,   2, nbdy,nbdy, halo_vs)
      call xctilr(vbflx_mn, 1,   2, nbdy,nbdy, halo_vs)
      call xctilr(vbflxs, 1,   3, nbdy,nbdy, halo_vs)
      call xctilr(vbflxs_p, 1, 2, nbdy,nbdy, halo_vs)
      call xctilr(vbcors_p, 1, 1, nbdy,nbdy, halo_vs)
      call xctilr(vtotm,  1,   1, nbdy,nbdy, halo_vs)
      call xctilr(vflux,  1,   1, nbdy,nbdy, halo_vs)
      call xctilr(vflux2, 1,   1, nbdy,nbdy, halo_vs)
      call xctilr(vflux3, 1,   1, nbdy,nbdy, halo_vs)
      call xctilr(via,    1,   1, nbdy,nbdy, halo_vs)
      call xctilr(vib,    1,   1, nbdy,nbdy, halo_vs)
      call xctilr(vflx,   1,2*kk, nbdy,nbdy, halo_vs)
      call xctilr(vtflx,  1,2*kk, nbdy,nbdy, halo_vs)
      call xctilr(vsflx,  1,2*kk, nbdy,nbdy, halo_vs)
      call xctilr(vmfltd, 1,2*kk, nbdy,nbdy, halo_vs)
      call xctilr(v   ,   1,2*kk, nbdy,nbdy, halo_vs)
#ifdef TRC
c$OMP PARALLEL DO PRIVATE(i,nt)
      do j=1-nbdy,jdm+nbdy
        do i=1-nbdy,idm+nbdy
          do nt=1,ntr
            vflxtr(nt,i,j)=vflux(i,j)
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
#endif
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'inivar:'
        endif
c       call chksummsk(p,ip,kk+1,'p')
c       call chksummsk(pu,iu,kk+1,'pu')
c       call chksummsk(pv,iv,kk+1,'pv')
c       call chksummsk(pb,ip,2,'pb')
c       call chksummsk(pb_mn,ip,2,'pb')
c       call chksummsk(pb_p,ip,1,'pb_p')
c       call chksummsk(ub,iu,2,'ub')
c       call chksummsk(vb,iv,2,'vb')
c       call chksummsk(pbu,iu,2,'pbu')
c       call chksummsk(pbu_p,iu,1,'pbu_p')
c       call chksummsk(pbv,iv,2,'pbv')
c       call chksummsk(pbv_p,iv,1,'pbv_p')
c       call chksummsk(ubflx,iu,2,'ubflx')
c       call chksummsk(ubflx_mn,iu,2,'ubflx')
c       call chksummsk(vbflx,iv,2,'vbflx')
c       call chksummsk(vbflx_mn,iv,2,'vbflx')
c       call chksummsk(ubflxs,iu,3,'ubflxs')
c       call chksummsk(ubflxs_p,iu,2,'ubflxs_p')
c       call chksummsk(vbflxs,iv,3,'vbflxs')
c       call chksummsk(vbflxs_p,iv,2,'vbflxs_p')
c       call chksummsk(ubcors_p,iu,1,'ubcors_p')
c       call chksummsk(vbcors_p,iv,1,'vbcors_p')
c       call chksummsk(utotm,iu,1,'utotm')
c       call chksummsk(vtotm,iv,1,'vtotm')
c       call chksummsk(utotn,iu,1,'utotn')
c       call chksummsk(vtotn,iv,1,'vtotn')
c       call chksummsk(uflux,iu,1,'uflux')
c       call chksummsk(vflux,iv,1,'vflux')
c       call chksummsk(uflux2,iu,1,'uflux2')
c       call chksummsk(vflux2,iv,1,'vflux2')
c       call chksummsk(uflux3,iu,1,'uflux3')
c       call chksummsk(vflux3,iv,1,'vflux3')
c       call chksummsk(uja,iu,1,'uja')
c       call chksummsk(ujb,iu,1,'ujb')
c       call chksummsk(via,iv,1,'via')
c       call chksummsk(vib,iv,1,'vib')
c       call chksummsk(u,iu,2*kk,'u')
c       call chksummsk(v,iv,2*kk,'v')
c       call chksummsk(uflx,iu,2*kk,'uflx')
c       call chksummsk(vflx,iv,2*kk,'vflx')
c       call chksummsk(dp,ip,2*kk,'dp')
c       call chksummsk(dpu,iu,2*kk,'dpu')
c       call chksummsk(dpuold,iu,kk,'dpuold')
c       call chksummsk(dpv,iv,2*kk,'dpv')
c       call chksummsk(dpvold,iv,kk,'dpvold')
      endif
c
      return
      end
