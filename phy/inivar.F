! ------------------------------------------------------------------------------
! Copyright (C) 2015-2020 Mats Bentsen, Jerry Tjiputra
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
          p(i,j,1)=spval
          pu(i,j,1)=spval
          pv(i,j,1)=spval
          phi(i,j,1)=spval
          do k=1,3
            ubflxs(i,j,k)=spval
            vbflxs(i,j,k)=spval
          enddo
          do k=1,2
            pb(i,j,k)=spval
            ubflx(i,j,k)=spval
            vbflx(i,j,k)=spval
            pb_mn(i,j,k)=spval
            ubflx_mn(i,j,k)=spval
            vbflx_mn(i,j,k)=spval
            ub(i,j,k)=spval
            vb(i,j,k)=spval
            ubflxs_p(i,j,k)=spval
            vbflxs_p(i,j,k)=spval
            pbu(i,j,k)=spval
            pbv(i,j,k)=spval
            pvtrop(i,j,k)=spval
          enddo
          pb_p(i,j)=spval
          pbu_p(i,j)=spval
          pbv_p(i,j)=spval
          ubcors_p(i,j)=spval
          vbcors_p(i,j)=spval
          defor1(i,j)=spval
          defor2(i,j)=spval
          pvtrop_o(i,j)=spval
          utotm(i,j)=spval
          vtotm(i,j)=spval
          utotn(i,j)=spval
          vtotn(i,j)=spval
          uflux (i,j)=spval
          vflux (i,j)=spval
          uflux2(i,j)=spval
          vflux2(i,j)=spval
          uflux3(i,j)=spval
          vflux3(i,j)=spval
          uja(i,j)=spval
          ujb(i,j)=spval
          via(i,j)=spval
          vib(i,j)=spval
          pgfxm(i,j,1)=spval
          pgfxm(i,j,2)=spval
          pgfxm_o(i,j)=spval
          xixp(i,j,1)=spval
          xixp(i,j,2)=spval
          xixp_o(i,j)=spval
          xixm(i,j,1)=spval
          xixm(i,j,2)=spval
          xixm_o(i,j)=spval
          pgfym(i,j,1)=spval
          pgfym(i,j,2)=spval
          pgfym_o(i,j)=spval
          xiyp(i,j,1)=spval
          xiyp(i,j,2)=spval
          xiyp_o(i,j)=spval
          xiym(i,j,1)=spval
          xiym(i,j,2)=spval
          xiym_o(i,j)=spval
          difwgt(i,j)=spval
          sealv(i,j)=spval
          twedon(i,j)=spval
          taux(i,j)=spval
          tauy(i,j)=spval
          swa(i,j)=spval
          nsf(i,j)=spval
          hmltfz(i,j)=spval
          hmlt(i,j)=spval
          dfl(i,j)=spval
          lip(i,j)=spval
          sop(i,j)=spval
          eva(i,j)=spval
          rnf(i,j)=spval
          rfi(i,j)=spval
          fmltfz(i,j)=spval
          sfl(i,j)=spval
          ztx(i,j)=spval
          mty(i,j)=spval
          ustarw(i,j)=spval
          slp(i,j)=spval
          abswnd(i,j)=spval
          albw(i,j)=spval
          frzpot(i,j)=spval
          mltpot(i,j)=spval
          atmco2(i,j)=spval
          flxco2(i,j)=spval
          flxdms(i,j)=spval
          alb(i,j)=spval
          rnfres(i,j)=spval
          rnfflx(i,j)=spval
          rfiflx(i,j)=spval
          eiacc(i,j)=spval
          pracc(i,j)=spval
          ficem(i,j)=spval
          hicem(i,j)=spval
          hsnwm(i,j)=spval
          ustari(i,j)=spval
          tauxice(i,j)=spval
          tauyice(i,j)=spval
          uicem(i,j)=spval
          vicem(i,j)=spval
          iagem(i,j)=spval
          do k=1,kk
            u  (i,j,k   )=spval
            u  (i,j,k+kk)=spval
            v  (i,j,k   )=spval
            v  (i,j,k+kk)=spval
            uflx(i,j,k   )=spval
            uflx(i,j,k+kk)=spval
            utflx(i,j,k   )=spval
            utflx(i,j,k+kk)=spval
            usflx(i,j,k   )=spval
            usflx(i,j,k+kk)=spval
            vflx(i,j,k   )=spval
            vflx(i,j,k+kk)=spval
            vtflx(i,j,k   )=spval
            vtflx(i,j,k+kk)=spval
            vsflx(i,j,k   )=spval
            vsflx(i,j,k+kk)=spval
            umfltd(i,j,k   )=spval
            umfltd(i,j,k+kk)=spval
            vmfltd(i,j,k   )=spval
            vmfltd(i,j,k+kk)=spval
            dp (i,j,k   )=spval
            dp (i,j,k+kk)=spval
            dpu(i,j,k   )=spval
            dpu(i,j,k+kk)=spval
            dpuold(i,j,k)=spval
            dpv(i,j,k   )=spval
            dpv(i,j,k+kk)=spval
            dpvold(i,j,k)=spval
            p (i,j,k+1)=spval
            pu(i,j,k+1)=spval
            pv(i,j,k+1)=spval
            phi(i,j,k+1)=spval
c
            temp(i,j,k)=spval
            saln(i,j,k)=spval
            sigma(i,j,k   )=spval
            sigma(i,j,k+kk)=spval
            pgfx(i,j,k   )=spval
            pgfx(i,j,k+kk)=spval
            pgfy(i,j,k   )=spval
            pgfy(i,j,k+kk)=spval
            pgfxo(i,j,k)=spval
            pgfyo(i,j,k)=spval
#ifdef TRC
            do nt=1,ntr
              trc(i,j,k   ,nt)=spval
              trc(i,j,k+kk,nt)=spval
            enddo
#endif
          enddo
          do k=1,4
            uml(i,j,k)=spval
            vml(i,j,k)=spval
          enddo
          do k=1,2
            umlres(i,j,k)=spval
            vmlres(i,j,k)=spval
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
