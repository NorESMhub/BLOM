      subroutine advect(m,n,mm,nn,k1m,k1n)
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
c
      integer i,j,k,l,km,kn,iw,ie,js,jn,isw,jsw,ise,jse,inw,jnw,ine,jne
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: pbmin
      real, dimension(1-nbdy:idm+nbdy) :: pbfac
      real flxhi,q,dpo,dpni
      logical isopyc
c
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      integer nt
#endif
c
      real dpeps
      parameter (dpeps=1.e-4)
c
c --- ---------------------------------------------------------------------
c --- estimate diffusive flux by biharmonic thickness diffusion
c --- (literally, interface depth diffusion)
c --- ---------------------------------------------------------------------
c
c$OMP PARALLEL DO
      do j=-1,jj+2
        do l=1,isu(j)
        do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
          uflxdf(i,j,1)=0.
          uflxdf(i,j,2)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
      do j=0,jj+2
        do l=1,isv(j)
        do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
          vflxdf(i,j,1)=0.
          vflxdf(i,j,2)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      do k=3,kk
        kn=k+nn
c
c$OMP PARALLEL DO
        do j=-1,jj+2
          do l=1,isp(j)
          do i=max(-1,ifp(j,l)),min(ii+2,ilp(j,l))
            if ((p(i,j,kk+1)-p(i,j,k)).lt.onem.or.
     .          (p(i,j,k)-p(i,j,2)).lt.onem) then
              util1(i,j)=0.
              util2(i,j)=0.
            else
              util1(i,j)=p(i,j,k)-.5*(p(i-1,j,k)+p(i+1,j,k))
              util2(i,j)=p(i,j,k)-.5*(p(i,j-1,k)+p(i,j+1,k))
              if (util1(i,j).gt.0..and.
     .            ((p(i-1,j,kk+1)-p(i-1,j,k)).lt.onem.or.
     .             (p(i+1,j,kk+1)-p(i+1,j,k)).lt.onem)) util1(i,j)=0.
              if (util2(i,j).gt.0..and.
     .            ((p(i,j-1,kk+1)-p(i,j-1,k)).lt.onem.or.
     .             (p(i,j+1,kk+1)-p(i,j+1,k)).lt.onem)) util2(i,j)=0.
              if (util1(i,j).lt.0..and.
     .            ((p(i-1,j,k)-p(i-1,j,2)).lt.onem.or.
     .             (p(i+1,j,k)-p(i+1,j,2)).lt.onem)) util1(i,j)=0.
              if (util2(i,j).lt.0..and.
     .            ((p(i,j-1,k)-p(i,j-1,2)).lt.onem.or.
     .             (p(i,j+1,k)-p(i,j+1,2)).lt.onem)) util2(i,j)=0.
            endif
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(flxhi,q)
        do j=-1,jj+2
          do l=1,isu(j)
          do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
            flxhi=.125*min(scp2(i-1,j),scp2(i,j))
     .                *min(dpu(i,j,kn-1),dpu(i,j,kn))
            q=min(flxhi,max(-flxhi,
     .            delt1*thkdff*(util1(i-1,j)-util1(i,j))*scuy(i,j)))
            uflxdf(i,j,k-1)=uflxdf(i,j,k-1)+q
            uflxdf(i,j,k  )=-q
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(flxhi,q)
        do j=0,jj+2
          do l=1,isv(j)
          do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            flxhi=.125*min(scp2(i,j-1),scp2(i,j))
     .                *min(dpv(i,j,kn-1),dpv(i,j,kn))
            q=min(flxhi,max(-flxhi,
     .            delt1*thkdff*(util2(i,j-1)-util2(i,j))*scvx(i,j)))
            vflxdf(i,j,k-1)=vflxdf(i,j,k-1)+q
            vflxdf(i,j,k  )=-q
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      enddo
c
c$OMP PARALLEL DO PRIVATE(iw,ie,js,jn,isw,jsw,ise,jse,inw,jnw,ine,jne)
      do j=-1,jj+2
        do l=1,isp(j)
        do i=max(-1,ifp(j,l)),min(ii+2,ilp(j,l))
          iw=i-iu(i  ,j)
          ie=i+iu(i+1,j)
          js=j-iv(i,j  )
          jn=j+iv(i,j+1)
          isw=i*(1-ip(iw,js))+iw*ip(iw,js)
          jsw=j*(1-ip(iw,js))+js*ip(iw,js)
          ise=i*(1-ip(ie,js))+ie*ip(ie,js)
          jse=j*(1-ip(ie,js))+js*ip(ie,js)
          inw=i*(1-ip(iw,jn))+iw*ip(iw,jn)
          jnw=j*(1-ip(iw,jn))+jn*ip(iw,jn)
          ine=i*(1-ip(ie,jn))+ie*ip(ie,jn)
          jne=j*(1-ip(ie,jn))+jn*ip(ie,jn)
          pbmin(i,j)=
     .      min(p(isw,jsw,kk+1),p(i  ,js ,kk+1),p(ise,jse,kk+1),
     .          p(iw ,j  ,kk+1),p(i  ,j  ,kk+1),p(ie ,j  ,kk+1),
     .          p(inw,jnw,kk+1),p(i  ,jn ,kk+1),p(ine,jne,kk+1))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call xctilr(u, 1,2*kk, 2,2, halo_uv)
      call xctilr(v, 1,2*kk, 2,2, halo_vv)
      call xctilr(ubflxs_p, 1,1, 2,2, halo_uv)
      call xctilr(vbflxs_p, 1,1, 2,2, halo_vv)
      call xctilr(temp(1-nbdy,1-nbdy,k1n), 1,kk, 3,3, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1n), 1,kk, 3,3, halo_ps)
#ifdef TRC
      do nt=1,ntrc
        call xctilr(trc(1-nbdy,1-nbdy,1,nt), 1,kk, 3,3, halo_ps)
      enddo
#endif
c
c$OMP PARALLEL DO PRIVATE(km,kn,isopyc)
c$OMP+FIRSTPRIVATE(utotm,vtotm)
      do k=1,kk
        km=k+mm
        kn=k+nn
c
c --- - advective and diffusive velocity at mid time level
c
        do j=-1,jj+2
          do l=1,isu(j)
          do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
            utotm(i,j)=u(i,j,km)
     .        +((ubflxs(i,j,m)+ubflxs_p(i,j))*dlt/pbu(i,j,m)
     .         +uflxdf(i,j,k)/max(onemm,dpu(i,j,kn)))/(delt1*scuy(i,j))
          enddo
          enddo
        enddo
        do j=0,jj+2
          do l=1,isv(j)
          do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            vtotm(i,j)=v(i,j,km)
     .        +((vbflxs(i,j,m)+vbflxs_p(i,j))*dlt/pbv(i,j,m)
     .         +vflxdf(i,j,k)/max(onemm,dpv(i,j,kn)))/(delt1*scvx(i,j))
          enddo
          enddo
        enddo
c
        if (k.eq.1) then
          isopyc=.false.
        else
          isopyc=.true.
        endif
c
        call remap(scuy,scvx,scp2i,scp2,pbmin,
     .             pbu(1-nbdy,1-nbdy,n),pbv(1-nbdy,1-nbdy,n),
     .             p(1-nbdy,1-nbdy,k+1),utotm,vtotm,delt1,isopyc,1,
     .             dp(1-nbdy,1-nbdy,kn),temp(1-nbdy,1-nbdy,kn),
     .             saln(1-nbdy,1-nbdy,kn),uflx(1-nbdy,1-nbdy,k),
     .             vflx(1-nbdy,1-nbdy,k)
#ifdef TRC
     .            ,k,ntrc,trc
#  ifdef ATRC
     .            ,natrc
#  endif
#endif
     .            )
c
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- by utilizing the time averaged barotropic mass flux, correct the
c --- layer thicknesses to better match the predictive bottom pressure
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(kn)
      do j=0,jj+1
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
            p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
          utotm(i,j)=dlt*(ubflxs(i,j,m)+ubflxs_p(i,j))
        enddo
        enddo
        do k=1,kk
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            utotm(i,j)=utotm(i,j)-uflx(i,j,k)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
      do j=1,jj+1
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          vtotm(i,j)=dlt*(vbflxs(i,j,m)+vbflxs_p(i,j))
        enddo
        enddo
        do k=1,kk
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            vtotm(i,j)=vtotm(i,j)-vflx(i,j,k)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      do k=1,kk
        kn=k+nn
c
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            if (utotm(i,j).gt.0.) then
              uflux(i,j)=utotm(i,j)*dp(i-1,j,kn)/p(i-1,j,kk+1)
              uflux2(i,j)=uflux(i,j)*saln(i-1,j,kn)
              uflux3(i,j)=uflux(i,j)*temp(i-1,j,kn)
#ifdef TRC
              do nt=1,ntrc
                uflxtr(nt,i,j)=uflux(i,j)*trc(i-1,j,k,nt)
              enddo
#endif
            else
              uflux(i,j)=utotm(i,j)*dp(i  ,j,kn)/p(i  ,j,kk+1)
              uflux2(i,j)=uflux(i,j)*saln(i  ,j,kn)
              uflux3(i,j)=uflux(i,j)*temp(i  ,j,kn)
#ifdef TRC
              do nt=1,ntrc
                uflxtr(nt,i,j)=uflux(i,j)*trc(i  ,j,k,nt)
              enddo
#endif
            endif
            uflx(i,j,k)=uflx(i,j,k)+uflux(i,j)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
        do j=1,jj+1
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            if (vtotm(i,j).gt.0.) then
              vflux(i,j)=vtotm(i,j)*dp(i,j-1,kn)/p(i,j-1,kk+1)
              vflux2(i,j)=vflux(i,j)*saln(i,j-1,kn)
              vflux3(i,j)=vflux(i,j)*temp(i,j-1,kn)
#ifdef TRC
              do nt=1,ntrc
                vflxtr(nt,i,j)=vflux(i,j)*trc(i,j-1,k,nt)
              enddo
#endif
            else
              vflux(i,j)=vtotm(i,j)*dp(i,j  ,kn)/p(i,j  ,kk+1)
              vflux2(i,j)=vflux(i,j)*saln(i,j  ,kn)
              vflux3(i,j)=vflux(i,j)*temp(i,j  ,kn)
#ifdef TRC
              do nt=1,ntrc
                vflxtr(nt,i,j)=vflux(i,j)*trc(i,j  ,k,nt)
              enddo
#endif
            endif
            vflx(i,j,k)=vflx(i,j,k)+vflux(i,j)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(dpo,dpni)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            dpo=dp(i,j,kn)
            dp(i,j,kn)=max(0.,dpo-(uflux(i+1,j)-uflux(i,j)
     .                            +vflux(i,j+1)-vflux(i,j))*scp2i(i,j))
            dpo=dpo+dpeps
            dpni=1./(dp(i,j,kn)+dpeps)
            saln(i,j,kn)=(dpo*saln(i,j,kn)
     .                   -(uflux2(i+1,j)-uflux2(i,j)
     .                    +vflux2(i,j+1)-vflux2(i,j))*scp2i(i,j))*dpni
            temp(i,j,kn)=(dpo*temp(i,j,kn)
     .                   -(uflux3(i+1,j)-uflux3(i,j)
     .                    +vflux3(i,j+1)-vflux3(i,j))*scp2i(i,j))*dpni
#ifdef TRC
            do nt=1,ntrc
              trc(i,j,k,nt)=(dpo*trc(i,j,k,nt)
     .          -(uflxtr(nt,i+1,j)-uflxtr(nt,i,j)
     .           +vflxtr(nt,i,j+1)-vflxtr(nt,i,j))*scp2i(i,j))*dpni
            enddo
#endif
            p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      enddo
c
c$OMP PARALLEL DO PRIVATE(kn,pbfac)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          pbfac(i)=pb_p(i,j)/p(i,j,kk+1)
        enddo
        enddo
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            dp(i,j,kn)=dp(i,j,kn)*pbfac(i)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'advect:'
        endif
        call chksummsk(dp,ip,2*kk,'dp')
        call chksummsk(temp,ip,2*kk,'temp')
        call chksummsk(saln,ip,2*kk,'saln')
        call chksummsk(uflx,iu,kk,'uflx')
        call chksummsk(vflx,iv,kk,'vflx')
#ifdef TRC
        do nt=1,ntrc
          call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,kk,'trc')
        enddo
#endif
      endif
c
      return
      end
