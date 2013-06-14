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
      real dpo,dpni
c
#if defined(TKE) && defined(TKEADV)
#  include "common_tke.h"
#endif
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      integer nt
#endif
c
      real dpeps
      parameter (dpeps=1.e-4)
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
      call xctilr(ubflxs_p, 1,2, 2,2, halo_uv)
      call xctilr(vbflxs_p, 1,2, 2,2, halo_vv)
      call xctilr(pbu, 1,2, 2,2, halo_us)
      call xctilr(pbv, 1,2, 2,2, halo_vs)
      call xctilr(temp(1-nbdy,1-nbdy,k1n), 1,kk, 3,3, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1n), 1,kk, 3,3, halo_ps)
#if defined(TKE) && defined(TKEADV)
      call xctilr(tke(1-nbdy,1-nbdy,k1n), 1,kk, 3,3, halo_ps)
#  ifdef GLS
      call xctilr(gls_psi(1-nbdy,1-nbdy,k1n), 1,kk, 3,3, halo_ps)
#  endif
#endif
#ifdef TRC
      do nt=1,ntrc
        call xctilr(trc(1-nbdy,1-nbdy,1,nt), 1,kk, 3,3, halo_ps)
      enddo
#endif
c
c$OMP PARALLEL DO PRIVATE(km,kn)
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
     .        +(ubflxs_p(i,j,m)*dlt/pbu(i,j,m)
     .         +umfltd(i,j,km)/max(onemm,dpu(i,j,kn)))/(delt1*scuy(i,j))
            utotm(i,j)=max(-umax(i,j),min(umax(i,j),utotm(i,j)))
          enddo
          enddo
        enddo
        do j=0,jj+2
          do l=1,isv(j)
          do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            vtotm(i,j)=v(i,j,km)
     .        +(vbflxs_p(i,j,m)*dlt/pbv(i,j,m)
     .         +vmfltd(i,j,km)/max(onemm,dpv(i,j,kn)))/(delt1*scvx(i,j))
            vtotm(i,j)=max(-vmax(i,j),min(vmax(i,j),vtotm(i,j)))
          enddo
          enddo
        enddo
c
        call remap(scuy,scvx,scp2i,scp2,pbmin,
     .             pbu(1-nbdy,1-nbdy,n),pbv(1-nbdy,1-nbdy,n),
     .             p(1-nbdy,1-nbdy,k+1),utotm,vtotm,delt1,1,
     .             dp(1-nbdy,1-nbdy,kn),temp(1-nbdy,1-nbdy,kn),
     .             saln(1-nbdy,1-nbdy,kn),uflx(1-nbdy,1-nbdy,km),
     .             vflx(1-nbdy,1-nbdy,km)
     .            ,utflx(1-nbdy,1-nbdy,km),vtflx(1-nbdy,1-nbdy,km)
     .            ,usflx(1-nbdy,1-nbdy,km),vsflx(1-nbdy,1-nbdy,km)
#if defined(TKE) && defined(TKEADV)
     .            ,tke(1-nbdy,1-nbdy,kn)

#  ifdef GLS
     .            ,gls_psi(1-nbdy,1-nbdy,kn)
#  endif
#endif
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
c$OMP PARALLEL DO PRIVATE(km)
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
          utotm(i,j)=dlt*ubflxs_p(i,j,m)
        enddo
        enddo
        do k=1,kk
          km=k+mm
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            utotm(i,j)=utotm(i,j)-uflx(i,j,km)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(km)
      do j=1,jj+1
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          vtotm(i,j)=dlt*vbflxs_p(i,j,m)
        enddo
        enddo
        do k=1,kk
          km=k+mm
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            vtotm(i,j)=vtotm(i,j)-vflx(i,j,km)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      do k=1,kk
        km=k+mm
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
#if defined(TKE) && defined(TKEADV)
              uflxtke(i,j)=uflux(i,j)*tke(i-1,j,kn)
#  ifdef GLS
              uflxgls_psi(i,j)=uflux(i,j)*gls_psi(i-1,j,kn)
#  endif
#endif
#ifdef TRC
              do nt=1,ntrc
                uflxtr(nt,i,j)=uflux(i,j)*trc(i-1,j,k,nt)
              enddo
#endif
            else
              uflux(i,j)=utotm(i,j)*dp(i  ,j,kn)/p(i  ,j,kk+1)
              uflux2(i,j)=uflux(i,j)*saln(i  ,j,kn)
              uflux3(i,j)=uflux(i,j)*temp(i  ,j,kn)
#if defined(TKE) && defined(TKEADV)
              uflxtke(i,j)=uflux(i,j)*tke(i  ,j,kn)
#  ifdef GLS
              uflxgls_psi(i,j)=uflux(i,j)*gls_psi(i  ,j,kn)
#  endif
#endif
#ifdef TRC
              do nt=1,ntrc
                uflxtr(nt,i,j)=uflux(i,j)*trc(i  ,j,k,nt)
              enddo
#endif
            endif
            uflx(i,j,km)=uflx(i,j,km)+uflux(i,j)
            usflx(i,j,km)=usflx(i,j,km)+uflux2(i,j)
            utflx(i,j,km)=utflx(i,j,km)+uflux3(i,j)
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
#if defined(TKE) && defined(TKEADV)
              vflxtke(i,j)=vflux(i,j)*tke(i,j-1,kn)
#  ifdef GLS
              vflxgls_psi(i,j)=vflux(i,j)*gls_psi(i,j-1,kn)
#  endif
#endif
#ifdef TRC
              do nt=1,ntrc
                vflxtr(nt,i,j)=vflux(i,j)*trc(i,j-1,k,nt)
              enddo
#endif
            else
              vflux(i,j)=vtotm(i,j)*dp(i,j  ,kn)/p(i,j  ,kk+1)
              vflux2(i,j)=vflux(i,j)*saln(i,j  ,kn)
              vflux3(i,j)=vflux(i,j)*temp(i,j  ,kn)
#if defined(TKE) && defined(TKEADV)
              vflxtke(i,j)=vflux(i,j)*tke(i,j  ,kn)
#  ifdef GLS
              vflxgls_psi(i,j)=vflux(i,j)*gls_psi(i,j  ,kn)
#  endif
#endif
#ifdef TRC
              do nt=1,ntrc
                vflxtr(nt,i,j)=vflux(i,j)*trc(i,j  ,k,nt)
              enddo
#endif
            endif
            vflx(i,j,km)=vflx(i,j,km)+vflux(i,j)
            vsflx(i,j,km)=vsflx(i,j,km)+vflux2(i,j)
            vtflx(i,j,km)=vtflx(i,j,km)+vflux3(i,j)
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
#if defined(TKE) && defined(TKEADV)
            tke(i,j,kn)=(dpo*tke(i,j,kn)
     .                  -(uflxtke(i+1,j)-uflxtke(i,j)
     .                   +vflxtke(i,j+1)-vflxtke(i,j))*scp2i(i,j))*dpni
#  ifdef GLS
            gls_psi(i,j,kn)=(dpo*gls_psi(i,j,kn)
     .                      -(uflxgls_psi(i+1,j)-uflxgls_psi(i,j)
     .                       +vflxgls_psi(i,j+1)-vflxgls_psi(i,j))
     .                       *scp2i(i,j))*dpni
#  endif
#endif
#ifdef TRC
            do nt=1,ntrc
              trc(i,j,k,nt)=(dpo*trc(i,j,k,nt)
     .          -(uflxtr(nt,i+1,j)-uflxtr(nt,i,j)
     .           +vflxtr(nt,i,j+1)-vflxtr(nt,i,j))*scp2i(i,j))*dpni
            enddo
#endif
            if (dp(i,j,kn).lt.1.e-6) dp(i,j,kn)=0.
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
        call chksummsk(uflx,iu,2*kk,'uflx')
        call chksummsk(vflx,iv,2*kk,'vflx')
#if defined(TKE) && defined(TKEADV)
        call chksummsk(tke,ip,2*kk,'tke')
#  ifdef GLS
        call chksummsk(gls_psi,ip,2*kk,'gls_psi')
#  endif
#endif
#ifdef TRC
        do nt=1,ntrc
          call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,kk,'trc')
        enddo
#endif
      endif
c
      return
      end
