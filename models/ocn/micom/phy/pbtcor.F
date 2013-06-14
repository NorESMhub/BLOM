      subroutine pbtcor(m,n,mm,nn,k1m,k1n)
c
c --- correct the layer thicknesses to better match the bottom pressure
c --- from the barotropic solution
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
c
      real, dimension(1-nbdy:idm+nbdy) :: pbfac
      real dpo,dpni
      integer i,j,k,l,kn,km
c
#ifdef TKE
#  include "common_tke.h"
#endif
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      integer nt
#endif
c
      real sig
      external sig
c
      call xctilr(ubflxs(1-nbdy,1-nbdy,n), 1,1, 1,1, halo_uv)
      call xctilr(vbflxs(1-nbdy,1-nbdy,n), 1,1, 1,1, halo_vv)
#ifdef TKE
      call xctilr(tke(1-nbdy,1-nbdy,k1m), 1,kk, 1,1, halo_ps)
#  ifdef GLS
      call xctilr(gls_psi(1-nbdy,1-nbdy,k1m), 1,kk, 1,1, halo_ps)
#  endif
#endif
#ifdef TRC
      do nt=1,ntr-ntrc
        call xctilr(trc(1-nbdy,1-nbdy,1,nt), 1,kk, 1,1, halo_ps)
      enddo
#endif
c
c$OMP PARALLEL DO PRIVATE(kn)
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
          utotn(i,j)=dlt*ubflxs(i,j,n)
        enddo
        enddo
        do k=1,kk
          kn=k+nn
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            utotn(i,j)=utotn(i,j)-uflx(i,j,kn)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(kn)
      do j=1,jj+1
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          vtotn(i,j)=dlt*vbflxs(i,j,n)
        enddo
        enddo
        do k=1,kk
          kn=k+nn
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            vtotn(i,j)=vtotn(i,j)-vflx(i,j,kn)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(km)
      do j=0,jj+1
        do k=1,kk
          km=k+mm
          do l=1,isp(j)
          do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
            dp(i,j,km)=max(0.,dp(i,j,km))+epsil
            p(i,j,k+1)=p(i,j,k)+dp(i,j,km)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      do k=1,kk
        kn=k+nn
        km=k+mm
c
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isu(j)
          do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            if (utotn(i,j).gt.0.) then
              uflux(i,j)=utotn(i,j)*dp(i-1,j,km)/p(i-1,j,kk+1)
              uflux2(i,j)=uflux(i,j)*saln(i-1,j,km)
              uflux3(i,j)=uflux(i,j)*temp(i-1,j,km)
#ifdef TKE
              uflxtke(i,j)=uflux(i,j)*tke(i-1,j,km)
#  ifdef GLS
              uflxgls_psi(i,j)=uflux(i,j)*gls_psi(i-1,j,km)
#  endif
#endif
#ifdef TRC
              do nt=1,ntr-ntrc
                uflxtr(nt,i,j)=uflux(i,j)*trc(i-1,j,k,nt)
              enddo
#endif
            else
              uflux(i,j)=utotn(i,j)*dp(i  ,j,km)/p(i  ,j,kk+1)
              uflux2(i,j)=uflux(i,j)*saln(i  ,j,km)
              uflux3(i,j)=uflux(i,j)*temp(i  ,j,km)
#ifdef TKE
              uflxtke(i,j)=uflux(i,j)*tke(i  ,j,km)
#  ifdef GLS
              uflxgls_psi(i,j)=uflux(i,j)*gls_psi(i  ,j,km)
#  endif
#endif
#ifdef TRC
              do nt=1,ntr-ntrc
                uflxtr(nt,i,j)=uflux(i,j)*trc(i  ,j,k,nt)
              enddo
#endif
            endif
            uflx(i,j,kn)=uflx(i,j,kn)+uflux(i,j)
            usflx(i,j,kn)=usflx(i,j,kn)+uflux2(i,j)
            utflx(i,j,kn)=utflx(i,j,kn)+uflux3(i,j)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
        do j=1,jj+1
          do l=1,isv(j)
          do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
            if (vtotn(i,j).gt.0.) then
              vflux(i,j)=vtotn(i,j)*dp(i,j-1,km)/p(i,j-1,kk+1)
              vflux2(i,j)=vflux(i,j)*saln(i,j-1,km)
              vflux3(i,j)=vflux(i,j)*temp(i,j-1,km)
#ifdef TKE
              vflxtke(i,j)=vflux(i,j)*tke(i,j-1,km)                 
#  ifdef GLS
              vflxgls_psi(i,j)=vflux(i,j)*gls_psi(i,j-1,km)
#  endif
#endif
#ifdef TRC
              do nt=1,ntr-ntrc
                vflxtr(nt,i,j)=vflux(i,j)*trc(i,j-1,k,nt)
              enddo
#endif
            else
              vflux(i,j)=vtotn(i,j)*dp(i,j  ,km)/p(i,j  ,kk+1)
              vflux2(i,j)=vflux(i,j)*saln(i,j  ,km)
              vflux3(i,j)=vflux(i,j)*temp(i,j  ,km)
#ifdef TKE
              vflxtke(i,j)=vflux(i,j)*tke(i,j  ,km)                 
#  ifdef GLS
              vflxgls_psi(i,j)=vflux(i,j)*gls_psi(i,j  ,km)
#  endif
#endif
#ifdef TRC
              do nt=1,ntr-ntrc
                vflxtr(nt,i,j)=vflux(i,j)*trc(i,j  ,k,nt)
              enddo
#endif
            endif
            vflx(i,j,kn)=vflx(i,j,kn)+vflux(i,j)
            vsflx(i,j,kn)=vsflx(i,j,kn)+vflux2(i,j)
            vtflx(i,j,kn)=vtflx(i,j,kn)+vflux3(i,j)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(dpo,dpni)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            dpo=dp(i,j,km)
            dp(i,j,km)=dpo-scp2i(i,j)*(uflux(i+1,j)-uflux(i,j)
     .                                +vflux(i,j+1)-vflux(i,j))
            dpni=1./dp(i,j,km)
            saln(i,j,km)=(dpo*saln(i,j,km)
     .                   -scp2i(i,j)*(uflux2(i+1,j)-uflux2(i,j)
     .                               +vflux2(i,j+1)-vflux2(i,j)))*dpni
            temp(i,j,km)=(dpo*temp(i,j,km)
     .                   -scp2i(i,j)*(uflux3(i+1,j)-uflux3(i,j)
     .                               +vflux3(i,j+1)-vflux3(i,j)))*dpni
#ifdef TKE
            tke(i,j,km)=(dpo*tke(i,j,km)
     .                  -scp2i(i,j)*(uflxtke(i+1,j)-uflxtke(i,j)
     .                              +vflxtke(i,j+1)-vflxtke(i,j)))*dpni
#  ifdef GLS
            gls_psi(i,j,km)=(dpo*gls_psi(i,j,km)
     .        -scp2i(i,j)*(uflxgls_psi(i+1,j)-uflxgls_psi(i,j)
     .                    +vflxgls_psi(i,j+1)-vflxgls_psi(i,j)))*dpni
#  endif
#endif
#ifdef TRC
            do nt=1,ntr-ntrc
              trc(i,j,k,nt)=(dpo*trc(i,j,k,nt)
     .          -(uflxtr(nt,i+1,j)-uflxtr(nt,i,j)
     .           +vflxtr(nt,i,j+1)-vflxtr(nt,i,j))*scp2i(i,j))*dpni
            enddo
#endif
            sigma(i,j,km)=sig(temp(i,j,km),saln(i,j,km))
            dp(i,j,km)=dp(i,j,km)-epsil
            if (dp(i,j,km).lt.1.e-6) dp(i,j,km)=0.
            p(i,j,k+1)=p(i,j,k)+dp(i,j,km)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      enddo
c
c$OMP PARALLEL DO PRIVATE(km,pbfac)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          pbfac(i)=pb(i,j,m)/p(i,j,kk+1)
        enddo
        enddo
        do k=1,kk
          km=k+mm
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            dp(i,j,km)=dp(i,j,km)*pbfac(i)
            p(i,j,k+1)=p(i,j,k)+dp(i,j,km)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'pbtcor:'
        endif
        call chksummsk(dp,ip,2*kk,'dp')
        call chksummsk(temp,ip,2*kk,'temp')
        call chksummsk(saln,ip,2*kk,'saln')
        call chksummsk(p,ip,kk+1,'p')
        call chksummsk(sigma,ip,2*kk,'sigma')
        call chksummsk(uflx,iu,2*kk,'uflx')
        call chksummsk(vflx,iv,2*kk,'vflx')
#ifdef TKE
        call chksummsk(tke,ip,2*kk,'tke')
#  ifdef GLS
        call chksummsk(gls_psi,ip,2*kk,'gls_psi')
#  endif
#endif
#ifdef TRC
        do nt=1,ntr-ntrc
          call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,kk,'trc')
        enddo
#endif
      endif
c
      return
      end
