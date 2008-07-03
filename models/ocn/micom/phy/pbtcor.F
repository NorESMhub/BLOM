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
      integer i,j,k,l,kn
c
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      integer nt
#endif
c
      real sig
      external sig
c
      call xctilr(ubflxs(1-nbdy,1-nbdy,n), 1,1, 2,2, halo_uv)
      call xctilr(vbflxs(1-nbdy,1-nbdy,n), 1,1, 2,2, halo_vv)
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
          utotm(i,j)=dlt*(ubflxs(i,j,m)+ubflxs(i,j,n))
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
          vtotm(i,j)=dlt*(vbflxs(i,j,m)+vbflxs(i,j,n))
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
c$OMP PARALLEL DO PRIVATE(kn)
      do j=0,jj+1
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
            dp(i,j,kn)=max(0.,dp(i,j,kn))+epsil
            p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
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
            dp(i,j,kn)=dpo-scp2i(i,j)*(uflux(i+1,j)-uflux(i,j)
     .                                +vflux(i,j+1)-vflux(i,j))
            dpni=1./dp(i,j,kn)
            saln(i,j,kn)=(dpo*saln(i,j,kn)
     .                   -scp2i(i,j)*(uflux2(i+1,j)-uflux2(i,j)
     .                               +vflux2(i,j+1)-vflux2(i,j)))*dpni
            temp(i,j,kn)=(dpo*temp(i,j,kn)
     .                   -scp2i(i,j)*(uflux3(i+1,j)-uflux3(i,j)
     .                               +vflux3(i,j+1)-vflux3(i,j)))*dpni
#ifdef TRC
            do nt=1,ntrc
              trc(i,j,k,nt)=(dpo*trc(i,j,k,nt)
     .          -(uflxtr(nt,i+1,j)-uflxtr(nt,i,j)
     .           +vflxtr(nt,i,j+1)-vflxtr(nt,i,j))*scp2i(i,j))*dpni
            enddo
#endif
            sigma(i,j,kn)=sig(temp(i,j,kn),saln(i,j,kn))
            dp(i,j,kn)=dp(i,j,kn)-epsil
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
          pbfac(i)=pb(i,j,n)/p(i,j,kk+1)
        enddo
        enddo
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            dp(i,j,kn)=dp(i,j,kn)*pbfac(i)
            p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
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
