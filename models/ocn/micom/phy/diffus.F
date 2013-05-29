      subroutine diffus(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- diffusion of tracers
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
      integer i,j,k,l,kn,km
      real q
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
      real dpeps
      parameter (dpeps=1.e-4)
c
      call xctilr(dp  (1-nbdy,1-nbdy,k1n), 1,kk, 2,2, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,k1n), 1,kk, 2,2, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1n), 1,kk, 2,2, halo_ps)
#ifdef TRC
      do nt=1,ntrc
        call xctilr(trc(1-nbdy,1-nbdy,1,nt), 1,kk, 2,2, halo_ps)
      enddo
#endif
      call xctilr(difiso, 1,kk, 2,2, halo_ps)
c
      do k=1,kk
        kn=k+nn
        km=k+mm
c
c$OMP PARALLEL DO PRIVATE(q)
        do j=0,jj+1
          do l=1,isu(j)
          do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
            q=delt1*.5*(difiso(i-1,j,k)+difiso(i,j,k))
     .        *scuy(i,j)*scuxi(i,j)
     .        *max(min(dp(i-1,j,kn),dp(i,j,kn)),dpeps)
            usflld(i,j,km)=q*(saln(i-1,j,kn)-saln(i,j,kn))
            utflld(i,j,km)=q*(temp(i-1,j,kn)-temp(i,j,kn))
#ifdef TRC
            do nt=1,ntrc
              uflxtr(nt,i,j)=q*(trc(i-1,j,k,nt)-trc(i,j,k,nt))
            enddo
#endif
            usflx(i,j,km)=usflx(i,j,km)+usflld(i,j,km)
            utflx(i,j,km)=utflx(i,j,km)+utflld(i,j,km)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(q)
        do j=0,jj+2
          do l=1,isv(j)
          do i=max(0,ifv(j,l)),min(ii+1,ilv(j,l))
            q=delt1*.5*(difiso(i,j-1,k)+difiso(i,j,k))
     .        *scvx(i,j)*scvyi(i,j)
     .        *max(min(dp(i,j-1,kn),dp(i,j,kn)),dpeps)
            vsflld(i,j,km)=q*(saln(i,j-1,kn)-saln(i,j,kn))
            vtflld(i,j,km)=q*(temp(i,j-1,kn)-temp(i,j,kn))
#ifdef TRC
            do nt=1,ntrc
              vflxtr(nt,i,j)=q*(trc(i,j-1,k,nt)-trc(i,j,k,nt))
            enddo
#endif
            vsflx(i,j,km)=vsflx(i,j,km)+vsflld(i,j,km)
            vtflx(i,j,km)=vtflx(i,j,km)+vtflld(i,j,km)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(q)
        do j=0,jj+1
          do l=1,isp(j)
          do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
            q=1./(scp2(i,j)*max(dp(i,j,kn),dpeps))
            saln(i,j,kn)=saln(i,j,kn)
     .                  -q*(usflld(i+1,j,km)-usflld(i,j,km)
     .                     +vsflld(i,j+1,km)-vsflld(i,j,km))
            temp(i,j,kn)=temp(i,j,kn)
     .                  -q*(utflld(i+1,j,km)-utflld(i,j,km)
     .                     +vtflld(i,j+1,km)-vtflld(i,j,km))
#ifdef TRC
            do nt=1,ntrc
              trc(i,j,k,nt)=trc(i,j,k,nt)
     .                     -q*(uflxtr(nt,i+1,j)-uflxtr(nt,i,j)
     .                        +vflxtr(nt,i,j+1)-vflxtr(nt,i,j))
            enddo
#endif
            sigma(i,j,kn)=sig(temp(i,j,kn),saln(i,j,kn))
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      enddo
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'diffus:'
        endif
        call chksummsk(temp,ip,2*kk,'temp')
        call chksummsk(saln,ip,2*kk,'saln')
#ifdef TRC
        do nt=1,ntrc
          call chksummsk(trc(1-nbdy,1-nbdy,1,nt),ip,kk,'trc')
        enddo
#endif
      endif
c
      return
      end
