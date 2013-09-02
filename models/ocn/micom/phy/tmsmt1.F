      subroutine tmsmt1(m,n,mm,nn,k1m,k1n)
c
c --- save old layer thickness, temperature and salinity for time
c --- smoothing
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
#ifdef TKE
#  include "common_tke.h"
#endif
c
      integer i,j,k,l,kn
c
c$OMP PARALLEL DO PRIVATE(kn)
      do j=1,jj
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            dpold(i,j,kn)=dp(i,j,kn)
            told(i,j,k)=temp(i,j,kn)
            sold(i,j,k)=saln(i,j,kn)
#ifdef TKE
            tke_old(i,j,k)=tke(i,j,kn)
#  ifdef GLS
            gls_psi_old(i,j,k)=gls_psi(i,j,kn)
#  endif
#endif
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'tmsmt1:'
        endif
        call chksummsk(dpold,ip,2*kk,'dpold')
        call chksummsk(told,ip,kk,'told')
        call chksummsk(sold,ip,kk,'sold')
#ifdef TKE
        call chksummsk(tke_old,ip,kk,'tke_old')
#  ifdef GLS
        call chksummsk(gls_psi_old,ip,kk,'gls_psi_old')
#  endif
#endif
      endif
c
      return
      end
