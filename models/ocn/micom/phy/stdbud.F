      subroutine stdbud(ncall,n,nn)
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "param_bud.h"
#include "common_bud.h"
c
      integer ncall,n,nn
c
      integer i,j,k,l,kn
      real q
c
#ifdef TKE
#  include "common_tke.h"
#endif
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
#endif
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=0.
          util2(i,j)=0.
#ifdef TKE
          util3(i,j)=0.
#  ifdef GLS
          util4(i,j)=0.
#  endif
#endif
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(kn,q)
      do j=1,jj
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            q=dp(i,j,kn)*scp2(i,j)
            util1(i,j)=util1(i,j)+saln(i,j,kn)*q
            util2(i,j)=util2(i,j)+temp(i,j,kn)*q
#ifdef TKE
            util3(i,j)=util3(i,j)+tke(i,j,kn)*q
#  ifdef GLS
            util4(i,j)=util4(i,j)+gls_psi(i,j,kn)*q
#  endif
#endif
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call xcsum(sdp(ncall,n),util1,ips)
      call xcsum(tdp(ncall,n),util2,ips)
#ifdef TKE
      call xcsum(tkedp(ncall,n),util3,ips)
#  ifdef GLS
      call xcsum(glsdp(ncall,n),util4,ips)
#  endif
#endif
c
#ifdef TRC
      if (n.eq.1) then
c
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j)=0.
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(kn,q)
        do j=1,jj
          do k=1,kk
            kn=k+nn
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              q=dp(i,j,kn)*scp2(i,j)
              util1(i,j)=util1(i,j)+trc(i,j,k,1)*q
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
        call xcsum(trdp(ncall),util1,ips)
c
      endif
#endif
c
      return
      end
