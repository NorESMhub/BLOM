      subroutine stdbud(ncall,n,nn,k1n)
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "param_bud.h"
#include "common_bud.h"
c
      integer ncall,n,nn,k1n
c
      integer i,j,k,l,kn
      real q,vol
c
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
#endif
c
      real sig
      external sig
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=0.
          util2(i,j)=0.
          util3(i,j)=0.
          util4(i,j)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(kn,q)
#ifdef ARCTIC
      do j=1,min(jj,jtdm-j0-1)
#else
      do j=1,jj
#endif
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          q=dp(i,j,k1n)*scp2(i,j)
          util1(i,j)=util1(i,j)+saln(i,j,k1n)*q
          util2(i,j)=util2(i,j)+temp(i,j,k1n)*q
        enddo
        enddo
        do k=2,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            q=dp(i,j,kn)*scp2(i,j)
            util1(i,j)=util1(i,j)+saln(i,j,kn)*q
            util2(i,j)=util2(i,j)+temp(i,j,kn)*q
            util3(i,j)=util3(i,j)
     .        +(sig(temp(i,j,kn),saln(i,j,kn))-sigmar(i,j,k))*q
            util4(i,j)=util4(i,j)+q
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      call xcsum(sdp(ncall,n),util1,ip)
      call xcsum(tdp(ncall,n),util2,ip)
      call xcsum(ddp(ncall,n),util3,ip)
      call xcsum(vol         ,util4,ip)
      ddp(ncall,n)=ddp(ncall,n)/vol
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
#ifdef ARCTIC
        do j=1,min(jj,jtdm-j0-1)
#else
        do j=1,jj
#endif
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
        call xcsum(trdp(ncall),util1,ip)
c
      endif
#endif
c
      return
      end
