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
      real sdpj(jdm),tdpj(jdm),ddpj(jdm),volj(jdm)
      real vol
c
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      real trdpj(jdm)
#endif
c
      real sig
      external sig
c
      do j=1,jj
        sdpj(j)=0.
        tdpj(j)=0.
        ddpj(j)=0.
        volj(j)=0.
      enddo
c$OMP PARALLEL DO PRIVATE(kn)
      do j=1,jj
        do l=1,isp(j)
        do i=ifp(j,l),ilp(j,l)
          sdpj(j)=sdpj(j)+saln(i,j,k1n)*dp(i,j,k1n)*scp2(i,j)
          tdpj(j)=tdpj(j)+temp(i,j,k1n)*dp(i,j,k1n)*scp2(i,j)
        enddo
        enddo
        do k=2,kk
          kn=k+nn
          do l=1,isp(j)
          do i=ifp(j,l),ilp(j,l)
            sdpj(j)=sdpj(j)+saln(i,j,kn)*dp(i,j,kn)*scp2(i,j)
            tdpj(j)=tdpj(j)+temp(i,j,kn)*dp(i,j,kn)*scp2(i,j)
            ddpj(j)=ddpj(j)
     .             +(sig(temp(i,j,kn),saln(i,j,kn))-sigmar(k))
     .              *dp(i,j,kn)*scp2(i,j)
            volj(j)=volj(j)+dp(i,j,kn)*scp2(i,j)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      sdp(ncall,n)=0.
      tdp(ncall,n)=0.
      ddp(ncall,n)=0.
      vol=0.
      do j=1,jj
        sdp(ncall,n)=sdp(ncall,n)+sdpj(j)
        tdp(ncall,n)=tdp(ncall,n)+tdpj(j)
        ddp(ncall,n)=ddp(ncall,n)+ddpj(j)
        vol=vol+volj(j)
      enddo
      ddp(ncall,n)=ddp(ncall,n)/vol
c
#ifdef TRC
      if (n.eq.1) then
        do j=1,jj
          trdpj(j)=0.
        enddo
c$OMP PARALLEL DO PRIVATE(kn)
        do j=1,jj
          do k=1,kk
            kn=k+nn
            do l=1,isp(j)
            do i=ifp(j,l),ilp(j,l)
              trdpj(j)=trdpj(j)+trc(i,j,k,ntr)*dp(i,j,kn)*scp2(i,j)
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
        trdp(ncall)=0.
        do j=1,jj
          trdp(ncall)=trdp(ncall)+trdpj(j)
        enddo
      endif
#endif
c
      return
      end
