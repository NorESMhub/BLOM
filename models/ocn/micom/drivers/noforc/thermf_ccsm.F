      subroutine thermf_ccsm(m,n,mm,nn,k1m,k1n)
c
c --- NERSC version of thermf. To be used when coupled to CCSM
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
#include "common_forc.h"
c
      integer i,j,k,l
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      integer nt
#endif
c
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
          rnfflx(i,j)=0.
          rfiflx(i,j)=0.
          brnflx(i,j)=0.
          frzpot(i,j)=0.
          mltpot(i,j)=0.
          hmltfz(i,j)=0.
          surflx(i,j)=0.
          sswflx(i,j)=0.
          surrlx(i,j)=0.
          salrlx(i,j)=0.
          salflx(i,j)=0.
#ifdef TRC
          do nt=1,ntrc
            trflx(nt,i,j)=0.
          enddo
#endif
          ustar(i,j)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
