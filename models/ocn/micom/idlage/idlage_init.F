      subroutine idlage_init
c
c --- ------------------------------------------------------------------
c --- initialization of ideal age tracer
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "param_trc.h"
#include "common_trc.h"
c
      integer i,j,k,l
c
c$OMP PARALLEL DO
      do j=1-nbdy,jj+nbdy
        do k=1,kk
          do l=1,isp(j)
          do i=ifp(j,l),ilp(j,l)
            trc(i,j,k,itriag)=0.
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
