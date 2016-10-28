      subroutine icaux_cesm
c
c --- ------------------------------------------------------------------
c --- Define initial conditions specifically when coupled to CESM.
c --- ------------------------------------------------------------------
c
c
      use mod_xc
c
      implicit none
c
#include "common_forc.h"
c
      integer i,j,l
c
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          frzpot(i,j)=0.
          mltpot(i,j)=0.
          flxco2(i,j)=0.
          flxdms(i,j)=0.
        enddo
        enddo
      enddo
c
      return
      end
