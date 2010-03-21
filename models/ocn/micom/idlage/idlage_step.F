      subroutine idlage_step
c
c --- ------------------------------------------------------------------
c --- update ideal age tracer
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "param_trc.h"
#include "common_trc.h"
#include "common_clndr.h"
c
      integer i,j,k,l
      real q
c
c --- set ideal age to zero in the surface layer and increment the age
c --- in the subsurface layers
c$OMP PARALLEL DO
      do j=1-nbdy,jj+nbdy
        do l=1,isp(j)
        do i=ifp(j,l),ilp(j,l)
          trc(i,j,1,itriag)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      q=delt1/(86400.*nday_in_year)
c$OMP PARALLEL DO
      do j=1-nbdy,jj+nbdy
        do k=2,kk
          do l=1,isp(j)
          do i=ifp(j,l),ilp(j,l)
            trc(i,j,k,itriag)=trc(i,j,k,itriag)+q
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
