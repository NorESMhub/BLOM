      subroutine trc_limitc(nn)
c***********************************************************************
c
c**** *SUBROUTINE trc_limitc* - remove negative tracer values.
c
c     J. Schwinger      *GFI, UiB        initial version, 2014-06-17
c      - 
c
c  
c     Purpose
c     -------
c      Remove negative tracer values in the first layer in a mass
c      conservative fashion (i.e. the mass deficit removed is 
c      transfered to non-negative points by a multiplicative 
c      correction). This is done since the virtual tracer fluxes
c      (applied in mxlayr.F directly before HAMOCC is called) can
c      cause negative tracer values in regions with low concentration
c      and strong precipitation.
c
c***********************************************************************
c
      use mod_xc
c
      implicit none
c
#include "param_trc.h"
#include "common_trc.h"
#include "common_blocks.h"
c
      integer       :: nn
      integer       :: i,j,l,nt,kn
      real          :: trbudo(ntrbgc),trbudn,q
c
c --- ------------------------------------------------------------------
c --- - compute tracer budgets before removing negative values
c --- ------------------------------------------------------------------
c
      kn=1+nn
c
      do nt=1,ntrbgc
c
        util1(:,:)=0.
c
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j)=util1(i,j)+trc(i,j,kn,itrbgc+nt-1)*dp(i,j,kn)*scp2(i,j)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
        call xcsum(trbudo(nt),util1,ips)
c
      enddo
c
c
c --- ------------------------------------------------------------------
c --- - remove negative tracer values in the surface layer
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO
      do nt=itrbgc,itrbgc+ntrbgc-1
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              trc(i,j,kn,nt) = max(trc(i,j,kn,nt),0.0)
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c
c --- ------------------------------------------------------------------
c --- - recalculate and correct tracer budgets 
c --- ------------------------------------------------------------------
c
      do nt=1,ntrbgc
c
        util1(:,:)=0.
c
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j)=util1(i,j)+trc(i,j,kn,itrbgc+nt-1)*dp(i,j,kn)*scp2(i,j)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
        call xcsum(trbudn,util1,ips)
        q=trbudo(nt)/max(1.e-14,trbudn)
c
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            trc(i,j,kn,itrbgc+nt-1) = trc(i,j,kn,itrbgc+nt-1)*q
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      enddo
c
      return
      end
