      subroutine ocntrc_init
c
c --- ------------------------------------------------------------------
c --- initialization of ocean tracers
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "param_trc.h"
#include "common_trc.h"
#include "common_geo.h"
c
      integer i,j,k,l,nt
#ifdef ATRC
      integer nat
#endif
c
c --- ------------------------------------------------------------------
c --- if no ocean tracers are defined, return
c --- ------------------------------------------------------------------
c
      if (ntrocn.eq.0) return
c
#ifdef ATRC
c --- ------------------------------------------------------------------
c --- initialization of tracers
c --- ------------------------------------------------------------------
c
      do nt=1,ntrocn-natr
        do k=1,kk
          do j=1,jj
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              trc(i,j,k,nt)=(mod(k,5)+1)*(plat(i,j)+90.)/(5.*180.)+nt
              trc(i,j,k+kk,nt)=trc(i,j,k,nt)
            enddo
            enddo
          enddo
        enddo
      enddo
c
c --- ------------------------------------------------------------------
c --- initialization of age tracers
c --- ------------------------------------------------------------------
c
      do nt=1,natr
        nat=ntr-natr+nt
        do k=1,kk
          do j=1,jj
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              trc(i,j,k,nat)=(mod(k,5)*(plon(i,j)+180.)/(4.*360.)+nat)
     .                       *trc(i,j,k,nt)
              trc(i,j,k+kk,nat)=trc(i,j,k,nat)
            enddo
            enddo
          enddo
        enddo
      enddo
c
#else
c --- ------------------------------------------------------------------
c --- initialization of tracers
c --- ------------------------------------------------------------------
c
      do nt=1,ntrocn
        do k=1,kk
          do j=1,jj
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              trc(i,j,k,nt)=(mod(k,5)+1)*(plat(i,j)+90.)/(5.*180.)
              trc(i,j,k+kk,nt)=trc(i,j,k,nt)
            enddo
            enddo
          enddo
        enddo
      enddo
#endif
c
      return
      end
