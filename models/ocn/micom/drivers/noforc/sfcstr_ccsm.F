      subroutine sfcstr_ccsm(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Compute the surface stress. To be used when coupled to CCSM
c --- ------------------------------------------------------------------
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
      integer i,j,l
c
c$OMP PARALLEL DO
      do 200 j=1,jj
c
      do 210 l=1,isu(j)
      do 210 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
 210  taux(i,j)=0.
c
      do 200 l=1,isv(j)
      do 200 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
 200  tauy(i,j)=0.
c$OMP END PARALLEL DO
c
      return
      end
