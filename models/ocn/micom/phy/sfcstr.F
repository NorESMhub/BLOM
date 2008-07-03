      subroutine sfcstr(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- compute the surface stress
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
#include "common_ice.h"
c
      integer i,j,l
      real facice
c
      call xctilr(ficem, 1,1, 1,1, halo_ps)
      call xctilr(hicem, 1,1, 1,1, halo_ps)
c
c$OMP PARALLEL DO PRIVATE(facice)
      do 200 j=1,jj
c
      do 210 l=1,isu(j)
      do 210 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
      facice=(ficem(i,j)+ficem(i-1,j))
     .       *min(2.,hicem(i,j)+hicem(i-1,j))*.25
 210  taux(i,j)=10.*(ztx(i,j)*(1.-facice)+tauxice(i,j)*facice)
c
      do 200 l=1,isv(j)
      do 200 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
      facice=(ficem(i,j)+ficem(i,j-1))
     .       *min(2.,hicem(i,j)+hicem(i,j-1))*.25
 200  tauy(i,j)=10.*(mty(i,j)*(1.-facice)+tauyice(i,j)*facice)
c$OMP END PARALLEL DO
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'sfcstr:'
        endif
        call chksummsk(taux,iu,1,'taux')
        call chksummsk(tauy,iv,1,'tauy')
      endif
c
      return
      end
