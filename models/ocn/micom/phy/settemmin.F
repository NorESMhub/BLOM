      subroutine settemmin
c
c --- ------------------------------------------------------------------
c --- Set minimum physical temperature values in isopycnic layers
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "common_eos.h"
c
      integer i,j,k,l
      real gam,salfrz,a,b,c
c
#ifdef CCSMCOUPLED
c
c --- Set temmin to a constant freezing temperature for all layers
      do k=2,kk
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            temmin(i,j,k)=-1.8
          enddo
          enddo
        enddo
      enddo
#else
c --- Let temmin be the freezing temperature of a given potential
c --- density. This can be achieved by using potential density given in
c --- the function sig and the salinity dependent freezing temperature
c --- given in the function swtfrz.
c
      gam=-.0547
c
      do k=2,kk
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            a=((ap14-ap24*sigmar(i,j,k))*gam
     .        + ap15-ap25*sigmar(i,j,k) )*gam
     .       +ap16-ap26*sigmar(i,j,k)
            b=(ap12-ap22*sigmar(i,j,k))*gam+ap13-ap23*sigmar(i,j,k)
            c=ap11-ap21*sigmar(i,j,k)
            salfrz=(-b+sqrt(b*b-4.*a*c))/(2.*a)
            temmin(i,j,k)=gam*salfrz
          enddo
          enddo
        enddo
      enddo
#endif
c
      if (mnproc.eq.ptest) then
        write (lp,'(a/(6(i5,f8.3)))') 'minimum temperature values:',
     .    (k,temmin(itest,jtest,k),k=2,kk)
      endif
c
      return
      end
