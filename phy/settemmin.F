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
      real salfrz,a,b,c
c
      if     (expcnf.eq.'cesm') then
c
c --- - Set temmin to a constant freezing temperature for all layers
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=2,kk
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              temmin(i,j,k)=-1.8
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      elseif (expcnf.eq.'ben02clim'.or.expcnf.eq.'ben02syn') then
c
c --- - Let temmin be the freezing temperature of a given potential
c --- - density. This can be achieved by using potential density given
c --- - in the function sig and the salinity dependent freezing
c --- - temperature given in the function swtfrz.
c
c$OMP PARALLEL DO PRIVATE(k,l,i,a,b,c,salfrz)
        do j=1,jj
          do k=2,kk
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              a=((ap14-ap24*sigmar(i,j,k))*atf
     .          + ap15-ap25*sigmar(i,j,k) )*atf
     .         +ap16-ap26*sigmar(i,j,k)
              b=(ap12-ap22*sigmar(i,j,k))*atf+ap13-ap23*sigmar(i,j,k)
              c=ap11-ap21*sigmar(i,j,k)
              salfrz=(-b+sqrt(b*b-4.*a*c))/(2.*a)
              temmin(i,j,k)=atf*salfrz
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      elseif (expcnf.eq.'isomip1'.or.expcnf.eq.'isomip2') then
c
c --- - Set temmin to a low value.
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=2,kk
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              temmin(i,j,k)=-10.
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      else
        if (mnproc.eq.1) then
          write (lp,'(3a)') ' expcnf=',trim(expcnf),' is unsupported!'
        endif
        call xcstop('(settemmin)')
               stop '(settemmin)'
      endif
c
      if (mnproc.eq.ptest) then
        write (lp,'(a/(6(i5,f8.3)))') 'minimum temperature values:',
     .    (k,temmin(itest,jtest,k),k=2,kk)
      endif
c
      return
      end
