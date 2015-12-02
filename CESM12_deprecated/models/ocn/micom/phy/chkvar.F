      subroutine chkvar(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- check for NaN of Inf in layer thickness
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
c
      integer i,j,k,l,kn
c
      logical isnan,isinf
      external isnan,isinf
c
c$OMP PARALLEL DO PRIVATE(kn)
      do j=1,jj
        do k=1,kk
          kn=k+nn
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (isnan(dp(i,j,kn))) then
              write (lp,'(a,i4,a,i4,a,i4,a)')
     .          ' chkvar: dp is NaN at (i =',i0+i,', j =',j0+j,
     .          ', k =',k,')'
             call xchalt('(chkvar') 
                    stop '(chkvar)'
            endif
            if (isinf(dp(i,j,kn))) then
              write (lp,'(a,i4,a,i4,a,i4,a)')
     .          ' chkvar: dp is Inf at (i =',i0+i,', j =',j0+j,
     .          ', k =',k,')'
             call xchalt('(chkvar') 
                    stop '(chkvar)'
            endif
            if (isnan(temp(i,j,kn))) then
              write (lp,'(a,i4,a,i4,a,i4,a)')
     .          ' chkvar: temp is NaN at (i =',i0+i,', j =',j0+j,
     .          ', k =',k,')'
             call xchalt('(chkvar') 
                    stop '(chkvar)'
            endif
            if (isinf(temp(i,j,kn))) then
              write (lp,'(a,i4,a,i4,a,i4,a)')
     .          ' chkvar: temp is Inf at (i =',i0+i,', j =',j0+j,
     .          ', k =',k,')'
             call xchalt('(chkvar') 
                    stop '(chkvar)'
            endif
            if (isnan(saln(i,j,kn))) then
              write (lp,'(a,i4,a,i4,a,i4,a)')
     .          ' chkvar: saln is NaN at (i =',i0+i,', j =',j0+j,
     .          ', k =',k,')'
             call xchalt('(chkvar') 
                    stop '(chkvar)'
            endif
            if (isinf(saln(i,j,kn))) then
              write (lp,'(a,i4,a,i4,a,i4,a)')
     .          ' chkvar: saln is Inf at (i =',i0+i,', j =',j0+j,
     .          ', k =',k,')'
             call xchalt('(chkvar') 
                    stop '(chkvar)'
            endif
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
