! ------------------------------------------------------------------------------
! Copyright (C) 2009-2018 Mats Bentsen
!
! This file is part of BLOM.
!
! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

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
c$OMP PARALLEL DO PRIVATE(k,kn,l,i)
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
