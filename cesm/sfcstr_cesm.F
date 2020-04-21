! ------------------------------------------------------------------------------
! Copyright (C) 2015-2018 Mats Bentsen
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

      subroutine sfcstr_cesm(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Compute the surface stress. To be used when coupled to CESM
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
c$OMP PARALLEL DO PRIVATE(l,i)
      do 200 j=1,jj
c
      do 210 l=1,isu(j)
      do 210 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
 210  taux(i,j)=10.*ztx(i,j)
c
      do 200 l=1,isv(j)
      do 200 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
 200  tauy(i,j)=10.*mty(i,j)
c$OMP END PARALLEL DO
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'sfcstr_cesm:'
        endif
        call chksummsk(taux,iu,1,'taux')
        call chksummsk(tauy,iv,1,'tauy')
      endif
c
      return
      end
