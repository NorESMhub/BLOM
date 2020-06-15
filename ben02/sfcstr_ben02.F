! ------------------------------------------------------------------------------
! Copyright (C) 2004-2015 Mats Bentsen
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

      subroutine sfcstr_ben02(m,n,mm,nn,k1m,k1n)
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
#include "common_seaice.h"
c
      integer i,j,l
      real facice
c
      call xctilr(ficem, 1,1, 1,1, halo_ps)
      call xctilr(hicem, 1,1, 1,1, halo_ps)
c
c$OMP PARALLEL DO PRIVATE(l,i,facice)
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          facice=(ficem(i,j)+ficem(i-1,j))
     .           *min(2.,hicem(i,j)+hicem(i-1,j))*.25
          taux(i,j)=10.*(ztx(i,j)*(1.-facice)+tauxice(i,j)*facice)
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          facice=(ficem(i,j)+ficem(i,j-1))
     .           *min(2.,hicem(i,j)+hicem(i,j-1))*.25
          tauy(i,j)=10.*(mty(i,j)*(1.-facice)+tauyice(i,j)*facice)
        enddo
        enddo
      enddo
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
