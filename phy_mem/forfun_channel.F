! ------------------------------------------------------------------------------
! Copyright (C) 2020 Aleksi Nummelin, Mats Bentsen
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

      subroutine forfun_channel

c
c --- ------------------------------------------------------------------
c --- Specify forcing functions for channel
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_forc.h"
c
      integer i,j,l,k
c
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          do k=1,12
            sstclm(i,j,k)=5.+15.*j/jj    
            sssclm(i,j,k)=35.+1.*j/jj
          enddo
        enddo
        enddo
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          ztx(i,j)=1.0e-2
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          mty(i,j)=0.
        enddo
        enddo
      enddo
c
c relaxation timescales
      trxday=0.
      srxday=0.
c
      return
      end
