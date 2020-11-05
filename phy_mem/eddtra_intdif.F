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

      subroutine eddtra_intdif(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Estimate eddy-induced transport by interface diffusion.
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
      integer i,j,k,l,km,kn
      real flxhi,flxlo,q
c
      call xctilr(difint, 1,kk, 2,2, halo_ps)
c
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=-1,jj+2
        do l=1,isu(j)
        do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
          umfltd(i,j,1+mm)=0.
          umfltd(i,j,2+mm)=0.
          umfltd(i,j,3+mm)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=0,jj+2
        do l=1,isv(j)
        do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
          vmfltd(i,j,1+mm)=0.
          vmfltd(i,j,2+mm)=0.
          vmfltd(i,j,3+mm)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      do k=4,kk
        km=k+mm
        kn=k+nn
c
c$OMP PARALLEL DO PRIVATE(l,i,flxhi,flxlo,q)
        do j=-1,jj+2
          do l=1,isu(j)
          do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
            flxhi= .125*min(dp(i-1,j,kn-1)*scp2(i-1,j),
     .                      dp(i  ,j,kn  )*scp2(i  ,j))
            flxlo=-.125*min(dp(i  ,j,kn-1)*scp2(i  ,j),
     .                      dp(i-1,j,kn  )*scp2(i-1,j))
            q=.25*(difint(i-1,j,k-1)+difint(i,j,k-1)
     .            +difint(i-1,j,k  )+difint(i,j,k  ))
            q=min(flxhi,max(flxlo,
     .            delt1*q*(p(i-1,j,k)-p(i,j,k))*scuy(i,j)*scuxi(i,j)))
            umfltd(i,j,km-1)=umfltd(i,j,km-1)+q
            umfltd(i,j,km  )=-q
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(l,i,flxhi,flxlo,q)
        do j=0,jj+2
          do l=1,isv(j)
          do i=max(-1,ifv(j,l)),min(ii+2,ilv(j,l))
            flxhi= .125*min(dp(i,j-1,kn-1)*scp2(i,j-1),
     .                      dp(i,j  ,kn  )*scp2(i,j  ))
            flxlo=-.125*min(dp(i,j  ,kn-1)*scp2(i,j  ),
     .                      dp(i,j-1,kn  )*scp2(i,j-1))
            q=.25*(difint(i,j-1,k-1)+difint(i,j,k-1)
     .            +difint(i,j-1,k  )+difint(i,j,k  ))
            q=min(flxhi,max(flxlo,
     .            delt1*q*(p(i,j-1,k)-p(i,j,k))*scvx(i,j)*scvyi(i,j)))
            vmfltd(i,j,km-1)=vmfltd(i,j,km-1)+q
            vmfltd(i,j,km  )=-q
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      enddo
c
      return
      end
