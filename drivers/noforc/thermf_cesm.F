! ------------------------------------------------------------------------------
! Copyright (C) 2009-2020 Mats Bentsen
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

      subroutine thermf_cesm(m,n,mm,nn,k1m,k1n)
c
c --- NERSC version of thermf. To be used when coupled to CESM
c
      use mod_xc
      use mod_forcing, only: hmltfz, surflx, surrlx, sswflx, salflx,
     .                       brnflx, salrlx, ustar
      use mod_cesm, only: frzpot, mltpot
#ifdef TRC
      use mod_tracers, only: ntrc, trflx
#endif
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
      integer i,j,k,l,nt
#ifdef TRC
      integer nt
#endif
c
c$OMP PARALLEL DO PRIVATE(l,i
#ifdef TRC
c$OMP+ ,nt
#endif
c$OMP+ )
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
          surflx(i,j)=0.
          sswflx(i,j)=0.
          salflx(i,j)=0.
          brnflx(i,j)=0.
          surrlx(i,j)=0.
          salrlx(i,j)=0.
          hmltfz(i,j)=0.
          ustar(i,j)=0.
          frzpot(i,j)=0.
          mltpot(i,j)=0.
#ifdef TRC
          do nt=1,ntrc
            trflx(nt,i,j)=0.
          enddo
#endif
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
