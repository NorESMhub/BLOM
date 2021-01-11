! ------------------------------------------------------------------------------
! Copyright (C) 2015-2020 Mats Bentsen
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

      subroutine ocntrc_init
c
c --- ------------------------------------------------------------------
c --- initialization of ocean tracers
c --- ------------------------------------------------------------------
c
      use mod_types, only: r8
      use mod_xc
      use mod_grid, only: plon, plat
      use mod_tracers, only: ntrocn, natr, ntr, trc
c
      implicit none
c
      integer i,j,k,l,nt,nat
c
c --- ------------------------------------------------------------------
c --- if no ocean tracers are defined, return
c --- ------------------------------------------------------------------
c
      if (ntrocn.eq.0) return
c
c --- ------------------------------------------------------------------
c --- if number of age tracers is greater than zero, CPP flag ATRC must
c --- be defined
c --- ------------------------------------------------------------------
#ifndef ATRC
      if (natr.gt.0) then
        if (mnproc.eq.1) then
          write (lp,'(2a)') ' Since number of age tracers is greater ',
     .      'than zero, ATRC must be defined!'
        endif
        call xcstop('(ocntrc_init)')
               stop '(ocntrc_init)'
      endif
#endif
c
c --- ------------------------------------------------------------------
c --- check number of age tracers
c --- ------------------------------------------------------------------
c
      if (natr.lt.0.or.2*natr.gt.ntrocn) then
        if (mnproc.eq.1) then
          write (lp,'(3a)') ' Number of age tracers must be greater ',
     .      'than zero and less or equal half the total number of ',
     .      'ocean tracers!'
        endif
        call xcstop('(ocntrc_init)')
               stop '(ocntrc_init)'
      endif
c
c --- ------------------------------------------------------------------
c --- initialization of tracers
c --- ------------------------------------------------------------------
c
      do nt=1,ntrocn-natr
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=1,kk
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              trc(i,j,k,nt)=
     .          (mod(k,5)+1)*(plat(i,j)+90._r8)/(5._r8*180._r8)+nt
              trc(i,j,k+kk,nt)=trc(i,j,k,nt)
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c
c --- ------------------------------------------------------------------
c --- initialization of age tracers
c --- ------------------------------------------------------------------
c
      do nt=1,natr
        nat=ntr-natr+nt
c$OMP PARALLEL DO PRIVATE(k,l,i)
        do j=1,jj
          do k=1,kk
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              trc(i,j,k,nat)=
     .          (mod(k,5)*(plon(i,j)+180._r8)/(4._r8*360._r8)+nat)
     .          *trc(i,j,k,nt)
              trc(i,j,k+kk,nat)=trc(i,j,k,nat)
            enddo
            enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c
      end
