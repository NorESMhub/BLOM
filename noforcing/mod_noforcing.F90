! ------------------------------------------------------------------------------
! Copyright (C) 2025 Mats Bentsen
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

module mod_noforcing
! ------------------------------------------------------------------------------
! This module contains variables and routines related to experiment
! configuration with no surface forcing.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_xc
   use mod_forcing, only: surflx, surrlx, sswflx, salflx, brnflx, salrlx, &
                          ustar, taux, tauy

   implicit none

   private

   public :: inifrc_noforcing

contains

   ! ---------------------------------------------------------------------------
   ! Public procedures.
   ! ---------------------------------------------------------------------------

   subroutine inifrc_noforcing
   ! ---------------------------------------------------------------------------
   ! Set all forcing arrays to zero at ocean points.
   ! ---------------------------------------------------------------------------

      integer :: i, j, l

      !$omp parallel do private(l, i)
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
            surflx(i,j) = 0._r8
            sswflx(i,j) = 0._r8
            salflx(i,j) = 0._r8
            brnflx(i,j) = 0._r8
            surrlx(i,j) = 0._r8
            salrlx(i,j) = 0._r8
            ustar (i,j) = 0._r8
         enddo
         enddo
         do l = 1, isu(j)
         do i = max(1, ifu(j,l)), min(ii, ilu(j,l))
            taux(i,j) = 0._r8
         enddo
         enddo
         do l = 1, isv(j)
         do i = max(1, ifv(j,l)), min(ii, ilv(j,l))
            tauy(i,j) = 0._r8
         enddo
         enddo
       enddo
       !$omp end parallel do

   end subroutine inifrc_noforcing

end module mod_noforcing
