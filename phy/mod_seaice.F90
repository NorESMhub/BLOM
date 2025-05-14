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

module mod_seaice
! ------------------------------------------------------------------------------
! This module contains variables related to sea ice.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_constants, only: spval
   use mod_xc

   implicit none

   private

   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy) :: &
      ficem, &   ! Sea-ice concentration [].
      hicem, &   ! Sea-ice thickness [m].
      hsnwm, &   ! Snow thickness on sea ice [m].
      ustari, &  ! Friction velocity at ocean and sea ice interface [m s-1].
      tauxice, & ! x-component of momentum flux from sea-ice to ocean [N m-2].
      tauyice, & ! y-component of momentum flux from sea-ice to ocean [N m-2].
      uicem, &   ! x-component of sea-ice velocity [m s-1].
      vicem, &   ! y-component of sea-ice velocity [m s-1].
      iagem      ! Age of sea ice [days].

   public :: ficem, hicem, hsnwm, ustari, tauxice, tauyice, uicem, vicem, &
             iagem, inivar_seaice

contains

   subroutine inivar_seaice
   ! ---------------------------------------------------------------------------
   ! Initialize variables related to sea ice.
   ! ---------------------------------------------------------------------------

      integer :: i, j, l

      ficem(:,:) = spval
      hicem(:,:) = spval
      hsnwm(:,:) = spval
      ustari(:,:) = spval
      tauxice(:,:) = spval
      tauyice(:,:) = spval
      uicem(:,:) = spval
      vicem(:,:) = spval
      iagem(:,:) = spval

   !$omp parallel do private(l, i)
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
            ustari(i, j) = 0._r8
            tauxice(i, j) = 0._r8
            tauyice(i, j) = 0._r8
            uicem(i, j) = 0._r8
            vicem(i, j) = 0._r8
         enddo
         enddo
      enddo
   !$omp end parallel do

   end subroutine inivar_seaice

end module mod_seaice
