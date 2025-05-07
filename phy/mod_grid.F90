! ------------------------------------------------------------------------------
! Copyright (C) 2020-2025 Mats Bentsen, Aleksi Nummmelin
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

module mod_grid
! ------------------------------------------------------------------------------
! This module contains declaration of grid variables.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_constants, only: spval
   use mod_xc, only: idm, jdm, kdm, nbdy, ii, jj, kk
   use mod_utility, only: fnmlen

   implicit none

   private

   ! Variable to be set in namelist:
   character(len = fnmlen) :: &
      grfile     ! Name of file containing grid specification.

   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, 4) :: &
      qclon, &   ! Longitude of q-cell corners [degrees].
      qclat, &   ! Latitude of q-cell corners [degrees].
      pclon, &   ! Longitude of p-cell corners [degrees].
      pclat, &   ! Latitude of p-cell corners [degrees].
      uclon, &   ! Longitude of u-cell corners [degrees].
      uclat, &   ! Latitude of u-cell corners [degrees].
      vclon, &   ! Longitude of v-cell corners [degrees].
      vclat      ! Latitude of v-cell corners [degrees].

   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy) :: &
      scqx, &    ! Grid size in x-direction centered at q-point [m].
      scqy, &    ! Grid size in y-direction centered at q-point [m].
      scpx, &    ! Grid size in x-direction centered at p-point [m].
      scpy, &    ! Grid size in y-direction centered at p-point [m].
      scux, &    ! Grid size in x-direction centered at u-point [m].
      scuy, &    ! Grid size in y-direction centered at u-point [m].
      scvx, &    ! Grid size in x-direction centered at v-point [m].
      scvy, &    ! Grid size in y-direction centered at v-point [m].
      scq2, &    ! Area of grid cell centered at q-point [m2].
      scp2, &    ! Area of grid cell centered at p-point [m2].
      scu2, &    ! Area of grid cell centered at u-point [m2].
      scv2, &    ! Area of grid cell centered at v-point [m2].
      scq2i, &   ! Multiplicative inverse of scq2 [m-2].
      scp2i, &   ! Multiplicative inverse of scp2 [m-2].
      scuxi, &   ! Multiplicative inverse of scux [m-1].
      scuyi, &   ! Multiplicative inverse of scuy [m-1].
      scvxi, &   ! Multiplicative inverse of scvx [m-1].
      scvyi, &   ! Multiplicative inverse of scvy [m-1].
      qlon, &    ! Longitude of q-point [degrees].
      qlat, &    ! Latitude of q-point [degrees].
      plon, &    ! Longitude of p-point [degrees].
      plat, &    ! Latitude of p-point [degrees].
      ulon, &    ! Longitude of u-point [degrees].
      ulat, &    ! Latitude of u-point [degrees].
      vlon, &    ! Longitude of v-point [degrees].
      vlat, &    ! Latitude of v-point [degrees].
      depths, &  ! Water depth [m].
      corioq, &  ! Coriolis parameter at q-point [s-1].
      coriop, &  ! Coriolis parameter at p-point [s-1].
      betafp, &  ! Derivative of Coriolis parameter with respect to meridional
                 ! distance at p-point [m-1 s-1].
      betatp, &  ! Topographic Rhines scale [m-1 s-1].
      angle, &   ! Local angle between x-direction and eastward direction at
                 ! p-points [radians].
      cosang, &  ! Cosine of local angle between x-direction and eastward
                 ! direction at p-points [].
      sinang, &  ! Sine of local angle between x-direction and eastward
                 ! direction at p-points [].
      hangle     ! Angle between the bottom slope vector and local
                 ! x-direction [radians]

   real(r8) :: &
      area       ! Total grid area [m2].

   integer :: &
      nwp        ! Number of wet grid cells.

   public :: grfile, &
             qclon, qclat, pclon, pclat, uclon, uclat, vclon, vclat, &
             scqx, scqy, scpx, scpy, scux, scuy, scvx, scvy, &
             scq2, scp2, scu2, scv2, scq2i, scp2i, &
             scuxi, scuyi, scvxi, scvyi, &
             qlon, qlat, plon, plat, ulon, ulat, vlon, vlat, &
             depths, corioq, coriop, betafp, betatp, &
             angle, cosang, sinang, hangle, &
             area, nwp, &
             inivar_grid

contains

   subroutine inivar_grid
   ! ---------------------------------------------------------------------------
   ! Initialize arrays.
   ! ---------------------------------------------------------------------------

      scqx(:,:) = spval
      scqy(:,:) = spval
      scpx(:,:) = spval
      scpy(:,:) = spval
      scux(:,:) = spval
      scuy(:,:) = spval
      scvx(:,:) = spval
      scvy(:,:) = spval
      scq2(:,:) = spval
      scp2(:,:) = spval
      scu2(:,:) = spval
      scv2(:,:) = spval
      scq2i(:,:) = spval
      scp2i(:,:) = spval
      scuxi(:,:) = spval
      scuyi(:,:) = spval
      scvxi(:,:) = spval
      scvyi(:,:) = spval
      qlon(:,:) = spval
      qlat(:,:) = spval
      plon(:,:) = spval
      plat(:,:) = spval
      ulon(:,:) = spval
      ulat(:,:) = spval
      vlon(:,:) = spval
      vlat(:,:) = spval
      depths(:,:) = spval
      corioq(:,:) = spval
      coriop(:,:) = spval
      betafp(:,:) = spval
      betatp(:,:) = spval
      angle(:,:) = spval
      cosang(:,:) = spval
      sinang(:,:) = spval
      hangle(:,:) = spval
      qclon(:,:,:) = spval
      qclat(:,:,:) = spval
      pclon(:,:,:) = spval
      pclat(:,:,:) = spval
      uclon(:,:,:) = spval
      uclat(:,:,:) = spval
      vclon(:,:,:) = spval
      vclat(:,:,:) = spval

   end subroutine inivar_grid

end module mod_grid
