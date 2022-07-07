! ------------------------------------------------------------------------------
! Copyright (C) 2021 Mehmet Ilicak, Mats Bentsen
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

 module mod_single_column
! ----------------------------------------------------------------------
! This module contains routines for generating geometry information for
! single column test case.
! ----------------------------------------------------------------------

   use mod_types, only: r8
   use mod_xc
   use mod_vcoord, only: sigmar
   use mod_grid, only: qclon, qclat, pclon, pclat, uclon, uclat, vclon, vclat, &
                       scqx, scqy, scpx, scpy, scux, scuy, scvx, scvy, &
                       scq2, scp2, scu2, scv2, &
                       qlon, qlat, plon, plat, ulon, ulat, vlon, vlat, &
                       depths, corioq, coriop, betafp, angle, cosang, sinang, &
                       nwp
   implicit none
   public :: geoenv_single_column

contains

   ! ---------------------------------------------------------------------------
   ! Public procedures.
   ! ---------------------------------------------------------------------------
      subroutine geoenv_single_column
! --- ------------------------------------------------------------------
! --- Define bathymetry, grid specification and Coriolis parameter for
! --- single column case [cm s]
! --- ------------------------------------------------------------------
      depths=1000._r8
!
      nwp=jtdm*itdm
      qlon = -165.6_r8
      qlat = 0._r8
      plon = -165.5_r8
      plat = 0._r8
      ulon = -165.6_r8
      ulat = -0.9_r8
      vlon = -165.5_r8
      vlat = -1.1_r8
      qclon = 0._r8
      qclat = 0._r8
      pclon = 0._r8
      pclat = 0._r8
      uclon = 0._r8
      uclat = 0._r8
      vclon = 0._r8
      vclat = 0._r8
      scqx = 1100000.0_r8
      scqy = 1100000.0_r8
      scpx = 1100000.0_r8
      scpy = 1100000.0_r8
      scux = 1100000.0_r8
      scuy = 1100000.0_r8
      scvx = 1100000.0_r8
      scvy = scuy
      scq2 = scqx*scqy
      scp2 = scpx*scpy
      scu2 = scux*scuy
      scv2 = scvx*scvy
      corioq = 0._r8
      coriop = 0._r8
      betafp = 0._r8
      angle = 0._r8
      cosang = 1._r8
      sinang = 0._r8

      return
      end subroutine geoenv_single_column

! --- ------------------------------------------------------------------

end module mod_single_column
