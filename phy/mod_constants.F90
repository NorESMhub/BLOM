! ------------------------------------------------------------------------------
! Copyright (C) 2020-2025 Mats Bentsen, Mehmet Ilicak
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

module mod_constants
! ------------------------------------------------------------------------------
! This module contains physical constants and various other constant parameters.
! ------------------------------------------------------------------------------

   use mod_types, only: r8

   implicit none

   private

   real(r8), parameter :: &
      grav   = 9.806_r8, &        ! Gravitational acceleration [m s-2].
      rearth = 6.37122e6_r8, &    ! Radius of the Earth [m].
      spcifh = 3990._r8, &        ! Specific heat capacity of sea water
                                  ! [J kg-1 K-1].
      t0deg  = 273.15_r8, &       ! Zero degrees Celsius in Kelvin [K].
      alpha0 = 1.e-3_r8, &        ! Reference value of specific volume
                                  ! [m3 kg-1].
      rho0   = 1.e3_r8, &         ! Reference value of density [kg m-3].
      pi     = 3.1415926536_r8, & ! pi [].
      radian = 57.295779513_r8, & ! 180/pi [].
      epsilpl = 1.e-14_r8, &      ! Small value for pressure*dx [].
      epsilp  = 1.e-12_r8, &      ! Small value for pressure [].
      epsilz  = 1.e-9_r8, &       ! Small value for depth [].
      epsilt  = 1.e-11_r8, &      ! Small value for time [].
      epsilk  = 1.e-15_r8, &      ! Small value for kappa [].
      spval  = 1.e33_r8, &        ! Large value [].
      tenm   = 98060._r8, &       ! 10 m in units of pressure [kg m-1 s-2].
      onem   = 9806._r8, &        ! 1 m in units of pressure [kg m-1 s-2].
      tencm  = 980.6_r8, &        ! 10 cm in units of pressure [kg m-1 s-2].
      onecm  = 98.06_r8, &        ! 1 cm in units of pressure [kg m-1 s-2].
      onemm  = 9.806_r8, &        ! 1 mm in units of pressure [kg m-1 s-2].
      onemu  = .009806_r8, &      ! 1 micrometer in units of pressure
                                  ! [kg m-1 s-2].
      g2kg   = 1.e-3_r8, &        ! convert g to kg coeff
      kg2g   = 1.e3_r8            ! convert kg to g coeff

   public :: grav, rearth, spcifh, t0deg, alpha0, rho0, pi, radian, &
             epsilpl, epsilp, epsilz, epsilt, epsilk, spval, &
             tenm, onem, tencm, onecm, onemm, onemu, g2kg, kg2g

end module mod_constants
