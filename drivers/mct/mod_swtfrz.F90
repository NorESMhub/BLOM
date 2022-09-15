! ------------------------------------------------------------------------------
! Copyright (C) 2018-2022 Mats Bentsen
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

module mod_swtfrz
! ------------------------------------------------------------------------------
! This module contains routines for computing the freezing point of sea water.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use shr_frz_mod, only: shr_frz_freezetemp

   implicit none

   private

   public :: swtfrz

   interface swtfrz
      module procedure swtfrz_0d
      module procedure swtfrz_1d
      module procedure swtfrz_2d
   end interface swtfrz

contains

   function swtfrz_0d(p,s) result(swtfrz)
   ! ---------------------------------------------------------------------------
   ! Retrieve freezing temperature from shared CESM function.
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: p      ! Pressure [g cm-1 s-2]
      real(r8), intent(in) :: s      ! Salinity [g kg-1]
      real(r8) :: swtfrz

      swtfrz = shr_frz_freezetemp(s)

   end function swtfrz_0d

   function swtfrz_1d(p,s) result(swtfrz)
   ! ---------------------------------------------------------------------------
   ! Retrieve freezing temperature from shared CESM function.
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: p(:)   ! Pressure [g cm-1 s-2]
      real(r8), intent(in) :: s(:)   ! Salinity [g kg-1]
      real(r8) :: swtfrz(size(s))

      swtfrz(:) = shr_frz_freezetemp(s(:))

   end function swtfrz_1d

   function swtfrz_2d(p,s) result(swtfrz)
   ! ---------------------------------------------------------------------------
   ! Retrieve freezing temperature from shared CESM function.
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: p(:,:) ! Pressure [g cm-1 s-2]
      real(r8), intent(in) :: s(:,:) ! Salinity [g kg-1]
      real(r8) :: swtfrz(size(s,1),size(s,2))

      swtfrz(:,:) = shr_frz_freezetemp(s(:,:))

   end function swtfrz_2d

end module mod_swtfrz
