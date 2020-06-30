! ------------------------------------------------------------------------------
! Copyright (C) 2018 Mats Bentsen
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
c
c --- ------------------------------------------------------------------
c --- This module contains routines for computing the freezing point of
c --- sea water.
c --- ------------------------------------------------------------------
c
      implicit none
c
      private
c
      public :: swtfrz
c
      interface swtfrz
         module procedure swtfrz_0d
         module procedure swtfrz_1d
         module procedure swtfrz_2d
      end interface swtfrz
c
      contains
c
c --- ------------------------------------------------------------------
c
      function swtfrz_0d(p,s) result(swtfrz)
c
c --- ------------------------------------------------------------------
c --- Freezing temperature as function of pressure and salinity.
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_eos.h"
c
      real, intent(in) :: p      ! Pressure [g cm-1 s-2]
      real, intent(in) :: s      ! Salinity [g kg-1]
      real :: swtfrz
c
      swtfrz=atf*s+btf+ctf*p
c
      end function swtfrz_0d
c
c --- ------------------------------------------------------------------
c
      function swtfrz_1d(p,s) result(swtfrz)
c
c --- ------------------------------------------------------------------
c --- Freezing temperature as function of pressure and salinity.
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_eos.h"
c
      real, intent(in) :: p(:)   ! Pressure [g cm-1 s-2]
      real, intent(in) :: s(:)   ! Salinity [g kg-1]
      real :: swtfrz(size(s))
c
      swtfrz=atf*s+btf+ctf*p
c
      end function swtfrz_1d
c
c --- ------------------------------------------------------------------
c
      function swtfrz_2d(p,s) result(swtfrz)
c
c --- ------------------------------------------------------------------
c --- Freezing temperature as function of pressure and salinity.
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_eos.h"
c
      real, intent(in) :: p(:,:) ! Pressure [g cm-1 s-2]
      real, intent(in) :: s(:,:) ! Salinity [g kg-1]
      real :: swtfrz(size(s,1),size(s,2))
c
      swtfrz=atf*s+btf+ctf*p
c
      end function swtfrz_2d
c
c --- ------------------------------------------------------------------
c
      end module mod_swtfrz
