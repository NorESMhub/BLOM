! ------------------------------------------------------------------------------
! Copyright (C) 2008-2020 Mats Bentsen
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

module mod_types
! ------------------------------------------------------------------------------
! This module defines numeric data types.
! ------------------------------------------------------------------------------

   use, intrinsic :: iso_fortran_env, only: &
      int8, int16, int32, int64, real32, real64

   implicit none

   private

   integer, parameter :: &
      i1 = int8, &
      i2 = int16, &
      i4 = int32, &
      i8 = int64, &
      r4 = real32, &
      r8 = real64
   integer,parameter :: BLOM_CHAR_S = 80      ! short char
   integer,parameter :: BLOM_CHAR_M = 160     ! mid-sized char
   integer,parameter :: BLOM_CHAR_L = 256     ! long char
   integer,parameter :: BLOM_CHAR_X = 512     ! extra-long char
   integer,parameter :: BLOM_CHAR_XX= 4096    ! extra-extra-long char

   public :: i1, i2, i4, i8, r4, r8,  &
        BLOM_CHAR_S,                  &
        BLOM_CHAR_M,                  &
        BLOM_CHAR_L,                  &
        BLOM_CHAR_X,                  &
        BLOM_CHAR_XX

end module mod_types
