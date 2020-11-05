! ------------------------------------------------------------------------------
! Copyright (C) 2008 Mats Bentsen
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

module types

   ! -------------------------------------------------------------------
   ! This module defines numerical data types
   ! -------------------------------------------------------------------

   implicit none

   integer, parameter, public :: &
      i1 = selected_int_kind(2), &
      i2 = selected_int_kind(4), &
      i4 = selected_int_kind(9), &
      i8 = selected_int_kind(18), &
      r4 = selected_real_kind(6,37), &
      r8 = selected_real_kind(13,307)

end module types
