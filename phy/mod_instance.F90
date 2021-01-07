! ------------------------------------------------------------------------------
! Copyright (C) 2020 Ping-Gin Chiu and Mats Bentsen
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

module mod_instance

   implicit none

   private

   integer             :: inst_index  = 0
   character(len = 16) :: inst_name   = ''
   character(len = 16) :: inst_suffix = ''

   public :: inst_index, inst_name, inst_suffix

end module mod_instance
