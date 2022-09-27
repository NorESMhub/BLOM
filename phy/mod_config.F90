! ------------------------------------------------------------------------------
! Copyright (C) 2020 Mats Bentsen, Ping-Gin Chiu
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

module mod_config
! ------------------------------------------------------------------------------
! This module contains configuration variables.
! ------------------------------------------------------------------------------

   implicit none

   private

   character(len = 256) :: &
      expcnf, &           ! Experiment configuration.
      runid               ! Experiment name.
   character(len = 16) :: &
      inst_name   = '', & ! Instance name.
      inst_suffix = ''    ! Instance suffix.
   integer :: &
      inst_index = 0      ! Instance index.
   logical :: &
      resume_flag = .false.    ! resume flag, use at ocn_run_mct()

   public ::  expcnf, runid, inst_name, inst_suffix, inst_index, resume_flag

end module mod_config
