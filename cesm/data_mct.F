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

      module data_mct

      ! Data set by mct to be used by blom

      implicit none

      private

      character (len=256) :: runid_mct  ! case name
      character (len=256) :: runtyp_mct ! run type
      integer :: ocn_cpl_dt_mct         ! Coupling frequency

      public :: runid_mct, runtyp_mct, ocn_cpl_dt_mct

      end module data_mct
