! ------------------------------------------------------------------------------
! Copyright (C) 2010-2020 Ingo Bethke, Mats Bentsen, Mehmet Ilicak,
!                         Alok Kumar Gupta, Jörg Schwinger
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

use seq_comm_mct, only: seq_comm_suffix, seq_comm_inst, seq_comm_name

implicit none
private
save

public :: blom_instance_init

integer,           public :: ocn_id
integer,           public :: inst_index
character(len=16), public :: inst_name
character(len=16), public :: inst_suffix

!===============================================================================
CONTAINS
!===============================================================================

subroutine blom_instance_init(in_ocn_id)

   integer, intent(in) :: in_ocn_id

   ocn_id      = in_ocn_id
   inst_name   = seq_comm_name(ocn_id)
   inst_index  = seq_comm_inst(ocn_id)
   inst_suffix = seq_comm_suffix(ocn_id)

end subroutine blom_instance_init

end module mod_instance
