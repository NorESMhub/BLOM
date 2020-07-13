!!  This file is modified from cam_instance.F90 in CAM component.
!!  Intent to provide multiple-instance support of BLOM component in NorESM2.
!!                                                      Ping-Gin Jul 2020 

module blom_instance

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
   print*,'ocn_id,name,index,suffix: ',ocn_id,trim(inst_name),inst_index,trim(inst_suffix)

end subroutine blom_instance_init

end module blom_instance
