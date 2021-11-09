! ------------------------------------------------------------------------------
! Copyright (C) 2021 Tomas Torsvik
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

subroutine restart_getfile(rstfnm_in, rstlabel, rstfnm_out, rstfnm_err)
  !-----------------------------------------------------------------------------
  ! Generate filename for restart files to read or write tracer fields
  !-----------------------------------------------------------------------------

  use mod_config, only: expcnf

  implicit none

  character(len=*), intent(in)  :: rstfnm_in     ! Original restart file name
  character(len=*), intent(in)  :: rstlabel      ! Label to insert in new file
  character(len=*), intent(out) :: rstfnm_out    ! New restart file name
  logical,          intent(out) :: rstfnm_err    ! Error flag

  integer :: i_suffix, i_time, i_restart
  character(len=:), allocatable  :: str_suffix, str_time, str_restart

  rstfnm_err = .false.

  if (expcnf.eq.'cesm') then
     ! Assume file format: <str_restart.>'r'<.str_timestamp><.str_suffix>
     ! Search for '.' starting from end of "rstfnm_in" filename
     !-- File suffix
     i_suffix = index(rstfnm_in, '.', back=.true.)
     str_suffix = trim(rstfnm_in(i_suffix:))
     !-- File timestamp
     i_time = index(rstfnm_in(:(i_suffix-1)), '.', back=.true.)
     str_time = rstfnm_in(i_time:(i_suffix-1))
     !-- File without original restart label
     i_restart = index(rstfnm_in(:(i_time-1)), '.', back=.true.)
     str_restart = rstfnm_in(:i_restart)

     if (i_suffix == 0 .or. i_time == 0 .or. i_restart == 0) then
        rstfnm_err = .true.
     else
        rstfnm_out = str_restart // trim(rstlabel) // str_time // str_suffix
     end if
  else
     ! Assume file format: <str_restart_>'restphy'<_str_suffix>
     ! Search for '_' starting from end of "rstfnm_in" filename
     !-- File suffix
     i_suffix = index(rstfnm_in, '_', back=.true.)
     str_suffix = trim(rstfnm_in(i_suffix:))
     !-- File without original restart label
     i_restart = index(rstfnm_in(:(i_suffix-1)), '_', back=.true.)
     str_restart = rstfnm_in(:i_restart)

     if (i_suffix == 0 .or. i_restart == 0) then
        rstfnm_err = .true.
     else
        rstfnm_out = str_restart // trim(rstlabel) // str_suffix
     end if
  end if
end subroutine restart_getfile
