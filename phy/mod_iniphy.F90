! ------------------------------------------------------------------------------
! Copyright (C) 2015-2021 Mats Bentsen, Mehmet Ilicak, Aleksi Nummelin

! This file is part of BLOM.

! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.

! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.

! You should have received a copy of the GNU Lesser General Public License
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

module mod_iniphy

  use mod_config,      only: expcnf
  use mod_xc,          only: lp, mnproc, xcstop
  use mod_vcoord,      only: vcoord_type_tag, cntiso_hybrid
  use mod_tidaldissip, only: read_tidaldissip
  use mod_difest,      only: init_difest

  implicit none
  private

  public :: iniphy

contains

  subroutine iniphy

    ! --- ------------------------------------------------------------------
    ! --- Initialize physical parameterizations
    ! --- ------------------------------------------------------------------

    if (expcnf == 'cesm'.or. &
       expcnf == 'ben02clim'.or.expcnf == 'ben02syn') then
      call read_tidaldissip
    else if (expcnf == 'channel') then
    else if (expcnf == 'fuk95') then
    else if (expcnf == 'single_column') then
    else if (expcnf == 'isomip1') then
    else if (expcnf == 'isomip2') then
    else
      if (mnproc == 1) then
        write (lp,'(3a)') ' expcnf = ',trim(expcnf),' is unsupported!'
      end if
      call xcstop('(iniphy)')
      stop '(iniphy)'
    end if

    if (vcoord_type_tag == cntiso_hybrid) then
      call init_difest
    end if

  end subroutine iniphy

end module mod_iniphy
