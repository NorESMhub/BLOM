! ------------------------------------------------------------------------------
! Copyright (C) 2015-2025 Mats Bentsen, Mehmet Ilicak, Aleksi Nummelin,
!                         Mariana Vertenstein
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

module mod_inifrc

  use mod_config,    only: expcnf
  use mod_xc,        only: lp, mnproc, xcstop
  use mod_cesm,      only: inifrc_cesm
  use mod_ben02,     only: inifrc_ben02clim, inifrc_ben02syn
  use mod_noforcing, only: inifrc_noforcing
  use mod_fuk95,     only: inifrc_fuk95
  use mod_channel,   only: inifrc_channel

  implicit none
  private

  public :: inifrc

contains

  subroutine inifrc
  ! ---------------------------------------------------------------------------
  ! Initialize forcing.
  ! ---------------------------------------------------------------------------

    select case (trim(expcnf))
      case ('cesm')
        call inifrc_cesm
      case ('ben02clim')
        call inifrc_ben02clim
      case ('ben02syn')
        call inifrc_ben02syn
      case ('noforcing')
        call inifrc_noforcing
      case ('fuk95')
        call inifrc_fuk95
      case ('channel')
        call inifrc_channel
      case ('isomip1')
        !        call inifrc_isomip1
      case ('isomip2')
        !        call inifrc_isomip2
      case ('single_column')
        call inifrc_ben02clim
      case default
        if (mnproc == 1) then
          write (lp,'(3a)') ' inifrc: expcnf = ', trim(expcnf), &
               ' is unsupported!'
        endif
        call xcstop('(inifrc)')
        stop '(inifrc)'
    end select

  end subroutine inifrc

end module mod_inifrc
