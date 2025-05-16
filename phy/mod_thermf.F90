! ------------------------------------------------------------------------------
! Copyright (C) 2015-2025 Mats Bentsen, Mehmet Ilicak, Aleksi Nummelin

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

module mod_thermf

  use mod_config,         only: expcnf
  use mod_xc,             only: lp, mnproc, xcstop
  use mod_thermf_cesm,    only: thermf_cesm
  use mod_thermf_ben02,   only: thermf_ben02
  use mod_thermf_channel, only: thermf_channel

  implicit none
  private

  public :: thermf

contains

  subroutine thermf(m,n,mm,nn,k1m,k1n)
  ! ----------------------------------------------------------------------------
  ! Get surface forcing functions.
  ! ----------------------------------------------------------------------------

    ! Arguments
    integer, intent(in) :: m,n,mm,nn,k1m,k1n

    select case (trim(expcnf))
      case ('cesm')
        call thermf_cesm(m,n,mm,nn,k1m,k1n)
      case ('ben02clim', 'ben02syn', 'single_column')
        call thermf_ben02(m,n,mm,nn,k1m,k1n)
      case ('noforcing')
      case ('channel')
        call thermf_channel(m,n,mm,nn,k1m,k1n)
      case ('fuk95')
      case ('isomip1')
        ! call thermf_isomip1(m,n,mm,nn,k1m,k1n)
      case ('isomip2')
        ! call thermf_isomip2(m,n,mm,nn,k1m,k1n)
      case default
        if (mnproc == 1) then
          write (lp,'(3a)') ' thermf: expcnf = ', trim(expcnf), &
                            ' is unsupported!'
        end if
        call xcstop('(thermf)')
        stop '(thermf)'
    end select

  end subroutine thermf

end module mod_thermf
