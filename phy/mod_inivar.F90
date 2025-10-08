! ------------------------------------------------------------------------------
! Copyright (C) 2015-2024 Mats Bentsen, Jerry Tjiputra, Mariana Vertenstein
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

module mod_inivar

  use mod_vcoord,      only: inivar_vcoord
  use mod_state,       only: inivar_state
  use mod_pgforc,      only: inivar_pgforc
  use mod_momtum,      only: inivar_momtum
  use mod_barotp,      only: inivar_barotp
  use mod_tmsmt,       only: inivar_tmsmt
  use mod_diffusion,   only: inivar_diffusion
  use mod_difest,      only: inivar_difest
  use mod_eddtra,      only: inivar_eddtra
  use mod_utility,     only: inivar_utility
  use mod_mxlayr,      only: inivar_mxlayr
  use mod_seaice,      only: inivar_seaice
  use mod_forcing,     only: inivar_forcing
  use mod_cmnfld,      only: inivar_cmnfld
  use mod_niw,         only: inivar_niw
  use mod_swabs,       only: inivar_swabs
  use mod_tidaldissip, only: inivar_tidaldissip
  use mod_tracers,     only: inivar_tracers
  use mod_cesm,        only: inivar_cesm

  implicit none
  private

  public :: inivar

contains

  subroutine inivar

    ! -------------------------------------------------
    ! Call initialization routines for various modules.
    ! -------------------------------------------------

    call inivar_tracers
    call inivar_vcoord
    call inivar_state
    call inivar_pgforc
    call inivar_momtum
    call inivar_barotp
    call inivar_tmsmt
    call inivar_diffusion
    call inivar_difest
    call inivar_eddtra
    call inivar_utility
    call inivar_mxlayr
    call inivar_seaice
    call inivar_forcing
    call inivar_cmnfld
    call inivar_niw
    call inivar_swabs
    call inivar_tidaldissip
    call inivar_cesm

  end subroutine inivar

end module mod_inivar
