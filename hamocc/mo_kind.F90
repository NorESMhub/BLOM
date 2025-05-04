! Copyright (C) 2025  J. Maerz, J. Schwinger, T. Torsvik
!
! This file is part of BLOM/iHAMOCC.
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
! along with BLOM. If not, see https://www.gnu.org/licenses/.

module mo_kind

  use, intrinsic :: iso_fortran_env, only: int8, int16, int32, int64, real32, real64

  implicit none
  !----------------------------------------------------------------------------
  ! precision/kind constants add data public
  !----------------------------------------------------------------------------

  public
  integer,parameter :: i1 = int8               ! 8-bit integer
  integer,parameter :: i2 = int16              ! 16-bit integer
  integer,parameter :: i4 = int32              ! 32-bit integer
  integer,parameter :: i8 = int64              ! 64-bit integer
  integer,parameter :: r4 = real32             ! 4 byte real
  integer,parameter :: r8 = real64             ! 8 byte real
  integer,parameter :: HAMOCC_KIND_CS = 80     ! short char
  integer,parameter :: HAMOCC_KIND_CM = 160    ! mid-sized char
  integer,parameter :: HAMOCC_KIND_CL = 256    ! long char
  integer,parameter :: HAMOCC_KIND_CX = 512    ! extra-long char
  integer,parameter :: HAMOCC_KIND_CXX= 4096   ! extra-extra-long char

  integer,parameter :: bgc_fnmlen = 512        ! default filename length
end module mo_kind
