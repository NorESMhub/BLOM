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

  use mod_types,   only: blom_i1 => i1, blom_i2 => i2, blom_i4 => i4, blom_i8 => i8, &
                         blom_r4 => r4, blom_r8 => r8, blom_kind_cs, blom_kind_cm,   &
                         blom_kind_cl, blom_kind_cx, blom_kind_cxx
  use mod_utility, only: fnmlen

  implicit none
  !----------------------------------------------------------------------------
  ! precision/kind constants add data public
  !----------------------------------------------------------------------------

  public
  integer,parameter :: i1 = blom_i1                     ! 8-bit integer
  integer,parameter :: i2 = blom_i2                     ! 16-bit integer
  integer,parameter :: i4 = blom_i4                     ! 32-bit integer
  integer,parameter :: i8 = blom_i8                     ! 64-bit integer
  integer,parameter :: r4 = blom_r4                     ! 4 byte real
  integer,parameter :: r8 = blom_r8                     ! 8 byte real
  integer,parameter :: HAMOCC_KIND_CS  = blom_kind_cs   ! short char
  integer,parameter :: HAMOCC_KIND_CM  = blom_kind_cm   ! mid-sized char
  integer,parameter :: HAMOCC_KIND_CL  = blom_kind_cl   ! long char
  integer,parameter :: HAMOCC_KIND_CX  = blom_kind_cx   ! extra-long char
  integer,parameter :: HAMOCC_KIND_CXX = blom_kind_cxx  ! extra-extra-long char

  integer,parameter :: bgc_fnmlen = fnmlen              ! default filename length
end module mo_kind
