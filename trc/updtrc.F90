! ------------------------------------------------------------------------------
! Copyright (C) 2007-2020 Mats Bentsen
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

subroutine updtrc(m,n,mm,nn,k1m,k1n)
  !
  ! --- ------------------------------------------------------------------
  ! --- update tracers due to non-passive processes
  ! --- ------------------------------------------------------------------
  !
  use mod_xc
  use mo_hamocc_step
  !
  implicit none
  !
  integer m,n,mm,nn,k1m,k1n
  !
#ifdef HAMOCC
  call hamocc_step(m,n,mm,nn,k1m,k1n)
#endif
#ifdef IDLAGE
  call idlage_step(m,n,mm,nn,k1m,k1n)
#endif
  !
  return
end subroutine updtrc
