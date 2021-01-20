! ------------------------------------------------------------------------------
! Copyright (C) 2020 Mats Bentsen
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

module mod_pointtest
! ------------------------------------------------------------------------------
! This module contains variables and procedures related to detailed point
! diagnostics.
! ------------------------------------------------------------------------------

   use mod_xc, only: i0, j0, ii, jj, lp, mnproc
   implicit none

   private

   integer :: &
      itest = 60, & ! Global i-index of point diagnostics.
      jtest = 60, & ! Global j-index of point diagnostics.
      ptest = 0     ! Processor containing the test point.

   public :: itest, jtest, ptest, init_ptest

contains

   subroutine init_ptest
   ! ---------------------------------------------------------------------------
   ! Identify processor and local horizontal indexes where detailed diagnostics
   ! are desired.
   ! ---------------------------------------------------------------------------

       if (itest > i0 .and. itest <= i0 + ii .and. &
           jtest > j0 .and. jtest <= j0 + jj) then
          write (lp, '(a,i4,a,i4,a,i5)') &
             ' itest = ', itest,', jtest = ', jtest, &
             ' found on processor ', mnproc
          call flush(lp)
          ptest = mnproc
          itest = itest - i0
          jtest = jtest - j0
       endif

   end subroutine init_ptest

end module mod_pointtest
