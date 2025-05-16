! ------------------------------------------------------------------------------
! Copyright (C) 2015-2025 Mats Bentsen, Mehmet Ilicak
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

module mod_cmnfld
! ------------------------------------------------------------------------------
! This module contains variables and procedures related to common fields used by
! several subsequent routines.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_constants, only: spval, onem
   use mod_xc

   implicit none

   private

   ! Parameters:
   real(r8) :: &
      sls0 = 10._r8*onem, &      ! Minimum smoothing length scale in the
                                 ! computation of filtered BFSQ [kg m-1 s-2].
      slsmfq = 2._r8, &          ! Factor to be multiplied with the mixed
                                 ! layer depth to find the smoothing length
                                 ! scale at the base of the mixed layer in the
                                 ! computation of filtered BFSQ [].
      slsels = 2._r8, &          ! Factor to be multiplied with the mixed
                                 ! layer depth to find the e-folding length
                                 ! scale of the smoothing length scale in the
                                 ! computation of filtered BFSQ [].
      bfsqmn = 1.e-7_r8, &       ! Minimum value of BFSQ used in the
                                 ! computation of neutral slope [s-2].
      dbcrit = .0003_r8          ! Critical buoyancy difference used in the
                                 ! mixed layer thickness estimation (Levitus,
                                 ! 1982) [m s-2].

   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, kdm + 1) :: &
      bfsqi, &                   ! Interface buoyancy frequency squared [s-2].
      bfsqf, &                   ! Filtered interface buoyancy frequency
                                 ! squared [s-2].
      z                          ! Interface depth [m].
   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, kdm) :: &
      bfsql, &                   ! Layer buoyancy frequency squared [s-2].
      nslpx, &                   ! x-component of local neutral slope [].
      nslpy, &                   ! y-component of local neutral slope [].
      nnslpx, &                  ! x-component of local neutral slope times
                                 ! buoyancy frequency [s-1].
      nnslpy, &                  ! y-component of local neutral slope times
                                 ! buoyancy frequency [s-1].
      dz                         ! Layer thickness [m].
   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy) :: &
      mlts                       ! Mixed layer depth defined by density
                                 ! criterion [m].

   public :: sls0, slsmfq, slsels, bfsqmn, dbcrit, &
             bfsqi, bfsqf, z, bfsql, nslpx, nslpy, nnslpx, nnslpy, dz, mlts, &
             inivar_cmnfld

   contains

   ! ---------------------------------------------------------------------------
   ! Public procedures.
   ! ---------------------------------------------------------------------------

   subroutine inivar_cmnfld
   ! ---------------------------------------------------------------------------
   ! Initialize arrays.
   ! ---------------------------------------------------------------------------

      bfsqi (:,:,:) = spval
      bfsqf (:,:,:) = spval
      z     (:,:,:) = spval
      bfsql (:,:,:) = spval
      nslpx (:,:,:) = spval
      nslpy (:,:,:) = spval
      nnslpx(:,:,:) = spval
      nnslpy(:,:,:) = spval
      dz    (:,:,:) = spval
      mlts  (:,:) = spval

   end subroutine inivar_cmnfld

end module mod_cmnfld
