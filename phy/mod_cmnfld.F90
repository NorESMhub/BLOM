! ------------------------------------------------------------------------------
! Copyright (C) 2015-2022 Mats Bentsen, Mehmet Ilicak
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
   use mod_constants, only: spval, onem, L_mks2cgs
   use mod_xc

   implicit none

   private

   ! Parameters:
   real(r8) :: &
      sls0 = 10._r8*onem, &      ! Minimum smoothing length scale in the
                                 ! computation of filtered BFSQ [g cm-1 s-2].
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
      dbcrit = .0003_r8*L_mks2cgs! Critical buoyancy difference used in the
                                 ! mixed layer thickness estimation (Levitus,
                                 ! 1982) [cm s-2].

   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, kdm + 1) :: &
      bfsqi, &                   ! Interface buoyancy frequency squared [s-2].
      bfsqf, &                   ! Filtered interface buoyancy frequency
                                 ! squared [s-2].
      z                          ! Interface depth [cm].
   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, kdm) :: &
      bfsql, &                   ! Layer buoyancy frequency squared [s-2].
      nslpx, &                   ! x-component of local neutral slope [].
      nslpy, &                   ! y-component of local neutral slope [].
      nnslpx, &                  ! x-component of local neutral slope times
                                 ! buoyancy frequency [s-1].
      nnslpy, &                  ! y-component of local neutral slope times
                                 ! buoyancy frequency [s-1].
      dz                         ! Layer thickness [cm].
   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy) :: &
      mlts                       ! Mixed layer depth defined by density
                                 ! criterion [cm].

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

      integer :: i,j,k

      !$omp parallel do private(k, i)
      do j = 1 - nbdy, jj + nbdy
         do k = 1, kk + 1
            do i = 1 - nbdy, ii + nbdy
               bfsqi(i, j, k) = spval
               bfsqf(i, j, k) = spval
               z    (i, j, k) = spval
            enddo
         enddo
         do k = 1, kk
            do i = 1 - nbdy, ii + nbdy
               bfsql (i, j, k) = spval
               nslpx (i, j, k) = spval
               nslpy (i, j, k) = spval
               nnslpx(i, j, k) = spval
               nnslpy(i, j, k) = spval
               dz    (i, j, k) = spval
            enddo
         enddo
         do i = 1 - nbdy, ii + nbdy
            mlts(i, j) = spval
         enddo
      enddo
      !$omp end parallel do

   end subroutine inivar_cmnfld

end module mod_cmnfld
