! ------------------------------------------------------------------------------
! Copyright (C) 2020-2025 Mats Bentsen
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

module mod_utility
! ------------------------------------------------------------------------------
! This module contains various utility variables.
! ------------------------------------------------------------------------------

   use mod_types, only: r8, blom_char_x
   use mod_constants, only: spval
   use mod_xc
   use mod_checksum, only: csdiag, chksummsk

   implicit none

   private

   ! Total (barotropic + baroclinic) velocity components at two time levels.
   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy) :: &
      utotm, vtotm, utotn, vtotn

   ! Horizontal mass flux components.
   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy) :: &
      uflux, vflux, uflux2, vflux2, uflux3, vflux3 

   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy) :: &
      umax, & ! u-component of maximum allowable velocity ensuring stability of
              ! the upwind scheme [m s-1].
      vmax    ! v-component of maximum allowable velocity ensuring stability of
              ! the upwind scheme [m s-1].

   ! Arrays for temporary storage.
   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy) :: &
      util1, util2, util3, util4

   ! Default length of character variables for filenames
   integer,parameter :: fnmlen = blom_char_x

   public :: utotm, vtotm, utotn, vtotn, &
             uflux, vflux, uflux2, vflux2, uflux3, vflux3, &
             umax, vmax, &
             util1, util2, util3, util4, &
             inivar_utility, fnmlen
   
contains

   subroutine inivar_utility
   ! ---------------------------------------------------------------------------
   ! Initialize state variables.
   ! ---------------------------------------------------------------------------

      integer :: i, j, l

   !$omp parallel do private(i)
      do j = 1 - nbdy, jj + nbdy
         do i = 1 - nbdy, ii + nbdy
            utotm (i, j) = spval
            vtotm (i, j) = spval
            utotn (i, j) = spval
            vtotn (i, j) = spval
            uflux (i, j) = spval
            vflux (i, j) = spval
            uflux2(i, j) = spval
            vflux2(i, j) = spval
            uflux3(i, j) = spval
            vflux3(i, j) = spval
            umax  (i, j) = spval
            vmax  (i, j) = spval
            util1 (i, j) = spval
            util2 (i, j) = spval
            util3 (i, j) = spval
            util4 (i, j) = spval
         enddo
      enddo
   !$omp end parallel do

      ! Initialize 'utotm', 'uflux', 'uflux2', 'uflux3' at points located
      ! upstream and downstream (in i-direction) of p-points. Initialize 'utotn'
      ! upstream and downstream of p-points as well as at lateral neighbors of
      ! interior u-points.
   !$omp parallel do private(l, i)
      do j = 0, jj + 1
         do l = 1, isu(j)
         do i = max(1, ifu(j, l)), min(ii, ilu(j, l))
            utotn(i, j - 1) = 0._r8
            utotn(i, j + 1) = 0._r8
         enddo
         enddo
      enddo
   !$omp end parallel do
   !$omp parallel do private(l, i)
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j, l)), min(ii, ilp(j, l) + 1)
            utotn (i, j) = 0._r8
            utotm (i, j) = 0._r8
            uflux (i, j) = 0._r8
            uflux2(i, j) = 0._r8
            uflux3(i, j) = 0._r8
         enddo
         enddo
      enddo
   !$omp end parallel do
      call xctilr(utotn,  1, 1, nbdy, nbdy, halo_us)
      call xctilr(utotm,  1, 1, nbdy, nbdy, halo_us)
      call xctilr(uflux,  1, 1, nbdy, nbdy, halo_us)
      call xctilr(uflux2, 1, 1, nbdy, nbdy, halo_us)
      call xctilr(uflux3, 1, 1, nbdy, nbdy, halo_us)

      ! Initialize 'vtotm', 'vflux', 'vflux2', 'vflux3' at points located
      ! upstream and downstream (in j-direction) of p-points. Initialize 'vtotn'
      ! upstream and downstream of p-points as well as at lateral neighbors of
      ! interior v-points.
   !$omp parallel do private(l, j)
      do i = 0, ii + 1
         do l = 1, jsv(i)
         do j = max(1, jfv(i, l)), min(jj, jlv(i, l))
            vtotn(i - 1, j) = 0._r8
            vtotn(i + 1, j) = 0._r8
         enddo
         enddo
      enddo
   !$omp end parallel do
   !$omp parallel do private(l, j)
      do i = 1, ii
         do l = 1, jsp(i)
         do j = max(1, jfp(i, l)), min(jj, jlp(i, l) + 1)
            vtotn (i, j) = 0._r8
            vtotm (i, j) = 0._r8
            vflux (i, j) = 0._r8
            vflux2(i, j) = 0._r8
            vflux3(i, j) = 0._r8
         enddo
         enddo
      enddo
   !$omp end parallel do
      call xctilr(vtotn,  1, 1, nbdy, nbdy, halo_vs)
      call xctilr(vtotm,  1, 1, nbdy, nbdy, halo_vs)
      call xctilr(vflux,  1, 1, nbdy, nbdy, halo_vs)
      call xctilr(vflux2, 1, 1, nbdy, nbdy, halo_vs)
      call xctilr(vflux3, 1, 1, nbdy, nbdy, halo_vs)

      if (csdiag) then
         if (mnproc == 1) then
            write (lp, *) 'inivar_utility:'
         endif
!        call chksummsk(utotm, iu, 1, 'utotm')
!        call chksummsk(vtotm, iv, 1, 'vtotm')
!        call chksummsk(utotn, iu, 1, 'utotn')
!        call chksummsk(vtotn, iv, 1, 'vtotn')
!        call chksummsk(uflux, iu, 1, 'uflux')
!        call chksummsk(vflux, iv, 1, 'vflux')
!        call chksummsk(uflux2, iu, 1, 'uflux2')
!        call chksummsk(vflux2, iv, 1, 'vflux2')
!        call chksummsk(uflux3, iu, 1, 'uflux3')
!        call chksummsk(vflux3, iv, 1, 'vflux3')
      endif

   end subroutine inivar_utility

end module mod_utility
