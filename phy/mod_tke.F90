! ------------------------------------------------------------------------------
! Copyright (C) 2013-2020 Mehmet Ilicak, Mats Bentsen
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

module mod_tke
! ------------------------------------------------------------------------------
! This module contains constants, variables and routines related to a k-epsilon
! model using a second-order turbulence closure.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_xc
   use mod_diffusion, only: difdia
   use mod_forcing, only: ustarb

   implicit none

   private

   real(r8), parameter :: &
      gls_cmu0 = .527_r8, &      ! cmu0
      Pr_t = 1._r8, &            ! Turbulent Prandtl number [non-dimensional].
      zos = .0002_r8, &          !
      gls_p = 3._r8, &           !
      gls_m = 1.5_r8, &          !
      gls_n = -1._r8, &          !
      gls_c1 = 1.44_r8, &        !
      gls_c2 = 1.92_r8, &        !
      gls_c3plus = 1._r8, &      !
      gls_c3minus = -.63_r8, &   !
      gls_L1 = .107_r8, &        !
      gls_L2 = .0032_r8, &       !
      gls_L3 = .0864_r8, &       !
      gls_L4 = .12_r8, &         !
      gls_L5 = 11.9_r8, &        !
      gls_L6 = .4_r8, &          !
      gls_L7 = .0_r8, &          !
      gls_L8 = .48_r8, &         !
      gls_Gh0 = .0329_r8, &      !
      gls_Ghmin = -.28_r8, &     !
      gls_Ghcri = .03_r8, &      !
      vonKar = .4_r8             !

#if defined(CGS)
   real(r8), parameter :: &
      tke_min = 7.6e-4_r8, &     ! Minimum TKE value [cm2/s2].
      gls_psi_min = 1.e-10_r8, & ! Minimum GLS value [cm2/s3].
      Ls_unlmt_min = 1.e-6_r8    ! [cm]
#endif
#if defined(MKS)
   real(r8), parameter :: &
      tke_min = 7.6e-8_r8, &     ! Minimum TKE value [m2/s2].
      gls_psi_min = 1.e-14_r8, & ! Minimum GLS value [m2/s3].
      Ls_unlmt_min = 1.e-8_r8    ! [m]
#endif

   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, kdm) :: &
      Prod, &                    ! Shear production [?].
      Buoy, &                    ! Buoyancy production [?].
      Shear2, &                  ! Square of the shear frequency [?].
      L_scale                    ! Dissipative length scale [?].

   ! Various coefficients.
   real(r8) :: gls_s0, gls_s1, gls_s2, gls_s4, gls_s5, gls_s6, &
               gls_b0, gls_b1, gls_b2, gls_b3, gls_b4, gls_b5, &
               sqrt2, cmu_fac1, cmu_fac2, cmu_fac3, tke_exp1, gls_exp1, gls_fac6

   public :: gls_cmu0, Pr_t, tke_min, zos, gls_psi_min, gls_p, gls_m, gls_n, &
             gls_c1, gls_c2, gls_c3plus, gls_c3minus, gls_Gh0, gls_Ghmin, &
             gls_Ghcri, vonKar, Ls_unlmt_min, Prod, Buoy, Shear2, L_scale, &
             gls_s0, gls_s1, gls_s2, gls_s4, gls_s5, gls_s6, gls_b0, gls_b1, &
             gls_b2, gls_b3, gls_b4, gls_b5, sqrt2, cmu_fac1, cmu_fac2, &
             cmu_fac3, tke_exp1, gls_exp1, gls_fac6, &
             initke

contains

   subroutine initke
   ! ---------------------------------------------------------------------------
   ! Initialization of second order turbulence closure.
   ! ---------------------------------------------------------------------------

#if defined(TRC) && defined(TKE)
      use mod_xc
      use mod_tracers, only: itrtke, itrgls, trc

      integer :: i, j, k, l

      ! Initialize fields holding turbulent kinetic energy, generic length
      ! scale, and other fields used in the turbulence closure.
   !$omp parallel do private(k, l, i)
      do j = 1 - nbdy, jj + nbdy
         do k = 1, 2*kdm
            do l = 1, isp(j)
            do i = ifp(j, l), ilp(j, l)
               trc(i, j, k, itrtke) = tke_min
               trc(i, j, k, itrgls) = gls_psi_min
            enddo
            enddo
         enddo
         do k = 1, kk
            do l = 1, isp(j)
            do i = ifp(j, l), ilp(j, l)
               difdia(i, j, k) = 0._r8
               L_scale(i, j, k) = Ls_unlmt_min
            enddo
            enddo
         enddo
         do l = 1, isp(j)
         do i = ifp(j, l), ilp(j, l)
            ustarb(i, j) = 0._r8
         enddo
         enddo
      enddo
   !$omp end parallel do

      ! Precompute various coefficients
      sqrt2 = sqrt(2._r8)
      cmu_fac1 = gls_cmu0**(- gls_p/gls_n)
      cmu_fac2 = gls_cmu0**(3._r8 + gls_p/gls_n)
      cmu_fac3 = sqrt2
      tke_exp1 = gls_m/gls_n
      gls_exp1 = 1._r8/gls_n
      gls_fac6 = 8._r8/gls_cmu0**6
      gls_s0 = 1.5_r8*gls_L1*gls_L5**2
      gls_s1 = - gls_L4*(gls_L6 + gls_L7) &
               + 2._r8*gls_L4*gls_L5*(gls_L1 - 1._r8/3._r8*gls_L2 - gls_L3) &
               + 1.5_r8*gls_L1*gls_L5*gls_L8
      gls_s2 = - 3._r8/8._r8*gls_L1*(gls_L6**2 - gls_L7**2)
      gls_s4 = 2._r8*gls_L5
      gls_s5 = 2._r8*gls_L4
      gls_s6 = 2._r8/3._r8*gls_L5*(3._r8*gls_L3**2 - gls_L2**2) &
             - .5_r8*gls_L5*gls_L1*(3._r8*gls_L3 - gls_L2) &
             + .75_r8*gls_L1*(gls_L6 - gls_L7)
      gls_b0 = 3._r8*gls_L5**2
      gls_b1 = gls_L5*(7._r8*gls_L4 + 3._r8*gls_L8)
      gls_b2 = gls_L5**2*(3._r8*gls_L3**2 - gls_L2**2) &
             - .75_r8*(gls_L6**2 - gls_L7**2)
      gls_b3 = gls_L4*(4._r8*gls_L4 + 3._r8*gls_L8)
      gls_b4 = gls_L4*(gls_L2*gls_L6 - 3._r8*gls_L3*gls_L7 &
             - gls_L5*(gls_L2**2 - gls_L3**2)) &
             + gls_L5*gls_L8*(3._r8*gls_L3**2 - gls_L2**2)
      gls_b5 = .25_r8*(gls_L2**2 - 3._r8*gls_L3**2)*(gls_L6**2 - gls_L7**2)

#endif

      end subroutine initke

end module mod_tke
