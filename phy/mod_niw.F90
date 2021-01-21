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

module mod_niw
! ------------------------------------------------------------------------------
! This module contains variables and procedures related near-intertial waves.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_constants, only: g, alpha0, pi, spval
   use mod_time, only: delt1, dlt
   use mod_grid, only: scuy, scvx, coriop
   use mod_xc
   use mod_state, only: u, v, dpu, dpv, pbu, pbv, ubflxs_p, vbflxs_p
   use mod_utility, only: util1, util2

   implicit none

   private

   ! Variables to be set by namelist.
   real(r8) :: &
      niwgf, &   ! Global factor applied to the energy input by near-inertial
                 ! waves [].
      niwbf, &   ! Fraction of near-inertial energy dissipated in the boundary
                 ! layer [].
      niwlf      ! Fraction of near-inertial energy dissipated locally beneath
                 ! the boundary layer [].

   ! Parameters:
   real(r8), parameter :: &
      ipfac = 2._r8, &      ! Inertial period factor for the velocity averaging
                            ! time scale [].
      cori10 = 2.5256e-5_r8 ! Coriolis parameter at 10N [1/s].

   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, 4) :: &
      uml, &    ! u-component of mixed layer velocities [cm s-1].
      vml       ! v-component of mixed layer velocities [cm s-1].

   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, 2) :: &
      umlres, & ! u-component of mixed layer velocity reservoar for temporal
                ! smoothing [cm s-1].
      vmlres    ! v-component of mixed layer velocity reservoar for temporal
                ! smoothing [cm s-1].

   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy) :: &
      idkedt  ! Vertically integrated inertial kinetic energy tendency
              ! [cm3 s-3].


   public :: niwgf, niwbf, niwlf, &
             uml, vml, umlres, vmlres, idkedt, &
             inivar_niw, niw_ke_tendency

contains

   subroutine inivar_niw
   ! ---------------------------------------------------------------------------
   ! Initialize arrays.
   ! ---------------------------------------------------------------------------

      integer :: i, j, k, l

   !$omp parallel do private(i, k)
      do j = 1 - nbdy, jj + nbdy
         do k = 1, 4
            do i = 1 - nbdy, ii + nbdy
               uml(i, j, k) = spval
               vml(i, j, k) = spval
            enddo
         enddo
         do k = 1, 2
            do i = 1 - nbdy, ii + nbdy
               umlres(i, j, k) = spval
               vmlres(i, j, k) = spval
            enddo
         enddo
         do i = 1 - nbdy, ii + nbdy
            idkedt(i, j) = spval
         enddo
      enddo
   !$omp end parallel do

   !$omp parallel do private(k, l, i)
      do j = 1, jj
        do k = 1, 4
          do l = 1, isu(j)
          do i = max(1, ifu(j, l)), min(ii, ilu(j, l))
            uml(i, j, k) = 0._r8
          enddo
          enddo
          do l = 1, isv(j)
          do i = max(1, ifv(j, l)), min(ii, ilv(j, l))
            vml(i, j, k) = 0._r8
          enddo
          enddo
        enddo
        do k = 1, 2
          do l = 1, isu(j)
          do i = max(1, ifu(j, l)), min(ii+1, ilu(j, l))
            umlres(i, j, k) = 0._r8
          enddo
          enddo
          do l = 1, isv(j)
          do i = max(1, ifv(j, l)), min(ii, ilv(j, l))
            vmlres(i, j, k) = 0._r8
          enddo
          enddo
        enddo
      enddo
   !$omp end parallel do

   end subroutine inivar_niw

   subroutine niw_ke_tendency(m, n, mm, nn, k1m, k1n)
   ! ---------------------------------------------------------------------------
   ! Estimate energy input by near-inertial waves by considering the mixed layer
   ! kinetic energy change over a model time step of velocity deviations from
   ! time averaged velocities. The velocity averaging time scale is the inertial
   ! period times 'ipfac' with the inertial period limited to the one at 10N/S.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      integer :: mmm, i, j, l
      real(r8) :: ubt, uml1t, uml2t, q, uml1a, uml2a, vbt, vml1t, vml2t, &
                  vml1a, vml2a

      mmm = (m - 1)*2

   !$omp parallel do private(l, i, ubt, uml1t, uml2t, q, uml1a, uml2a)
      do j = 1, jj
         do l = 1, isu(j)
         do i = max(1, ifu(j, l)), min(ii + 1, ilu(j, l))

            ubt = ubflxs_p(i, j, m)*dlt/(delt1*scuy(i, j)*pbu(i, j, m))
            uml1t = u(i, j, 1 + mm) + ubt
            uml2t = u(i, j, 2 + mm) + ubt

            q = delt1*max(cori10, &
                          abs(.5_r8*(coriop(i - 1, j) + coriop(i, j)))) &
                /(ipfac*2._r8*pi)
            umlres(i, j, 1) = umlres(i, j, 1) + uml1t
            uml1a = umlres(i, j, 1)*q
            umlres(i, j, 1) = umlres(i, j, 1)*(1._r8 - q)
            umlres(i, j, 2) = umlres(i, j, 2) + uml2t
            uml2a = umlres(i, j, 2)*q
            umlres(i, j, 2) = umlres(i, j, 2)*(1._r8 - q)

            util1(i, j) = ( (uml1t              - uml1a)**2 &
                          - (uml(i, j, 1 + mmm) - uml1a)**2) &
                          *dpu(i, j, 1 + mm) &
                        + ( (uml2t              - uml2a)**2 &
                          - (uml(i, j, 2 + mmm) - uml2a)**2) &
                          *dpu(i, j, 2 + mm)

            uml(i, j, 1 + mmm) = uml1t
            uml(i, j, 2 + mmm) = uml2t

         enddo
         enddo
      enddo
   !$omp end parallel do
   !$omp parallel do private(l, i, vbt, vml1t, vml2t, q, vml1a, vml2a)
      do j = 1, jj + 1
         do l = 1, isv(j)
         do i = max(1, ifv(j, l)), min(ii, ilv(j, l))

            vbt = vbflxs_p(i, j, m)*dlt/(delt1*scvx(i, j)*pbv(i, j, m))
            vml1t = v(i, j, 1 + mm) + vbt
            vml2t = v(i, j, 2 + mm) + vbt

            q = delt1*max(cori10, &
                          abs(.5_r8*(coriop(i, j - 1) + coriop(i, j)))) &
                /(ipfac*2._r8*pi)
            vmlres(i, j, 1) = vmlres(i, j, 1) + vml1t
            vml1a = vmlres(i, j, 1)*q
            vmlres(i, j, 1) = vmlres(i, j, 1)*(1._r8 - q)
            vmlres(i, j, 2) = vmlres(i, j, 2) + vml2t
            vml2a = vmlres(i, j, 2)*q
            vmlres(i, j, 2) = vmlres(i, j, 2)*(1._r8 - q)

            util2(i, j) = ( (vml1t              - vml1a)**2 &
                          - (vml(i, j, 1 + mmm) - vml1a)**2) &
                          *dpv(i, j, 1 + mm) &
                        + ( (vml2t              - vml2a)**2 &
                          - (vml(i, j, 2 + mmm) - vml2a)**2) &
                          *dpv(i, j, 2 + mm)

            vml(i, j, 1 + mmm) = vml1t
            vml(i, j, 2 + mmm) = vml2t

         enddo
         enddo
      enddo
   !$omp end parallel do
   !$omp parallel do private(l, i)
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
            idkedt(i, j) = abs(  ( util1(i    , j)*iu(i    , j) &
                                 + util1(i + 1, j)*iu(i + 1, j)) &
                                 /max(1, iu(i, j) + iu(i + 1, j)) &
                               + ( util2(i, j    )*iv(i, j    ) &
                                 + util2(i, j + 1)*iv(i, j + 1)) &
                                 /max(1, iv(i, j) + iv(i, j + 1))) &
                           *alpha0/(2._r8*g*delt1)
         enddo
         enddo
      enddo
   !$omp end parallel do

   end subroutine niw_ke_tendency

end module mod_niw
