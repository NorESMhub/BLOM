! ------------------------------------------------------------------------------
! Copyright (C) 2015-2020 Mats Bentsen
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
   use mod_constants, only: g, alpha0, epsil, spval, onem, onemm
   use mod_xc
   use mod_vcoord, only: vcoord_type_tag, isopyc_bulkml, cntiso_hybrid
   use mod_grid, only: scuxi, scvyi
   use mod_eos, only: rho, p_alpha
   use mod_state, only: dp, temp, saln, p, phi, kfpla
!  use mod_dia, only : nphy, ACC_BFSQ
   use mod_diffusion, only: eitmth, edritp 
   use mod_utility, only: util1
   use mod_checksum, only: csdiag, chksummsk

   implicit none

   private

   ! Parameters:
   real(r8) :: &
      sls0 = 10._r8*98060._r8, & ! Minimum smoothing length scale in the
                                 ! computation of filtered BFSQ [g cm-1 s-2].
      slsmfq = 2._r8, &          ! Factor to be multiplied with the mixed
                                 ! layer depth to find the smoothing length
                                 ! scale at the base of the mixed layer in the
                                 ! computation of filtered BFSQ [].
      slsels = 2._r8, &          ! Factor to be multiplied with the mixed
                                 ! layer depth to find the e-folding length
                                 ! scale of the smoothing length scale in the
                                 ! computation of filtered BFSQ [].
      bfsqmn = 1.e-7_r8          ! Minimum value of BFSQ used in the
                                 ! computation of neutral slope [s-2].

   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, kdm + 1) :: &
      bfsqi, &                   ! Interface buoyancy frequency squared [s-2].
      bfsqf                      ! Filtered interface buoyancy frequency
                                 ! squared [s-2].
   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, kdm) :: &
      bfsql, &                   ! Layer buoyancy frequency squared [s-2].
      nslpx, &                   ! x-component of local neutral slope [].
      nslpy, &                   ! y-component of local neutral slope [].
      nnslpx, &                  ! x-component of local neutral slope times
                                 ! buoyancy frequency [s-1].
      nnslpy                     ! y-component of local neutral slope times
                                 ! buoyancy frequency [s-1].

   public :: bfsql, nslpx, nslpy, nnslpx, nnslpy, inivar_cmnfld, cmnfld
   public :: bfsqi

   contains

   ! ---------------------------------------------------------------------------
   ! Private procedures.
   ! ---------------------------------------------------------------------------

   subroutine cmnfld_bfsqf_isopyc_bulkml(m, n, mm, nn, k1m, k1n)
   ! ---------------------------------------------------------------------------
   ! Compute buoyancy frequency squared (BFSQ) on layer interfaces and
   ! representative of the layer itself. Also compute a filtered BFSQ on
   ! interfaces.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      real(r8), dimension(kdm) :: delp, bfsq, sls2, atd, btd, ctd, rtd, gam
      real(r8) :: pml, q, pup, tup, sup, plo, tlo, slo, bei
      integer :: i, j, k, l, km, kfpl

      ! ------------------------------------------------------------------------
      ! The BFSQ is estimated locally at layer interfaces. The filtered BFSQ is
      ! smoothed in the vertical direction by solving a diffusion equation. At
      ! the mixed layer base the diffusion length scale is set to the maximum of
      ! sls0 and mixed layer depth (MLD) times slsmfq. Below the mixed layer,
      ! the diffusion length scale approaches sls0 with an e-folding length
      ! scale of MLD times slsels.
      ! ------------------------------------------------------------------------

   !$omp parallel do private(l, i, kfpl, k, pml, delp, bfsq, q, sls2, &
   !$omp                     pup, tup, sup, km, plo, tlo, slo, &
   !$omp                     ctd, btd, rtd, atd, bei, gam)
      do j = - 1, jj + 2
         do l = 1, isp(j)
         do i = max(- 1, ifp(j, l)), min(ii + 2, ilp(j, l))

            ! Compute BFSQ in the mixed layer.
            bfsqi(i, j, 1) = &
               .5_r8*g*g*( rho(p(i, j, 2), &
                               temp(i, j, 2 + mm), saln(i, j, 2 + mm)) &
                         - rho(p(i, j, 2), &
                               temp(i, j, 1 + mm), saln(i, j, 1 + mm))) &
               /(dp(i, j, 1 + mm) + dp(i, j, 2 + mm))
            bfsqi(i, j, 2) = bfsqi(i, j, 1)
            bfsql(i, j, 1) = bfsqi(i, j, 1)
            bfsql(i, j, 2) = bfsqi(i, j, 1)

            kfpl = kfpla(i, j, m)

            if (kfpl > kk) then

               ! If the mixed layer extends to the bottom, propagate the
               ! interface and layer BFSQ of the mixed layer downwards while the
               ! filtered BFSQ is set to a minimum value.
               do k = 3, kk
                  bfsqi(i, j, k) = bfsqi(i, j, 1)
                  bfsql(i, j, k) = bfsqi(i, j, 1)
               enddo
               bfsqi(i, j, kk + 1) = bfsqi(i, j, 1)
               do k = 1, kk + 1
                  bfsqf(i, j, k) = bfsqmn
               enddo

            else

               ! At layer interfaces, compute BFSQ and length scale for the
               ! subsequent smoothing.
               pml = max(.5_r8*(p(i, j, 3) + p(i, j, 1)), &
                         .5_r8*(3._r8*p(i, j, 3) - p(i, j, kfpl + 1)))
               delp(kfpl - 1) = pml - p(i, j, 1)
               bfsqi(i, j, kfpl - 1) = bfsqi(i, j, 2)
               bfsq(kfpl - 1) = bfsqmn
               q = max(sls0, delp(kfpl - 1)*slsmfq)
               sls2(kfpl - 1) = q*q
               pup = pml
               tup = temp(i, j, 2 + mm)
               sup = saln(i, j, 2 + mm)
               do k = kfpl, kk
                  km = k + mm
                  if (p(i, j, kk + 1) - p(i, j, k) < epsil) then
                     delp(k) = onemm
                     bfsqi(i, j, k) = bfsqi(i, j, k - 1)
                     bfsq(k) = bfsqmn
                     q = exp(- (p(i, j, kk + 1) - pml)/(slsels*delp(kfpl - 1)))
                     q = max(sls0, delp(kfpl - 1)*slsmfq*q + sls0*(1._r8 - q))
                     sls2(k) = q*q
                  else
                     if (p(i, j, kk + 1) - p(i, j, k + 1) < epsil) then
                        plo = p(i, j, kk + 1)
                     else
                        plo = .5_r8*(p(i, j, k) + p(i, j, k + 1))
                     endif
                     tlo = temp(i, j, km)
                     slo = saln(i, j, km)
                     delp(k) = max(onemm, plo - pup)
                     bfsqi(i, j, k) = g*g*( rho(p(i, j, k), tlo, slo) &
                                          - rho(p(i, j, k), tup, sup))/delp(k)
                     bfsq(k) = max(bfsqmn, bfsqi(i, j, k))
                     bfsqi(i, j, k) = bfsqi(i, j, k)*delp(k)/max(onem, delp(k))
                     if (p(i, j, kk + 1) - p(i, j, k) < onem) then
                        bfsqi(i, j, k) = bfsqi(i, j, k - 1)
                     endif
                     q = exp(- (p(i, j, k) - pml)/(slsels*delp(kfpl - 1)))
                     q = max(sls0, delp(kfpl - 1)*slsmfq*q + sls0*(1._r8 - q))
                     sls2(k) = q*q
                     pup = plo
                     tup = tlo
                     sup = slo
                  endif
               enddo

               ! Compute the layer BFSQ as the arithmetic mean of the layer
               ! interface BFSQ.
               do k = kfpl, kk - 1
                  bfsql(i, j, k) = .5_r8*(bfsqi(i, j, k) + bfsqi(i, j, k + 1))
               enddo
               bfsql(i, j, kk) = bfsqi(i, j, kk)
               do k = 3, kfpl - 1
                  bfsqi(i, j, k) = bfsqi(i, j, kfpl)
                  bfsql(i, j, k) = bfsql(i, j, kfpl)
               enddo

               ! For the filtered BFSQ, compute the coefficients for the
               ! tridiagonal set of equations arising from the implicit backward
               ! discretization.
               k = kfpl - 1
               ctd(k) = - 2._r8*sls2(k    ) &
                          /(delp(k)*(delp(k    ) + delp(k + 1)))
               btd(k) = 1._r8 - ctd(k)
               rtd(k) = bfsq(k)
               do k = kfpl, kk - 1
                  atd(k) = - 2._r8*sls2(k - 1) &
                             /(delp(k)*(delp(k - 1) + delp(k    )))
                  ctd(k) = - 2._r8*sls2(k    ) &
                             /(delp(k)*(delp(k    ) + delp(k + 1)))
                  btd(k) = 1._r8 - atd(k) - ctd(k)
                  rtd(k) = bfsq(k)
               enddo
               k = kk
               atd(k) = - 2._r8*sls2(k - 1) &
                          /(delp(k)*(delp(k - 1) + delp(k    )))
               btd(k) = 1._r8 - atd(k)
               rtd(k) = bfsq(k)

               ! Solve the tridiagonal set of equations.
               bei = 1._r8/btd(kfpl - 1)
               bfsqf(i, j, kfpl - 1) = rtd(kfpl - 1)*bei
               do k = kfpl, kk
                  gam(k) = ctd(k - 1)*bei
                  bei = 1._r8/(btd(k) - atd(k)*gam(k))
                  bfsqf(i, j, k) = (rtd(k) - atd(k)*bfsqf(i, j, k - 1))*bei
               enddo
               do k = kk - 1, kfpl - 1, - 1
                  bfsqf(i, j, k) = bfsqf(i, j, k) &
                                 - gam(k + 1)*bfsqf(i, j, k + 1)
               enddo
               do k = 1, kfpl - 2
                  bfsqf(i, j, k) = bfsqf(i, j, kfpl - 1)
               enddo

               ! Extrapolate to the bottom interface.
               bfsqi(i, j, kk + 1) = bfsqi(i, j, kk)
               bfsqf(i, j, kk + 1) = bfsqf(i, j, kk)

            endif

         enddo
         enddo
      enddo
   !$omp end parallel do

      if (csdiag) then
         if (mnproc == 1) then
            write(lp,*) 'cmnfld_bfsqf_isopyc_bulkml:'
         endif
         call chksummsk(bfsqi, ip, kk + 1, 'bfsqi')
         call chksummsk(bfsql, ip, kk, 'bfsql')
         call chksummsk(bfsqf, ip, kk + 1, 'bfsqf')
      endif

   end subroutine cmnfld_bfsqf_isopyc_bulkml

   subroutine cmnfld_bfsqf_cntiso_hybrid(m, n, mm, nn, k1m, k1n)
   ! ---------------------------------------------------------------------------
   ! Compute buoyancy frequency squared (BFSQ) on layer interfaces and
   ! representative of the layer itself. Also compute a filtered BFSQ on
   ! interfaces.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      real(r8), dimension(kdm) :: delp, bfsq, sls2, atd, btd, ctd, rtd, gam
      real(r8) :: q, pup, tup, sup, plo, tlo, slo, bei
      integer :: i, j, k, l, km

      ! ------------------------------------------------------------------------
      ! The BFSQ is estimated locally at layer interfaces. The filtered BFSQ is
      ! smoothed in the vertical direction by solving a diffusion equation. 
      ! ------------------------------------------------------------------------

      bfsqi = 0.0_r8
      bfsql = 0.0_r8
   !$omp parallel do private(l, i, k, delp, bfsq, q, sls2, pup, tup, sup, km, &
   !$omp                     plo, tlo, slo, ctd, btd, rtd, atd, bei, gam)
      do j = - 1, jj + 2
         do l = 1, isp(j)
         do i = max(- 1, ifp(j, l)), min(ii + 2, ilp(j, l))

            ! At layer interfaces, compute BFSQ and length scale for the
            ! subsequent smoothing.
            bfsqi(i, j, 1) = bfsqmn
            pup = .5_r8*(p(i, j, 1) + p(i, j, 2))
            tup = temp(i, j, 1 + mm)
            sup = saln(i, j, 1 + mm)
            do k = 2, kk
               km = k + mm
               if (p(i, j, kk + 1) - p(i, j, k) < epsil) then
                  delp(k) = onemm
                  bfsqi(i, j, k) = bfsqi(i, j, k - 1)
                  bfsq(k) = bfsqmn
                  sls2(k) = sls0*sls0
               else
                  if (p(i, j, kk + 1) - p(i, j, k + 1) < epsil) then
                     plo = p(i, j, kk + 1)
                  else
                     plo = .5_r8*(p(i, j, k) + p(i, j, k + 1))
                  endif
                  tlo = temp(i, j, km)
                  slo = saln(i, j, km)
                  delp(k) = max(onemm, plo - pup)
                  bfsqi(i, j, k) = g*g*( rho(p(i, j, k), tlo, slo) &
                                       - rho(p(i, j, k), tup, sup))/delp(k)
                  bfsq(k) = max(bfsqmn, bfsqi(i, j, k))
                  bfsqi(i, j, k) = bfsqi(i, j, k)*delp(k)/max(onem, delp(k))
                  if (p(i, j, kk + 1) - p(i, j, k) < onem) then
                     bfsqi(i, j, k) = bfsqi(i, j, k - 1)
                  endif
                  sls2(k) = sls0*sls0
                  pup = plo
                  tup = tlo
                  sup = slo
               endif
            enddo
            delp(1) = dp(i, j, 1 + mm)
            bfsqi(i, j, 1) = bfsqi(i, j, 2)
            bfsq(1) = max(bfsqmn, bfsqi(i, j, 1))
            sls2(1) = sls0*sls0

            ! Compute the layer BFSQ as the arithmetic mean of the layer
            ! interface BFSQ.
            do k = 1, kk - 1
               bfsql(i, j, k) = .5_r8*(bfsqi(i, j, k) + bfsqi(i, j, k + 1))
            enddo
            bfsql(i, j, kk) = bfsqi(i, j, kk)

            ! For the filtered BFSQ, compute the coefficients for the
            ! tridiagonal set of equations arising from the implicit backward
            ! discretization.
            k = 1
            ctd(k) = - 2._r8*sls2(k    ) &
                       /(delp(k)*(delp(k    ) + delp(k + 1)))
            btd(k) = 1._r8 - ctd(k)
            rtd(k) = bfsq(k)
            do k = 2, kk - 1
               atd(k) = - 2._r8*sls2(k - 1) &
                          /(delp(k)*(delp(k - 1) + delp(k    )))
               ctd(k) = - 2._r8*sls2(k    ) &
                          /(delp(k)*(delp(k    ) + delp(k + 1)))
               btd(k) = 1._r8 - atd(k) - ctd(k)
               rtd(k) = bfsq(k)
            enddo
            k = kk
            atd(k) = - 2._r8*sls2(k - 1) &
                       /(delp(k)*(delp(k - 1) + delp(k    )))
            btd(k) = 1._r8 - atd(k)
            rtd(k) = bfsq(k)

            ! Solve the tridiagonal set of equations.
            bei = 1._r8/btd(1)
            bfsqf(i, j, 1) = rtd(1)*bei
            do k = 2, kk
               gam(k) = ctd(k - 1)*bei
               bei = 1._r8/(btd(k) - atd(k)*gam(k))
               bfsqf(i, j, k) = (rtd(k) - atd(k)*bfsqf(i, j, k - 1))*bei
            enddo
            do k = kk - 1, 1, - 1
               bfsqf(i, j, k) = bfsqf(i, j, k) - gam(k + 1)*bfsqf(i, j, k + 1)
            enddo

            ! Extrapolate to the bottom interface.
            bfsqi(i, j, kk + 1) = bfsqi(i, j, kk)
            bfsqf(i, j, kk + 1) = bfsqf(i, j, kk)

         enddo
         enddo
      enddo
   !$omp end parallel do

      if (csdiag) then
         if (mnproc == 1) then
            write(lp,*) 'cmnfld_bfsqf_cntiso_hybrid:'
         endif
         call chksummsk(bfsqi, ip, kk + 1, 'bfsqi')
         call chksummsk(bfsql, ip, kk, 'bfsql')
         call chksummsk(bfsqf, ip, kk + 1, 'bfsqf')
      endif

   end subroutine cmnfld_bfsqf_cntiso_hybrid

   subroutine cmnfld_nslope_isopyc_bulkml(m, n, mm, nn, k1m, k1n)
   ! ---------------------------------------------------------------------------
   ! Estimate slope of local neutral surface.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      real(r8) :: rho0, pm, rho_x, phi_x, bfsqm, rho_y, phi_y
      integer :: i, j, k, l, km, kintr, kmax, knnsl

      ! ------------------------------------------------------------------------
      ! Compute geopotential at layer interfaces.
      ! ------------------------------------------------------------------------

   !$omp parallel do private(k, km, l, i)
      do j = - 1, jj + 2
         do k = kk, 1, - 1
            km = k + mm
            do l = 1, isp(j)
            do i = max(- 1, ifp(j, l)), min(ii + 2, ilp(j, l))
               if (dp(i, j, km) < epsil) then
                  phi(i, j, k) = phi(i, j, k + 1)
               else
                  phi(i, j, k) = phi(i, j, k + 1) &
                               - p_alpha(p(i, j, k + 1), p(i, j, k), &
                                         temp(i, j, km), saln(i, j, km))
               endif
            enddo
            enddo
         enddo
      enddo
   !$omp end parallel do

      ! ------------------------------------------------------------------------
      ! Compute slope vector of local neutral surfaces and also slope vector
      ! times Brunt-Vaisala frequency (optionally used in the computation of
      ! eddy growth rate). The latter is not computed when the gradient of the
      ! geopotential is expected to be influenced by the gradient of the
      ! bathymetry and in this case values are extrapolated from above.
      ! ------------------------------------------------------------------------

      rho0 = 1._r8/alpha0

   !$omp parallel do private(l, i, k, kmax, km, kintr, knnsl, pm, rho_x, &
   !$omp                     phi_x, bfsqm)
      do j = - 1, jj + 2
         do l = 1, isu(j)
         do i = max(0, ifu(j, l)), min(ii + 2, ilu(j, l))

            ! Set the x-component of the slope vector to zero initially.
            do k = 1, kk
               nslpx(i, j, k) = 0._r8
               nnslpx(i, j, k) = 0._r8
            enddo

            if     (kfpla(i - 1, j, m) <= kk .or. kfpla(i, j, m) <= kk) then

               ! Index of last layer containing mass at either of the scalar
               ! points adjacent to the velocity point.
               kmax = 1
               do k = 3, kk
                  km = k + mm
                  if (dp(i - 1, j, km) > epsil .or. dp(i, j, km) > epsil) &
                     kmax = k
               enddo

               ! The first interior interface where the x-component of the slope
               ! vector is estimated is at index kintr + 1.
               kintr = max(kfpla(i - 1, j, m), kfpla(i, j, m))

               ! Index of last interface where slope vector times Brunt-Vaisala
               ! frequency is computed.
               knnsl = 2

               ! Compute the x-component of the slope vector at the mixed layer
               ! base.
               pm = .5_r8*(p(i - 1, j, 3) + p(i, j, 3))
               rho_x = rho(pm, temp(i    , j, 2 + mm), saln(i    , j, 2 + mm)) &
                     - rho(pm, temp(i - 1, j, 2 + mm), saln(i - 1, j, 2 + mm))
               phi_x = phi(i, j, 3) - phi(i - 1, j, 3)
               bfsqm = .5_r8*(bfsqf(i - 1, j, 3) + bfsqf(i, j, 3))
               nslpx(i, j, 3) = (g*rho_x/(rho0*bfsqm) + phi_x/g)*scuxi(i, j)
               if (phi(i    , j, 3) > phi(i - 1, j, kk + 1) .and. &
                   phi(i - 1, j, 3) > phi(i    , j, kk + 1)) then
                  nnslpx(i, j, 3) = sqrt(bfsqm)*nslpx(i, j, 3)
                  knnsl = 3
               endif

               ! Compute the x-component of the slope vector at interior
               ! interfaces.
               do k = kintr + 1, kmax
                  km = k + mm
                  pm = .5_r8*(p(i - 1, j, k) + p(i, j, k))
                  rho_x = .5_r8*( rho(pm, temp(i    , j, km - 1), &
                                          saln(i    , j, km - 1)) &
                                - rho(pm, temp(i - 1, j, km - 1), &
                                          saln(i - 1, j, km - 1)) &
                                + rho(pm, temp(i    , j, km    ), &
                                          saln(i    , j, km    )) &
                                - rho(pm, temp(i - 1, j, km    ), &
                                          saln(i - 1, j, km    )))
                  phi_x = phi(i, j, k) - phi(i - 1, j, k)
                  bfsqm = .5_r8*(bfsqf(i - 1, j, k) + bfsqf(i, j, k))
                  nslpx(i, j, k) = (g*rho_x/(rho0*bfsqm) + phi_x/g)*scuxi(i, j)
                  if (phi(i    , j, k) > phi(i - 1, j, kk + 1) .and. &
                      phi(i - 1, j, k) > phi(i    , j, kk + 1)) then
                     nnslpx(i, j, k) = sqrt(bfsqm)*nslpx(i, j, k)
                     knnsl = k
                  endif
               enddo
               do k = knnsl + 1, kmax
                  nnslpx(i, j, k) = nnslpx(i, j, knnsl)
               enddo
               if (kintr < kmax) then
                  do k = 4, kintr
                     nslpx(i, j, k) = nslpx(i, j, kintr + 1)
                     nnslpx(i, j, k) = nnslpx(i, j, kintr + 1)
                  enddo
               else
                  do k = 4, kmax
                     nslpx(i, j, k) = nslpx(i, j, 3)
                     nnslpx(i, j, k) = nnslpx(i, j, 3)
                  enddo
               endif

            endif

         enddo
         enddo
      enddo
   !$omp end parallel do

   !$omp parallel do private(l, i, k, kmax, km, kintr, knnsl, pm, rho_y, &
   !$omp                     phi_y, bfsqm)
      do j = 0, jj + 2
         do l = 1, isv(j)
         do i = max(- 1, ifv(j, l)), min(ii + 2, ilv(j, l))

            ! Set the y-component of the slope vector to zero initially.
            do k = 1, kk
               nslpy(i, j, k) = 0._r8
               nnslpy(i, j, k) = 0._r8
            enddo

            if     (kfpla(i, j - 1, m) <= kk .or. kfpla(i, j, m) <= kk) then

               ! Index of last layer containing mass at either of the scalar
               ! points adjacent to the velocity point.
               kmax = 1
               do k = 3, kk
                  km = k + mm
                  if (dp(i, j - 1, km) > epsil .or. dp(i, j, km) > epsil) &
                     kmax = k
               enddo

               ! The first interior interface where the y-component of the slope
               ! vector is estimated is at index kintr + 1.
               kintr = max(kfpla(i, j - 1, m), kfpla(i, j, m))

               ! Index of last interface where slope vector times Brunt-Vaisala
               ! frequency is computed.
               knnsl = 2

               ! Compute the y-component of the slope vector at the mixed layer
               ! base.
               pm = .5_r8*(p(i, j - 1, 3) + p(i, j, 3))
               rho_y = rho(pm, temp(i, j    , 2 + mm), saln(i, j    , 2 + mm)) &
                     - rho(pm, temp(i, j - 1, 2 + mm), saln(i, j - 1, 2 + mm))
               phi_y = phi(i, j, 3) - phi(i, j - 1, 3)
               bfsqm = .5_r8*(bfsqf(i, j - 1, 3) + bfsqf(i, j, 3))
               nslpy(i, j, 3) = (g*rho_y/(rho0*bfsqm) + phi_y/g)*scvyi(i, j)
               if (phi(i, j    , 3) > phi(i, j - 1, kk + 1) .and. &
                   phi(i, j - 1, 3) > phi(i, j    , kk + 1)) then
                  nnslpy(i, j, 3) = sqrt(bfsqm)*nslpy(i, j, 3)
                  knnsl = 3
               endif

               ! Compute the y-component of the slope vector at interior
               ! interfaces.
               do k = kintr + 1, kmax
                  km = k + mm
                  pm = .5_r8*(p(i, j - 1, k) + p(i, j, k))
                  rho_y = .5_r8*( rho(pm, temp(i, j    , km - 1), &
                                          saln(i, j    , km - 1)) &
                                - rho(pm, temp(i, j - 1, km - 1), &
                                          saln(i, j - 1, km - 1)) &
                                + rho(pm, temp(i, j    , km    ), &
                                          saln(i, j    , km    )) &
                                - rho(pm, temp(i, j - 1, km    ), &
                                          saln(i, j - 1, km    )))
                  phi_y = phi(i, j, k) - phi(i, j - 1, k)
                  bfsqm = .5_r8*(bfsqf(i, j - 1, k) + bfsqf(i, j, k))
                  nslpy(i, j, k) = (g*rho_y/(rho0*bfsqm) + phi_y/g)*scvyi(i, j)
                  if (phi(i, j    , k) > phi(i, j - 1, kk + 1) .and. &
                      phi(i, j - 1, k) > phi(i, j    , kk + 1)) then
                     nnslpy(i, j, k) = sqrt(bfsqm)*nslpy(i, j, k)
                     knnsl = k
                  endif
               enddo
               do k = knnsl + 1, kmax
                  nnslpy(i, j, k) = nnslpy(i, j, knnsl)
               enddo
               if (kintr < kmax) then
                  do k = 4, kintr
                     nslpy(i, j, k) = nslpy(i, j, kintr + 1)
                     nnslpy(i, j, k) = nnslpy(i, j, kintr + 1)
                  enddo
               else
                  do k = 4, kmax
                     nslpy(i, j, k) = nslpy(i, j, 3)
                     nnslpy(i, j, k) = nnslpy(i, j, 3)
                  enddo
               endif

            endif

         enddo
         enddo
      enddo
   !$omp end parallel do

      if (csdiag) then
         if (mnproc == 1) then
           write (lp,*) 'cmnfld_nslope_isopyc_bulkml:'
         endif
         call chksummsk(nslpx, iu, kk, 'nslpx')
         call chksummsk(nslpy, iv, kk, 'nslpy')
         call chksummsk(nnslpx, iu, kk, 'nnslpx')
         call chksummsk(nnslpy, iv, kk, 'nnslpy')
      endif

   end subroutine cmnfld_nslope_isopyc_bulkml

   subroutine cmnfld_nslope_cntiso_hybrid(m, n, mm, nn, k1m, k1n)
   ! ---------------------------------------------------------------------------
   ! Estimate slope of local neutral surface.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      real(r8) :: rho0, pm, rho_x, phi_x, bfsqm, rho_y, phi_y
      integer :: i, j, k, l, km, kmax, knnsl

      ! ------------------------------------------------------------------------
      ! Compute geopotential at layer interfaces.
      ! ------------------------------------------------------------------------

   !$omp parallel do private(k, km, l, i)
      do j = - 1, jj + 2
         do k = kk, 1, - 1
            km = k + mm
            do l = 1, isp(j)
            do i = max(- 1, ifp(j, l)), min(ii + 2, ilp(j, l))
               if (dp(i, j, km) < epsil) then
                  phi(i, j, k) = phi(i, j, k + 1)
               else
                  phi(i, j, k) = phi(i, j, k + 1) &
                               - p_alpha(p(i, j, k + 1), p(i, j, k), &
                                         temp(i, j, km), saln(i, j, km))
               endif
            enddo
            enddo
         enddo
      enddo
   !$omp end parallel do

      ! ------------------------------------------------------------------------
      ! Compute slope vector of local neutral surfaces and also slope vector
      ! times Brunt-Vaisala frequency (optionally used in the computation of
      ! eddy growth rate). The latter is not computed when the gradient of the
      ! geopotential is expected to be influenced by the gradient of the
      ! bathymetry and in this case values are extrapolated from above.
      ! ------------------------------------------------------------------------

      rho0 = 1._r8/alpha0

   !$omp parallel do private(l, i, k, kmax, km, knnsl, pm, rho_x, phi_x, bfsqm)
      do j = - 1, jj + 2
         do l = 1, isu(j)
         do i = max(0, ifu(j, l)), min(ii + 2, ilu(j, l))

            ! Set the x-component of the slope vector to zero initially.
            do k = 1, kk
               nslpx(i, j, k) = 0._r8
               nnslpx(i, j, k) = 0._r8
            enddo

            ! Index of last layer containing mass at either of the scalar
            ! points adjacent to the velocity point.
            kmax = 1
            do k = 2, kk
               km = k + mm
               if (dp(i - 1, j, km) > epsil .or. dp(i, j, km) > epsil) kmax = k
            enddo

            ! Index of last interface where slope vector times Brunt-Vaisala
            ! frequency is computed.
            knnsl = 2

            ! Compute the x-component of the slope vector at interfaces.
            do k = 2, kmax
               km = k + mm
               pm = .5_r8*(p(i - 1, j, k) + p(i, j, k))
               rho_x = .5_r8*( rho(pm, temp(i    , j, km - 1), &
                                       saln(i    , j, km - 1)) &
                             - rho(pm, temp(i - 1, j, km - 1), &
                                       saln(i - 1, j, km - 1)) &
                             + rho(pm, temp(i    , j, km    ), &
                                       saln(i    , j, km    )) &
                             - rho(pm, temp(i - 1, j, km    ), &
                                       saln(i - 1, j, km    )))
               phi_x = phi(i, j, k) - phi(i - 1, j, k)
               bfsqm = .5_r8*(bfsqf(i - 1, j, k) + bfsqf(i, j, k))
               nslpx(i, j, k) = (g*rho_x/(rho0*bfsqm) + phi_x/g)*scuxi(i, j)
               if (phi(i    , j, k) > phi(i - 1, j, kk + 1) .and. &
                   phi(i - 1, j, k) > phi(i    , j, kk + 1)) then
                  nnslpx(i, j, k) = sqrt(bfsqm)*nslpx(i, j, k)
                  knnsl = k
               endif
            enddo
            do k = knnsl + 1, kmax
               nnslpx(i, j, k) = nnslpx(i, j, knnsl)
            enddo

         enddo
         enddo
      enddo
   !$omp end parallel do

   !$omp parallel do private(l, i, k, kmax, km, knnsl, pm, rho_y, phi_y, bfsqm)
      do j = 0, jj + 2
         do l = 1, isv(j)
         do i = max(- 1, ifv(j, l)), min(ii + 2, ilv(j, l))

            ! Set the y-component of the slope vector to zero initially.
            do k = 1, kk
               nslpy(i, j, k) = 0._r8
               nnslpy(i, j, k) = 0._r8
            enddo

            ! Index of last layer containing mass at either of the scalar
            ! points adjacent to the velocity point.
            kmax = 1
            do k = 2, kk
               km = k + mm
               if (dp(i, j - 1, km) > epsil .or. dp(i, j, km) > epsil) kmax = k
            enddo

            ! Index of last interface where slope vector times Brunt-Vaisala
            ! frequency is computed.
            knnsl = 2

            ! Compute the y-component of the slope vector at interfaces.
            do k = 2, kmax
               km = k + mm
               pm = .5_r8*(p(i, j - 1, k) + p(i, j, k))
               rho_y = .5_r8*( rho(pm, temp(i, j    , km - 1), &
                                       saln(i, j    , km - 1)) &
                             - rho(pm, temp(i, j - 1, km - 1), &
                                       saln(i, j - 1, km - 1)) &
                             + rho(pm, temp(i, j    , km    ), &
                                       saln(i, j    , km    )) &
                             - rho(pm, temp(i, j - 1, km    ), &
                                       saln(i, j - 1, km    )))
               phi_y = phi(i, j, k) - phi(i, j - 1, k)
               bfsqm = .5_r8*(bfsqf(i, j - 1, k) + bfsqf(i, j, k))
               nslpy(i, j, k) = (g*rho_y/(rho0*bfsqm) + phi_y/g)*scvyi(i, j)
               if (phi(i, j    , k) > phi(i, j - 1, kk + 1) .and. &
                   phi(i, j - 1, k) > phi(i, j    , kk + 1)) then
                  nnslpy(i, j, k) = sqrt(bfsqm)*nslpy(i, j, k)
                  knnsl = k
               endif
            enddo
            do k = knnsl + 1, kmax
               nnslpy(i, j, k) = nnslpy(i, j, knnsl)
            enddo

         enddo
         enddo
      enddo
   !$omp end parallel do

      if (csdiag) then
         if (mnproc == 1) then
           write (lp,*) 'cmnfld_nslope_cntiso_hybrid:'
         endif
         call chksummsk(nslpx, iu, kk, 'nslpx')
         call chksummsk(nslpy, iv, kk, 'nslpy')
         call chksummsk(nnslpx, iu, kk, 'nnslpx')
         call chksummsk(nnslpy, iv, kk, 'nnslpy')
      endif

   end subroutine cmnfld_nslope_cntiso_hybrid

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
         do k = 1, kk
            do i = 1 - nbdy, ii + nbdy
               bfsqi (i, j, k) = spval
               bfsql (i, j, k) = spval
               bfsqf (i, j, k) = spval
               nslpx (i, j, k) = spval
               nslpy (i, j, k) = spval
               nnslpx(i, j, k) = spval
               nnslpy(i, j, k) = spval
            enddo
         enddo
      enddo
   !$omp end parallel do

   end subroutine inivar_cmnfld

   subroutine cmnfld(m, n, mm, nn, k1m, k1n)
   ! ---------------------------------------------------------------------------
   ! Compute fields that are used by several subsequent routines
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      integer :: i, j, l

      ! ------------------------------------------------------------------------
      ! Update halos of various fields.
      ! ------------------------------------------------------------------------

      call xctilr(temp, 1, 2*kk, 3, 3, halo_ps)
      call xctilr(saln, 1, 2*kk, 3, 3, halo_ps)
!     call xctilr(temp(1 - nbdy, 1 - nbdy, k1n), 1, kk, 3, 3, halo_ps)
!     call xctilr(saln(1 - nbdy, 1 - nbdy, k1n), 1, kk, 3, 3, halo_ps)
   !$omp parallel do private(l, i)
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
            util1(i, j) = kfpla(i, j, m)
         enddo
         enddo
      enddo
   !$omp end parallel do
      call xctilr(util1, 1, 1, 2, 2, halo_ps)
   !$omp parallel do private(l, i)
      do j = - 1, jj + 2
         do l = 1, isp(j)
         do i = max(- 1, ifp(j, l)), min(ii + 2, ilp(j, l))
            kfpla(i, j, m) = nint(util1(i, j))
         enddo
         enddo
      enddo
   !$omp end parallel do

      ! ------------------------------------------------------------------------
      ! Compute fields depending on selection of physics.
      ! ------------------------------------------------------------------------

!     if (edritp == 'large scale' .or. eitmth == 'gm' .or.
!    .    sum(ACC_BFSQ(1:nphy)).ne.0) then
      if (edritp == 'large scale' .or. eitmth == 'gm') then

         ! ---------------------------------------------------------------------
         ! Compute filtered buoyancy frequency squared.
         ! ---------------------------------------------------------------------

         if (vcoord_type_tag == isopyc_bulkml) then
            call cmnfld_bfsqf_isopyc_bulkml(m, n, mm, nn, k1m, k1n)
         else
            call cmnfld_bfsqf_cntiso_hybrid(m, n, mm, nn, k1m, k1n)
         endif

      endif

      if (edritp == 'large scale' .or. eitmth == 'gm') then

         ! ---------------------------------------------------------------------
         ! Estimate slope of local neutral surface.
         ! ---------------------------------------------------------------------

         if (vcoord_type_tag == isopyc_bulkml) then
            call cmnfld_nslope_isopyc_bulkml(m, n, mm, nn, k1m, k1n)
         else
            call cmnfld_nslope_cntiso_hybrid(m, n, mm, nn, k1m, k1n)
         endif

      endif

      end subroutine cmnfld

end module mod_cmnfld
