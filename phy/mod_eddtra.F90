! ------------------------------------------------------------------------------
! Copyright (C) 2015-2021 Mats Bentsen
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

module mod_eddtra
! ------------------------------------------------------------------------------
! This module contains procedures related to the computation of eddy-induced
! transport.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_constants, only: g, alpha0, epsil, onemm
   use mod_time, only: delt1
   use mod_xc
   use mod_grid, only: scuy, scvx, scp2, scu2, scv2, scuxi, scvyi
   use mod_eos, only: rho
   use mod_state, only: dp, dpu, dpv, temp, saln, p, pbu, pbv, kfpla
   use mod_diffusion, only: eitmth, difint, umfltd, vmfltd, &
                            utfltd, vtfltd, usfltd, vsfltd
   use mod_cmnfld, only: nslpx, nslpy, nnslpx, nnslpy
   use mod_checksum, only: csdiag, chksummsk

   implicit none

   private

   public :: eddtra

contains

   ! ---------------------------------------------------------------------------
   ! Private procedures.
   ! ---------------------------------------------------------------------------

   subroutine eddtra_intdif(m, n, mm, nn, k1m, k1n)
   ! ---------------------------------------------------------------------------
   ! Estimate eddy-induced transport by interface diffusion.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      real(r8) :: flxhi, flxlo, q
      integer :: i, j, k, l, km, kn

      call xctilr(difint, 1,kk, 2,2, halo_ps)

   !$omp parallel do private(l, i)
      do j = - 1, jj + 2
         do l = 1, isu(j)
         do i = max(0, ifu(j, l)), min(ii + 2, ilu(j, l))
            umfltd(i, j, 1 + mm) = 0._r8
            umfltd(i, j, 2 + mm) = 0._r8
            umfltd(i, j, 3 + mm) = 0._r8
         enddo
         enddo
      enddo
   !$omp end parallel do
   !$omp parallel do private(l, i)
      do j = 0, jj + 2
         do l = 1, isv(j)
         do i = max(- 1, ifv(j, l)), min(ii + 2, ilv(j, l))
            vmfltd(i, j, 1 + mm) = 0._r8
            vmfltd(i, j, 2 + mm) = 0._r8
            vmfltd(i, j, 3 + mm) = 0._r8
         enddo
         enddo
      enddo
   !$omp end parallel do

      do k = 4, kk
         km = k + mm
         kn = k + nn

      !$omp parallel do private(l, i, flxhi, flxlo, q)
         do j = - 1, jj + 2
            do l = 1, isu(j)
            do i = max(0, ifu(j, l)), min(ii + 2, ilu(j, l))
               flxhi =   .125_r8*min(dp(i - 1, j, kn - 1)*scp2(i - 1, j), &
                                     dp(i    , j, kn    )*scp2(i    , j))
               flxlo = - .125_r8*min(dp(i    , j, kn - 1)*scp2(i    , j), &
                                     dp(i - 1, j, kn    )*scp2(i - 1, j))
               q = .25_r8*( difint(i - 1, j, k - 1) + difint(i, j, k - 1) &
                          + difint(i - 1, j, k    ) + difint(i, j, k    ))
               q = min(flxhi, max(flxlo, &
                       delt1*q*(p(i - 1, j, k) - p(i, j, k)) &
                       *scuy(i, j)*scuxi(i, j)))
               umfltd(i, j, km - 1) = umfltd(i, j, km - 1) + q
               umfltd(i, j, km    ) = - q
            enddo
            enddo
         enddo
      !$omp end parallel do

      !$omp parallel do private(l, i, flxhi, flxlo, q)
         do j = 0, jj + 2
            do l = 1, isv(j)
            do i = max(- 1, ifv(j, l)), min(ii + 2, ilv(j, l))
               flxhi =   .125_r8*min(dp(i, j - 1, kn - 1)*scp2(i, j - 1), &
                                     dp(i, j    , kn    )*scp2(i, j    ))
               flxlo = - .125_r8*min(dp(i, j    , kn - 1)*scp2(i, j    ), &
                                     dp(i, j - 1, kn    )*scp2(i, j - 1))
               q = .25_r8*( difint(i, j - 1, k - 1) + difint(i, j, k - 1) &
                          + difint(i, j - 1, k    ) + difint(i, j, k    ))
               q = min(flxhi, max(flxlo, &
                       delt1*q*(p(i, j - 1, k) - p(i, j, k)) &
                       *scvx(i, j)*scvyi(i, j)))
               vmfltd(i, j, km - 1) = vmfltd(i, j, km - 1) + q
               vmfltd(i, j, km    ) = - q
            enddo
            enddo
         enddo
      !$omp end parallel do

      enddo

   end subroutine eddtra_intdif

   subroutine eddtra_gm(m, n, mm, nn, k1m, k1n)
   ! ---------------------------------------------------------------------------
   ! Estimate eddy-induced transport following the Gent-McWilliams
   ! parameterization.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      ! Parameters:
      real(r8), parameter :: &
         ffac  = .0625_r8, & ! Fraction of the mass of a grid cell a mass flux
                             ! is allowed to deplete [].
         fface = .99_r8*ffac ! (1-epsilon)*ffac [].

      real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ptu, ptv
      real(r8), dimension(kdm+1) :: upsilon, mfl
      real(r8), dimension(kdm) :: dlm, dlp
      real(r8) :: rho0, q, et2mf, kappa, fhi, flo
      integer :: i, j, k, l, km, kn, kintr, kmax, kmin, niter, kdir
      logical :: changed

      rho0 = 1._r8/alpha0

      call xctilr(difint, 1, kk, 2, 2, halo_ps)
      call xctilr(pbu, 1, 2, 2, 2, halo_us)
      call xctilr(pbv, 1, 2, 2, 2, halo_vs)


      ! Compute top pressure at velocity points.
   !$omp parallel do private(l, i)
      do j= - 1, jj + 2
         do l = 1, isu(j)
         do i = max(0, ifu(j, l)), min(ii + 2, ilu(j, l))
            ptu(i, j) = max(p(i - 1, j, 1), p(i, j, 1))
         enddo
         enddo
      enddo
   !$omp end parallel do
   !$omp parallel do private(l, i)
      do j = 0, jj + 2
         do l = 1, isv(j)
         do i = max( - 1, ifv(j, l)), min(ii + 2, ilv(j, l))
            ptv(i, j) = max(p(i, j - 1, 1), p(i, j, 1))
         enddo
         enddo
      enddo
   !$omp end parallel do

     ! -------------------------------------------------------------------------
     ! Compute u-component of eddy-induced mass fluxes.
     ! -------------------------------------------------------------------------

   !$omp parallel do private(l, i, k, km, et2mf, kmax, kn, kintr, kappa, &
   !$omp                     upsilon, kmin, mfl, dlm, dlp, fhi, flo, changed, &
   !$omp                     niter, kdir, q)
      do j = - 1, jj + 2
         do l = 1, isu(j)
         do i = max(0, ifu(j, l)), min(ii + 2, ilu(j, l))

            ! Set eddy-induced mass fluxes to zero initially.
            do k = 1, kk
               km = k + mm
               umfltd(i, j, km) = 0._r8
            enddo

            ! Eddy transport to mass flux conversion factor.
            et2mf = - g*rho0*delt1*scuy(i, j)

            ! Index of last layer containing mass at either of the scalar points
            ! adjacent to the velocity point.
            kmax = 1
            do k = 3, kk
               kn = k + nn
               if (dp(i - 1, j, kn) > epsil .or. dp(i, j, kn) > epsil) kmax = k
            enddo

            ! ------------------------------------------------------------------
            ! Proceed with mass flux computation if at least one of the adjacent
            ! scalar points to the velocity point has a mass containing interior
            ! layer. Mass fluxes will be assigned at layer interface
            ! corresponding to the eddy induced transport. The final layer mass
            ! flux will be the lower minus the upper interface flux. The mass
            ! fluxes are limited to keep interfaces within the water column.
            ! There are 3 cases to consider:
            !   Case 1: The mixed layer extends to the bottom at both adjacent
            !           scalar points to the velocity point
            !   Case 2: The mixed layer extends to the bottom at scalar point
            !           (i, j).
            !   Case 3: The mixed layer extends to the bottom at scalar point
            !           (i - 1, j).
            !   Case 4: The mixed layer does not reach the bottom at neither of
            !           the scalar points adjacent to the velocity point.
            ! ------------------------------------------------------------------

            if     (kfpla(i - 1, j, n) >  kk .and. kfpla(i, j, n) > kk) then
               ! ---------------------------------------------------------------
               ! Case 1:
               ! ---------------------------------------------------------------

               ! Keep the initial zero mass fluxes for this column.
               cycle

            elseif (kfpla(i - 1, j, n) <= kk .and. kfpla(i, j, n) > kk) then
               ! ---------------------------------------------------------------
               ! Case 2:
               ! ---------------------------------------------------------------

               ! Find the index of the first layer at (i - 1, j) that is
               ! hydrostatically stable at the mixed layer base at (i, j).
               km = 2 + nn
               kintr = kfpla(i - 1, j, n)
               kn = kintr + nn
               do while (rho(p(i    , j, 3), &
                             temp(i - 1, j, kn), saln(i - 1, j, kn)) < &
                         rho(p(i    , j, 3), &
                             temp(i    , j, km), saln(i    , j, km)) .or. &
                         dp(i - 1, j, kn) < epsil)
                  kintr = kintr + 1
                  if (kintr == kmax + 1) exit
                  kn = kintr + nn
               enddo

               ! If a physical layer cannot be found, keep the initial zero mass
               ! fluxes for this column.
               if (kintr == kmax + 1) cycle

               ! Compute the eddy induced transport (upsilon) at the mixed layer
               ! base.
               kappa = .5_r8*(difint(i - 1, j, 2) + difint(i, j, 2))
               upsilon(3) = - kappa*nslpx(i, j, 3)

               ! If the eddy-induced transport at the base of the mixed layer
               ! would cause a negative mass flux below the mixed layer, keep
               ! the initial zero mass fluxes for this column.
               if (upsilon(3) <= 0._r8) cycle

               ! Assign interface mass fluxes.
               kmin = kintr - 1
               mfl(kmin) = 0._r8
               mfl(kintr) = et2mf*upsilon(3)
               do k = kintr + 1, kmax + 1
                  mfl(k) = 0._r8
               enddo

            elseif (kfpla(i - 1, j, n) >  kk .and. kfpla(i, j, n) <= kk) then
               ! ---------------------------------------------------------------
               ! Case 3:
               ! ---------------------------------------------------------------

               ! Find the index of the first layer at (i, j) that is
               ! hydrostatically stable at the mixed layer base at (i - 1, j).
               km = 2 + nn
               kintr = kfpla(i    , j, n)
               kn = kintr + nn
               do while (rho(p(i - 1, j, 3), &
                             temp(i    , j, kn), saln(i    , j, kn)) < &
                         rho(p(i - 1, j, 3), &
                             temp(i - 1, j, km), saln(i - 1, j, km)) .or. &
                         dp(i    , j, kn) < epsil)
                  kintr = kintr + 1
                  if (kintr == kmax + 1) exit
                  kn = kintr + nn
               enddo

               ! If a physical layer cannot be found, keep the initial zero mass
               ! fluxes for this column.
               if (kintr == kmax + 1) cycle

               ! Compute the eddy induced transport (upsilon) at the mixed layer
               ! base.
               kappa = .5_r8*(difint(i - 1, j, 2) + difint(i, j, 2))
               upsilon(3) = - kappa*nslpx(i, j, 3)

               ! If the eddy-induced transport at the base of the mixed layer
               ! would cause a positive mass flux below the mixed layer, keep
               ! the initial zero mass fluxes for this column.
               if (upsilon(3) >= 0._r8) cycle

               ! Assign interface mass fluxes.
               kmin = kintr - 1
               mfl(kmin) = 0._r8
               mfl(kintr) = et2mf*upsilon(3)
               do k = kintr + 1, kmax + 1
                  mfl(k) = 0._r8
               enddo

            else
               ! ---------------------------------------------------------------
               ! Case 4:
               ! ---------------------------------------------------------------

               ! The first interior interface where the eddy induced transport
               ! is estimated is at index kintr + 1.
               kintr = max(kfpla(i - 1, j, n), kfpla(i, j, n))

               ! Compute the eddy induced transport (upsilon) at the mixed layer
               ! base.
               kappa = .5_r8*(difint(i - 1, j, 2) + difint(i, j, 2))
               upsilon(3) = - kappa*nslpx(i, j, 3)

               ! Compute the eddy induced transport at interior interfaces.
               do k = kintr + 1, kmax
                  kn = k + nn
                  kappa = .25_r8*( difint(i - 1, j, k - 1) &
                                 + difint(i    , j, k - 1) &
                                 + difint(i - 1, j, k    ) &
                                 + difint(i    , j, k    ))
                  upsilon(k) = - kappa*nslpx(i, j, k)
               enddo
               upsilon(kmax + 1) = 0._r8

               ! If the layer kintr - 1 is a physical layer at either of the
               ! adjacent scalar points to the velocity point, then apply an
               ! upper interface mass flux corresponding to the eddy induced
               ! transport at the mixed layer base and a lower interface mass
               ! flux corresponding to the eddy induced transport at the
               ! kintr + 1 interface if this would lead to a hydrostatically
               ! stable layer arrangement.
               km = 2 + nn
               kn = kintr - 1 + nn
               if ((kfpla(i - 1, j, n) < kintr .and. &
                    upsilon(3) - upsilon(kintr + 1) > 0._r8 .and. &
                    rho(p(i    , j, 3), &
                        temp(i - 1, j, kn), saln(i - 1, j, kn)) > &
                    rho(p(i    , j, 3), &
                        temp(i    , j, km), saln(i    , j, km))) .or. &
                   (kfpla(i    , j, n) < kintr .and. &
                    upsilon(3) - upsilon(kintr + 1) < 0._r8 .and. &
                    rho(p(i - 1, j, 3), &
                        temp(i    , j, kn), saln(i    , j, kn)) > &
                    rho(p(i - 1, j, 3), &
                        temp(i - 1, j, km), saln(i - 1, j, km)))) then
                  kintr = kintr - 1
                  upsilon(kintr + 1) = upsilon(kintr + 2)
               endif

               ! Assign interface mass fluxes.
               kmin = kintr - 1
               mfl(kmin) = 0._r8
               mfl(kintr) = et2mf*upsilon(3)
               do k = kintr + 1, kmax
                  mfl(k) = et2mf*upsilon(k)
               enddo
               mfl(kmax + 1) = 0._r8

            endif

            ! ------------------------------------------------------------------
            ! Ensure that mass fluxes do not create negative layer thicknesses.
            ! ------------------------------------------------------------------

            ! Compute the layer thicknesses available to be depleted by mass
            ! fluxes at the scalar points adjacent to the velocity point. These
            ! bounded layer thicknesses are consistent with the transport
            ! algorithm.
            dlm(kmin) = max(0._r8, min(p(i - 1, j, 3), pbu(i, j, n)) &
                                 - max(p(i - 1, j, 1), ptu(i, j)))
            dlp(kmin) = max(0._r8, min(p(i    , j, 3), pbu(i, j, n)) &
                                 - max(p(i    , j, 1), ptu(i, j)))
            do k = kintr, kmax
               dlm(k) = max(0._r8, min(p(i - 1, j, k + 1), pbu(i, j, n)) &
                                 - max(p(i - 1, j, k    ), ptu(i, j)))
               dlp(k) = max(0._r8, min(p(i    , j, k + 1), pbu(i, j, n)) &
                                 - max(p(i    , j, k    ), ptu(i, j)))
            enddo

            ! If excessive depletion of layers occur beneath the mixed layer
            ! base, try to adjust interface fluxes other that the mixed layer
            ! base interface flux.
            fhi =   fface*max(0._r8, min((p(i - 1, j, 3) - ptu(i, j)) &
                                         *scp2(i - 1, j), &
                                         (pbu(i, j, n) - p(i    , j, kintr)) &
                                         *scp2(i    , j)))
            flo = - fface*max(0._r8, min((p(i    , j, 3) - ptu(i, j)) &
                                         *scp2(i    , j), &
                                         (pbu(i, j, n) - p(i - 1, j, kintr)) &
                                         *scp2(i - 1, j)))
            mfl(kmin + 1) = min(fhi, max(flo, mfl(kmin + 1)))
            do k = kmin + 1, kmax - 1
               if     (mfl(k + 1) - mfl(k) > &
                        ffac*max(epsil, dlm(k))*scp2(i - 1, j)) then
                  mfl(k + 1) = mfl(k) + fface*dlm(k)*scp2(i - 1, j)
               elseif (mfl(k + 1) - mfl(k) < &
                      - ffac*max(epsil, dlp(k))*scp2(i    , j)) then
                  mfl(k + 1) = mfl(k) - fface*dlp(k)*scp2(i    , j)
               else
                  exit
               endif
            enddo

            ! Apply an iterative procedure for flux limiting by alternate upward
            ! and downward propagation through the layers.

            changed = .true.
            niter = 0
            kdir = 1

            do while (changed)

               niter = niter + 1
               if (niter == 1000) then
                  k = kmin
                  write(lp,*)
                  write(lp,'(i3,3e16.8)') &
                     1, mfl(k + 1), mfl(k), &
                     (mfl(k + 1) - mfl(k)) &
                     /(max(onemm, dpu(i, j, 1 + nn) + dpu(i, j, 2 + nn)) &
                       *delt1*scuy(i, j))
                  do k = kintr, kmax
                     kn = k + nn
                     write(lp,'(i3,3e16.8)') &
                        k, mfl(k + 1), mfl(k), &
                        (mfl(k + 1) - mfl(k)) &
                        /(max(onemm, dpu(i, j, kn))*delt1*scuy(i, j))
                  enddo
                  write(lp,*) 'no convergence u', i + i0, j + j0
                  call xchalt('(eddtra_gm)')
                         stop '(eddtra_gm)'
               endif

               changed = .false.
               kdir = - kdir

               do k = ((1 - kdir)*kmax + (1 + kdir)*kmin)/2, &
                      ((1 - kdir)*kmin + (1 + kdir)*kmax)/2, kdir

                  ! Proceed with flux limiting of this layer if the mass flux
                  ! difference between lower and upper interface is beyond the
                  ! floating point accuracy limitation.
                  if (abs(mfl(k + 1) - mfl(k)) > &
                      1.e-14_r8*max(epsil*scu2(i, j), &
                                    abs(mfl(k + 1) + mfl(k)))) then

                     if     (mfl(k + 1) - mfl(k) > &
                              ffac*max(epsil, dlm(k))*scp2(i - 1, j)) then
                        ! In this case, the mass fluxes are removing too much
                        ! mass from the grid cell at (i - 1, j, k). Limit the
                        ! dominating interface flux.
                        q = fface*dlm(k)*scp2(i - 1, j)
                        if (mfl(k + 1) > - mfl(k)) then
                           if (mfl(k    ) > - .5_r8*q) then
                              mfl(k + 1) =   mfl(k    ) + q
                           else
                              mfl(k + 1) =   .5_r8*q
                              mfl(k    ) = - mfl(k + 1)
                           endif
                        else
                           if (mfl(k + 1) <   .5_r8*q) then
                              mfl(k    ) =   mfl(k + 1) - q
                           else
                              mfl(k    ) = - .5_r8*q
                              mfl(k + 1) = - mfl(k    )
                           endif
                        endif
                        changed = .true.
                     elseif (mfl(k + 1) - mfl(k) < &
                            - ffac*max(epsil, dlp(k))*scp2(i    , j)) then
                        ! In this case, the mass fluxes are removing too much
                        ! mass from the grid cell at (i, j, k). Limit the
                        ! dominating interface flux.
                        q = fface*dlp(k)*scp2(i    , j)
                        if (mfl(k + 1) < - mfl(k)) then
                           if (mfl(k    ) <   .5_r8*q) then
                              mfl(k + 1) =   mfl(k    ) - q
                           else
                              mfl(k + 1) = - .5_r8*q
                              mfl(k    ) = - mfl(k + 1)
                           endif
                        else
                           if (mfl(k + 1) > - .5_r8*q) then
                              mfl(k    ) =   mfl(k + 1) + q
                           else
                              mfl(k    ) =   .5_r8*q
                              mfl(k + 1) = - mfl(k    )
                           endif
                        endif
                        changed = .true.
                     endif
                  endif

               enddo

            enddo

            ! ------------------------------------------------------------------
            ! Compute the final mass fluxes.
            ! ------------------------------------------------------------------

            k = kmin
            if (abs(mfl(k + 1) - mfl(k)) > &
                1.e-14_r8*max(epsil*scu2(i, j), &
                              abs(mfl(k + 1) + mfl(k)))) then
               umfltd(i, j, 2 + mm) = mfl(k + 1) - mfl(k)
               umfltd(i, j, 1 + mm) = umfltd(i, j, 2 + mm) &
                                      *dpu(i, j, 1 + nn)/( dpu(i, j, 1 + nn) &
                                                         + dpu(i, j, 2 + nn))
               umfltd(i, j, 2 + mm) = umfltd(i, j, 2 + mm) &
                                    - umfltd(i, j, 1 + mm)
            else
               umfltd(i, j, 1 + mm) = 0._r8
               umfltd(i, j, 2 + mm) = 0._r8
            endif
            do k = kintr, kmax
               km = k + mm
               if (abs(mfl(k + 1) - mfl(k)) > &
                   1.e-14_r8*max(epsil*scu2(i, j), &
                                 abs(mfl(k + 1) + mfl(k)))) then
                  umfltd(i, j, km) = mfl(k + 1) - mfl(k)
               else
                  umfltd(i, j, km) = 0._r8
               endif
               if (umfltd(i, j, km) > &
                    ffac*max(epsil, dlm(k))*scp2(i - 1, j)) then
                  write(lp,*) 'eddtra_gm u >', &
                              i + i0, j + j0, k, umfltd(i, j, km), &
                              ffac*max(epsil, dlm(k))*scp2(i - 1, j)
                  call xchalt('(eddtra_gm)')
                         stop '(eddtra_gm)'
               endif
               if (umfltd(i, j, km) < &
                  - ffac*max(epsil, dlp(k))*scp2(i    , j)) then
                  write(lp,*) 'eddtra_gm u <', &
                              i + i0, j + j0, k, umfltd(i, j, km), &
                            - ffac*max(epsil, dlp(k))*scp2(i    , j)
                  call xchalt('(eddtra_gm)')
                         stop '(eddtra_gm)'
               endif
            enddo

         enddo
         enddo
      enddo
   !$omp end parallel do

     ! -------------------------------------------------------------------------
     ! Compute v-component of eddy-induced mass fluxes.
     ! -------------------------------------------------------------------------

   !$omp parallel do private(l, i, k, km, et2mf, kmax, kn, kintr, kappa, &
   !$omp                     upsilon, kmin, mfl, dlm, dlp, fhi, flo, changed, &
   !$omp                     niter, kdir, q)
      do j = 0, jj + 2
         do l = 1, isv(j)
         do i = max( - 1, ifv(j, l)), min(ii + 2, ilv(j, l))

            ! Set eddy-induced mass fluxes to zero initially.
            do k = 1, kk
               km = k + mm
               vmfltd(i, j, km) = 0._r8
            enddo

            ! Eddy transport to mass flux conversion factor.
            et2mf = - g*rho0*delt1*scvx(i, j)

            ! Index of last layer containing mass at either of the scalar points
            ! adjacent to the velocity point.
            kmax = 1
            do k = 3, kk
               kn = k + nn
               if (dp(i, j - 1, kn) > epsil .or. dp(i, j, kn) > epsil) kmax = k
            enddo

            ! ------------------------------------------------------------------
            ! Proceed with mass flux computation if at least one of the adjacent
            ! scalar points to the velocity point has a mass containing interior
            ! layer. Mass fluxes will be assigned at layer interface
            ! corresponding to the eddy induced transport. The final layer mass
            ! flux will be the lower minus the upper interface flux. The mass
            ! fluxes are limited to keep interfaces within the water column.
            ! There are 3 cases to consider:
            !   Case 1: The mixed layer extends to the bottom at both adjacent
            !           scalar points to the velocity point
            !   Case 2: The mixed layer extends to the bottom at scalar point
            !           (i, j).
            !   Case 3: The mixed layer extends to the bottom at scalar point
            !           (i, j - 1).
            !   Case 4: The mixed layer does not reach the bottom at neither of
            !           the scalar points adjacent to the velocity point.
            ! ------------------------------------------------------------------

            if     (kfpla(i, j - 1, n) > kk .and. kfpla(i, j, n) > kk) then
               ! ---------------------------------------------------------------
               ! Case 1:
               ! ---------------------------------------------------------------

               ! Keep the initial zero mass fluxes for this column.
               cycle

            elseif (kfpla(i, j - 1, n) <= kk .and. kfpla(i, j, n) > kk) then
               ! ---------------------------------------------------------------
               ! Case 2:
               ! ---------------------------------------------------------------

               ! Find the index of the first layer at (i, j - 1) that is
               ! hydrostatically stable at the mixed layer base at (i, j).
               km = 2 + nn
               kintr = kfpla(i, j - 1, n)
               kn = kintr + nn
               do while (rho(p(i, j    , 3), &
                             temp(i, j - 1, kn), saln(i, j - 1, kn)) < &
                         rho(p(i, j    , 3), &
                             temp(i, j    , km), saln(i, j    , km)) .or. &
                         dp(i, j - 1, kn) < epsil)
                  kintr = kintr + 1
                  if (kintr == kmax + 1) exit
                  kn = kintr + nn
               enddo

               ! If a physical layer cannot be found, keep the initial zero mass
               ! fluxes for this column.
               if (kintr == kmax + 1) cycle

               ! Compute the eddy induced transport (upsilon) at the mixed layer
               ! base.
               kappa = .5_r8*(difint(i, j - 1, 2) + difint(i, j, 2))
               upsilon(3) = - kappa*nslpy(i, j, 3)

               ! If the eddy-induced transport at the base of the mixed layer
               ! would cause a negative mass flux below the mixed layer, keep
               ! the initial zero mass fluxes for this column.
               if (upsilon(3) <= 0._r8) cycle

               ! Assign interface mass fluxes.
               kmin = kintr - 1
               mfl(kmin) = 0._r8
               mfl(kintr) = et2mf*upsilon(3)
               do k = kintr + 1, kmax + 1
                  mfl(k) = 0._r8
               enddo

            elseif (kfpla(i, j - 1, n) > kk .and. kfpla(i, j, n) <= kk) then
               ! ---------------------------------------------------------------
               ! Case 3:
               ! ---------------------------------------------------------------

               ! Find the index of the first layer at (i, j) that is
               ! hydrostatically stable at the mixed layer base at (i, j - 1).
               km = 2 + nn
               kintr = kfpla(i, j    , n)
               kn = kintr + nn
               do while (rho(p(i, j - 1, 3), &
                             temp(i, j    , kn), saln(i, j    , kn)) < &
                         rho(p(i, j - 1, 3), &
                             temp(i, j - 1, km), saln(i, j - 1, km)) .or. &
                         dp(i, j    , kn) < epsil)
                  kintr = kintr + 1
                  if (kintr == kmax + 1) exit
                  kn = kintr + nn
               enddo

               ! If a physical layer cannot be found, keep the initial zero mass
               ! fluxes for this column.
               if (kintr == kmax + 1) cycle

               ! Compute the eddy induced transport (upsilon) at the mixed layer
               ! base.
               kappa = .5_r8*(difint(i, j - 1, 2) + difint(i, j, 2))
               upsilon(3) = - kappa*nslpy(i, j, 3)

               ! If the eddy-induced transport at the base of the mixed layer
               ! would cause a positive mass flux below the mixed layer, keep
               ! the initial zero mass fluxes for this column.
               if (upsilon(3) >= 0._r8) cycle

               ! Assign interface mass fluxes.
               kmin = kintr - 1
               mfl(kmin) = 0._r8
               mfl(kintr) = et2mf*upsilon(3)
               do k = kintr + 1, kmax + 1
                  mfl(k) = 0._r8
               enddo

            else
               ! ---------------------------------------------------------------
               ! Case 4:
               ! ---------------------------------------------------------------

               ! The first interior interface where the eddy induced transport
               ! is estimated is at index kintr + 1.
               kintr = max(kfpla(i, j - 1, n), kfpla(i, j, n))

               ! Compute the eddy induced transport (upsilon) at the mixed layer
               ! base.
               kappa = .5_r8*(difint(i, j - 1, 2) + difint(i, j, 2))
               upsilon(3) = - kappa*nslpy(i, j, 3)

               ! Compute the eddy induced transport at interior interfaces.
               do k = kintr + 1, kmax
                  kn = k + nn
                  kappa = .25_r8*( difint(i, j - 1, k - 1) &
                                 + difint(i, j    , k - 1) &
                                 + difint(i, j - 1, k    ) &
                                 + difint(i, j    , k    ))
                  upsilon(k) = - kappa*nslpy(i, j, k)
               enddo
               upsilon(kmax + 1) = 0._r8

               ! If the layer kintr - 1 is a physical layer at either of the
               ! adjacent scalar points to the velocity point, then apply an
               ! upper interface mass flux corresponding to the eddy induced
               ! transport at the mixed layer base and a lower interface mass
               ! flux corresponding to the eddy induced transport at the
               ! kintr + 1 interface if this would lead to a hydrostatically
               ! stable layer arrangement.
               km = 2 + nn
               kn = kintr - 1 + nn
               if ((kfpla(i, j - 1, n) < kintr .and. &
                    upsilon(3) - upsilon(kintr + 1) > 0._r8 .and. &
                    rho(p(i, j    , 3), &
                        temp(i, j - 1, kn), saln(i, j - 1, kn)) > &
                    rho(p(i, j    , 3), &
                        temp(i, j    , km), saln(i, j    , km))) .or. &
                   (kfpla(i, j    , n) < kintr .and. &
                    upsilon(3) - upsilon(kintr + 1) < 0._r8 .and. &
                    rho(p(i, j - 1, 3), &
                        temp(i, j    , kn), saln(i, j    , kn)) > &
                    rho(p(i, j - 1, 3), &
                        temp(i, j - 1, km), saln(i, j - 1, km)))) then
                  kintr = kintr - 1
                  upsilon(kintr + 1) = upsilon(kintr + 2)
               endif

               ! Assign interface mass fluxes.
               kmin = kintr - 1
               mfl(kmin) = 0._r8
               mfl(kintr) = et2mf*upsilon(3)
               do k = kintr + 1, kmax
                  mfl(k) = et2mf*upsilon(k)
               enddo
               mfl(kmax + 1) = 0._r8

            endif

            ! ------------------------------------------------------------------
            ! Ensure that mass fluxes do not create negative layer thicknesses.
            ! ------------------------------------------------------------------

            ! Compute the layer thicknesses available to be depleted by mass
            ! fluxes at the scalar points adjacent to the velocity point. These
            ! bounded layer thicknesses are consistent with the transport
            ! algorithm.
            dlm(kmin) = max(0._r8, min(p(i, j - 1, 3), pbv(i, j, n)) &
                                 - max(p(i, j - 1, 1), ptv(i, j)))
            dlp(kmin) = max(0._r8, min(p(i, j    , 3), pbv(i, j, n)) &
                                 - max(p(i, j    , 1), ptv(i, j)))
            do k = kintr, kmax
               dlm(k) = max(0._r8, min(p(i, j - 1, k + 1), pbv(i, j, n)) &
                                 - max(p(i, j - 1, k    ), ptv(i, j)))
               dlp(k) = max(0._r8, min(p(i, j    , k + 1), pbv(i, j, n)) &
                                 - max(p(i, j    , k    ), ptv(i, j)))
            enddo

            ! If excessive depletion of layers occur beneath the mixed layer
            ! base, try to adjust interface fluxes other that the mixed layer
            ! base interface flux.
            fhi =   fface*max(0._r8, min((p(i, j - 1, 3) - ptv(i, j)) &
                                         *scp2(i, j - 1), &
                                         (pbv(i, j, n) - p(i, j    , kintr)) &
                                         *scp2(i, j    )))
            flo = - fface*max(0._r8, min((p(i, j   , 3) - ptv(i, j)) &
                                         *scp2(i, j    ), &
                                         (pbv(i, j, n) - p(i, j - 1, kintr)) &
                                         *scp2(i, j - 1)))
            mfl(kmin + 1) = min(fhi, max(flo, mfl(kmin + 1)))
            do k = kmin + 1, kmax - 1
               if     (mfl(k + 1) - mfl(k) > &
                        ffac*max(epsil, dlm(k))*scp2(i, j - 1)) then
                  mfl(k + 1) = mfl(k) + fface*dlm(k)*scp2(i, j - 1)
               elseif (mfl(k + 1) - mfl(k) < &
                      - ffac*max(epsil, dlp(k))*scp2(i, j    )) then
                  mfl(k + 1) = mfl(k) - fface*dlp(k)*scp2(i, j    )
               else
                  exit
               endif
            enddo

            ! Apply an iterative procedure for flux limiting by alternate upward
            ! and downward propagation through the layers.

            changed = .true.
            niter = 0
            kdir = 1

            do while (changed)

               niter = niter + 1
               if (niter == 1000) then
                  k = kmin
                  write(lp,*)
                  write(lp,'(i3,3e16.8)') &
                     1, mfl(k + 1), mfl(k), &
                     (mfl(k + 1) - mfl(k)) &
                     /(max(onemm, dpv(i, j, 1 + nn) + dpv(i, j, 2 + nn)) &
                       *delt1*scvx(i, j))
                  do k = kintr, kmax
                     kn = k + nn
                     write(lp,'(i3,3e16.8)') &
                        k, mfl(k + 1), mfl(k), &
                        (mfl(k + 1) - mfl(k)) &
                        /(max(onemm, dpv(i, j, kn))*delt1*scvx(i, j))
                  enddo
                  write(lp,*) 'no convergence v', i + i0, j + j0
                  call xchalt('(eddtra_gm)')
                         stop '(eddtra_gm)'
               endif

               changed = .false.
               kdir = - kdir

               do k = ((1 - kdir)*kmax + (1 + kdir)*kmin)/2, &
                      ((1 - kdir)*kmin + (1 + kdir)*kmax)/2, kdir

                  ! Proceed with flux limiting of this layer if the mass flux
                  ! difference between lower and upper interface is beyond the
                  ! floating point accuracy limitation.
                  if (abs(mfl(k + 1) - mfl(k)) > &
                      1.e-14_r8*max(epsil*scv2(i, j), &
                                    abs(mfl(k + 1) + mfl(k)))) then

                     if     (mfl(k + 1) - mfl(k) > &
                              ffac*max(epsil, dlm(k))*scp2(i, j - 1)) then
                        ! In this case, the mass fluxes are removing too much
                        ! mass from the grid cell at (i, j - 1, k). Limit the
                        ! dominating interface flux.
                        q = fface*dlm(k)*scp2(i, j - 1)
                        if (mfl(k + 1) > - mfl(k)) then
                           if (mfl(k    ) > - .5_r8*q) then
                              mfl(k + 1) =   mfl(k    ) + q
                           else
                              mfl(k + 1) =   .5_r8*q
                              mfl(k    ) = - mfl(k + 1)
                           endif
                        else
                           if (mfl(k + 1) <   .5_r8*q) then
                              mfl(k    ) =   mfl(k + 1) - q
                           else
                              mfl(k    ) = - .5_r8*q
                              mfl(k + 1) = - mfl(k    )
                           endif
                        endif
                        changed = .true.
                     elseif (mfl(k + 1) - mfl(k) < &
                            - ffac*max(epsil, dlp(k))*scp2(i, j    )) then
                        ! In this case, the mass fluxes are removing too much
                        ! mass from the grid cell at (i, j, k). Limit the
                        ! dominating interface flux.
                        q = fface*dlp(k)*scp2(i, j    )
                        if (mfl(k + 1) < - mfl(k)) then
                           if (mfl(k    ) <   .5_r8*q) then
                              mfl(k + 1) =   mfl(k    ) - q
                           else
                              mfl(k + 1) = - .5_r8*q
                              mfl(k    ) = - mfl(k + 1)
                           endif
                        else
                           if (mfl(k + 1) > - .5_r8*q) then
                              mfl(k    ) =   mfl(k + 1) + q
                           else
                              mfl(k    ) =   .5_r8*q
                              mfl(k + 1) = - mfl(k    )
                           endif
                        endif
                        changed = .true.
                     endif
                  endif

               enddo

            enddo

            ! ------------------------------------------------------------------
            ! Compute the final mass fluxes.
            ! ------------------------------------------------------------------

            k = kmin
            if (abs(mfl(k + 1) - mfl(k)) > &
                1.e-14_r8*max(epsil*scv2(i, j), &
                              abs(mfl(k + 1) + mfl(k)))) then
               vmfltd(i, j, 2 + mm) = mfl(k + 1) - mfl(k)
               vmfltd(i, j, 1 + mm) = vmfltd(i, j, 2 + mm) &
                                      *dpv(i, j, 1 + nn)/( dpv(i, j, 1 + nn) &
                                                         + dpv(i, j, 2 + nn))
               vmfltd(i, j, 2 + mm) = vmfltd(i, j, 2 + mm) &
                                    - vmfltd(i, j, 1 + mm)
            else
               vmfltd(i, j, 1 + mm) = 0._r8
               vmfltd(i, j, 2 + mm) = 0._r8
            endif
            do k = kintr, kmax
               km = k + mm
               if (abs(mfl(k + 1) - mfl(k)) > &
                   1.e-14_r8*max(epsil*scv2(i, j), &
                                 abs(mfl(k + 1) + mfl(k)))) then
                  vmfltd(i, j, km) = mfl(k + 1) - mfl(k)
               else
                  vmfltd(i, j, km) = 0._r8
               endif
               if (vmfltd(i, j, km) > &
                    ffac*max(epsil, dlm(k))*scp2(i, j - 1)) then
                  write (lp,*) 'eddtra_gm v >', &
                               i + i0, j + j0, k, vmfltd(i, j, km), &
                               ffac*max(epsil, dlm(k))*scp2(i, j - 1)
                  call xchalt('(eddtra_gm)')
                         stop '(eddtra_gm)'
               endif
               if (vmfltd(i, j, km) < &
                  - ffac*max(epsil, dlp(k))*scp2(i, j    )) then
                  write (lp,*) 'eddtra_gm v <', &
                               i + i0, j + j0, k, vmfltd(i, j, km), &
                             - ffac*max(epsil, dlp(k))*scp2(i, j    )
                  call xchalt('(eddtra_gm)')
                         stop '(eddtra_gm)'
               endif
            enddo

         enddo
         enddo
      enddo
   !$omp end parallel do

   end subroutine eddtra_gm

   ! ---------------------------------------------------------------------------
   ! Public procedures.
   ! ---------------------------------------------------------------------------

   subroutine eddtra(m, n, mm, nn, k1m, k1n)
   ! ---------------------------------------------------------------------------
   ! Compute eddy-induced transport.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      integer :: i, j, k, l, km

      ! Compute eddy-induced transport of mass.
      if     (eitmth == 'intdif') then
         call eddtra_intdif(m, n, mm, nn, k1m, k1n)
      elseif (eitmth == 'gm') then
         call eddtra_gm(m, n, mm, nn, k1m, k1n)
      else
         if (mnproc == 1) then
            write(lp,'(3a)') ' eitmth=', trim(eitmth), ' is unsupported!'
         endif
         call xcstop('(eddtra)')
                stop '(eddtra)'
      endif

      ! Diagnose eddy-induced transport components of heat and salt.
   !$omp parallel do private(k,km,l,i)
      do j = 1, jj
         do k = 1, kk
            km = k + mm
            do l = 1, isu(j)
            do i = max(1, ifu(j, l)), min(ii, ilu(j, l))
               utfltd(i, j, km) = .5_r8*umfltd(i, j, km) &
                                  *(temp(i - 1, j, km) + temp(i, j, km))
               usfltd(i, j, km) = .5_r8*umfltd(i, j, km) &
                                  *(saln(i - 1, j, km) + saln(i, j, km))
            enddo
            enddo
            do l = 1, isv(j)
            do i = max(1, ifv(j, l)), min(ii, ilv(j, l))
               vtfltd(i, j, km) = .5_r8*vmfltd(i, j, km) &
                                  *(temp(i, j - 1, km) + temp(i, j, km))
               vsfltd(i, j, km) = .5_r8*vmfltd(i, j, km) &
                                  *(saln(i, j - 1, km) + saln(i, j, km))
            enddo
            enddo
         enddo
      enddo
   !$omp end parallel do

      if (csdiag) then
         if (mnproc == 1) then
            write(lp,*) 'eddtra:'
         endif
         call chksummsk(umfltd(1 - nbdy, 1 - nbdy, k1m), iu, kk, 'umfltd')
         call chksummsk(vmfltd(1 - nbdy, 1 - nbdy, k1m), iv, kk, 'vmfltd')
         call chksummsk(utfltd(1 - nbdy, 1 - nbdy, k1m), iu, kk, 'utfltd')
         call chksummsk(vtfltd(1 - nbdy, 1 - nbdy, k1m), iv, kk, 'vtfltd')
         call chksummsk(usfltd(1 - nbdy, 1 - nbdy, k1m), iu, kk, 'usfltd')
         call chksummsk(vsfltd(1 - nbdy, 1 - nbdy, k1m), iv, kk, 'vsfltd')
      endif

   end subroutine eddtra

end module mod_eddtra
