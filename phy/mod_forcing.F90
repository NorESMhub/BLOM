! ------------------------------------------------------------------------------
! Copyright (C) 2002-2020 Mats Bentsen, Jerry Tjiputra
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

module mod_forcing
! ------------------------------------------------------------------------------
! This module contains variables related to the application forcing fields that
! is shared among the various sources of forcing.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_constants, only: spval
   use mod_time, only: nday_of_year, nstep, nstep_in_day, baclin
   use mod_xc
   use mod_grid, only: scp2
   use mod_checksum, only: csdiag, chksummsk

   implicit none

   private

   ! Variables to be set in namelist:
   logical :: &
      aptflx, &       ! If true, apply diagnosed heat flux.
      apsflx, &       ! If true, apply diagnosed freshwater flux.
      ditflx, &       ! If true, diagnose heat flux.
      disflx, &       ! If true, diagnose freshwater flux.
      srxbal, &       ! If true, globally balance the SSS relaxation flux.
      sprfac          ! If true, apply factor to precipitation and runoff for
                      ! balancing the freshwater forcing budget. In case of
                      ! coupling to CESM, this implies sending the factor to the
                      ! coupler for application.
   real(r8) :: &
      trxday, &       ! e-folding relaxation time scale for SST [days].
      srxday, &       ! e-folding relaxation time scale for SSS [days].
      trxdpt, &       ! Maximum mixed layer depth nudged by SST relaxation [m].
      srxdpt, &       ! Maximum mixed layer depth nudged by SSS relaxation [m].
      trxlim, &       ! Maximum absolute value of SST difference in relaxation
                      ! [deg C].
      srxlim          ! Maximum absolute value of SSS difference in relaxation
                      ! [g kg-1].
   character(len = 256) :: &
      scfile          ! Name of file containing monthly SSS climatology.


   ! Constants used in forcing computations.
   real(r8) :: &
      sref = 34.65_r8 ! Global reference sea-surface salinity [g kg-1].

   ! Variables for diagnosed relaxation fluxes.
   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, 48) :: &
      tflxap, &       ! Heat flux to be applied [W cm-2].
      sflxap, &       ! Salt flux to be applied [10e-3 g cm-2 s-1].
      tflxdi, &       ! Diagnosed heat flux [W cm-2].
      sflxdi          ! Diagnosed salt flux [10e-3 g cm-2 s-1].
   integer, dimension(48) :: &
      nflxdi          ! Accumulation counter for diagnosed relaxation fluxes.

   ! Monthly climatological fields used in the computation of climatological
   ! fluxes and relaxation of sea surface temperature and salinity.
   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, 12) :: &
      sstclm, &       ! Sea-surface temperature [deg C].
      ricclm, &       ! Sea-ice concentration [].
      sssclm          ! Sea-surface salinity [g kg-1].

   ! Variables related to balancing the freshwater forcing budget.
   real(r8) :: &
      prfac           ! Correction factor for precipitation and runoff [].
   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy) :: &
      eiacc, &        ! Accumulation of freshwater fluxes related to
                      ! evaporation and sea-ice melt ing and freezing.
      pracc           ! Accumulation of freshwaterfluxes related to
                      ! precipitation and runoff.

   ! Various surface state and flux fields.

   real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy) :: &
      swa, &          ! Solar heat flux [W m-2].
      nsf, &          ! Non-solar heat flux [W m-2].
      hmltfz, &       ! Heat flux due to melting and freezing [W m-2].
      lip, &          ! Liquid water flux [kg m-2 s-1].
      sop, &          ! Solid precipitation [kg m-2 s-1].
      eva, &          ! Evaporation [kg m-2 s-1].
      rnf, &          ! Liquid runoff [kg m-2 s-1].
      rfi, &          ! Frozen runoff [kg m-2 s-1].
      fmltfz, &       ! Fresh water flux due to melting and freezing
                      ! [kg m-2 s-1].
      sfl, &          ! Salt flux [kg m-2 s-1].
      ztx, &          ! u-component of wind stress [kg m-1 s-2].
      mty, &          ! v-component of wind stress [kg m-1 s-2].
      ustarw, &       ! Friction velocity for open water [m s-1].
      slp, &          ! Sea-level pressure [kg m-1 s-2].
      abswnd, &       ! Wind speed at measurement height (zu) [m s-1].
      atmco2, &       ! Atmospheric CO2 concentration [ppm].
      flxco2, &       ! Air-sea CO2 flux [kg m-2 s-1].
      flxdms, &       ! Sea-air DMS flux [kg m-2 s-1].
      flxbrf, &       ! sea-air bromoform flux
      atmbrf          ! atmospheric bromoform concentration


   real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy) :: &
      surflx, &       ! Surface thermal energy flux [W cm-2].
      surrlx, &       ! Surface relaxation thermal energy flux [W cm-2].
      sswflx, &       ! Surface solar energy flux [W cm-2].
      salflx, &       ! Surface salinity flux [10e-3 g cm-2 s-1].
      brnflx, &       ! Surface brine flux [10e-3 g cm-2 s-1].
      salrlx, &       ! Surface relaxation salinity flux [10e-3 g cm-2 s-1].
      taux, &         ! u-component of surface stress [g cm-1 s-2].
      tauy, &         ! v-component of surface stress [g cm-1 s-2].
      ustar, &        ! Surface friction velocity [cm s-1].
      ustarb, &       ! Bottom friction velocity [cm s-1].
      ustar3, &       ! Friction velocity cubed [cm3 s-3].
      buoyfl          ! Surface buoyancy flux [cm2 s-3].

   public :: aptflx, apsflx, ditflx, disflx, srxbal, sprfac, &
             trxday, srxday, trxdpt, srxdpt, trxlim, srxlim, scfile, &
             sref, tflxap, sflxap, tflxdi, sflxdi, nflxdi, &
             sstclm, ricclm, sssclm, prfac, eiacc, pracc, &
             swa, nsf, hmltfz, lip, sop, eva, rnf, rfi, fmltfz, sfl, ztx, mty, &
             ustarw, slp, abswnd, atmco2, flxco2, flxdms, flxbrf, atmbrf, &
             surflx, surrlx, sswflx, salflx, brnflx, salrlx, taux, tauy, &
             ustar, ustarb, ustar3, buoyfl, &
             inivar_forcing, fwbbal

contains

   subroutine inivar_forcing
   ! ---------------------------------------------------------------------------
   ! Initialize variables related to forcing.
   ! ---------------------------------------------------------------------------

      integer :: i, j, l

   !$omp parallel do private(i)
      do j = 1 - nbdy, jj + nbdy
         do i = 1 - nbdy, ii + nbdy
            eiacc(i, j) = spval
            pracc(i, j) = spval
            swa(i, j) = spval
            nsf(i, j) = spval
            hmltfz(i, j) = spval
            lip(i, j) = spval
            sop(i, j) = spval
            eva(i, j) = spval
            rnf(i, j) = spval
            rfi(i, j) = spval
            fmltfz(i, j) = spval
            sfl(i, j) = spval
            ztx(i, j) = spval
            mty(i, j) = spval
            ustarw(i, j) = spval
            slp(i, j) = spval
            abswnd(i, j) = spval
            atmco2(i, j) = spval
            flxco2(i, j) = spval
            flxdms(i, j) = spval
            atmbrf(i, j) = spval
            flxbrf(i, j) = spval
            surflx(i, j) = spval
            surrlx(i, j) = spval
            sswflx(i, j) = spval
            salflx(i, j) = spval
            brnflx(i, j) = spval
            salrlx(i, j) = spval
            taux(i, j) = spval
            tauy(i, j) = spval
            ustar(i, j) = spval
            ustarb(i, j) = spval
            ustar3(i, j) = spval
            buoyfl(i, j) = spval
         enddo
      enddo
   !$omp end parallel do

   !$omp parallel do private(l, i)
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
            flxco2(i, j) = 0._r8
            flxdms(i, j) = 0._r8
            ustar (i, j) = 0._r8
            ustarb(i, j) = 0._r8
            buoyfl(i, j) = 0._r8
#ifdef BROMO
            flxbrf(i, j) = 0._r8
#endif
         enddo
         enddo
      enddo
   !$omp end parallel do

      if (sprfac) then
         prfac = 1._r8
      !$omp parallel do private(l, i)
         do j = 1, jj
            do l = 1, isp(j)
            do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
               eiacc(i, j) = 0._r8
               pracc(i, j) = 0._r8
            enddo
            enddo
         enddo
      !$omp end parallel do
      endif

   end subroutine inivar_forcing

   subroutine fwbbal(m, n, mm, nn, k1m, k1n)
   ! ---------------------------------------------------------------------------
   ! Balance the freshwater budget by computing a correcting factor to be
   ! applied to precipitation and runoff. The correction factor is based on the
   ! fresh water budget for the previous year.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      real(r8) :: totei, totpr
      integer :: i, j, l

      if (.not. sprfac) return

      ! Accumulate two groups of fresh water fluxes. One is evaporation and
      ! sea-ice melting/freezing and the other is precipitation and runoff. The
      ! fresh water fluxes are weighted with the time step in case it varies
      ! during the accumulation period.
   !$omp parallel do private(l, i)
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
            eiacc(i, j) = eiacc(i, j) &
                        + (eva(i, j) + fmltfz(i, j))*baclin
            pracc(i, j) = pracc(i, j) &
                        + (lip(i, j) + sop(i, j) + rnf(i, j) + rfi(i, j))*baclin
         enddo
         enddo
      enddo
   !$omp end parallel do

      ! Compute new correction factor at the end of a year and reset
      ! accumulation arrays.
      if (nday_of_year == 1 .and. mod(nstep, nstep_in_day) == 0) then

         ! Weight the accumulated fluxes with grid cell area and do global sums,
         ! but only including grid cells connected to the world ocean.
      !$omp parallel do private(l, i)
         do j = 1, jj
            do l = 1, isp(j)
            do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
               eiacc(i, j) = eiacc(i, j)*scp2(i, j)
               pracc(i, j) = pracc(i, j)*scp2(i, j)
            enddo
            enddo
         enddo
     !$omp end parallel do
         call xcsum(totei, eiacc, ipwocn)
         call xcsum(totpr, pracc, ipwocn)

         ! Update correction factor.
         prfac = - prfac*totei/totpr
         if (mnproc == 1) then
            write (lp, *) &
               'new correction factor for precipitation/runoff:', prfac
            call flush(lp)
         endif

         ! Reset accumulation arrays.
      !$omp parallel do private(l, i)
         do j = 1, jj
            do l = 1, isp(j)
            do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
               eiacc(i, j) = 0._r8
               pracc(i, j) = 0._r8
            enddo
            enddo
         enddo
      !$omp end parallel do

      endif

      if (csdiag) then
         if (mnproc == 1) then
            write (lp, *) 'fwbbal:'
         endif
         call chksummsk(eiacc, ip, 1, 'eiacc')
         call chksummsk(pracc, ip, 1, 'pracc')
      endif

   end subroutine fwbbal

end module mod_forcing
