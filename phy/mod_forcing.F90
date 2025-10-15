! ------------------------------------------------------------------------------
! Copyright (C) 2002-2025 Mats Bentsen, Jerry Tjiputra, Jörg Schwinger,
!                         Mariana Vertenstein, Joeran Maerz
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
   use mod_tracers, only: ntr
   use mod_ifdefs, only: use_TRC
   use mod_checksum, only: csdiag, chksum
   use mod_utility, only: fnmlen

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

   character(len = fnmlen) :: &
      scfile, &       ! Name of file containing monthly SSS climatology.
      wavsrc          ! Source of wave fields. Valid source: 'none', 'param',
                      ! 'extern'.

   ! Options derived from string options.
   integer :: &
      wavsrc_opt

   ! Parameters:
   integer, parameter :: &
      ! Wave source options:
      wavsrc_none        = 0, & ! No wave fields.
      wavsrc_param       = 1, & ! Parameterized wave fields.
      wavsrc_extern      = 2    ! Receive external wave fields.

   ! Constants used in forcing computations.
   real(r8) :: &
      sref = 34.65_r8 ! Global reference sea-surface salinity [g kg-1].

   ! Variables for diagnosed relaxation fluxes.
   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, 48) :: &
      tflxap, &       ! Heat flux to be applied [W m-2].
      sflxap, &       ! Salt flux to be applied [g m-2 s-1].
      tflxdi, &       ! Diagnosed heat flux [W m-2].
      sflxdi          ! Diagnosed salt flux [g m-2 s-1].
   integer, dimension(48) :: &
      nflxdi          ! Accumulation counter for diagnosed relaxation fluxes.

   ! Monthly climatological fields used in the computation of climatological
   ! fluxes and relaxation of sea surface temperature and salinity.
   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, 12) :: &
      sstclm, &       ! Sea-surface temperature [deg C].
      ricclm, &       ! Sea-ice concentration [].
      sssclm          ! Sea-surface salinity [g kg-1].

   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy) :: &
      sst_stream, &    ! Sea-surface temperature [deg C] from stream data.
      ice_stream, &    ! Sea-ice concentration [] from stream data.
      sss_stream       ! Sea-surface salinity [g kg-1] from stream data.

   logical :: use_stream_relaxation ! If true, use nuopc stream relaxation capability

   ! Allocate dust_stream in ocn_stream_dust.F90
   real(r8), allocatable  :: dust_stream(:,:,:) ! iron dust deposition flux (hamocc)
   logical :: use_stream_dust    ! If true, use nuopc stream dust capability (hamocc only)

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
      hmat, &         ! Surface material enthalpy flux [W m-2].
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
      lamult, &       ! Langmuir enhancement factor [].
      lasl, &         ! Surface layer averaged Langmuir number [].
      ustokes, &      ! u-component of surface Stokes drift [m s-1].
      vstokes, &      ! v-component of surface Stokes drift [m s-1].
      atmco2, &       ! Atmospheric CO2 concentration [ppm].
      flxco2, &       ! Air-sea CO2 flux [kg m-2 s-1].
      flxdms, &       ! Sea-air DMS flux [kg m-2 s-1].
      flxbrf, &       ! sea-air bromoform flux
      atmbrf, &       ! atmospheric bromoform concentration
      flxn2o, &       ! sea-air nitrous oxide flux [kg N2O m-2 s-1]
      atmn2o, &       ! atmospheric nitrous oxide concentration [pptv]
      flxnh3, &       ! sea-air ammonia flux [kg NH3 m-2 s-1]
      atmnh3, &       ! atmospheric ammonia concentration [pptv]
      atmnhxdep,&     ! atmospheric nhx deposition [kgN/m2/s]
      atmnoydep       ! atmospheric noy deposition [kgN/m2/s]

   real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy) :: &
      surflx, &       ! Surface thermal energy flux [W m-2].
      surrlx, &       ! Surface relaxation thermal energy flux [W m-2].
      sswflx, &       ! Surface solar energy flux [W m-2].
      salflx, &       ! Surface salinity flux [g m-2 s-1].
      brnflx, &       ! Surface brine flux [g m-2 s-1].
      salrlx, &       ! Surface relaxation salinity flux [g m-2 s-1].
      taux, &         ! u-component of surface stress [kg m-1 s-2].
      tauy, &         ! v-component of surface stress [kg m-1 s-2].
      ustar, &        ! Surface friction velocity [m s-1].
      ustarb, &       ! Bottom friction velocity [m s-1].
      ustar3, &       ! Friction velocity cubed [m3 s-3].
      wstar3          ! Convective velocity cubed [m3 s-3].

   ! Correction fields for keeping tracers above zero value (expectation that
   ! this will only occure due to the application of tracer virtual surface
   ! fluxes.

   real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
      salt_corr       ! Correction of salt [g m-2].
   real(r8), allocatable, dimension(:,:,:) :: &
      trc_corr        ! Correction of tracers (unit depending on tracer).

   ! Flux fields at model interfaces.

   real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy, kk + 1) :: &
      buoyfl, &       ! Buoyancy flux [m2 s-3].
      t_sw_nonloc, &  ! Non-local transport term that is the fraction of
                      ! shortwave flux passing a layer interface [].
      t_rs_nonloc, &  ! Non-local transport term that is the fraction of
                      ! restoring heat flux passing a layer interface [].
      s_br_nonloc, &  ! Non-local transport term that is the fraction of
                      ! brine flux passing a layer interface [].
      s_rs_nonloc     ! Non-local transport term that is the fraction of
                      ! restoring salt flux passing a layer interface [].

   public :: aptflx, apsflx, ditflx, disflx, srxbal, sprfac, &
             trxday, srxday, trxdpt, srxdpt, trxlim, srxlim, scfile, &
             wavsrc, wavsrc_opt, wavsrc_none, wavsrc_param, wavsrc_extern, &
             sref, tflxap, sflxap, tflxdi, sflxdi, nflxdi, &
             sstclm, ricclm, sssclm, prfac, eiacc, pracc, &
             swa, nsf, hmltfz, hmat, lip, sop, eva, rnf, rfi, fmltfz, sfl, &
             ztx, mty, ustarw, slp, abswnd, lamult, lasl, ustokes, vstokes, &
             atmco2, flxco2, flxdms, flxbrf, atmbrf, &
             atmn2o,flxn2o,atmnh3,flxnh3, atmnhxdep,atmnoydep, &
             surflx, surrlx, sswflx, salflx, brnflx, salrlx, taux, tauy, &
             ustar, ustarb, ustar3, wstar3, buoyfl, salt_corr, trc_corr, &
             t_sw_nonloc, t_rs_nonloc, s_br_nonloc, s_rs_nonloc, &
             inivar_forcing, fwbbal, sss_stream, sst_stream, ice_stream, &
             dust_stream, &
             use_stream_relaxation, use_stream_dust

contains

   subroutine inivar_forcing
   ! ---------------------------------------------------------------------------
   ! Initialize variables related to forcing.
   ! ---------------------------------------------------------------------------

      integer :: i, j, k, l, errstat, nt

      eiacc(:,:) = spval
      pracc(:,:) = spval
      swa(:,:) = spval
      nsf(:,:) = spval
      hmltfz(:,:) = spval
      hmat(:,:) = spval
      lip(:,:) = spval
      sop(:,:) = spval
      eva(:,:) = spval
      rnf(:,:) = spval
      rfi(:,:) = spval
      fmltfz(:,:) = spval
      sfl(:,:) = spval
      ztx(:,:) = spval
      mty(:,:) = spval
      ustarw(:,:) = spval
      slp(:,:) = spval
      abswnd(:,:) = spval
      lamult(:,:) = spval
      lasl(:,:) = spval
      ustokes(:,:) = spval
      vstokes(:,:) = spval
      atmco2(:,:) = spval
      flxco2(:,:) = spval
      flxdms(:,:) = spval
      atmbrf(:,:) = spval
      flxbrf(:,:) = spval
      atmn2o(:,:) = -spval
      flxn2o(:,:) = spval
      atmnh3(:,:) = -spval
      flxnh3(:,:) = spval
      atmnhxdep(:,:) = spval
      atmnoydep(:,:) = spval
      surflx(:,:) = spval
      surrlx(:,:) = spval
      sswflx(:,:) = spval
      salflx(:,:) = spval
      brnflx(:,:) = spval
      salrlx(:,:) = spval
      taux(:,:) = spval
      tauy(:,:) = spval
      ustar(:,:) = spval
      ustarb(:,:) = spval
      ustar3(:,:) = spval
      wstar3(:,:) = spval
      salt_corr(:,:) = spval
      buoyfl(:,:,:) = spval
      t_sw_nonloc(:,:,:) = spval
      t_rs_nonloc(:,:,:) = spval
      s_br_nonloc(:,:,:) = spval
      s_rs_nonloc(:,:,:) = spval

      if (use_TRC) then
         allocate(trc_corr(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ntr), stat = errstat)
         if (errstat /= 0) then
            write(lp,*) 'Failed to allocate trc_corr!'
            call xchalt('(inivar_forcing)')
            stop '(inivar_forcing)'
         endif
         trc_corr(:,:,:) = spval
      endif

   !$omp parallel do private(l, i, k)
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
            t_rs_nonloc(i, j, 1) = 1._r8
            s_rs_nonloc(i, j, 1) = 1._r8
         enddo
         enddo
         do k = 2, kk+1
           do l = 1, isp(j)
           do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
              t_rs_nonloc(i, j, k) = 0._r8
              s_rs_nonloc(i, j, k) = 0._r8
           enddo
           enddo
         enddo
      enddo
   !$omp end parallel do
   !$omp parallel do private(l, i)
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
            flxco2(i, j) = 0._r8
            flxdms(i, j) = 0._r8
            flxbrf(i, j) = 0._r8
            flxn2o(i, j) = 0._r8
            flxnh3(i, j) = 0._r8
            atmnhxdep(i, j) = 0._r8
            atmnoydep(i, j) = 0._r8
            ustar (i, j) = 0._r8
            ustarb(i, j) = 0._r8
            wstar3(i, j) = 0._r8
            salt_corr(i, j) = 0._r8
         enddo
         enddo
      enddo
   !$omp end parallel do

      if (use_TRC) then
      !$omp parallel do private(nt, l, i)
         do j = 1, jj
            do nt = 1, ntr
               do l = 1, isp(j)
               do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
                  trc_corr(i, j, nt) = 0._r8
               enddo
               enddo
            enddo
         enddo
      !$omp end parallel do
      endif

   !$omp parallel do private(k, l, i)
      do j = 1, jj
         do k = 1, kk + 1
            do l = 1, isp(j)
            do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
               buoyfl(i, j, k) = 0._r8
            enddo
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
         call chksum(eiacc, 1, halo_ps, 'eiacc')
         call chksum(pracc, 1, halo_ps, 'pracc')
      endif

   end subroutine fwbbal

end module mod_forcing
