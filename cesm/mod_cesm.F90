! ------------------------------------------------------------------------------
! Copyright (C) 2011-2022 Mats Bentsen, Jerry Tjiputra, Jörg Schwinger
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

module mod_cesm
! ------------------------------------------------------------------------------
! This module contains variables and routines related to the coupling to CESM
! that must be available to BLOM routines.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_constants, only: pi
   use mod_time, only: nstep
   use mod_xc
   use mod_forcing, only: trxday, srxday, swa, nsf, lip, sop, eva, rnf, rfi, &
                          fmltfz, sfl, ztx, mty, ustarw, slp, abswnd, &
                          lamult, lasl, ustokes, vstokes, atmco2, atmbrf, flxdms
   use mod_ben02, only: initai, rdcsic, rdctsf, fnlzai
   use mod_seaice, only: ficem
   use mod_checksum, only: csdiag, chksummsk

   implicit none

   private

   character(len = 256) :: &
      runid_cesm, &      ! Case name received from CESM.
      runtyp_cesm        ! Run type received from CESM.
   integer :: &
      ocn_cpl_dt_cesm, & ! Coupling time interval.
      nstep_in_cpl       ! Number of ocean time steps in a coupling interval.

   ! Heat flux due to melting received from CESM amd freezing and melting
   ! potentials sent to CESM.
   real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy) :: &
      hmlt, &            ! Heat flux due to melting [W m-2].
      frzpot, &          ! Freezing potential [J m-2].
      mltpot             ! Melting potential [J m-2].

   ! Forcing arrays, keeping two forcing intervals for time smoothing.
   real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy, 2) :: &
      swa_da, &          ! Solar heat flux [W m-2].
      nsf_da, &          ! Non-solar heat flux [W m-2].
      hmlt_da, &         ! Heat flux due to melting [W m-2].
      lip_da, &          ! Liquid water flux [kg m-2 s-1].
      sop_da, &          ! Solid precipitation [kg m-2 s-1].
      eva_da, &          ! Evaporation [kg m-2 s-1].
      rnf_da, &          ! Liquid runoff [kg m-2 s-1].
      rfi_da, &          ! Frozen Runoff [kg m-2 s-1].
      fmltfz_da, &       ! Fresh water flux due to melting and freezing
                         ! [kg m-2 s-1].
      sfl_da, &          ! Salt flux [kg m-2 s-1].
      ztx_da, &          ! u-component of wind stress [kg m-1 s-2].
      mty_da, &          ! v-component of wind stress [kg m-1 s-2].
      ustarw_da, &       ! Friction velocity for open water [m s-1].
      slp_da, &          ! Sea-level pressure [kg m-1 s-2].
      abswnd_da, &       ! Wind speed at measurement height (zu) [m s-1].
      ficem_da, &        ! Ice concentration [].
      lamult_da, &       ! Langmuir enhancement factor [].
      lasl_da, &         ! Surface layer averaged Langmuir number [].
      ustokes_da, &      ! u-component of surface Stokes drift [m s-1].
      vstokes_da, &      ! v-component of surface Stokes drift [m s-1].
      atmco2_da, &       ! Atmospheric CO2 concentration [ppm].
      atmbrf_da, &       ! Atmospheric bromoform concentration [ppt].
      flxdms_da          ! dms surface flux computed by mediator [kg m-2 s-1]

   logical :: &
      smtfrc             ! If true, time smooth CESM forcing fields.

   integer :: &
      l1ci, l2ci         ! Time-level indices for time smoothing of CESM fields.

   logical :: get_flxdms_from_med

   public :: runid_cesm, runtyp_cesm, ocn_cpl_dt_cesm, nstep_in_cpl, hmlt, &
             frzpot, mltpot, swa_da, nsf_da, hmlt_da, lip_da, sop_da, eva_da, &
             rnf_da, rfi_da, fmltfz_da, sfl_da, ztx_da, mty_da, ustarw_da, &
             slp_da, abswnd_da, ficem_da, lamult_da, lasl_da, flxdms_da, &
             ustokes_da, vstokes_da, atmco2_da, atmbrf_da, smtfrc, l1ci, l2ci, &
             inicon_cesm, inifrc_cesm, getfrc_cesm, get_flxdms_from_med
contains

   subroutine inicon_cesm
   ! ---------------------------------------------------------------------------
   ! Set initial conditions for variables specifically when coupled to CESM.
   ! ---------------------------------------------------------------------------

      integer :: i, j

   !$omp parallel do private(i)
      do j = 1, jj
         do i = 1, ii
            frzpot(i, j) = 0._r8
            mltpot(i, j) = 0._r8
         enddo
      enddo
   !$omp end parallel do

   end subroutine inicon_cesm

   subroutine inifrc_cesm
   ! ---------------------------------------------------------------------------
   ! Initialize climatological fields for surface restoring and interpolation of
   ! CESM forcing fields.
   ! ---------------------------------------------------------------------------

      ! If SST restoring is requested, prepare interpolation of surface fields
      ! and read climatological sea-ice concentration and surface temperature.
      if (trxday > 0._r8) then
        call initai
        call rdcsic
        call rdctsf
      endif

      ! If SSS restoring is requested, read climatological sea surface salinity.
      if (srxday > 0._r8) call rdcsss

      ! Initialize diagnosing/application of relaxation fluxes.
      call idarlx

      ! Deallocate memory used for interpolation of surface fields.
      if (trxday > 0._r8) then
        call fnlzai
      endif

      ! Initialize time level indexes
      l1ci = 1
      l2ci = 1

   end subroutine inifrc_cesm

   subroutine getfrc_cesm
   ! ---------------------------------------------------------------------------
   ! Interpolate CESM forcing fields.
   ! ---------------------------------------------------------------------------

#define DIAG
#undef DIAG
#ifdef DIAG
      use mod_nctools
      use mod_dia, only : iotype
#endif

      integer :: i, j, l
      real(r8) :: w1, w2

      if (smtfrc) then
         w1 = .5_r8*( 1._r8 &
                    + cos((mod(nstep - 1, nstep_in_cpl) + 1)*pi/nstep_in_cpl))
      else
         w1 = 0._r8
      endif
      w2 = 1._r8 - w1

   !$omp parallel do private(l, i)
      do j = 1, jj
        do l = 1, isp(j)
        do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
           ustarw(i, j)  = w1*ustarw_da(i, j, l1ci)  + w2*ustarw_da(i, j, l2ci)
           lip(i, j)     = w1*lip_da(i, j, l1ci)     + w2*lip_da(i, j, l2ci)
           sop(i, j)     = w1*sop_da(i, j, l1ci)     + w2*sop_da(i, j, l2ci)
           eva(i, j)     = w1*eva_da(i, j, l1ci)     + w2*eva_da(i, j, l2ci)
           rnf(i, j)     = w1*rnf_da(i, j, l1ci)     + w2*rnf_da(i, j, l2ci)
           rfi(i, j)     = w1*rfi_da(i, j, l1ci)     + w2*rfi_da(i, j, l2ci)
           fmltfz(i, j)  = w1*fmltfz_da(i, j, l1ci)  + w2*fmltfz_da(i, j, l2ci)
           sfl(i, j)     = w1*sfl_da(i, j, l1ci)     + w2*sfl_da(i, j, l2ci)
           swa(i, j)     = w1*swa_da(i, j, l1ci)     + w2*swa_da(i, j, l2ci)
           nsf(i, j)     = w1*nsf_da(i, j, l1ci)     + w2*nsf_da(i, j, l2ci)
           hmlt(i, j)    = w1*hmlt_da(i, j, l1ci)    + w2*hmlt_da(i, j, l2ci)
           slp(i, j)     = w1*slp_da(i, j, l1ci)     + w2*slp_da(i, j, l2ci)
           abswnd(i, j)  = w1*abswnd_da(i, j, l1ci)  + w2*abswnd_da(i, j, l2ci)
           ficem(i, j)   = w1*ficem_da(i, j, l1ci)   + w2*ficem_da(i, j, l2ci)
           lamult(i, j)  = w1*lamult_da(i, j, l1ci)  + w2*lamult_da(i, j, l2ci)
           lasl(i, j)    = w1*lasl_da(i, j, l1ci)    + w2*lasl_da(i, j, l2ci)
           ustokes(i, j) = w1*ustokes_da(i, j, l1ci) + w2*ustokes_da(i, j, l2ci)
           vstokes(i, j) = w1*vstokes_da(i, j, l1ci) + w2*vstokes_da(i, j, l2ci)
           atmco2(i, j)  = w1*atmco2_da(i, j, l1ci)  + w2*atmco2_da(i, j, l2ci)
           atmbrf(i, j)  = w1*atmbrf_da(i, j, l1ci)  + w2*atmbrf_da(i, j, l2ci)
           if (get_flxdms_from_med) then
              flxdms(i, j)  = w1*flxdms_da(i, j, l1ci)  + w2*flxdms_da(i, j, l2ci) 
           end if
        enddo
        enddo
        do l = 1, isu(j)
        do i = max(1, ifu(j, l)), min(ii, ilu(j, l))
           ztx(i, j) = w1*ztx_da(i, j, l1ci) + w2*ztx_da(i, j, l2ci)
        enddo
        enddo
        do l = 1, isv(j)
        do i = max(1, ifv(j, l)), min(ii, ilv(j, l))
           mty(i, j) = w1*mty_da(i, j, l1ci) + w2*mty_da(i, j, l2ci)
        enddo
        enddo
      enddo
   !$omp end parallel do

#ifdef DIAG
      call ncfopn('getfrc_cesm.nc', 'w', 'c', 1, iotype)
      call ncdims('x', itdm)
      call ncdims('y', jtdm)
      call ncdefvar('ustarw_da', 'x y', ndouble, 8) 
      call ncdefvar('lip_da', 'x y', ndouble, 8)
      call ncdefvar('sop_da', 'x y', ndouble, 8)
      call ncdefvar('eva_da', 'x y', ndouble, 8)
      call ncdefvar('rnf_da', 'x y', ndouble, 8)
      call ncdefvar('rfi_da', 'x y', ndouble, 8)
      call ncdefvar('fmltfz_da', 'x y', ndouble, 8)
      call ncdefvar('sfl_da', 'x y', ndouble, 8)
      call ncdefvar('swa_da', 'x y', ndouble, 8)
      call ncdefvar('nsf_da', 'x y', ndouble, 8)
      call ncdefvar('hmlt_da', 'x y', ndouble, 8)
      call ncdefvar('slp_da', 'x y', ndouble, 8)
      call ncdefvar('abswnd_da', 'x y', ndouble, 8)
      call ncdefvar('ficem_da', 'x y', ndouble, 8)
      call ncdefvar('lamult_da', 'x y', ndouble, 8)
      call ncdefvar('lasl_da', 'x y', ndouble, 8)
      call ncdefvar('ustokes_da', 'x y', ndouble, 8)
      call ncdefvar('vstokes_da', 'x y', ndouble, 8)
      call ncdefvar('atmco2_da', 'x y', ndouble, 8)
      call ncdefvar('atmbrf_da', 'x y', ndouble, 8)
      call ncdefvar('ztx_da', 'x y', ndouble, 8)
      call ncdefvar('mty_da', 'x y', ndouble, 8)
      call ncedef

      call ncwrtr('ustarw_da', 'x y', ustarw_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('lip_da', 'x y', lip_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('sop_da', 'x y', sop_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('eva_da', 'x y', eva_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('rnf_da', 'x y', rnf_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('rfi_da', 'x y', rfi_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('fmltfz_da', 'x y', fmltfz_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('sfl_da', 'x y', sfl_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('swa_da', 'x y', swa_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('nsf_da', 'x y', nsf_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('hmlt_da', 'x y', hmlt_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('slp_da', 'x y', slp_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('abswnd_da', 'x y', abswnd_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('ficem_da', 'x y', ficem_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('lamult_da', 'x y', lamult_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('lasl_da', 'x y', lasl_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('ustokes_da', 'x y', ustokes_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('vstokes_da', 'x y', vstokes_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('atmco2_da', 'x y', atmco2_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('atmbrf_da', 'x y', atmbrf_da(1 - nbdy, 1 - nbdy, l2ci), &
                  ip, 1, 1._r8, 0._r8, 8)
      call ncwrtr('ztx_da', 'x y', ztx_da(1 - nbdy, 1 - nbdy, l2ci), &
                  iu, 1, 1._r8, 0._r8, 8)
      call ncwrtr('mty_da', 'x y', mty_da(1 - nbdy, 1 - nbdy, l2ci), &
                  iv, 1, 1._r8, 0._r8, 8)
      call ncfcls
      call xcstop('(getfrc_cesm)')
             stop '(getfrc_cesm)'
#endif

      if (csdiag) then
         if (mnproc == 1) then
            write (lp, *) 'getfrc_cesm:'
         endif
         call chksummsk(ustarw, ip, 1, 'ustarw')
         call chksummsk(ztx, iu, 1, 'ztx')
         call chksummsk(mty, iv, 1, 'mty')
         call chksummsk(lip, ip, 1, 'lip')
         call chksummsk(sop, ip, 1, 'sop')
         call chksummsk(eva, ip, 1, 'eva')
         call chksummsk(rnf, ip, 1, 'rnf')
         call chksummsk(rfi, ip, 1, 'rfi')
         call chksummsk(fmltfz, ip, 1, 'fmltfz')
         call chksummsk(sfl, ip, 1, 'sfl')
         call chksummsk(swa, ip, 1, 'swa')
         call chksummsk(nsf, ip, 1, 'nsf')
         call chksummsk(hmlt, ip, 1, 'hmlt')
         call chksummsk(slp, ip, 1, 'slp')
         call chksummsk(abswnd, ip, 1, 'abswnd')
         call chksummsk(ficem, ip, 1, 'ficem')
         call chksummsk(lamult, ip, 1, 'lamult')
         call chksummsk(lasl, ip, 1, 'lasl')
         call chksummsk(ustokes, ip, 1, 'ustokes')
         call chksummsk(vstokes, ip, 1, 'vstokes')
         call chksummsk(atmco2, ip, 1, 'atmco2')
         call chksummsk(atmbrf, ip, 1, 'atmbrf')
      endif

   end subroutine getfrc_cesm

end module mod_cesm
