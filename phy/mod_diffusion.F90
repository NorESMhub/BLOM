! ------------------------------------------------------------------------------
! Copyright (C) 2020-2025 Mats Bentsen, Mehmet Ilicak, Aleksi Nummelin
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

module mod_diffusion
! ------------------------------------------------------------------------------
! This module contains variables related to diffusion.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_config, only: inst_suffix
   use mod_constants, only: spval, epsilk
   use mod_xc
   use mod_forcing, only: wavsrc_opt, wavsrc_none, wavsrc_param, wavsrc_extern
   use mod_utility, only: fnmlen

   implicit none
   private

   ! Variables to be set by namelist.
   real(r8) :: &
      egc, &    ! The parameter 'c' in the Eden and Greatbatch (2008)
                ! parameterization [].
      eggam, &  ! The parameter 'gamma' in the Eden and Greatbatch (2008)
                ! parameterization [].
      eglsmn, & ! Minimum eddy length scale in the Eden and Greatbatch (2008)
                ! parameterization [m].
      egmndf, & ! Minimum diffusivity in the Eden and Greatbatch (2008)
                ! parameterization [m2 s-1].
      egmxdf, & ! Maximum diffusivity in the Eden and Greatbatch (2008)
                ! parameterization [m2 s-1].
      egidfq, & ! Factor relating the isopycnal diffusivity to the layer
                ! interface diffusivity in the Eden and Greatbatch (2008)
                ! parameterization (egidfq = difint/difiso) [].
      rhiscf, & ! Linear scaling parameter for topographic rhines scale [].
      ri0, &    ! Critical gradient richardson number for shear driven vertical
                ! mixing [].
      bdmc1, &  ! Background diapycnal diffusivity times buoyancy frequency
                ! [m2 s-2].
      bdmc2, &  ! Background diapycnal diffusivity [m2 s-1].
      tkepf     ! Fraction of surface TKE that penetrates beneath mixed layer
                ! [].
   integer :: &
      bdmtyp    ! Type of background diapycnal mixing. If bdmtyp = 1 the
                ! background diffusivity is a constant divided by the
                ! Brunt-Vaisala frequency, if bdmtyp = 2 the background
                ! diffusivity is constant [].
   logical :: &
      eddf2d, & ! If true, eddy diffusivity has a 2d structure.
      edsprs, & ! If true, apply eddy mixing suppression away from steering
                ! level.
      edanis, & ! If true, apply anisotropy correction to diffusivity.
      redi3d, & ! If true, then isopycnal/neutral diffusion will have 3D
                ! structure based in the 3D structure of anisotropy.
      rhsctp, & ! If true, use the minimum of planetary and topographic beta
                ! to define the Rhines scale.
      bdmldp, & ! If true, make the background mixing latitude dependent
                ! according to Gregg et al. (2003).
      smobld    ! If true, apply lateral smoothing of CVMix estimated boundary
                ! layer depth.
   character(len = fnmlen) :: &
      tbfile    ! Name of file containing topographic beta parameter.
   character(len = 80) :: &
      lngmtp, & ! Type of Langmuir turbulence parameterization. Valid types:
                ! 'none', 'vr12-ma', 'lf17'
      eitmth, & ! Eddy-induced transport parameterization method. Valid
                ! methods: 'intdif', 'gm'.
      edritp, & ! Type of Richardson number used in eddy diffusivity
                ! computation. Valid types: 'shear', 'large scale'.
      edwmth, & ! Method to estimate eddy diffusivity weight as a function of
                ! the ration of Rossby radius of deformation to the horizontal
                ! grid spacing. Valid methods: 'smooth', 'step'.
      ltedtp    ! Type of lateral tracer eddy diffusion: Valid methods: 'layer',
                ! 'neutral'.

   ! Options derived from string options.
   integer :: &
      eitmth_opt, &
      edritp_opt, &
      edwmth_opt, &
      ltedtp_opt

   ! Parameters:
   integer, parameter :: &
      ! Eddy-induced transport parameterization methods:
      eitmth_intdif      = 1, & ! Interface diffusion.
      eitmth_gm          = 2, & ! Gent-McWilliams.
      ! Type of Richardson number used in eddy diffusivity computation:
      edritp_shear       = 1, & ! Using local vertical velocity shear.
      edritp_large_scale = 2, & ! Using large scale variables.
      ! Method to estimate eddy diffusivity weight:
      edwmth_smooth      = 1, & ! Smooth function of Rossby radius over grid
                                ! spacing.
      edwmth_step        = 2, & ! Step function of Rossby radius over grid
                                ! spacing.
      ! Lateral tracer eddy diffusion type:
      ltedtp_layer       = 1, & ! Diffusion along model layers.
      ltedtp_neutral     = 2    ! Diffusion along neutral sublayers.

   real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy, kdm) :: &
      difint, & ! Layer interface diffusivity [m2 s-1].
      difiso, & ! Isopycnal diffusivity [m2 s-1].
      difdia    ! Diapycnal diffusivity [m2 s-1].

   real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy, kdm+1) :: &
      Kvisc_m, &      ! momentum eddy viscosity [m2 s-1].
      Kdiff_t, &      ! temperature eddy diffusivity [m2 s-1].
      Kdiff_s, &      ! salinity eddy  diffusivity [m2 s-1].
      t_ns_nonloc, &  ! Non-local transport term that is the fraction of
                      ! non-shortwave flux passing a layer interface [].
      s_nb_nonloc, &  ! Non-local transport term that is the fraction of
                      ! non-brine material tracer flux passing a layer interface
                      ! [].
      mu_nonloc, &    ! Non-local transport term that is the fraction of
                      ! u-component momentum flux passing a layer interface [].
      mv_nonloc       ! Non-local transport term that is the fraction of
                      ! u-component momentum flux passing a layer interface [].

   real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy) :: &
      difmxp, & ! Maximum lateral diffusivity at p-points [m2 s-1].
      difmxq, & ! Maximum lateral diffusivity at q-points [m2 s-1].
      difwgt    ! Eddy diffusivity weight [].

   real(r8), dimension(1 - nbdy:idm + nbdy,1 - nbdy:jdm + nbdy, 2*kdm) :: &
      umfltd, & ! u-component of horizontal mass flux due to thickness diffusion
                ! [kg m s-2].
      vmfltd, & ! v-component of horizontal mass flux due to thickness diffusion
                ! [kg m s-2].
      umflsm, & ! u-component of horizontal mass flux due to submesoscale
                ! eddy-induced transport [kg m s-2].
      vmflsm, & ! v-component of horizontal mass flux due to submesoscale
                ! eddy-induced transport [kg m s-2].
      utfltd, & ! u-component of horizontal heat flux due to thickness diffusion
                ! [K kg m s-2].
      vtfltd, & ! v-component of horizontal heat flux due to thickness diffusion
                ! [K kg m s-2].
      utflsm, & ! u-component of horizontal heat flux due to submesoscale
                ! eddy-induced transport [K kg m s-2].
      vtflsm, & ! v-component of horizontal heat flux due to submesoscale
                ! eddy-induced transport [K kg m s-2].
      utflld, & ! u-component of horizontal heat flux due to lateral diffusion
                ! [K kg m s-2].
      vtflld, & ! v-component of horizontal heat flux due to lateral diffusion
                ! [K kg m s-2].
      usfltd, & ! u-component of horizontal salt flux due to thickness diffusion
                ! [g m s-2].
      vsfltd, & ! v-component of horizontal salt flux due to thickness diffusion
                ! [g m s-2].
      usflsm, & ! u-component of horizontal salt flux due to submesoscale
                ! eddy-induced transport [g m s-2].
      vsflsm, & ! v-component of horizontal salt flux due to submesoscale
                ! eddy-induced transport [g m s-2].
      usflld, & ! u-component of horizontal salt flux due to lateral diffusion
                ! [g m s-2].
      vsflld    ! v-component of horizontal salt flux due to lateral diffusion
                ! [g m s-2].

   ! Public variables
   public :: egc, eggam, eglsmn, egmndf, egmxdf, egidfq, &
             rhiscf, ri0, bdmc1, bdmc2, bdmldp, tkepf, bdmtyp, &
             eddf2d, edsprs, edanis, redi3d, rhsctp, tbfile, smobld, lngmtp, &
             eitmth_opt, eitmth_intdif, eitmth_gm, edritp_opt, edritp_shear, &
             edritp_large_scale, edwmth_opt, edwmth_smooth, edwmth_step, &
             ltedtp_opt, ltedtp_layer, ltedtp_neutral, &
             difint, difiso, difdia, difmxp, difmxq, difwgt, &
             umfltd, vmfltd, umflsm, vmflsm, utfltd, vtfltd, &
             utflsm, vtflsm, utflld, vtflld, usfltd, vsfltd, &
             usflsm, vsflsm, usflld, vsflld, &
             Kvisc_m, Kdiff_t, Kdiff_s, &
             t_ns_nonloc, s_nb_nonloc, mu_nonloc, mv_nonloc

   ! Public routines
   public :: readnml_diffusion, inivar_diffusion

contains

   subroutine readnml_diffusion
      ! ---------------------------------------------------------------------------
      ! Read variables in the namelist group 'diffusion' and resolve options.
      ! ---------------------------------------------------------------------------

      ! Local variables
      character(len = 80) :: nml_fname
      integer :: nfu, ios
      logical :: fexist

      namelist /diffusion/ &
         egc, eggam, eglsmn, egmndf, egmxdf, egidfq, rhiscf, ri0, &
         bdmc1, bdmc2, bdmldp, tkepf, bdmtyp, eddf2d, edsprs, edanis, redi3d, &
         rhsctp, tbfile, smobld, lngmtp, eitmth, edritp, edwmth, ltedtp

      ! Read variables in the namelist group 'diffusion'.
      if (mnproc == 1) then
         nml_fname = 'ocn_in'//trim(inst_suffix)
         inquire(file = nml_fname, exist = fexist)
         if (fexist) then
            open (newunit = nfu, file = nml_fname, status = 'old', &
                 action = 'read')
         else
            nml_fname = 'limits'//trim(inst_suffix)
            inquire(file = nml_fname, exist = fexist)
            if (fexist) then
               open (newunit = nfu, file = nml_fname, status = 'old', &
                     action = 'read')
            else
               write (lp,*) 'readnml_diffusion: could not find namelist file!'
               call xchalt('(readnml_diffusion)')
                      stop '(readnml_diffusion)'
            endif
         endif
         read (unit = nfu, nml = diffusion, iostat = ios)
         close (unit = nfu)
      endif
      call xcbcst(ios)
      if (ios /= 0) then
         if (mnproc == 1) &
            write (lp,*) 'readnml_diffusion: No diffusion variable '//  &
                         'group found in namelist. Using defaults.'
      else
        call xcbcst(egc)
        call xcbcst(eggam)
        call xcbcst(eglsmn)
        call xcbcst(egmndf)
        call xcbcst(egmxdf)
        call xcbcst(egidfq)
        call xcbcst(rhiscf)
        call xcbcst(ri0)
        call xcbcst(bdmc1)
        call xcbcst(bdmc2)
        call xcbcst(bdmldp)
        call xcbcst(tkepf)
        call xcbcst(bdmtyp)
        call xcbcst(eddf2d)
        call xcbcst(edsprs)
        call xcbcst(edanis)
        call xcbcst(redi3d)
        call xcbcst(rhsctp)
        call xcbcst(tbfile)
        call xcbcst(smobld)
        call xcbcst(lngmtp)
        call xcbcst(eitmth)
        call xcbcst(edritp)
        call xcbcst(edwmth)
        call xcbcst(ltedtp)
      endif
      if (mnproc == 1) then
         write (lp,*) 'readnml_diffusion: diffusion variables:'
         write (lp,*) '  egc    = ', egc
         write (lp,*) '  eggam  = ', eggam
         write (lp,*) '  eglsmn = ', eglsmn
         write (lp,*) '  egmndf = ', egmndf
         write (lp,*) '  egmxdf = ', egmxdf
         write (lp,*) '  egidfq = ', egidfq
         write (lp,*) '  rhiscf = ', rhiscf
         write (lp,*) '  ri0    = ', ri0
         write (lp,*) '  bdmc1  = ', bdmc1
         write (lp,*) '  bdmc2  = ', bdmc2
         write (lp,*) '  bdmldp = ', bdmldp
         write (lp,*) '  tkepf  = ', tkepf
         write (lp,*) '  bdmtyp = ', bdmtyp
         write (lp,*) '  eddf2d = ', eddf2d
         write (lp,*) '  edsprs = ', edsprs
         write (lp,*) '  edanis = ', edanis
         write (lp,*) '  redi3d = ', redi3d
         write (lp,*) '  rhsctp = ', rhsctp
         write (lp,*) '  smobld = ', smobld
         write (lp,*) '  tbfile = ', trim(tbfile)
         write (lp,*) '  lngmtp = ', trim(lngmtp)
         write (lp,*) '  eitmth = ', trim(eitmth)
         write (lp,*) '  edritp = ', trim(edritp)
         write (lp,*) '  edwmth = ', trim(edwmth)
         write (lp,*) '  ltedtp = ', trim(ltedtp)
      endif

      ! Resolve options.
      select case (trim(eitmth))
         case ('intdif')
            eitmth_opt = eitmth_intdif
         case ('gm')
            eitmth_opt = eitmth_gm
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') &
                  ' readnml_diffusion: eitmth = ', trim(eitmth), &
                  ' is unsupported!'
            call xcstop('(readnml_diffusion)')
                   stop '(readnml_diffusion)'
      end select

      select case (trim(edritp))
         case ('shear')
            edritp_opt = edritp_shear
         case ('large scale')
            edritp_opt = edritp_large_scale
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') &
                  ' readnml_diffusion: edritp = ', trim(edritp), &
                  ' is unsupported!'
            call xcstop('(readnml_diffusion)')
                   stop '(readnml_diffusion)'
      end select

      select case (trim(edwmth))
         case ('smooth')
            edwmth_opt = edwmth_smooth
         case ('step')
            edwmth_opt = edwmth_step
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') &
                  ' readnml_diffusion: edwmth = ', trim(edwmth), &
                  ' is unsupported!'
            call xcstop('(readnml_diffusion)')
                   stop '(readnml_diffusion)'
      end select

      select case (trim(ltedtp))
         case ('layer')
            ltedtp_opt = ltedtp_layer
         case ('neutral')
            ltedtp_opt = ltedtp_neutral
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') &
                  ' readnml_diffusion: ltedtp = ', trim(ltedtp), &
                  ' is unsupported!'
            call xcstop('(readnml_diffusion)')
                   stop '(readnml_diffusion)'
      end select

      select case (trim(lngmtp))
         case ('none')
         case ('vr12-ma')
            if (wavsrc_opt == wavsrc_none) then
               if (mnproc == 1) &
                  write (lp,'(3a)') &
                     ' readnml_diffusion: lngmtp = ', trim(lngmtp), &
                     ' requires wavsrc = param or wavsrc = extern!'
               call xcstop('(readnml_diffusion)')
                      stop '(readnml_diffusion)'
            endif
         case ('lf17')
            if (wavsrc_opt /= wavsrc_extern) then
               if (mnproc == 1) &
                  write (lp,'(3a)') &
                     ' readnml_diffusion: lngmtp = ', trim(lngmtp), &
                     ' requires wavsrc = extern!'
               call xcstop('(readnml_diffusion)')
                      stop '(readnml_diffusion)'
            endif
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') &
                  ' readnml_diffusion: lngmtp = ', trim(lngmtp), &
                  ' is unsupported!'
            call xcstop('(readnml_diffusion)')
                   stop '(readnml_diffusion)'
      end select

   end subroutine readnml_diffusion

   subroutine inivar_diffusion
   ! ---------------------------------------------------------------------------
   ! Initialize arrays.
   ! ---------------------------------------------------------------------------

      integer :: i, j, k, l

      !$omp parallel do private(i, k)
      do j = 1 - nbdy, jj + nbdy
         do i = 1 - nbdy, ii + nbdy
            difmxp(i, j) = spval
            difmxq(i, j) = spval
            difwgt(i, j) = spval
         enddo
         do k = 1, kk
            do i = 1 - nbdy, ii + nbdy
               difint(i, j, k) = spval
               difiso(i, j, k) = spval
               difdia(i, j, k) = spval
            enddo
         enddo
         do k = 1, 2*kk
            do i = 1 - nbdy, ii + nbdy
               umfltd(i, j, k) = spval
               vmfltd(i, j, k) = spval
               umflsm(i, j, k) = spval
               vmflsm(i, j, k) = spval
               utfltd(i, j, k) = spval
               vtfltd(i, j, k) = spval
               utflsm(i, j, k) = spval
               vtflsm(i, j, k) = spval
               utflld(i, j, k) = spval
               vtflld(i, j, k) = spval
               usfltd(i, j, k) = spval
               vsfltd(i, j, k) = spval
               usflsm(i, j, k) = spval
               vsflsm(i, j, k) = spval
               usflld(i, j, k) = spval
               vsflld(i, j, k) = spval
            enddo
         enddo
         do k = 1, kk+1
            do i = 1 - nbdy, ii + nbdy
               Kvisc_m(i, j, k) = epsilk
               Kdiff_t(i, j, k) = epsilk
               Kdiff_s(i, j, k) = epsilk
               t_ns_nonloc(i, j, k) = spval
               s_nb_nonloc(i, j, k) = spval
               mu_nonloc(i, j, k) = spval
               mv_nonloc(i, j, k) = spval
            enddo
         enddo
      enddo
      !$omp end parallel do

      ! Initialize isopycnal diffusivity.
      !$omp parallel do private(k, l, i)
      do j = 1, jj
         do k = 1, kk
            do l = 1, isp(j)
            do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
               difiso(i, j, k) = 0._r8
            enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
      call xctilr(difiso, 1, kk, nbdy, nbdy, halo_ps)

      ! Initialize diffusive fluxes at points located upstream and downstream
      ! (in i-direction) of p-points.
      !$omp parallel do private(k, l, i)
      do j = 1, jj
         do k = 1, 2*kk
            do l = 1, isp(j)
            do i = max(1, ifp(j, l)), min(ii, ilp(j, l) + 1)
               umfltd(i, j, k) = 0._r8
               umflsm(i, j, k) = 0._r8
               utfltd(i, j, k) = 0._r8
               utflsm(i, j, k) = 0._r8
               utflld(i, j, k) = 0._r8
               usfltd(i, j, k) = 0._r8
               usflsm(i, j, k) = 0._r8
               usflld(i, j, k) = 0._r8
            enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
      call xctilr(umfltd, 1, 2*kk, nbdy, nbdy, halo_us)
      call xctilr(umflsm, 1, 2*kk, nbdy, nbdy, halo_us)
      call xctilr(utflld, 1, 2*kk, nbdy, nbdy, halo_us)
      call xctilr(usflld, 1, 2*kk, nbdy, nbdy, halo_us)

      ! Initialize diffusive fluxes at points located upstream and downstream
      ! (in j-direction) of p-points.
      !$omp parallel do private(k, l, j)
      do i = 1, ii
         do k = 1, 2*kk
            do l = 1, jsp(i)
            do j = max(1, jfp(i, l)), min(jj, jlp(i, l) + 1)
               vmfltd(i, j, k) = 0._r8
               vmflsm(i, j, k) = 0._r8
               vtfltd(i, j, k) = 0._r8
               vtflsm(i, j, k) = 0._r8
               vtflld(i, j, k) = 0._r8
               vsfltd(i, j, k) = 0._r8
               vsflsm(i, j, k) = 0._r8
               vsflld(i, j, k) = 0._r8
            enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
      call xctilr(vmfltd, 1, 2*kk, nbdy, nbdy, halo_vs)
      call xctilr(vmflsm, 1, 2*kk, nbdy, nbdy, halo_vs)
      call xctilr(vtflld, 1, 2*kk, nbdy, nbdy, halo_vs)
      call xctilr(vsflld, 1, 2*kk, nbdy, nbdy, halo_vs)

      ! Initialize non-local transport.
      !$omp parallel do private(k, l, i)
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
            t_ns_nonloc(i, j, 1) = 1._r8
            s_nb_nonloc(i, j, 1) = 1._r8
         enddo
         enddo
         do l = 1, isu(j)
         do i = max(1, ifu(j, l)), min(ii, ilu(j, l))
            mu_nonloc(i, j, 1) = 1._r8
         enddo
         enddo
         do l = 1, isv(j)
         do i = max(1, ifv(j, l)), min(ii, ilv(j, l))
            mv_nonloc(i, j, 1) = 1._r8
         enddo
         enddo
         do k = 2, kk + 1
            do l = 1, isp(j)
            do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
               t_ns_nonloc(i, j, k) = 0._r8
               s_nb_nonloc(i, j, k) = 0._r8
            enddo
            enddo
            do l = 1, isu(j)
            do i = max(1, ifu(j, l)), min(ii, ilu(j, l))
               mu_nonloc(i, j, k) = 0._r8
            enddo
            enddo
            do l = 1, isv(j)
            do i = max(1, ifv(j, l)), min(ii, ilv(j, l))
               mv_nonloc(i, j, k) = 0._r8
            enddo
            enddo
         enddo
      enddo
      !$omp end parallel do

   end subroutine inivar_diffusion

end module mod_diffusion
