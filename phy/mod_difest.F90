! ------------------------------------------------------------------------------
! Copyright (C) 2009-2025 Mats Bentsen, Mehmet Ilicak, Aleksi Nummelin,
!                         Mariana Vertenstein
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
! ----------------------------------------------------------------------------

module mod_difest

  use mod_types,             only: r8
  use mod_constants,         only: grav, alpha0, pi, epsilp, spval, onem, &
                                   onecm
  use mod_time,              only: delt1, dlt
  use mod_xc
  use mod_vcoord,            only: sigmar
  use mod_grid,              only: scpx, scpy, scp2, scuyi, scvxi, plat, &
                                   coriop, betafp, betatp, cosang, sinang, hangle
  use mod_eos,               only: rho
  use mod_state,             only: u, v, dp, dpu, dpv, temp, saln, sigma, p, &
                                   pbu, pbv, ubflxs_p, vbflxs_p, kfpla
  use mod_diffusion,         only: egc, eggam, eglsmn, egmndf, egmxdf, &
                                   egidfq, rhiscf, ri0, bdmc1, bdmc2, &
                                   bdmldp, tkepf, bdmtyp, eddf2d, edsprs, &
                                   edanis, redi3d, rhsctp, smobld, lngmtp, &
                                   edritp_opt, edritp_shear, &
                                   edritp_large_scale, edwmth_opt, &
                                   edwmth_smooth, edwmth_step, &
                                   difint, difiso, difdia, difmxp, difwgt, &
                                   Kvisc_m, Kdiff_t, Kdiff_s, &
                                   t_ns_nonloc, s_nb_nonloc, &
                                   mu_nonloc, mv_nonloc
  use mod_cmnfld,            only: bfsqi, nnslpx, nnslpy, mlts
  use mod_forcing,           only: wavsrc_opt, wavsrc_param, &
                                   abswnd, lamult, lasl, &
                                   ustar, ustarb, ustar3, wstar3, &
                                   buoyfl, t_sw_nonloc, surflx, sswflx, salflx
  use mod_tidaldissip,       only: twedon
  use mod_niw,               only: niwgf, niwbf, niwlf, idkedt, niw_ke_tendency
  use mod_seaice,            only: ficem
  use mod_utility,           only: util1
  use mod_checksum,          only: csdiag, chksummsk
  use CVMix_kpp,             only: CVMix_coeffs_kpp
  use CVMix_kpp,             only: CVMix_kpp_compute_turbulent_scales
  use CVMix_kpp,             only: CVMix_kpp_compute_bulk_Richardson
  use CVMix_kpp,             only: CVMix_kpp_compute_OBL_depth
  use CVMix_kpp,             only: CVmix_kpp_compute_unresolved_shear
  use CVMix_kpp,             only: CVMix_kpp_compute_kOBL_depth
  use CVMix_kpp,             only: cvmix_kpp_EFactor_model
  use CVMix_shear,           only: CVMix_init_shear, CVMix_coeffs_shear
  use CVMix_background,      only: CVMix_init_bkgnd, CVMix_coeffs_bkgnd
  use CVMix_convection,      only: CVMix_init_conv, CVMix_coeffs_conv
  use CVMix_tidal,           only: CVMix_init_tidal, CVMix_coeffs_tidal
  use CVMix_tidal,           only: CVMix_compute_Simmons_invariant, &
                                   CVMix_tidal_params_type
  use CVMix_kinds_and_types, only: CVMix_global_params_type
  use CVMix_kpp,             only: CVMix_kpp_params_type
  use CVMix_kpp,             only: CVMix_put_kpp
  use CVMix_kpp,             only: CVMix_init_kpp
  use CVMix_put_get,         only: CVMix_put
  use mod_tracers,           only: itrtke, itrgls, trc
  use mod_tke,               only: gls_cmu0, Pr_t, tke_min, gls_psi_min, gls_p, &
                                   gls_m, gls_n, gls_c1, gls_c2, gls_c3plus, &
                                   gls_c3minus, gls_Gh0, gls_Ghmin, gls_Ghcri, &
                                   Ls_unlmt_min, Prod, Buoy, Shear2, L_scale, &
                                   gls_s0, gls_s1, gls_s2, gls_s4, gls_s5, gls_s6, &
                                   gls_b0, gls_b1, gls_b2, gls_b3, gls_b4, gls_b5, &
                                   sqrt2, cmu_fac1, cmu_fac2, cmu_fac3, tke_exp1, &
                                   gls_exp1, gls_fac6
  use mod_ifdefs,            only: use_TRC, use_TKE, use_GLS

  implicit none
  private

  !     Initialize hOBL with hOBL_static = 3. for consistency with bulk
  !     mixed layer formulation in iHAMOCC: kmle = nint(hOBL) - 1 = 2
  real, PARAMETER :: hOBL_static = 3.

  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) :: &
       rig
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: &
       du2l,drhol,up,vp
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
       OBLdepth
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
       hOBL
  integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: &
       mskv,msku
  integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
       kmax,kfil

  type(CVMix_tidal_params_type)   :: CVMix_tidal_params
  type(CVMix_global_params_type)  :: CVMix_glb_params   !<CVMix-specific for Prandtl number only
  type(CVMix_kpp_params_type)     :: KPP_params
  !      type(CVMix_kpp_params_type), pointer :: CVmix_kpp_params_in
  !      type(CVMix_kpp_params_type) ::  CVmix_kpp_params_in

  ! parameters:
  !   iidtyp - type of interface and isopycnal diffusivities. If
  !            iidtyp=1 the diffusivities are diffusive velocities
  !            multiplied by the local horizontal grid scale, if
  !            iidtyp=2 the diffusivities are parameterized according
  !            to Eden and Greatbatch (2008).
  !   tdmflg - If tdmflg=1, apply tidally driven diapycnal mixing.
  !   iwdflg - If iwdflg=1, reduce background diapycnal diffusivity
  !            due to internal wave damping under sea-ice.
  !   dpbmin - smallest layer thickness allowed in evaluating
  !            local gradient richardson number [kg/m/s^2].
  !   drhomn - minimum density difference in evaluations the
  !            Brunt-Vaisala frequency and the local gradient
  !            Richardson number [kg/m*3].
  !   thkdff - diffusive velocity for thickness diffusion [m/s].
  !   temdff - diffusive velocity for tracer isopycnal diffusion [m/s].
  !   nu0    - diapycnal diffusivity when range of isopycnic physical
  !            layers is restricted [m^2/s].
  !   nus0   - maximum shear driven diapycnal diffusivity [m^2/s].
  !   nug0   - maximum gravity current diapycnal diffusivity [m^2/s].
  !   drho0  - critical local interface density difference [kg/m^3]
  !   nuls0  - maximum diapycnal diffusivity applied when local
  !            stability is weak [m^2/s].
  !   iwdfac - internal wave dissipation factor under sea ice [].
  !   dmxeff - diapycnal mixing efficiency [].
  !   tdmq   - tidal dissipation efficiency [].
  !   tdmls0 - tidal driven mixing length scale below critical
  !            latitude [kg/m/s^2].
  !   tdmls1 - tidal driven mixing length scale above critical
  !            latitude [kg/m/s^2].
  !   tdclat - critical latitude for tide M2 propagation [].
  !   tddlat - latitudinal transition zone for different tidal driven
  !            mixing length scales near the critical latitude.
  !   tkepls - length scale of surface TKE penetration beneath the
  !            mixed layer [kg/m/s^2]
  !   niwls  - near-inertial waves driven mixing length scale
  !            beneath the mixed layer [kg/m/s^2].
  !   cori30 - coriolis parameter at 30N [1/s].
  !   bvf0   - reference stratification in the parameterization of
  !            latitude dependent background diapycnal mixing [1/s].
  !   nubmin - minimum background diapycnal diffusivity [m^2/s].
  !   dpgc   - thickness of region near the bottom where the maximum
  !            diffusivity is increased due to gravity current mixing
  !            processes [kg/m/s^2].
  !   dpgrav - thickness of region below the non-isopycnic surface
  !            layers used to estimate upper ocean Eady growth rate
  !            [kg/m/s^2].
  !   dpdiav - thickness of region below the non-isopycnic surface
  !            layers used to estimate lateral diffusivities in the
  !            non-isopycnic layers [kg/m/s^2].
  !   dpddav - thickness of region below the non-isopycnic surface
  !            layers used to estimate diapycnal diffusivities in the
  !            non-isopycnic layers [kg/m/s^2].
  !   dpnbav - thickness of region near the bottom used to estimate
  !            bottom Brunt-Vaisala frequency [kg/m/s^2].
  !  cpsemin - Lower bound of zonal eddy phase speed minus zonal
  !            barotropic velocity [m/s].
  !  urmsemin- Lower bound of absolute value of RMS eddy velocity
  !            [m/s].
  integer , parameter :: iidtyp=2
  integer , parameter :: tdmflg=1
  integer , parameter :: iwdflg=1
  integer , parameter :: dptmin = onem
  real    , parameter ::  dpbmin=onecm
  real    , parameter :: drhomn = 6.e-3
  real    , parameter :: thkdff=5.e-3
  real    , parameter :: temdff = 3.5e-3
  real    , parameter :: nu0 = 1.e-5
  real    , parameter :: nus0 = 5.e-3
  real    , parameter :: nug0 = 2.5e-1
  real    , parameter :: drho0 = 6.e-3
  real    , parameter :: nuls0=5.e-2
  real    , parameter :: iwdfac = .06
  real    , parameter :: dmxeff=.2
  real    , parameter :: tdmq=1./3.
  real    , parameter :: tdmls0 = 500.*onem
  real    , parameter :: tdmls1=100.*onem
  real    , parameter :: tdclat=74.5
  real    , parameter :: tddlat = 3.
  real    , parameter :: tkepls=20.*onem
  real    , parameter :: niwls=300.*onem
  real    , parameter :: cori30 = 7.2722e-5
  real    , parameter :: bvf0=5.24e-3
  real    , parameter :: nubmin = 1.e-6
  real    , parameter :: dpgc=300.*onem
  real    , parameter :: dpgrav=100.*onem
  real    , parameter :: dpdiav = 100.*onem
  real    , parameter :: dpddav=10.*onem
  real    , parameter :: dpnbav=250.*onem
  real    , parameter :: ustmin = .001
  real    , parameter :: kappa=.4
  real    , parameter :: bfeps=1.e-16
  real    , parameter :: sleps=.1
  real    , parameter :: zetas = -1.
  real    , parameter :: cpsemin=-0.2
  real    , parameter :: urmsemin = 0.05
  real    , parameter :: as=-28.86
  real    , parameter :: cs=98.96
  real    , parameter :: minOBLdepth = 1.0

  character(len=16) :: langmuir_mixing_opt, langmuir_entrainment_opt

  public :: OBLdepth, inivar_difest, init_difest, difest_isobml, &
            difest_lateral_hybrid, difest_vertical_hybrid, hOBL

contains

  subroutine inivar_difest()

    !-----------------------------------------------------------
    ! Initialize arrays.
    !-----------------------------------------------------------

    ! Local variables
    integer :: i,j,k

    !$omp parallel do private(i,k)
    do j = 1-nbdy,jj+nbdy
      do k = 1,kk+1
        do i = 1-nbdy,ii+nbdy
          rig(i,j,k) = spval
        end do
      end do
      do k = 1,kk
        do i = 1-nbdy,ii+nbdy
          du2l(i,j,k) = spval
          drhol(i,j,k) = spval
          up(i,j,k) = spval
          vp(i,j,k) = spval
        end do
      end do
      do i = 1-nbdy,ii+nbdy
        OBLdepth(i,j) = spval
        hOBL(i,j) = hOBL_static
      end do
    end do
    !$omp end parallel do

  end subroutine inivar_difest

  subroutine init_difest

    !-----------------------------------------------------------
    ! Initialize CVmix variables.
    !-----------------------------------------------------------

    ! Local variables
    integer :: i,j,l

    select case (lngmtp)
      case ('none')
        if (mnproc == 1) then
          write (lp,'(a38)') 'no Langmuir turbulence parameterization'
        end if
        langmuir_mixing_opt      = 'NONE'
        langmuir_entrainment_opt = 'NONE'
      case ('vr12-ma')
        if (mnproc == 1) then
          write (lp,'(a38)') 'Langmuir param. type: vr12-ma    '
        end if
        langmuir_mixing_opt      = 'LWF16'
        langmuir_entrainment_opt = 'LWF16'
      case ('lf17')
        if (mnproc == 1) then
          write (lp,'(a38)') 'Langmuir param. type: lf17       '
        end if
        langmuir_mixing_opt      = 'LWF16'
        langmuir_entrainment_opt = 'LF17'
      case default
        if (mnproc == 1) then
          write (lp,'(3a)') "Error: '", trim(lngmtp), &
               "' is not a valid type of Langmuir turbulence parameterization"
        end if
        call xcstop('(init_difest)')
        stop '(init_difest)'
    end select

    ! -- ------- Background diapycnal mixing.
    !         The Bryan-Lewis parameterization is based on the following:
    !         \begin{eqnarray*}
    !         \kappa_{BL} &=& \textrm{bl1} +
    !                   \frac{\textrm{bl2}}{\pi}\tan^{-1}\bigg(
    !                        \textrm{bl3}(|z|-\textrm{bl4})\bigg)\\
    !         \nu_{BL} &=& \textrm{Pr}\cdot\kappa_{BL}
    !         \end{eqnarray*}
    !- Diapycnal mixing when local stability is weak
    !- convection routine based on N2 not rho
    !- if lBruntVaisala is TRUE, otherwise based on rho
    !- convert nuls0 to m2/s
    call CVMix_init_conv(convect_diff = 20.0*nuls0, &
         convect_visc = 20.0*nuls0, &
         lBruntVaisala = .true., &
         BVsqr_convect = 0.0)
    call CVMix_put(CVMix_glb_params,'max_nlev',kk)
    call CVMix_put(CVMix_glb_params,'Prandtl',1.0)
    call CVMix_put(CVMix_glb_params,'FreshWaterDensity',1000.0)
    call CVMix_put(CVMix_glb_params,'SaltWaterDensity',1025.0)
    call CVMix_put(CVMix_glb_params,'Gravity',grav)
    call cvmix_init_shear(mix_scheme = 'KPP', &
         KPP_nu_zero = nus0, &
         KPP_Ri_zero = ri0, &
         KPP_exp = 3.0)
    !  CVmix_kpp_params_in => CVmix_kpp_params_user
    !       call CVMix_init_kpp(Ri_crit=0.3, &
    !             minOBLdepth=minOBLdepth, &
    !             minVtsqr=1e-10, &
    !             vonKarman=0.4, &
    !             surf_layer_ext=0.1, &
    !             interp_type='quadratic', &
    !             interp_type2='LMD94', &
    !             lEkman=.false., &
    !             lMonOb=.false., &
    !             MatchTechnique='SimpleShapes', &
    !             lenhanced_diff=.true., &
    !             lnonzero_surf_nonlocal=.false.  , &
    !             lnoDGat1=.true.                 , &
    !             Langmuir_mixing_str=Langmuir_mixing_str, &
    !             Langmuir_entrainment_str=Langmuir_entrainment_str, &
    !             CVMix_kpp_params_user=KPP_params )
    !       call CVMix_init_kpp(Ri_crit=0.3, &
    !             minOBLdepth=minOBLdepth, &
    !             minVtsqr=1e-10, &
    !             vonKarman=0.4, &
    !             surf_layer_ext=0.1, &
    !             interp_type='quadratic', &
    !             interp_type2='LMD94', &
    !             lEkman=.false., &
    !             lMonOb=.false., &
    !             MatchTechnique='MatchGradient', &
    !             lenhanced_diff=.true., &
    !             lnonzero_surf_nonlocal=.false.  , &
    !             lnoDGat1=.false.                 , &
    !             Langmuir_mixing_str=Langmuir_mixing_str, &
    !             Langmuir_entrainment_str=Langmuir_entrainment_str, &
    !             CVMix_kpp_params_user=KPP_params )
    call CVMix_init_kpp(Ri_crit = 0.3, &
         minOBLdepth = minOBLdepth, &
         minVtsqr = 1e-10, &
         vonKarman = 0.4, &
         surf_layer_ext = 0.1, &
         interp_type = 'quadratic', &
         interp_type2 = 'LMD94', &
         lEkman = .false., &
         lMonOb = .false., &
         MatchTechnique = 'ParabolicNonLocal', &
         lenhanced_diff = .true., &
         lnonzero_surf_nonlocal = .true.  , &
         lnoDGat1 = .true.                 , &
         Langmuir_mixing_str = langmuir_mixing_opt, &
         Langmuir_entrainment_str = langmuir_entrainment_opt, &
         CVMix_kpp_params_user = KPP_params )

    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          OBLdepth(i,j) = 10.
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine init_difest

  subroutine difest_common_iso(m,n,mm,nn,k1m,k1n)

    !-----------------------------------------------------------
    ! Obtain common fields for the estimation of lateral and vertical
    ! diffusivities diapycnal diffusivities when vcoord == 'isopyc_bulkml'.
    !-----------------------------------------------------------

    integer :: m,n,mm,nn,k1m,k1n

    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: dv2
    real, dimension(1-nbdy:idm+nbdy,kdm) :: du2
    real, dimension(1-nbdy:idm+nbdy) :: tup
    integer, dimension(1-nbdy:idm+nbdy) :: kfpl,klpl
    integer :: i,j,k,l,kn
    real :: q

    ! Locate the range of layers to be considered in the computation of
    ! diffusivities.
    do j = 0,jj+1
      do i = 0,ii+1
        kmax(i,j) = 0
      end do
      do l = 1,isp(j)
        do i = max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          kmax(i,j) = 1
          do k = 3,kk
            kn = k+nn
            if (dp(i,j,kn) > dpbmin) kmax(i,j) = k
          end do
          if (kfpla(i,j,n) >= kmax(i,j)) then
            kfil(i,j) = kfpla(i,j,n)+1
          else
            if (sigma(i,j,kfpla(i,j,n)+nn) < &
            .5*(sigmar(i,j,kfpla(i,j,n)) + &
                sigmar(i,j,kfpla(i,j,n)+1))) then
              kfil(i,j) = kfpla(i,j,n)+1
            else
              kfil(i,j) = kfpla(i,j,n)+2
            end if
          end if
        end do
      end do
    end do

    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j) = kfil(i,j)
        end do
      end do
    end do
    !$omp end parallel do
    call xctilr(util1, 1,1, 1,1, halo_ps)
    !$omp parallel do private(l,i)
    do j = 0,jj+1
      do l = 1,isp(j)
        do i = max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          kfil(i,j) = nint(util1(i,j))
        end do
      end do
    end do
    !$omp end parallel do

    ! Compute squared vertical velocity difference of v-component
    !$omp parallel do private(l,i,kfpl,klpl,k,kn,q,tup)
    do j = 1,jj+1
      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          kfpl(i) = kk+1
          klpl(i) = 1
        end do
      end do
      do k = 3,kk
        kn = k+nn
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            if (dpv(i,j,kn) > dpbmin) then
              klpl(i) = k
            end if
          end do
        end do
      end do
      do k = kk,4,-1
        kn = k+nn
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            if (k >= max(kfil(i,j-1),kfil(i,j)).and.dpv(i,j,kn) > dptmin) then
              kfpl(i) = k
            end if
          end do
        end do
      end do
      do k = 1,kk
        kn = k+nn
        do i = 1,ii
          dv2(i,j,k) = 0.
          mskv(i,j,k) = 0
        end do
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            if (k >= kfpl(i).and.k <= klpl(i).and.klpl(i)-kfpl(i) >= 1) then
              if     (k == kfpl(i)) then
                q = v(i,j,kn+1)-v(i,j,kn)
                q = q*q
                dv2(i,j,k) = q
                tup(i) = q
              else if (k < klpl(i)) then
                q = v(i,j,kn+1)-v(i,j,kn)
                q = q*q
                dv2(i,j,k) = .5*(tup(i)+q)
                tup(i) = q
              else
                dv2(i,j,k) = tup(i)
              end if
              mskv(i,j,k) = 1
            end if
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(l,i,kfpl,klpl,k,kn,du2,q,tup)
    do j = 1,jj

      ! Compute squared vertical velocity difference of u-component
      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
          kfpl(i) = kk+1
          klpl(i) = 1
        end do
      end do
      do k = 3,kk
        kn = k+nn
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            if (dpu(i,j,kn) > dpbmin) klpl(i) = k
          end do
        end do
      end do
      do k = kk,4,-1
        kn = k+nn
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            if (k >= min(kfil(i-1,j),kfil(i,j)).and.dpu(i,j,kn) > dptmin) then
              kfpl(i) = k
            end if
          end do
        end do
      end do
      do k = 1,kk
        kn = k+nn
        do i = 1,ii+1
          du2(i,k) = 0.
          msku(i,j,k) = 0
        end do
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            if (k >= kfpl(i).and.k <= klpl(i).and.klpl(i)-kfpl(i) >= 1) then
              if     (k == kfpl(i)) then
                q = u(i,j,kn+1)-u(i,j,kn)
                q = q*q
                du2(i,k) = q
                tup(i) = q
              else if (k < klpl(i)) then
                q = u(i,j,kn+1)-u(i,j,kn)
                q = q*q
                du2(i,k) = .5*(tup(i)+q)
                tup(i) = q
              else
                du2(i,k) = tup(i)
              end if
              msku(i,j,k) = 1
            end if
          end do
        end do
      end do

      ! - Centered at layers, compute vertical in-situ density difference,
      ! - squared vertical velocity difference and local gradient
      ! - Richardson number.
      do k = 4,kk
        kn = k+nn
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k >= kfil(i,j).and.k <= kmax(i,j).and. &
                 kmax(i,j)-kfil(i,j) >= 1) then

              ! Vertical in-situ density difference.
              if     (k == kfil(i,j)) then
                q = max(0.,rho(p(i,j,k+1),temp(i,j,kn+1),saln(i,j,kn+1)) &
                          -rho(p(i,j,k+1),temp(i,j,kn  ),saln(i,j,kn  )))
                drhol(i,j,k) = q
                tup(i) = q
              else if (k < kmax(i,j)) then
                q = max(0.,rho(p(i,j,k+1),temp(i,j,kn+1),saln(i,j,kn+1)) &
                          -rho(p(i,j,k+1),temp(i,j,kn  ),saln(i,j,kn  )))
                drhol(i,j,k) = 2.*tup(i)*q/max(1.e-11,tup(i)+q)
                tup(i) = q
              else
                drhol(i,j,k) = tup(i)
              end if

              ! Vertical squared velocity difference.
              du2l(i,j,k) = (msku(i  ,j,k)*du2(i  ,k) &
                            +msku(i+1,j,k)*du2(i+1,k)) &
                            /max(1,msku(i,j,k)+msku(i+1,j,k)) &
                            +(mskv(i,j  ,k)*dv2(i,j  ,k) &
                             +mskv(i,j+1,k)*dv2(i,j+1,k)) &
                            /max(1,mskv(i,j,k)+mskv(i,j+1,k))

              ! Local gradient Richardson number.
              rig(i,j,k) = alpha0*alpha0 &
                   *max(drhomn,drhol(i,j,k))*dp(i,j,kn) &
                   /max(1.e-13,du2l(i,j,k))

            end if
          end do
        end do
      end do

    end do

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'difest_common_iso:'
      end if
      call chksummsk(drhol,ip,kk,'drhol')
      call chksummsk(du2l,ip,kk,'du2l')
      call chksummsk(rig,ip,kk,'rig')
    end if

  end subroutine difest_common_iso

  subroutine difest_common_hyb(m,n,mm,nn,k1m,k1n)

    !-----------------------------------------------------------
    ! Obtain common fields for the estimation of lateral and vertical
    ! diffusivities diapycnal diffusivities when vcoord == 'isopyc_bulkml'.
    !-----------------------------------------------------------

    integer :: m,n,mm,nn,k1m,k1n

    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: dv2
    real, dimension(1-nbdy:idm+nbdy,kdm) :: du2
    integer, dimension(1-nbdy:idm+nbdy) :: klpl
    integer :: i,j,k,l,kn
    real :: q,dz

    ! Compute squared vertical velocity difference of v-component
    !$omp parallel do private(l,i,klpl,k,kn,q)
    do j = 1,jj+1
      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          klpl(i) = 1
        end do
      end do
      do k = 2,kk
        kn = k+nn
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            if (dpv(i,j,kn) > dpbmin) klpl(i) = k
          end do
        end do
      end do
      do i = 1,ii
        mskv(i,j,1) = iv(i,j)
      end do
      do k = 2,kk
        do i = 1,ii
          dv2(i,j,k) = 0.
          mskv(i,j,k) = 0
        end do
        kn = k+nn
        do l = 1,isv(j)
          do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
            if (klpl(i) >= 2.and.k <= klpl(i)) then
              q = v(i,j,kn)-v(i,j,kn-1)
              dv2(i,j,k) = q*q
              mskv(i,j,k) = 1
            end if
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(l,i,klpl,k,kn,du2,q,dz)
    do j = 1,jj

      ! Compute squared vertical velocity difference of u-component
      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
          klpl(i) = 1
        end do
      end do
      do k = 2,kk
        kn = k+nn
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            if (dpu(i,j,kn) > dpbmin) klpl(i) = k
          end do
        end do
      end do
      do i = 1,ii+1
        msku(i,j,1) = iu(i,j)
      end do
      do k = 2,kk
        do i = 1,ii+1
          du2(i,k) = 0.
          msku(i,j,k) = 0
        end do
        kn = k+nn
        do l = 1,isu(j)
          do i = max(1,ifu(j,l)),min(ii+1,ilu(j,l))
            if (klpl(i) >= 2.and.k <= klpl(i)) then
              q = u(i,j,kn)-u(i,j,kn-1)
              du2(i,k) = q*q
              msku(i,j,k) = 1
            end if
          end do
        end do
      end do

      ! - Compute local gradient Richardson number at interfaces.
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          rig(i,j,1) = 0.
        end do
      end do
      do k = 2,kk
        kn = k+nn
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (msku(i,j,k)+msku(i+1,j,k)+mskv(i,j,k)+mskv(i,j+1,k) > 0) then
              q = (msku(i,j,k)*du2(i,k)  +msku(i+1,j,k)*du2(i+1,k)) &
                   /max(1,msku(i,j,k)+msku(i+1,j,k)) &
                   +(mskv(i,j,k)*dv2(i,j,k)+mskv(i,j+1,k)*dv2(i,j+1,k)) &
                   /max(1,mskv(i,j,k)+mskv(i,j+1,k))
              dz = .5*(dp(i,j,kn-1)+dp(i,j,kn))*alpha0/grav
              rig(i,j,k) = max(0.,bfsqi(i,j,k)*dz*dz) &
                          /max(1.e-13,q)
            else
              rig(i,j,k) = rig(i,j,k-1)
            end if
          end do
        end do
      end do
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          rig(i,j,1) = rig(i,j,2)
          rig(i,j,kk+1) = rig(i,j,kk)
        end do
      end do

      ! - Compute velocity components at p-points.
      do k = 1,kk
        kn = k+nn
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            up(i,j,k) = (msku(i,j,k)*u(i,j,kn)+msku(i+1,j,k)*u(i+1,j,kn)) &
                        /max(1,msku(i,j,k)+msku(i+1,j,k))
            vp(i,j,k) = (mskv(i,j,k)*v(i,j,kn)+mskv(i,j+1,k)*v(i,j+1,kn)) &
                        /max(1,mskv(i,j,k)+mskv(i,j+1,k))
          end do
        end do
      end do

    end do

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'difest_common_hyb:'
      end if
      call chksummsk(rig,ip,kk+1,'rig')
      call chksummsk(up,ip,kk,'up')
      call chksummsk(vp,ip,kk,'vp')
    end if

  end subroutine difest_common_hyb

  subroutine difest_isobml(m,n,mm,nn,k1m,k1n)

    !-----------------------------------------------------------
    ! estimate diffusivities for eddy-induced transport, layer-wise
    ! diffusion and vertical diffusion
    !-----------------------------------------------------------

    integer :: m,n,mm,nn,k1m,k1n

    integer :: i,j,k,l,kn

    !-----------------------------------------------------------
    ! update halos of various fields
    !-----------------------------------------------------------

    call xctilr(u, 1,2*kk, 2,2, halo_uv)
    call xctilr(v, 1,2*kk, 2,2, halo_vv)
    call xctilr(ubflxs_p, 1,2, 2,2, halo_uv)
    call xctilr(vbflxs_p, 1,2, 2,2, halo_vv)
    call xctilr(pbu, 1,2, 2,2, halo_us)
    call xctilr(pbv, 1,2, 2,2, halo_vs)

    !-----------------------------------------------------------
    ! Update layer interface pressure.
    !-----------------------------------------------------------

    !$omp parallel do private(k,kn,l,i)
    do j = -2,jj+3
      do k = 1,kk
        kn = k+nn
        do l = 1,isp(j)
          do i = max(-2,ifp(j,l)),min(ii+3,ilp(j,l))
            p(i,j,k+1) = p(i,j,k)+dp(i,j,kn)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !-----------------------------------------------------------
    ! Estimate friction velocity cubed.
    !-----------------------------------------------------------

    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          ustar3(i,j) = ustar(i,j)**3
        end do
      end do
    end do
    !$omp end parallel do

    ! Estimate energy input by near-inertial waves.
    call niw_ke_tendency(m,n,mm,nn,k1m,k1n)

    ! Obtain common fields for the estimation of lateral and vertical
    ! diffusivities diapycnal diffusivities.
    call difest_common_iso(m,n,mm,nn,k1m,k1n)

    ! Estimate vertical diffusivity.
    call difest_vertical_iso(m,n,mm,nn,k1m,k1n)

    ! Estimate diffusivities for eddy-induced transport and layer-wise
    ! diffusion.
    call difest_lateral_iso(m,n,mm,nn,k1m,k1n)

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'difest_isobml:'
      end if
      call chksummsk(ustar3,ip,1,'ustar3')
    end if

  end subroutine difest_isobml

  subroutine difest_lateral_hybrid(m,n,mm,nn,k1m,k1n)

    !-----------------------------------------------------------
    ! estimate diffusivities for eddy-induced transport, layer-wise
    ! diffusion and vertical diffusion
    !-----------------------------------------------------------

    integer :: m,n,mm,nn,k1m,k1n

    integer :: i,j,k,l,kn

    !-----------------------------------------------------------
    ! update halos of various fields
    !-----------------------------------------------------------

    call xctilr(u, 1,2*kk, 2,2, halo_uv)
    call xctilr(v, 1,2*kk, 2,2, halo_vv)
    call xctilr(ubflxs_p, 1,2, 2,2, halo_uv)
    call xctilr(vbflxs_p, 1,2, 2,2, halo_vv)
    call xctilr(pbu, 1,2, 2,2, halo_us)
    call xctilr(pbv, 1,2, 2,2, halo_vs)

    !-----------------------------------------------------------
    ! Estimate friction velocity cubed.
    !-----------------------------------------------------------

    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          ustar3(i,j) = ustar(i,j)**3
        end do
      end do
    end do
    !$omp end parallel do

    ! Obtain common fields for the estimation of lateral and vertical
    ! diffusivities diapycnal diffusivities.
    call difest_common_hyb(m,n,mm,nn,k1m,k1n)

    ! Estimate diffusivities for eddy-induced transport and layer-wise
    ! diffusion.
    call difest_lateral_hyb(m,n,mm,nn,k1m,k1n)

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'difest_lateral_hybrid:'
      end if
      call chksummsk(ustar3,ip,1,'ustar3')
    end if

  end subroutine difest_lateral_hybrid

  subroutine difest_vertical_hybrid(m,n,mm,nn,k1m,k1n)

    ! -----------------------------------------------------------
    ! estimate diffusivities for eddy-induced transport, layer-wise
    ! diffusion and vertical diffusion
    ! -----------------------------------------------------------

    integer :: m,n,mm,nn,k1m,k1n

    ! -----------------------------------------------------------
    ! update halos of various fields
    ! -----------------------------------------------------------

    call xctilr(u(1-nbdy,1-nbdy,k1n), 1,kk, 1,1, halo_uv)
    call xctilr(v(1-nbdy,1-nbdy,k1n), 1,kk, 1,1, halo_vv)

    ! Obtain common fields for the estimation of lateral and vertical
    ! diffusivities diapycnal diffusivities.
    call difest_common_hyb(m,n,mm,nn,k1m,k1n)

    ! Estimate vertical diffusivities..
    call difest_vertical_hyb(m,n,mm,nn,k1m,k1n)

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'difest_vertical_hybrid:'
      end if
    end if

  end subroutine difest_vertical_hybrid

  subroutine difest_vertical_hyb(m,n,mm,nn,k1m,k1n)

    ! -----------------------------------------------------------
    ! estimate layer diapycnal, diffusivities for hybrid
    ! coordinates
    ! -----------------------------------------------------------

    integer :: m,n,mm,nn,k1m,k1n

    real, dimension(kdm+1) :: rig_i, rig_i_prev, bvfsq_i_prev
    integer :: i,j,k,l,kn
    real :: q

    type(CVMix_tidal_params_type)   :: CVMix_tidal_params
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) :: z_int
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: z_mid
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: rig_lf, bfsqi_lf
    real, dimension(kdm+1) :: depth_int
    real, dimension(kdm+1) :: Kv_col, Kd_col ! background visc/diff
    real, dimension(kdm+1) :: Kv_shr, Kd_shr ! shear driven visc/diff
    real, dimension(kdm+1) :: Kv_conv, Kd_conv ! convection visc/diff
    real, dimension(kdm+1) :: vert_dep       ! vertical deposition
    real, dimension(kdm+1) :: Kv_tidal, Kd_tidal ! tidal viscosity,diffusivity
    real, dimension(kdm+1) :: Kv_kpp, Kt_kpp, Ks_kpp ! vertical viscosity,diffusivity temp/salt
    real, dimension(kdm+1) :: iFaceHeight    ! Height of interfaces [m]
    real, dimension(kdm+1) :: bvfsq_i, bvf_i ! N2, N at interfaces
    real, dimension(kdm) :: cellHeight       ! Height of cell centers [m]
    real, dimension(kdm) :: rho_zeros, rho_lwr  ! dummy vars for convection
    real, dimension(kdm) :: rho_1d           ! 1D density at the layer center
    real, dimension(kdm) :: Ws_1d            ! Profile of vertical velocity scale for scalars [m s-1]
    real, dimension(kdm) :: surfBuoyFlux2
    real, dimension(kdm) :: BulkRi_1d        ! Bulk Richardson number for each layer
    real, dimension(kdm) :: deltaU2          ! square of delta U (shear) in denominator of Bulk Ri [m2 s-2]
    real, dimension(kdm) :: VT2              ! unresolved shear used for  Bulk Ri
    real, dimension(kdm) :: deltaRho         ! delta Rho [kg/m3] in numerator of Bulk Ri number
    real, dimension(kdm+1,2) :: nonLocalTrans  ! Non-local transport for heat/salt at interfaces [nondim]
    real :: surf_layer_ext, surfFricVel
    real :: surfBuoyFlux
    real :: delH, bvfbot, dps
    real :: dh, hcorr
    real :: Uk, Vk
    real :: surfU, surfV
    real :: surfHu, surfHv
    real :: surfTemp, surfSalt, surfRho
    real :: surfHtemp, surfHsalt
    real :: SLdepth_0d, hTot
    real :: Simmons_coeff, zBottomMinusOffset
    real :: bl1, bl2, bl3, bl4
    real :: ws, ww, we, wn, wc
    integer :: ki, ksfc, ktmp, kobl, kn1, n_iter

    surf_layer_ext = 0.1
    bl1 = 8e-5
    bl2 = 1.05e-4
    bl3 = 4.5e-3
    bl4 = 2500.0

    ! Obtain interface and cell heights.
    do j = 0,jj+1
      do l = 1,isp(j)
        do i = max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          z_int(i,j,1) = -p(i,j,1)/onem
          hcorr = 0.
          do k = 1,kk
            kn = k + nn
            dh = dp(i,j,kn)/onem
            dh = dh + hcorr ! Take away the accumulated error (could temporarily make dh<0)
            hcorr = min(dh - 1e-10, 0.) ! If inflating then hcorr<0
            dh = max(dh, 1e-10) ! Limit increment  dh>=min_thicknes
            z_mid(i,j,k) = z_int(i,j,k) - 0.5 * dh
            z_int(i,j,k+1) = z_int(i,j,k) - dh
          end do
        end do
      end do
    end do

    call xctilr(rig, 1,kk, 1,1, halo_ps)

    do k = 1,kk
      kn = k+nn
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            ws = .125*ip(i  ,j-1)*min(onem,dp(i  ,j-1,kn))/onem
            ww = .125*ip(i-1,j  )*min(onem,dp(i-1,j  ,kn))/onem
            we = .125*ip(i+1,j  )*min(onem,dp(i+1,j  ,kn))/onem
            wn = .125*ip(i  ,j+1)*min(onem,dp(i  ,j+1,kn))/onem
            wc = 1. - (ws + ww + we + wn)
            rig_lf(i,j,k) =   ws*rig(i  ,j-1,k) &
                            + ww*rig(i-1,j  ,k) &
                            + wc*rig(i  ,j  ,k) &
                            + we*rig(i+1,j  ,k) &
                            + wn*rig(i  ,j+1,k)
            bfsqi_lf(i,j,k) =   ws*bfsqi(i  ,j-1,k) &
                              + ww*bfsqi(i-1,j  ,k) &
                              + wc*bfsqi(i  ,j  ,k) &
                              + we*bfsqi(i+1,j  ,k) &
                              + wn*bfsqi(i  ,j+1,k)
          end do
        end do
      end do
    end do

    ! single column diffusivity
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))

          ! -- ------- CVMix variables computed below

          surfBuoyFlux2(:) = 0.0
          surfBuoyFlux = 0.0
          iFaceHeight(:) = 0.
          cellHeight(:) = 0.
          rho_1d(:)    = 0.
          rig_i(:) = 1.e8 !Initialize w/ large Richardson value

          iFaceHeight(1) = z_int(i,j,1)
          surfFricVel = ustar(i,j)
          surfBuoyFlux = buoyfl(i,j,2) - buoyfl(i,j,1)
          do k = 1,kk
            kn = k + nn
            kn1 = max(nn+1,kn-1)
            cellHeight(k) = z_mid(i,j,k)
            iFaceHeight(k+1) = z_int(i,j,k+1)

            ! compute rho_1d at the interfaces
            rho_1d(k) = rho(p(i,j,k),temp(i,j,kn),saln(i,j,kn))

            ! find ksfc for cell where "surface layer" sits
            SLdepth_0d = surf_layer_ext* &
                 max(max(-cellHeight(k),-iFaceHeight(2)), &
                 minOBLdepth)
            ksfc = k
            do ki = 1,k
              if (-1.0*iFaceHeight(ki+1) >= SLdepth_0d) then
                ksfc = ki
                exit
              end if
            end do
            surfHu    = 0.0
            surfHv    = 0.0
            surfHtemp = 0.0
            surfHsalt = 0.0
            hTot      = 0.0
            do ki = 1,ksfc
              ktmp = ki+nn
              ! SLdepth_0d can be between cell interfaces
              delH = min( max(0.0, SLdepth_0d - hTot), &
                   dp(i,j,ktmp)/onem )
              ! surface layer thickness
              hTot = hTot + delH
              ! surface averaged fields
              surfHtemp = surfHtemp + temp(i,j,ktmp)*delH
              surfHsalt = surfHsalt + saln(i,j,ktmp)*delH
              surfHu = surfHu+up(i,j,ki)*delH
              surfHv = surfHv+vp(i,j,ki)*delH
            end do
            surfTemp = surfHtemp / hTot
            surfSalt = surfHsalt / hTot
            surfU    = surfHu / hTot
            surfV    = surfHv / hTot
            surfRho  = rho(p(i,j,k),surfTemp,surfSalt)
            if (p(i,j,kk+1)-p(i,j,k) < epsilp) then
              deltaRho(k) = deltaRho(k-1)
            else
              deltaRho(k) = rho_1d(k) - surfRho
            end if
            ! vertical shear between present layer and
            ! surface layer averaged surfU,surfV.
            ! C-grid average to get Uk and Vk on T-points.
            Uk = up(i,j,k) - surfU
            Vk = vp(i,j,k) - surfV
            deltaU2(k) = (Uk**2 + Vk**2)

            ! XXX: Temporary de-scaling of N2_int(i,:) into a
            ! temporary variable
            bvfsq_i(k) = bfsqi_lf(i,j,k)

            ! Local gradient Richardson number
            rig_i(k) = rig_lf(i,j,k)

            surfBuoyFlux2(k) = ( buoyfl(i,j,k+1) &
                 - buoyfl(i,j,1  ))

          end do  ! k

          ! bottom values for the Ri and N
          rig_i(kk+1) = rig_i(kk)
          bvfsq_i(kk+1) = bvfsq_i(kk)

          rig_i_prev(:) = rig_i(:)
          bvfsq_i_prev(:) = bvfsq_i(:)
          do k = 2, kk
             rig_i(k) = .5_r8*rig_i_prev(k) + .25*(rig_i_prev(k-1) + rig_i_prev(k+1))
             bvfsq_i(k) = .5_r8*bvfsq_i_prev(k) + .25*(bvfsq_i_prev(k-1) + bvfsq_i_prev(k+1))
          enddo
          bvf_i(:) = sqrt( max( bvfsq_i(:), 0.) )

          !- turbulent velocity scales w_s and w_m computed at the cell
          !- centers.
          call CVMix_kpp_compute_turbulent_scales( &
               surf_layer_ext,   & ! (in)  Normalized surface layer Cdepth; sigma = CS%surf_layer_ext
               -cellHeight(:),   & ! (in)  Assume here that OBL depth [m] = -cellHeight(k)
               surfBuoyFlux2(:), & ! (in)  Buoyancy flux at surface [m2 s-3] &
               surfFricVel,      & ! (in)  Turbulent friction  velocity at surface [m s-1] &
               w_s = Ws_1d(:),   & ! (out) Turbulent velocity scale profile [m s-1] &
               CVMix_kpp_params_user = KPP_params)

          if (wavsrc_opt == wavsrc_param) then
            lamult(i,j) = cvmix_kpp_EFactor_model( &
                          abswnd(i,j), &
                          surfFricVel, &
                          OBLdepth(i,j), &
                          CVMix_glb_params)
            lamult(i,j) = lamult(i,j)*(1. - ficem(i,j)) + ficem(i,j)
          end if

          ! Compute unresolved shear for CVMix
          VT2(:) = CVmix_kpp_compute_unresolved_shear( &
               zt_cntr = cellHeight(:), & ! Depth ofcell center [m]
               ws_cntr = Ws_1d(:),      & ! Turbulent velocity scale profile, at centers [m s-1]
               N_iface = bvf_i(:),      & ! Buoyancy frequency at the interface [s-1]
               EFactor = lamult(i,j),   & ! Langmuir enhancement factor []
               LaSL = lasl(i,j),        & ! Surface layer averaged Langmuir number []
               bfsfc = surfBuoyFlux,    & ! Surface buoyancy flux [m2 s-3]
               ustar = surfFricVel,     & ! Friction velocity []
               CVMix_kpp_params_user = KPP_params)

          ! Calculate Bulk Richardson number from eq (21) of LMD94
          BulkRi_1d(:) = CVmix_kpp_compute_bulk_Richardson(       &
               zt_cntr = cellHeight(:),                           & ! Depth of cell  center [m]
               delta_buoy_cntr = grav*alpha0*deltaRho(:),         & ! Bulk buoyancy difference, Br-B(z) [m s-2]
               delta_Vsqr_cntr = deltaU2(:),                      & ! Square of resolved velocity difference [m2 s-2]
               Vt_sqr_cntr = VT2(:),                              & ! Unresolved shear [m2 s-2]
               ws_cntr = Ws_1d(:),                                & ! Turbulent velocity scale profile [m s-1]
               N_iface = bvf_i(:),                                & ! Buoyancy frequency at the interface [s-1]
               CVMix_kpp_params_user = KPP_params )                 ! KPP parameters

          ! Compute OBL depth for KPP
          call CVMix_kpp_compute_OBL_depth( &
               BulkRi_1d(:),                & ! (in) Bulk Richardson number
               iFaceHeight(:),              & ! (in)  Height of interfaces [m]
               OBLdepth(i,j),               & ! (out) OBL depth [m]
               hOBL(i,j),                   & ! (out) level (+fraction) of OBL extent
               zt_cntr = cellHeight(:),     & ! Depth of cell  center [m]
               surf_fric = surfFricVel,     & ! (in)  Turbulent friction  velocity at surface [m s-1]
               surf_buoy = surfBuoyFlux,    & ! (in) Buoyancy flux at surface [m2 s-3]
               Coriolis = coriop(i,j),      & ! (in) Coriolis parameter [s-1]
               CVMix_kpp_params_user = KPP_params ) ! KPP parameters

          ! Avoid KPP reaching bottom
          zBottomMinusOffset = iFaceHeight(kk+1) &
               + min(1.0,-0.1*iFaceHeight(kk+1))
          OBLdepth(i,j) = min(OBLdepth(i,j), -zBottomMinusOffset)
          ! no shallower than top layer
          OBLdepth(i,j) = max(OBLdepth(i,j), -iFaceHeight(2))
          ! no deeper than bottom
          OBLdepth(i,j) = min(OBLdepth(i,j), -iFaceHeight(kk+1))

        end do
      end do
    end do

    if (smobld) then

      ! -- Smooth OBLdepth
      call xctilr(OBLdepth, 1,1, 2,2, halo_ps)
      util1(:,:) = OBLdepth(:,:)
      do j = 0,jj+1
        do l = 1,isp(j)
          do i = max(0,ifp(j,l)),min(ii+1,ilp(j,l))
            ws = .125*ip(i  ,j-1)
            ww = .125*ip(i-1,j  )
            we = .125*ip(i+1,j  )
            wn = .125*ip(i  ,j+1)
            wc = 1. - (ws + ww + we + wn)
            OBLdepth(i,j) = ws*util1(i  ,j-1) &
                          + ww*util1(i-1,j  ) &
                          + wc*util1(i  ,j  ) &
                          + we*util1(i+1,j  ) &
                          + wn*util1(i  ,j+1)
            OBLdepth(i,j) = min(OBLdepth(i,j), -z_int(i,j,kk+1))
          end do
        end do
      end do
    else
      call xctilr(OBLdepth, 1,1, 1,1, halo_ps)
    end if

    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))

          ! -- ------- CVMix variables computed below

          surfBuoyFlux = 0.0
          bvfbot    = 0.
          dps       = 0.
          depth_int(:) = 0.
          Kv_col(:)    = 0.
          Kd_col(:)    = 0.
          vert_dep(:)  = 0.
          Kv_tidal(:)  = 0.
          Kd_tidal(:)  = 0.
          Kv_conv(:)   = 0.
          Kd_conv(:)   = 0.
          Kv_shr(:)    = 0.
          Kd_shr(:)    = 0.
          iFaceHeight(:) = 0.
          cellHeight(:) = 0.
          bvfsq_i(:)   = 0.
          rho_lwr(:)= drho0
          rho_zeros(:)= 0.
          nonLocalTrans(:,:) = 0.0
          rig_i(:) = 1.e8 !Initialize w/ large Richardson value
          Kv_kpp(:) = 0.0
          Kt_kpp(:) = 0.0
          Ks_kpp(:) = 0.0
          do k = 1,kk+1
            Kv_kpp(k) = Kvisc_m(i,j,k)
            Kt_kpp(k) = Kdiff_t(i,j,k)
            Ks_kpp(k) = Kdiff_s(i,j,k)
          end do
          iFaceHeight(1) = z_int(i,j,1)
          depth_int(1) = -iFaceHeight(1)
          surfFricVel = ustar(i,j)

          do k = 1,kk
            kn = k+nn
            cellHeight(k) = z_mid(i,j,k)
            iFaceHeight(k+1) = z_int(i,j,k+1)
            depth_int(k+1) = -iFaceHeight(k+1)

            ! XXX: Temporary de-scaling of N2_int(i,:) into a
            ! temporary variable
            bvfsq_i(k) = bfsqi_lf(i,j,k)
            bvf_i(k) = sqrt( max( bvfsq_i(k), 0.) )
            ! Accumulate Brunt-Vaisala frequency in a region near the
            ! bottom
            q = max(0.,p(i,j,k+1)-max(p(i,j,kk+1)-dpnbav,p(i,j,k)))
            if (q > 0.) then
              bvfbot = bvfbot+bvf_i(k)*q
              dps = dps+q
            end if

            ! Local gradient Richardson number
            rig_i(k) = rig_lf(i,j,k)

          end do  ! k

          if(dps > 0.) bvfbot = bvfbot/dps

          ! bottom values for the Ri and N2
          rig_i(kk+1) = rig_i(kk)
          bvfsq_i(kk+1) = bfsqi(i,j,kk+1)

          rig_i_prev(:) = rig_i(:)
          bvfsq_i_prev(:) = bvfsq_i(:)
          do k = 2, kk
             rig_i(k) = .5_r8*rig_i_prev(k) + .25*(rig_i_prev(k-1) + rig_i_prev(k+1))
             bvfsq_i(k) = .5_r8*bvfsq_i_prev(k) + .25*(bvfsq_i_prev(k-1) + bvfsq_i_prev(k+1))
          enddo

          ! gets index of the level and interface above hbl
          hOBL(i,j) = CVMix_kpp_compute_kOBL_depth(iFaceHeight(:), &
               cellHeight(:),OBLdepth(i,j))

          ! -- ------- Background diapycnal mixing.
          if (bdmtyp == 1) then
            ! zw interface depths relative to the surface in m, must be positive.
            call CVMix_init_bkgnd(max_nlev = kk, zw = depth_int(:), &
                 bl1 = bl1, bl2 = bl2, bl3 = bl3, bl4 = bl4, &
                 prandtl = CVMix_glb_params%Prandtl)
            call CVMix_coeffs_bkgnd(Mdiff_out = Kv_col(:), &
                 Tdiff_out=Kd_col(:), nlev=kk, max_nlev = kk)
          else if (bdmtyp == 2) then
            ! Type 2: Background diffusivity is a constant
            Kv_col(:) = bdmc2
            Kd_col(:) = bdmc2
          else
            Kv_col(:) = 0.
            Kd_col(:) = 0.
          end if
          if (iwdflg == 1) then
            Kv_col = Kv_col*(1.+(iwdfac-1.)*ficem(i,j))
            Kd_col = Kd_col*(1.+(iwdfac-1.)*ficem(i,j))
          end if

          !- Latitude dependency of background diapycnal mixing
          if (bdmldp) then
            q = max(1.e-9,abs(coriop(i,j)))
            Kv_col = Kv_col*q/cori30*log(2.*bvf0/q)/log(2.*bvf0/cori30)
            Kd_col = Kd_col*q/cori30*log(2.*bvf0/q)/log(2.*bvf0/cori30)
          end if

          !- Tidally driven diapycnal mixing

          if (tdmflg == 1) then
            call CVMix_init_tidal( &
                 CVmix_tidal_params_user = CVMix_tidal_params, &
                 mix_scheme = 'Simmons', &
                 efficiency=dmxeff, local_mixing_frac = tdmq)

            call CVMix_compute_Simmons_invariant(nlev = kk, &
                 energy_flux = twedon(i,j)*bvfbot, &
                 rho = CVMix_glb_params%FreshWaterDensity, &
                 SimmonsCoeff = Simmons_coeff, &
                 VertDep = vert_dep(:), &
                 zw = iFaceHeight(:), zt = cellHeight(:), &
                 CVmix_tidal_params_user = CVMix_tidal_params)

            call CVMix_coeffs_tidal(Mdiff_out = Kv_tidal(:), &
                 Tdiff_out = Kd_tidal(:), Nsqr = bvfsq_i(:), &
                 OceanDepth = -iFaceHeight(kk+1), &
                 SimmonsCoeff = Simmons_coeff, &
                 vert_dep = vert_dep(:), &
                 nlev=kk, max_nlev = kk, &
                 cvmix_params = CVMix_glb_params, &
                 CVmix_tidal_params_user = CVMix_tidal_params)
          else
            Kd_tidal(:) = 0.
          end if

          ! Call to CVMix wrapper for computing interior mixing coefficients.
          call  CVMix_coeffs_shear(Mdiff_out = Kv_shr(:), &
               Tdiff_out = Kd_shr(:), &
               RICH = rig_i(:), &
               nlev = kk, &
               max_nlev = kk)

          ! gets index of the level and interface above hbl
          kOBL = int(hOBL(i,j))    ! index of interface above OBL depth

          !- Diapycnal mixing when local stability is weak
          !- convection routine based on N2 not rho
          !- make sure it is in metrics if stability depends on rho
          call CVMix_coeffs_conv(Mdiff_out = Kv_conv(:), &
               Tdiff_out = Kd_conv(:), Nsqr = bvfsq_i(:), &
               dens=rho_zeros(:),dens_lwr = rho_lwr(:), &
               nlev=kk, max_nlev = kk, &
               OBL_ind = kOBL)
          ! Do not apply mixing due to convection within the boundary layer
          do k = 1,kOBL
            Kv_conv(k) = 0.0
            Kd_conv(k) = 0.0
          end do

          ! total diffusivities without KPP
          Kv_kpp(:) = Kv_col(:)+Kv_conv(:)+Kv_shr(:)
          Kt_kpp(:) = Kd_col(:)+Kd_conv(:)+Kd_shr(:)+Kd_tidal(:)
          Ks_kpp(:) = Kd_col(:)+Kd_conv(:)+Kd_shr(:)+Kd_tidal(:)

          ! Buoyancy flux acting on the OBL
          surfBuoyFlux = buoyfl(i,j,kOBL+1) - buoyfl(i,j,1  )

          ! Compute KPP using CVMix
          call CVMix_coeffs_kpp(&
               Kv_kpp(:),                         &  ! (inout) Total viscosity [m2 s-1]
               Kt_kpp(:),                         &  ! (inout) Total temp diffusivity [m2 s-1]
               Ks_kpp(:),                         &  ! (inout) Total salt  diffusivity [m2 s-1]
               iFaceHeight(:),                    &  ! (in)  Height of interfaces [m]
               cellHeight(:),                     &  ! (in)   Height of level centers [m]
               Kv_kpp(:),                         &  ! (in) Original viscosity [m2 s-1]
               Kt_kpp(:),                         &  ! (in) Original temp diffusivity [m2 s-1]
               Ks_kpp(:),                         &  ! (in) Original salt  diffusivity [m2 s-1]
               OBLdepth(i,j),                     &  ! (in) OBL depth [m]
               hOBL(i,j),                         &  ! (in) level (+fraction) of OBL extent
               nonLocalTrans(:,1),                &  ! (out) Non-local heat transport [nondim]
               nonLocalTrans(:,2),                &  ! (out) Non-local salt transport [nondim]
               surfFricVel,                       &  ! (in)  Turbulent friction  velocity at surface [m s-1]
               surfBuoyFlux,                      &  ! (in) Buoyancy flux at surface [m2 s-3]
               kk,                                &  ! (in) Number of levels to compute coeffs for
               kk,                                &  ! (in) Number of levels in array shape
               lamult(i,j),                       &  ! (in) Langmuir enhancement factor []
               CVMix_kpp_params_user = KPP_params )  ! KPP parameters

          !- ccc -------
          Kv_kpp = Kv_kpp
          Kt_kpp = Kt_kpp
          Ks_kpp = Ks_kpp
          Kv_kpp = max(nubmin,Kv_kpp)
          Kt_kpp = max(nubmin,Kt_kpp)
          Ks_kpp = max(nubmin,Ks_kpp)
          Kvisc_m(i,j,:) = Kv_kpp(:)
          Kdiff_t(i,j,:) = Kt_kpp(:)
          Kdiff_s(i,j,:) = Ks_kpp(:)
          t_ns_nonloc(i,j,:) = nonLocalTrans(:,1)
          s_nb_nonloc(i,j,:) = nonLocalTrans(:,2)
          do k = 1, kk+1
            t_sw_nonloc(i,j,k) = max(t_sw_nonloc(i,j,k), &
                 nonLocalTrans(k,1))
          end do

          ! Compute convective velocity cubed [m3 s-3].
          wstar3(i,j) = max(0.,-surfBuoyFlux)*OBLdepth(i,j)

        end do
      end do
      ! end of single column

    end do ! j-index

    do j = 1, jj
      do l = 1, isu(j)
        do i = max(1, ifu(j,l)), min(ii, ilu(j,l))
          dps = 0.
          do k = 2, kk+1
            dps = dps + dpu(i,j,k-1+mm)
            q = 2.*dps/((OBLdepth(i-1,j) + OBLdepth(i,j))*onem)
            if (q < 1.) then
              mu_nonloc(i,j,k) = (1. - q)**2
            else
              mu_nonloc(i,j,k) = 0.
            end if
          end do
        end do
      end do
      do l = 1, isv(j)
        do i = max(1, ifv(j,l)), min(ii, ilv(j,l))
          dps = 0.
          do k = 2, kk+1
            dps = dps + dpv(i,j,k-1+mm)
            q = 2.*dps/((OBLdepth(i,j-1) + OBLdepth(i,j))*onem)
            if (q < 1.) then
              mv_nonloc(i,j,k) = (1. - q)**2
            else
              mv_nonloc(i,j,k) = 0.
            end if
          end do
        end do
      end do
    end do

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'difest_vertical_hyb:'
      end if
      call chksummsk(Kvisc_m,ip,kk+1,'Kvisc_m')
      call chksummsk(Kdiff_t,ip,kk+1,'Kdiff_t')
      call chksummsk(Kdiff_s,ip,kk+1,'Kdiff_s')
      call chksummsk(t_ns_nonloc,ip,kk+1,'t_ns_nonloc')
      call chksummsk(s_nb_nonloc,ip,kk+1,'s_nb_nonloc')
      call chksummsk(mu_nonloc,iu,kk+1,'mu_nonloc')
      call chksummsk(mv_nonloc,iv,kk+1,'mv_nonloc')
    end if

  end subroutine difest_vertical_hyb

  subroutine difest_lateral_hyb(m,n,mm,nn,k1m,k1n)

    ! -----------------------------------------------------------
    ! estimate layer interface, isopycnal, diffusivities for hybrid
    ! coordinates
    ! -----------------------------------------------------------

    integer :: m,n,mm,nn,k1m,k1n


    real, dimension(1-nbdy:idm+nbdy,kdm) :: egr,anisok
    real, dimension(1-nbdy:idm+nbdy) :: &
         tup,pup,sup,cr,bcrrd,afeql,dps,egrs,egrup,dfints,anisos,ubt,vbt, &
         umls,vmls,urmse,cpse
    integer :: i,j,k,l,kn
    real :: q,plo,tlo,slo,tsfac,rhisc,ubc,vbc,speed,rhisct,falign, &
         els,egrlo,umnsc,esfac

    ! Locate the range of layers to be considered in the computation of
    ! diffusivities.
    do j = 0,jj+1
      do i = 0,ii+1
        kmax(i,j) = 0
      end do
      do l = 1,isp(j)
        do i = max(0,ifp(j,l)),min(ii+1,ilp(j,l))
          kmax(i,j) = 1
          do k = 2,kk
            kn = k+nn
            if (dp(i,j,kn) > dpbmin) kmax(i,j) = k
          end do
        end do
      end do
    end do
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          kfil(i,j) = kk+1
          do k = kk,2,-1
            if (p(i,j,k) > mlts(i,j)*(onem)) kfil(i,j) = k
          end do
        end do
      end do
    end do

    !$omp parallel do private( &
    !$omp l,i,k,kn,q,tup,pup,sup,cr,plo,tlo,slo,bcrrd,afeql,dps,egrs,egr, &
    !$omp egrup,egrlo,dfints,anisos,ubt,vbt,rhisc,ubc,vbc,speed,rhisct, &
    !$omp falign,els,anisok,umls,vmls,urmse,cpse,umnsc,esfac)
    do j = 1,jj

      ! Compute the first baroclinic rossby radius of deformation using
      ! the WKB approximation by Chelton at al. (1998).
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          pup(i) = p(i,j,1)
          kn = 1+nn
          tup(i) = temp(i,j,kn)
          sup(i) = saln(i,j,kn)
          cr(i) = 0.
        end do
      end do
      do k = 2,kk
        kn = k+nn
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (p(i,j,kk+1)-p(i,j,k+1) < epsilp) then
              plo = p(i,j,kk+1)
            else
              plo = .5*(p(i,j,k)+p(i,j,k+1))
            end if
            tlo = temp(i,j,kn)
            slo = saln(i,j,kn)
            cr(i) = cr(i) &
                 +sqrt(max(0.,(rho(p(i,j,k),tlo,slo) &
                              -rho(p(i,j,k),tup(i),sup(i))) &
                 *(plo-pup(i))))
            pup(i) = plo
            tup(i) = tlo
            sup(i) = slo
          end do
        end do
      end do
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          cr(i) = alpha0*cr(i)/pi
          bcrrd(i)= &
               sqrt(cr(i)*cr(i) &
               /max(coriop(i,j)*coriop(i,j)+2.*betafp(i,j)*cr(i), &
               1.e-24))
          afeql(i) = max(abs(coriop(i,j)),sqrt(2.*betafp(i,j)*cr(i)))
        end do
      end do

      ! - Compute diffusivity weigth to reduce eddy diffusivity when the
      ! - Rossby radius is resolved by the grid.
      if     (edwmth_opt == edwmth_smooth) then
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            q = bcrrd(i)/sqrt(.5*(scpx(i,j)*scpx(i,j) &
                 +scpy(i,j)*scpy(i,j)))
            difwgt(i,j) = 1./(1.+.25*q**4)
          end do
        end do
      else if (edwmth_opt == edwmth_step) then
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            q = bcrrd(i)/sqrt(.5*(scpx(i,j)*scpx(i,j) + scpy(i,j)*scpy(i,j)))
            if (q <= 2.) then
              difwgt(i,j) = 1.
            else
              difwgt(i,j) = 0.
            end if
          end do
        end do
      end if

      ! -----------------------------------------------------------
      ! Compute layer interface and isopycnal diffusivities
      ! -----------------------------------------------------------

      if (iidtyp == 1) then

        ! Type 1: Diffusivities are diffusive velocities multiplied by
        ! the local horizontal grid scale.
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            q = sqrt(scp2(i,j))
            difint(i,j,1) = thkdff*q
            difiso(i,j,1) = temdff*q
          end do
        end do
        do k = 2,kk
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              difint(i,j,k) = difint(i,j,1)
              difiso(i,j,k) = difiso(i,j,1)
            end do
          end do
        end do

      else

        ! Type 2: Diffusivities are parameterized according to Eden and
        ! Greatbatch (2008), extended with options for vertically averaged
        ! diffusivities and various eddy supression approaches.


        ! Eady growth rate.
        if (edsprs.or.edanis) then
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              egrs(i) = 0.
              dps(i) = 0.
            end do
          end do
        end if
        if     (edritp_opt == edritp_shear) then
          do k = 2,kk
            do l = 1,isp(j)
              do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                if (k >= kfil(i,j).and.k <= kmax(i,j).and. &
                     kmax(i,j)-kfil(i,j) >= 1) then
                  egr(i,k) = afeql(i) &
                       /sqrt(.5*(rig(i,j,k)+rig(i,j,k+1))+eggam)
                  if (edsprs.or.edanis) then
                    if (eddf2d) then
                      q = max(0.,p(i,j,k+1)-p(i,j,k))
                    else
                      q = max(0.,min(p(i,j,kfil(i,j))+dpgrav, &
                           p(i,j,k+1))-p(i,j,k))
                    end if
                    dps(i) = dps(i)+q
                    egrs(i) = egrs(i)+egr(i,k)*q
                  end if
                end if
              end do
            end do
          end do
        else if (edritp_opt == edritp_large_scale) then
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (kmax(i,j)-kfil(i,j) >= 1) then
                k = kfil(i,j)
                if     (kmax(i-1,j) >= k.and.kmax(i+1,j) >= k) then
                  q = .25*(nnslpx(i,j,k)+nnslpx(i+1,j,k))**2
                else if (kmax(i-1,j) >= k) then
                  q = nnslpx(i,j,k)**2
                else if (kmax(i+1,j) >= k) then
                  q = nnslpx(i+1,j,k)**2
                else
                  q = 0.
                end if
                if     (kmax(i,j-1) >= k.and.kmax(i,j+1) >= k) then
                  q = q+.25*(nnslpy(i,j,k)+nnslpy(i,j+1,k))**2
                else if (kmax(i,j-1) >= k) then
                  q = q+nnslpy(i,j,k)**2
                else if (kmax(i,j+1) >= k) then
                  q = q+nnslpy(i,j+1,k)**2
                end if
                egrup(i) = sqrt(q)
              end if
            end do
          end do
          do k = 2,kk
            do l = 1,isp(j)
              do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                if (kmax(i,j)-kfil(i,j) >= 1) then
                  if     (k >= kfil(i,j).and.k < kmax(i,j)) then
                    if     (kmax(i-1,j) > k.and.kmax(i+1,j) > k) then
                      q = .25*(nnslpx(i,j,k+1)+nnslpx(i+1,j,k+1))**2
                    else if (kmax(i-1,j) > k) then
                      q = nnslpx(i,j,k+1)**2
                    else if (kmax(i+1,j) > k) then
                      q = nnslpx(i+1,j,k+1)**2
                    else
                      q = 0.
                    end if
                    if     (kmax(i,j-1) > k.and.kmax(i,j+1) > k) then
                      q = q+.25*(nnslpy(i,j,k+1)+nnslpy(i,j+1,k+1))**2
                    else if (kmax(i,j-1) > k) then
                      q = q+nnslpy(i,j,k+1)**2
                    else if (kmax(i,j+1) > k) then
                      q = q+nnslpy(i,j+1,k+1)**2
                    end if
                    egrlo = sqrt(q)
                    egr(i,k) = .5*(egrup(i)+egrlo)
                    egrup(i) = egrlo
                    if (edsprs.or.edanis) then
                      if (eddf2d) then
                        q = max(0.,p(i,j,k+1)-p(i,j,k))
                      else
                        q = max(0.,min(p(i,j,kfil(i,j))+dpgrav, &
                                       p(i,j,k+1))-p(i,j,k))
                      end if
                      dps(i) = dps(i)+q
                      egrs(i) = egrs(i)+egr(i,k)*q
                    end if
                  else if (k == kmax(i,j)) then
                    egr(i,k) = egr(i,k-1)
                  end if
                end if
              end do
            end do
          end do
        end if
        if (edsprs.or.edanis) then
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (dps(i) > 0.) then
                egrs(i) = egrs(i)/dps(i)
              else
                egrs(i) = 0.
              end if
            end do
          end do
        end if

        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            difint(i,j,1) = egmndf
            dfints(i) = 0.
            dps(i) = 0.
            anisos(i) = 0.
          end do
        end do

        if (rhsctp) then

          ! Obtain barotropic velocities at p-points. TODO: check if
          ! weighting by bottom pressure is appropriate.
          tsfac = dlt/delt1
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              ubt(i) = (ubflxs_p(i  ,j,n)*scuyi(i  ,j) &
                       +ubflxs_p(i+1,j,n)*scuyi(i+1,j))*tsfac &
                       /max(epsilp,pbu(i,j,n)+pbu(i+1,j,n))
              vbt(i) = (vbflxs_p(i,j  ,n)*scvxi(i,j  ) &
                       +vbflxs_p(i,j+1,n)*scvxi(i,j+1))*tsfac &
                       /max(epsilp,pbv(i,j,n)+pbv(i,j+1,n))
            end do
          end do
        end if

        do k = 2,kk
          kn = k+nn
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (k >= kfil(i,j).and.k <= kmax(i,j).and. &
                   kmax(i,j)-kfil(i,j) >= 1) then

                ! Planetary Rhines scale.
                rhisc = egr(i,k)/max(1.e-22,betafp(i,j))

                if (edanis.or.rhsctp) then

                  ! Baroclinic velocities at p-points. TODO: check if
                  ! weighting by layer thickness is appropriate.
                  ubc = (u(i,j,kn)*dpu(i,j,kn)+u(i+1,j,kn)*dpu(i+1,j,kn)) &
                       /max(epsilp,dpu(i,j,kn)+dpu(i+1,j,kn))
                  vbc = (v(i,j,kn)*dpv(i,j,kn)+v(i,j+1,kn)*dpv(i,j+1,kn)) &
                       /max(epsilp,dpv(i,j,kn)+dpv(i,j+1,kn))

                  ! Current speed.
                  speed = max(1.e-22,sqrt(ubc*ubc+vbc*vbc))

                  if (rhsctp) then

                    ! Topographic Rhines scale.
                    rhisct = egr(i,k)/max(1.e-22,betatp(i,j))

                    ! Mask rhsctp if flow is not along bottom topography
                    !  1) option 1 with cosine to power of 10 will be >10
                    !     from +/- 35 degrees if rhiscf=1.
                    !  2) option 2 with hyberpolic tangent to power of
                    !     will be >10 from +/- 22.5 degrees if rhiscf=0.1
                    !     falign=1.
                    falign = 1./max(sin(atan2(vbc+vbt(i),ubc+ubt(i))&
                         - hangle(i,j))**10,1.e-10)
                    ! falign=1./max((1.-tanh(abs( &
                    !       abs(atan((vbc+vbt(i))/(ubc+ubt(i)))-hangle(i,j)) &
                    !       -pi/2.)))**6,1.e-24)

                    rhisc = min(rhisc,falign*rhiscf*rhisct)

                  end if

                end if

                ! Eddy length scale.
                els = max(eglsmn,min(bcrrd(i),rhisc))

                ! Temporary layer interface diffusivity.
                difint(i,j,k) = egc*egr(i,k)*els*els

                ! Accumulate diffusivities and estimate and accumulate
                ! anisotrophy if requested.

                if (eddf2d) then
                  q = max(0.,p(i,j,k+1)-p(i,j,k))
                else

                  ! Only consider a region below the first physical layer
                  ! with 3d structure of diffusivity.
                  q = max(0.,min(p(i,j,kfil(i,j))+dpdiav, &
                                 p(i,j,k+1))-p(i,j,k))
                end if

                dps(i) = dps(i)+q
                dfints(i) = dfints(i)+difint(i,j,k)*q

                if (edanis) then
                  anisok(i,k) = 1./(1.+(speed/max(1.e-22,&
                                                  egr(i,k)*els))**2)
                  anisos(i) = anisos(i)+anisok(i,k)*q
                end if

              else
                difint(i,j,k) = difint(i,j,k-1)
              end if
            end do
          end do
        end do

        ! Apply eddy diffusivity limiting, suppression when the Rossby
        ! radius is resolved by the grid, and suppression due to
        ! anisotropy or away from steering levels.

        ! Eddy diffusivity modification of surface non-isopycnic
        ! layers.

        if (edsprs.or.edanis) then

          ! Mixed layer velocity.
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              umls(i) = 0.
              vmls(i) = 0.
            end do
          end do
          do k = 1,kk
            do l = 1,isp(j)
              do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                q = max(0.,min(p(i,j,k+1),OBLdepth(i,j)*onem)-p(i,j,k))
                umls(i) = umls(i)+up(i,j,k)*q
                vmls(i) = vmls(i)+vp(i,j,k)*q
              end do
            end do
          end do
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              q = 1./(OBLdepth(i,j)*onem)
              umls(i) = umls(i)*q
              vmls(i) = vmls(i)*q
            end do
          end do
        end if

        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))

            if (edsprs.or.edanis) then

              ! Rhines scale using vertically averaged Eady growth rate.

              ! Planetary Rhines scale.
              rhisc = egrs(i)/max(1.e-22,betafp(i,j))

              if (rhsctp) then

                ! Topographic Rhines scale.
                rhisct = egrs(i)/max(1.e-22,betatp(i,j))
                rhisc = min(rhisc,rhiscf*rhisct)

              end if

              ! Eddy length scale.
              els = max(eglsmn,min(bcrrd(i),rhisc))

              if (edsprs) then

                ! RMS eddy velocity estimated from K = Gamma*u_rms*L,
                ! where a mixing efficiency of Gamma = 0.35 is used
                ! (Klocker and Abernathey, 2014). Note that 1/0.35 = 2.86.
                urmse(i) = 2.86*egc*egrs(i)*els

                ! Zonal eddy phase speed minus zonal barotropic velocity
                ! with a lower bound of -0.20 m s-1.
                cpse(i) = max(cpsemin,-betafp(i,j)*bcrrd(i)**2)

              end if

            end if

            if (dps(i) > 0.) then

              dfints(i) = dfints(i)/dps(i)

              if (edsprs) then

                ! Zonal mixed layer velocity minus eddy phase speed. Note
                ! that only the baroclinic component is used since the
                ! barotropic velocity is subtracted from the estimate of
                ! eddy phase speed.
                umnsc = umls(i)*cosang(i,j)-vmls(i)*sinang(i,j)-cpse(i)

                ! Eddy mixing suppresion factor where lower bounds of
                ! zonal velocity minus eddy phase speed and absolute value
                ! of RMS eddy velocity is set to -0.20 m s-1 and 0.05 m s-1,
                ! respectively.
                esfac = 1./(1.+4.*(umnsc/max(urmsemin, &
                                             abs(urmse(i))))**2)

              else if (edanis) then
                anisos(i) = anisos(i)/dps(i)
                speed = max(1.e-22, &
                            sqrt(umls(i)*umls(i)+vmls(i)*vmls(i)))
                esfac = 1./(1.+(speed/max(1.e-22,egrs(i)*els))**2)
              else
                esfac = 1.
              end if

              if (eddf2d) then

                ! GM is always 2D, will feel the influence of mean
                ! anisotropy.
                if (edanis) then
                  difint(i,j,1) = anisos(i)*dfints(i)*difwgt(i,j)
                else
                  difint(i,j,1) = dfints(i)*difwgt(i,j)
                end if
                if (redi3d) then

                  ! Allow Redi diffusivity to have a 3D profile.
                  difiso(i,j,1)= &
                       min(difmxp(i,j),egmxdf, &
                       max(egmndf,esfac*dfints(i)*egidfq*difwgt(i,j)))
                else
                  difiso(i,j,1)= &
                       min(difmxp(i,j),egmxdf, &
                       max(egmndf,difint(i,j,1)*egidfq))
                end if
                difint(i,j,1)= &
                     min(difmxp(i,j),egmxdf,max(egmndf,difint(i,j,1)))
              else
                difint(i,j,1)= &
                     min(difmxp(i,j),egmxdf, &
                     max(egmndf,dfints(i)*difwgt(i,j)*esfac))
                difiso(i,j,1) = difint(i,j,1)*egidfq
              end if

            else
              dfints(i) = egmndf
              difiso(i,j,1) = difint(i,j,1)*egidfq
            end if
          end do
        end do

        ! Eddy diffusivity modification of isopycnic layers.
        do k = 2,kk
          kn = k+nn
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (k >= kfil(i,j).and.k <= kmax(i,j).and. &
                   kmax(i,j)-kfil(i,j) >= 1) then

                if (edsprs) then

                  ! Zonal velocity minus eddy phase speed.
                  umnsc = up(i,j,k)*cosang(i,j)-vp(i,j,k)*sinang(i,j) &
                       -cpse(i)

                  ! Eddy mixing suppresion factor.
                  esfac = 1./ &
                       (1.+4.*(umnsc/max(urmsemin,abs(urmse(i))))**2)

                else if (edanis) then
                  esfac = anisok(i,k)
                else
                  esfac = 1.
                end if

                if (eddf2d) then

                  ! GM is always 2D, will feel the influence of mean
                  ! anisotropy.
                  difint(i,j,k) = difint(i,j,1)

                  if (redi3d) then

                    ! Allow Redi diffusivity to have a 3D profile.
                    difiso(i,j,k)= &
                         min(difmxp(i,j),egmxdf, &
                         max(egmndf,esfac*dfints(i)*egidfq &
                         *difwgt(i,j)))
                  else
                    difiso(i,j,k) = difiso(i,j,1)
                  end if
                else
                  difint(i,j,k)= &
                       min(difmxp(i,j),egmxdf, &
                       max(egmndf,difint(i,j,k)*difwgt(i,j)*esfac))
                  difiso(i,j,k) = difint(i,j,k)*egidfq
                end if
              else
                difint(i,j,k) = difint(i,j,k-1)
                difiso(i,j,k) = difiso(i,j,k-1)
              end if
            end do
          end do
        end do

      end if
    end do
    !$omp end parallel do

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'difest_lateral_hyb:'
      end if
      call chksummsk(difint,ip,kk,'difint')
      call chksummsk(difiso,ip,kk,'difiso')
    end if

  end subroutine difest_lateral_hyb

  subroutine difest_lateral_iso(m,n,mm,nn,k1m,k1n)

    !-----------------------------------------------------------
    ! estimate layer interface, isopycnal, diffusivities for isopycnal
    ! coordinates
    !-----------------------------------------------------------

    integer :: m,n,mm,nn,k1m,k1n

    real, dimension(1-nbdy:idm+nbdy,kdm) :: egr,anisok
    real, dimension(1-nbdy:idm+nbdy) :: &
         tup,pup,sup,cr,bcrrd,afeql,dps,egrs,egrup,dfints,anisos,ubt,vbt, &
         urmse,cpse
    integer :: i,j,k,l,kn
    real :: q,plo,tlo,slo,tsfac,rhisc,ubc,vbc,speed,rhisct,falign, &
         els,egrlo,umnsc,esfac

    !$omp parallel do private( &
    !$omp l,i,k,kn,q,tup,pup,sup,cr,plo,tlo,slo,bcrrd,afeql,dps,egrs,egr, &
    !$omp egrup,egrlo,dfints,anisos,ubt,vbt,rhisc,ubc,vbc,speed,rhisct, &
    !$omp falign,els,anisok,urmse,cpse,umnsc,esfac)
    do j = 1,jj

      ! Compute the first baroclinic rossby radius of deformation using
      ! the WKB approximation by Chelton at al. (1998).
      ! !!! Could include top layer in computation !!!
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          pup(i) = .5*(3.*p(i,j,3)-p(i,j,min(kk,kfpla(i,j,n))+1))
          kn = 2+nn
          tup(i) = temp(i,j,kn)
          sup(i) = saln(i,j,kn)
          cr(i) = 0.
        end do
      end do
      do k = 3,kk
        kn = k+nn
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k >= kfpla(i,j,n)) then
              if (p(i,j,kk+1)-p(i,j,k+1) < epsilp) then
                plo = p(i,j,kk+1)
              else
                plo = .5*(p(i,j,k)+p(i,j,k+1))
              end if
              tlo = temp(i,j,kn)
              slo = saln(i,j,kn)
              cr(i) = cr(i) &
                   +sqrt(max(0.,(rho(p(i,j,k),tlo,slo) &
                   -rho(p(i,j,k),tup(i),sup(i))) &
                   *(plo-pup(i))))
              pup(i) = plo
              tup(i) = tlo
              sup(i) = slo
            end if
          end do
        end do
      end do
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          cr(i) = alpha0*cr(i)/pi
          bcrrd(i)= &
               sqrt(cr(i)*cr(i) &
               /max(coriop(i,j)*coriop(i,j)+2.*betafp(i,j)*cr(i), &
               1.e-24))
          afeql(i) = max(abs(coriop(i,j)),sqrt(2.*betafp(i,j)*cr(i)))
        end do
      end do

      ! - Compute diffusivity weigth to reduce eddy diffusivity when the
      ! - Rossby radius is resolved by the grid.
      if     (edwmth_opt == edwmth_smooth) then
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            q = bcrrd(i)/sqrt(.5*(scpx(i,j)*scpx(i,j) &
                                 +scpy(i,j)*scpy(i,j)))
            difwgt(i,j) = 1./(1.+.25*q**4)
          end do
        end do
      else if (edwmth_opt == edwmth_step) then
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            q = bcrrd(i)/sqrt(.5*(scpx(i,j)*scpx(i,j) &
                                 +scpy(i,j)*scpy(i,j)))
            if (q <= 2.) then
              difwgt(i,j) = 1.
            else
              difwgt(i,j) = 0.
            end if
          end do
        end do
      end if

      !-----------------------------------------------------------
      ! Compute layer interface and isopycnal diffusivities
      !-----------------------------------------------------------

      if (iidtyp == 1) then

        ! Type 1: Diffusivities are diffusive velocities multiplied by
        ! the local horizontal grid scale.
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            q = sqrt(scp2(i,j))
            difint(i,j,1) = thkdff*q
            difiso(i,j,1) = temdff*q
          end do
        end do
        do k = 2,kk
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              difint(i,j,k) = difint(i,j,1)
              difiso(i,j,k) = difiso(i,j,1)
            end do
          end do
        end do

      else

        ! Type 2: Diffusivities are parameterized according to Eden and
        ! Greatbatch (2008), extended with options for vertically averaged
        ! diffusivities and various eddy supression approaches.

        ! Eady growth rate.
        if (edsprs.or.edanis) then
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              egrs(i) = 0.
              dps(i) = 0.
            end do
          end do
        end if
        if     (edritp_opt == edritp_shear) then
          do k = 2,kk
            do l = 1,isp(j)
              do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                if (k >= kfil(i,j).and.k <= kmax(i,j).and. &
                     kmax(i,j)-kfil(i,j) >= 1) then
                  egr(i,k) = afeql(i)/sqrt(rig(i,j,k)+eggam)
                  if (edsprs.or.edanis) then
                    if (eddf2d) then
                      q = max(0.,p(i,j,k+1)-p(i,j,k))
                    else
                      q = max(0.,min(p(i,j,kfil(i,j))+dpgrav, &
                                     p(i,j,k+1))-p(i,j,k))
                    end if
                    dps(i) = dps(i)+q
                    egrs(i) = egrs(i)+egr(i,k)*q
                  end if
                end if
              end do
            end do
          end do
        else if (edritp_opt == edritp_large_scale) then
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (kmax(i,j)-kfil(i,j) >= 1) then
                k = kfil(i,j)
                if     (kmax(i-1,j) >= k.and.kmax(i+1,j) >= k) then
                  q = .25*(nnslpx(i,j,k)+nnslpx(i+1,j,k))**2
                else if (kmax(i-1,j) >= k) then
                  q = nnslpx(i,j,k)**2
                else if (kmax(i+1,j) >= k) then
                  q = nnslpx(i+1,j,k)**2
                else
                  q = 0.
                end if
                if     (kmax(i,j-1) >= k.and.kmax(i,j+1) >= k) then
                  q = q+.25*(nnslpy(i,j,k)+nnslpy(i,j+1,k))**2
                else if (kmax(i,j-1) >= k) then
                  q = q+nnslpy(i,j,k)**2
                else if (kmax(i,j+1) >= k) then
                  q = q+nnslpy(i,j+1,k)**2
                end if
                egrup(i) = sqrt(q)
              end if
            end do
          end do
          do k = 2,kk
            do l = 1,isp(j)
              do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
                if (kmax(i,j)-kfil(i,j) >= 1) then
                  if     (k >= kfil(i,j).and.k < kmax(i,j)) then
                    if     (kmax(i-1,j) > k.and.kmax(i+1,j) > k) then
                      q = .25*(nnslpx(i,j,k+1)+nnslpx(i+1,j,k+1))**2
                    else if (kmax(i-1,j) > k) then
                      q = nnslpx(i,j,k+1)**2
                    else if (kmax(i+1,j) > k) then
                      q = nnslpx(i+1,j,k+1)**2
                    else
                      q = 0.
                    end if
                    if     (kmax(i,j-1) > k.and.kmax(i,j+1) > k) then
                      q = q+.25*(nnslpy(i,j,k+1)+nnslpy(i,j+1,k+1))**2
                    else if (kmax(i,j-1) > k) then
                      q = q+nnslpy(i,j,k+1)**2
                    else if (kmax(i,j+1) > k) then
                      q = q+nnslpy(i,j+1,k+1)**2
                    end if
                    egrlo = sqrt(q)
                    egr(i,k) = .5*(egrup(i)+egrlo)
                    egrup(i) = egrlo
                    if (edsprs.or.edanis) then
                      if (eddf2d) then
                        q = max(0.,p(i,j,k+1)-p(i,j,k))
                      else
                        q = max(0.,min(p(i,j,kfil(i,j))+dpgrav, &
                                       p(i,j,k+1))-p(i,j,k))
                      end if
                      dps(i) = dps(i)+q
                      egrs(i) = egrs(i)+egr(i,k)*q
                    end if
                  else if (k == kmax(i,j)) then
                    egr(i,k) = egr(i,k-1)
                  end if
                end if
              end do
            end do
          end do
        end if
        if (edsprs.or.edanis) then
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (dps(i) > 0.) then
                egrs(i) = egrs(i)/dps(i)
              else
                egrs(i) = 0.
              end if
            end do
          end do
        end if

        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            difint(i,j,1) = egmndf
            dfints(i) = 0.
            dps(i) = 0.
            anisos(i) = 0.
          end do
        end do

        if (rhsctp) then

          ! Obtain barotropic velocities at p-points. TODO: check if
          ! weighting by bottom pressure is appropriate.
          tsfac = dlt/delt1
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              ubt(i) = (ubflxs_p(i  ,j,n)*scuyi(i  ,j) &
                       +ubflxs_p(i+1,j,n)*scuyi(i+1,j))*tsfac &
                       /max(epsilp,pbu(i,j,n)+pbu(i+1,j,n))
              vbt(i) = (vbflxs_p(i,j  ,n)*scvxi(i,j  ) &
                       +vbflxs_p(i,j+1,n)*scvxi(i,j+1))*tsfac &
                       /max(epsilp,pbv(i,j,n)+pbv(i,j+1,n))
            end do
          end do
        end if

        do k = 2,kk
          kn = k+nn
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (k >= kfil(i,j).and.k <= kmax(i,j).and. &
                   kmax(i,j)-kfil(i,j) >= 1) then

                ! Planetary Rhines scale.
                rhisc = egr(i,k)/max(1.e-22,betafp(i,j))

                if (edanis.or.rhsctp) then

                  ! Baroclinic velocities at p-points. TODO: check if
                  ! weighting by layer thickness is appropriate.
                  ubc = (u(i,j,kn)*dpu(i,j,kn)+u(i+1,j,kn)*dpu(i+1,j,kn)) &
                       /max(epsilp,dpu(i,j,kn)+dpu(i+1,j,kn))
                  vbc = (v(i,j,kn)*dpv(i,j,kn)+v(i,j+1,kn)*dpv(i,j+1,kn)) &
                       /max(epsilp,dpv(i,j,kn)+dpv(i,j+1,kn))

                  ! Current speed.
                  speed = max(1.e-22,sqrt(ubc*ubc+vbc*vbc))

                  if (rhsctp) then

                    !- -------- Topographic Rhines scale.
                    rhisct = egr(i,k)/max(1.e-22,betatp(i,j))

                    !- --- ---- Mask rhsctp if flow is not along bottom topography
                    !                   1) option 1 with cosine to power of 10 will be >10
                    !                      from +/- 35 degrees if rhiscf=1.
                    !                   2) option 2 with hyberpolic tangent to power of
                    !                      will be >10 from +/- 22.5 degrees if rhiscf=0.1
                    !                      falign=1.
                    falign = 1./max(sin(atan2(vbc+vbt(i),ubc+ubt(i)) &
                         -hangle(i,j))**10,1.e-10)
                    !                   falign=1./max((1.-tanh(abs(
                    !     .                    abs(atan((vbc+vbt(i))/(ubc+ubt(i)))-hangle(i,j))
                    !     .                    -pi/2.)))**6,1.e-24)

                    rhisc = min(rhisc,falign*rhiscf*rhisct)

                  end if

                end if

                ! Eddy length scale.
                els = max(eglsmn,min(bcrrd(i),rhisc))

                ! Temporary layer interface diffusivity.
                difint(i,j,k) = egc*egr(i,k)*els*els

                ! Accumulate diffusivities and estimate and accumulate
                ! anisotrophy if requested.

                if (eddf2d) then
                  q = max(0.,p(i,j,k+1)-p(i,j,k))
                else

                  ! Only consider a region below the first physical layer
                  ! with 3d structure of diffusivity.
                  q = max(0.,min(p(i,j,kfil(i,j))+dpdiav, &
                                 p(i,j,k+1))-p(i,j,k))
                end if

                dps(i) = dps(i)+q
                dfints(i) = dfints(i)+difint(i,j,k)*q

                if (edanis) then
                  anisok(i,k) = 1./(1.+(speed/max(1.e-22,egr(i,k)*els))**2)
                  anisos(i) = anisos(i)+anisok(i,k)*q
                end if

              else
                difint(i,j,k) = difint(i,j,k-1)
              end if
            end do
          end do
        end do

        ! Apply eddy diffusivity limiting, suppression when the Rossby
        ! radius is resolved by the grid, and suppression due to
        ! anisotropy or away from steering levels.

        ! Eddy diffusivity modification of surface non-isopycnic
        ! layers.
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))

            if (edsprs.or.edanis) then

              ! Rhines scale using vertically averaged Eady growth rate.

              ! Planetary Rhines scale.
              rhisc = egrs(i)/max(1.e-22,betafp(i,j))

              if (rhsctp) then

                ! Topographic Rhines scale.
                rhisct = egrs(i)/max(1.e-22,betatp(i,j))
                rhisc = min(rhisc,rhiscf*rhisct)

              end if

              ! Eddy length scale.
              els = max(eglsmn,min(bcrrd(i),rhisc))

              if (edsprs) then

                ! RMS eddy velocity estimated from K = Gamma*u_rms*L,
                ! where a mixing efficiency of Gamma = 0.35 is used
                ! (Klocker and Abernathey, 2014). Note that 1/0.35 = 2.86.
                urmse(i) = 2.86*egc*egrs(i)*els

                ! Zonal eddy phase speed minus zonal barotropic velocity
                ! with a lower bound of -0.20 m s-1.
                cpse(i) = max(cpsemin,-betafp(i,j)*bcrrd(i)**2)

              end if

            end if

            if (dps(i) > 0.) then

              dfints(i) = dfints(i)/dps(i)

              if (edsprs.or.edanis) then

                ! Baroclinic velocities at p-points.
                if     (ip(i-1,j)+ip(i+1,j) == 2) then
                  ubc = .5*((u(i  ,j,1+nn)*dpu(i  ,j,1+nn) &
                            +u(i  ,j,2+nn)*dpu(i  ,j,2+nn)) &
                         /(dpu(i  ,j,1+nn)+dpu(i  ,j,2+nn)) &
                           +(u(i+1,j,1+nn)*dpu(i+1,j,1+nn) &
                            +u(i+1,j,2+nn)*dpu(i+1,j,2+nn)) &
                         /(dpu(i+1,j,1+nn)+dpu(i+1,j,2+nn)))
                else if (ip(i-1,j) == 1) then
                  ubc = (u(i  ,j,1+nn)*dpu(i  ,j,1+nn) &
                        +u(i  ,j,2+nn)*dpu(i  ,j,2+nn)) &
                     /(dpu(i  ,j,1+nn)+dpu(i  ,j,2+nn))
                else if (ip(i+1,j) == 1) then
                  ubc = (u(i+1,j,1+nn)*dpu(i+1,j,1+nn) &
                        +u(i+1,j,2+nn)*dpu(i+1,j,2+nn)) &
                     /(dpu(i+1,j,1+nn)+dpu(i+1,j,2+nn))
                else
                  ubc = 0.
                end if
                if     (ip(i,j-1)+ip(i,j+1) == 2) then
                  vbc = .5*((v(i,j  ,1+nn)*dpv(i,j  ,1+nn) &
                            +v(i,j  ,2+nn)*dpv(i,j  ,2+nn)) &
                         /(dpv(i,j  ,1+nn)+dpv(i,j  ,2+nn)) &
                           +(v(i,j+1,1+nn)*dpv(i,j+1,1+nn) &
                            +v(i,j+1,2+nn)*dpv(i,j+1,2+nn)) &
                         /(dpv(i,j+1,1+nn)+dpv(i,j+1,2+nn)))
                else if (ip(i,j-1) == 1) then
                  vbc = (v(i,j  ,1+nn)*dpv(i,j  ,1+nn) &
                        +v(i,j  ,2+nn)*dpv(i,j  ,2+nn)) &
                     /(dpv(i,j  ,1+nn)+dpv(i,j  ,2+nn))
                else if (ip(i,j+1) == 1) then
                  vbc = (v(i,j+1,1+nn)*dpv(i,j+1,1+nn) &
                        +v(i,j+1,2+nn)*dpv(i,j+1,2+nn)) &
                     /(dpv(i,j+1,1+nn)+dpv(i,j+1,2+nn))
                else
                  vbc = 0.
                end if

                if (edanis) then
                  anisos(i) = anisos(i)/dps(i)
                  speed = max(1.e-22,sqrt(ubc*ubc+vbc*vbc))
                  esfac = 1./(1.+(speed/max(1.e-22, egrs(i)*els))**2)
                else

                  ! Zonal mixed layer velocity minus eddy phase speed.
                  ! Note that only the baroclinic component is used since
                  ! the barotropic velocity is subtracted from the
                  ! estimate of eddy phase speed.
                  umnsc = ubc*cosang(i,j)
                  umnsc = umnsc-vbc*sinang(i,j)-cpse(i)

                  ! Eddy mixing suppression factor where lower bounds of
                  ! zonal velocity minus eddy phase speed and absolute
                  ! value of RMS eddy velocity is set to -0.20 m s-1 and
                  ! 0.05 m s-1, respectively.
                  esfac = 1./(1.+4.*(umnsc/max(urmsemin, abs(urmse(i))))**2)
                end if

              else
                esfac = 1.
              end if

              if (eddf2d) then

                ! GM is always 2D, will feel the influence of mean
                ! anisotropy.
                if (edanis) then
                  difint(i,j,1) = anisos(i)*dfints(i)*difwgt(i,j)
                else
                  difint(i,j,1) = dfints(i)*difwgt(i,j)
                end if
                if (redi3d) then

                  ! Allow Redi diffusivity to have a 3D profile.
                  difiso(i,j,1)= &
                       min(difmxp(i,j),egmxdf, &
                       max(egmndf,esfac*dfints(i)*egidfq*difwgt(i,j)))
                else
                  difiso(i,j,1)= &
                       min(difmxp(i,j),egmxdf, &
                       max(egmndf,difint(i,j,1)*egidfq))
                end if
                difint(i,j,1)= &
                     min(difmxp(i,j),egmxdf,max(egmndf,difint(i,j,1)))
              else
                difint(i,j,1)= &
                     min(difmxp(i,j),egmxdf, &
                     max(egmndf,dfints(i)*difwgt(i,j)*esfac))
                difiso(i,j,1) = difint(i,j,1)*egidfq
              end if

            else
              dfints(i) = egmndf
              difiso(i,j,1) = difint(i,j,1)*egidfq
            end if
          end do
        end do

        ! Eddy diffusivity modification of isopycnic layers.
        do k = 2,kk
          kn = k+nn
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (k >= kfil(i,j).and.k <= kmax(i,j).and. &
                   kmax(i,j)-kfil(i,j) >= 1) then

                if (edsprs) then

                  ! Zonal velocity minus eddy phase speed.
                  umnsc= &
                       (msku(i,j,k)*u(i,j,kn)+msku(i+1,j,k)*u(i+1,j,kn)) &
                       /max(1,msku(i,j,k)+msku(i+1,j,k))*cosang(i,j) &
                       -(mskv(i,j,k)*v(i,j,kn)+mskv(i,j+1,k)*v(i,j+1,kn)) &
                       /max(1,mskv(i,j,k)+mskv(i,j+1,k))*sinang(i,j) &
                       -cpse(i)

                  ! Eddy mixing suppression factor.
                  esfac = 1./(1.+4.*(umnsc/max(urmsemin, abs(urmse(i))))**2)

                else if (edanis) then
                  esfac = anisok(i,k)
                else
                  esfac = 1.
                end if

                if (eddf2d) then

                  ! GM is always 2D, will feel the influence of mean
                  ! anisotropy.
                  difint(i,j,k) = difint(i,j,1)

                  if (redi3d) then

                    ! Allow Redi diffusivity to have a 3D profile.
                    difiso(i,j,k)= &
                         min(difmxp(i,j),egmxdf, &
                         max(egmndf,esfac*dfints(i)*egidfq &
                         *difwgt(i,j)))
                  else
                    difiso(i,j,k) = difiso(i,j,1)
                  end if
                else
                  difint(i,j,k)= &
                       min(difmxp(i,j),egmxdf, &
                       max(egmndf,difint(i,j,k)*difwgt(i,j)*esfac))
                  difiso(i,j,k) = difint(i,j,k)*egidfq
                end if
              else
                difint(i,j,k) = difint(i,j,k-1)
                difiso(i,j,k) = difiso(i,j,k-1)
              end if
            end do
          end do
        end do

      end if
    end do
    !$omp end parallel do


    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'difest_lateral_iso:'
      end if
      call chksummsk(difint,ip,kk,'difint')
      call chksummsk(difiso,ip,kk,'difiso')
    end if

  end subroutine difest_lateral_iso

  subroutine difest_vertical_iso(m,n,mm,nn,k1m,k1n)

    !-----------------------------------------------------------
    ! estimate diapycnal diffusivities for isopycnal model
    !-----------------------------------------------------------

    integer :: m,n,mm,nn,k1m,k1n

    real, dimension(1-nbdy:idm+nbdy,kdm) :: bvfsq,bvf
    real, dimension(1-nbdy:idm+nbdy) :: bvfbot,dps,dfddsu,dfddsl
    integer :: i,j,k,l,kn
    real :: q,nus,nub,nut,nuls,vsf,nusm,ust,mols,h,sg,zeta,phis,ws
    real :: gls_c3,tke_prod,tke_buoy,tke_epsilon,ls_unlmt,ls_lmt,tke_q
    real :: Gm,Gh,Sm,Sh,cff,ql
    real :: gls_prod,gls_buoy,gls_diss,gls_q

    !$omp parallel do private( &
    !$omp l,i,k,kn,q,bvfbot,dps,bvfsq,bvf,dfddsu,dfddsl,nub,nus,ust,vsf, &
    !$omp nut,nuls,nusm,mols,h,sg,zeta,phis,ws, &
    !$omp gls_c3,tke_epsilon,tke_prod,tke_buoy,tke_q,ls_unlmt,ls_lmt,gh, &
    !$omp gm,cff,sm,sh,ql, &
    !$omp gls_prod,gls_buoy,gls_diss,gls_q)
    do j = 1,jj

      ! Compute Brunt-Vaisala frequency.
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          bvfbot(i) = 0.
          dps(i) = 0.
        end do
      end do
      do k = 4,kk
        kn = k+nn
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k >= kfil(i,j).and.k <= kmax(i,j).and. &
                 kmax(i,j)-kfil(i,j) >= 1) then

              ! Brunt-Vaisala frequency squared
              bvfsq(i,k) = grav*grav*max(drhomn,drhol(i,j,k)) &
                           /max(epsilp,dp(i,j,kn))

              ! Brunt-Vaisala frequency
              bvf(i,k) = sqrt(bvfsq(i,k))

              if (use_TRC .and. use_TKE) then
                if (dp(i,j,kn) > dpbmin) then
                  Buoy(i,j,k) = -difdia(i,j,k)*bvfsq(i,k)
                  h = max(onem,dp(i,j,kn))*alpha0/grav
                  !               h=max(onem*1e-8,dp(i,j,kn))*alpha0/grav
                  !               h=max(onemm,dp(i,j,kn))*alpha0/grav
                  Shear2(i,j,k) = max(1.e-13,du2l(i,j,k))/(h*h)
                  Prod(i,j,k) = difdia(i,j,k)*Pr_t*Shear2(i,j,k)
                else
                  Buoy(i,j,k) = 0.
                  Shear2(i,j,k) = 1.e-9
                  Prod(i,j,k) = 0.
                end if
              end if

              ! Accumulate Brunt-Vaisala frequency in a region near the
              ! bottom
              q = max(0.,p(i,j,k+1)-max(p(i,j,kk+1)-dpnbav,p(i,j,k)))
              if (q > 0.) then
                bvfbot(i) = bvfbot(i)+bvf(i,k)*q
                dps(i) = dps(i)+q
              end if
            end if
          end do
        end do
      end do
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          if (dps(i) > 0.) then
            bvfbot(i) = bvfbot(i)/dps(i)
          end if
        end do
      end do

      !-----------------------------------------------------------
      ! Compute diapycnal diffusivity.
      !-----------------------------------------------------------

      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          difdia(i,j,1) = nu0
          dfddsu(i) = 0.
          dfddsl(i) = 0.
          dps(i) = 0.
        end do
      end do
      do k = 2,kk
        kn = k+nn
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k >= kfil(i,j).and.k <= kmax(i,j).and. &
                 kmax(i,j)-kfil(i,j) >= 1) then

              ! Background diapycnal mixing.
              if     (bdmtyp == 1) then

                ! Type 1: Background diffusivity is a constant divided by
                ! Brunt-Vaisala frequency.
                nub = bdmc1/bvf(i,k)
              else if (bdmtyp == 2) then

                ! Type 2: Background diffusivity is a constant
                nub = bdmc2
              else
                nub = 0.
              end if
              if (iwdflg == 1) then
                nub = nub*(1.+(iwdfac-1.)*ficem(i,j))
              end if

              ! Latitude dependency of background diapycnal mixing
              if (bdmldp) then
                q = max(1.e-9,abs(coriop(i,j)))
                nub = nub*q/cori30*log(2.*bvf0/q)/log(2.*bvf0/cori30)
              end if

              nub = max(nubmin,nub)

              if (.not. use_TRC .or. .not. use_TKE) then
                ! Shear driven diapycnal mixing.
                if (rig(i,j,k) < ri0) then

                  ! Maximum diffusivity is increased near the bottom to
                  ! provide additional mixing of gravity currents.
                  q = (p(i,j,kk+1)-p(i,j,k)+.5*dp(i,j,kn)) &
                       /min(dpgc,.5*p(i,j,kk+1))
                  q = max(0.,1.-q*q)
                  q = q*q*q
                  nus = q*nug0+(1.-q)*nus0

                  ! Parameterization of diffusivity as a function of local
                  ! gradient richardson number.
                  q = rig(i,j,k)/ri0
                  q = max(0.,1.-q*q)
                  nus = nus*q*q*q
                else
                  nus = 0.
                end if
              else
                if (bvfsq(i,k) > 0.) then  ! stable stratification
                  gls_c3 = gls_c3minus
                else                        ! unstable stratification
                  gls_c3 = gls_c3plus
                end if
                if (.not. use_GLS) then
                  trc(i,j,kn,itrgls) = max((gls_c1*Prod(i,j,k) &
                       +gls_c3*Buoy(i,j,k))/gls_c2, &
                       gls_psi_min)
                end if
                tke_epsilon = cmu_fac2*trc(i,j,kn,itrtke)**(1.5+gls_m/gls_n) &
                     *trc(i,j,kn,itrgls)**(-1./gls_n)
                tke_prod = Prod(i,j,k)
                tke_buoy = Buoy(i,j,k)
                tke_Q = tke_epsilon/trc(i,j,kn,itrtke)
                if (use_GLS) then
                  gls_prod = (trc(i,j,kn,itrgls)/trc(i,j,kn,itrtke)) &
                       *gls_c1*Prod(i,j,k)
                  gls_buoy = (trc(i,j,kn,itrgls)/trc(i,j,kn,itrtke)) &
                       *gls_c3*Buoy(i,j,k)
                  gls_diss = (trc(i,j,kn,itrgls)/trc(i,j,kn,itrtke)) &
                       *gls_c2*tke_epsilon
                  gls_Q = gls_diss/trc(i,j,kn,itrgls)
                  if (gls_prod+gls_buoy >= 0.) then
                    trc(i,j,kn,itrgls)= &
                         (trc(i,j,kn,itrgls)+delt1*(gls_prod+gls_buoy)) &
                         /(1.+delt1*gls_Q)
                  else
                    trc(i,j,kn,itrgls)= &
                         (trc(i,j,kn,itrgls)+delt1*gls_prod) &
                         /(1.+delt1*(gls_Q-(gls_buoy/trc(i,j,kn,itrgls))))
                  end if
                  trc(i,j,kn,itrgls) = max(trc(i,j,kn,itrgls),gls_psi_min)
                  q = .56**(.5*gls_n)*gls_cmu0**gls_p &
                       *trc(i,j,kn,itrtke)**(gls_m+.5*gls_n) &
                       *bvf(i,k)**(-gls_n)
                  if (gls_n < 0.) then
                    trc(i,j,kn,itrgls) = max(trc(i,j,kn,itrgls),q)
                  else
                    trc(i,j,kn,itrgls) = min(trc(i,j,kn,itrgls),q)
                  end if
                end if ! use_GLS

                tke_epsilon = cmu_fac2*trc(i,j,kn,itrtke)**(1.5+gls_m/gls_n) &
                     *trc(i,j,kn,itrgls)**(-1./gls_n)
                tke_Q = tke_epsilon/trc(i,j,kn,itrtke)

                if (tke_prod+tke_buoy >= 0.) then
                  trc(i,j,kn,itrtke)= &
                       (trc(i,j,kn,itrtke)+delt1*(tke_prod+tke_buoy)) &
                       /(1.+delt1*tke_Q)
                else
                  trc(i,j,kn,itrtke)= &
                       (trc(i,j,kn,itrtke)+delt1*tke_prod) &
                       /(1.+delt1*(tke_Q-(tke_buoy/trc(i,j,kn,itrtke))))
                  trc(i,j,kn,itrtke) = max(trc(i,j,kn,itrtke),tke_min)
                end if

                ! Penetration of surface TKE below mixed layer.
                if (tkepf > 0.) then
                  if (dp(i,j,kn) < epsilp) then
                    q = exp(-p(i,j,k)/tkepls)
                  else
                    q = tkepls*(exp(-p(i,j,k  )/tkepls) &
                         -exp(-p(i,j,k+1)/tkepls))/dp(i,j,kn)
                  end if
                  trc(i,j,kn,itrtke) = trc(i,j,kn,itrtke) &
                       +67.83*tkepf*q*ustar(i,j)**2
                end if

                ! Set TKE and GLS to prescribed minimum values in surface
                ! mixed layers and thin layers
                if (dp(i,j,kn) < epsilp) then
                  trc(i,j,kn,itrtke) = tke_min
                  trc(i,j,kn,itrgls) = gls_psi_min
                end if
                trc(i,j,1+nn,itrtke) = tke_min
                trc(i,j,2+nn,itrtke) = tke_min
                trc(i,j,1+nn,itrgls) = gls_psi_min
                trc(i,j,2+nn,itrgls) = gls_psi_min

                ! Bottom Boundary Conditions
                if (k == kmax(i,j)) then
                  ust = max(ustarb(i,j),ustmin)
                  trc(i,j,kn,itrtke) = max(tke_min,(ust/gls_cmu0)**2)
                  if (use_GLS) then
                    trc(i,j,kn,itrgls) = max(gls_psi_min, &
                         (gls_cmu0**(gls_p-2.*gls_m)) &
                         *(ust**(2.*gls_m)) &
                         *(kappa)**gls_n)
                  end if
                end if

                Ls_unlmt = max(Ls_unlmt_min, &
                     cmu_fac1*trc(i,j,kn,itrgls)**(gls_exp1) &
                     *trc(i,j,kn,itrtke)**(-tke_exp1))

                if (bvfsq(i,k) > 0.) then  ! stable stratification
                  !               Ls_lmt=min(Ls_unlmt,
                  !    .                     sqrt(.56*trc(i,j,kn,itrtke)
                  !    .                          /max(bvfsq(i,k),1.e-10)))

                  Ls_lmt = min(Ls_unlmt,trc(i,j,kn,itrtke)**(-gls_m/gls_n) &
                       *trc(i,j,kn,itrgls)**gls_n)
                  !               Ls_lmt=Ls_unlmt
                else                        ! unstable stratification
                  Ls_lmt = Ls_unlmt
                end if

                ! Compute nondimensional stability functions for tracers
                ! (Sh) and momentum (Sm). Canuto-A
                Gh = min(gls_Gh0,-bvfsq(i,k)*Ls_lmt*Ls_lmt &
                     /(2.*trc(i,j,kn,itrtke)))
                Gh = min(Gh,(Gh-(Gh-gls_Ghcri)**2) &
                     /(Gh+gls_Gh0-2.*gls_Ghcri))
                Gh = max(Gh,gls_Ghmin)
                Gh = min(Gh,gls_Gh0)

                ! Compute shear number.
                Gm = (gls_b0/gls_fac6-gls_b1*Gh+gls_b3*gls_fac6*(Gh**2)) &
                     /(gls_b2-gls_b4*gls_fac6*Gh)
                Gm = min(Gm,Shear2(i,j,k)*Ls_lmt*Ls_lmt &
                     /(2.*trc(i,j,kn,itrtke)))

                ! Compute stability functions
                cff = gls_b0-gls_b1*gls_fac6*Gh+gls_b2*gls_fac6*Gm &
                     +gls_b3*gls_fac6**2*Gh**2-gls_b4*gls_fac6**2*Gh*Gm &
                     +gls_b5*gls_fac6**2*Gm*Gm
                Sm = (gls_s0-gls_s1*gls_fac6*Gh+gls_s2*gls_fac6*Gm)/cff
                Sh = (gls_s4-gls_s5*gls_fac6*Gh+gls_s6*gls_fac6*Gm)/cff
                Sm = max(Sm,0.)
                Sh = max(Sh,0.)

                ! Relate Canuto stability to BLOM notation
                Sm = Sm*cmu_fac3/gls_cmu0**3
                Sh = Sh*cmu_fac3/gls_cmu0**3

                ql = sqrt2*(Ls_lmt) &
                     *sqrt(trc(i,j,kn,itrtke))
                !             ql=sqrt2*.5*(Ls_lmt+L_scale(i,j,k))
                !    .           *sqrt(trc(i,j,kn,itrtke))

                !             nus=Sh*ql
                !             nus=min(0.1*ql,4.05*nug0)
                nus = min(Sh*ql,4.05*nug0)
                !             nus=Sh*(trc(i,j,k,itrtke)*trc(i,j,k,itrtke))
                !    .            /trc(i,j,k,itrgls)
                L_scale(i,j,k) = max(Ls_lmt,Ls_unlmt_min)
                if (use_GLS) then
                  ! Recompute gls based on limited length scale
                  trc(i,j,kn,itrgls)= &
                       max(gls_cmu0**gls_p*trc(i,j,kn,itrtke)**gls_m &
                       *L_scale(i,j,k)**gls_n,gls_psi_min)
                end if
              end if

              ! Tidally driven diapycnal mixing
              if (tdmflg == 1) then
                q = .5*(tanh(4.*(abs(plat(i,j))-tdclat)/tddlat-2.)+1.)
                q = (1.-q)*tdmls0+q*tdmls1
                if (dp(i,j,kn) < epsilp) then
                  vsf = exp(p(i,j,k)/q)/(q*(exp(p(i,j,kk+1)/q)-1.))
                else
                  vsf = (exp(p(i,j,k+1)/q)-exp(p(i,j,k)/q)) &
                       /(dp(i,j,kn)*(exp(p(i,j,kk+1)/q)-1.))
                end if
                nut = grav*tdmq*dmxeff*twedon(i,j)*bvfbot(i)*vsf/bvfsq(i,k)
              else
                nut = 0.
              end if

              ! Diapycnal mixing when local stability is weak
              if (drhol(i,j,k) < drho0) then
                q = drhol(i,j,k)/drho0
                q = max(0.,1.-q*q)
                nuls = nuls0*q*q*q
              else
                nuls = 0.
              end if

              ! Total diapycnal diffusivity.
              difdia(i,j,k) = nub+nus+nut+nuls

              ! Accumulate diffusivities in a region below the first
              ! physical layer
              q = max(0.,min(p(i,j,kfil(i,j))+dpddav,p(i,j,k+1))-p(i,j,k))
              dps(i) = dps(i)+q
              dfddsu(i) = dfddsu(i)+nub*q
              dfddsl(i) = dfddsl(i)+difdia(i,j,k)*q

            else
              difdia(i,j,k) = difdia(i,j,k-1)
              if (use_TRC .and. use_TKE) then
                !             trc(i,j,kn,itrtke)=tke_min
                !             L_scale(i,j,k)=Ls_unlmt_min
                trc(i,j,kn,itrtke) = trc(i,j,kn-1,itrtke)
                L_scale(i,j,k) = L_scale(i,j,k-1)
                if (use_GLS) then
                  !             trc(i,j,kn,itrgls)=gls_psi_min
                  trc(i,j,kn,itrgls) = trc(i,j,kn-1,itrgls)
                end if
              end if
            end if
          end do
        end do
      end do
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          if (dps(i) > 0.) then
            dfddsu(i) = dfddsu(i)/dps(i)
            dfddsl(i) = dfddsl(i)/dps(i)
          else
            dfddsu(i) = nu0
            dfddsl(i) = nu0
          end if
        end do
      end do
      do k = 2,kk-1
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k < kfil(i,j)) then
              if (k > 2.and.kfil(i,j) <= kk.and. &
                   p(i,j,min(kk,kfil(i,j)))-p(i,j,3) > epsilp) then
                q = .5*(p(i,j,k+1)+p(i,j,k))
                difdia(i,j,k) = ((q-p(i,j,3))*dfddsl(i) &
                                +(p(i,j,kfil(i,j))-q)*dfddsu(i)) &
                                /(p(i,j,kfil(i,j))-p(i,j,3))
              else
                difdia(i,j,k) = dfddsu(i)
              end if
            end if
          end do
        end do
      end do

      ! - Diapycnal diffusivity beneath mixed layer by dissipation of
      ! - energy originating from near-inertial waves.
      do k = 2,kk-1
        kn = k+nn
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (k <= kmax(i,j).and.kmax(i,j)-kfil(i,j) >= 1) then
              q = niwls
              if (k == 2.or.dp(i,j,kn) < epsilp) then
                vsf = exp((p(i,j,3)-p(i,j,k+1))/q) &
                     /(q*(1.-exp((p(i,j,3)-p(i,j,kk+1))/q)))
              else
                vsf = (exp((p(i,j,3)-p(i,j,k  ))/q) &
                      -exp((p(i,j,3)-p(i,j,k+1))/q)) &
                     /(dp(i,j,kn)*(1.-exp((p(i,j,3)-p(i,j,kk+1))/q)))
              end if
              nusm = grav*niwgf*(1.-niwbf)*niwlf*dmxeff*idkedt(i,j)*vsf &
                     /(alpha0*bvfsq(i,max(k,kfil(i,j))))
              difdia(i,j,k) = difdia(i,j,k)+nusm
            end if
          end do
        end do
      end do

      ! - Diffusivity at the lower interface of the top layer
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))

          ! Lower bounded friction velocity
          ust = max(ustmin,ustar(i,j))

          ! Monin-Obukhov length scale
          mols = ust**3/(kappa*sign(max(abs(buoyfl(i,j,1)),bfeps), &
                                    -buoyfl(i,j,1)))

          ! Mixed layer thickness
          h = (p(i,j,3)-p(i,j,1))/onem

          ! Dimensionless vertical coordinate in the boundary layer
          sg = (p(i,j,2)-p(i,j,1))/(p(i,j,3)-p(i,j,1))

          ! Velocity scale
          if (mols < 0.) then
            zeta = min(sleps,sg)*h/mols
            if (zeta > zetas) then
              phis = (1.-16.*zeta)**(-1./2.)
            else
              phis = (as-cs*zeta)**(-1./3.)
            end if
          else
            zeta = sg*h/mols
            phis = 1.+5.*zeta
          end if
          ws = kappa*ust/phis

          difdia(i,j,1) = h*ws*sg*(1.-sg)**2
        end do
      end do

    end do
    !$omp end parallel do

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'difest_vertical_iso:'
      end if
      call chksummsk(idkedt,ip,1,'idkedt')
      call chksummsk(difdia,ip,kk,'difdia')
      if (use_TRC .and. use_TKE) then
        call chksummsk(trc(1-nbdy,1-nbdy,1,itrtke),ip,2*kk,'tke')
        if (use_GLS) then
          call chksummsk(trc(1-nbdy,1-nbdy,1,itrgls),ip,2*kk,'gls_psi')
        end if
      end if
    end if

  end subroutine difest_vertical_iso

end module mod_difest
