!>
!! @par (c) Copyright
!! This software is provided under:
!!
!! The 3-Clause BSD License
!! SPDX short identifier: BSD-3-Clause
!! See https://opensource.org/licenses/BSD-3-Clause
!!
!! (c) Copyright 2016-2021 MPI-M, Joeran Maerz, Irene Stemmler;
!!     first published 2020
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions are met:
!!
!! 1. Redistributions of source code must retain the above copyright notice,
!!    this list of conditions and the following disclaimer.
!! 2. Redistributions in binary form must reproduce the above copyright notice,
!!    this list of conditions and the following disclaimer in the documentation
!!    and/or other materials provided with the distribution.
!! 3. Neither the name of the copyright holder nor the names of its contributors
!!    may be used to endorse or promote products derived from this software
!!    without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
!! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.[7]
!!
!!
!! -----------------------------------------------------------------------------
!! -----------------------------------------------------------------------------
!! @file mo_m4ago.F90
!! @brief Module for Marine Aggregates:
!!        The Microstructure, Multiscale, Mechanistic, Marine Aggregates
!!        in the Global Ocean (M4AGO) sinking scheme
!!
!! The mo_aggregates module contains routines to calculate:
!!      - aggregate properties
!!      - mean sinking velocity of aggregates
!!
!! See:
!! Maerz et al. 2020: Microstructure and composition of marine aggregates
!!                    as co-determinants for vertical particulate organic
!!                    carbon transfer in the global ocean.
!!                    Biogeosciences, 17, 1765-1803,
!!                    https://doi.org/10.5194/bg-17-1765-2020
!!
!! This module is written within the project:
!! Multiscale Approach on the Role of Marine Aggregates (MARMA)
!! funded by the Max Planck Society (MPG)
!!
!! @author: joeran maerz (joeran.maerz@mpimet.mpg.de), MPI-M, HH
!! 2019, June, revised by Irene Stemmler (refactoring, cleaning), MPI-M, HH
!!
!! 2023 adopted to iHAMOCC by joeran maerz, UiB, Bergen
!!
!! -----------------------------------------------------------------------------
!! -----------------------------------------------------------------------------
!!
!!



module mo_m4ago
  use mo_vgrid,       only: dp_min
  use mo_control_bgc, only: dtb, dtbgc,io_stdo_bgc
  use mo_param_bgc,   only: calcdens, claydens, opaldens, calcwei, opalwei, ropal
  use mo_carbch,      only: ocetra
  use mo_param1_bgc,  only: iopal, ifdust, icalc, idet

  implicit none

  private

  ! Public subroutines
  public :: mean_aggregate_sinking_speed, init_m4ago_nml_params, init_m4ago_params, alloc_mem_m4ago

  ! Public fields and parameters
  public :: ws_agg,&
          & aggregate_diagnostics,kav_dp,kav_rho_p,kav_d_C,kws_agg,kdf_agg,kstickiness_agg,kb_agg, &
          & kstickiness_frustule,kLmax_agg,kdynvis,kav_rhof_V,kav_por_V

  integer  :: i,j,k


  ! model parameters
  ! primary particle diameter for POM & PIM species involved in parametrized aggregation (m)
  real :: dp_dust ! primary particle diameter dust
  real :: dp_det  ! primary particle diameter detritus
  real :: dp_calc ! primary particle diameter calc
  real :: dp_opal ! primary particle diameter opal
  real :: stickiness_TEP  ! stickiness of TEP (related to opal frustules)
  real :: stickiness_det  ! normal detritus stickiness
  real :: stickiness_opal ! stickiness of opal (without TEP - just normal coating)
  real :: stickiness_calc ! stickiness of calc particles (coated with organics)
  real :: stickiness_dust ! stickiness of dust particles (coated with organics)
  real :: agg_df_max      ! maximum fractal dimension of aggregates (~2.5)
  real :: agg_df_min      ! minimum fractal dimension of aggregates (~1.2 - 1.6)
  real :: rho_TEP         ! density of TEP particles
  real :: agg_org_dens    ! organic detritus density (alternative to orgdens to avoid negative ws)

  real :: agg_Re_crit ! critical particle Reynolds number for nr-distribution limiting
  real :: POM_remin_q10 ! Q10 factor for organic remineralization (POC)
  real :: POM_remin_Tref
  real :: opal_remin_q10 ! Q10 factor for silicate remineralization (OPAL)
  real :: opal_remin_Tref

  real,allocatable ::  av_dp(:,:,:),              &  ! mean primary particle diameter
                    &  av_rho_p(:,:,:),           &  ! mean primary particle density
                    &  df_agg(:,:,:),             &  ! fractal dimension of aggregates
                    &  b_agg(:,:,:),              &  ! aggregate number distribution slope
                    &  Lmax_agg(:,:,:),           &  ! maximum diameter of aggregates
                    &  ws_agg(:,:,:),             &  ! aggregate mean sinking velocity
                    &  stickiness_agg(:,:,:),     &  ! mean aggregate stickiness
                    &  stickiness_frustule(:,:,:),&  ! frustule stickiness
                    &  N_agg(:,:,:),              &  ! Number of aggregates
                    &  av_d_C(:,:,:),             &  ! concentration-weighted mean diameter of aggs
                    &  dyn_vis(:,:,:),            &  ! molecular dynamic viscosity
                    &  m4ago_ppo(:,:,:)              ! pressure

  integer, parameter :: &
                       kav_dp               =  1, &
                       kav_rho_p            =  2, &
                       kav_d_C              =  3, &
                       kws_agg              =  4, &
                       kdf_agg              =  5, &
                       kstickiness_agg      =  6, &
                       kb_agg               =  7, &
                       kstickiness_frustule =  8, &
                       kLmax_agg            =  9, &
                       kdynvis              = 10, &
                       kav_rhof_V           = 11, &
                       kav_por_V            = 12, &
                       naggdiag             = 12

  real, dimension (:,:,:,:), allocatable, target :: aggregate_diagnostics    ! 3d concentration EU



  ! Internally used parameters and values
  real, parameter :: ONE_SIXTH = 1./6.
  real, parameter :: PI = 3.141592654
  real, parameter :: NUM_FAC = 1.e9             ! factor to avoid numerical precision problems
  real, parameter :: EPS_ONE = EPSILON(1.)

  real :: det_mol2mass ! mol detritus P/m^3 to kg POM /m^3 (according to stoichiometry)
  real :: AJ1, AJ2, AJ3, BJ1, BJ2, BJ3 ! constants for CD
  real :: grav_acc_const  ! gravitational acceleration constant
  real :: rho_aq          ! water reference density  (1025 kg/m^3)
  real :: n_det,n_opal,n_calc,n_dust               ! total primary particle number (#)
  real :: mf                                       ! mass factor for aggregates
  real :: V_dp_dust,V_dp_det,V_dp_calc,V_dp_opal   ! volumes of primary particles (L^3)
  real :: A_dp_dust,A_dp_det,A_dp_calc,A_dp_opal   ! surface areas of primary particles (L^2)
  real :: A_dust,A_det,A_calc,A_opal,A_total       ! total surface area of primary particles per unit volume (L^2/L^3)
  real :: stickiness_min, stickiness_max           ! minimum and maximum stickiness of primary particles
  real :: stickiness_mapped                        ! mapped mean stickiness of particles on range (0,1)
  real :: df_slope                                 ! slope for stickiness to fractal dimension mapping
  real :: rho_V_dp_dust,rho_V_dp_det,rho_V_dp_calc ! rho_V_dp_opal ! mass of primary particles (M)
  real :: V_det,V_opal,V_calc,V_dust,V_solid       ! total volume of primary particles in a unit volume (L^3/L^3)
  real :: Rm_SiP                                   ! molar mass ratio opal (SiO_2) to POM
  real :: thick_shell                              ! diatom frustule shell thickness (L)
  real :: d_frustule_inner                         ! diameter of hollow part in diatom frustule (L)
  real :: V_frustule_inner                         ! volume of hollow part in diatom frustule (L^3)
  real :: V_frustule_opal                          ! volume of opal shell material (L^3)
  real :: rho_V_frustule_opal                      ! mass of frustule material (M)
  real :: cell_det_mass                            ! mass of detritus material in diatoms
  real :: cell_pot_det_mass                        ! potential (max) mass detritus material in diatoms
  real :: free_detritus                            ! freely available detritus mass outside the frustule
  real :: V_POM_cell                               ! volume of POM in frustule
  real :: V_aq                                     ! volume of water space in frustule
  real :: rho_frustule                             ! density of diatom frustule incl. opal, detritus and water
  real :: rho_diatom                               ! density of either hollow frustule

contains

  !===================================================================================== m4ago_init_params
  subroutine init_m4ago_nml_params
    !>
    !! Initialization of namelist parameters
    !!
    implicit none
    ! Primary particle sizes
    dp_dust = 2.e-6      ! following the classical HAMOCC parametrization
    dp_det  = 4.e-6   ! not well defined
    dp_calc = 3.e-6   ! following Henderiks 2008, Henderiks & Pagani 2008
    dp_opal = 20.e-6  ! rough guestimate - literature search required

    ! Stickiness values
    stickiness_TEP    = 0.19
    stickiness_det    = 0.1
    stickiness_opal   = 0.08
    stickiness_calc   = 0.09
    stickiness_dust   = 0.07

    ! minimum and maximum aggregate fractal dimension
    agg_df_min        = 1.6
    agg_df_max        = 2.4

    ! Density of primary particles
    rho_TEP           = 800. ! 700.-840. kg/m^3 Azetsu-Scott & Passow 2004
    agg_org_dens      = 1100. ! detritus density - don't use orgdens to avoid negative ws

    agg_Re_crit       = 20.  ! critical particle Reynolds number for limiting nr-distribution

  end subroutine init_m4ago_nml_params

  subroutine init_m4ago_params
    !>
    !! Initilization of parameters
    !!

    implicit none
    det_mol2mass   = 3166.  ! mol detritus P/m^3 to kg POM /m^3 (according to stoichiometry)
    grav_acc_const = 9.81   ! gravitational acceleration constant
    rho_aq         = 1025.  ! water reference density  (1025 kg/m^3)

    ! CD parameters (formula 16)
    AJ1 = 24.00
    AJ2 = 29.03
    AJ3 = 14.15
    BJ1 = 1.0
    BJ2 = 0.871
    BJ3 = 0.547

    V_dp_dust = ONE_SIXTH*PI*dp_dust**3.*NUM_FAC
    V_dp_det  = ONE_SIXTH*PI*dp_det**3.*NUM_FAC
    V_dp_calc = ONE_SIXTH*PI*dp_calc**3.*NUM_FAC
    V_dp_opal = ONE_SIXTH*PI*dp_opal**3.*NUM_FAC
    A_dp_dust = PI*dp_dust**2.*NUM_FAC
    A_dp_det  = PI*dp_det**2.*NUM_FAC
    A_dp_calc = PI*dp_calc**2.*NUM_FAC
    A_dp_opal = PI*dp_opal**2.*NUM_FAC

    rho_V_dp_dust = V_dp_dust*claydens
    rho_V_dp_det  = V_dp_det*agg_org_dens
    rho_V_dp_calc = V_dp_calc*calcdens

    Rm_SiP = ropal*opalwei/det_mol2mass
    ! shell thickness
    thick_shell = 0.5*dp_opal*(1. - (opaldens/(Rm_SiP*agg_org_dens+opaldens))**(1./3.))
    d_frustule_inner = dp_opal - 2.*thick_shell
    ! volume of hollow part of frustule
    V_frustule_inner = ONE_SIXTH* PI*d_frustule_inner**3.*NUM_FAC
    ! volume of opal part of frustule
    V_frustule_opal = ONE_SIXTH*PI*(dp_opal**3. - d_frustule_inner**3.)*NUM_FAC
    rho_V_frustule_opal = V_frustule_opal*opaldens

    stickiness_min = min(stickiness_TEP, stickiness_det, stickiness_opal, stickiness_calc, stickiness_dust)
    stickiness_max = max(stickiness_TEP, stickiness_det, stickiness_opal, stickiness_calc, stickiness_dust)
    df_slope = log(agg_df_min / agg_df_max)
  end subroutine init_m4ago_params


  subroutine alloc_mem_m4ago(kpie, kpje, kpke)
    !-----------------------------------------------------------------------
    !>
    !! Initialization/allocation fields
    !! Called in ini_bgc after read_namelist
    !!

    implicit none

    integer, intent(in)  :: kpie                  !< 1st real of model grid.
    integer, intent(in)  :: kpje                  !< 2nd real of model grid.
    integer, intent(in)  :: kpke                  !< 3rd (vertical) real of model grid.

    ! allocate memory space for aggregate properties
    allocate(av_dp(kpie,kpje,kpke))
    allocate(av_rho_p(kpie,kpje,kpke))
    allocate(df_agg(kpie,kpje,kpke))
    allocate(b_agg(kpie,kpje,kpke))
    allocate(Lmax_agg(kpie,kpje,kpke))
    allocate(av_d_C(kpie,kpje,kpke))
    allocate(stickiness_agg(kpie,kpje,kpke))
    allocate(stickiness_frustule(kpie,kpje,kpke))
    allocate(aggregate_diagnostics(kpie, kpje, kpke, naggdiag))

    ! mean sinking velocity
    allocate(ws_agg(kpie,kpje,kpke))

    ! molecular dynamic viscosity
    allocate(dyn_vis(kpie, kpje, kpke))
    allocate(m4ago_ppo(kpie,kpje,kpke))

    av_dp    = 0.
    av_rho_p = 0.
    df_agg   = 0.
    b_agg    = 0.
    Lmax_agg = 0.
    av_d_C   = 0.
    stickiness_agg = 0.
    stickiness_frustule = 0.
    aggregate_diagnostics = 0.
    m4ago_ppo = 0.

  end subroutine alloc_mem_m4ago

  subroutine cleanup_mem_m4ago
    deallocate(av_dp)
    deallocate(av_rho_p)
    deallocate(df_agg)
    deallocate(b_agg)
    deallocate(Lmax_agg)
    deallocate(av_d_C)
    deallocate(stickiness_agg)
    deallocate(stickiness_frustule)
    deallocate(aggregate_diagnostics)
    deallocate(ws_agg)
    deallocate(dyn_vis)
    deallocate(m4ago_ppo)
  end subroutine cleanup_mem_m4ago

  !===================================================================================== pressure
  subroutine calc_pressure(kpie, kpje, kpke,kbnd, pddpo,omask)

    use mo_vgrid, only: ptiestu

    implicit none

    integer, intent(in)  :: kpie                  !< 1st real of model grid.
    integer, intent(in)  :: kpje                  !< 2nd real of model grid.
    integer, intent(in)  :: kpke                  !< 3rd (vertical) real of model grid.
    integer, intent(in)  :: kbnd
    real, intent(in) :: pddpo(kpie,kpje,kpke)     !< size of scalar grid cell (3rd dimension) [m]
    real, intent(in) :: omask(kpie,kpje)          !< mask

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = 1,kpke
      do j = 1,kpje
        do i = 1,kpie
          if(omask(i,j) > 0.5 .and. pddpo(i,j,k).gt.dp_min) then
            m4ago_ppo(i,j,k) = 1e5 * ptiestu(i,j,k)*98060.*1.027e-6 ! pressure in unit Pa, 98060 = onem
          endif
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine calc_pressure

  !===================================================================================== mean_agg_ws
  subroutine mean_aggregate_sinking_speed(kpie, kpje, kpke, kbnd, pddpo, omask, ptho, psao, ppao, prho)
    !-----------------------------------------------------------------------
    !>
    !! calculates the mass concentration-weighted mean sinking velocity of marine
    !! aggregates
    !!

    implicit none

    integer, intent(in)  :: kpie                  !< 1st real of model grid.
    integer, intent(in)  :: kpje                  !< 2nd real of model grid.
    integer, intent(in)  :: kpke                  !< 3rd (vertical) real of model grid.
    integer, intent(in)  :: kbnd
    real, intent(in) :: pddpo(kpie,kpje,kpke) !< size of scalar grid cell (3rd dimension) [m]
    real, intent(in) :: omask(kpie,kpje)
    real, intent(in) :: ptho (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) !< potential temperature [deg C]
    real, intent(in) :: psao (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) !< salinity [psu.].
    real, intent(in) :: ppao (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd) !< pressure at sea level [Pa].
    real, intent(in) :: prho (kpie,kpje,kpke) !< density [g/cm3]

    call calc_pressure(kpie, kpje, kpke,kbnd, pddpo, omask)

    ! molecular dynamic viscosity
    call dynvis(kpie, kpje, kpke, kbnd, pddpo, omask, ptho, psao, m4ago_ppo)

    ! ======== calculate the mean sinking velocity of aggregates =======
    call aggregate_properties(kpie, kpje, kpke, kbnd, pddpo, omask, ptho)
    call ws_Re_approx(kpie, kpje, kpke, pddpo, omask)

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do j = 1,kpje
      do i = 1,kpie
        do k = 1,kpke
          if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then
            ! Limit settling velocity wrt CFL:
            ws_agg(i,j,k) = min(ws_agg(i,j,k), 0.99*pddpo(i,j,k))

            ! ============================== Write general diagnostics ============
            ! ----- settling velocity-related -----
            aggregate_diagnostics(i,j,k,kws_agg) = ws_agg(i,j,k)/dtb  ! applied ws conversion  m/time_step  to  m/d for output

            ! ----- settling environment -----
            aggregate_diagnostics(i,j,k,kdynvis) = dyn_vis(i,j,k)     ! dynamic viscosity

            ! ----- aggregate properties -----
            av_d_C(i,j,k) = (1. + df_agg(i,j,k) - b_agg(i,j,k))                      &
                          & /(2. + df_agg(i,j,k) - b_agg(i,j,k))                     &
                          & *(Lmax_agg(i,j,k)**(2. + df_agg(i,j,k) - b_agg(i,j,k))   &
                          & - av_dp(i,j,k)**(2. + df_agg(i,j,k) - b_agg(i,j,k)))     &
                          & / (Lmax_agg(i,j,k)**(1.+df_agg(i,j,k)-b_agg(i,j,k))      &
                          & - av_dp(i,j,k)**(1. + df_agg(i,j,k)-b_agg(i,j,k)))

            aggregate_diagnostics(i,j,k,kstickiness_agg) = stickiness_agg(i,j,k)           ! aggre. stickiness
            aggregate_diagnostics(i,j,k,kstickiness_frustule) = stickiness_frustule(i,j,k) ! frustule stickiness

            aggregate_diagnostics(i,j,k,kLmax_agg)  = Lmax_agg(i,j,k)   ! applied max. diameter
            aggregate_diagnostics(i,j,k,kav_dp)     = av_dp(i,j,k)      ! mean primary particle diameter
            aggregate_diagnostics(i,j,k,kav_rho_p)  = av_rho_p(i,j,k)   ! mean primary particle density
            aggregate_diagnostics(i,j,k,kav_d_C)    = av_d_C(i,j,k)     ! conc-weighted mean agg. diameter
            aggregate_diagnostics(i,j,k,kdf_agg)    = df_agg(i,j,k)     ! aggregate fractal dim
            aggregate_diagnostics(i,j,k,kb_agg)     = b_agg(i,j,k)      ! aggre number distr. slope

            ! volume-weighted aggregate density
            aggregate_diagnostics(i,j,k,kav_rhof_V) = (av_rho_p(i,j,k)-rho_aq)*av_dp(i,j,k)**(3.-df_agg(i,j,k)) &
                        & *(4.-b_agg(i,j,k))*(Lmax_agg(i,j,k)**(1.+df_agg(i,j,k)-b_agg(i,j,k))                  &
                        &         - av_dp(i,j,k)**(1.+df_agg(i,j,k)-b_agg(i,j,k)))                              &
                        &  / ((1.+df_agg(i,j,k)-b_agg(i,j,k))                                                   &
                        & *(Lmax_agg(i,j,k)**(4.-b_agg(i,j,k)) - av_dp(i,j,k)**(4.-b_agg(i,j,k)))) + rho_aq

            ! volume-weighted aggregate porosity
            aggregate_diagnostics(i,j,k,kav_por_V)  =  1. - ((4.-b_agg(i,j,k))                                  &
                        & *av_dp(i,j,k)**(3.-df_agg(i,j,k))                                                     &
                        & *(Lmax_agg(i,j,k)**(1.+df_agg(i,j,k)-b_agg(i,j,k))                                    &
                        & - av_dp(i,j,k)**(1.+df_agg(i,j,k)-b_agg(i,j,k))))                                     &
                        & / ((1.+df_agg(i,j,k)-b_agg(i,j,k))                                                    &
                        & *(Lmax_agg(i,j,k)**(4.-b_agg(i,j,k)) - av_dp(i,j,k)**(4.-b_agg(i,j,k))))
          endif
        enddo
      enddo
    enddo
  end subroutine mean_aggregate_sinking_speed

  !===================================================================================== aggregate_properties
  subroutine aggregate_properties(kpie, kpje, kpke, kbnd, pddpo, omask, ptho)
    !-----------------------------------------------------------------------
    !>
    !! aggregate_properties calculates
    !!   - mean stickiness/aggrega
    !!   - fractal dimension
    !!   - slope of aggregate spectrum
    !!   - mean primary particle diameter
    !!   - mean primary particle density
    !!   - maximum aggregate diameter
    !!

    implicit none

    integer, intent(in)  :: kpie                  !< 1st real of model grid.
    integer, intent(in)  :: kpje                  !< 2nd real of model grid.
    integer, intent(in)  :: kpke                  !< 3rd (vertical) real of model grid.
    integer, intent(in)  :: kbnd
    real, intent(in) :: pddpo(kpie,kpje,kpke) !< size of scalar grid cell (3rd dimension) [m]
    real, intent(in) :: omask(kpie,kpje)
    real, intent(in) :: ptho (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) !< potential temperature [deg C]

    real :: C_det,C_opal,C_calc,C_dust
    !$OMP PARALLEL DO PRIVATE(i,j,k,C_det,C_opal,C_calc,C_dust,n_det,n_opal,n_dust,n_calc,mf,V_det,&
    !$OMP                     V_opal,V_calc,V_dust,V_solid,free_detritus,rho_diatom,cell_det_mass, &
    !$OMP                     cell_pot_det_mass,V_POM_cell,V_aq,rho_frustule,A_det,A_opal,         &
    !$OMP                     A_calc,A_dust,A_total,stickiness_mapped)
    do j = 1,kpje
      do i = 1,kpie
        do k = 1,kpke
          if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then
            C_det  = 0.
            C_opal = 0.
            C_calc = 0.
            C_dust = 0.

            C_det  = abs(ocetra(i,j,k,idet))
            C_opal = abs(ocetra(i,j,k,iopal))
            C_calc = abs(ocetra(i,j,k,icalc))
            C_dust = abs(ocetra(i,j,k,ifdust))

            n_det   = 0. ! number of primary particles
            n_opal  = 0.
            n_dust  = 0.
            n_calc  = 0.
            mf      = 0.

            V_det   = 0. ! total volume of primary particles in a unit volume
            V_opal  = 0.
            V_calc  = 0.
            V_dust  = 0.
            V_solid = 0.

            free_detritus = 0.
            rho_diatom    = 0.
            ! n_det are detritus primary particle that are
            ! NOT linked to any diatom frustule
            ! n_opal are number of frustule-like primary particles possessing
            ! a density i) different from pure opal ii) due to a mixture of
            ! opal frustule, detritus inside the frustule and potentially water
            ! inside the frustule

            ! describing diatom frustule as hollow sphere
            ! that is completely or partially filled with detritus
            ! and water
            cell_det_mass     = 0.
            cell_pot_det_mass = 0.
            V_POM_cell        = 0.
            V_aq              = 0.
            rho_frustule      = 0.

            ! number of opal frustules (/NUM_FAC)
            n_opal = C_opal*opalwei/rho_V_frustule_opal
            ! maximum mass of detritus inside a frustule
            cell_pot_det_mass = n_opal*V_frustule_inner*agg_org_dens

            ! detritus mass inside frustules
            cell_det_mass = min(cell_pot_det_mass, C_det*det_mol2mass - EPS_ONE)

            ! volume of detritus component in cell
            V_POM_cell = (cell_det_mass/n_opal)/agg_org_dens

            ! if not detritus is available, water is added
            V_aq = V_frustule_inner -  V_POM_cell

            ! density of the diatom frsutules incl. opal, detritus and water
            rho_frustule = (rho_V_frustule_opal + cell_det_mass/n_opal + V_aq*rho_aq)/V_dp_opal

            ! mass of extra cellular detritus particles
            free_detritus = C_det*det_mol2mass  - cell_det_mass
            rho_diatom = (rho_frustule + cell_det_mass/cell_pot_det_mass*rho_TEP)                  &
                           /(1. + cell_det_mass/cell_pot_det_mass)

            ! number of primary particles
            n_det  = free_detritus/rho_V_dp_det  ! includes NUM_FAC
            n_calc = C_calc*calcwei/rho_V_dp_calc
            n_dust = C_dust/rho_V_dp_dust     ! dust is in kg/m3

            ! primary particles surface weighted stickiness is mapped
            ! on range between 0 and 1
            ! fractal dimension of aggregates is based on that mapped df
            ! number distribution slope b is based on df

            ! calc total areas
            A_det   = n_det*A_dp_det
            A_opal  = n_opal*A_dp_opal
            A_calc  = n_calc*A_dp_calc
            A_dust  = n_dust*A_dp_dust
            A_total = A_det + A_opal + A_calc + A_dust

            ! calc frustule stickiness
            stickiness_frustule(i,j,k) = cell_det_mass/(cell_pot_det_mass +EPS_ONE)*stickiness_TEP &
                                       & + (1. - cell_det_mass/(cell_pot_det_mass + EPS_ONE))      &
                                       &   *stickiness_opal

            ! calc mean stickiness
            stickiness_agg(i,j,k) = stickiness_frustule(i,j,k)*A_opal                              &
                                  & + stickiness_det*A_det                                         &
                                  & + stickiness_calc*A_calc                                       &
                                  & + stickiness_dust*A_dust

            stickiness_agg(i,j,k) = stickiness_agg(i,j,k)/(A_total+EPS_ONE)

            stickiness_mapped = (stickiness_agg(i,j,k) - stickiness_min)                           &
                               & /(stickiness_max - stickiness_min)

            df_agg(i,j,k) = agg_df_max*exp(df_slope*stickiness_mapped)

            ! Slope is here positive defined (as n(d)~d^-b), so *-1 of
            ! Jiang & Logan 1991: Fractal dimensions of aggregates
            ! determined from steady-state size distributions.
            ! Environ. Sci. Technol. 25, 2031-2038.
            !
            ! See also:
            ! Hunt 1980: Prediction of oceanic particle size distributions
            !            from coagulation and sedimentation mechanisms.
            !
            ! Additional assumptions made here:
            ! b in Jiang & Logan     (used for       Re <   0.1: b=1
            !                              for 0.1 < Re <  10  : b=0.871
            !                              for 10  < Re < 100  : b=0.547)
            ! is set to 0.871 as an 'average for our range of 0<Re<Re_crit'
            ! D2=min(2,df(3d)) (Meakin 1988)
            !
            ! => Formulation in Jiang & Logan 1991:
            ! slope = -0.5*(3+df+(2+df-D2)/(2-b)) reduces to:

            b_agg(i,j,k) = 0.5*(3. + df_agg(i,j,k)                                                 &
                        & + (2. + df_agg(i,j,k) - min(2., df_agg(i,j,k)))/(2. - BJ2))
 
            ! careful: for df=1.5904: b_agg=2*df where w_s is undefined.

            ! total volume of primary particles
            V_det   = n_det*V_dp_det*NUM_FAC
            V_opal  = n_opal*V_dp_opal*NUM_FAC
            V_calc  = n_calc*V_dp_calc*NUM_FAC
            V_dust  = n_dust*V_dp_dust*NUM_FAC
            V_solid = V_det + V_opal + V_calc + V_dust
 
            ! primary particle mean diameter according to Bushell & Amal 1998, 2000
            ! sum(n_i) not changing - can be pulled out and thus cancels out
            av_dp(i,j,k) = (n_calc*dp_calc**3. + n_dust*dp_dust**3. + n_opal*dp_opal**3.           &
                         &  + n_det*dp_det**3.)
            av_dp(i,j,k) = av_dp(i,j,k)/(n_calc*dp_calc**df_agg(i,j,k)                             &
                         & + n_dust*dp_dust**df_agg(i,j,k) &
                         & + n_opal*dp_opal**df_agg(i,j,k) + n_det*dp_det**df_agg(i,j,k))
            av_dp(i,j,k) = av_dp(i,j,k)**(1./(3. - df_agg(i,j,k)))

            ! density of mean primary particles
            av_rho_p(i,j,k) = (V_det*agg_org_dens + V_opal*rho_diatom + V_calc*calcdens            &
                            &  + V_dust*claydens)/V_solid
          endif
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    ! calculate the maximum diameter of aggregates based on agg props
    call max_agg_diam(kpie, kpje, kpke, pddpo, omask)

  end subroutine aggregate_properties


  !================================== Reynolds number based on diameter
  real function Re_fun(ws,d,mu,rho)
    !-----------------------------------------------------------------------
    !>
    !! Reynolds number for settling particles
    !!

    implicit none

    real,intent(in) :: ws,d,mu,rho

    Re_fun = abs(ws*d*rho/mu)

  end function Re_fun


  !==================================================================================================
  !===================================================================================== ws_Re_approx
  subroutine ws_Re_approx(kpie, kpje, kpke, pddpo, omask)
    !-----------------------------------------------------------------------
    !>
    !! ws_Re_approx:  distribution integrated to Lmax (Re crit dependent maximum agg size)
    !! Renolds number-dependent sinking velocity.
    !! Approximation for c_D-value taken from Jiang & Logan 1991:
    !! c_D=a*Re^-b
    !!

    implicit none

    integer, intent(in)  :: kpie                  !< 1st real of model grid.
    integer, intent(in)  :: kpje                  !< 2nd real of model grid.
    integer, intent(in)  :: kpke                  !< 3rd (vertical) real of model grid.
    real, intent(in) :: pddpo(kpie,kpje,kpke) !< size of scalar grid cell (3rd dimension) [m]
    real, intent(in) :: omask(kpie,kpje)

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do j = 1,kpje
      do i = 1,kpie
        do k = 1,kpke
          if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then
            ws_agg(i,j,k) = ws_Re(i,j,k,Lmax_agg(i,j,k))
          endif
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine ws_Re_approx

  real function get_dRe(i, j, k, AJ, BJ, Re)
    implicit none
    ! Arguments
    integer, intent(in)  :: i                  !< 1st real of model grid.
    integer, intent(in)  :: j                  !< 2nd real of model grid.
    integer, intent(in)  :: k                  !< 3rd (vertical) real of model grid.
    real, intent(in) :: AJ
    real, intent(in) :: BJ
    real, intent(in) :: Re

    ! Local variables

    real :: nu_vis

    nu_vis =  dyn_vis(i,j,k)/rho_aq

    get_dRe = (Re*nu_vis)**((2. - BJ)/df_agg(i,j,k))/(4./3.*(av_rho_p(i,j,k) - rho_aq)/rho_aq &
           *av_dp(i,j,k)**(3. - df_agg(i,j,k))*grav_acc_const/(AJ*nu_vis**(BJ)))**(1./df_agg(i,j,k))

  end function get_dRe

  real function get_ws_agg_integral(i, j, k, AJ, BJ, lower_bound, upper_bound)
    implicit none

    integer, intent(in)  :: i                  !< 1st real of model grid.
    integer, intent(in)  :: j                  !< 2nd real of model grid.
    integer, intent(in)  :: k                  !< 3rd (vertical) real of model grid.

    real, intent(in) :: AJ
    real, intent(in) :: BJ
    real, intent(in) :: upper_bound
    real, intent(in) :: lower_bound

    ! Local variables
    real :: nu_vis

    nu_vis =  dyn_vis(i,j,k)/rho_aq
    get_ws_agg_integral = (4./3.*(av_rho_p(i,j,k) - rho_aq)/rho_aq                                 &
                     & *av_dp(i,j,k)**(3. - df_agg(i,j,k))*grav_acc_const                          &
                     & /(AJ*nu_vis**BJ))**(1./(2. - BJ))                                           &
                     & *(upper_bound**(1. - b_agg(i,j,k) + df_agg(i,j,k)                           &
                     & + (BJ + df_agg(i,j,k) - 2.)/(2. - BJ))                                      &
                     & /(1. - b_agg(i,j,k) + df_agg(i,j,k) + (BJ + df_agg(i,j,k) - 2.)/(2. - BJ))  &
                     & - lower_bound**(1. - b_agg(i,j,k) + df_agg(i,j,k) + (BJ + df_agg(i,j,k) -2.)&
                     & /(2. - BJ))                                                                 &
                     & /(1. - b_agg(i,j,k) + df_agg(i,j,k) + (BJ + df_agg(i,j,k) - 2.)/(2. - BJ)))

  end function get_ws_agg_integral

  !===================================================================================== ws_Re
  real function ws_Re(i,j,k,dmax_agg)
    !-----------------------------------------------------------------------
    !>
    !! ws_Re:  distribution integrated to Lmax (Re crit dependent maximum agg size)
    !! Reynolds number-dependent sinking velocity.
    !! Approximation for c_D-value taken from Jiang & Logan 1991:
    !! c_D=a*Re^-b
    !! written in such a way that we check the critical Reynolds
    !! number (in case that we extend the maximum size by shear-
    !! driven break-up).
    !!

    implicit none

    integer, intent(in)  :: i                  !< 1st real of model grid.
    integer, intent(in)  :: j                  !< 2nd real of model grid.
    integer, intent(in)  :: k                  !< 3rd (vertical) real of model grid.
    real, intent(in) :: dmax_agg

    ! Local
    real :: d_Re01, d_Re10, d_low, ws_agg_ints

    ! for Re-dependent, it should always be agg_Re_crit>10
    ! for shear-driven break-up, check against integration bounds
    ! calc integration limits for Re-dependent sinking:
    ! Re=0.1
    d_Re01 = get_dRe(i,j,k, AJ1, BJ1, 0.1)
    ! Re=10
    d_Re10 = get_dRe(i,j,k, AJ2, BJ2, 10.)
    d_low = av_dp(i,j,k)

    ws_agg_ints = 0.
    if(dmax_agg >= d_Re01)then ! Re > 0.1
                                       ! - collect full range up to
                                       ! 0.1, (dp->d_Re1) and set lower bound to
                                       ! Re=0.1 val
                                       ! aj=AJ1, bj=1
        ws_agg_ints = get_ws_agg_integral(i, j, k, AJ1, BJ1, av_dp(i,j,k), d_Re01)
        d_low = d_Re01
    endif

    if(dmax_agg >= d_Re10)then ! Re > 10
                                         ! - collect full range Re=0.1-10 (d_Re1-> d_Re2)
                                         ! and set lower bound to
                                         ! Re=10 val
                                         ! aj=AJ2, bj=0.871
        ws_agg_ints = ws_agg_ints  + get_ws_agg_integral(i, j, k, AJ2, BJ2, d_Re01, d_Re10)
        d_low = d_Re10
    endif

    if(d_low < d_Re01)then ! Re<0.1 and Lmax < d_Re1
        ws_agg_ints = get_ws_agg_integral(i, j, k, AJ1, BJ1, av_dp(i,j,k), dmax_agg)
    else ! Re > 10, aj=AJ3, bj=BJ3
        ws_agg_ints = ws_agg_ints + get_ws_agg_integral(i, j, k, AJ3, BJ3, d_low, dmax_agg)
    endif

    ! concentration-weighted mean sinking velocity
    ws_Re = (ws_agg_ints &
            & /((dmax_agg**(1. + df_agg(i,j,k) - b_agg(i,j,k))  &
            & - av_dp(i,j,k)**(1. + df_agg(i,j,k) - b_agg(i,j,k)))  &
            & / (1. + df_agg(i,j,k) - b_agg(i,j,k))))*dtbgc   ! (m/s -> m/d)  *dtb

  end function ws_Re


  subroutine max_agg_diam(kpie, kpje, kpke, pddpo, omask)
    !-----------------------------------------------------------------------
    !>
    !! max_agg_diam calculates the maximum aggregate diameter of the aggregate
    !! number distribution, assumes Re_crit > 10
    !!
    integer, intent(in)  :: kpie                  !< 1st real of model grid.
    integer, intent(in)  :: kpje                  !< 2nd real of model grid.
    integer, intent(in)  :: kpke                  !< 3rd (vertical) real of model grid.
    real, intent(in) :: pddpo(kpie,kpje,kpke) !< size of scalar grid cell (3rd dimension) [m]
    real, intent(in) :: omask(kpie,kpje)

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    ! base on analytical Jiang approximation
    do j = 1,kpje
      do i = 1,kpie
        do k = 1,kpke
          if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then
            Lmax_agg(i,j,k)   = max_agg_diam_white(i,j,k)
          endif
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine max_agg_diam

  !================================================ maximum diameter of agg in non-stratified fluid
  real  function max_agg_diam_white(i,j,k)
    !-------------------------------------------------------------------------
    !>
    !! maximum aggregate diameter in a non-stratified fluid - following the
    !! White drag approaximation by Jiang & Logan 1991, assuming agg_re_crit > 10
    !! (otherwise AJX,BJX needs to be adjusted)
    !!

    implicit none

    integer,intent(in) :: i,j,k
    real        :: nu_vis

    nu_vis  =  dyn_vis(i,j,k)/rho_aq
    max_agg_diam_white = (agg_Re_crit*nu_vis)**((2. - BJ3)/df_agg(i,j,k))                          &
                        & /((4./3.)*(av_rho_p(i,j,k) - rho_aq)/rho_aq                              &
                        & *av_dp(i,j,k)**(3. - df_agg(i,j,k))*grav_acc_const                       &
                        & /(AJ3*nu_vis**BJ3))**(1./df_agg(i,j,k))

  end function max_agg_diam_white

  !===================================================================================== mass factor
  real  function mass_factor(dp,df,rhop)
    !-----------------------------------------------------------------------
    !>
    !! mass_factor calculates the mass factor for the mass of a single
    !! aggregate
    !!
    implicit none

    real, intent(in) :: dp
    real, intent(in) :: df
    real, intent(in) :: rhop

    ! mass factor
    mass_factor = ONE_SIXTH * PI * dp**(3. - df) * rhop

  end function mass_factor


  !===================================================================================== rho_agg
  real function rho_agg(d,rhop,dp,df,rho)
    !-----------------------------------------------------------------------
    !>
    !! rho_agg provides the aggregate density
    !!

    implicit none

    real, intent(in) :: d
    real, intent(in) :: rhop
    real, intent(in) :: dp
    real, intent(in) :: df
    real, intent(in) :: rho

    rho_agg =  (rhop - rho)*(dp/d)**(3. - df) + rho

  end function rho_agg

  !===================================================================================== dynvis
  subroutine dynvis(kpie, kpje, kpke, kbnd, pddpo, omask, ptho, psao, ppo)
    !-----------------------------------------------------------------------
    !>
    !! dynvis calculates the molecular dynamic viscosity according to
    !! Richards 1998: The effect of temperature, pressure, and salinity
    !! on sound attenuation in turbid seawater. J. Acoust. Soc. Am. 103 (1),
    !! originally published by  Matthaeus, W. (1972): Die Viskositaet des
    !! Meerwassers. Beitraege zur Meereskunde, Heft 29 (in German).
    !!

    implicit none

    integer, intent(in)  :: kpie                  !< 1st real of model grid.
    integer, intent(in)  :: kpje                  !< 2nd real of model grid.
    integer, intent(in)  :: kpke                  !< 3rd (vertical) real of model grid.
    integer, intent(in)  :: kbnd

    real, intent(in) :: pddpo(kpie,kpje,kpke) !< size of scalar grid cell (3rd dimension) [m]
    real, intent(in) :: omask(kpie,kpje)
    real, intent(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) !< potential temperature [deg C]
    real, intent(in) :: psao(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)  !< salinity [psu.].
    real, intent(in) :: ppo(kpie,kpje,kpke)  !< pressure [Pa].

    ! Local variables
    real:: press_val  ! Pascal/rho -> dbar
    real:: ptho_val,psao_val
    integer :: kch
    kch = 0
    !$OMP PARALLEL DO PRIVATE(i,j,k,press_val,ptho_val,psao_val,kch)
    do j = 1,kpje
      do i = 1,kpie
        do k = 1,kpke
          if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then
            kch = merge(k+1,k,k<kpke)
            if(pddpo(i,j,kch) > 0.5) then
              press_val    = 0.5*(ppo(i,j,k)  + ppo(i,j,kch))*1.e-5 ! Pascal -> dbar
              ptho_val     = 0.5*(ptho(i,j,k) + ptho(i,j,kch))
              psao_val     = 0.5*(psao(i,j,k) + ptho(i,j,kch))
            else
              press_val    = ppo(i,j,k)*1.e-5 ! Pascal -> dbar
              ptho_val     = ptho(i,j,k)
              psao_val     = psao(i,j,k)
            endif

            ! molecular dynamic viscosity
            dyn_vis(i,j,k) = 0.1    & ! Unit: g / (cm*s) -> kg / (m*s)
                           &     *(1.79e-2                                                         &
                           &     - 6.1299e-4*ptho_val + 1.4467e-5*ptho_val**2.                     &
                           &     - 1.6826e-7*ptho_val**3.                                          &
                           &     - 1.8266e-7*press_val  + 9.8972e-12*press_val**2.                 &
                           &     + 2.4727e-5*psao_val                                              &
                           &     + psao_val*(4.8429e-7*ptho_val - 4.7172e-8*ptho_val**2.           &
                           &     + 7.5986e-10*ptho_val**3.)                                        &
                           &     + press_val*(1.3817e-8*ptho_val - 2.6363e-10*ptho_val**2.)        &
                           &     - press_val**2.*(6.3255e-13*ptho_val - 1.2116e-14*ptho_val**2.))
          endif
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine dynvis

end module mo_m4ago

