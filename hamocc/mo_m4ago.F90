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



MODULE mo_m4ago
     USE mo_vgrid,       ONLY: dp_min
     USE mo_control_bgc, ONLY: dtb, dtbgc,io_stdo_bgc
     USE mo_sedmnt,      ONLY: calcdens, claydens, opaldens, calcwei, opalwei                  
     USE mo_carbch,      ONLY: ocetra
     USE mo_param1_bgc,  ONLY: iopal, ifdust, icalc, idet
     USE mo_biomod,      ONLY: ropal

     IMPLICIT NONE

     PRIVATE
     
     ! Public subroutines
     PUBLIC :: mean_aggregate_sinking_speed, init_m4ago_nml_params, init_m4ago_params, alloc_mem_m4ago
     
     ! Public fields and parameters
     PUBLIC :: ws_agg, POM_remin_q10, POM_remin_Tref, opal_remin_q10, opal_remin_Tref, &
             & aggregate_diagnostics,kav_dp,kav_rho_p,kav_d_C,kws_agg,kdf_agg,kstickiness_agg,kb_agg,kstickiness_frustule, &
             & kLmax_agg,kdynvis,kav_rhof_V,kav_por_V   

     INTEGER  :: i,j,k


     ! model parameters 
     ! primary particle diameter for POM & PIM species involved in parametrized aggregation (m) 
     REAL :: dp_dust ! primary particle diameter dust
     REAL :: dp_det  ! primary particle diameter detritus
     REAL :: dp_calc ! primary particle diameter calc
     REAL :: dp_opal ! primary particle diameter opal
     REAL :: stickiness_TEP  ! stickiness of TEP (related to opal frustules)
     REAL :: stickiness_det  ! normal detritus stickiness
     REAL :: stickiness_opal ! stickiness of opal (without TEP - just normal coating)
     REAL :: stickiness_calc ! stickiness of calc particles (coated with organics)
     REAL :: stickiness_dust ! stickiness of dust particles (coated with organics)
     REAL :: agg_df_max      ! maximum fractal dimension of aggregates (~2.5)
     REAL :: agg_df_min      ! minimum fractal dimension of aggregates (~1.2 - 1.6)        
     REAL :: rho_TEP         ! density of TEP particles 
     REAL :: agg_org_dens    ! organic detritus density (alternative to orgdens to avoid negative ws)

     REAL :: agg_Re_crit ! critical particle Reynolds number for nr-distribution limiting
     REAL :: POM_remin_q10 ! Q10 factor for organic remineralization (POC) 
     REAL :: POM_remin_Tref
     REAL :: opal_remin_q10 ! Q10 factor for silicate remineralization (OPAL)
     REAL :: opal_remin_Tref

     REAL,ALLOCATABLE ::  av_dp(:,:,:),              &  ! mean primary particle diameter
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
                       &  m4ago_ppo(:,:,:)              ! in situ pressure - potentially to replace by BLOM pressure

     INTEGER, PARAMETER :: &
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

     REAL, DIMENSION (:,:,:,:), ALLOCATABLE, TARGET :: aggregate_diagnostics    ! 3d concentration EU



     ! Internally used parameters and values
     REAL, PARAMETER :: ONE_SIXTH = 1./6.
     REAL, PARAMETER :: PI = 3.141592654
     REAL, PARAMETER :: NUM_FAC = 1.e9             ! factor to avoid numerical precision problems
     REAL, PARAMETER :: EPS_ONE = EPSILON(1.)
     
     REAL :: det_mol2mass ! mol detritus P/m^3 to kg POM /m^3 (according to stoichiometry)
     REAL :: AJ1, AJ2, AJ3, BJ1, BJ2, BJ3 ! constants for CD 
     REAL :: grav_acc_const  ! gravitational acceleration constant 
     REAL :: rho_aq          ! water reference density  (1025 kg/m^3)
     REAL :: n_det,n_opal,n_calc,n_dust               ! total primary particle number (#)
     REAL :: mf                                       ! mass factor for aggregates 
     REAL :: V_dp_dust,V_dp_det,V_dp_calc,V_dp_opal   ! volumes of primary particles (L^3)
     REAL :: A_dp_dust,A_dp_det,A_dp_calc,A_dp_opal   ! surface areas of primary particles (L^2)
     REAL :: A_dust,A_det,A_calc,A_opal,A_total       ! total surface area of primary particles per unit volume (L^2/L^3)
     REAL :: stickiness_min, stickiness_max           ! minimum and maximum stickiness of primary particles
     REAL :: stickiness_mapped                        ! mapped mean stickiness of particles on range (0,1)
     REAL :: df_slope                                 ! slope for stickiness to fractal dimension mapping
     REAL :: rho_V_dp_dust,rho_V_dp_det,rho_V_dp_calc ! rho_V_dp_opal ! mass of primary particles (M)
     REAL :: V_det,V_opal,V_calc,V_dust,V_solid       ! total volume of primary particles in a unit volume (L^3/L^3) 
     REAL :: Rm_SiP                                   ! molar mass ratio opal (SiO_2) to POM   
     REAL :: thick_shell                              ! diatom frustule shell thickness (L)
     REAL :: d_frustule_inner                         ! diameter of hollow part in diatom frustule (L)
     REAL :: V_frustule_inner                         ! volume of hollow part in diatom frustule (L^3)
     REAL :: V_frustule_opal                          ! volume of opal shell material (L^3)
     REAL :: rho_V_frustule_opal                      ! mass of frustule material (M) 
     REAL :: cell_det_mass                            ! mass of detritus material in diatoms
     REAL :: cell_pot_det_mass                        ! potential (max) mass detritus material in diatoms
     REAL :: free_detritus                            ! freely available detritus mass outside the frustule
     REAL :: V_POM_cell                               ! volume of POM in frustule
     REAL :: V_aq                                     ! volume of water space in frustule
     REAL :: rho_frustule                             ! density of diatom frustule incl. opal, detritus and water
     REAL :: rho_diatom                               ! density of either hollow frustule 

   CONTAINS

  !===================================================================================== m4ago_init_params
  SUBROUTINE init_m4ago_nml_params
     !>
     !! Initialization of namelist parameters
     !!
     IMPLICIT NONE
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

     POM_remin_q10     = 2.1 ! Bidle et al. 2002: Regulation of Oceanic Silicon...
     opal_remin_q10    = 2.6 ! Bidle et al. 2002: Regulation of Oceanic Silicon...
     POM_remin_Tref    = 10.
     opal_remin_Tref   = 10.
  END SUBROUTINE init_m4ago_nml_params

  SUBROUTINE init_m4ago_params
     !> 
     !! Initilization of parameters 
     !!
 
     IMPLICIT NONE
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

     stickiness_min = MIN(stickiness_TEP, stickiness_det, stickiness_opal, stickiness_calc, stickiness_dust)
     stickiness_max = MAX(stickiness_TEP, stickiness_det, stickiness_opal, stickiness_calc, stickiness_dust)
     df_slope = LOG( agg_df_min / agg_df_max)
  END SUBROUTINE init_m4ago_params


  SUBROUTINE alloc_mem_m4ago(kpie, kpje, kpke)
     !-----------------------------------------------------------------------
     !> 
     !! Initialization/allocation fields 
     !! Called in ini_bgc after read_namelist
     !!

     IMPLICIT NONE

     INTEGER, INTENT(in)  :: kpie                  !< 1st REAL of model grid.
     INTEGER, INTENT(in)  :: kpje                  !< 2nd REAL of model grid.
     INTEGER, INTENT(in)  :: kpke                  !< 3rd (vertical) REAL of model grid.

     ! allocate memory space for aggregate properties 
     ALLOCATE(av_dp(kpie,kpje,kpke))  
     ALLOCATE(av_rho_p(kpie,kpje,kpke))
     ALLOCATE(df_agg(kpie,kpje,kpke))
     ALLOCATE(b_agg(kpie,kpje,kpke))
     ALLOCATE(Lmax_agg(kpie,kpje,kpke))
     ALLOCATE(av_d_C(kpie,kpje,kpke))
     ALLOCATE(stickiness_agg(kpie,kpje,kpke))
     ALLOCATE(stickiness_frustule(kpie,kpje,kpke))
     ALLOCATE(aggregate_diagnostics(kpie, kpje, kpke, naggdiag))
     
     ! mean sinking velocity
     ALLOCATE(ws_agg(kpie,kpje,kpke))

     ! molecular dynamic viscosity
     ALLOCATE(dyn_vis(kpie, kpje, kpke)) 
     ALLOCATE(m4ago_ppo(kpie,kpje,kpke))

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

  END SUBROUTINE alloc_mem_m4ago

  SUBROUTINE cleanup_mem_m4ago
  
     DEALLOCATE(av_dp)  
     DEALLOCATE(av_rho_p)
     DEALLOCATE(df_agg)
     DEALLOCATE(b_agg)
     DEALLOCATE(Lmax_agg)
     DEALLOCATE(av_d_C)
     DEALLOCATE(stickiness_agg)
     DEALLOCATE(stickiness_frustule)
     DEALLOCATE(aggregate_diagnostics)
     DEALLOCATE(ws_agg)
     DEALLOCATE(dyn_vis) 
     DEALLOCATE(m4ago_ppo)
  END SUBROUTINE cleanup_mem_m4ago

  !===================================================================================== pressure
  SUBROUTINE calc_pressure(kpie, kpje, kpke,kbnd, pddpo,omask, ppao, prho)
     IMPLICIT NONE

     INTEGER, INTENT(in)  :: kpie                  !< 1st REAL of model grid.
     INTEGER, INTENT(in)  :: kpje                  !< 2nd REAL of model grid.
     INTEGER, INTENT(in)  :: kpke                  !< 3rd (vertical) REAL of model grid.
     INTEGER, INTENT(in)  :: kbnd
     REAL, INTENT(in) :: pddpo(kpie,kpje,kpke) !< size of scalar grid cell (3rd dimension) [m]
     REAL, INTENT(in) :: omask(kpie,kpje) 
     REAL, INTENT(in) :: ppao (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd) !< pressure at sea level [Pa].
     REAL, INTENT(in) :: prho (kpie,kpje,kpke) !< salinity [psu.].

     !$OMP PARALLEL DO PRIVATE(i,j,k)
     do j = 1,kpje
     do i = 1,kpie
        if(omask(i,j) > 0.5) then
           m4ago_ppo(i,j,1) = ppao(i,j) + prho(i,j,1)*grav_acc_const*pddpo(i,j,1)   
           do k = 2,kpke
            if(pddpo(i,j,k) > dp_min) then
               m4ago_ppo(i,j,k) = m4ago_ppo(i,j,k-1) + prho(i,j,k)*grav_acc_const*pddpo(i,j,k)  
            endif
           enddo
        endif
      enddo
      enddo 
      !$OMP END PARALLEL DO
  END SUBROUTINE calc_pressure

  !===================================================================================== mean_agg_ws
  SUBROUTINE mean_aggregate_sinking_speed(kpie, kpje, kpke, kbnd, pddpo, omask, ptho, psao, ppao, prho)
     !-----------------------------------------------------------------------
     !> 
     !! calculates the mass concentration-weighted mean sinking velocity of marine
     !! aggregates 
     !!
    
     IMPLICIT NONE
     
     INTEGER, INTENT(in)  :: kpie                  !< 1st REAL of model grid.
     INTEGER, INTENT(in)  :: kpje                  !< 2nd REAL of model grid.
     INTEGER, INTENT(in)  :: kpke                  !< 3rd (vertical) REAL of model grid.
     INTEGER, INTENT(in)  :: kbnd
     REAL, INTENT(in) :: pddpo(kpie,kpje,kpke) !< size of scalar grid cell (3rd dimension) [m]
     REAL, INTENT(in) :: omask(kpie,kpje) 
     REAL, INTENT(in) :: ptho (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) !< potential temperature [deg C]
     REAL, INTENT(in) :: psao (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) !< salinity [psu.].
     REAL, INTENT(in) :: ppao (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd) !< pressure at sea level [Pa].
     REAL, INTENT(in) :: prho (kpie,kpje,kpke) !< salinity [psu.].

     CALL calc_pressure(kpie, kpje, kpke,kbnd, pddpo, omask, ppao, prho)

     ! molecular dynamic viscosity 
     CALL dynvis(kpie, kpje, kpke, kbnd, pddpo, omask, ptho, psao, m4ago_ppo)

     ! ======== calculate the mean sinking velocity of aggregates =======
     CALL aggregate_properties(kpie, kpje, kpke, kbnd, pddpo, omask, ptho) 
     CALL ws_Re_approx(kpie, kpje, kpke, pddpo, omask)

     !$OMP PARALLEL DO PRIVATE(i,j,k)
     DO j = 1,kpje
     DO i = 1,kpie
     DO k = 1,kpke
        IF(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) THEN
            ! Limit settling velocity wrt CFL:
            ws_agg(i,j,k) = MIN(ws_agg(i,j,k), 0.99*pddpo(i,j,k))

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
        END IF
     END DO
     END DO
     END DO 

  END SUBROUTINE mean_aggregate_sinking_speed

  !===================================================================================== aggregate_properties
  SUBROUTINE aggregate_properties(kpie, kpje, kpke, kbnd, pddpo, omask, ptho) 
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
  
     IMPLICIT NONE

     INTEGER, INTENT(in)  :: kpie                  !< 1st REAL of model grid.
     INTEGER, INTENT(in)  :: kpje                  !< 2nd REAL of model grid.
     INTEGER, INTENT(in)  :: kpke                  !< 3rd (vertical) REAL of model grid.
     INTEGER, INTENT(in)  :: kbnd
     REAL, INTENT(in) :: pddpo(kpie,kpje,kpke) !< size of scalar grid cell (3rd dimension) [m]
     REAL, INTENT(in) :: omask(kpie,kpje) 
     REAL, INTENT(in) :: ptho (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) !< potential temperature [deg C]

     REAL :: C_det,C_opal,C_calc,C_dust
     !$OMP PARALLEL DO PRIVATE(i,j,k,C_det,C_opal,C_calc,C_dust,n_det,n_opal,n_dust,n_calc,mf,V_det,V_opal,V_calc,V_dust,V_solid,  &
     !$OMP                    free_detritus,rho_diatom,cell_det_mass,cell_pot_det_mass,V_POM_cell,V_aq,rho_frustule,A_det,A_opal,  &
     !$OMP                    A_calc,A_dust,A_total,stickiness_mapped)
     DO j = 1,kpje
     DO i = 1,kpie
     DO k = 1,kpke
        IF(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) THEN
           C_det  = 0.
           C_opal = 0.
           C_calc = 0.
           C_dust = 0.
              
           C_det  = ABS(ocetra(i,j,k,idet))
           C_opal = ABS(ocetra(i,j,k,iopal))
           C_calc = ABS(ocetra(i,j,k,icalc))
           C_dust = ABS(ocetra(i,j,k,ifdust))

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
           cell_det_mass = MIN(cell_pot_det_mass, C_det*det_mol2mass - EPS_ONE)

           ! volume of detritus component in cell
           V_POM_cell = (cell_det_mass/n_opal)/agg_org_dens
   
           ! if not detritus is available, water is added
           V_aq = V_frustule_inner -  V_POM_cell                
                 
           ! density of the diatom frsutules incl. opal, detritus and water 
           rho_frustule = (rho_V_frustule_opal + cell_det_mass/n_opal + V_aq*rho_aq)/V_dp_opal               
                 
           ! mass of extra cellular detritus particles
           free_detritus = C_det*det_mol2mass  - cell_det_mass
           rho_diatom = (rho_frustule + cell_det_mass/cell_pot_det_mass*rho_TEP) &
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
           stickiness_frustule(i,j,k) = cell_det_mass/(cell_pot_det_mass + EPS_ONE)*stickiness_TEP &
                                      & + (1. - cell_det_mass/(cell_pot_det_mass + EPS_ONE))*stickiness_opal

           ! calc mean stickiness
           stickiness_agg(i,j,k) = stickiness_frustule(i,j,k)*A_opal  &
                                 & + stickiness_det*A_det              &
                                 & + stickiness_calc*A_calc             &
                                 & + stickiness_dust*A_dust

           stickiness_agg(i,j,k) = stickiness_agg(i,j,k)/(A_total+EPS_ONE)

           stickiness_mapped = (stickiness_agg(i,j,k) - stickiness_min) & 
                               & /(stickiness_max - stickiness_min)

           df_agg(i,j,k) = agg_df_max*EXP(df_slope*stickiness_mapped)

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
                    
           b_agg(i,j,k) = 0.5*(3. + df_agg(i,j,k) &  
                        & + (2. + df_agg(i,j,k) - MIN(2., df_agg(i,j,k)))/(2. - BJ2))
 
           ! careful: for df=1.5904: b_agg=2*df where w_s is undefined.

           ! total volume of primary particles
           V_det   = n_det*V_dp_det*NUM_FAC
           V_opal  = n_opal*V_dp_opal*NUM_FAC
           V_calc  = n_calc*V_dp_calc*NUM_FAC
           V_dust  = n_dust*V_dp_dust*NUM_FAC
           V_solid = V_det + V_opal + V_calc + V_dust
 
           ! primary particle mean diameter according to Bushell & Amal 1998, 2000
           ! sum(n_i) not changing - can be pulled out and thus cancels out
           av_dp(i,j,k) = (n_calc*dp_calc**3. + n_dust*dp_dust**3. + n_opal*dp_opal**3. + n_det*dp_det**3.)
           av_dp(i,j,k) = av_dp(i,j,k)/(n_calc*dp_calc**df_agg(i,j,k) + n_dust*dp_dust**df_agg(i,j,k) &
                        & + n_opal*dp_opal**df_agg(i,j,k) + n_det*dp_det**df_agg(i,j,k))
           av_dp(i,j,k) = av_dp(i,j,k)**(1./(3. - df_agg(i,j,k)))

           ! density of mean primary particles
           av_rho_p(i,j,k) = (V_det*agg_org_dens + V_opal*rho_diatom + V_calc*calcdens + V_dust*claydens)/V_solid  
        END IF
     END DO
     END DO
     END DO
     !$OMP END PARALLEL DO
 
     ! calculate the maximum diameter of aggregates based on agg props     
     CALL max_agg_diam(kpie, kpje, kpke, pddpo, omask) 

  END SUBROUTINE aggregate_properties
  

  !================================== Reynolds number based on diameter
  REAL FUNCTION Re_fun(ws,d,mu,rho)
     !-----------------------------------------------------------------------
     !>
     !! Reynolds number for settling particles
     !!
     
     IMPLICIT NONE
 
     REAL,INTENT(in) :: ws,d,mu,rho

     Re_fun = ABS(ws*d*rho/mu)
 
   END FUNCTION Re_fun


  !==================================================================================================
  !===================================================================================== ws_Re_approx
  SUBROUTINE ws_Re_approx(kpie, kpje, kpke, pddpo, omask) 
     !-----------------------------------------------------------------------
     !>
     !! ws_Re_approx:  distribution integrated to Lmax (Re crit dependent maximum agg size)
     !! Renolds number-dependent sinking velocity. 
     !! Approximation for c_D-value taken from Jiang & Logan 1991:
     !! c_D=a*Re^-b
     !! 
     
     IMPLICIT NONE

     INTEGER, INTENT(in)  :: kpie                  !< 1st REAL of model grid.
     INTEGER, INTENT(in)  :: kpje                  !< 2nd REAL of model grid.
     INTEGER, INTENT(in)  :: kpke                  !< 3rd (vertical) REAL of model grid.
     REAL, INTENT(in) :: pddpo(kpie,kpje,kpke) !< size of scalar grid cell (3rd dimension) [m]
     REAL, INTENT(in) :: omask(kpie,kpje) 

     !$OMP PARALLEL DO PRIVATE(i,j,k)
     DO j = 1,kpje
     DO i = 1,kpie
     DO k = 1,kpke
       IF(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) THEN
            ws_agg(i,j,k) = ws_Re(i,j,k,Lmax_agg(i,j,k))
       END IF
     END DO
     END DO
     END DO
     !$OMP END PARALLEL DO

  END SUBROUTINE ws_Re_approx

  REAL FUNCTION get_dRe(i, j, k, AJ, BJ, Re) 
     IMPLICIT NONE
     ! Arguments
     INTEGER, INTENT(in)  :: i                  !< 1st REAL of model grid.
     INTEGER, INTENT(in)  :: j                  !< 2nd REAL of model grid.
     INTEGER, INTENT(in)  :: k                  !< 3rd (vertical) REAL of model grid.
     REAL, INTENT(in) :: AJ
     REAL, INTENT(in) :: BJ
     REAL, INTENT(in) :: Re
    
    ! Local variables

    REAL :: nu_vis 
     
    nu_vis =  dyn_vis(i,j,k)/rho_aq 

    get_dRe = (Re*nu_vis)**((2. - BJ)/df_agg(i,j,k))/(4./3.*(av_rho_p(i,j,k) - rho_aq)/rho_aq &
           *av_dp(i,j,k)**(3. - df_agg(i,j,k))*grav_acc_const/(AJ*nu_vis**(BJ)))**(1./df_agg(i,j,k))

  END FUNCTION get_dRe
  
  REAL FUNCTION get_ws_agg_integral(i, j, k, AJ, BJ, lower_bound, upper_bound)
    IMPLICIT NONE
      
    INTEGER, INTENT(in)  :: i                  !< 1st REAL of model grid.
    INTEGER, INTENT(in)  :: j                  !< 2nd REAL of model grid.
    INTEGER, INTENT(in)  :: k                  !< 3rd (vertical) REAL of model grid.
    
    REAL, INTENT(in) :: AJ
    REAL, INTENT(in) :: BJ
    REAL, INTENT(in) :: upper_bound
    REAL, INTENT(in) :: lower_bound
    
    ! Local variables
    REAL :: nu_vis 
     
    nu_vis =  dyn_vis(i,j,k)/rho_aq 
    get_ws_agg_integral = (4./3.*(av_rho_p(i,j,k) - rho_aq)/rho_aq &
                     & *av_dp(i,j,k)**(3. - df_agg(i,j,k))*grav_acc_const  &
                     & /(AJ*nu_vis**BJ))**(1./(2. - BJ)) &
                     & *(upper_bound**(1. - b_agg(i,j,k) + df_agg(i,j,k)    &
                     & + (BJ + df_agg(i,j,k) - 2.)/(2. - BJ)) &
                     & /(1. - b_agg(i,j,k) + df_agg(i,j,k) + (BJ + df_agg(i,j,k) - 2.)/(2. - BJ)) &
                     & - lower_bound**(1. - b_agg(i,j,k) + df_agg(i,j,k) + (BJ + df_agg(i,j,k) - 2.)   &
                     & /(2. - BJ)) &
                     & /(1. - b_agg(i,j,k) + df_agg(i,j,k) + (BJ + df_agg(i,j,k) - 2.)/(2. - BJ)))
    
  END FUNCTION get_ws_agg_integral 
  
  !===================================================================================== ws_Re
  REAL FUNCTION ws_Re(i,j,k,dmax_agg)
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

    IMPLICIT NONE

    INTEGER, INTENT(in)  :: i                  !< 1st REAL of model grid.
    INTEGER, INTENT(in)  :: j                  !< 2nd REAL of model grid.
    INTEGER, INTENT(in)  :: k                  !< 3rd (vertical) REAL of model grid.
    REAL, INTENT(in) :: dmax_agg

    ! Local
    REAL :: d_Re01, d_Re10, d_low, ws_agg_ints

    ! for Re-dependent, it should always be agg_Re_crit>10
    ! for shear-driven break-up, check against integration bounds
    ! calc integration limits for Re-dependent sinking:
    ! Re=0.1
    d_Re01 = get_dRe(i,j,k, AJ1, BJ1, 0.1)
    ! Re=10
    d_Re10 = get_dRe(i,j,k, AJ2, BJ2, 10.)
    d_low = av_dp(i,j,k) 

    ws_agg_ints = 0.
    IF(dmax_agg >= d_Re01)THEN ! Re > 0.1  
                                       ! - collect full range up to
                                       ! 0.1, (dp->d_Re1) and set lower bound to
                                       ! Re=0.1 val
                                       ! aj=AJ1, bj=1
        ws_agg_ints = get_ws_agg_integral(i, j, k, AJ1, BJ1, av_dp(i,j,k), d_Re01)
        d_low = d_Re01
    ENDIF

    IF(dmax_agg >= d_Re10)THEN ! Re > 10
                                         ! - collect full range Re=0.1-10 (d_Re1-> d_Re2)
                                         ! and set lower bound to
                                         ! Re=10 val
                                         ! aj=AJ2, bj=0.871
        ws_agg_ints = ws_agg_ints  + get_ws_agg_integral(i, j, k, AJ2, BJ2, d_Re01, d_Re10)
        d_low = d_Re10
    ENDIF

    IF(d_low < d_Re01)THEN ! Re<0.1 and Lmax < d_Re1 
        ws_agg_ints = get_ws_agg_integral(i, j, k, AJ1, BJ1, av_dp(i,j,k), dmax_agg)
    ELSE ! Re > 10, aj=AJ3, bj=BJ3
        ws_agg_ints = ws_agg_ints + get_ws_agg_integral(i, j, k, AJ3, BJ3, d_low, dmax_agg) 
    ENDIF

    ! concentration-weighted mean sinking velocity
    ws_Re = (ws_agg_ints &
            & /((dmax_agg**(1. + df_agg(i,j,k) - b_agg(i,j,k))  &
            & - av_dp(i,j,k)**(1. + df_agg(i,j,k) - b_agg(i,j,k)))  &
            & / (1. + df_agg(i,j,k) - b_agg(i,j,k))))*dtbgc   ! (m/s -> m/d)  *dtb

  END FUNCTION ws_Re


  SUBROUTINE max_agg_diam(kpie, kpje, kpke, pddpo, omask)
     !-----------------------------------------------------------------------
     !>
     !! max_agg_diam calculates the maximum aggregate diameter of the aggregate
     !! number distribution, assumes Re_crit > 10 
     !!
     INTEGER, INTENT(in)  :: kpie                  !< 1st REAL of model grid.
     INTEGER, INTENT(in)  :: kpje                  !< 2nd REAL of model grid.
     INTEGER, INTENT(in)  :: kpke                  !< 3rd (vertical) REAL of model grid.
     REAL, INTENT(in) :: pddpo(kpie,kpje,kpke) !< size of scalar grid cell (3rd dimension) [m]
     REAL, INTENT(in) :: omask(kpie,kpje) 

     !$OMP PARALLEL DO PRIVATE(i,j,k)
     ! base on analytical Jiang approximation
     DO j = 1,kpje
     DO i = 1,kpie
     DO k = 1,kpke
        IF(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) THEN
             Lmax_agg(i,j,k)   = max_agg_diam_white(i,j,k)
        END IF
     END DO
     END DO
     END DO
     !$OMP END PARALLEL DO
  END SUBROUTINE max_agg_diam
  
  !================================================ maximum diameter of agg in non-stratified fluid
  REAL  FUNCTION max_agg_diam_white(i,j,k)
     !-------------------------------------------------------------------------
     !>
     !! maximum aggregate diameter in a non-stratified fluid - following the 
     !! White drag approaximation by Jiang & Logan 1991, assuming agg_re_crit > 10
     !! (otherwise AJX,BJX needs to be adjusted)
     !!

     IMPLICIT NONE

     INTEGER,INTENT(in) :: i,j,k
     REAL        :: nu_vis

     nu_vis  =  dyn_vis(i,j,k)/rho_aq
     max_agg_diam_white = (agg_Re_crit*nu_vis)**((2. - BJ3)/df_agg(i,j,k))       &
                         & /((4./3.)*(av_rho_p(i,j,k) - rho_aq)/rho_aq        &
                         & *av_dp(i,j,k)**(3. - df_agg(i,j,k))*grav_acc_const    &
                         & /(AJ3*nu_vis**BJ3))**(1./df_agg(i,j,k)) 

  END FUNCTION max_agg_diam_white

  !===================================================================================== mass factor
  REAL  FUNCTION mass_factor(dp,df,rhop)
     !-----------------------------------------------------------------------
     !> 
     !! mass_factor calculates the mass factor for the mass of a single
     !! aggregate
     !!
    IMPLICIT NONE

    REAL, INTENT(in) :: dp
    REAL, INTENT(in) :: df
    REAL, INTENT(in) :: rhop
    
    ! mass factor 
    mass_factor = ONE_SIXTH * PI * dp**(3. - df) * rhop

  END FUNCTION mass_factor


  !===================================================================================== rho_agg
  REAL FUNCTION rho_agg(d,rhop,dp,df,rho)
     !-----------------------------------------------------------------------
     !> 
     !! rho_agg provides the aggregate density
     !!
 
     IMPLICIT NONE

     REAL, INTENT(in) :: d
     REAL, INTENT(in) :: rhop
     REAL, INTENT(in) :: dp
     REAL, INTENT(in) :: df
     REAL, INTENT(in) :: rho

     rho_agg =  (rhop - rho)*(dp/d)**(3. - df) + rho

  END FUNCTION rho_agg

  !===================================================================================== dynvis
  SUBROUTINE dynvis(kpie, kpje, kpke, kbnd, pddpo, omask, ptho, psao, ppo)
     !-----------------------------------------------------------------------
     !> 
     !! dynvis calculates the molecular dynamic viscosity according to
     !! Richards 1998: The effect of temperature, pressure, and salinity 
     !! on sound attenuation in turbid seawater. J. Acoust. Soc. Am. 103 (1), 
     !! originally published by  Matthaeus, W. (1972): Die Viskositaet des 
     !! Meerwassers. Beitraege zur Meereskunde, Heft 29 (in German).
     !!

     IMPLICIT NONE

     INTEGER, INTENT(in)  :: kpie                  !< 1st REAL of model grid.
     INTEGER, INTENT(in)  :: kpje                  !< 2nd REAL of model grid.
     INTEGER, INTENT(in)  :: kpke                  !< 3rd (vertical) REAL of model grid.
     INTEGER, INTENT(in)  :: kbnd

     REAL, INTENT(in) :: pddpo(kpie,kpje,kpke) !< size of scalar grid cell (3rd dimension) [m]
     REAL, INTENT(in) :: omask(kpie,kpje) 
     REAL, INTENT(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) !< potential temperature [deg C]
     REAL, INTENT(in) :: psao(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)  !< salinity [psu.].
     REAL, INTENT(in) :: ppo(kpie,kpje,kpke)  !< pressure [Pa].

     ! Local variables
     REAL:: press_val  ! Pascal/rho -> dbar 
     REAL:: ptho_val,psao_val
     INTEGER :: kch
     kch = 0   
     !$OMP PARALLEL DO PRIVATE(i,j,k,press_val,ptho_val,psao_val,kch)
     DO j = 1,kpje
     DO i = 1,kpie
     DO k = 1,kpke
       IF(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) THEN
         kch = MERGE(k+1,k,k<kpke)
         IF(pddpo(i,j,kch) > 0.5) THEN
            press_val    = 0.5*(ppo(i,j,k)  + ppo(i,j,kch))*1.e-5 ! Pascal -> dbar
            ptho_val     = 0.5*(ptho(i,j,k) + ptho(i,j,kch))
            psao_val     = 0.5*(psao(i,j,k) + ptho(i,j,kch))  
         ELSE
            press_val    = ppo(i,j,k)*1.e-5 ! Pascal -> dbar 
            ptho_val     = ptho(i,j,k)
            psao_val     = psao(i,j,k)  
         END IF

     
         ! molecular dynamic viscosity
         dyn_vis(i,j,k) = 0.1    & ! Unit: g / (cm*s) -> kg / (m*s)
             &     *(1.79e-2                                                                &
             &     - 6.1299e-4*ptho_val + 1.4467e-5*ptho_val**2.                      &
             &     - 1.6826e-7*ptho_val**3.                                              &
             &     - 1.8266e-7*press_val  + 9.8972e-12*press_val**2.                  &
             &     + 2.4727e-5*psao_val                                                     &
             &     + psao_val*(4.8429e-7*ptho_val - 4.7172e-8*ptho_val**2.            &
             &     + 7.5986e-10*ptho_val**3.)                                            &
             &     + press_val*(1.3817e-8*ptho_val - 2.6363e-10*ptho_val**2.)         &
             &     - press_val**2.*(6.3255e-13*ptho_val - 1.2116e-14*ptho_val**2.))
       END IF
     END DO
     END DO
     END DO
     !$OMP END PARALLEL DO
  END SUBROUTINE dynvis
 END MODULE mo_m4ago

