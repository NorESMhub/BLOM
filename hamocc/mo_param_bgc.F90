! Copyright (C) 2001  Ernst Maier-Reimer, S. Legutke
! Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger, I. Kriest,
!                     A. Moree, C. Heinze
!
! This file is part of BLOM/iHAMOCC.
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
! along with BLOM. If not, see https://www.gnu.org/licenses/.

module mo_param_bgc

  !*************************************************************************************************
  ! BELEG_PARM - now mo_param_bgc - initialize bgc parameters.
  !
  ! Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
  !
  ! Modified
  ! J.Schwinger,        *NORCE Climate, Bergen*    2020-05-19
  !  - split the original BELEG_BGC in two parts, BELEG_PARM and BELEG_VARS
  ! jmaerz
  !  - rename beleg_parm to mo_param_bgc
  !  T. Bourgeois,     *NORCE climate, Bergen*   2025-04-14
  !  - implement R2OMIP protocol
  !*************************************************************************************************

  use mo_kind,        only: rp
  use mo_carbch,      only: atm_co2
  use mo_control_bgc, only: io_stdo_bgc,bgc_namelist,use_AGG,use_natDIC,                           &
                            use_BROMO,use_cisonew,use_WLIN,use_FB_BGC_OCE,                         &
                            do_ndep,do_oalk,do_rivinpt,do_sedspinup,l_3Dvarsedpor,                 &
                            use_BOXATM,use_CFC,use_PBGC_CK_TIMESTEP,                               &
                            use_sedbypass,with_dmsph,use_PBGC_OCNP_TIMESTEP,ocn_co2_type,use_M4AGO,&
                            do_n2o_coupled,do_nh3_coupled,use_extNcycle,                           &
                            lkwrbioz_off,lTO2depremin,use_shelfsea_res_time,use_sediment_quality,  &
                            use_pref_tracers,use_coupler_ndep,use_river2omip,use_DOMclasses,       &
                            linit_DOMclasses_sim,ldyn_sed_age,sedspin_yr_s,sedspin_yr_e,           &
                            sedspin_ncyc,ldtbgc
  use mod_xc,         only: mnproc,xchalt

  implicit none
  private

  ! For writing parameters to netcdf file
  integer, parameter,public                   :: nlength = 30
  type,public :: pinfo
    character(len=nlength) :: pname
    real(rp)                   :: pvalue
  end type
  integer, protected, public                  :: nentries = 0
  type(pinfo), allocatable, protected, public :: param4nc(:)

  type, public :: ctrinfo
    character(len=nlength) :: cname
    integer                :: cvalue
  end type
  integer, protected, public                    :: centries = 0
  type(ctrinfo), allocatable, protected, public :: controls4nc(:)

  ! Routines
  public  :: ini_parambgc
  public  :: ini_bgctimes
  private :: ini_aggregation
  private :: read_bgcnamelist
  private :: calc_param_atm
  private :: calc_param_biol
  private :: rates_2_timestep

  ! Module variables set by bgcparams namelist
  public :: wpoc_const,wcal_const,wopal_const,wdust_const
  public :: bkopal,bkphy,bluefix,bkzoo
  public :: drempoc,dremopal,dremcalc,dremn2o,dremsul
  public :: drempoc_anaerob,bkox_drempoc
  public :: grazra,gammap,gammaz,spemor
  public :: ecan,epsher,fetune
  public :: relaxfe,rcalc,ropal
  public :: wmin,wmax,wlin,zinges

  ! Other module variables
  public :: ro2ut,rcar,rnit,rnoi,riron,rdnit0,rdnit1,rdnit2,rdn2o1,rdn2o2
  public :: rcar_tdoclc,rhyd_tdoclc,roxy_tdoclc,rnit_tdoclc,ro2ut_tdoclc,ro2utammo_tdoclc
  public :: rcar_tdochc,rhyd_tdochc,roxy_tdochc,rnit_tdochc,ro2ut_tdochc,ro2utammo_tdochc
  public :: atm_n2,atm_o2,atm_co2_nat,atm_bromo,re1312,atm_n2o,atm_nh3
  public :: srfdic_min,re14to,prei13,prei14,ctochl
  public :: atten_w,atten_c,atten_uv,atten_f
  public :: frac_ironindust,frac_soliron,fesoly,phytomi,pi_alpha
  public :: dyphy,tf2,tf1,tf0,tff,bifr13_ini,bifr14_ini,c14_t_half
  public :: rbro,fbro1,fbro2,grami
  public :: calmax,remido,rem_tdoclc,rem_tdochc
  public :: dustd1,dustd2,dustd3,dustsink
  public :: SinkExp, FractDim, Stick, cellmass, cellsink
  public :: fsh,fse,alow1, alow2,alow3
  public :: alar1,alar2,alar3,TSFac,TMFac
  public :: vsmall,safe,pupper,plower,zdis,nmldmin
  public :: beta13,alpha14,atm_c13,atm_c14,c14fac,c14dec
  public :: sedict,silsat,disso_poc,disso_sil,disso_caco3
  public :: sed_denit,sed_sulf,calcwei,opalwei,orgwei
  public :: calcdens,opaldens,orgdens,claydens
  public :: dmsp1,dmsp2,dmsp3,dmsp4,dmsp5,dmsp6,dms_gamma
  public :: gammapsl,gammazsl,alphasl,alphasr,docl_remin,docsl_remin,docsr_remin,docr_remin
  public :: POM_remin_q10,opal_remin_q10,POM_remin_Tref,opal_remin_Tref
  public :: O2thresh_aerob,O2thresh_hypoxic,NO3thresh_sulf
  public :: sed_O2thresh_sulf,sed_O2thresh_hypoxic,sed_NO3thresh_sulf
  public :: shelfbreak_depth
  public :: sed_alpha_poc,sed_qual_sc

  ! extended nitrogen cycle
  public :: q10ano3denit,sc_ano3denit,Trefano3denit,rano3denit,bkano3denit,      &
          & rano2anmx,q10anmx,Trefanmx,alphaanmx,bkoxanmx,bkano2anmx,bkanh4anmx, &
          & rano2denit,q10ano2denit,Trefano2denit,bkoxano2denit,bkano2denit,     &
          & ran2odenit,q10an2odenit,Trefan2odenit,bkoxan2odenit,bkan2odenit,     &
          & rdnra,q10dnra,Trefdnra,bkoxdnra,bkdnra,ranh4nitr,q10anh4nitr,        &
          & Trefanh4nitr,bkoxamox,bkanh4nitr,bkamoxn2o,bkyamox,                  &
          & rano2nitr,q10ano2nitr,Trefano2nitr,bkoxnitr,bkano2nitr,n2omaxy,      &
          & n2oybeta,NOB2AOAy,bn2o,mufn2o,                                       &
          & rc2n,ro2nnit,rnoxp,rnoxpi,rno2anmx,rno2anmxi,rnh4anmx,               &
          & rnh4anmxi,rno2dnra,rno2dnrai,rnh4dnra,rnh4dnrai,rnm1,                &
          & bkphyanh4,bkphyano3,bkphosph,bkiron,ro2utammo,                       &
          & q10ano3denit_sed,sc_ano3denit_sed,Trefano3denit_sed,rano3denit_sed,  &
          & bkano3denit_sed,rano2anmx_sed,q10anmx_sed,Trefanmx_sed,alphaanmx_sed,&
          & bkoxanmx_sed,bkano2anmx_sed,bkanh4anmx_sed,rano2denit_sed,           &
          & q10ano2denit_sed,Trefano2denit_sed,bkoxano2denit_sed,bkano2denit_sed,&
          & ran2odenit_sed,q10an2odenit_sed,Trefan2odenit_sed,bkoxan2odenit_sed, &
          & bkan2odenit_sed,rdnra_sed,q10dnra_sed,Trefdnra_sed,bkoxdnra_sed,     &
          & bkdnra_sed,ranh4nitr_sed,q10anh4nitr_sed,Trefanh4nitr_sed,           &
          & bkoxamox_sed,bkanh4nitr_sed,bkamoxn2o_sed,bkyamox_sed,               &
          & rano2nitr_sed,q10ano2nitr_sed,Trefano2nitr_sed,bkoxnitr_sed,         &
          & bkano2nitr_sed,n2omaxy_sed,n2oybeta_sed,NOB2AOAy_sed,bn2o_sed,       &
          & mufn2o_sed,POM_remin_q10_sed, POM_remin_Tref_sed,bkox_drempoc_sed,   &
          & max_limiter

  ! Time variables
  public :: sec_per_year,sec_per_day,days_per_year
  !********************************************************************
  ! Time parameters
  !********************************************************************
  ! NOTE: days_per_year and sec_per_year are updated in ini_bgctimes(!)
  real(rp), parameter :: sec_per_day   = 86400._rp           ! [s/d]  seconds per day
  real(rp), protected :: days_per_year = 365._rp             ! [d/yr] days per year
  real(rp), protected :: sec_per_year  = 365._rp*sec_per_day ! [s/yr] seconds per year
  logical             :: lini=.true.
  !********************************************************************
  ! Stoichiometry and fixed parameters
  !********************************************************************

  ! extended redfield ratio declaration
  ! Note: stoichiometric ratios are based on Takahashi et al. (1985)
  ! P:N:C:-O2 + 1:16:122:172
  real(rp), parameter :: ro2ut  = 172._rp                ! Oxygen utilization
  real(rp), parameter :: rcar   = 122._rp                ! mol C per mol P
  real(rp), parameter :: rnit   = 16._rp                 ! mol N per mol P
  real(rp), parameter :: rnoi   = 1._rp/rnit             ! mol P per mol N
  real(rp), parameter :: riron  = 5._rp*rcar*1.e-6_rp    ! fe to P ratio in organic matter

  ! stoichiometric ratios for denitrification from Paulmier et al. 2009, Table 1 and
  ! equation 18. Note that their R_0=ro2ut-2*rnit.
  real(rp), parameter :: rdnit0 = 0.8_rp*ro2ut           ! moles nitrate lost for remineralisation of 1 mole P
  real(rp), parameter :: rdnit1 = 0.8_rp*ro2ut - rnit    ! moles nitrate net  for remineralisation of 1 mole P
  real(rp), parameter :: rdnit2 = 0.4_rp*ro2ut           ! moles N2 released  for remineralisation of 1 mole P

  ! stoichiometric ratios for N2O loss by "intermediate dinitrification". Note that there
  ! is no nitrate created by this process, organic N is released as N2
  real(rp), parameter :: rdn2o1 = 2._rp*ro2ut - 2.5_rp*rnit  ! moles N2O used for remineralisation of 1 mole P
  real(rp), parameter :: rdn2o2 = 2._rp*ro2ut - 2._rp*rnit   ! moles N2 released  for remineralisation of 1 mole P

  ! Decay parameter for C14, HalfLive = 5700 years
  real(rp), parameter :: c14_t_half = 5700._rp*365._rp       ! Half life of 14C [days]

  ! Minimum surface DIC concentration for gas-exchange parameterization
  real(rp), parameter :: srfdic_min = 1.0e-5_rp              ! kmol C m-3

  ! Extended nitrogen cycle
  real(rp), parameter :: max_limiter   = 0.9999_rp          ! maximum in concentrations that can be consumed at once
  real(rp), parameter :: rc2n          = rcar/rnit          ! iHAMOCC C:N ratio
  real(rp), parameter :: ro2utammo     = 140._rp            ! Oxygen utilization per mol detritus during ammonification
  real(rp), parameter :: ro2nnit       = ro2utammo/rnit     !
  real(rp), parameter :: rnoxp         = 280._rp            ! consumption of NOx per mol detritus during denitrification
  real(rp), parameter :: rnoxpi        = 1._rp/rnoxp        ! inverse
  real(rp), parameter :: rno2anmx      = 1144._rp           ! consumption of NO2 per mol organic production by anammox
  real(rp), parameter :: rno2anmxi     = 1._rp/rno2anmx     ! inverse
  real(rp), parameter :: rnh4anmx      = 880._rp            ! consumption of NH4 per mol organic production by anammox
  real(rp), parameter :: rnh4anmxi     = 1._rp/rnh4anmx     ! inverse
  real(rp), parameter :: rno2dnra      = 93._rp + 1._rp/3._rp ! consumption of NO2 per mol OM degradation during DNRA
  real(rp), parameter :: rno2dnrai     = 1._rp/rno2dnra     ! inverse
  real(rp), parameter :: rnh4dnra      = rno2dnra + rnit    ! production of NH4 per mol OM during DNRA
  real(rp), parameter :: rnh4dnrai     = 1._rp/rnh4dnra     ! inverse
  real(rp), parameter :: rnm1          = rnit - 1._rp

  ! Terrestrial dissolved organic carbon tDOC (river2oceanmip)
  ! Low-carbon tDOC
  real(rp), parameter :: rcar_tdoclc  = 276._rp                             ! mol C per mol P
  real(rp), parameter :: rnit_tdoclc  = 25._rp                              ! mol N per mol P
  real(rp), parameter :: rhyd_tdoclc  = 2._rp*rcar_tdoclc+3._rp*rnit_tdoclc+3._rp ! =630. mol H per mol P in low-C tDOC
  real(rp), parameter :: roxy_tdoclc  = rcar_tdoclc + 4._rp                 ! =280. mol O per mol P in low-C tDOC
  real(rp), parameter :: ro2ut_tdoclc = (4._rp*rcar_tdoclc+rhyd_tdoclc-2._rp*roxy_tdoclc+5._rp*rnit_tdoclc+5._rp)  &
                                  & /4._rp                              ! =326. Oxygen utilization per mol tDOClc
                                                                        ! during remineralisation
  real(rp), parameter :: ro2utammo_tdoclc = (4._rp*rcar_tdoclc+rhyd_tdoclc-2._rp*roxy_tdoclc           &
                                  & -3._rp*rnit_tdoclc+5._rp)/4._rp     ! =276. Oxygen utilization per mol tDOClc
                                                                        ! during ammonification
  ! High-carbon tDOC
  real(rp), parameter :: rcar_tdochc  = 2583._rp                            ! mol C per mol P
  real(rp), parameter :: rnit_tdochc  = 103._rp                             ! mol N per mol P
  real(rp), parameter :: rhyd_tdochc  = 2._rp*rcar_tdochc+3._rp*rnit_tdochc+3._rp ! =5478. mol H per mol P
  real(rp), parameter :: roxy_tdochc  = rcar_tdochc + 4._rp                 ! =2587. mol O per mol P
  real(rp), parameter :: ro2ut_tdochc = (4._rp*rcar_tdochc+rhyd_tdochc-2._rp*roxy_tdochc+5._rp*rnit_tdochc+5._rp)  &
                                      & /4._rp                              ! =2789. Oxygen utilization per mol tDOChc
                                                                            ! during remineralisation
  real(rp), parameter :: ro2utammo_tdochc = (4._rp*rcar_tdochc+rhyd_tdochc-2._rp*roxy_tdochc           &
                                      & -3._rp*rnit_tdochc+5._rp)/4._rp     ! =2583. Oxygen utilization per mol tDOChc
                                                                            ! during ammonification

  !********************************************************************
  ! Atmosphere:
  !********************************************************************

  real(rp), protected :: atm_n2      = 802000._rp ! atmosphere dinitrogen concentration
  real(rp), protected :: atm_n2o     = 270.1e3_rp ! atmosphere N2O conc. pre-industrial: 270.1 (+-6ppb) IPCC 2021, p708, provided in ppt,300ppb = 300e3ppt = 3e-7 mol/mol
  real(rp), protected :: atm_nh3     = 0._rp      ! Six & Mikolajewicz 2022: less than 1nmol m-3
  real(rp), protected :: atm_o2      = 196800._rp ! atmosphere oxygen concentration
  real(rp), protected :: atm_co2_nat = 284.32_rp  ! atmosphere CO2 concentration CMIP6 pre-industrial reference
  real(rp), protected :: atm_bromo   = 3.4_rp     ! atmosphere bromophorme concentration
                                           ! For now use 3.4ppt from Hense and Quack (2009; Biogeosciences) NEED TO
                                           ! BE UPDATED WITH Ziska et al. (2013) climatology database
  ! set standard carbon isotope ratios
  real(rp), protected :: re1312     = 0.0112372_rp
  real(rp), protected :: re14to     = 1.170e-12_rp    ! Karlen et al. 1965 / Orr et al. 2017
  ! set preindustr. d13c and bigd14C in atmosphere
  real(rp), protected :: prei13     = -6.5_rp
  real(rp), protected :: prei14     = 0._rp
  real(rp), protected :: c14fac                       ! factor for normalizing 14C tracers (~1e-12)
  real(rp), protected :: c14dec
  ! calculate atm_c13 and atm_c14
  real(rp), protected :: atm_c13
  real(rp), protected :: atm_c14                      ! absolute 14c concentration in preindustrial atmosphere

  !********************************************************************
  ! Water column biogeochemistry
  !********************************************************************

  ! Initialize default biogeochemistry parameters
  ! Note that rates are initialized here in /d or equivalent and
  ! time step adjustment is done after reading the BGCPARAMS namelist

  !********************************************************************
  ! SW-radiation and attenuation parameters
  !********************************************************************

  ! Analog to Moore et al., Deep-Sea Research II 49 (2002), 403-462
  ! 1 kmolP = (122*12/60)*10^6 mg[Chlorophyl]
  real(rp), parameter :: ctochl     = 60._rp             ! C to Chlorophyl ratio
  real(rp), protected :: atten_w    = 0.04_rp            ! yellow substances attenuation in 1/m
  real(rp), protected :: atten_c    = 0.03_rp*rcar*(12._rp/ctochl)*1.e6_rp  ! phytoplankton attenuation in 1/m
  real(rp), protected :: atten_uv   = 0.33_rp
  real(rp), protected :: atten_f    = 0.4_rp             ! fraction of sw-radiation directly absorbed in surface layer
                                                  ! (only if FB_BGC_OCE) [feedback bgc-ocean]

  !********************************************************************
  ! Dust deposition and iron solubility parameters
  !********************************************************************
  !ik weight percent iron in dust deposition times Fe solubility
  ! the latter three values come from Johnson et al., 1997
  real(rp), protected :: fetune          = 1.0_rp        ! factor introduced to tune deposition/solubility
  real(rp), protected :: frac_ironindust = 0.035_rp      ! fraction of total iron in dust (typically 3.5%)
  real(rp), protected :: frac_soliron    = 0.01_rp       ! fraction of total iron that is immediately soluble
  real(rp), protected :: fesoly          = 0.55e-9_rp    ! max. diss. iron concentration in deep water
  real(rp), protected :: relaxfe         = 0.05_rp/365._rp  ! 1/d complexation rate to relax iron concentration to fesoly

  !********************************************************************
  ! Phytoplankton parameters (incl. cyanobacteria)
  !********************************************************************
  real(rp), protected :: phytomi    = 1.e-11_rp          ! kmol/m3 - i.e. 1e-5 mmol P/m3 minimum concentration of phyto plankton (?js)
  real(rp), protected :: pi_alpha   = 0.02_rp*0.4_rp     ! initial slope of production vs irradiance curve (alpha) (0.002 for 10 steps per day)
  real(rp), protected :: bkphy      = 4.e-8_rp           ! kmol/m3 - i.e. 0.04 mmol P/m3 half saturation constant
  real(rp), protected :: dyphy      = 0.004_rp           ! 1/d -mortality rate of phytoplankton

  ! Initial fractionation during photosynthesis
  real(rp), protected :: bifr13_ini = 0.98_rp
  real(rp), protected :: bifr14_ini

  ! N2-Fixation following the parameterization in Kriest and Oschlies, 2015.
  ! Factors tf2, tf1 and tf0 are a polynomial (2nd order)
  ! approximation to the functional relationship by Breitbarth et al. (2007),
  ! for temperature dependence of Trichodesmium growth, their eq. (2)
  ! The relation will be scaled to their max. growth rate, tff.
  ! Note that the second order approx. is basically similar to their
  ! function 2 for T-dependent nitrogen fixation multiplied by 4
  ! (2 [N atoms per mole] * 12 [light hrs per day]/6 [C-atoms per N-atoms])
  real(rp), protected :: bluefix    =  0.005_rp          ! 1/d  ! nitrogen fixation rate by blue green algae (cyanobacteria)
  real(rp), protected :: tf2        = -0.0042_rp
  real(rp), protected :: tf1        =  0.2253_rp
  real(rp), protected :: tf0        = -2.7819_rp
  real(rp), protected :: tff        =  0.2395_rp

  !********************************************************************
  ! Zooplankton parameters
  !********************************************************************
  real(rp), protected :: grami      = 1.e-10_rp       ! kmol/m3 - i.e. 1e-5 mmol P/m3 minimum concentration of zooplankton
  real(rp), protected :: bkzoo      = 1.e-7_rp        ! kmol/m3 - i.e. 0.08 mmol P/m3 half saturation constant
  real(rp), protected :: grazra     = 1.5_rp          ! 1/d - grazing rate
  real(rp), protected :: spemor     = 3._rp*1.e6_rp   ! 1/d - mortality rate
  real(rp), protected :: gammap     = 0.03_rp         ! 1/d - exudation rate
  real(rp), protected :: gammaz     = 0.06_rp         ! 1/d - excretion rate
  real(rp), protected :: ecan       = 0.95_rp         ! fraction of mortality as PO_4
  real(rp), protected :: zinges                       ! dimensionless fraction - assimilation efficiency
  real(rp), protected :: epsher                       ! dimensionless fraction - fraction of grazing egested

  !*************************************************
  ! Extended DOM parameters
  !*************************************************
  real(rp), protected :: gammapsl  = 0.02_rp        ! DOC_sl exudation rate [day-1]
  real(rp), protected :: gammazsl  = 0.03_rp        ! DOC_sl excretion rate [day-1]
  real(rp), protected :: alphasl   = 0.18_rp        ! fraction of DOC_sl converted to DOC_sr []
  real(rp), protected :: alphasr   = 0.19_rp        ! fraction of DOC_sr converted to DOC_r  []
  real(rp), protected :: docl_remin  = 1.7e6_rp     ! theor. remineralization rate of labile DOC [d-1]
  real(rp), protected :: docsl_remin = 5.0e7_rp     ! theor. remineralization rate of semi-labile DOC [d-1]
  real(rp), protected :: docsr_remin = 1.7e17_rp    ! theor. remineralization rate of semi-refractory DOC [d-1]
  real(rp), protected :: docr_remin  = 5.0e26_rp    ! theor. remineralization rate of refractory DOC [d-1]


  !********************************************************************
  ! Shell production (CaCO3 and opal) parameters
  !********************************************************************
  real(rp), protected :: bkopal     = 1.e-5_rp        ! kmol/m3 - i.e. 10 mmol Si/m3 half saturation constant
  real(rp), protected :: rcalc                        ! calcium carbonate to organic phosphorous production ratio
  real(rp), protected :: ropal                        ! opal to organic phosphorous production ratio
  real(rp), protected :: calmax                       ! maximum CaCO3 production fraction

  !********************************************************************
  ! Remineralization and dissolution parameters
  !********************************************************************
  real(rp), parameter :: O2thresh_aerob   = 5.e-8_rp   ! Above O2thresh_aerob aerob remineralization takes place
  real(rp), parameter :: O2thresh_hypoxic = 5.e-7_rp   ! Below O2thresh_hypoxic denitrification and sulfate reduction takes place (default model version)
  real(rp), parameter :: NO3thresh_sulf   = 3.e-6_rp   ! Below NO3thresh_sulf 'sufate reduction' takes place
  real(rp), protected :: remido     = 0.004_rp         ! 1/d - remineralization rate (of DOM)
  real(rp), protected :: rem_tdoclc = 1._rp/(1.5_rp*365._rp) ! 1/d Degradation time scale of low-C tDOC (1.5 yr)
  real(rp), protected :: rem_tdochc = 1._rp/(1.5_rp*365._rp) ! 1/d Degradation time scale of high-C tDOC (1.5 yr)
  ! deep sea remineralisation constants
  real(rp), protected :: drempoc         = 0.025_rp    ! 1/d Aerob remineralization rate detritus
  real(rp), protected :: drempoc_anaerob = 1.25e-3_rp  ! =0.05*drempoc - remin in sub-/anoxic environm. - not be overwritten by M4AGO
  real(rp), protected :: bkox_drempoc    = 1.e-5_rp    ! half-saturation constant for oxygen for ammonification (aerobic remin via drempoc)
  real(rp), protected :: dremopal        = 0.003_rp    ! 1/d Dissolution rate for opal
  real(rp), protected :: dremcalc        = 0.02_rp     ! 1/d Dissolution rate for CaCO3 (applied if Omega_c < 1)
  real(rp), protected :: dremn2o         = 0.01_rp     ! 1/d Remineralization rate of detritus on N2O
  real(rp), protected :: dremsul         = 0.005_rp    ! 1/d Remineralization rate for sulphate reduction
  real(rp), protected :: POM_remin_q10   = 2.1_rp      ! Bidle et al. 2002: Regulation of Oceanic Silicon...
  real(rp), protected :: opal_remin_q10  = 2.6_rp      ! Bidle et al. 2002: Regulation of Oceanic Silicon...
  real(rp), protected :: POM_remin_Tref  = 10._rp      ! [deg C] reference temperatue for Q10-dep. POC remin
  real(rp), protected :: opal_remin_Tref = 10._rp      ! [deg C] reference temperature for Q10-dep. opal dissolution

  !********************************************************************
  ! Extended nitrogen cycle
  !********************************************************************
  ! WATER COLUMN
  ! Phytoplankton growth
  real(rp), protected :: bkphyanh4     = 0.10e-6_rp  ! Half-saturation constant for NH4 uptake by bulk phytoplankton (kmol/m3)
  real(rp), protected :: bkphyano3     = 0.16e-6_rp  ! Half-saturation constant for NO3 uptake by bulk phytoplankton (kmol/m3)
  real(rp), protected :: bkphosph      = 0.01e-6_rp  ! Half-saturation constant for PO4 uptake by bulk phytoplankton (kmol/m3)
  real(rp), protected :: bkiron                      ! = bkphosph*riron - Half-saturation constant for Fe uptake by bulk phytoplankton (kmol/m3)

  ! === Nitrification on NH4
  real(rp), protected :: ranh4nitr     = 0.6_rp      ! Maximum growth rate nitrification on NH4 at reference T (1/d -> 1/dt)
  real(rp), protected :: q10anh4nitr   = 3.3_rp      ! Q10 factor for nitrification on NH4 (-)
  real(rp), protected :: Trefanh4nitr  = 20._rp      ! Reference temperature for nitrification on NH4 (degr C)
  real(rp), protected :: bkoxamox      = 0.333e-6_rp ! Half-saturation constant for oxygen limitation of nitrification on NH4 (kmol/m3)
  real(rp), protected :: bkanh4nitr    = 0.133e-6_rp ! Half-saturation constant for nitrification on NH4 (kmol/m3)
  real(rp), protected :: bkamoxn2o     = 0.1e-6_rp   ! Half saturation constant for NH4 in pathway splitting function N2O for nitrification on NH4 (kmol/m3)
  real(rp), protected :: yield_n2o_inf = 0.077_rp    ! Santoro et al. a_2 = 0.077+-0.07
  real(rp), protected :: mufn2o                      !       = 0.11/(50.*1e6*bkoxamox) !=6.61e-3  0.11/(50*1e6)=2.2e-9 - ~Santoro et al. 2011 with simple MM,
  real(rp), protected :: bn2o                        !       = 0.077/(50.*mufn2o)  !=0.2331 - before set to 0.3 - base fraction entering N2O
  real(rp), protected :: n2omaxy       = 0.003_rp    ! Maximum yield of OM on NH4 nitrification (-)
  real(rp), protected :: n2oybeta      = 18._rp      ! Decay factor for inhibition function for yield during nitrification on NH4 (kmol/m3)
  real(rp), protected :: bkyamox       = 0.333e-6_rp ! Half saturation constant for pathway splitting function OM-yield for nitrification on NH4 (kmol/m3)

  ! === Nitrification on NO2
  real(rp), protected :: rano2nitr     = 0.75_rp     ! Maximum growth rate nitrification on NO2 at reference T (1/d -> 1/dt)
  real(rp), protected :: q10ano2nitr   = 2.7_rp      ! Q10 factor for nitrification on NO2 (-)
  real(rp), protected :: Trefano2nitr  = 20._rp      ! Reference temperature for nitrification on NO2 (degr C)
  real(rp), protected :: bkoxnitr      = 0.788e-6_rp ! Half-saturation constant for oxygen limitation of nitrification on NO2 (kmol/m3)
  real(rp), protected :: bkano2nitr    = 0.287e-6_rp ! Half-saturation constant for NO2 for nitrification on NO2 (kmol/m3)
  real(rp), protected :: NOB2AOAy      = 0.44_rp     ! Ratio of NOB versus AOA yield per energy source ~0.043/0.098 according to Zakem et al. 2022

  ! === Denitrification step NO3 -> NO2:
  real(rp), protected :: rano3denit    = 0.0002_rp   ! Maximum growth rate denitrification on NO3 at reference T (1/d -> 1/dt)
  real(rp), protected :: q10ano3denit  = 2._rp       ! Q10 factor for denitrification on NO3 (-)
  real(rp), protected :: Trefano3denit = 10._rp      ! Reference temperature for denitrification on NO3 (degr C)
  real(rp), protected :: sc_ano3denit  = 0.12e6_rp   ! Shape factor for NO3 denitrification oxygen inhibition function (m3/kmol)
  real(rp), protected :: bkano3denit   = 5.e-6_rp    ! Half-saturation constant for NO3 denitrification (kmol/m3)

  ! === Anammox
  real(rp), protected :: rano2anmx     = 0.001_rp    ! Maximum growth rate for anammox at reference T (1/d -> 1/dt)
  real(rp), protected :: q10anmx       = 1.6_rp      ! Q10 factor for anammox (-)
  real(rp), protected :: Trefanmx      = 10._rp      ! Reference temperature for anammox (degr C)
  real(rp), protected :: alphaanmx     = 0.45e6_rp   ! Shape factor for anammox oxygen inhibition function (m3/kmol)
  real(rp), protected :: bkoxanmx      = 11.3e-6_rp  ! Half-saturation constant for oxygen inhibition function (kmol/m3)
  real(rp), protected :: bkano2anmx    = 5.e-6_rp    ! Half-saturation constant for NO2 limitation (kmol/m3)
  real(rp), protected :: bkanh4anmx                  !   = bkano2anmx * rnh4anmx/rno2anmx !Half-saturation constant for NH4 limitation of anammox (kmol/m3)

  ! === Denitrification step NO2 -> N2O
  real(rp), protected :: rano2denit    = 0.002_rp    ! Maximum growth rate denitrification on NO2 at reference T (1/d -> 1/dt)
  real(rp), protected :: q10ano2denit  = 2.0_rp      ! Q10 factor for denitrification on NO2 (-)
  real(rp), protected :: Trefano2denit = 10._rp      ! Reference temperature for denitrification on NO2 (degr C)
  real(rp), protected :: bkoxano2denit = 2.e-6_rp    ! Half-saturation constant for (quadratic) oxygen inhibition function of denitrification on NO2 (kmol/m3)
  real(rp), protected :: bkano2denit   = 5.6e-6_rp   ! Half-saturation constant for denitrification on NO2 (kmol/m3)

  ! === DNRA NO2 -> NH4
  real(rp), protected :: rdnra         = 0.0002_rp   ! Maximum growth rate DNRA on NO2 at reference T (1/d -> 1/dt)
  real(rp), protected :: q10dnra       = 2._rp       ! Q10 factor for DNRA on NO2 (-)
  real(rp), protected :: Trefdnra      = 10._rp      ! Reference temperature for DNRA (degr C)
  real(rp), protected :: bkoxdnra      = 2.5e-6_rp   ! Half saturation constant for (quadratic) oxygen inhibition function of DNRA on NO2 (kmol/m3)
  real(rp), protected :: bkdnra        = 0.05e-6_rp  ! Half-saturation constant for DNRA on NO2 (kmol/m3)

  ! === Denitrification step N2O -> N2
  real(rp), protected :: ran2odenit    = 0.0002_rp  ! Maximum growth rate denitrification on N2O at reference T (1/d -> 1/dt)
  real(rp), protected :: q10an2odenit  = 3._rp       ! Q10 factor for denitrificationj on N2O (-)
  real(rp), protected :: Trefan2odenit = 10._rp      ! Reference temperature for denitrification on N2O (degr C)
  real(rp), protected :: bkoxan2odenit = 10.e-6_rp   ! Half-saturation constant for (quadratic) oxygen inhibition function of denitrification on N2O (kmol/m3)
  real(rp), protected :: bkan2odenit   = 0.1e-6_rp   ! Half-saturation constant for denitrification on N2O (kmol/m3)

  !SEDIMENT
      ! === Ammonification in the sediment
  real(rp), protected :: POM_remin_q10_sed  = 2.1_rp     ! ammonification Q10 in sediment
  real(rp), protected :: POM_remin_Tref_sed = 10._rp     ! ammonification Tref in sediment
  real(rp), protected :: bkox_drempoc_sed   = 1.e-5_rp   ! half saturation constant for O2 limitation of ammonification in sediment (kmol/m3)

      ! === Nitrification on NH4
  real(rp), protected :: ranh4nitr_sed      = 20._rp      ! Maximum growth rate nitrification on NH4 at reference T (1/d -> 1/dt)
  real(rp), protected :: q10anh4nitr_sed    = 3.3_rp      ! Q10 factor for nitrification on NH4 (-)
  real(rp), protected :: Trefanh4nitr_sed   = 20._rp      ! Reference temperature for nitrification on NH4 (degr C)
  real(rp), protected :: bkoxamox_sed       = 0.333e-6_rp ! Half-saturation constant for oxygen limitation of nitrification on NH4 (kmol/m3)
  real(rp), protected :: bkanh4nitr_sed     = 0.133e-6_rp ! Half-saturation constant for nitrification on NH4 (kmol/m3)
  real(rp), protected :: bkamoxn2o_sed      = 0.5e-6_rp   ! Half saturation constant for NH4 in pathway splitting function N2O for nitrification on NH4 (kmol/m3)
  real(rp), protected :: mufn2o_sed                       ! = 0.11/(50.*1e6*bkoxamox_sed)  !=6.61e-3  0.11/(50*1e6)=2.2e-9 - ~Santoro et al. 2011 with simple MM
  real(rp), protected :: bn2o_sed                         ! = 0.077/(50.*mufn2o_sed)       !=0.2331 - before set to 0.3 - base fraction entering N2O
  real(rp), protected :: n2omaxy_sed        = 0.003_rp    ! Maximum yield of OM on NH4 nitrification (-)
  real(rp), protected :: n2oybeta_sed       = 18._rp      ! Decay factor for inhibition function for yield during nitrification on NH4 (kmol/m3)
  real(rp), protected :: bkyamox_sed        = 0.333e-6_rp ! Half saturation constant for pathway splitting function OM-yield for nitrification on NH4 (kmol/m3)

      ! === Nitrification on NO2
  real(rp), protected :: rano2nitr_sed      = 20._rp      ! Maximum growth rate nitrification on NO2 at reference T (1/d -> 1/dt)
  real(rp), protected :: q10ano2nitr_sed    = 2.7_rp      ! Q10 factor for nitrification on NO2 (-)
  real(rp), protected :: Trefano2nitr_sed   = 20._rp      ! Reference temperature for nitrification on NO2 (degr C)
  real(rp), protected :: bkoxnitr_sed       = 0.788e-6_rp ! Half-saturation constant for oxygen limitation of nitrification on NO2 (kmol/m3)
  real(rp), protected :: bkano2nitr_sed     = 0.287e-6_rp ! Half-saturation constant for NO2 for nitrification on NO2 (kmol/m3)
  real(rp), protected :: NOB2AOAy_sed       = 0.44_rp     ! Ratio of NOB versus AOA yield per energy source ~0.043/0.098 according to Zakem et al. 2022

      ! === Denitrification step NO3 -> NO2:
  real(rp), protected :: rano3denit_sed     = 0.3_rp      ! Maximum growth rate denitrification on NO3 at reference T (1/d -> 1/dt)
  real(rp), protected :: q10ano3denit_sed   = 2.57_rp     ! Q10 factor for denitrification on NO3 (-)
  real(rp), protected :: Trefano3denit_sed  = 10._rp      ! Reference temperature for denitrification on NO3 (degr C)
  real(rp), protected :: sc_ano3denit_sed   = 0.12e6_rp   ! Shape factor for NO3 denitrification oxygen inhibition function (m3/kmol)
  real(rp), protected :: bkano3denit_sed    = 5.e-6_rp    ! Half-saturation constant for NO3 denitrification (kmol/m3)

      ! === Anammox
  real(rp), protected :: rano2anmx_sed      = 0.84_rp     ! Maximum growth rate for anammox at reference T (1/d -> 1/dt)
  real(rp), protected :: q10anmx_sed        = 2.12_rp     ! Q10 factor for anammox (-)
  real(rp), protected :: Trefanmx_sed       = 10._rp      ! Reference temperature for anammox (degr C)
  real(rp), protected :: alphaanmx_sed      = 0.45e6_rp   ! Shape factor for anammox oxygen inhibition function (m3/kmol)
  real(rp), protected :: bkoxanmx_sed       = 11.3e-6_rp  ! Half-saturation constant for oxygen inhibition function (kmol/m3)
  real(rp), protected :: bkano2anmx_sed     = 5.e-6_rp    ! Half-saturation constant for NO2 limitation (kmol/m3)
  real(rp), protected :: bkanh4anmx_sed                   ! = bkano2anmx_sed * rnh4anmx/rno2anmx !Half-saturation constant for NH4 limitation of anammox (kmol/m3)

      ! === Denitrification step NO2 -> N2O
  real(rp), protected :: rano2denit_sed     = 2.2_rp      ! Maximum growth rate denitrification on NO2 at reference T (1/d -> 1/dt)
  real(rp), protected :: q10ano2denit_sed   = 2.97_rp     ! Q10 factor for denitrification on NO2 (-)
  real(rp), protected :: Trefano2denit_sed  = 10._rp      ! Reference temperature for denitrification on NO2 (degr C)
  real(rp), protected :: bkoxano2denit_sed  = 2.e-6_rp    ! Half-saturation constant for (quadratic) oxygen inhibition function of denitrification on NO2 (kmol/m3)
  real(rp), protected :: bkano2denit_sed    = 5.6e-6_rp   ! Half-saturation constant for denitrification on NO2 (kmol/m3)

      ! === DNRA NO2 -> NH4
  real(rp), protected :: rdnra_sed          = 0.5_rp      ! Maximum growth rate DNRA on NO2 at reference T (1/d -> 1/dt)
  real(rp), protected :: q10dnra_sed        = 2._rp       ! Q10 factor for DNRA on NO2 (-)
  real(rp), protected :: Trefdnra_sed       = 10._rp      ! Reference temperature for DNRA (degr C)
  real(rp), protected :: bkoxdnra_sed       = 2.5e-6_rp   ! Half saturation constant for (quadratic) oxygen inhibition function of DNRA on NO2 (kmol/m3)
  real(rp), protected :: bkdnra_sed         = 0.05e-6_rp  ! Half-saturation constant for DNRA on NO2 (kmol/m3)

      ! === Denitrification step N2O -> N2
  real(rp), protected :: ran2odenit_sed     = 2.8_rp      ! Maximum growth rate denitrification on N2O at reference T (1/d -> 1/dt)
  real(rp), protected :: q10an2odenit_sed   = 2.37_rp     ! Q10 factor for denitrification on N2O (-)
  real(rp), protected :: Trefan2odenit_sed  = 10._rp      ! Reference temperature for denitrification on N2O (degr C)
  real(rp), protected :: bkoxan2odenit_sed  = 5.e-6_rp    ! Half-saturation constant for (quadratic) oxygen inhibition function of denitrification on N2O (kmol/m3)
  real(rp), protected :: bkan2odenit_sed    = 1.e-6_rp    ! Half-saturation constant for denitrification on N2O (kmol/m3)

  !********************************************************************
  ! Parameters for DMS and BrO schemes
  !********************************************************************
  ! Set constants for dms scheme following Kloster et al. (2006), Table 1
  real(rp), protected :: dmsp1 = 10._rp             ! 2*5. production with temp
  real(rp), protected :: dmsp2 = 0.0011_rp
  real(rp), protected :: dmsp3 = 0.0864_rp          ! bacterial removal, but reduced 50% to increase dms emissions
  real(rp), protected :: dmsp4 = 1.25_rp*0.10_rp    ! production with delcar, but reduced by ~7%
  real(rp), protected :: dmsp5 = 1.25_rp*0.02_rp    ! production with delsil, but increased by a factor of ~2
  real(rp), protected :: dmsp6 = 0.1e-07_rp         ! half saturation microbial

  ! Scaling factor for pH dependency (used if with_dmsph=.true.)
  real(rp), parameter :: dms_gamma = 0.87_rp

  !Bromoform to phosphate ratio (Hense and Quack, 2009)
  !JT: too little production: 0.25Gmol/yr     rbro=6.72e-7*rnit
  !      rbro=2.*6.72e-7*rnit
  !JT Following discussion with B. Quack and D. Booge (01.07.2021), we agree to use 2.4e-6
  real(rp), protected :: rbro       = 2.4e-6_rp*rnit
  real(rp), protected :: fbro1      = 1.0_rp
  real(rp), protected :: fbro2      = 1.0_rp

  !********************************************************************
  ! Sinking parameters
  !********************************************************************
  real(rp), protected :: wpoc_const  =  5._rp           ! m/d   Sinking speed of detritus iris : 5.
  real(rp), protected :: wcal_const  = 30._rp           ! m/d   Sinking speed of CaCO3 shell material
  real(rp), protected :: wopal_const = 30._rp           ! m/d   Sinking speed of opal iris : 60
  real(rp), protected :: wdust_const                    ! m/d   Sinking speed of dust
  real(rp), protected :: wmin        =  5._rp           ! m/d   minimum sinking speed
  real(rp), protected :: wmax        = 60._rp           ! m/d   maximum sinking speed
  real(rp), protected :: wlin        = 60._rp/2400._rp  ! m/d/m constant describing incr. with depth, r/a=1.0
  real(rp), protected :: dustd1      = 0.0001_rp        ! cm = 1 um, boundary between clay and silt
  real(rp), protected :: dustd2                         ! dust diameter squared
  real(rp), protected :: dustd3                         ! dust diameter cubed
  real(rp), protected :: dustsink                       ! sinking speed of dust (used use_AGG)

  real(rp), protected :: SinkExp, FractDim, Stick, cellmass
  real(rp), protected :: fsh, fse,alow1, alow2,alow3,alar1,alar2,alar3,TSFac,TMFac
  real(rp), protected :: vsmall,safe,pupper,plower,zdis,nmldmin
  real(rp), protected :: cellsink = 9999._rp

  !********************************************************************
  ! Shelfsea water residence time
  !********************************************************************
  real(rp), protected :: shelfbreak_depth = 200._rp ! [m] shelf-break depth fall-back value, if no shelfseaa mask file provided

  !********************************************************************
  ! Sediment biogeochemistry
  !********************************************************************
  ! Note that the rates in the sediment are given in per second here!
  !
  real(rp), protected :: sed_O2thresh_hypoxic = 1.e-6_rp   ! Below sed_O2thresh_hypoxic denitrification takes place (default model version)
  real(rp), protected :: sed_O2thresh_sulf    = 3.e-6_rp   ! Below sed_O2thresh_sulf 'sulfate reduction' takes place
  real(rp), protected :: sed_NO3thresh_sulf   = 3.e-6_rp   ! Below sed_NO3thresh_sulf 'sufate reduction' takes place
  real(rp), protected :: sedict      = 1.e-9_rp            ! m2/s Molecular diffusion coefficient
  real(rp), protected :: silsat      = 0.001_rp            ! kmol/m3 Silicate saturation concentration is 1 mol/m3
  real(rp), protected :: disso_poc   = 1.3e-6_rp           ! 1/(kmol O2/m3 s)      Degradation rate constant of POP
! real(rp), protected :: disso_poc   = 0.19_rp/sec_per_day ! 1/(kmol O2/m3 s)      Degradation rate constant of POP
  real(rp), protected :: disso_sil   = 1.4e-7_rp           ! 1/(kmol Si(OH)4/m3 s) Dissolution rate constant of opal
  real(rp), protected :: disso_caco3 = 1.e-7_rp            ! 1/(kmol CO3--/m3 s) Dissolution rate constant of CaCO3
  real(rp), protected :: sed_denit   = 0.01_rp/sec_per_day ! 1/s Denitrification rate constant of POP
  real(rp), protected :: sed_sulf    = 0.01_rp/sec_per_day ! 1/s "Sulfate reduction" rate constant of POP
  real(rp), protected :: sed_alpha_poc = 1._rp/90._rp         ! 1/d 1/decay time for sediment moving average - assuming ~3 month memory here
  real(rp), protected :: sed_qual_sc = 1._rp               ! scaling factor for sediment quality-based remineralization
  !********************************************************************
  ! Densities etc. for SEDIMENT SHIFTING
  !********************************************************************
  ! define weight of calcium carbonate, opal, and poc [kg/kmol]
  real(rp), parameter :: calcwei = 100._rp    ! 40+12+3*16 kg/kmol C
  real(rp), parameter :: opalwei = 60._rp     ! 28 + 2*16  kg/kmol Si
  real(rp), parameter :: orgwei  = 30._rp     ! from 12 kg/kmol * 2.5 POC[kg]/DW[kg]
                                              ! after Alldredge, 1998:
                                              ! POC(g)/DW(g) = 0.4 of diatom marine snow, size 1mm3

  ! define densities of opal, caco3, poc [kg/m3]
  real(rp), parameter :: calcdens = 2600._rp
  real(rp), parameter :: opaldens = 2200._rp
  real(rp), parameter :: orgdens  = 1000._rp
  real(rp), parameter :: claydens = 2600._rp  ! quartz

  !********************************************************************
  ! Module-wide variables used in more than one subroutine
  !********************************************************************
  real(rp) :: beta13, alpha14, d14cat, d13c_atm

contains

  !********************************************************************
  subroutine ini_bgctimes(nday_in_year)

    use mod_config, only: expcnf
    ! NOTE: called also at run time after initialization
    integer,intent(in) :: nday_in_year

    days_per_year = real(nday_in_year)

    if ((nday_in_year /= 365) .and. (mnproc == 1) .and. (lini .eqv. .true.)) then
      lini=.false.
      if (.not. (expcnf == 'single_column' .or. expcnf == 'fuk95' .or. expcnf == 'channel' .or. expcnf == 'noforcing')) then
        ! for production runs, issue an error and stop
        write (io_stdo_bgc,*) 'Error: Init iHAMOCC time variables: non-standard calendar selected with [days] ',days_per_year
        call xchalt('(ini_bgctimes)')
        stop '(ini_bgctimes)'
      else
        ! for test cases, pass, but issue a warning
        write (io_stdo_bgc,*) 'WARNING: Init iHAMOCC time variables: non-standard calendar selected with [days] ',days_per_year
      endif
    endif

    sec_per_year = days_per_year*sec_per_day
  end subroutine ini_bgctimes

  !********************************************************************
  subroutine ini_parambgc()
    !
    ! First, Initialze parameters of individual components with default values.
    ! The order of initialization can matter due to interdependcies.
    ! Then read the namelist and initialize dependent parameter values
    ! adjust rates to 'per time step'
    ! Eventually write out the used parameters to the log file
    !

    call ini_param_biol()    ! initialize biological parameters
    if (use_AGG) then
      call ini_aggregation() ! Initialize aggregation module of Iris Kriest (no NML read thus far)
    endif

    call read_bgcnamelist()  ! read the BGCPARAMS namelist
    call calc_param_atm()    ! calculate atmospheric parameters after updating parameters via nml
    call calc_param_biol()   ! potentially readjust namlist parameter-dependent parameters
    call rates_2_timestep()  ! Converting rates from /d... to /dtb

    call write_parambgc()    ! write out used parameters and calculate back rates from /dtb to /d..
  end subroutine ini_parambgc

  !********************************************************************
  subroutine calc_param_atm()
    !
    ! AFTER having read the namelist:
    ! calculate parameters for atmosphere from given parameters
    !
    if (use_cisonew) then
      beta13   = (prei13/1000._rp)+1._rp
      alpha14  = 2._rp*(prei13+25._rp)
      d14cat   = (prei14+alpha14)/(1._rp-alpha14/1000._rp)
      ! calculate atm_c13 and atm_c14
      atm_c13  = beta13*re1312*atm_co2/(1._rp+beta13*re1312)
      d13C_atm = (((atm_c13/(atm_co2-atm_c13))/re1312)-1._rp)*1000._rp
      ! absolute 14c concentration in preindustrial atmosphere
      atm_c14  = ((d14cat/1000._rp)+1._rp)*re14to*atm_co2
      ! factor for normalizing 14C tracers (~1e-12)
      c14fac   = atm_c14/atm_co2
    endif
  end subroutine calc_param_atm

  !********************************************************************
  subroutine ini_param_biol()
    !
    ! BEFORE reading the namelist:
    ! Default parameters that depend on use case
    !
    !*************************************************
    !     Zooplankton parameters
    !*************************************************
    if (use_AGG) then
      zinges  = 0.5_rp        ! dimensionless fraction -assimilation efficiency
      epsher  = 0.9_rp        ! dimensionless fraction -fraction of grazing egested
    else if ((use_WLIN .eqv. .true.) .or. (use_M4AGO .eqv. .true.))  then
      zinges  = 0.7_rp        ! dimensionless fraction -assimilation efficiency
      epsher  = 0.75_rp       ! dimensionless fraction -fraction of grazing egested
    else
      zinges  = 0.6_rp        ! dimensionless fraction -assimilation efficiency
      epsher  = 0.8_rp        ! dimensionless fraction -fraction of grazing egest
    endif

    !*************************************************
    ! Shell production (CaCO3 and opal) parameters
    !*************************************************
    if (use_AGG) then
      rcalc  = 14._rp         ! calcium carbonate to organic phosphorous production ratio
      ropal  = 10.5_rp        ! opal to organic phosphorous production ratio
      calmax = 0.20_rp
    else if ((use_WLIN .eqv. .true.) .or. (use_M4AGO .eqv. .true.)) then
      rcalc  =  8._rp         ! calcium carbonate to organic phosphorous production ratio
      ropal  = 70._rp         ! opal to organic phosphorous production ratio
    else
      rcalc  = 40._rp         ! iris 40 !calcium carbonate to organic phosphorous production ratio
      ropal  = 30._rp         ! iris 25 !opal to organic phosphorous production ratio
    endif

    if (use_extNcycle) then
      sed_sulf = 1.e-9_rp
    endif
    if (use_M4AGO) then
      ! reset drempoc and dremopal for Q10 T-dep remin/dissolution
      drempoc  = 0.10_rp
      dremopal = 0.023_rp
    endif

  end subroutine ini_param_biol

  !********************************************************************
  subroutine read_bgcnamelist()
    !
    ! READ the bgcparams namelist for parameter tuning.
    ! Note that afterward, i) rates need to be adjusted for timestep
    ! and some depending parameters need re-calculation
    !
    integer  :: iounit

    namelist /bgcparams/ bkphy,dyphy,bluefix,bkzoo,grazra,spemor,gammap,gammaz,  &
                         ecan,zinges,epsher,bkopal,rcalc,ropal,                  &
                         remido,drempoc,dremopal,dremcalc,dremn2o,dremsul,       &
                         fetune,fesoly,relaxfe,wmin,wmax,wlin,wpoc_const,        &
                         wcal_const,wopal_const,disso_poc,disso_sil,disso_caco3, &
                         rano3denit,rano2anmx,rano2denit,ran2odenit,rdnra,       &
                         ranh4nitr,rano2nitr,rano3denit_sed,rano2anmx_sed,       &
                         rano2denit_sed,ran2odenit_sed,rdnra_sed,ranh4nitr_sed,  &
                         rano2nitr_sed,atm_nh3,atm_n2o,bkphyanh4,bkphyano3,      &
                         bkphosph,rem_tdoclc,rem_tdochc,                         &
                         q10ano3denit,sc_ano3denit,bkano3denit,q10anmx,alphaanmx,&
                         bkoxanmx,bkano2anmx,q10ano2denit,                       &
                         bkoxano2denit,bkano2denit,q10an2odenit,bkoxan2odenit,   &
                         bkan2odenit,q10dnra,bkoxdnra,bkdnra,q10anh4nitr,        &
                         bkoxamox,bkanh4nitr,q10ano2nitr,bkoxnitr,bkano2nitr,    &
                         q10ano3denit_sed,sc_ano3denit_sed,bkano3denit_sed,      &
                         q10anmx_sed,alphaanmx_sed,bkox_drempoc_sed,             &
                         bkoxanmx_sed,bkano2anmx_sed,q10ano2denit_sed,           &
                         bkoxano2denit_sed,bkano2denit_sed,q10an2odenit_sed,     &
                         bkoxan2odenit_sed,bkan2odenit_sed,q10dnra_sed,          &
                         bkoxdnra_sed,bkdnra_sed,q10anh4nitr_sed,                &
                         bkoxamox_sed,bkanh4nitr_sed,q10ano2nitr_sed,            &
                         bkoxnitr_sed,bkano2nitr_sed,sed_alpha_poc,sed_qual_sc,  &
                         sed_denit,sed_sulf,                                     &
                         sed_O2thresh_hypoxic,sed_O2thresh_sulf,sed_NO3thresh_sulf,&
                         gammapsl,gammazsl,alphasl,alphasr,docl_remin,docsl_remin, &
                         docsr_remin,docr_remin,yield_n2o_inf,bkamoxn2o

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)
      write(io_stdo_bgc,*)'********************************************'
      write(io_stdo_bgc,*) 'iHAMOCC: read namelist bgcparams'
    endif

    open (newunit=iounit, file=bgc_namelist, status='old',action='read')
    read (unit=iounit, nml=BGCPARAMS)
    close(unit=iounit)

    if (mnproc.eq.1) then
      write(io_stdo_bgc,nml=BGCPARAMS)
    endif

  end subroutine read_bgcnamelist

  !********************************************************************
  subroutine calc_param_biol()
    !
    ! AFTER reading the namelist:
    ! calulate parameters that depend on other tunable parameters
    !
    bifr14_ini  = bifr13_ini**2

    dustd2   = dustd1*dustd1
    dustsink = (9.81_rp * sec_per_day / 18._rp               & ! g * sec per day / 18.
               * (claydens - 1025._rp) / 1.567_rp * 1000._rp & ! excess density / dyn. visc.
               * dustd2 * 1.e-4_rp)                            ! m/d

    if (use_extNcycle) then
      bkiron        = bkphosph*riron                 ! Half-saturation constant for Fe uptake by bulk phytoplankton (kmol/m3)
      bkanh4anmx    = bkano2anmx * rnh4anmx/rno2anmx ! Half-saturation constant for NH4 limitation of anammox (kmol/m3)
      mufn2o        = 0.11_rp/(50._rp*1.e6_rp*bkoxamox)   ! =6.61e-3  0.11/(50*1e6)=2.2e-9 - ~Santoro et al. 2011 with simple MM,
      bn2o          = yield_n2o_inf/(50._rp*mufn2o)            ! =0.2331 - before set to 0.3 - base fraction entering N2O
      bkanh4anmx_sed = bkano2anmx_sed * rnh4anmx/rno2anmx !Half-saturation constant for NH4 limitation of anammox (kmol/m3)
      mufn2o_sed     = 0.11_rp/(50._rp*1.e6_rp*bkoxamox_sed) !=6.61e-3  0.11/(50*1e6)=2.2e-9 - ~Santoro et al. 2011 with simple MM
      bn2o_sed       = yield_n2o_inf/(50._rp*mufn2o_sed)       !=0.2331 - before set to 0.3 - base fraction entering N2O
      lTO2depremin   = .true.
    endif
    if (use_M4AGO) lTO2depremin = .true.
  end subroutine calc_param_biol

  !********************************************************************
  subroutine rates_2_timestep()
    !
    ! AFTER potential update of rates, convert them to rates per timestep
    !
    use mo_control_bgc, only: dtb,dtbgc

    !********************************************************************
    !     Phytoplankton parameters (incl. cyanobacteria)
    !********************************************************************
    dyphy    = dyphy*dtb       ! 1/d -mortality rate of phytoplankton

    ! nitrogen fixation by blue green algae
    bluefix  = bluefix*dtb     ! 1/d

    if (use_cisonew) then
      c14dec = 1._rp-(log(2._rp)/c14_t_half)*dtb   ! lambda [1/day]; c14dec[-]
    endif

    !********************************************************************
    !     Zooplankton parameters
    !********************************************************************
    grazra   = grazra*dtb      ! 1/d to 1/time step - grazing rate
    spemor   = spemor*dtb      ! 1/d to 1/time step - mortality rate
    gammap   = gammap*dtb      ! 1/d to 1/time step - exudation rate
    gammaz   = gammaz*dtb      ! 1/d to 1/time step - excretion rate
    if (use_DOMclasses) then
      gammapsl    = gammapsl*dtb      ! 1/d to 1/time step - exudation rate
      gammazsl    = gammazsl*dtb      ! 1/d to 1/time step - exudation rate
      docl_remin  = docl_remin*dtb
      docsl_remin = docsl_remin*dtb
      docsr_remin = docsr_remin*dtb
      docr_remin  = docr_remin*dtb
    endif

    !********************************************************************
    !     Remineralization and dissolution parameters
    !********************************************************************
    remido   = remido*dtb     ! 1/d to 1/time step - remineralization rate (of DOM)
    ! deep sea remineralisation constants
    drempoc  = drempoc*dtb    ! 1/d to 1/time step  Aerob remineralization rate of detritus
    drempoc_anaerob = drempoc_anaerob*dtb ! 1/d Anaerob remin rate of detritus
    dremopal = dremopal*dtb   ! 1/d to 1/time step  Dissolution rate of opal
    dremcalc = dremcalc*dtb   ! 1/d to 1/time step  Dissolution rate of CaCO3
    dremn2o  = dremn2o*dtb    ! 1/d to 1/time step  Remineralization rate of detritus on N2O
    dremsul  = dremsul*dtb    ! 1/d to 1/time step  Remineralization rate for sulphate reduction
    rem_tdoclc = rem_tdoclc*dtb ! 1/d to 1/time step - remineralisation time scale of terrestrial DOC
    rem_tdochc = rem_tdochc*dtb ! 1/d to 1/time step - remineralisation time scale of terrestrial DOC

    !********************************************************************
    !     Parameters for DMS and BrO schemes
    !********************************************************************
    dmsp2 = dmsp2*dtb
    dmsp3 = dmsp3*dtb         ! 1/d to 1/time-step bacterial removal

    !********************************************************************
    !     Dust deposition and iron solubility parameters
    !********************************************************************
    relaxfe  = relaxfe*dtb    ! 1/d to 1/time step  iron complexation rate

    !********************************************************************
    !     Sinking parameters
    !********************************************************************
    wpoc_const  = wpoc_const*dtb  ! m/d to m/time step      Sinking speed detritusiris : 5.
    wcal_const  = wcal_const*dtb  ! m/d to m/time step      Sinking speed CaCO3
    wopal_const = wopal_const*dtb ! m/d to m/time step      Sinking speed opal iris : 60
    wmin        = wmin*dtb        ! m/d to m/time step      minimum sinking speed
    wmax        = wmax*dtb        ! m/d to m/time step      maximum sinking speed
    wlin        = wlin*dtb        ! m/d/m to m/time step/m  constant describing incr. with depth, r/a=1.0
    dustsink    = dustsink*dtb    ! m/d to m/time step      Sinking speed dust (used in use_AGG)
    wdust_const = dustsink        ! m/d to m/time step      Sinking speed dust

    if(dustsink.gt.cellsink .and. use_AGG) then
      if (mnproc.eq.1)then
        write(io_stdo_bgc,*) ' dust sinking speed greater than cellsink'
        write(io_stdo_bgc,*) ' set dust sinking speed to cellsink'
      endif
      dustsink = cellsink
    endif

    !********************************************************************
    !     Sediment rates
    !********************************************************************
    sedict      = sedict      * dtbgc ! m2/time step                  Molecular diffusion coefficient
    disso_sil   = disso_sil   * dtbgc ! 1/(kmol Si(OH)4/m3 time step) Dissolution rate constant of opal
    disso_poc   = disso_poc   * dtbgc ! 1/(kmol O2/m3 time step)      Degradation rate constant of POP
    disso_caco3 = disso_caco3 * dtbgc ! 1/(kmol CO3--/m3 time step)   Dissolution rate constant of CaCO3
    sed_denit   = sed_denit   * dtbgc ! 1/time step                   Denitrification rate constant of POP
    sed_sulf    = sed_sulf    * dtbgc ! 1/time step                   "Sulfate reduction" rate constant of POP

    if (use_extNcycle) then
      rano3denit = rano3denit *dtb ! Maximum growth rate denitrification on NO3 at reference T (1/d -> 1/dt)
      rano2anmx  = rano2anmx  *dtb ! Maximum growth rate for anammox at reference T (1/d -> 1/dt)
      rano2denit = rano2denit *dtb ! Maximum growth rate denitrification on NO2 at reference T (1/d -> 1/dt)
      ran2odenit = ran2odenit *dtb ! Maximum growth rate denitrification on N2O at reference T (1/d -> 1/dt)
      rdnra      = rdnra      *dtb ! Maximum growth rate DNRA on NO2 at reference T (1/d -> 1/dt)
      ranh4nitr  = ranh4nitr  *dtb ! Maximum growth rate nitrification on NH4 at reference T (1/d -> 1/dt)
      rano2nitr  = rano2nitr  *dtb ! Maximum growth rate nitrification on NO2 at reference T (1/d -> 1/dt)

      rano3denit_sed = rano3denit_sed *dtb ! Maximum growth rate denitrification on NO3 at reference T (1/d -> 1/dt)
      rano2anmx_sed  = rano2anmx_sed  *dtb ! Maximum growth rate for anammox at reference T (1/d -> 1/dt)
      rano2denit_sed = rano2denit_sed *dtb ! Maximum growth rate denitrification on NO2 at reference T (1/d -> 1/dt)
      ran2odenit_sed = ran2odenit_sed *dtb ! Maximum growth rate denitrification on N2O at reference T (1/d -> 1/dt)
      rdnra_sed      = rdnra_sed      *dtb ! Maximum growth rate DNRA on NO2 at reference T (1/d -> 1/dt)
      ranh4nitr_sed  = ranh4nitr_sed  *dtb ! Maximum growth rate nitrification on NH4 at reference T (1/d -> 1/dt)
      rano2nitr_sed  = rano2nitr_sed  *dtb ! Maximum growth rate nitrification on NO2 at reference T (1/d -> 1/dt)
    endif
  end subroutine rates_2_timestep

  !********************************************************************
  subroutine ini_aggregation()
    !
    ! parameters needed for the aggregation module
    !
    use mo_control_bgc, only: dtb

    real(rp) :: shear
    
    dustd3   = dustd1**3
    SinkExp  = 0.62_rp
    FractDim = 1.62_rp
    !Stick = 0.40_rp
    !Stick = 0.25_rp
    Stick    = 0.15_rp
    cellmass = 0.012_rp / rnit ![nmol P]
    !ik      cellmass = 0.0039_rp/ rnit ![nmol P] for a 10 um diameter
    cellsink = 1.40_rp *dtb! [m/d]
    !ik      cellsink = 0.911_rp *dtb! [m/d]  for a 10 um diameter
    !shear = 86400._rp !shear in the mixed layer,   1.0  d-1
    !shear = 64800._rp !shear in the mixed layer,   0.75 d-1
    shear   = 43200._rp !shear in the mixed layer,   0.5  d-1
    fsh     = 0.163_rp * shear *dtb
    fse     = 0.125_rp * 3.1415927_rp * cellsink * 100._rp
    alow1   = 0.002_rp !diameter of smallest particle [cm]
    !ik      alow1 = 0.001_rp !diameter of smallest particle [cm]
    alow2   = alow1 * alow1
    alow3   = alow2 * alow1
    !alar1 = 1.0_rp !diameter of largest particle for size dependend aggregation and sinking [cm]
    !alar1 = 0.75_rp !diameter of largest particle for size dependend aggregation and sinking [cm]
    alar1   = 0.5_rp !diameter of largest particle for size dependend aggregation and sinking [cm]
    vsmall  = 1.e-9_rp
    safe    = 1.e-6_rp
    pupper  = safe/((FractDim+safe)*cellmass)
    plower  = 1._rp/(1.1_rp*cellmass)
    zdis    = 0.01_rp / ((FractDim + 0.01_rp)*cellmass)
    nmldmin = 0.1_rp ! minimum particle number in mixed layer

    alar2   = alar1 * alar1
    alar3   = alar2 * alar1
    TSFac   = (alar1/alow1)**SinkExp
    TMFac   = (alar1/alow1)**FractDim

  end subroutine ini_aggregation

  !********************************************************************
  subroutine write_parambgc()
    !
    ! Write parameters
    !
    use mo_control_bgc, only: dtb,dtbgc

    ! Local variables
    real(rp) :: dtbinv,dtbgcinv

    dtbinv    = 1._rp/dtb
    dtbgcinv  = 1._rp/dtbgc

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)
      write(io_stdo_bgc,*) '********************************************'
      write(io_stdo_bgc,*) '* iHAMOCC configuration: '
      call cinfo_add_entry('use_BROMO',              use_BROMO)
      call cinfo_add_entry('use_AGG',                use_AGG)
      call cinfo_add_entry('use_WLIN',               use_WLIN)
      call cinfo_add_entry('use_natDIC',             use_natDIC)
      call cinfo_add_entry('use_CFC',                use_CFC)
      call cinfo_add_entry('use_cisonew',            use_cisonew)
      call cinfo_add_entry('use_extNcycle',          use_extNcycle)
      call cinfo_add_entry('use_PBGC_OCNP_TIMESTEP', use_PBGC_OCNP_TIMESTEP)
      call cinfo_add_entry('use_PBGC_CK_TIMESTEP',   use_PBGC_CK_TIMESTEP)
      call cinfo_add_entry('use_FB_BGC_OCE',         use_FB_BGC_OCE)
      call cinfo_add_entry('use_BOXATM',             use_BOXATM)
      call cinfo_add_entry('use_sedbypass',          use_sedbypass)
      write(io_stdo_bgc,*) '*   ocn_co2_type           = ',ocn_co2_type
      call cinfo_add_entry('do_ndep',                do_ndep)
      call cinfo_add_entry('do_rivinpt',             do_rivinpt)
      call cinfo_add_entry('do_oalk',                do_oalk)
      call cinfo_add_entry('with_dmsph',             with_dmsph)
      call cinfo_add_entry('do_sedspinup',           do_sedspinup)
      call pinfo_add_entry('sedspin_yr_s',           real(sedspin_yr_s))
      call pinfo_add_entry('sedspin_yr_e',           real(sedspin_yr_e))
      call pinfo_add_entry('sedspin_ncyc',           real(sedspin_ncyc))
      call cinfo_add_entry('l_3Dvarsedpor',          l_3Dvarsedpor)
      call cinfo_add_entry('lkwrbioz_off',           lkwrbioz_off)
      call cinfo_add_entry('lTO2depremin',           lTO2depremin)
      call cinfo_add_entry('use_shelfsea_res_time',  use_shelfsea_res_time)
      call cinfo_add_entry('use_sediment_quality',   use_sediment_quality)
      call cinfo_add_entry('ldyn_sed_age',           ldyn_sed_age)
      call cinfo_add_entry('use_M4AGO',              use_M4AGO)
      call cinfo_add_entry('use_pref_tracers',       use_pref_tracers)
      call cinfo_add_entry('use_coupler_ndep',       use_coupler_ndep)
      call cinfo_add_entry('use_river2omip',         use_river2omip)
      call cinfo_add_entry('use_DOMclasses',         use_DOMclasses)
      call cinfo_add_entry('linit_DOMclasses_sim',   linit_DOMclasses_sim)
      call pinfo_add_entry('dtb',                    dtb)
      call pinfo_add_entry('dtbgc',                  dtbgc)
      call pinfo_add_entry('ldtbgc',                 real(ldtbgc))
      if (use_extNcycle) then
        call cinfo_add_entry('do_n2o_coupled',       do_n2o_coupled)
        call cinfo_add_entry('do_nh3_coupled',       do_nh3_coupled)
      endif
      write(io_stdo_bgc,*) '* '
      write(io_stdo_bgc,*) '* Values of MO_PARAM_BGC variables : '
      call pinfo_add_entry('atm_co2',     atm_co2)
      call pinfo_add_entry('srfdic_min',  srfdic_min)
      if (use_cisonew) then
        call pinfo_add_entry('atm_c13',     atm_c13)
        call pinfo_add_entry('d13C_atm',    d13C_atm)
        call pinfo_add_entry('atm_c14',     atm_c14)
        call pinfo_add_entry('bifr13_ini',  bifr13_ini)
        call pinfo_add_entry('bifr14_ini',  bifr14_ini)
        call pinfo_add_entry('c14fac',      c14fac)
        call pinfo_add_entry('prei13',      prei13)
        call pinfo_add_entry('prei14',      prei14)
        call pinfo_add_entry('re1312',      re1312)
        call pinfo_add_entry('re14to',      re14to)
        call pinfo_add_entry('c14_t_half',  c14_t_half)
        call pinfo_add_entry('c14dec',      c14dec)
        call pinfo_add_entry('beta13',      beta13)
        call pinfo_add_entry('alpha14',     alpha14)
        call pinfo_add_entry('d14cat',      d14cat)
        call pinfo_add_entry('c14fac',      c14fac)
      endif
      call pinfo_add_entry('atm_o2',      atm_o2)
      call pinfo_add_entry('atm_n2',      atm_n2)
      call pinfo_add_entry('atm_n2o',     atm_n2o)
      if (use_natDIC) then
        call pinfo_add_entry('atm_co2_nat', atm_co2_nat)
      endif
      if (use_extNcycle) then
        call pinfo_add_entry('atm_nh3',atm_nh3)
      endif
      call pinfo_add_entry('phytomi',     phytomi)
      call pinfo_add_entry('grami',       grami)
      call pinfo_add_entry('remido',      remido*dtbinv)
      call pinfo_add_entry('dyphy',       dyphy*dtbinv)
      call pinfo_add_entry('zinges',      zinges)
      call pinfo_add_entry('epsher',      epsher)
      call pinfo_add_entry('grazra',      grazra*dtbinv)
      call pinfo_add_entry('spemor',      spemor*dtbinv)
      call pinfo_add_entry('gammap',      gammap*dtbinv)
      call pinfo_add_entry('gammaz',      gammaz*dtbinv)
      if (use_DOMclasses) then
        call pinfo_add_entry('gammapsl',      gammapsl*dtbinv)
        call pinfo_add_entry('gammazsl',      gammazsl*dtbinv)
        call pinfo_add_entry('alphasl',       alphasl)
        call pinfo_add_entry('alphasr',       alphasr)
        call pinfo_add_entry('docl_remin',    docl_remin*dtbinv)
        call pinfo_add_entry('docsl_remin',   docsl_remin*dtbinv)
        call pinfo_add_entry('docsr_remin',   docsr_remin*dtbinv)
        call pinfo_add_entry('docr_remin',    docr_remin*dtbinv)
      endif
      call pinfo_add_entry('ecan',        ecan)
      call pinfo_add_entry('pi_alpha',    pi_alpha)
      call pinfo_add_entry('bkphy',       bkphy)
      call pinfo_add_entry('bkzoo',       bkzoo)
      call pinfo_add_entry('bkopal',      bkopal)
      call pinfo_add_entry('wpoc_const',  wpoc_const*dtbinv)
      call pinfo_add_entry('wcal_const',  wcal_const*dtbinv)
      call pinfo_add_entry('wopal_const', wopal_const*dtbinv)
      call pinfo_add_entry('O2thresh_aerob', O2thresh_aerob)
      call pinfo_add_entry('O2thresh_hypoxic',O2thresh_hypoxic)
      call pinfo_add_entry('NO3thresh_sulf',  NO3thresh_sulf)
      call pinfo_add_entry('drempoc',         drempoc*dtbinv)
      if (lTO2depremin) then
        call pinfo_add_entry('POM_remin_q10',   POM_remin_q10)
        call pinfo_add_entry('POM_remin_Tref',  POM_remin_Tref)
        call pinfo_add_entry('bkox_drempoc',    bkox_drempoc)
      endif
      call pinfo_add_entry('drempoc_anaerob', drempoc_anaerob*dtbinv)
      call pinfo_add_entry('dremopal',    dremopal*dtbinv)
      call pinfo_add_entry('dremcalc',    dremcalc*dtbinv)
      call pinfo_add_entry('dremn2o',     dremn2o*dtbinv)
      call pinfo_add_entry('dremsul',     dremsul*dtbinv)
      call pinfo_add_entry('bluefix',     bluefix*dtbinv)
      call pinfo_add_entry('tf0',         tf0)
      call pinfo_add_entry('tf1',         tf1)
      call pinfo_add_entry('tf2',         tf2)
      call pinfo_add_entry('tff',         tff)
      call pinfo_add_entry('ro2ut',       ro2ut)
      call pinfo_add_entry('rcar',        rcar)
      call pinfo_add_entry('rnit',        rnit)
      call pinfo_add_entry('rnoi',        rnoi)
      call pinfo_add_entry('rdnit0',      rdnit0)
      call pinfo_add_entry('rdnit1',      rdnit1)
      call pinfo_add_entry('rdnit2',      rdnit2)
      call pinfo_add_entry('rdn2o1',      rdn2o1)
      call pinfo_add_entry('rdn2o2',      rdn2o2)
      call pinfo_add_entry('rcalc',       rcalc)
      call pinfo_add_entry('ropal',       ropal)
      call pinfo_add_entry('ctochl',      ctochl)
      call pinfo_add_entry('atten_w',     atten_w)
      call pinfo_add_entry('atten_c',     atten_c)
      call pinfo_add_entry('atten_f',     atten_f)
      call pinfo_add_entry('atten_uv',    atten_uv)
      call pinfo_add_entry('fetune',      fetune)
      call pinfo_add_entry('frac_ironindust',frac_ironindust)
      call pinfo_add_entry('frac_soliron',frac_soliron )
      call pinfo_add_entry('riron',       riron)
      call pinfo_add_entry('fesoly',      fesoly)
      call pinfo_add_entry('relaxfe',     relaxfe*dtbinv)
      call pinfo_add_entry('dmsp1',       dmsp1)
      call pinfo_add_entry('dmsp2',       dmsp2*dtbinv)
      call pinfo_add_entry('dmsp3',       dmsp3*dtbinv)
      call pinfo_add_entry('dmsp4',       dmsp4)
      call pinfo_add_entry('dmsp5',       dmsp5)
      call pinfo_add_entry('dmsp6',       dmsp6)
      if (with_dmsph) then
        call pinfo_add_entry('dms_gamma',   dms_gamma)
      endif
      if (use_BROMO) then
        call pinfo_add_entry('rbro',        rbro)
        call pinfo_add_entry('atm_bromo',   atm_bromo)
        call pinfo_add_entry('fbro1',       fbro1)
        call pinfo_add_entry('fbro2',       fbro2)
      endif
      if (use_WLIN .and. .not. use_AGG) then
        call pinfo_add_entry('wmin',        wmin*dtbinv)
        call pinfo_add_entry('wmax',        wmax*dtbinv)
        call pinfo_add_entry('wlin',        wlin*dtbinv)
      endif
      if (.not. use_AGG) then
        call pinfo_add_entry('dustd1',      dustd1)
        call pinfo_add_entry('dustd2',      dustd2)
        call pinfo_add_entry('dustsink',    dustsink*dtbinv)
        call pinfo_add_entry('wdust_const', wdust_const*dtbinv)
      else
        write(io_stdo_bgc,*)
        write(io_stdo_bgc,*) '********************************************'
        write(io_stdo_bgc,*) '* iHAMOCC aggregate sinking scheme:'
        call pinfo_add_entry('alar1',       alar1)
        call pinfo_add_entry('alar2',       alar2)
        call pinfo_add_entry('alar3',       alar3)
        call pinfo_add_entry('alow1',       alow1)
        call pinfo_add_entry('alow2',       alow2)
        call pinfo_add_entry('alow3',       alow3)
        call pinfo_add_entry('calmax',      calmax)
        call pinfo_add_entry('cellmass',    cellmass)
        call pinfo_add_entry('cellsink',    cellsink)
        call pinfo_add_entry('dustd1',      dustd1)
        call pinfo_add_entry('dustd2',      dustd2)
        call pinfo_add_entry('dustd3',      dustd3)
        call pinfo_add_entry('fractdim',    fractdim)
        call pinfo_add_entry('fse',         fse)
        call pinfo_add_entry('fsh',         fsh)
        call pinfo_add_entry('nmldmin',     nmldmin)
        call pinfo_add_entry('plower',      plower)
        call pinfo_add_entry('pupper',      pupper)
        call pinfo_add_entry('safe',        safe)
        call pinfo_add_entry('sinkexp',     sinkexp)
        call pinfo_add_entry('stick',       stick)
        call pinfo_add_entry('tmfac',       tmfac)
        call pinfo_add_entry('tsfac',       tsfac)
        call pinfo_add_entry('vsmall',      vsmall)
        call pinfo_add_entry('zdis',        zdis)
        write(io_stdo_bgc,*) '* Maximum sinking speed for aggregates of '
        write(io_stdo_bgc,*) '* maximum size ', alar1, ' cm is '
        write(io_stdo_bgc,*)   cellsink/dtb*(alar1/alow1)**SinkExp, ' m/day'
        call pinfo_add_entry('dust diameter (cm)',       dustd1)
        call pinfo_add_entry('dust sinking speed (m d^-1)', dustsink / dtb)
      endif
      if (use_M4AGO) then
        call pinfo_add_entry('opal_remin_q10', opal_remin_q10)
        call pinfo_add_entry('opal_remin_Tref',opal_remin_Tref)
      endif
      write(io_stdo_bgc,*)
      write(io_stdo_bgc,*) '********************************************'
      write(io_stdo_bgc,*) '* Values of MO_PARAM_BGC sediment variables : '
      call pinfo_add_entry('sed_O2thresh_hypoxic',sed_O2thresh_hypoxic)
      call pinfo_add_entry('sed_O2thresh_sulf',   sed_O2thresh_sulf)
      call pinfo_add_entry('sed_NO3thresh_sulf',  sed_NO3thresh_sulf)
      call pinfo_add_entry('sedict',      sedict      * dtbgcinv)
      call pinfo_add_entry('disso_poc',   disso_poc   * dtbgcinv)
      call pinfo_add_entry('disso_sil',   disso_sil   * dtbgcinv)
      call pinfo_add_entry('disso_caco3', disso_caco3 * dtbgcinv)
      call pinfo_add_entry('sed_denit',   sed_denit   * dtbgcinv)
      call pinfo_add_entry('sed_sulf',    sed_sulf    * dtbgcinv)
      call pinfo_add_entry('silsat',      silsat)
      call pinfo_add_entry('orgwei',      orgwei)
      call pinfo_add_entry('opalwei',     opalwei)
      call pinfo_add_entry('calcwei',     calcwei)
      call pinfo_add_entry('orgdens',     orgdens)
      call pinfo_add_entry('opaldens',    opaldens)
      call pinfo_add_entry('calcdens',    calcdens)
      call pinfo_add_entry('claydens',    claydens)
      if (use_sediment_quality) then
        call pinfo_add_entry('sed_alpha_poc',    sed_alpha_poc)
        call pinfo_add_entry('sed_qual_sc',      sed_qual_sc)
      endif
      if (lTO2depremin) then
        call pinfo_add_entry('POM_remin_q10_sed', POM_remin_q10_sed)
        call pinfo_add_entry('POM_remin_Tref_sed',POM_remin_Tref_sed)
        call pinfo_add_entry('bkox_drempoc_sed',  bkox_drempoc_sed)
      endif
      if (use_extNcycle) then
        write(io_stdo_bgc,*) '*********************************************************'
        write(io_stdo_bgc,*) '* HAMOCC extended nitrogen cycle parameters water column:'
        call pinfo_add_entry('rc2n',          rc2n)
        call pinfo_add_entry('ro2utammo',     ro2utammo)
        call pinfo_add_entry('ro2nnit',       ro2nnit)
        call pinfo_add_entry('rnoxp',         rnoxp)
        call pinfo_add_entry('rnoxpi',        rnoxpi)
        call pinfo_add_entry('rno2anmx',      rno2anmx)
        call pinfo_add_entry('rno2anmxi',     rno2anmxi)
        call pinfo_add_entry('rnh4anmx',      rnh4anmx)
        call pinfo_add_entry('rnh4anmxi',     rnh4anmxi)
        call pinfo_add_entry('rno2dnra',      rno2dnra)
        call pinfo_add_entry('rno2dnrai',     rno2dnrai)
        call pinfo_add_entry('rnh4dnra',      rnh4dnra)
        call pinfo_add_entry('rnh4dnrai',     rnh4dnrai)
        call pinfo_add_entry('rnm1',          rnm1)
        call pinfo_add_entry('bkphyanh4',     bkphyanh4)
        call pinfo_add_entry('bkphyano3',     bkphyano3)
        call pinfo_add_entry('bkphosph',      bkphosph)
        call pinfo_add_entry('bkiron',        bkiron)
        call pinfo_add_entry('rano3denit',    rano3denit    *dtbinv)
        call pinfo_add_entry('q10ano3denit',  q10ano3denit)
        call pinfo_add_entry('Trefano3denit', Trefano3denit)
        call pinfo_add_entry('sc_ano3denit',  sc_ano3denit)
        call pinfo_add_entry('bkano3denit',   bkano3denit)
        call pinfo_add_entry('rano2anmx',     rano2anmx     *dtbinv)
        call pinfo_add_entry('q10anmx',       q10anmx)
        call pinfo_add_entry('Trefanmx',      Trefanmx)
        call pinfo_add_entry('alphaanmx',     alphaanmx)
        call pinfo_add_entry('bkoxanmx',      bkoxanmx)
        call pinfo_add_entry('bkano2anmx',    bkano2anmx)
        call pinfo_add_entry('yield_n2o_inf', yield_n2o_inf)
        call pinfo_add_entry('bkanh4anmx',    bkanh4anmx)
        call pinfo_add_entry('rano2denit',    rano2denit    *dtbinv)
        call pinfo_add_entry('q10ano2denit',  q10ano2denit)
        call pinfo_add_entry('Trefano2denit', Trefano2denit)
        call pinfo_add_entry('bkoxano2denit', bkoxano2denit)
        call pinfo_add_entry('bkano2denit',   bkano2denit)
        call pinfo_add_entry('ran2odenit',    ran2odenit    *dtbinv)
        call pinfo_add_entry('q10an2odenit',  q10an2odenit)
        call pinfo_add_entry('Trefan2odenit', Trefan2odenit)
        call pinfo_add_entry('bkoxan2odenit', bkoxan2odenit)
        call pinfo_add_entry('bkan2odenit',   bkan2odenit)
        call pinfo_add_entry('rdnra',         rdnra         *dtbinv)
        call pinfo_add_entry('q10dnra',       q10dnra)
        call pinfo_add_entry('Trefdnra',      Trefdnra)
        call pinfo_add_entry('bkoxdnra',      bkoxdnra)
        call pinfo_add_entry('bkdnra',        bkdnra)
        call pinfo_add_entry('ranh4nitr',     ranh4nitr     *dtbinv)
        call pinfo_add_entry('q10anh4nitr',   q10anh4nitr)
        call pinfo_add_entry('Trefanh4nitr',  Trefanh4nitr)
        call pinfo_add_entry('bkoxamox',      bkoxamox)
        call pinfo_add_entry('bkanh4nitr',    bkanh4nitr)
        call pinfo_add_entry('bkamoxn2o',     bkamoxn2o)
        call pinfo_add_entry('mufn2o',        mufn2o)
        call pinfo_add_entry('bn2o',          bn2o)
        call pinfo_add_entry('n2omaxy',       n2omaxy)
        call pinfo_add_entry('n2oybeta',      n2oybeta)
        call pinfo_add_entry('bkyamox',       bkyamox)
        call pinfo_add_entry('rano2nitr',     rano2nitr     *dtbinv)
        call pinfo_add_entry('q10ano2nitr',   q10ano2nitr)
        call pinfo_add_entry('Trefano2nitr',  Trefano2nitr)
        call pinfo_add_entry('bkoxnitr',      bkoxnitr)
        call pinfo_add_entry('bkano2nitr',    bkano2nitr)
        call pinfo_add_entry('NOB2AOAy',      NOB2AOAy)

        write(io_stdo_bgc,*) '****************************************************************'
        write(io_stdo_bgc,*) '* HAMOCC extended nitrogen cycle parameters sediment:'
        call pinfo_add_entry('rano3denit_sed',    rano3denit_sed    *dtbinv)
        call pinfo_add_entry('q10ano3denit_sed',  q10ano3denit_sed)
        call pinfo_add_entry('Trefano3denit_sed', Trefano3denit_sed)
        call pinfo_add_entry('sc_ano3denit_sed',  sc_ano3denit_sed)
        call pinfo_add_entry('bkano3denit_sed',   bkano3denit_sed)
        call pinfo_add_entry('rano2anmx_sed',     rano2anmx_sed     *dtbinv)
        call pinfo_add_entry('q10anmx_sed',       q10anmx_sed)
        call pinfo_add_entry('Trefanmx_sed',      Trefanmx_sed)
        call pinfo_add_entry('alphaanmx_sed',     alphaanmx_sed)
        call pinfo_add_entry('bkoxanmx_sed',      bkoxanmx_sed)
        call pinfo_add_entry('bkano2anmx_sed',    bkano2anmx_sed)
        call pinfo_add_entry('bkanh4anmx_sed',    bkanh4anmx_sed)
        call pinfo_add_entry('rano2denit_sed',    rano2denit_sed    *dtbinv)
        call pinfo_add_entry('q10ano2denit_sed',  q10ano2denit_sed)
        call pinfo_add_entry('Trefano2denit_sed', Trefano2denit_sed)
        call pinfo_add_entry('bkoxano2denit_sed', bkoxano2denit_sed)
        call pinfo_add_entry('bkano2denit_sed',   bkano2denit_sed)
        call pinfo_add_entry('ran2odenit_sed',    ran2odenit_sed    *dtbinv)
        call pinfo_add_entry('q10an2odenit_sed',  q10an2odenit_sed)
        call pinfo_add_entry('Trefan2odenit_sed', Trefan2odenit_sed)
        call pinfo_add_entry('bkoxan2odenit_sed', bkoxan2odenit_sed)
        call pinfo_add_entry('bkan2odenit_sed',   bkan2odenit_sed)
        call pinfo_add_entry('rdnra_sed',         rdnra_sed         *dtbinv)
        call pinfo_add_entry('q10dnra_sed',       q10dnra_sed)
        call pinfo_add_entry('Trefdnra_sed',      Trefdnra_sed)
        call pinfo_add_entry('bkoxdnra_sed',      bkoxdnra_sed)
        call pinfo_add_entry('bkdnra_sed',        bkdnra_sed)
        call pinfo_add_entry('ranh4nitr_sed',     ranh4nitr_sed     *dtbinv)
        call pinfo_add_entry('q10anh4nitr_sed',   q10anh4nitr_sed)
        call pinfo_add_entry('Trefanh4nitr_sed',  Trefanh4nitr_sed)
        call pinfo_add_entry('bkoxamox_sed',      bkoxamox_sed)
        call pinfo_add_entry('bkanh4nitr_sed',    bkanh4nitr_sed)
        call pinfo_add_entry('bkamoxn2o_sed',     bkamoxn2o_sed)
        call pinfo_add_entry('mufn2o_sed',        mufn2o_sed)
        call pinfo_add_entry('bn2o_sed',          bn2o_sed)
        call pinfo_add_entry('n2omaxy_sed',       n2omaxy_sed)
        call pinfo_add_entry('n2oybeta_sed',      n2oybeta_sed)
        call pinfo_add_entry('bkyamox_sed',       bkyamox_sed)
        call pinfo_add_entry('rano2nitr_sed',     rano2nitr_sed     *dtbinv)
        call pinfo_add_entry('q10ano2nitr_sed',   q10ano2nitr_sed)
        call pinfo_add_entry('Trefano2nitr_sed',  Trefano2nitr_sed)
        call pinfo_add_entry('bkoxnitr_sed',      bkoxnitr_sed)
        call pinfo_add_entry('bkano2nitr_sed',    bkano2nitr_sed)
        call pinfo_add_entry('NOB2AOAy_sed',      NOB2AOAy_sed)
      endif ! N-cycle
      if (use_river2omip) then
        write(io_stdo_bgc,*) '****************************************************************'
        write(io_stdo_bgc,*) '* HAMOCC river2omip parameters'
        call pinfo_add_entry('rcar_tdoclc', rcar_tdoclc)
        call pinfo_add_entry('rhyd_tdoclc', rhyd_tdoclc)
        call pinfo_add_entry('roxy_tdoclc', roxy_tdoclc)
        call pinfo_add_entry('rnit_tdoclc', rnit_tdoclc)
        call pinfo_add_entry('ro2ut_tdoclc', ro2ut_tdoclc)
        call pinfo_add_entry('ro2utammo_tdoclc', ro2utammo_tdoclc)
        call pinfo_add_entry('rem_tdoclc',  rem_tdoclc*dtbinv)
        call pinfo_add_entry('rcar_tdochc', rcar_tdochc)
        call pinfo_add_entry('rhyd_tdochc', rhyd_tdochc)
        call pinfo_add_entry('roxy_tdochc', roxy_tdochc)
        call pinfo_add_entry('rnit_tdochc', rnit_tdochc)
        call pinfo_add_entry('ro2ut_tdochc', ro2ut_tdochc)
        call pinfo_add_entry('ro2utammo_tdochc', ro2utammo_tdochc)
        call pinfo_add_entry('rem_tdochc',  rem_tdochc*dtbinv)
      endif
    endif ! mnproc
  end subroutine write_parambgc

  subroutine cinfo_add_entry(controlname,controlval)
    character(len=*), intent(in)  :: controlname
    logical,          intent(in)  :: controlval

    ! local:
    type(ctrinfo),    allocatable :: temp(:)
    integer :: sw

    if (controlval) then
      sw = 1
    else
      sw = 0
    endif

    !Write to log-file
    write(io_stdo_bgc,*) '*   ',controlname,'     = ',controlval

    if (centries==0) then
      allocate(controls4nc(1))
      controls4nc(1)%cname  = controlname
      controls4nc(1)%cvalue = sw
      centries = centries + 1
    else
      allocate(temp(centries+1))
      temp(1:centries) = controls4nc
      call move_alloc(temp,controls4nc)

      centries = centries + 1
      controls4nc(centries)%cname  = controlname
      controls4nc(centries)%cvalue = sw
    endif
  end subroutine

  subroutine pinfo_add_entry(parname,parvalue)
    character(len=*), intent(in)  :: parname
    real(rp),         intent(in)  :: parvalue

    ! local:
    type(pinfo),      allocatable :: temp(:)

    ! write to logfile:
    write(io_stdo_bgc,*) '*   ',parname,'     = ',parvalue

    if (nentries==0) then
      allocate(param4nc(1))
      param4nc(1)%pname  = parname
      param4nc(1)%pvalue = parvalue
      nentries = nentries + 1
    else
      allocate(temp(nentries+1))
      temp(1:nentries) = param4nc
      call move_alloc(temp,param4nc)

      nentries =nentries + 1
      param4nc(nentries)%pname  = parname
      param4nc(nentries)%pvalue = parvalue
    endif
  end subroutine

end module mo_param_bgc
