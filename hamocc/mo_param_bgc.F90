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
  !*************************************************************************************************

  use mo_carbch,      only: atm_co2
  use mo_control_bgc, only: io_stdo_bgc,bgc_namelist,use_AGG,use_natDIC,                           &
                            use_BROMO,use_cisonew,use_WLIN,use_FB_BGC_OCE,                         &
                            do_ndep,do_oalk,do_rivinpt,do_sedspinup,l_3Dvarsedpor,                 &
                            use_BOXATM,use_CFC,use_PBGC_CK_TIMESTEP,                               &
                            use_sedbypass,with_dmsph,use_PBGC_OCNP_TIMESTEP,ocn_co2_type,lm4ago,   &
                            leuphotic_cya,do_ndep_coupled,do_n2onh3_coupled,use_extNcycle
  use mod_xc,         only: mnproc
!  use mo_m4ago,       only: init_m4ago_nml_params, init_m4ago_params
  use mo_param1_bgc,  only: iatmnh3,iatmn2o
!      use mo_extNwatercol,only: extNwatercol_param_init,extNwatercol_param_update,extNwatercol_param_write,                        &
!                                rano3denit,rano2anmx,rano2denit,ran2odenit,rdnra,ranh4nitr,rano2nitr
!      use mo_extNsediment,only: extNsediment_param_init,extNsediment_param_update,extNsediment_param_write,                        &
!                                rano3denit_sed,rano2anmx_sed,rano2denit_sed,ran2odenit_sed,rdnra_sed,ranh4nitr_sed,rano2nitr_sed

  implicit none
  private

  ! Routines
  public  :: ini_parambgc
  private :: ini_aggregation
  private :: read_bgcnamelist
  private :: calc_param_atm
  private :: calc_param_biol
  private :: rates_2_timestep

  ! Module variables set by bgcparams namelist
  public :: wpoc_const,wcal_const,wopal_const,wdust_const
  public :: bkopal,bkphy,bluefix,bkzoo
  public :: drempoc,dremopal,dremn2o,dremsul
  public :: drempoc_anaerob,bkox_drempoc
  public :: grazra,gammap,gammaz,spemor
  public :: ecan,epsher,fetune
  public :: relaxfe,rcalc,ropal
  public :: wmin,wmax,wlin,zinges

  ! Other module variables
  public :: ro2ut,rcar,rnit,rnoi,riron,rdnit0,rdnit1,rdnit2,rdn2o1,rdn2o2
  public :: atm_n2,atm_o2,atm_co2_nat,atm_bromo,re1312,atm_n2o,atm_nh3
  public :: re14to,prei13,prei14,ctochl
  public :: atten_w,atten_c,atten_uv,atten_f
  public :: perc_diron,fesoly,phytomi,pi_alpha
  public :: dyphy,tf2,tf1,tf0,tff,bifr13_ini,bifr14_ini,c14_t_half
  public :: rbro,fbro1,fbro2,grami
  public :: calmax,remido
  public :: dustd1,dustd2,dustd3,dustsink
  public :: SinkExp, FractDim, Stick, cellmass, cellsink
  public :: fsh,fse,alow1, alow2,alow3
  public :: alar1,alar2,alar3,TSFac,TMFac
  public :: vsmall,safe,pupper,plower,zdis,nmldmin
  public :: beta13,alpha14,atm_c13,atm_c14,c14fac,c14dec
  public :: sedict,silsat,disso_poc,disso_sil,disso_caco3
  public :: sed_denit,calcwei,opalwei,orgwei
  public :: calcdens,opaldens,orgdens,claydens
  public :: dmsp1,dmsp2,dmsp3,dmsp4,dmsp5,dmsp6,dms_gamma
  public :: POM_remin_q10,opal_remin_q10,POM_remin_Tref,opal_remin_Tref


  !********************************************************************
  ! Stoichiometry and fixed parameters
  !********************************************************************

  ! extended redfield ratio declaration
  ! Note: stoichiometric ratios are based on Takahashi etal. (1985)
  ! P:N:C:-O2 + 1:16:122:172
  real, parameter :: ro2ut  = 172.                ! Oxygen utilization
  real, parameter :: rcar   = 122.                ! mol C per mol P
  real, parameter :: rnit   = 16.                 ! mol N per mol P
  real, parameter :: rnoi   = 1./rnit             ! mol P per mol N
  real, parameter :: riron  = 5.*rcar*1.e-6       ! fe to P ratio in organic matter

  ! stoichiometric ratios for denitrification from Paulmier et al. 2009, Table 1 and
  ! equation 18. Note that their R_0=ro2ut-2*rnit.
  real, parameter :: rdnit0 = 0.8*ro2ut           ! moles nitrate lost for remineralisation of 1 mole P
  real, parameter :: rdnit1 = 0.8*ro2ut - rnit    ! moles nitrate net  for remineralisation of 1 mole P
  real, parameter :: rdnit2 = 0.4*ro2ut           ! moles N2 released  for remineralisation of 1 mole P

  ! stoichiometric ratios for N2O loss by "intermediate dinitrification". Note that there
  ! is no nitrate created by this process, organic N is released as N2
  real, parameter :: rdn2o1 = 2*ro2ut - 2.5*rnit  ! moles N2O used for remineralisation of 1 mole P
  real, parameter :: rdn2o2 = 2*ro2ut - 2*rnit    ! moles N2 released  for remineralisation of 1 mole P

  ! Decay parameter for C14, HalfLive = 5700 years
  real, parameter :: c14_t_half = 5700.*365.      ! Half life of 14C [days]

  !********************************************************************
  ! Atmosphere:
  !********************************************************************

  real, protected :: atm_n2      = 802000. ! atmosphere dinitrogen concentration
  real, protected :: atm_n2o     = 300e3   ! atmosphere laughing gas mixing ratio around 1980: 300 ppb,provided in ppt,300ppb = 300e3ppt = 3e-7 mol/mol
  real, protected :: atm_nh3     = 0.      ! Six & Mikolajewicz 2022: less than 1nmol m-3
  real, protected :: atm_o2      = 196800. ! atmosphere oxygen concentration
  real, protected :: atm_co2_nat = 284.32  ! atmosphere CO2 concentration CMIP6 pre-industrial reference
  real, protected :: atm_bromo   = 3.4     ! atmosphere bromophorme concentration
                                           ! For now use 3.4ppt from Hense and Quack (2009; Biogeosciences) NEED TO
                                           ! BE UPDATED WITH Ziska et al. (2013) climatology database
  ! set standard carbon isotope ratios
  real, protected :: re1312     = 0.0112372
  real, protected :: re14to     = 1.170e-12       ! Karlen et al. 1965 / Orr et al. 2017
  ! set preindustr. d13c and bigd14C in atmosphere
  real, protected :: prei13     = -6.5
  real, protected :: prei14     = 0.
  real, protected :: c14fac                       ! factor for normalizing 14C tracers (~1e-12)
  real, protected :: c14dec
  ! calculate atm_c13 and atm_c14
  real, protected :: atm_c13
  real, protected :: atm_c14                      ! absolute 14c concentration in preindustrial atmosphere

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
  real, parameter :: ctochl     = 60.             ! C to Chlorophyl ratio
  real, protected :: atten_w    = 0.04            ! yellow substances attenuation in 1/m
  real, protected :: atten_c    = 0.03*rcar*(12./ctochl)*1.e6  ! phytoplankton attenuation in 1/m
  real, protected :: atten_uv   = 0.33
  real, protected :: atten_f    = 0.4             ! fraction of sw-radiation directly absorbed in surface layer
                                                  ! (only if FB_BGC_OCE) [feedback bgc-ocean]

  !********************************************************************
  ! Dust deposition and iron solubility parameters
  !********************************************************************
  !ik weight percent iron in dust deposition times Fe solubility
  ! the latter three values come from Johnson et al., 1997
  real, protected :: fetune     = 0.6             ! factor introduced to tune deposition/solubility
  real, protected :: perc_diron
  real, protected :: fesoly     = 0.5*1.e-9       ! max. diss. iron concentration in deep water
  real, protected :: relaxfe    = 0.05/365.       ! 1/d complexation rate to relax iron concentration to fesoly

  !********************************************************************
  ! Phytoplankton parameters (incl. cyanobacteria)
  !********************************************************************
  real, protected :: phytomi    = 1.e-11          ! kmol/m3 - i.e. 1e-5 mmol P/m3 minimum concentration of phyto plankton (?js)
  real, protected :: pi_alpha   = 0.02*0.4        ! initial slope of production vs irradiance curve (alpha) (0.002 for 10 steps per day)
  real, protected :: bkphy      = 4.e-8           ! kmol/m3 - i.e. 0.04 mmol P/m3 half saturation constant
  real, protected :: dyphy      = 0.004           ! 1/d -mortality rate of phytoplankton

  ! Initial fractionation during photosynthesis
  real, protected :: bifr13_ini = 0.98
  real, protected :: bifr14_ini

  ! N2-Fixation following the parameterization in Kriest and Oschlies, 2015.
  ! Factors tf2, tf1 and tf0 are a polynomial (2nd order)
  ! approximation to the functional relationship by Breitbarth et al. (2007),
  ! for temperature dependence of Trichodesmium growth, their eq. (2)
  ! The relation will be scaled to their max. growth rate, tff.
  ! Note that the second order approx. is basically similar to their
  ! function 2 for T-dependent nitrogen fixation multiplied by 4
  ! (2 [N atoms per mole] * 12 [light hrs per day]/6 [C-atoms per N-atoms])
  real, protected :: bluefix    =  0.005          ! 1/d  ! nitrogen fixation rate by blue green algae (cyanobacteria)
  real, protected :: tf2        = -0.0042
  real, protected :: tf1        =  0.2253
  real, protected :: tf0        = -2.7819
  real, protected :: tff        =  0.2395

  !********************************************************************
  ! Zooplankton parameters
  !********************************************************************
  real, protected :: grami      = 1.e-10          ! kmol/m3 - i.e. 1e-5 mmol P/m3 minimum concentration of zooplankton
  real, protected :: bkzoo      = 8.e-8           ! kmol/m3 - i.e. 0.08 mmol P/m3 half saturation constant

  !ik addded parameter definition; taken from OCPROD.F
  real, protected :: grazra     = 1.2             ! 1/d - grazing rate
  real, protected :: spemor     = 3.*1.e6         ! 1/d - mortality rate
  real, protected :: gammap     = 0.04            ! 1/d - exudation rate
  real, protected :: gammaz     = 0.06            ! 1/d - excretion rate
  real, protected :: ecan       = 0.95            ! fraction of mortality as PO_4
  real, protected :: zinges                       ! dimensionless fraction - assimilation efficiency
  real, protected :: epsher                       ! dimensionless fraction - fraction of grazing egested

  !********************************************************************
  ! Shell production (CaCO3 and opal) parameters
  !********************************************************************
  real, protected :: bkopal     = 5.e-6           ! kmol/m3 - i.e. 4.0 mmol Si/m3 half saturation constant
  real, protected :: rcalc                        ! calcium carbonate to organic phosphorous production ratio
  real, protected :: ropal                        ! opal to organic phosphorous production ratio
  real, protected :: calmax                       ! maximum CaCO3 production fraction

  !********************************************************************
  ! Remineralization and dissolution parameters
  !********************************************************************
  real, protected :: remido     = 0.004           ! 1/d - remineralization rate (of DOM)
  ! deep sea remineralisation constants
  real, protected :: drempoc    = 0.025           ! 1/d Aerob remineralization rate detritus
  real, protected :: drempoc_anaerob = 1.25e-3    ! =0.05*drempoc - remin in sub-/anoxic environm. - not be overwritten by lm4ago
  real, protected :: bkox_drempoc = 1e-7 ! half-saturation constant for oxygen for ammonification (aerobic remin via drempoc)
  real, protected :: dremopal   = 0.003           ! 1/d Dissolution rate for opal
  real, protected :: dremn2o    = 0.01            ! 1/d Remineralization rate of detritus on N2O
  real, protected :: dremsul    = 0.005           ! 1/d Remineralization rate for sulphate reduction
  real, protected :: POM_remin_q10   = 2.1        ! Bidle et al. 2002: Regulation of Oceanic Silicon...
  real, protected :: opal_remin_q10  = 2.6        ! Bidle et al. 2002: Regulation of Oceanic Silicon...
  real, protected :: POM_remin_Tref  = 10.        ! [deg C] reference temperatue for Q10-dep. POC remin
  real, protected :: opal_remin_Tref = 10.        ! [deg C] reference temperature for Q10-dep. opal dissolution


  !********************************************************************
  ! Parameters for DMS and BrO schemes
  !********************************************************************
  ! Set constants for dms scheme following Kloster et al. (2006), Table 1
  real, protected :: dmsp1 = 10.             ! 2*5. production with temp
  real, protected :: dmsp2 = 0.0011
  real, protected :: dmsp3 = 0.0864          ! bacterial removal, but reduced 50% to increase dms emissions
  real, protected :: dmsp4 = 1.25*0.10       ! production with delcar, but reduced by ~7%
  real, protected :: dmsp5 = 1.25*0.02       ! production with delsil, but increased by a factor of ~2
  real, protected :: dmsp6 = 0.100000000E-07 ! half saturation microbial

  ! Scaling factor for pH dependency (used if with_dmsph=.true.)
  real, parameter :: dms_gamma = 0.87

  !Bromoform to phosphate ratio (Hense and Quack, 2009)
  !JT: too little production: 0.25Gmol/yr     rbro=6.72e-7*rnit
  !      rbro=2.*6.72e-7*rnit
  !JT Following discussion with B. Quack and D. Booge (01.07.2021), we agree to use 2.4e-6
  real, protected :: rbro       = 2.4e-6*rnit
  real, protected :: fbro1      = 1.0
  real, protected :: fbro2      = 1.0

  !********************************************************************
  ! Sinking parameters
  !********************************************************************
  real, protected :: wpoc_const  =  5.             ! m/d   Sinking speed of detritus iris : 5.
  real, protected :: wcal_const  = 30.             ! m/d   Sinking speed of CaCO3 shell material
  real, protected :: wopal_const = 30.             ! m/d   Sinking speed of opal iris : 60
  real, protected :: wdust_const                             ! m/d   Sinking speed of dust
  real, protected :: wmin        =  1.             ! m/d   minimum sinking speed
  real, protected :: wmax        = 60.             ! m/d   maximum sinking speed
  real, protected :: wlin        = 60./2400.       ! m/d/m constant describing incr. with depth, r/a=1.0
  real, protected :: dustd1      = 0.0001          ! cm = 1 um, boundary between clay and silt
  real, protected :: dustd2                       ! dust diameter squared
  real, protected :: dustd3                       ! dust diameter cubed
  real, protected :: dustsink                     ! sinking speed of dust (used use_AGG)

  real, protected :: SinkExp, FractDim, Stick, cellmass
  real, protected :: fsh, fse,alow1, alow2,alow3,alar1,alar2,alar3,TSFac,TMFac
  real, protected :: vsmall,safe,pupper,plower,zdis,nmldmin
  real, protected :: cellsink = 9999.


  !********************************************************************
  ! Sediment biogeochemistry
  !********************************************************************
  ! Note that the rates in the sediment are given in per second here!
  !
  real, protected :: sedict      = 1.e-9          ! m2/s Molecular diffusion coefficient
  real, protected :: silsat      = 0.001          ! kmol/m3 Silicate saturation concentration is 1 mol/m3
  real, protected :: disso_poc   = 0.01 / 86400.  ! 1/(kmol O2/m3 s) disso=3.e-5 was quite high - Degradation rate constant of POP
  real, protected :: disso_sil   = 1.e-6          ! 1/(kmol Si(OH)4/m3 s) Dissolution rate constant of opal
  ! THIS NEEDS TO BE CHANGED TO disso=3.e-8! THIS IS ONLY KEPT FOR THE MOMENT
  ! FOR BACKWARDS COMPATIBILITY
  ! disso_sil = 3.e-8*dtbgc  ! (2011-01-04) EMR
  ! disso_sil = 1.e-6*dtbgc  ! test vom 03.03.04 half live sil ca. 20.000 yr
  real, protected :: disso_caco3 = 1.e-7          ! 1/(kmol CO3--/m3 s) Dissolution rate constant of CaCO3
  real, protected :: sed_denit   = 0.01/86400.    ! 1/s Denitrification rate constant of POP

  !********************************************************************
  ! Densities etc. for SEDIMENT SHIFTING
  !********************************************************************
  ! define weight of calcium carbonate, opal, and poc [kg/kmol]
  real, parameter :: calcwei = 100.    ! 40+12+3*16 kg/kmol C
  real, parameter :: opalwei = 60.     ! 28 + 2*16  kg/kmol Si
  real, parameter :: orgwei  = 30.     ! from 12 kg/kmol * 2.5 POC[kg]/DW[kg]
                                       ! after Alldredge, 1998:
                                       ! POC(g)/DW(g) = 0.4 of diatom marine snow, size 1mm3

  ! define densities of opal, caco3, poc [kg/m3]
  real, parameter :: calcdens = 2600.
  real, parameter :: opaldens = 2200.
  real, parameter :: orgdens  = 1000.
  real, parameter :: claydens = 2600.  ! quartz

  !********************************************************************
  ! Module-wide variables used in more than one subroutine
  !********************************************************************
  real :: beta13, alpha14, d14cat, d13c_atm

contains

  !********************************************************************
  subroutine ini_parambgc(kpie,kpje)
    !
    ! First, Initialze parameters of individual components with default values.
    ! The order of initialization can matter due to interdependcies.
    ! Then read the namelist and initialize dependent parameter values
    ! adjust rates to 'per time step'
    ! Eventually write out the used parameters to the log file
    !
    ! Arguments
    integer, intent(in) :: kpie  ! 1st dimension of model grid.
    integer, intent(in) :: kpje  ! 2nd dimension of model grid.

    call ini_param_biol()    ! initialize biological parameters
    if (use_AGG) then
      call ini_aggregation() ! Initialize aggregation module of Iris Kriest (no NML read thus far)
    endif

    call read_bgcnamelist()  ! read the BGCPARAMS namelist
    call calc_param_atm()    ! calculate atmospheric parameters after updating parameters via nml
    call calc_param_biol()   ! potentially readjust namlist parameter-dependent parameters
    call rates_2_timestep()  ! Converting rates from /d... to /dtb
!    call init_m4ago_params()       ! Initialize M4AGO parameters relying on nml parameters
    call write_parambgc()    ! write out used parameters and calculate back rates from /dtb to /d..
  end subroutine ini_parambgc

  !********************************************************************
  subroutine calc_param_atm()
    !
    ! AFTER having read the namelist:
    ! calculate parameters for atmosphere from given parameters
    !
    if (use_cisonew) then
      beta13   = (prei13/1000.)+1.
      alpha14  = 2.*(prei13+25.)
      d14cat   = (prei14+alpha14)/(1.-alpha14/1000.)
      ! calculate atm_c13 and atm_c14
      atm_c13  = beta13*re1312*atm_co2/(1.+beta13*re1312)
      d13C_atm = (((atm_c13/(atm_co2-atm_c13))/re1312)-1.)*1000.
      ! absolute 14c concentration in preindustrial atmosphere
      atm_c14  = ((d14cat/1000.)+1.)*re14to*atm_co2
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
      zinges  = 0.5        ! dimensionless fraction -assimilation efficiency
      epsher  = 0.9        ! dimensionless fraction -fraction of grazing egested
    else if (use_WLIN) then
      zinges  = 0.7        ! dimensionless fraction -assimilation efficiency
      epsher  = 0.85       ! dimensionless fraction -fraction of grazing egested
    else
      zinges  = 0.6        ! dimensionless fraction -assimilation efficiency
      epsher  = 0.8        ! dimensionless fraction -fraction of grazing egest
    endif

    !*************************************************
    ! Shell production (CaCO3 and opal) parameters
    !*************************************************
    if (use_AGG) then
      rcalc  = 14.         ! calcium carbonate to organic phosphorous production ratio
      ropal  = 10.5        ! opal to organic phosphorous production ratio
      calmax = 0.20
    else if (use_WLIN) then
      rcalc  = 33.         ! calcium carbonate to organic phosphorous production ratio
      ropal  = 45.         ! opal to organic phosphorous production ratio
    else
      rcalc  = 40.         ! iris 40 !calcium carbonate to organic phosphorous production ratio
      ropal  = 30.         ! iris 25 !opal to organic phosphorous production ratio
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
                         remido,drempoc,dremopal,dremn2o,dremsul,fetune,relaxfe, &
                         wmin,wmax,wlin,wpoc_const,wcal_const,wopal_const

    open (newunit=iounit, file=bgc_namelist, status='old',action='read')
    read (unit=iounit, nml=BGCPARAMS)
    close(unit=iounit)

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*) 
      write(io_stdo_bgc,*)'********************************************'
      write(io_stdo_bgc,*) 'iHAMOCC: read namelist bgcparams'
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

    perc_diron = fetune * 0.035 * 0.01 / 55.85

    dustd2   = dustd1*dustd1
    dustsink = (9.81 * 86400. / 18.                 & ! g * sec per day / 18.
               * (claydens - 1025.) / 1.567 * 1000. & ! excess density / dyn. visc.
               * dustd2 * 1.e-4)                      ! m/d

    if(lm4ago)then
        ! reset drempoc and dremopal for T-dep remin/dissolution
        drempoc  = 0.12
        dremopal = 0.023
    endif
    ! M4AGO parameters
!    call init_m4ago_nml_params()
    if (use_extNcycle) then
    ! initialize the extended nitrogen cycle parameters - first water column, then sediment, 
    ! since sediment relies on water column parameters for the extended nitrogen cycle 
    ! Sediment also relies on M4AGO being initialized (POM_remin_q10 and POM_remin_Tref)
!    call extNwatercol_param_init()
!    call extNsediment_param_init()
    endif
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
      c14dec = 1.-(log(2.)/c14_t_half)*dtb   ! lambda [1/day]; c14dec[-]
    endif

    !********************************************************************
    !     Zooplankton parameters
    !********************************************************************
    grazra   = grazra*dtb      ! 1/d to 1/time step - grazing rate
    spemor   = spemor*dtb      ! 1/d to 1/time step - mortality rate
    gammap   = gammap*dtb      ! 1/d to 1/time step - exudation rate
    gammaz   = gammaz*dtb      ! 1/d to 1/time step - excretion rate

    !********************************************************************
    !     Remineralization and dissolution parameters
    !********************************************************************
    remido   = remido*dtb     ! 1/d to 1/time step - remineralization rate (of DOM)
    ! deep sea remineralisation constants
    drempoc  = drempoc*dtb    ! 1/d to 1/time step  Aerob remineralization rate of detritus
    drempoc_anaerob = drempoc_anaerob*dtb ! 1/d Anaerob remin rate of detritus
    dremopal = dremopal*dtb   ! 1/d to 1/time step  Dissolution rate of opal
    dremn2o  = dremn2o*dtb    ! 1/d to 1/time step  Remineralization rate of detritus on N2O
    dremsul  = dremsul*dtb    ! 1/d to 1/time step  Remineralization rate for sulphate reduction

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
    if (use_extNcycle) then
!    call extNwatercol_param_update()
!    call extNsediment_param_update()
    endif
  end subroutine rates_2_timestep

  !********************************************************************
  subroutine ini_aggregation()
    !
    ! parameters needed for the aggregation module
    !
    use mo_control_bgc, only: dtb

    real :: shear

    SinkExp  = 0.62
    FractDim = 1.62
    !Stick = 0.40
    !Stick = 0.25
    Stick    = 0.15
    cellmass = 0.012 / rnit ![nmol P]
    !ik      cellmass = 0.0039/ rnit ![nmol P] for a 10 um diameter
    cellsink = 1.40 *dtb! [m/d]
    !ik      cellsink = 0.911 *dtb! [m/d]  for a 10 um diameter
    !shear = 86400. !shear in the mixed layer,   1.0  d-1
    !shear = 64800. !shear in the mixed layer,   0.75 d-1
    shear   = 43200. !shear in the mixed layer,   0.5  d-1
    fsh     = 0.163 * shear *dtb
    fse     = 0.125 * 3.1415927 * cellsink * 100.
    alow1   = 0.002 !diameter of smallest particle [cm]
    !ik      alow1 = 0.001 !diameter of smallest particle [cm]
    alow2   = alow1 * alow1
    alow3   = alow2 * alow1
    !alar1 = 1.0 !diameter of largest particle for size dependend aggregation and sinking [cm]
    !alar1 = 0.75 !diameter of largest particle for size dependend aggregation and sinking [cm]
    alar1   = 0.5 !diameter of largest particle for size dependend aggregation and sinking [cm]
    vsmall  = 1.e-9
    safe    = 1.e-6
    pupper  = safe/((FractDim+safe)*cellmass)
    plower  = 1./(1.1*cellmass)
    zdis    = 0.01 / ((FractDim + 0.01)*cellmass)
    nmldmin = 0.1 ! minimum particle number in mixed layer

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
    real :: dtbinv,dtbgcinv

    dtbinv    = 1./dtb
    dtbgcinv  = 1./dtbgc

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)
      write(io_stdo_bgc,*) '********************************************'
      write(io_stdo_bgc,*) '* iHAMOCC configuration: '
      write(io_stdo_bgc,*) '*   use_BROMO              = ',use_BROMO
      write(io_stdo_bgc,*) '*   use_AGG                = ',use_AGG
      write(io_stdo_bgc,*) '*   use_WLIN               = ',use_WLIN
      write(io_stdo_bgc,*) '*   use_natDIC             = ',use_natDIC
      write(io_stdo_bgc,*) '*   use_CFC                = ',use_CFC
      write(io_stdo_bgc,*) '*   use_cisonew            = ',use_cisonew
      write(io_stdo_bgc,*) '*   use_PBGC_OCNP_TIMESTEP = ',use_PBGC_OCNP_TIMESTEP
      write(io_stdo_bgc,*) '*   use_PBGC_CK_TIMESTEP   = ',use_PBGC_CK_TIMESTEP
      write(io_stdo_bgc,*) '*   use_FB_BGC_OCE BROMO   = ',use_FB_BGC_OCE
      write(io_stdo_bgc,*) '*   use_BOXATM             = ',use_BOXATM
      write(io_stdo_bgc,*) '*   use_sedbypass          = ',use_sedbypass
      write(io_stdo_bgc,*) '*   ocn_co2_type           = ',ocn_co2_type
      write(io_stdo_bgc,*) '*   do_ndep                = ',do_ndep
      write(io_stdo_bgc,*) '*   do_rivinpt             = ',do_rivinpt
      write(io_stdo_bgc,*) '*   do_oalk                = ',do_oalk
      write(io_stdo_bgc,*) '*   with_dmsph             = ',with_dmsph
      write(io_stdo_bgc,*) '*   do_sedspinup           = ',do_sedspinup
      write(io_stdo_bgc,*) '*   l_3Dvarsedpor          = ',l_3Dvarsedpor
      write(io_stdo_bgc,*) '*   leuphotic_cya          = ',leuphotic_cya
      write(io_stdo_bgc,*) '*   lm4ago                 = ',lm4ago
      if (use_extNcycle) then
        write(io_stdo_bgc,*) '*   do_ndep_coupled        = ',do_ndep_coupled
        write(io_stdo_bgc,*) '*   do_n2onh3_coupled      = ',do_n2onh3_coupled
      endif
      write(io_stdo_bgc,*) '* '
      write(io_stdo_bgc,*) '* Values of MO_PARAM_BGC variables : '
      write(io_stdo_bgc,*) '*   atm_co2      = ',atm_co2
      if (use_cisonew) then
        write(io_stdo_bgc,*) '*   atm_c13      = ',atm_c13
        write(io_stdo_bgc,*) '*   d13C_atm     = ',d13C_atm
        write(io_stdo_bgc,*) '*   atm_c14      = ',atm_c14
        write(io_stdo_bgc,*) '*   bifr13_ini   = ',bifr13_ini
        write(io_stdo_bgc,*) '*   bifr14_ini   = ',bifr14_ini
        write(io_stdo_bgc,*) '*   c14fac       = ',c14fac
        write(io_stdo_bgc,*) '*   prei13       = ',prei13
        write(io_stdo_bgc,*) '*   prei14       = ',prei14
        write(io_stdo_bgc,*) '*   re1312       = ',re1312
        write(io_stdo_bgc,*) '*   re14to       = ',re14to
        write(io_stdo_bgc,*) '*   c14_t_half   = ',c14_t_half
        write(io_stdo_bgc,*) '*   c14dec       = ',c14dec
        write(io_stdo_bgc,*) '*   beta13       = ',beta13
        write(io_stdo_bgc,*) '*   alpha14      = ',alpha14
        write(io_stdo_bgc,*) '*   d14cat       = ',d14cat
        write(io_stdo_bgc,*) '*   c14fac       = ',c14fac
      endif
      write(io_stdo_bgc,*) '*   atm_o2       = ',atm_o2
      write(io_stdo_bgc,*) '*   atm_n2       = ',atm_n2
      WRITE(io_stdo_bgc,*) '*   atm_n2o      = ',atm_n2o
      if (use_extNcycle) then
        write(io_stdo_bgc,*) '*   atm_nh3      = ',atm_nh3
      endif
      write(io_stdo_bgc,*) '*   phytomi      = ',phytomi
      write(io_stdo_bgc,*) '*   grami        = ',grami
      write(io_stdo_bgc,*) '*   remido       = ',remido*dtbinv
      write(io_stdo_bgc,*) '*   dyphy        = ',dyphy*dtbinv
      write(io_stdo_bgc,*) '*   zinges       = ',zinges
      write(io_stdo_bgc,*) '*   epsher       = ',epsher
      write(io_stdo_bgc,*) '*   grazra       = ',grazra*dtbinv
      write(io_stdo_bgc,*) '*   spemor       = ',spemor*dtbinv
      write(io_stdo_bgc,*) '*   gammap       = ',gammap*dtbinv
      write(io_stdo_bgc,*) '*   gammaz       = ',gammaz*dtbinv
      write(io_stdo_bgc,*) '*   ecan         = ',ecan
      write(io_stdo_bgc,*) '*   pi_alpha     = ',pi_alpha
      write(io_stdo_bgc,*) '*   bkphy        = ',bkphy
      write(io_stdo_bgc,*) '*   bkzoo        = ',bkzoo
      write(io_stdo_bgc,*) '*   bkopal       = ',bkopal
      write(io_stdo_bgc,*) '*   wpoc_const   = ',wpoc_const*dtbinv
      write(io_stdo_bgc,*) '*   wcal_const   = ',wcal_const*dtbinv
      write(io_stdo_bgc,*) '*   wopal_const  = ',wopal_const*dtbinv
      write(io_stdo_bgc,*) '*   drempoc      = ',drempoc*dtbinv
      write(io_stdo_bgc,*) '*   dremopal     = ',dremopal*dtbinv
      write(io_stdo_bgc,*) '*   dremn2o      = ',dremn2o*dtbinv
      write(io_stdo_bgc,*) '*   dremsul      = ',dremsul*dtbinv
      write(io_stdo_bgc,*) '*   bluefix      = ',bluefix*dtbinv
      write(io_stdo_bgc,*) '*   tf0          = ',tf0
      write(io_stdo_bgc,*) '*   tf1          = ',tf1
      write(io_stdo_bgc,*) '*   tf2          = ',tf2
      write(io_stdo_bgc,*) '*   tff          = ',tff
      write(io_stdo_bgc,*) '*   ro2ut        = ',ro2ut
      write(io_stdo_bgc,*) '*   rcar         = ',rcar
      write(io_stdo_bgc,*) '*   rnit         = ',rnit
      write(io_stdo_bgc,*) '*   rnoi         = ',rnoi
      write(io_stdo_bgc,*) '*   rdnit0       = ',rdnit0
      write(io_stdo_bgc,*) '*   rdnit1       = ',rdnit1
      write(io_stdo_bgc,*) '*   rdnit2       = ',rdnit2
      write(io_stdo_bgc,*) '*   rdn2o1       = ',rdn2o1
      write(io_stdo_bgc,*) '*   rdn2o2       = ',rdn2o2
      write(io_stdo_bgc,*) '*   rcalc        = ',rcalc
      write(io_stdo_bgc,*) '*   ropal        = ',ropal
      write(io_stdo_bgc,*) '*   ctochl       = ',ctochl
      write(io_stdo_bgc,*) '*   atten_w      = ',atten_w
      write(io_stdo_bgc,*) '*   atten_c      = ',atten_c
      write(io_stdo_bgc,*) '*   atten_f      = ',atten_f
      write(io_stdo_bgc,*) '*   atten_uv     = ',atten_uv
      write(io_stdo_bgc,*) '*   fetune       = ',fetune
      write(io_stdo_bgc,*) '*   perc_diron   = ',perc_diron
      write(io_stdo_bgc,*) '*   riron        = ',riron
      write(io_stdo_bgc,*) '*   fesoly       = ',fesoly
      write(io_stdo_bgc,*) '*   relaxfe      = ',relaxfe*dtbinv
      write(io_stdo_bgc,*) '*   dmsp1        = ',dmsp1
      write(io_stdo_bgc,*) '*   dmsp2        = ',dmsp2
      write(io_stdo_bgc,*) '*   dmsp3        = ',dmsp3
      write(io_stdo_bgc,*) '*   dmsp4        = ',dmsp4
      write(io_stdo_bgc,*) '*   dmsp5        = ',dmsp5
      write(io_stdo_bgc,*) '*   dmsp6        = ',dmsp6
      if (use_BROMO) then
        write(io_stdo_bgc,*) '*   rbro         = ',rbro
        write(io_stdo_bgc,*) '*   atm_bromo    = ',atm_bromo
        write(io_stdo_bgc,*) '*   fbro1        = ',fbro1
        write(io_stdo_bgc,*) '*   fbro2        = ',fbro2
      endif
      if (use_WLIN .and. .not. use_AGG) then
        write(io_stdo_bgc,*) '*   wmin         = ',wmin*dtbinv
        write(io_stdo_bgc,*) '*   wmax         = ',wmax*dtbinv
        write(io_stdo_bgc,*) '*   wlin         = ',wlin*dtbinv
      endif
      if (.not. use_AGG) then
        write(io_stdo_bgc,*) '*   dustd1       = ',dustd1
        write(io_stdo_bgc,*) '*   dustd2       = ',dustd2
        write(io_stdo_bgc,*) '*   dustsink     = ',dustsink*dtbinv
        write(io_stdo_bgc,*) '*   wdust_const  = ',wdust_const*dtbinv
      else
        write(io_stdo_bgc,*)
        write(io_stdo_bgc,*) '********************************************'
        write(io_stdo_bgc,*) '* iHAMOCC aggregate sinking scheme:'
        write(io_stdo_bgc,*) '*   alar1      = ',alar1
        write(io_stdo_bgc,*) '*   alar2      = ',alar2
        write(io_stdo_bgc,*) '*   alar3      = ',alar3
        write(io_stdo_bgc,*) '*   alow1      = ',alow1
        write(io_stdo_bgc,*) '*   alow2      = ',alow2
        write(io_stdo_bgc,*) '*   alow3      = ',alow3
        write(io_stdo_bgc,*) '*   calmax     = ',calmax
        write(io_stdo_bgc,*) '*   cellmass   = ',cellmass
        write(io_stdo_bgc,*) '*   cellsink   = ',cellsink
        write(io_stdo_bgc,*) '*   dustd1     = ',dustd1
        write(io_stdo_bgc,*) '*   dustd2     = ',dustd2
        write(io_stdo_bgc,*) '*   dustd3     = ',dustd3
        write(io_stdo_bgc,*) '*   fractdim   = ',fractdim
        write(io_stdo_bgc,*) '*   fse        = ',fse
        write(io_stdo_bgc,*) '*   fsh        = ',fsh
        write(io_stdo_bgc,*) '*   nmldmin    = ',nmldmin
        write(io_stdo_bgc,*) '*   plower     = ',plower
        write(io_stdo_bgc,*) '*   pupper     = ',pupper
        write(io_stdo_bgc,*) '*   safe       = ',safe
        write(io_stdo_bgc,*) '*   sinkexp    = ',sinkexp
        write(io_stdo_bgc,*) '*   stick      = ',stick
        write(io_stdo_bgc,*) '*   tmfac      = ',tmfac
        write(io_stdo_bgc,*) '*   tsfac      = ',tsfac
        write(io_stdo_bgc,*) '*   vsmall     = ',vsmall
        write(io_stdo_bgc,*) '*   zdis       = ',zdis
        write(io_stdo_bgc,*) '* Maximum sinking speed for aggregates of '
        write(io_stdo_bgc,*) '* maximum size ', alar1, ' cm is '
        write(io_stdo_bgc,*)   cellsink/dtb*(alar1/alow1)**SinkExp, ' m/day'
        write(io_stdo_bgc,*) '* dust diameter (cm)', dustd1
        write(io_stdo_bgc,*) '* dust sinking speed (m/d)', dustsink / dtb
      endif
      write(io_stdo_bgc,*)
      write(io_stdo_bgc,*) '********************************************'
      write(io_stdo_bgc,*) '* Values of MO_PARAM_BGC sediment variables : '
      write(io_stdo_bgc,*) '*   sedict       = ',sedict      * dtbgcinv
      write(io_stdo_bgc,*) '*   disso_poc    = ',disso_poc   * dtbgcinv
      write(io_stdo_bgc,*) '*   disso_sil    = ',disso_sil   * dtbgcinv
      write(io_stdo_bgc,*) '*   disso_caco3  = ',disso_caco3 * dtbgcinv
      write(io_stdo_bgc,*) '*   sed_denit    = ',sed_denit   * dtbgcinv
      write(io_stdo_bgc,*) '*   silsat       = ',silsat
      write(io_stdo_bgc,*) '*   orgwei       = ',orgwei
      write(io_stdo_bgc,*) '*   opalwei      = ',opalwei
      write(io_stdo_bgc,*) '*   calcwei      = ',calcwei
      write(io_stdo_bgc,*) '*   orgdens      = ',orgdens
      write(io_stdo_bgc,*) '*   opaldens     = ',opaldens
      write(io_stdo_bgc,*) '*   calcdens     = ',calcdens
      write(io_stdo_bgc,*) '*   claydens     = ',claydens
    endif
    if (use_extNcycle) then
!    call extNwatercol_param_write()
!    call extNsediment_param_write()
    endif
  end subroutine write_parambgc

end module mo_param_bgc
