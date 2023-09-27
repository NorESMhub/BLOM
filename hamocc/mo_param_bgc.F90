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
!******************************************************************************
!
! BELEG_PARM - now mo_param_bgc - initialize bgc parameters.
!
!  Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!  Modified
!  --------
!  J.Schwinger,        *NORCE Climate, Bergen*    2020-05-19
!   -split the original BELEG_BGC in two parts, BELEG_PARM and BELEG_VARS
!  jmaerz 
!   - rename beleg_parm to mo_param_bgc 
!
!  Purpose
!  -------
!  - set bgc parameter values.
!     
!
!  Parameter list:
!  ---------------
!  *INTEGER*   *kpie*    - 1st dimension of model grid.
!  *INTEGER*   *kpje*    - 2nd dimension of model grid.
!
!******************************************************************************

  use mo_carbch,      only: atm,atm_co2,atm_n2,atm_n2o,atm_o2,dmspar
  use mo_biomod,      only: atten_c,atten_f,atten_uv,atten_w,bkopal,bkphy,bkopal,bkzoo,bluefix,ctochl,dremn2o,dremopal,            &
                          & drempoc,dremsul,dyphy,ecan,epsher,fesoly,fetune,gammap,gammaz,grami,grazra,perc_diron,phytomi,         &
                          & pi_alpha,rcalc,rcar, rdn2o1,rdn2o2,rdnit0,rdnit1,rdnit2,relaxfe,remido,riron,rnit,rnoi,ro2ut,          &
                          & ropal,spemor,tf0,tf1,tf2,tff,wcal,wdust,wopal,wpoc,zinges,drempoc_anaerob,bkox_drempoc 
  use mo_sedmnt,      only: claydens,o2ut,rno3
  use mo_control_bgc, only: io_stdo_bgc,bgc_namelist,lm4ago
  use mo_param1_bgc,  only: iatmco2,iatmnco2,iatmo2,iatmn2,iatmc13,iatmc14,iatmbromo
  use mod_xc,         only: mnproc
  use mo_m4ago,       only: init_m4ago_nml_params, init_m4ago_params
#ifdef AGG
  use mo_biomod,      only: alar1,alar2,alar3,alow1,alow2,alow3,calmax,cellmass,cellsink,dustd1,dustd2,dustd3,dustsink,            &
                          & fractdim,fse,fsh,nmldmin,plower,pupper,safe,sinkexp,stick,tmfac,tsfac,vsmall,zdis
#elif defined(WLIN)
  use mo_biomod,      only: wmin,wmax,wlin
#endif
#ifdef BROMO
  use mo_biomod,      only: rbro
  use mo_carbch,      only: atm_bromo,fbro1,fbro2
#endif
#ifdef cisonew
  use mo_biomod,      only: bifr13,bifr14,c14fac,prei13,prei14,re1312,re14to
  use mo_carbch,      only: atm_c13, atm_c14,c14_t_half,c14dec
#endif
#ifdef natDIC
  use mo_carbch,      only: atm_co2_nat
#endif
#ifdef extNcycle
      use mo_param1_bgc,  only: iatmnh3,iatmn2o
      use mo_carbch,      only: atm_nh3
      use mo_extNwatercol,only: extNwatercol_param_init,extNwatercol_param_update,extNwatercol_param_write,                        &
                                rano3denit,rano2anmx,rano2denit,ran2odenit,rdnra,ranh4nitr,rano2nitr
      use mo_extNsediment,only: extNsediment_param_init,extNsediment_param_update,extNsediment_param_write,                        &
                                rano3denit_sed,rano2anmx_sed,rano2denit_sed,ran2odenit_sed,rdnra_sed,ranh4nitr_sed,rano2nitr_sed
#endif

  implicit none

  private
   
  public :: ini_parambgc

  ! Module-wide parameters (used in more than one subroutine)
#ifndef AGG
  REAL :: dustd1, dustd2, dustsink
#endif
#ifdef cisonew 
  REAL :: beta13, alpha14, d14cat, d13c_atm  
#endif

  contains

  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine ini_parambgc(kpie,kpje)
    !
    ! First, Initialze parameters of individual components with default values. 
    ! The order of initialization can matter due to interdependcies.
    ! Then read the namelist and adjust rates to 'per time step'
    ! Re-adjust dependent parameter values
    ! Eventually write out the used parameters to the log file
    !
    implicit none      

    INTEGER, intent(in) :: kpie,kpje

    call ini_param_atm()           ! Initialize default atmospheric parameters
    call ini_stoichiometry()       ! Initialize fixed stoichiometric parameters 
    call ini_param_biol()          ! initialize biological parameters
#ifdef AGG
    call ini_aggregation()         ! Initialize aggregation module of Iris Kriest (no NML read thus far)
#endif

    call read_bgcnamelist()        ! read the BGCPARAMS namelist
    call calc_param_atm()          ! calculate atmospheric parameters after updating parameters via nml
    call ini_fields_atm(kpie,kpje) ! initialize atmospheric fields with (updated) parameter values 
    call readjust_param()          ! potentially readjust namlist parameter-dependent parameters
    call rates_2_timestep()        ! Converting rates from /d... to /dtb    
    call init_m4ago_params()       ! Initialize M4AGO parameters relying on nml parameters

    call write_parambgc()          ! write out used parameters and calculate back rates from /dtb to /d..
  end subroutine


  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine ini_param_atm()
    !
    ! Atmospheric concentrations (atm_co2 is set via BGCNML namelist).
    !
    atm_o2      = 196800.
    atm_n2      = 802000.
    atm_n2o     = 300e3   !Atmospheric mixing ratio of N2O around 1980 300 ppb,provided in ppt,300ppb = 300e3ppt = 3e-7 mol/mol

#ifdef natDIC
    atm_co2_nat = 284.32 ! CMIP6 pre-industrial reference
#endif
#ifdef BROMO
    !For now use 3.4ppt from Hense and Quack (2009; Biogeosciences) NEED TO
    !BE UPDATED WITH Ziska et al. (2013) climatology database
    atm_bromo   = 3.4
#endif
#ifdef cisonew
    ! set standard carbon isotope ratios
    re1312      = 0.0112372
    re14to      = 1.170e-12  ! Karlen et al. 1965 / Orr et al. 2017
    ! set preindustr. d13c and bigd14C in atmosphere
    prei13      = -6.5
    prei14      =  0.
#endif cisonew
#ifdef extNcycle
      ! Six & Mikolajewicz 2022: less than 1nmol m-3
      atm_nh3 = 0.
#endif
  end subroutine
  

  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine calc_param_atm()
    !
    ! AFTER having read the namelist:
    ! calculate parameters for atmosphere from given parameters 
    !
#ifdef cisonew
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
#endif
  end subroutine 

  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine ini_fields_atm(kpie,kpje)
    ! AFTER having read the nml: 
    ! Initialise atmosphere fields. We use a 2D representation of atmospheric
    ! fields for simplicity, even for cases where actually only a scalar value 
    ! is used. The overhead of this is small. If an atm-field is present in
    ! restart file (if BOXATM is activated), this will be overwritten later.

    INTEGER, intent(in) :: kpie,kpje

    ! local variables
    INTEGER :: i,j
  
    DO j=1,kpje
    DO i=1,kpie 
      atm(i,j,iatmco2)  = atm_co2
      atm(i,j,iatmo2)   = atm_o2
      atm(i,j,iatmn2)   = atm_n2
#ifdef natDIC
      atm(i,j,iatmnco2) = atm_co2_nat
#endif   
#ifdef cisonew
      atm(i,j,iatmc13)  = atm_c13
      atm(i,j,iatmc14)  = atm_c14/c14fac
#endif
#ifdef BROMO
      atm(i,j,iatmbromo)= atm_bromo
#endif
#ifdef extNcycle
        atm(i,j,iatmnh3)  = atm_nh3
        atm(i,j,iatmn2o)  = atm_n2o
#endif
    ENDDO
    ENDDO
  end subroutine

  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine ini_stoichiometry()
    !
    ! Initialize FIXED stoichiometric parameters for iHAMOCC
    !
    ! extended redfield ratio declaration
    ! Note: stoichiometric ratios are based on Takahashi etal. (1985)
    ! P:N:C:-O2 + 1:16:122:172
    ro2ut  = 172. 
    rcar   = 122.
    rnit   = 16.
    rnoi   = 1./rnit
    riron  = 5.*rcar*1.e-6       ! fe to P ratio in organic matter

    ! stoichiometric ratios for denitrification from Paulmier et al. 2009, Table 1 and
    ! equation 18. Note that their R_0=ro2ut-2*rnit.
    rdnit0 = 0.8*ro2ut           ! moles nitrate lost for remineralisation of 1 mole P
    rdnit1 = 0.8*ro2ut - rnit    ! moles nitrate net  for remineralisation of 1 mole P
    rdnit2 = 0.4*ro2ut           ! moles N2 released  for remineralisation of 1 mole P

    ! stoichiometric ratios for N2O loss by "intermediate dinitrification". Note that there
    ! is no nitrate created by this process, organic N is released as N2
    rdn2o1 = 2*ro2ut - 2.5*rnit  ! moles N2O used for remineralisation of 1 mole P
    rdn2o2 = 2*ro2ut - 2*rnit    ! moles N2 released  for remineralisation of 1 mole P

    !ik for interaction with sediment module
    o2ut   = 172.
    rno3   = 16.
  end subroutine

  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine ini_param_biol()
    !
    ! Initialize default biogeochemistry parameters.
    !
    ! Note that rates are initialized here in /d or equivalent and 
    ! time step adjustment is done after reading the BGCPARAMS namelist
    !
    !********************************************************************
    !     Phytoplankton parameters (incl. cyanobacteria)
    !********************************************************************
    phytomi = 1.e-11        !kmol/m3 - i.e. 1e-5 mmol P/m3 minimum concentration of phyto plankton (?js)
    pi_alpha= 0.02*0.4      ! initial slope of production vs irradiance curve (alpha) (0.002 for 10 steps per day)
    bkphy   = 4.e-8         !kmol/m3 - i.e. 0.04 mmol P/m3 half saturation constant
    dyphy   = 0.004         !1/d -mortality rate of phytoplankton 

    ! N2-Fixation following the parameterization in Kriest and Oschlies, 2015.
    ! Factors tf2, tf1 and tf0 are a polynomial (2nd order) 
    ! approximation to the functional relationship by Breitbarth et al. (2007),
    ! for temperature dependence of Trichodesmium growth, their eq. (2)
    ! The relation will be scaled to their max. growth rate, tff.
    ! Note that the second order approx. is basically similar to their
    ! function 2 for T-dependent nitrogen fixation multiplied by 4 
    ! (2 [N atoms per mole] * 12 [light hrs per day]/6 [C-atoms per N-atoms])
    bluefix =  0.005         !1/d  ! nitrogen fixation rate by blue green algae (cyanobacteria)
    tf2     = -0.0042
    tf1     =  0.2253
    tf0     = -2.7819
    tff     =  0.2395  

#ifdef cisonew
    ! Initial fractionation during photosynthesis
    bifr13     = 0.98
    bifr14     = bifr13**2
    ! Decay parameter for sco214, HalfLive = 5730 years
    c14_t_half = 5700.*365. ! Half life of 14C [days]	
#endif
#ifdef BROMO
    !Bromoform to phosphate ratio (Hense and Quack, 2009)
    !JT: too little production: 0.25Gmol/yr     rbro=6.72e-7*rnit
    !      rbro=2.*6.72e-7*rnit
    !JT Following discussion with B. Quack and D. Booge (01.07.2021), we agree to use 2.4e-6 
    rbro  = 2.4e-6*rnit
    fbro1 = 1.0
    fbro2 = 1.0
#endif

    !********************************************************************
    !     Zooplankton parameters
    !********************************************************************
    grami   = 1.e-10        !kmol/m3 - i.e. 1e-5 mmol P/m3 minimum concentration of zooplankton
    bkzoo   = 8.e-8         !kmol/m3 - i.e. 0.08 mmol P/m3 half saturation constant

    !ik addded parameter definition; taken from OCPROD.F
    grazra  = 1.2           !1/d -grazing rate
    spemor  = 3.*1.e6       !1/d -mortality rate
    gammap  = 0.04          !1/d -exudation rate
    gammaz  = 0.06          !1/d -excretion rate
    ecan    = 0.95          ! fraction of mortality as PO_4
#ifdef AGG
    zinges  = 0.5           !dimensionless fraction -assimilation efficiency
    epsher  = 0.9           !dimensionless fraction -fraction of grazing egested
#elif defined(WLIN)
    zinges  = 0.7           !dimensionless fraction -assimilation efficiency
    epsher  = 0.85          !dimensionless fraction -fraction of grazing egested
#else
    zinges  = 0.6           !dimensionless fraction -assimilation efficiency
    epsher  = 0.8           !dimensionless fraction -fraction of grazing egest      
#endif

    !********************************************************************
    !     Shell production (CaCO3 and opal) parameters
    !******************************************************************** 
    bkopal = 5.e-6     !kmol/m3 - i.e. 4.0 mmol Si/m3 half saturation constant
#ifdef AGG
    rcalc  = 14.  ! calcium carbonate to organic phosphorous production ratio
    ropal  = 10.5 ! opal to organic phosphorous production ratio      
    calmax = 0.20
#elif defined(WLIN)
    rcalc  = 33.  ! calcium carbonate to organic phosphorous production ratio
    ropal  = 45.  ! opal to organic phosphorous production ratio      
#else
    rcalc  = 40.  ! iris 40 !calcium carbonate to organic phosphorous production ratio
    ropal  = 30.  ! iris 25 !opal to organic phosphorous production ratio      
#endif
    
    !********************************************************************
    !     Remineralization and dissolution parameters (incl. DMS prod.)
    !********************************************************************
    remido   = 0.004    !1/d -remineralization rate (of DOM)
    ! deep sea remineralisation constants
    drempoc  = 0.025    !1/d Aerob remineralization rate detritus
    drempoc_anaerob = 0.05*drempoc ! remin in sub-/anoxic environm. - not be overwritten by lm4ago
    bkox_drempoc = 1e-7 ! half-saturation constant for oxygen for ammonification (aerobic remin via drempoc)
    dremopal = 0.003    !1/d Dissolution rate for opal 
    dremn2o  = 0.01     !1/d Remineralization rate of detritus on N2O
    dremsul  = 0.005    !1/d Remineralization rate for sulphate reduction 
    if(lm4ago)then
        ! reset drempoc and dremopal for T-dep remin/dissolution
        drempoc  = 0.12
        dremopal = 0.023
    endif

    ! M4AGO parameters
    call init_m4ago_nml_params()
#ifdef extNcycle
    ! initialize the extended nitrogen cycle parameters - first water column, then sediment, 
    ! since sediment relies on water column parameters for the extended nitrogen cycle 
    ! Sediment also relies on M4AGO being initialized (POM_remin_q10 and POM_remin_Tref)
    call extNwatercol_param_init()
    call extNsediment_param_init()
#endif   
    ! Set constants for calculation of dms ( mo_carbch )
    ! Parameters are a result from kettle optimisation 02.03.04
    dmspar(6)=0.100000000E-07  !0 half saturation microbial
    dmspar(5)=1.25*0.02  ! production with delsil, following Kloster et al., 06 Table 1, but increased by a factor of ~2
    dmspar(4)=1.25*0.10  ! production with delcar, following Kloster et al., 06 Table 1, but reduced by ~7%
    dmspar(3)=0.0864     ! following Kloster et al., 06 Table 1 with 50% reduction to reduce bacterial removal and increase dms emissions
    dmspar(2)=0.0011     ! following Kloster et al., 06 Table 1
    dmspar(1)=10.        ! 2*5. production with temp


    !********************************************************************
    !     Radiation and attenuation parameters
    !********************************************************************
    ! parameters for sw-radiation attenuation
    ! Analog to Moore et al., Deep-Sea Research II 49 (2002), 403-462
    ! 1 kmolP = (122*12/60)*10^6 mg[Chlorophyl] 
    ctochl  = 60.        ! C to Chlorophyl ratio
    atten_w = 0.04       ! yellow substances attenuation in 1/m
    atten_c = 0.03*rcar*(12./ctochl)*1.e6  ! phytoplankton attenuation in 1/m 
    atten_uv= 0.33       ! 
    atten_f = 0.4        ! fraction of sw-radiation directly absorbed in surface layer 
                         ! (only if FB_BGC_OCE) [feedback bgc-ocean]

    !********************************************************************
    !     Dust deposition and iron solubility parameters
    !********************************************************************
    !ik weight percent iron in dust deposition times Fe solubility
    ! the latter three values come from Johnson et al., 1997
    fetune     = 0.6                  ! factor introduced to tune deposition/solubility
    perc_diron = fetune * 0.035 * 0.01 / 55.85
    fesoly     = 0.5*1.e-9            ! max. diss. iron concentration in deep water 
    relaxfe    = 0.05/365.            ! 1/d complexation rate to relax iron concentration to fesoly

    !********************************************************************
    !     Sinking parameters
    !********************************************************************

    wpoc   =  5.       !m/d   Sinking speed of detritus iris : 5.
    wcal   = 30.       !m/d   Sinking speed of CaCO3 shell material
    wopal  = 30.       !m/d   Sinking speed of opal iris : 60
#if defined(WLIN) && ! defined(AGG)
    wmin   =  1.       !m/d   minimum sinking speed
    wmax   = 60.       !m/d   maximum sinking speed
    wlin   = 60./2400. !m/d/m constant describing incr. with depth, r/a=1.0
#endif
#ifndef AGG
    dustd1   = 0.0001 !cm = 1 um, boundary between clay and silt
    dustd2   = dustd1*dustd1
    dustsink = (9.81 * 86400. / 18.                    &  ! g * sec per day / 18.
     &         * (claydens - 1025.) / 1.567 * 1000.    &  !excess density / dyn. visc.
     &         * dustd2 * 1.e-4)
    wdust = dustsink
#endif    
  end subroutine

  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine read_bgcnamelist()
    !
    ! Read the bgcparams namelist for parameter tuning.
    ! Note that afterward, i) rates need to be adjusted for timestep
    ! and some depending parameters need re-calculation
    !

    integer  :: iounit

    namelist /bgcparams/ bkphy,dyphy,bluefix,bkzoo,grazra,spemor,gammap,gammaz,ecan,zinges,epsher,bkopal,rcalc,ropal,              &
                       & remido,drempoc,dremopal,dremn2o,dremsul,fetune,relaxfe,wpoc,                                              &
#if defined(WLIN) && ! defined(AGG)
                       & wmin,wmax,wlin,                                                                                           &
#endif 
#ifdef extNcycle
                       & rano3denit,rano2anmx,rano2denit,ran2odenit,rdnra,ranh4nitr,rano2nitr,                                     &
                       & rano3denit_sed,rano2anmx_sed,rano2denit_sed,ran2odenit_sed,rdnra_sed,ranh4nitr_sed,rano2nitr_sed,         &
                       & atm_nh3,                                                                                                  &
#endif
                       & wcal,wopal,atm_n2o

    open (newunit=iounit, file=bgc_namelist, status='old',action='read')
    read (unit=iounit, nml=BGCPARAMS)
    close (unit=iounit)

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*) '------------------------------------------'
      write(io_stdo_bgc,*) 'iHAMOCC: read namelist bgcparams'
      write(io_stdo_bgc,nml=BGCPARAMS)
      write(io_stdo_bgc,*) '------------------------------------------'
    endif

  end subroutine 
 
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine readjust_param()
    !
    !  AFTER reading the namelist:
    ! re-adjust parameters that depend on tuning parameters
    !
    perc_diron = fetune * 0.035 * 0.01 / 55.85

  end subroutine

  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine rates_2_timestep()
    !
    ! AFTER potential update of rates, convert them to rates per timestep
    !
    use mo_control_bgc, only: dtb

!    implicit none

    !********************************************************************
    !     Phytoplankton parameters (incl. cyanobacteria)
    !********************************************************************
    dyphy  = dyphy*dtb       !1/d -mortality rate of phytoplankton 

    ! nitrogen fixation by blue green algae
    bluefix = bluefix*dtb     !1/d

#ifdef cisonew
    c14dec = 1.-(log(2.)/c14_t_half)*dtb   ! lambda [1/day]; c14dec[-]
#endif

    !********************************************************************
    !     Zooplankton parameters
    !********************************************************************
    grazra = grazra*dtb      !1/d -grazing rate
    spemor = spemor*dtb      !1/d -mortality rate
    gammap = gammap*dtb      !1/d -exudation rate
    gammaz = gammaz*dtb      !1/d -excretion rate

    !********************************************************************
    !     Remineralization and dissolution parameters 
    !********************************************************************
    remido   = remido*dtb     !1/d -remineralization rate (of DOM)
    ! deep sea remineralisation constants
    drempoc  = drempoc*dtb    !1/d  Aerob remineralization rate of detritus
    drempoc_anaerob = drempoc_anaerob*dtb ! 1/d Anaerob remin rate of detritus
    dremopal = dremopal*dtb   !1/d  Dissolution rate of opal 
    dremn2o  = dremn2o*dtb    !1/d  Remineralization rate of detritus on N2O
    dremsul  = dremsul*dtb    !1/d  Remineralization rate for sulphate reduction 

    !********************************************************************
    !     Dust deposition and iron solubility parameters
    !********************************************************************
    relaxfe = relaxfe*dtb     !1/d  iron complexation rate

    
    !********************************************************************
    !     Sinking parameters
    !********************************************************************
    wpoc  = wpoc*dtb       !m/d  Sinking speed detritusiris : 5.
    wcal  = wcal*dtb       !m/d  Sinking speed CaCO3
    wopal = wopal*dtb      !m/d  Sinking speed opal iris : 60
#if defined(WLIN) && ! defined(AGG)
    wmin  = wmin*dtb       !m/d   minimum sinking speed
    wmax  = wmax*dtb       !m/d   maximum sinking speed
    wlin  = wlin*dtb       !m/d/m constant describing incr. with depth, r/a=1.0
#endif
#ifndef AGG
    wdust = wdust*dtb      !m/d   dust sinking speed
#endif 
#ifdef extNcycle
    call extNwatercol_param_update()
    call extNsediment_param_update()
#endif

  end subroutine

  !---------------------------------------------------------------------------------------------------------------------------------
#ifdef AGG
  subroutine ini_aggregation()
    !
    ! parameters needed for the aggregation module
    !
    use mo_control_bgc, only: dtb

    REAL :: shear

    SinkExp  = 0.62
    FractDim = 1.62
    !      Stick = 0.40
    !      Stick = 0.25
    Stick    = 0.15
    cellmass = 0.012 / rnit ![nmol P]
    !ik      cellmass = 0.0039/ rnit ![nmol P] for a 10 um diameter
    cellsink = 1.40 *dtb! [m/d]
    !ik      cellsink = 0.911 *dtb! [m/d]  for a 10 um diameter
    !      shear = 86400. !shear in the mixed layer,   1.0  d-1
    !      shear = 64800. !shear in the mixed layer,   0.75 d-1
    shear   = 43200. !shear in the mixed layer,   0.5  d-1
    fsh     = 0.163 * shear *dtb
    fse     = 0.125 * 3.1415927 * cellsink * 100.
    alow1   = 0.002 !diameter of smallest particle [cm]
    !ik      alow1 = 0.001 !diameter of smallest particle [cm]
    alow2   = alow1 * alow1
    alow3   = alow2 * alow1
    !      alar1 = 1.0 !diameter of largest particle for size dependend aggregation and sinking [cm]
    !      alar1 = 0.75 !diameter of largest particle for size dependend aggregation and sinking [cm]
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

    ! for shear aggregation of dust:
    dustd1  = 0.0001 !cm = 1 um, boundary between clay and silt
    dustd2  = dustd1*dustd1
    dustd3  = dustd2*dustd1
    dustsink = (9.81 * 86400. / 18.                & ! g * sec per day / 18.                 
     &         * (claydens - 1025.) / 1.567 * 1000.  & !excess density / dyn. visc.
     &         * dustd2 * 1.e-4)*dtb
    if(dustsink.gt.cellsink) then 
       if (mnproc.eq.1)then
         write(io_stdo_bgc,*) ' dust sinking speed greater than cellsink'
         write(io_stdo_bgc,*) ' set dust sinking speed to cellsink'
       endif
       dustsink = cellsink
    endif
  end subroutine 
#endif

  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine write_parambgc()
    !
    ! Write parameters
    !
    use mo_control_bgc, only: dtb

    REAL :: dtbinv
    dtbinv = 1./dtb

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*) '****************************************************************'
      WRITE(io_stdo_bgc,*) '* '
      WRITE(io_stdo_bgc,*) '* Values of MO_PARAM_BGC variables : '
      WRITE(io_stdo_bgc,*) '*          atm_co2      = ',atm_co2      
#ifdef cisonew
      WRITE(io_stdo_bgc,*) '*          atm_c13      = ',atm_c13      
      WRITE(io_stdo_bgc,*) '*          d13C_atm     = ',d13C_atm    
      WRITE(io_stdo_bgc,*) '*          atm_c14      = ',atm_c14  
      WRITE(io_stdo_bgc,*) '*          bifr13       = ',bifr13 
      WRITE(io_stdo_bgc,*) '*          bifr14       = ',bifr14
      WRITE(io_stdo_bgc,*) '*          c14fac       = ',c14fac
      WRITE(io_stdo_bgc,*) '*          prei13       = ',prei13
      WRITE(io_stdo_bgc,*) '*          prei14       = ',prei14
      WRITE(io_stdo_bgc,*) '*          re1312       = ',re1312
      WRITE(io_stdo_bgc,*) '*          re14to       = ',re14to
      WRITE(io_stdo_bgc,*) '*          c14_t_half   = ',c14_t_half
      WRITE(io_stdo_bgc,*) '*          c14dec       = ',c14dec
      WRITE(io_stdo_bgc,*) '*          beta13       = ',beta13
      WRITE(io_stdo_bgc,*) '*          alpha14      = ',alpha14
      WRITE(io_stdo_bgc,*) '*          d14cat       = ',d14cat
      WRITE(io_stdo_bgc,*) '*          c14fac       = ',c14fac
#endif
      WRITE(io_stdo_bgc,*) '*          atm_o2       = ',atm_o2           
      WRITE(io_stdo_bgc,*) '*          atm_n2       = ',atm_n2
#ifdef extNcycle
      WRITE(io_stdo_bgc,*) '*          atm_nh3      = ',atm_nh3
#endif
      WRITE(io_stdo_bgc,*) '*          atm_n2o      = ',atm_n2o
      WRITE(io_stdo_bgc,*) '*          phytomi      = ',phytomi
      WRITE(io_stdo_bgc,*) '*          grami        = ',grami
      WRITE(io_stdo_bgc,*) '*          remido       = ',remido*dtbinv
      WRITE(io_stdo_bgc,*) '*          dyphy        = ',dyphy*dtbinv
      WRITE(io_stdo_bgc,*) '*          zinges       = ',zinges
      WRITE(io_stdo_bgc,*) '*          epsher       = ',epsher
      WRITE(io_stdo_bgc,*) '*          grazra       = ',grazra*dtbinv
      WRITE(io_stdo_bgc,*) '*          spemor       = ',spemor*dtbinv
      WRITE(io_stdo_bgc,*) '*          gammap       = ',gammap*dtbinv
      WRITE(io_stdo_bgc,*) '*          gammaz       = ',gammaz*dtbinv
      WRITE(io_stdo_bgc,*) '*          ecan         = ',ecan    
      WRITE(io_stdo_bgc,*) '*          pi_alpha     = ',pi_alpha
      WRITE(io_stdo_bgc,*) '*          bkphy        = ',bkphy
      WRITE(io_stdo_bgc,*) '*          bkzoo        = ',bkzoo    
      WRITE(io_stdo_bgc,*) '*          bkopal       = ',bkopal    
      WRITE(io_stdo_bgc,*) '*          wpoc         = ',wpoc*dtbinv
      WRITE(io_stdo_bgc,*) '*          wcal         = ',wcal*dtbinv    
      WRITE(io_stdo_bgc,*) '*          wopal        = ',wopal*dtbinv   
      WRITE(io_stdo_bgc,*) '*          drempoc      = ',drempoc*dtbinv    
      WRITE(io_stdo_bgc,*) '*          dremopal     = ',dremopal*dtbinv   
      WRITE(io_stdo_bgc,*) '*          dremn2o      = ',dremn2o*dtbinv   
      WRITE(io_stdo_bgc,*) '*          dremsul      = ',dremsul*dtbinv   
      WRITE(io_stdo_bgc,*) '*          bluefix      = ',bluefix*dtbinv   
      WRITE(io_stdo_bgc,*) '*          tf0          = ',tf0   
      WRITE(io_stdo_bgc,*) '*          tf1          = ',tf1   
      WRITE(io_stdo_bgc,*) '*          tf2          = ',tf2   
      WRITE(io_stdo_bgc,*) '*          tff          = ',tff   
      WRITE(io_stdo_bgc,*) '*          ro2ut        = ',ro2ut   
      WRITE(io_stdo_bgc,*) '*          rcar         = ',rcar 
      WRITE(io_stdo_bgc,*) '*          rnit         = ',rnit
      WRITE(io_stdo_bgc,*) '*          rnoi         = ',rnoi
      WRITE(io_stdo_bgc,*) '*          rdnit0       = ',rdnit0
      WRITE(io_stdo_bgc,*) '*          rdnit1       = ',rdnit1
      WRITE(io_stdo_bgc,*) '*          rdnit2       = ',rdnit2
      WRITE(io_stdo_bgc,*) '*          rdn2o1       = ',rdn2o1
      WRITE(io_stdo_bgc,*) '*          rdn2o2       = ',rdn2o2
      WRITE(io_stdo_bgc,*) '*          rcalc        = ',rcalc
      WRITE(io_stdo_bgc,*) '*          ropal        = ',ropal
      WRITE(io_stdo_bgc,*) '*          ctochl       = ',ctochl
      WRITE(io_stdo_bgc,*) '*          atten_w      = ',atten_w
      WRITE(io_stdo_bgc,*) '*          atten_c      = ',atten_c
      WRITE(io_stdo_bgc,*) '*          atten_f      = ',atten_f
      WRITE(io_stdo_bgc,*) '*          atten_uv     = ',atten_uv
      WRITE(io_stdo_bgc,*) '*          o2ut         = ',o2ut
      WRITE(io_stdo_bgc,*) '*          rno3         = ',rno3
      WRITE(io_stdo_bgc,*) '*          fetune       = ',fetune
      WRITE(io_stdo_bgc,*) '*          perc_diron   = ',perc_diron
      WRITE(io_stdo_bgc,*) '*          riron        = ',riron
      WRITE(io_stdo_bgc,*) '*          fesoly       = ',fesoly
      WRITE(io_stdo_bgc,*) '*          relaxfe      = ',relaxfe*dtbinv
      WRITE(io_stdo_bgc,*) '*          dmspar(1)    = ',dmspar(1)
      WRITE(io_stdo_bgc,*) '*          dmspar(2)    = ',dmspar(2)
      WRITE(io_stdo_bgc,*) '*          dmspar(3)    = ',dmspar(3)
      WRITE(io_stdo_bgc,*) '*          dmspar(4)    = ',dmspar(4)
      WRITE(io_stdo_bgc,*) '*          dmspar(5)    = ',dmspar(5)
#ifdef BROMO
      WRITE(io_stdo_bgc,*) '*          rbro         = ',rbro
      WRITE(io_stdo_bgc,*) '*          atm_bromo    = ',atm_bromo
      WRITE(io_stdo_bgc,*) '*          fbro1        = ',fbro1
      WRITE(io_stdo_bgc,*) '*          fbro2        = ',fbro2
#endif
#if defined(WLIN) && ! defined(AGG)
      WRITE(io_stdo_bgc,*) '*          wmin         = ',wmin
      WRITE(io_stdo_bgc,*) '*          wmax         = ',wmax
      WRITE(io_stdo_bgc,*) '*          wlin         = ',wlin
#endif 
#ifndef AGG
      WRITE(io_stdo_bgc,*) '*          dustd1       = ',dustd1
      WRITE(io_stdo_bgc,*) '*          dustd2       = ',dustd2
      WRITE(io_stdo_bgc,*) '*          dustsink     = ',dustsink
      WRITE(io_stdo_bgc,*) '*          wdust        = ',wdust*dtbinv
#else
      write(io_stdo_bgc,*)
      write(io_stdo_bgc,*) '****************************************************************'
      write(io_stdo_bgc,*) 'HAMOCC aggregate sinking scheme:' 
      write(io_stdo_bgc,*) '        alar1      = ',alar1
      write(io_stdo_bgc,*) '        alar2      = ',alar2
      write(io_stdo_bgc,*) '        alar3      = ',alar3
      write(io_stdo_bgc,*) '        alow1      = ',alow1
      write(io_stdo_bgc,*) '        alow2      = ',alow2
      write(io_stdo_bgc,*) '        alow3      = ',alow3
      write(io_stdo_bgc,*) '        calmax     = ',calmax
      write(io_stdo_bgc,*) '        cellmass   = ',cellmass
      write(io_stdo_bgc,*) '        cellsink   = ',cellsink
      write(io_stdo_bgc,*) '        dustd1     = ',dustd1
      write(io_stdo_bgc,*) '        dustd2     = ',dustd2
      write(io_stdo_bgc,*) '        dustd3     = ',dustd3
      write(io_stdo_bgc,*) '        fractdim   = ',fractdim
      write(io_stdo_bgc,*) '        fse        = ',fse
      write(io_stdo_bgc,*) '        fsh        = ',fsh
      write(io_stdo_bgc,*) '        nmldmin    = ',nmldmin
      write(io_stdo_bgc,*) '        plower     = ',plower
      write(io_stdo_bgc,*) '        pupper     = ',pupper
      write(io_stdo_bgc,*) '        safe       = ',safe
      write(io_stdo_bgc,*) '        sinkexp    = ',sinkexp
      write(io_stdo_bgc,*) '        stick      = ',stick
      write(io_stdo_bgc,*) '        tmfac      = ',tmfac
      write(io_stdo_bgc,*) '        tsfac      = ',tsfac
      write(io_stdo_bgc,*) '        vsmall     = ',vsmall
      write(io_stdo_bgc,*) '        zdis       = ',zdis
      write(io_stdo_bgc,*) ' Maximum sinking speed for aggregates of '
      write(io_stdo_bgc,*) ' maximum size ', alar1, ' cm is '
      write(io_stdo_bgc,*)   cellsink/dtb*(alar1/alow1)**SinkExp, ' m/day'
      write(io_stdo_bgc,*) ' dust diameter (cm)', dustd1
      write(io_stdo_bgc,*) ' dust sinking speed (m/d)', dustsink / dtb
      write(io_stdo_bgc,*) '****************************************************************'
#endif 
      WRITE(io_stdo_bgc,*) '*          claydens     = ',claydens
#ifdef extNcycle
    call extNwatercol_param_write()
    call extNsediment_param_write()
#endif
      WRITE(io_stdo_bgc,*) '****************************************************************'
      ENDIF
  end subroutine
end module mo_param_bgc
