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


      SUBROUTINE BELEG_PARM(kpie,kpje)
!******************************************************************************
!
! BELEG_PARM - initialize bgc parameters.
!
!  Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!  Modified
!  --------
!  J.Schwinger,        *NORCE Climate, Bergen*    2020-05-19
!   -split the original BELEG_BGC in two parts, BELEG_PARM and BELEG_VARS
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
      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt,      only: claydens,o2ut,rno3
      USE mo_control_bgc, only: dtb,io_stdo_bgc
      use mo_param1_bgc,  only: iatmco2,iatmnco2,iatmo2,iatmn2
      USE mod_xc,         only: mnproc

      implicit none      

      INTEGER, intent(in) :: kpie,kpje

      ! local variables
      INTEGER :: i,j
#ifdef cisonew
      REAL :: rco213,rco214,alpha14,beta13,beta14,d13C_atm,d14cat
#endif
#ifdef AGG
      REAL :: shear,snow
#else
      REAL :: dustd1, dustd2, dustsink
#endif


!
! Atmospheric concentrations (atm_co2 is set via namelist).
!
      atm_o2  = 196800.
      atm_n2  = 802000.
#ifdef natDIC
      atm_co2_nat = 284.32 ! CMIP6 pre-industrial reference
#endif

#ifdef cisonew
! set standard carbon isotope ratios
      re1312=0.0112372
      re14to=1.176e-12
! set preindustr. d13c and bigd14C in atmosphere
      prei13  = -6.5
      prei14  = 0.
      beta13  = (prei13/1000.)+1.
      alpha14 = 2.*(prei13+25.)
      d14cat  = (prei14+alpha14)/(1.-alpha14/1000.)
! calculate atm_c13 and atm_c14
      atm_c13  = beta13*re1312*atm_co2/(1.+beta13*re1312)
      d13C_atm = (((atm_c13/(atm_co2-atm_c13))/re1312)-1.)*1000.
! absolute 14c concentration in preindustrial atmosphere
      atm_c14  = ((d14cat/1000.)+1.)*re14to*atm_co2
! factor for normalizing 14C tracers (~1e-12)
      c14fac   = atm_c14/atm_co2
#endif

! Initialise atmosphere fields. We use a 2D representation of atmospheric
! fields for simplicity, even for cases where actually only a scalar value 
! is used. The overhead of this is small. If an atm-field is present in
! restart file (if BOXATM is activated), this will be overwritten later.
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
      ENDDO
      ENDDO


!
! Biology
!
!ik note that the unit is kmol/m3!
      phytomi=1.e-11        !i.e. 1e-5 mmol P/m3 minimum concentration of phyto plankton (?js)
      grami=1.e-10          !i.e. 1e-5 mmol P/m3 minimum concentration of zooplankton

!ik addded parameter definition; taken from OCPROD.F
      remido=0.004*dtb      !1/d -remineralization rate (of DOM)
      dyphy=0.004*dtb       !1/d -mortality rate of phytoplankton 
      grazra=1.2*dtb        !1/d -grazing rate
      spemor=3.*1.e6*dtb    !1/d -mortality rate
      gammap=0.04*dtb       !1/d -exudation rate
      gammaz=0.06*dtb       !1/d -excretion rate
      ecan=0.95             ! fraction of mortality as PO_4
      pi_alpha=0.02*0.4     ! initial slope of production vs irradiance curve (alpha) (0.002 for 10 steps per day)
#ifdef AGG
      zinges = 0.5          !dimensionless fraction -assimilation efficiency
      epsher = 0.9          !dimensionless fraction -fraction of grazing egested
#elif defined(WLIN)
      zinges = 0.7          !dimensionless fraction -assimilation efficiency
      epsher = 0.85         !dimensionless fraction -fraction of grazing egested
#else
      zinges = 0.6          !dimensionless fraction -assimilation efficiency
      epsher = 0.8          !dimensionless fraction -fraction of grazing egest      
#endif

#ifdef cisonew
! Initial fractionation during photosynthesis
      bifr13=0.98
      bifr14=bifr13**2
! Decay parameter for sco214, HalfLive = 5730 years
      c14_t_half=5730.*365.                ! Half life of 14C [days]	
      c14dec=1.-(log(2.)/c14_t_half)*dtb   ! labda [1/day]; c14dec[-]
#endif

! half sat. constants, note that the units are kmol/m3 !
      bkphy  = 4.e-8    !i.e. 0.04 mmol P/m3
      bkzoo  = 8.e-8    !i.e. 0.08 mmol P/m3
      bkopal = 5.e-6    !i.e. 4.0 mmol Si/m3

!sinking speed
      wpoc  =  5.*dtb       !m/d  iris : 5.
      wcal  = 30.*dtb       !m/d 
      wopal = 30.*dtb       !m/d  iris : 60
#ifdef WLIN
      wmin  =  1.*dtb       !m/d   minimum sinking speed
      wmax  = 60.*dtb       !m/d   maximum sinking speed
      wlin  = 60./2400.*dtb !m/d/m constant describing incr. with depth, r/a=1.0
#endif

! deep see remineralisation constants
      drempoc  = 0.025*dtb    !1/d
      dremopal = 0.003*dtb    !1/d
      dremn2o  = 0.01*dtb     !1/d
      dremsul  = 0.005*dtb    ! remineralization rate for sulphate reduction 
      

! nirogen fixation by blue green algae
      bluefix=0.005*dtb       !1/d

! N2-Fixation following the parameterization in Kriest and Oschlies, 2015.
! Factors tf2, tf1 and tf0 are a polynomial (2nd order) 
! approximation to the functional relationship by Breitbarth et al. (2007),
! for temperature dependence of Trichodesmium growth, their eq. (2)
! The relation will be scaled to their max. growth rate, tff.
! Note that the second order approx. is basically similar to their
! function 2 for T-dependent nitrogen fixation multiplied by 4 
! (2 [N atoms per mole] * 12 [light hrs per day]/6 [C-atoms per N-atoms])
      tf2 = -0.0042
      tf1 = 0.2253
      tf0 = -2.7819
      tff = 0.2395  

! extended redfield ratio declaration
! Note: stoichiometric ratios are based on Takahashi etal. (1985)
! P:N:C:-O2 + 1:16:122:172
      ro2ut=172. 
      rcar=122.
      rnit=16.
      rnoi=1./rnit

! stoichiometric ratios for denitrification from Paulmier et al. 2009, Table 1 and
! equation 18. Note that their R_0=ro2ut-2*rnit.
      rdnit0=0.8*ro2ut           ! moles nitrate lost for remineralisation of 1 mole P
      rdnit1=0.8*ro2ut-rnit      ! moles nitrate net  for remineralisation of 1 mole P
      rdnit2=0.4*ro2ut           ! moles N2 released  for remineralisation of 1 mole P

! stoichiometric ratios for N2O loss by "intermediate dinitrification". Note that there
! is no nitrate created by this process, organic N is released as N2
      rdn2o1=2*ro2ut-2.5*rnit    ! moles N2O used for remineralisation of 1 mole P
      rdn2o2=2*ro2ut-2*rnit      ! moles N2 released  for remineralisation of 1 mole P


#ifdef AGG
      rcalc = 14.  ! calcium carbonate to organic phosphorous production ratio
      ropal = 10.5 ! opal to organic phosphorous production ratio      
      calmax= 0.20
#elif defined(WLIN)
      rcalc = 33.  ! calcium carbonate to organic phosphorous production ratio
      ropal = 45.  ! opal to organic phosphorous production ratio      
#else
      rcalc = 40.  ! iris 40 !calcium carbonate to organic phosphorous production ratio
      ropal = 30.  ! iris 25 !opal to organic phosphorous production ratio      
#endif

! parameters for sw-radiation attenuation
! Analog to Moore et al., Deep-Sea Research II 49 (2002), 403-462
! 1 kmolP = (122*12/60)*10^6 mg[Chlorophyl] 

      ctochl  = 60.        ! C to Chlorophyl ratio
      atten_w = 0.04       ! yellow substances attenuation in 1/m
      atten_c = 0.03*rcar*(12./ctochl)*1.e6  ! phytoplankton attenuation in 1/m 
      atten_f = 0.4        ! fraction of sw-radiation directly absorbed in surface layer 
                           ! (only if FB_BGC_OCE) [feedback bgc-ocean]
      	
!ik for interaction with sediment module
      o2ut=172.
      rno3=16.

!ik weight percent iron in dust deposition times Fe solubility
! the latter three values come from Johnson et al., 1997
      fetune=0.6                  ! factor introduced to tune deposistion/solubility
      perc_diron = fetune * 0.035 * 0.01 / 55.85
      riron= 5.*rcar*1.e-6        ! fe to P ratio in organic matter
      fesoly=0.5*1.e-9            ! max. diss. iron concentration in deep water 
      relaxfe = 0.05/365.*dtb

!                        
! Set constants for calculation of dms ( mo_carbch )
! Parameters are a result from kettle optimisation 02.03.04

       dmspar(6)=0.100000000E-07  !0 half saturation microbial
       dmspar(5)=1.25*0.02  ! production with delsil, following Kloster et al., 06 Table 1, but increased by a factor of ~2
       dmspar(4)=1.25*0.10  ! production with delcar, following Kloster et al., 06 Table 1, but reduced by ~7%
       dmspar(3)=0.0864     ! following Kloster et al., 06 Table 1 with 50% reduction to reduce bacterial removal and increase dms emissions
       dmspar(2)=0.0011     ! following Kloster et al., 06 Table 1
       dmspar(1)=10.        ! 2*5. production with temp


      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)                                             &
     &'****************************************************************'
      WRITE(io_stdo_bgc,*)                                             &
     &'* '
      WRITE(io_stdo_bgc,*)                                             &
     &'* Values of BELEG_BGC variables : '
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atm_co2      = ',atm_co2      
#ifdef cisonew
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atm_c13      = ',atm_c13      
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              d13C_atm     = ',d13C_atm    
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atm_c14      = ',atm_c14      
#endif
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atm_o2       = ',atm_o2           
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atm_n2       = ',atm_n2 
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              phytomi      = ',phytomi
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              grami        = ',grami
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              remido       = ',remido
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dyphy        = ',dyphy
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              zinges       = ',zinges
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              epsher       = ',epsher
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              grazra       = ',grazra
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              spemor       = ',spemor
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              gammap       = ',gammap
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              gammaz       = ',gammaz
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              ecan         = ',ecan    
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              bkphy        = ',bkphy
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              bkzoo        = ',bkzoo    
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              bkopal       = ',bkopal    
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              wpoc         = ',wpoc
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              wcal         = ',wcal    
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              wopal        = ',wopal   
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              drempoc      = ',drempoc    
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dremopal     = ',dremopal   
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              bluefix      = ',bluefix   
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              ro2ut        = ',ro2ut   
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rcar         = ',rcar 
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rnit         = ',rnit
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rnoi         = ',rnoi
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rdnit1       = ',rdnit1
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rdnit2       = ',rdnit2
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rcalc        = ',rcalc
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              ropal        = ',ropal
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              ctochl       = ',ctochl
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atten_w      = ',atten_w
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atten_c      = ',atten_c
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atten_f      = ',atten_f
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              o2ut         = ',o2ut
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rno3         = ',rno3
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              perc_diron   = ',perc_diron
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              riron        = ',riron
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              fesoly       = ',fesoly
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              relaxfe      = ',relaxfe
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dmspar(1)    = ',dmspar(1)
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dmspar(2)    = ',dmspar(2)
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dmspar(3)    = ',dmspar(3)
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dmspar(4)    = ',dmspar(4)
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dmspar(5)    = ',dmspar(5)
      ENDIF

#ifndef AGG
      dustd1 = 0.0001 !cm = 1 um, boundary between clay and silt
      dustd2=dustd1*dustd1
      dustsink = (9.81 * 86400. / 18.                  &  ! g * sec per day / 18.
     &         * (claydens - 1025.) / 1.567 * 1000.    &  !excess density / dyn. visc.
     &         * dustd2 * 1.e-4)*dtb
      wdust = dustsink

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dustd1       = ',dustd1
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dustd2       = ',dustd2
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dustsink     = ',dustsink
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              wdust        = ',wdust
      ENDIF
#endif

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)                                             &
     &'****************************************************************'
      ENDIF

#ifdef AGG
! parameters needed for the aggregation module

      SinkExp = 0.62
      FractDim = 1.62
!      Stick = 0.40
!      Stick = 0.25
      Stick = 0.15
      cellmass = 0.012 / rnit ![nmol P]
!ik      cellmass = 0.0039/ rnit ![nmol P] for a 10 um diameter
      cellsink = 1.40 *dtb! [m/d]
!ik      cellsink = 0.911 *dtb! [m/d]  for a 10 um diameter
!      shear = 86400. !shear in the mixed layer,   1.0  d-1
!      shear = 64800. !shear in the mixed layer,   0.75 d-1
      shear = 43200. !shear in the mixed layer,   0.5  d-1
      fsh = 0.163 * shear *dtb
      fse = 0.125 * 3.1415927 * cellsink * 100.
      alow1 = 0.002 !diameter of smallest particle [cm]
!ik      alow1 = 0.001 !diameter of smallest particle [cm]
      alow2 = alow1 * alow1
      alow3 = alow2 * alow1
!      alar1 = 1.0 !diameter of largest particle for size dependend aggregation and sinking [cm]
!      alar1 = 0.75 !diameter of largest particle for size dependend aggregation and sinking [cm]
      alar1 = 0.5 !diameter of largest particle for size dependend aggregation and sinking [cm]
      vsmall = 1.e-9 
      safe = 1.e-6     
      pupper = safe/((FractDim+safe)*cellmass)
      plower = 1./(1.1*cellmass)
      zdis = 0.01 / ((FractDim + 0.01)*cellmass)
      nmldmin = 0.1 ! minimum particle number in mixed layer

! Determine maximum sinking speed
      IF (mnproc.eq.1) THEN
      write(io_stdo_bgc,*)
      write(io_stdo_bgc,*)                                             &
     &'****************************************************************'
      write(io_stdo_bgc,*) 'HAMOCC aggregate sinking scheme:'
      write(io_stdo_bgc,*) ' Maximum sinking speed for aggregates of '
      write(io_stdo_bgc,*) ' maximum size ', alar1, ' cm is '
      write(io_stdo_bgc,*)   cellsink/dtb*(alar1/alow1)**SinkExp, ' m/day'
      endif

      alar2 = alar1 * alar1
      alar3 = alar2 * alar1
      TSFac = (alar1/alow1)**SinkExp
      TMFac = (alar1/alow1)**FractDim

! for shear aggregation of dust:
      dustd1 = 0.0001 !cm = 1 um, boundary between clay and silt
      dustd2=dustd1*dustd1
      dustd3=dustd2*dustd1
      dustsink = (9.81 * 86400. / 18.                & ! g * sec per day / 18.                 
     &         * (claydens - 1025.) / 1.567 * 1000.  & !excess density / dyn. visc.
     &         * dustd2 * 1.e-4)*dtb
      IF (mnproc.eq.1) THEN
      write(io_stdo_bgc,*) ' dust diameter (cm)', dustd1
      write(io_stdo_bgc,*) ' dust sinking speed (m/d)', dustsink / dtb
      if(dustsink.gt.cellsink) then 
        write(io_stdo_bgc,*) ' dust sinking speed greater than cellsink'
        dustsink=cellsink
        write(io_stdo_bgc,*) ' set dust sinking speed to cellsink'
      endif
      write(io_stdo_bgc,*)                                             &
     &'****************************************************************'
      endif

#endif /*AGG*/  



      RETURN
      END
