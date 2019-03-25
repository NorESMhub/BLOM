      SUBROUTINE BELEG_BGC(kpaufr,kpie,kpje,kpke,pddpo,ptiestw,prho,  &
     &                     omask,pglon,pglat,path)
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/beleg_bgc.f90,v $\\
!$Revision: 1.2 $\\
!$Date: 2004/11/12 15:37:21 $\\

!****************************************************************
!
!**** *BELEG_BGC* - initialize bgc variables.
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*     10.04.01
!
!     J.Schwinger,      *GFI, Bergen*     2013-10-21
!     - corrected units of tracer fields at initialisation to
!       mol/kg as tracers are passed to ocean model first,
!       unit conversion is done using in situ density
!     - code clean up
!
!     I. Kriest, GEOMAR, 11.08.2016
!     - included T-dependence of cyanobacteria growth
!     - modified stoichiometry for denitrification
! 
!     A.Moree,          *GFI, Bergen*   2018-04-12
!     - new version of carbon isotope code
!
!     J.Tjiputra,       *Uni Research, Bergen*   2018-04-12
!     - added preformed and saturated DIC tracers
!
!     J.Schwinger,      *Uni Research, Bergen*   2018-04-12
!     - moved reading of namelist and initialisation of dust to 
!       ini_hamocc.F90
!     - initialisation of tracer fields is skipped if the model is
!       restarted
!     - added sediment bypass preprocessor option
!
!
!     Purpose
!     -------
!     - set start values for bgc variables.
!
!     Method
!     -------
!     - bgc variables are initialized unless it is a restart run
!     - physical constants are defined
!     
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!     *INTEGER*   *kpaufr*  - 1/0 flag, 1 indicating a restart run
!     *INTEGER*   *kpie*    - 1st dimension of model grid.
!     *INTEGER*   *kpje*    - 2nd dimension of model grid.
!     *INTEGER*   *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*      *pddpo*   - size of grid cell (3rd dimension) [m].
!     *REAL*      *ptiestw* - depth of layer interfaces [m].
!     *REAL*      *prho*    - density [g/cm^3].
!     *REAL*      *omask*   - ocean mask.
!     *REAL*      *pglon*   - longitude of grid cell [deg].
!     *REAL*      *pglat*   - latitude  of grid cell [deg].
!     *CHARACTER* *path*    - path to input data files.
!
!
!**********************************************************************

      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      USE mo_control_bgc
      use mo_param1_bgc 
      USE mod_xc, only: mnproc

      implicit none      

      INTEGER :: kpaufr,kpie,kpje,kpke
      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: ptiestw(kpie,kpje,kpke+1)
      REAL :: prho (kpie,kpje,kpke)
      REAL :: omask(kpie,kpje)
      REAL :: pglon(kpie,kpje)
      REAL :: pglat(kpie,kpje)
      character(len=*) :: path

      ! local variables
      INTEGER :: i,j,k,l,ii,jj
      integer :: p_joff,p_ioff
      REAL :: north, south
#ifdef cisonew
      REAL :: rco213,rco214,alpha14,beta13,beta14,d13C_atm,d14cat
#endif
#ifdef AGG
      REAL :: shear,snow
#else
      REAL :: dustd1, dustd2, dustsink
#endif


!
! Initialize overall time step counter.
!
      ldtbgc = 0

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
! restart file (if BOXATM or DIFFAT is activated), this will be overwritten
! later.
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
      phytomi=1.e-11		!i.e. 1e-5 mmol P/m3 minimum concentration of phyto plankton (?js)
      grami=1.e-10              !i.e. 1e-5 mmol P/m3 minimum concentration of zooplankton

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
      epsher = 0.9          !dimensionless fraction -fraction of grazing egested
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
      rcalc = 35.  ! calcium carbonate to organic phosphorous production ratio
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
       dmspar(5)=1.25*0.109784522E-01  !2*0.02   production with delsil
       dmspar(4)=1.25*0.107638502E+00  !2*1.3e-5 production with delcar
       dmspar(3)=0.0864 ! following Kloster et al., 06 Table 1 with 50% reduction to reduce bacterial removal and increase dms emissions
       dmspar(2)=0.0011 ! following Kloster et al., 06 Table 1
       dmspar(1)=10.              !2*5.     production with temp


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


#ifdef ANTC14
      DO  i=1,kpie
        ii=1+(i+p_ioff-1)   ! global i-index for giph_g
        north=1.
        DO  j=1,kpje
          jj=1+(j+p_joff-1) ! global j-index for giph_g
          north=pglat(ii,jj)
          if (north .gt. 20.) then
            Rbomb(i,j) = D14C_north
          endif
          if (north .lt. -20.) then
            Rbomb(i,j) = D14C_south
          endif
          if ((north .le. 20.).and.(north .ge. -20.)) then
            Rbomb(i,j) = D14C_equator
          endif
!         WRITE(io_stdo_bgc,*)'Rbomb: ',i,j,north,Rbomb(i,j)
        ENDDO
      ENDDO
#endif


#ifdef FB_BGC_OCE
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
        abs_oce(i,j,k)=1.
      ENDDO
      ENDDO
      ENDDO
#endif


!
! Initialisation of ocean tracers and sediment
!

! Initialise ocean tracers with WOA and GLODAP data. This is done even in case
! of a restart since some tracers (e.g. C-isotopes) might not be in the restart 
! file and aufr.f90 instead expects an initialised field.
      call profile_gd(kpie,kpje,kpke,pglon,pglat,ptiestw,omask,TRIM(path))

! If this is a restart run initialisation is done in aufr.F90 
      IF(kpaufr.EQ.1) RETURN

      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
        IF (omask(i,j) .GT. 0.5 ) THEN
          ! convert WOA tracers kmol/m^3 -> mol/kg; GLODAP dic and alk
          ! are already in mol/kg. We need these units here, since after 
          ! initialisation the tracer field is passed to the ocean model
          ! first where units are mol/kg.
          ocetra(i,j,k,iphosph) = ocetra(i,j,k,iphosph)/prho(i,j,k)
          ocetra(i,j,k,ioxygen) = ocetra(i,j,k,ioxygen)/prho(i,j,k)
          ocetra(i,j,k,iano3)   = ocetra(i,j,k,iano3)  /prho(i,j,k)
          ocetra(i,j,k,isilica) = ocetra(i,j,k,isilica)/prho(i,j,k)
#ifdef cisonew
          ! d13C based on Eide data is read in above (profile_gd)                        
          ! Convert to 13C using model initial (ie GLODAP) total C
          ! If restarting, this is redone with model total C from restart in aufr_bgc.F90 
          beta13=ocetra(i,j,k,isco213)/1000.+1.
          ocetra(i,j,k,isco213) = ocetra(i,j,k,isco212)*beta13*re1312/(1.+beta13*re1312)

          ! 14C is read in as small delta14C (calculated from R. Key, 2003 and Eide et al. 2017)
          ! Convert to 14C using model total C, and normalize by c14fac to prevent numerical errors
          beta14=ocetra(i,j,k,isco214)/1000.+1.
          ocetra(i,j,k,isco214) = ocetra(i,j,k,isco212)*beta14*re14to/c14fac
#endif
        ENDIF
      ENDDO
      ENDDO
      ENDDO

! Initialise remaining ocean tracers 
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie

      IF(omask(i,j) .GT. 0.5) THEN
         ocetra(i,j,k,igasnit)=1.e-10
         ocetra(i,j,k,idoc)   =1.e-8
         ocetra(i,j,k,iphy)   =1.e-8 
         ocetra(i,j,k,izoo)   =1.e-8 
         ocetra(i,j,k,idet)   =1.e-8 
         ocetra(i,j,k,icalc)  =0. 
         ocetra(i,j,k,iopal)  =1.e-8 
         ocetra(i,j,k,ian2o)  =0. 
         ocetra(i,j,k,idms)   =0. 
         ocetra(i,j,k,ifdust) =0. 
         ocetra(i,j,k,iiron)  =fesoly
         ocetra(i,j,k,iprefo2)=0.
         ocetra(i,j,k,iprefpo4)=0.
         ocetra(i,j,k,iprefalk)=0.
         ocetra(i,j,k,iprefdic)=0.
         ocetra(i,j,k,idicsat)=1.e-8
         hi(i,j,k)            =1.e-8
         co3(i,j,k)           =0.
         co2star(i,j,k)       =20.e-6	   
#ifdef AGG
! calculate initial numbers from mass, to start with appropriate size distribution
         snow = (ocetra(i,j,k,iphy)+ocetra(i,j,k,idet))*1.e+6
         ocetra(i,j,k,inos)   = snow / cellmass / (FractDim+1.)
         ocetra(i,j,k,iadust) =0. 
#endif /*AGG*/
#ifdef ANTC14
         ocetra(i,j,k,iantc14)=ocetra(i,j,k,isco214)
#endif
#ifdef CFC
         ocetra(i,j,k,icfc11)   =0.
         ocetra(i,j,k,icfc12)   =0.
         ocetra(i,j,k,isf6)     =0.
#endif
#ifdef natDIC
         nathi(i,j,k)           =1.e-8
         natco3(i,j,k)          =0.
         ocetra(i,j,k,inatcalc) =0. 
#endif
#ifdef cisonew
         rco213=ocetra(i,j,k,isco213)/(ocetra(i,j,k,isco212)+safediv)
         rco214=ocetra(i,j,k,isco214)/(ocetra(i,j,k,isco212)+safediv)
         ocetra(i,j,k,iphy13) =ocetra(i,j,k,iphy)*rco213*bifr13
         ocetra(i,j,k,iphy14) =ocetra(i,j,k,iphy)*rco214*bifr14
         ocetra(i,j,k,izoo13) =ocetra(i,j,k,izoo)*rco213*bifr13
         ocetra(i,j,k,izoo14) =ocetra(i,j,k,izoo)*rco214*bifr14
         ocetra(i,j,k,idoc13) =ocetra(i,j,k,idoc)*rco213*bifr13
         ocetra(i,j,k,idoc14) =ocetra(i,j,k,idoc)*rco214*bifr14
         ocetra(i,j,k,idet13) =ocetra(i,j,k,idet)*rco213*bifr13
         ocetra(i,j,k,idet14) =ocetra(i,j,k,idet)*rco214*bifr14
         ocetra(i,j,k,icalc13)=ocetra(i,j,k,icalc)*rco213
         ocetra(i,j,k,icalc14)=ocetra(i,j,k,icalc)*rco214
#endif
      ENDIF ! omask > 0.5

      ENDDO
      ENDDO
      ENDDO

! Initialise preformed tracers in the mixed layer; note that the 
! whole field has been initialised to zero above
      DO k=1,kmle
      DO j=1,kpje
      DO i=1,kpie

      IF(omask(i,j) .GT. 0.5) THEN
         ocetra(i,j,k,iprefo2) =ocetra(i,j,k,ioxygen)
         ocetra(i,j,k,iprefpo4)=ocetra(i,j,k,iphosph)
         ocetra(i,j,k,iprefalk)=ocetra(i,j,k,ialkali)
         ocetra(i,j,k,iprefdic)=ocetra(i,j,k,isco212)
      ENDIF

      ENDDO
      ENDDO
      ENDDO


! Initial values for sediment
#ifndef sedbypass
      DO  k=1,ks
      DO  j=1,kpje
      DO  i=1,kpie 
      IF(omask(i,j) .GT. 0.5) THEN
         powtra(i,j,k,ipowaic)=ocetra(i,j,kbo(i,j),isco212)
         powtra(i,j,k,ipowaal)=ocetra(i,j,kbo(i,j),ialkali)
         powtra(i,j,k,ipowaph)=ocetra(i,j,kbo(i,j),iphosph)
         powtra(i,j,k,ipowaox)=ocetra(i,j,kbo(i,j),ioxygen)
         powtra(i,j,k,ipown2) =0.
         powtra(i,j,k,ipowno3)=ocetra(i,j,kbo(i,j),iano3)
         powtra(i,j,k,ipowasi)=ocetra(i,j,kbo(i,j),isilica)      
         sedlay(i,j,k,issso12)=1.e-8
         sedlay(i,j,k,isssc12)=1.e-8
         sedlay(i,j,k,issster)=30.
         sedlay(i,j,k,issssil)=1.e-8
         sedhpl(i,j,k)        =hi(i,j,kbo(i,j))
#ifdef cisonew
         rco213=ocetra(i,j,kbo(i,j),isco213)/(ocetra(i,j,kbo(i,j),isco212)+safediv)
         rco214=ocetra(i,j,kbo(i,j),isco214)/(ocetra(i,j,kbo(i,j),isco212)+safediv)
         powtra(i,j,k,ipowc13)=powtra(i,j,k,ipowaic)*rco213*bifr13
	 powtra(i,j,k,ipowc14)=powtra(i,j,k,ipowaic)*rco214*bifr14
         sedlay(i,j,k,issso13)=sedlay(i,j,k,issso12)*rco213*bifr13
         sedlay(i,j,k,issso14)=sedlay(i,j,k,issso12)*rco214*bifr14
         sedlay(i,j,k,isssc13)=sedlay(i,j,k,isssc12)*rco213
         sedlay(i,j,k,isssc14)=sedlay(i,j,k,isssc12)*rco214
#endif
      ELSE
         powtra(i,j,k,ipowno3)=rmasks
         powtra(i,j,k,ipown2) =rmasks
         powtra(i,j,k,ipowaic)=rmasks
         powtra(i,j,k,ipowaal)=rmasks
         powtra(i,j,k,ipowaph)=rmasks
         powtra(i,j,k,ipowaox)=rmasks
         powtra(i,j,k,ipowasi)=rmasks
         sedlay(i,j,k,issso12)=rmasks
         sedlay(i,j,k,isssc12)=rmasks
         sedlay(i,j,k,issssil)=rmasks
         sedlay(i,j,k,issster)=rmasks
         sedlay(i,j,k,issssil)=rmasks
         sedhpl(i,j,k)        =rmasks
#ifdef cisonew
         powtra(i,j,k,ipowc13)=rmasks
	 powtra(i,j,k,ipowc14)=rmasks
         sedlay(i,j,k,issso13)=rmasks
         sedlay(i,j,k,issso14)=rmasks
         sedlay(i,j,k,isssc13)=rmasks
         sedlay(i,j,k,isssc14)=rmasks
#endif
      ENDIF
      ENDDO
      ENDDO
      ENDDO

      ! last and final sediment layer
      DO  l=1,nsedtra
      DO  j=1,kpje
      DO  i=1,kpie
 	 burial(i,j,l)=0. 
      ENDDO
      ENDDO
      ENDDO
#endif


      RETURN
      END
