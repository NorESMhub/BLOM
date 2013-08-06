      SUBROUTINE BELEG_BGC(kpie,kpje,kpke,psao,ptho,pddpo,ptiestu,      &
     &                     omask,gila_g,giph_g,path,path_len)
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
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!     
!     Purpose
!     -------
!     - set start values for  bgc variables.
!
!     Method
!     -------
!     - bgc variables are initialized. They might be overwritten
!       if read from restart by call of AUFR_BGC.
!     - physical constants are defined
!     - fields used for budget calculations (should be in extra SBR!)
!       and as boundary conditions are set to zero.
!     
!
!**   Interface.
!     ----------
!
!     *CALL*       *BELEG_BGC*
!
!     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!     *COMMON*     *PARAM1_BGC.h* - declaration of ocean/sediment tracer.
!     *COMMON*     *COMMO1_BGC.h* - ocean/sediment tracer arrays.
!     *COMMON*     *UNITS_BGC.h*  - std I/O logical units.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of grid cell (3rd dimension) [m].
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************

      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      USE mo_control_bgc
      use mo_param1_bgc 
      
!      USE MO_COMMO1

#ifdef PDYNAMIC_BGC
      USE mo_dynamic
#endif /* PDYNAMIC_BGC */

      USE mod_xc, only: mnproc,xcminr


      implicit none      

      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: ptho (kpie,kpje,kpke)
      REAL :: psao (kpie,kpje,kpke)
      REAL :: ptiestu(kpie,kpje,kpke+1)
      REAL :: omask(kpie,kpje)
!      REAL :: gila_g(kpie,kpje)
!      REAL :: giph_g(kpie,kpje)
      REAL :: gila_g(kpie*2,kpje*2)
      REAL :: giph_g(kpie*2,kpje*2)
      REAL :: north, south, oarea
      REAL :: rho2,s,rrr
      INTEGER :: i,j,k,l,ii,jj,m,kpie,kpje,kpke,kmon
      character*(*) path
      integer path_len

#ifdef AGG
      REAL :: shear,zmini,talar1,snow, checksink
#endif 

#ifndef AGG
      REAL :: dustd1, dustd2, dustsink
#endif

      REAL :: xpi,rad,radi,rmissing
      integer p_joff,p_ioff

#ifdef CCSMCOUPLED
      namelist /bgcnml/ atm_co2
#endif

      xpi       = 4.*ATAN(1.)
      rad       = xpi/180.
      radi      = 1./rad

!
! Initialize overall time step counter.
!
      ldtbgc = 0
!
!
#ifndef DIFFAT            
#ifdef CCSMCOUPLED
!
! Obtain the CCSM value of atmospheric co2 concentration.
!
      open (unit=io_nml,file='ocn_in',status='old',action='read',      &
     &      recl=80)
      read (unit=io_nml,nml=BGCNML)
      close (unit=io_nml)
#else
         atm_co2 = 278.
#endif
      IF (mnproc.eq.1) THEN
        write(io_stdo_bgc,*) 'HAMOCC: atmospheric co2:',atm_co2
      ENDIF
#ifdef __c_isotopes
         atm_c13 = 276.2 ! preindustrial atmospheric d13c=-6.5 permil --> 276.2ppm? test js 15082006 ok.
                        ! piston velocity ~8m/yr -> equilibration time atmosphere/ocean ~500 years
         atm_c14 = 274.
#endif
         atm_o2  = 196800.
         atm_n2  = 802000.
#endif   
!
! Biology
!
!ik note that the unit is kmol/m3!
      phytomi=1.e-11		!i.e. 1e-5 mmol P/m3 minimum concentration of phyto plankton (?js)
      grami=1.e-10              !i.e. 1e-5 mmol P/m3 minimum concentration of zooplankton

!ik addded parameter definition; taken from OCPROD.F
      remido=0.025*dtb      !1/d -remineralization rate (of DOM)
      dyphy=0.008*dtb       !1/d -mortality rate of phytoplankton 
      zinges=0.5            !dimensionless fraction -assimilation efficiency
      epsher=0.8            !dimensionless fraction -fraction of grazing egested
      grazra=1.0*dtb        !1/d -grazing rate
      spemor=5.*1.e6*dtb      !1/d -mortality rate
      gammap=0.03*dtb       !1/d -exudation rate
      gammaz=0.03*dtb       !1/d -excretion rate
      ecan=0.95             ! fraction of mortality as PO_4
      pi_alpha=0.02*0.4     ! initial slope of production vs irradiance curve (alpha) (0.002 for 10 steps per day)

#ifdef __c_isotopes
! fractionation for photosynthesis plafr13 and plafr14 valid for particles
! for avoiding too many tracers, surface gross rates work with reduced
! values bifr13 and bifr14
       plafr13=1.                 ! planktonic fractionatio 13C   (never used) (HAMOCC3.1: 0.98) 
       plafr14=1.
       bifr13=0.98                ! biogenic fractionation ?
       bifr14=0.98
#endif

! half sat. constants, note that the units are kmol/m3 !
!JT      bkphy  = 1.e-7    !i.e. 0.04 mmol P/m3
      bkphy  = 2.e-7    !i.e. 0.04 mmol P/m3
      bkzoo  = 4.e-8    !i.e. 0.04 mmol P/m3
      bkopal = 1.5e-6    !i.e. 1.0 mmol Si/m3

!sinking speed
!JT      wpoc  = 10.*dtb       !m/d  iris : 5.
      wpoc  =  5.*dtb       !m/d  iris : 5.
      wcal  = 30.*dtb       !m/d 
!JT      wopal = 50.*dtb       !m/d  iris : 60
      wopal = 60.*dtb       !m/d  iris : 60

      
! deep see remineralisation constants

      drempoc  = 0.03*dtb     !1/d
      dremdoc  = 0.004*dtb      !1/d
      dphymor  = 0.1*dtb      !1/d
      dzoomor  = 0.02*dtb      !1/d
      dremopal = 0.01*dtb      !1/d
      dremn2o  = 0.01*dtb       !1/d
      dremsul  = 0.005*dtb            ! remineralization rate for sulphate reduction 
      dremcalc = 0.075 *dtb              ! 0.2 -> 0.02 js10072006 : slightly overdone --> 0.075
      
#ifdef AGG
      drempoc  = 0.05 *dtb      !1/d       re-introduced 09062006 -- too strong?? 0.1 -> 0.05 02072007
      dremopal = 0.03 *dtb   ! js 4.7.2006 0.0033 from .01/3. (60/20 m/d)
      dremcalc = 0.2*dtb
      dphymor  = 0.2 *dtb      !1/d
#endif

! nirogen fixation by blue green algae
   
      bluefix=0.005*dtb     !1/d

! extended redfield ratio declaration
! Note: stoichiometric ratios are based on Takahashi etal. (1985)
! P:N:C:-O2 + 1:16:122:172

      ro2ut=172. 
      rcar=122.
      rnit=16.
      rnoi=1./rnit
      rnit23 = ro2ut*2./3. !114 rno3 * 2 / 3 for nitrate reduction during denitrification
      rnit13 = ro2ut*1./3. !57   rno3 * 1 / 3 for n2 production during denitrification

      rcalc = 35. ! iris 40 !calcium carbonate to organic phosphorous production ratio
      ropal = 25. ! iris 25 !opal to organic phosphorous production ratio      
      calmax=0.20
      gutc = 0. 

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
      perc_diron = 0.035 * 0.01 / 55.85
!      riron= 5.*rcar*1.e-6  
      riron= 3.*rcar*1.e-6       ! 15/06/06 changed ka 
      fesoly=0.6*1.e-9 !max. diss. iron concentration in deep water 
!      relaxfe = 0.05/365.*dtb    ! set to this to get rid of sfc P in
      relaxfe = 0.005/365.*dtb    ! default

!ka surface input of nutrients, C, alk to balance sediment
       oarea=3.5969e14
       sco212_sfc=19.     *1.e9/oarea/365.
       alkali_sfc=29.47541*1.e9/oarea/365.
       phosph_sfc=4.0     /rcar*1.e9/oarea/365.
       ano3_sfc  =4.0     *rnit/rcar*1.e9/oarea/365.
       silica_sfc=5.5     *1.e9/oarea/365.

! decay parameter for sco214, HalfLive = 5730years
      c14dec=(alog(2.)/(5730*365))*dtb
      eins=1.
      c14ret=eins-c14dec
! Ratm: normalized atmospheric ratio of C-14/C-12, D14Catm: atmospheric Delta C-14
      D14Catm=0.             ! D14Catm=0. only for equilibrium runs
      Ratm=1+D14Catm/1000.      

!                        
! Set constants for calculation of dms ( mo_carbch )
! Parameters are a result from kettle optimisation 02.03.04

       dmspar(6)=0.100000000E-07  !0 half saturation microbial
!jt       dmspar(5)=1.25*0.107638502E+00  !2*1.3e-5 production with delsil
!jt       dmspar(4)=1.25*0.109784522E-01  !2*0.02   production with delcar
       dmspar(5)=1.25*0.109784522E-01  !2*0.02   production with delsil
       dmspar(4)=1.25*0.107638502E+00  !2*1.3e-5 production with delcar
       dmspar(3)=0.0096  !4.8e-5  !2*1.6e-3 microbial consumption
       dmspar(2)=0.0075  !0.0003  !2*0.005  photolysis
       dmspar(1)=10.              !2*5.     production with temp


      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)                                             &
     &'****************************************************************'
      WRITE(io_stdo_bgc,*)                                             &
     &'* '
      WRITE(io_stdo_bgc,*)                                             &
     &'* Values of BELEG_BGC variables : '
#ifndef DIFFAT      
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atm_co2      = ',atm_co2      
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atm_o2       = ',atm_o2           
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atm_n2       = ',atm_n2 
#endif  
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
     &'*                              dremdoc      = ',dremdoc   
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dremopal     = ',dremopal   
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dphymor      = ',dphymor   
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dzoomor      = ',dzoomor   
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
     &'*                              rnit23       = ',rnit23
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rnit13       = ',rnit13
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rcalc        = ',rcalc
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              ropal        = ',ropal
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              gutc         = ',gutc
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
     &'*                              c14dec       = ',c14dec
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              D14Catm      = ',D14Catm
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              Ratm         = ',Ratm
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
      Stick = 0.40
      cellmass = 0.012 / rnit ![nmol P]
!ik      cellmass = 0.0039/ rnit ![nmol P] for a 10 um diameter
      cellsink = 1.40 *dtb! [m/d]
!ik      cellsink = 0.911 *dtb! [m/d]  for a 10 um diameter
      shear = 86400. !shear in upper 5 layers, 1 d-1
      fsh = 0.163 * shear *dtb
      fse = 0.125 * 3.1415927 * cellsink * 100.
      alow1 = 0.002 !diameter of smallest particle [cm]
!ik      alow1 = 0.001 !diameter of smallest particle [cm]
      alow2 = alow1 * alow1
      alow3 = alow2 * alow1
      alar1 = 1.0 !diameter of largest particle for size dependend aggregation and sinking [cm]
      vsmall = 1.e-9 
      safe = 1.e-6     
      pupper = safe/((FractDim+safe)*cellmass)
      plower = 1./(1.1*cellmass)
      zdis = 0.01 / ((FractDim + 0.01)*cellmass)

!ik check max possible sinking speed in relation to min.
!ik layer thinkness and time step for all standard layers, except
!ik the bottom layer.
!ik if max possible sinking speed (per time step) is greater
!ik than min layer thickness, decrease max. length for sinking and
!ik aggregation

      zmini = 8000.
      DO  j=1,kpje
      DO  i=1,kpie
      DO  k=1,kbo(i,j)-1
        if(pddpo(i,j,k).gt.0.5) then
         zmini=min(pddpo(i,j,k),zmini)
        endif 
      ENDDO
      ENDDO
      ENDDO

      CALL xcminr(zmini)
            
      checksink =(zmini/cellsink)**(1./SinkExp)*alow1 
      if(checksink.lt.alar1) then 
      write(io_stdo_bgc,*) 'Allowed max. length for sinking               &
     &   with min. depth of '                                             &
     & , zmini, ' m for layers 1-(kbo-1) and time step of ',dtb           &
     & ,' days is' , checksink                                            &
     & ,'cm, which is smaller than prescribed value of', alar1, ' cm'
        talar1 = alar1
        alar1 = checksink
      write(io_stdo_bgc,*) 'Set max. length for sinking and aggregation   &
     &   from ',talar1,' to ', alar1
      endif

      alar2 = alar1 * alar1
      alar3 = alar2 * alar1
      TSFac = (alar1/alow1)**SinkExp
      TMFac = (alar1/alow1)**FractDim

!ik check the maximum possible sinking speed for the last layer (which
!ik may be smaller than zmini, and write to array alar1max, tsfmax, tmfmax

      DO j=1,kpje
      DO i=1,kpie
         alar1max(i,j) = alar1
         TSFmax(i,j) = TSFac
         TMFmax(i,j) = TMFac
         if(omask(i,j).gt.0.5) then

!ik evaluate safe length scale for size dependent sinking and
!ik aggregation, and the resulting sinking rate and aggregation rate.

          checksink = (pddpo(i,j,kbo(i,j))/cellsink)**(1./SinkExp)        &
     &                    *alow1
          if(checksink.lt.alar1) then
           alar1max(i,j) = checksink
           TSFmax(i,j) = (checksink/alow1)**SinkExp
           TMFmax(i,j) = (checksink/alow1)**FractDim
           write(io_stdo_bgc,*) 'resetting alar1 to',checksink,'at i =', &
     &     i,' j = ',j,' k = ', kbo(i,j), ' with dz = ',                 &
     &     pddpo(i,j,kbo(i,j))  
          endif
        ENDIF
      ENDDO
      ENDDO

! for shear aggregation of dust:
      dustd1 = 0.0001 !cm = 1 um, boundary between clay and silt
      dustd2=dustd1*dustd1
      dustd3=dustd2*dustd1
      dustsink = (9.81 * 86400. / 18.                & ! g * sec per day / 18.                 
     &         * (claydens - 1025.) / 1.567 * 1000.  & !excess density / dyn. visc.
     &         * dustd2 * 1.e-4)*dtb
      write(io_stdo_bgc,*) 'dust diameter (cm)', dustd1
      write(io_stdo_bgc,*) 'dust sinking speed (m/d)', dustsink / dtb
      if(dustsink.gt.cellsink) then 
        write(io_stdo_bgc,*) 'dust sinking speed greater than cellsink'
        dustsink=cellsink
        write(io_stdo_bgc,*) 'set dust sinking speed to cellsink'
      endif

#endif /*AGG*/  

      DO  j=1,kpje
      DO  i=1,kpie
      DO  k=1,8
      DO  kmon=1,12
         chemcm(i,j,k,kmon)=rmasko
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO  i=1,kpie
      DO  j=1,kpje
      DO  k=1,kpke
         aksp(i,j,k)=rmasko
      ENDDO
      ENDDO
      ENDDO

      DO  j=1,kpje
      DO  i=1,kpie
      DO  k=1,ks
      DO  l=1,nsedtra
         sedlay(i,j,k,l)=0.
	 burial(i,j,l)=0.   ! last and final sediment layer
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO  j=1,kpje
      DO  i=1,kpie
         expoor(i,j)=0.
         expoca(i,j)=0.
         exposi(i,j)=0.
      ENDDO
      ENDDO

!
!  Initial values for aquatic (advected) ocean tracers
! 
!ka calculate average profile for DIC,Alk,P,N,Si,O  
!      call profile(kpie,kpje,kpke,ptiestu,omask)

! Initialise ocean tracers with WOA and GLODAP data
      call profile_gd(kpie,kpje,kpke,gila_g,giph_g,ptiestu,omask,TRIM(path))


      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie

      IF(omask(i,j) .GT. 0.5) THEN
!         ocetra(i,j,k,isco212)=2.27e-3
!         ocetra(i,j,k,ialkali)=2.37e-3
!         ocetra(i,j,k,iphosph)=2.17e-6
!         ocetra(i,j,k,ioxygen)=3.e-4
         ocetra(i,j,k,igasnit)=0.0
!         ocetra(i,j,k,iano3)  =ocetra(i,j,k,iphosph)*rnit 
!         ocetra(i,j,k,iano3)  =2.17e-6*rnit ! old 32.e-6
!         ocetra(i,j,k,isilica)=1.2e-4
!JT         ocetra(i,j,k,idoc)   =1.e-10
         ocetra(i,j,k,idoc)   =100*1.e-10
         ocetra(i,j,k,iphy)   =1.e-8 
         ocetra(i,j,k,izoo)   =1.e-8 
         ocetra(i,j,k,idet)   =1.e-8 
         ocetra(i,j,k,icalc)  =0. 
         ocetra(i,j,k,iopal)  =1.e-8 
#ifdef __c_isotopes
	 ocetra(i,j,k,isco214)=0.75*2.27e-3 !Paddy: oldest water: 1600y --> X0.83 
#endif
         ocetra(i,j,k,ian2o)  =0. 
         ocetra(i,j,k,idms)   =0. 
         ocetra(i,j,k,ifdust) =0. 
         ocetra(i,j,k,iiron)  =0.6*1.e-9
         hi(i,j,k)            =3.e-9
         co3(i,j,k)           =0.
#ifdef __c_isotopes
         ocetra(i,j,k,isco213)=ocetra(i,j,k,isco212)     ! adjusted to reference ratio 13C/12C=1 (*100)!
         ocetra(i,j,k,isco214)=ocetra(i,j,k,isco212)
!         ocetra(i,j,k,isco213)=2.27e-3     ! adjusted to reference ratio 13C/12C=1 (*100)!
!         ocetra(i,j,k,isco214)=2.27e-3
         ocetra(i,j,k,idet13) =1.e-8
         ocetra(i,j,k,idet14) =1.e-8
         ocetra(i,j,k,icalc13)=0.
         ocetra(i,j,k,icalc14)=0.
#endif

#ifdef AGG
! calculate initial numbers from mass, to start with appropriate size distribution
         snow = (ocetra(i,j,k,iphy)+ocetra(i,j,k,idet))*1.e+6
         ocetra(i,j,k,inos)   = snow / cellmass / (FractDim+1.)
         ocetra(i,j,k,iadust) =0. 
#endif /*AGG*/

#ifdef ANTC14
         ocetra(i,j,k,iantc14)=ocetra(i,j,k,isco214)
#endif
#ifdef PCFC
         ocetra(i,j,k,icfc11)=0.
         ocetra(i,j,k,icfc12)=0.
#endif
!      ELSE
!         ocetra(i,j,k,iphosph)=0.
!         ocetra(i,j,k,isilica)=0.
!         ocetra(i,j,k,ioxygen)=0.
!         ocetra(i,j,k,ialkali)=0.
!         ocetra(i,j,k,isco212)=0.
!         ocetra(i,j,k,iano3)  =0.
!         ocetra(i,j,k,igasnit)=0.
!         ocetra(i,j,k,idoc)   =0.
!         ocetra(i,j,k,iphy)   =0.
!         ocetra(i,j,k,izoo)   =0.
!         ocetra(i,j,k,idet)   =0.
!         ocetra(i,j,k,icalc)  =0.
!         ocetra(i,j,k,iopal)  =0.
!	 ocetra(i,j,k,isco214)=0.
!         ocetra(i,j,k,ian2o)  =0.
!         ocetra(i,j,k,idms)   =0.
!         ocetra(i,j,k,ifdust) =0.
!         ocetra(i,j,k,iiron)  =0.
!         hi(i,j,k)            =0.
!         co3(i,j,k)           =0.
!#ifdef AGG
!         ocetra(i,j,k,inos)   =0.
!         ocetra(i,j,k,iadust) =0. 
!#endif /*AGG*/
!#ifdef ANTC14
!         ocetra(i,j,k,iantc14)  =0.
!#endif
!#ifdef PCFC
!         ocetra(i,j,k,icfc11)   =0.
!         ocetra(i,j,k,icfc12)   =0.
!#endif
      ENDIF

      ENDDO
      ENDDO
      ENDDO

!convert ocean tracer concentrations from kmol/kg to k

      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
      IF(omask(i,j) .GT. 0.5) THEN
         s=MAX(25.,psao(i,j,k))
         rrr=rho2(s,ptho(i,j,k),0.)*1.e-3
      do l=1,nocetra
         ocetra(i,j,k,l)=ocetra(i,j,k,l)*rrr
      enddo
      ENDIF
      ENDDO
      ENDDO
      ENDDO

!read in dust fields
     CALL GET_DUST(kpie,kpje,kpke,omask,path,path_len)
!      dusty(:,:,:)=1.e-12
!
!  Initial values for sediment pore water tracer.
      DO  k=1,ks
      DO  j=1,kpje
      DO  i=1,kpie 
!      IF(bolay(i,j) .GT. 0.) THEN
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
#ifdef __c_isotopes
         sedlay(i,j,k,issso13)=1.e-8
         sedlay(i,j,k,isssc13)=1.e-8
         sedlay(i,j,k,issso14)=1.e-8
         sedlay(i,j,k,isssc14)=1.e-8
#endif
         sedlay(i,j,k,issster)=30.
!JT         sedlay(i,j,k,issssil)=3.
         sedlay(i,j,k,issssil)=1.e-8
         sedhpl(i,j,k)        =hi(i,j,kbo(i,j))
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
#ifdef __c_isotopes
         sedlay(i,j,k,issso13)=rmasks
         sedlay(i,j,k,isssc13)=rmasks
         sedlay(i,j,k,issso14)=rmasks
         sedlay(i,j,k,isssc14)=rmasks
#endif
         sedlay(i,j,k,issssil)=rmasks
         sedlay(i,j,k,issster)=rmasks
         sedlay(i,j,k,issssil)=rmasks
         sedhpl(i,j,k)        =rmasks
      ENDIF
      ENDDO
      ENDDO
      ENDDO

!
! values for sediment fluxes
!
      
      DO j=1,kpje
      DO i=1,kpie
      DO k=1,npowtra
        sedfluxo(i,j,k)=0.     
      ENDDO
      ENDDO
      ENDDO
!
! values for sedimentation
!
      
      DO j=1,kpje
      DO i=1,kpie
        prorca(i,j)=0.
        prcaca(i,j)=0.
        silpro(i,j)=0.
        produs(i,j)=0.
      ENDDO
      ENDDO
      
!
!  atmospheric diffusion parameters.
!    
#ifdef DIFFAT
      DO  j=1,kpje
      DO  i=1,kpie 
         atm(i,j,iatmco2) = 278.
         atm(i,j,iatmo2)  = 196800.
         atm(i,j,iatmn2)  = 802000.
#ifdef __c_isotopes
         atm(i,j,iatmc13) = 278.-(278.*0.0065)
         atm(i,j,iatmc14) = 278.-(278.*0.0065)**2
#endif
         atdifv(i,j)=1.
      ENDDO
      ENDDO
#endif
   
#ifdef PCFC
      DO  i=1,kpie
        ii=1+(i+p_ioff-1)*2 ! global i-index for giph_g
        north=1.
        DO  j=1,kpje
          jj=1+(j+p_joff-1)*2 ! global j-index for giph_g
          IF(jj<=2) CYCLE
          north=giph_g(ii,jj)
          if (north .gt. 10.) then
            cfc_int(i,j) = 1.
          endif
          if (north .lt. -10.) then
            cfc_int(i,j) = 0.
          endif
          if ((north .le. 10.).and.(north .ge. -10.)) then
            cfc_int(i,j) = (north +10.)/20.
          endif
!         WRITE(io_stdo_bgc,*)'cfc_int: ',i,j,north,cfc_int(i,j)
        ENDDO
      ENDDO
#endif
#ifdef ANTC14
      DO  i=1,kpie
        ii=1+(i+p_ioff-1)*2 ! global i-index for giph_g
        north=1.
        DO  j=1,kpje
          jj=1+(j+p_joff-1)*2 ! global j-index for giph_g
          north=giph_g(ii,jj)
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
#ifdef PDYNAMIC_BGC
        DO k=1,nbgcdyn
        DO j=1,kpje
        DO i=1,kpie
          bgcdyntmp(i,j,k) = 0.
        ENDDO
        ENDDO
        ENDDO

        DO l=1,kdtot
        DO k=1,nbgcdyn
        DO j=1,kpje
        DO i=1,kpie
          bgcdyn(i,j,k,l) = 0.
        ENDDO
        ENDDO
        ENDDO
        ENDDO

        DO j=1,kpje
        DO i=1,kpie
          bgc_zmld(i,j) = 0.
          bgc_nmld(i,j) = 0.
          nbgc_mld(i,j) = 0
        ENDDO
        ENDDO
#endif /* PDYNAMIC_BGC */


#ifdef FB_BGC_OCE
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
        abs_oce(i,j,k)=1.
      ENDDO
      ENDDO
      ENDDO
#endif

      RETURN
      END
