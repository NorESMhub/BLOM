! Copyright (C) 2001  Ernst Maier-Reimer, S. Legutke
! Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger, A. Moree, 
!                     C. Heinze
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


      SUBROUTINE CARCHM(kpie,kpje,kpke,kbnd,                                  &
                        pdlxp,pdlyp,pddpo,prho,pglat,omask,                   &
                        psicomo,ppao,pfu10,ptho,psao)
!******************************************************************************
!
!**** *CARCHM* - .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,             *MPI-MaD, HH*    10.04.01
!     - rename: ssso12(i,j,k)=sedlay(i,j,k,issso12 ) etc.; no equivalence statements
!     - rename: powasi(i,j,k )=powtra(i,j,1,ipowasi) etc.; no equivalence statements
!     - interfacing with ocean model
!
!     J.Tjiputra,            *BCCR*           09.18.08
!     - modified all carbon chemistry formulations following the OCMIP protocols
!
!     J.Schwinger,           *GFI, UiB*       2013-04-22    
!     - Use density prho consistent with MICOM for conversion to mol/kg
!     - Calculate solubility of O2 and N2 every timestep, consistent with
!       what is done for carbon chemistry. Array chemcm not used any more.
!     - Added J.Tjiputras code for cfc- and sf6-fluxes
!     - Cautious code clean-up
!
!     J.Schwinger,           *UNI-RESEARCH*   2017-08-30
!      - Moved the accumulation of global fields for output to routine
!        hamocc4bgc.
!
!     A.Moree,          *GFI, Bergen*   2018-04-12
!     - new version of carbon isotope code
!
!     J.Tjiputra,       *Uni Research, Bergen*   2018-04-12
!     - added preformed and saturated DIC tracers
!
!     J.Schwinger,      *Uni Research, Bergen*   2018-04-12
!     - moved accumulation of all output fields to seperate subroutine,
!       related code-restructuring
!     - dissolution of CaCO3 moved into main loop
!     - added sediment bypass preprocessor option
!
!     Purpose
!     -------
!     Inorganic carbon cycle.
!
!     Method
!     -------
!     Surface fluxes of CO2 / N2O / dms
!     Dissolution of calcium
!
!
!**** Parameter list:
!     ---------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *INTEGER* *kbnd*    - nb of halo grid points
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st dimension) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (2nd dimension) [m].
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!     *REAL*    *prho*    - density [g/cm^3].
!     *REAL*    *pglat*   - latitude of grid cells [deg north].
!     *REAL*    *omask*   - ocean mask.
!     *REAL*    *psicomo* - sea ice.
!     *REAL*    *ppao*    - sea level presure [Pascal].
!     *REAL*    *pfu10*   - forcing field wind speed.
!     *REAL*    *ptho*    - potential temperature.
!     *REAL*    *psao*    - salinity [psu].
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************
      USE mo_carbch
      USE mo_chemcon
      USE mo_biomod
      USE mo_sedmnt
      USE mo_control_bgc
      USE mo_param1_bgc 
      use mo_vgrid, only: dp_min,kbo,ptiestu

      implicit none

      INTEGER, intent(in) :: kpie,kpje,kpke,kbnd     
      REAL,    intent(in) :: pdlxp(kpie,kpje)
      REAL,    intent(in) :: pdlyp(kpie,kpje)
      REAL,    intent(in) :: pddpo(kpie,kpje,kpke)
      REAL,    intent(in) :: prho(kpie,kpje,kpke)
      REAL,    intent(in) :: pglat(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in) :: omask(kpie,kpje)
      REAL,    intent(in) :: psicomo(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in) :: ppao(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in) :: pfu10(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)
      REAL,    intent(in) :: psao(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)

      ! Local variables
      INTEGER :: i,j,k,l,js
      INTEGER, parameter :: niter=20

      REAL    :: supsat, undsa, dissol
      REAL    :: rpp0,fluxd,fluxu
      REAL    :: kwco2,kwo2,kwn2,kwdms,kwn2o
      REAL    :: scco2,sco2,scn2,scdms,scn2o
      REAL    :: Xconvxa
      REAL    :: oxflux,niflux,dmsflux,n2oflux
      REAL    :: ato2,atn2,atco2,pco2
      REAL    :: oxy,ani,anisa 
      REAL    :: rrho,t,t2,t3,t4,tk,tk100,prb,s,rs
      REAL    :: Kh,Khd,K1,K2,Kb,K1p,K2p,K3p,Ksi,Kw,Ks1,Kf,Kspc,Kspa
      REAL    :: tc,ta,sit,pt,ah1,ac,cu,cb,cc,tc_sat
      REAL    :: omega
#ifdef CFC
      REAL    :: atm_cfc11,atm_cfc12,atm_sf6,fact
      REAL    :: sch_11,sch_12,sch_sf,kw_11,kw_12,kw_sf
      REAL    :: flx11,flx12,flxsf,a_11,a_12,a_sf
#endif
#ifdef natDIC
      REAL    :: natcu,natcb,natcc
      REAL    :: natpco2,natfluxd,natfluxu,natomega
      REAL    :: natsupsat,natundsa,natdissol
#endif
#ifdef cisonew
      REAL    :: rco213,rco214
      REAL    :: dissol13,dissol14
      REAL    :: flux14d,flux14u,flux13d,flux13u
      REAL    :: atco213,atco214,pco213,pco214      
      REAL    :: frac_k,frac_aqg,frac_dicg
#endif


! set variables for diagnostic output to zero
       atmflx (:,:,:)=0.
       co2fxd   (:,:)=0.
       co2fxu   (:,:)=0.
#ifdef cisonew
       co213fxd (:,:)=0.
       co213fxu (:,:)=0.
       co214fxd (:,:)=0.
       co214fxu (:,:)=0.
#endif
       pco2d    (:,:)=0. 
       kwco2sol (:,:)=0.
       co2star(:,:,:)=0.
       co3    (:,:,:)=0.
       satoxy (:,:,:)=0.
       omegaA (:,:,:)=0.
       omegaC (:,:,:)=0.
#ifdef natDIC
       natpco2d   (:,:)=0. 
       natco3   (:,:,:)=0.
       natomegaA(:,:,:)=0.
       natomegaC(:,:,:)=0.
#endif

!$OMP PARALLEL DO PRIVATE(t,tk,tk100,s,rs,prb,Kh,Khd,K1,K2,Kb,K1p,K2p &
!$OMP  ,K3p,Ksi,Kw,Ks1,Kf,Kspc,Kspa,tc,ta,sit,pt,ah1,ac,cu,cb,cc,pco2 &
!$OMP  ,rpp0,scco2,scdms,sco2,oxy,ani,anisa,Xconvxa,kwco2,kwdms,kwo2  &
!$OMP  ,atco2,ato2,atn2,fluxd,fluxu,oxflux,tc_sat,niflux,n2oflux      &
!$OMP  ,dmsflux,omega,supsat,undsa,dissol                             &
#ifdef CFC
!$OMP  ,sch_11,sch_12,sch_sf,kw_11,kw_12,kw_sf,a_11,a_12,a_sf,flx11   &
!$OMP  ,flx12,flxsf,atm_cfc11,atm_cfc12,atm_sf6                       &
#endif
#ifdef natDIC
!$OMP  ,natcu,natcb,natcc,natpco2,natfluxd,natfluxu,natomega          &
!$OMP  ,natsupsat,natundsa,natdissol                                  &
#endif
#ifdef cisonew
!$OMP  ,atco213,atco214,rco213,rco214,pco213,pco214,frac_aqg          &
!$OMP  ,frac_dicg,flux13d,flux13u,flux14d,flux14u,dissol13,dissol14   &
#endif
!$OMP  )
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie

      IF(omask(i,j).gt.0.5.and.pddpo(i,j,k).GT.dp_min) THEN

! Carbon chemistry: Caculate equilibrium constants and solve for [H+] and
! carbonate alkalinity (ac)
      t    = min(40.,max(-3.,ptho(i,j,k)))
      t2   = t**2
      t3   = t**3
      t4   = t**4
      tk   = t + tzero
      tk100= tk/100.0
      s    = min(40.,max( 25.,psao(i,j,k)))
      rrho = prho(i,j,k)                   ! seawater density [g/cm3]
      prb  = ptiestu(i,j,k)*98060*1.027e-6 ! pressure in unit bars, 98060 = onem

      tc   = ocetra(i,j,k,isco212) / rrho  ! convert to mol/kg
      ta   = ocetra(i,j,k,ialkali) / rrho
      sit  = ocetra(i,j,k,isilica) / rrho
      pt   = ocetra(i,j,k,iphosph) / rrho
      ah1  = hi(i,j,k)

      CALL CARCHM_KEQUI(t,s,prb,Kh,Khd,K1,K2,Kb,Kw,Ks1,Kf,Ksi,             &
                        K1p,K2p,K3p,Kspc,Kspa)

      CALL CARCHM_SOLVE(s,tc,ta,sit,pt,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p, &
                        ah1,ac,niter)

      if(ah1.gt.0.) then 
        hi(i,j,k)=max(1.e-20,ah1)
      endif

! Determine CO2*, HCO3- and CO3-- concentrations (in mol/kg soln)
      cu = ( 2. * tc - ac ) / ( 2. + K1 / ah1 )
      cb = K1 * cu / ah1
      cc = K2 * cb / ah1
      co2star(i,j,k)=cu

! Carbonate ion concentration, convert from mol/kg to kmol/m^3 
      co3(i,j,k)  = cc * rrho 

#ifdef natDIC
      tc   = ocetra(i,j,k,inatsco212) / rrho  ! convert to mol/kg
      ta   = ocetra(i,j,k,inatalkali) / rrho
      ah1  = nathi(i,j,k)

      CALL CARCHM_SOLVE(s,tc,ta,sit,pt,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p, &
                        ah1,ac,niter)

      if(ah1.gt.0.) then 
        nathi(i,j,k)=max(1.e-20,ah1)
      endif

! Determine natural CO2*, HCO3- and CO3-- concentrations (in mol/kg soln)
      natcu = ( 2. * tc - ac ) / ( 2. + K1 / ah1 )
      natcb = K1 * natcu / ah1
      natcc = K2 * natcb / ah1
! Natural carbonate ion concentration, convert from mol/kg to kmol/m^3 
      natco3(i,j,k) = natcc * rrho
#endif

! solubility of O2 (Weiss, R.F. 1970, Deep-Sea Res., 17, 721-735) for moist air
! at 1 atm; multiplication with oxyco converts to kmol/m^3/atm
      oxy=ox0+ox1/tk100+ox2*alog(tk100)+ox3*tk100+s*(ox4+ox5*tk100+ox6*tk100**2)
      satoxy(i,j,k)=exp(oxy)*oxyco

      if (k.eq.1) then
! Determine CO2 pressure and fugacity (in micoatm)
! NOTE: equation below for pCO2 needs requires CO2 in mol/kg
      pco2 = cu * 1.e6 / Kh
#ifdef natDIC
      natpco2 = natcu * 1.e6 / Kh
#endif


! Schmidt numbers according to Wanninkhof (2014), Table 1
      scco2 = 2116.8 - 136.25*t + 4.7353*t2 - 0.092307*t3 + 0.0007555 *t4
      sco2  = 1920.4 - 135.6 *t + 5.2122*t2 - 0.10939 *t3 + 0.00093777*t4
      scn2  = 2304.8 - 162.75*t + 6.2557*t2 - 0.13129 *t3 + 0.0011255 *t4
      scdms = 2855.7 - 177.63*t + 6.0438*t2 - 0.11645 *t3 + 0.00094743*t4 
      scn2o = 2356.2 - 166.38*t + 6.3952*t2 - 0.13422 *t3 + 0.0011506 *t4
#ifdef CFC
      sch_11= 3579.2 - 222.63*t + 7.5749*t2 - 0.14595 *t3 + 0.0011874 *t4
      sch_12= 3828.1 - 249.86*t + 8.7603*t2 - 0.1716  *t3 + 0.001408  *t4
      sch_sf= 3177.5 - 200.57*t + 6.8865*t2 - 0.13335 *t3 + 0.0010877 *t4
#endif

! solubility of N2 (Weiss, R.F. 1970, Deep-Sea Res., 17, 721-735) for moist air
! at 1 atm; multiplication with oxyco converts to kmol/m^3/atm
       ani=an0+an1/tk100+an2*alog(tk100)+an3*tk100+s*(an4+an5*tk100+an6*tk100**2)
       anisa=exp(ani)*oxyco

! solubility of laughing gas  (Weiss and Price 1980, Marine Chemistry, 8, 347-359) 
! for moist air at 1 atm in kmol/m^3/atm
       rs=al1+al2/tk100+al3*log(tk100)+al4*tk100**2+s*(bl1+bl2*tk100+bl3*tk100**2)
       satn2o(i,j)=exp(rs)

#ifdef CFC
! solubility of cfc11,12 (mol/(l*atm)) (Warner and Weiss 1985) and
! sf6 from eq. 6 of Bullister et al. (2002)
! These are the alpha in (1b) of the ocmpic2 howto
      a_11 = exp(-229.9261 + 319.6552*(100/tk) + 119.4471*log(tk100)  &
     &         -1.39165*(tk100)**2 + s*(-0.142382 + 0.091459*(tk100)  &
     &         -0.0157274*(tk100)**2)) 
      a_12 = exp(-218.0971 + 298.9702*(100/tk) + 113.8049*log(tk100)  &
     &         -1.39165*(tk100)**2 + s*(-0.143566 + 0.091015*(tk100)  &
     &         -0.0153924*(tk100)**2)) 
      a_sf = exp(-80.0343  + 117.232 *(100/tk) +  29.5817*log(tk100)  &
     &         +s*(0.033518-0.0373942*(tk100)+0.00774862*(tk100)**2)) 
! conversion from mol/(l * atm) to kmol/(m3 * pptv) 
      a_11 = 1e-12 * a_11 
      a_12 = 1e-12 * a_12
      a_sf = 1e-12 * a_sf
#endif

! Transfer (piston) velocity kw according to Wanninkhof (2014), in units of ms-1 
       Xconvxa = 6.97e-07   ! Wanninkhof's a=0.251 converted to ms-1/(ms-1)^2 
       kwco2 = (1.-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660./scco2)**0.5
       kwo2  = (1.-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660./sco2)**0.5 
       kwn2  = (1.-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660./scn2)**0.5 
       kwdms = (1.-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660./scdms)**0.5 
       kwn2o = (1.-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660./scn2o)**0.5 
#ifdef CFC
       kw_11 = (1.-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660./sch_11)**0.5
       kw_12 = (1.-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660./sch_12)**0.5
       kw_sf = (1.-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660./sch_sf)**0.5
#endif


       atco2 = atm(i,j,iatmco2)
       ato2  = atm(i,j,iatmo2)
       atn2  = atm(i,j,iatmn2)
#ifdef cisonew
       atco213 = atm(i,j,iatmc13)
       atco214 = atm(i,j,iatmc14)
#endif

! Ratio P/P_0, where P is the local SLP and P_0 is standard pressure (1 atm). This is
! used in all surface flux calculations where atmospheric concentration is given as a
! mixing ratio (i.e. partial presure = mixing ratio*SLP/P_0 [atm])
       rpp0 = ppao(i,j)/101325.0

       fluxd=atco2*rpp0*kwco2*dtbgc*Kh*1e-6*rrho ! Kh is in mol/kg/atm. Multiply by rrho (g/cm^3) 
       fluxu=pco2      *kwco2*dtbgc*Kh*1e-6*rrho ! to get fluxes in kmol/m^2      
!JT set limit for CO2 outgassing to avoid negative DIC concentration, set minimum DIC concentration to 1e-5 kmol/m3 
       fluxu=min(fluxu,fluxd-(1e-5 - ocetra(i,j,k,isco212))*pddpo(i,j,1))
#ifdef natDIC
       natfluxd=atm_co2_nat*rpp0*kwco2*dtbgc*Kh*1e-6*rrho
       natfluxu=natpco2         *kwco2*dtbgc*Kh*1e-6*rrho 
       natfluxu=min(natfluxu,natfluxd-(1e-5 - ocetra(i,j,k,inatsco212))*pddpo(i,j,1))
#endif


! Calculate saturation DIC concentration in mixed layer
       ta = ocetra(i,j,k,ialkali) / rrho
       CALL carchm_solve_DICsat(s,atco2*rpp0,ta,sit,pt,Kh,K1,K2,Kb,Kw,Ks1,Kf, &
                               Ksi,K1p,K2p,K3p,tc_sat,niter)
       ocetra(i,j,k,  idicsat)=tc_sat * rrho ! convert mol/kg to kmol/m^3 
       ocetra(i,j,k+1,idicsat)=tc_sat * rrho ! k+1 = the rest of the mixed layer

#ifdef cisonew 
! Ocean-Atmosphere fluxes for carbon isotopes
       rco213=ocetra(i,j,1,isco213)/(ocetra(i,j,1,isco212)+safediv) ! Fraction DIC13 over total DIC
       rco214=ocetra(i,j,1,isco214)/(ocetra(i,j,1,isco212)+safediv) ! Fraction DIC14 over total DIC

       pco213 = pco2 * rco213 ! Determine water CO213 pressure and fugacity (microatm)
       pco214 = pco2 * rco214 ! Determine water CO214 pressure and fugacity (microatm)

! fractionation factors for 13C during air-sea gas exchange (Zhang et al. 1995, Orr et al. 2017)
       frac_k    = 0.99912                                 !Constant kinetic fractionation
       frac_aqg  = (0.0049*t - 1.31)/1000. + 1.  !Gas dissolution fractionation
       frac_dicg = (0.0144*t*(cc/(cc+cu+cb)) - 0.107*t + 10.53)/1000. + 1. !DIC to CO2 frac
       flux13d=atco213*rpp0*kwco2*dtbgc*Kh*1.e-6*rrho*frac_aqg*frac_k         
       flux13u=pco213      *kwco2*dtbgc*Kh*1.e-6*rrho*frac_aqg*frac_k/frac_dicg   
       flux14d=atco214*rpp0*kwco2*dtbgc*Kh*1.e-6*rrho*(frac_aqg**2)*(frac_k**2)           
       flux14u=pco214      *kwco2*dtbgc*Kh*1.e-6*rrho*(frac_aqg**2)*(frac_k**2)/(frac_dicg**2)  
#endif

! Update DIC
       ocetra(i,j,1,isco212)=ocetra(i,j,1,isco212)+(fluxd-fluxu)/pddpo(i,j,1)
#ifdef natDIC
       ocetra(i,j,1,inatsco212)=ocetra(i,j,1,inatsco212)+(natfluxd-natfluxu)/pddpo(i,j,1)
#endif
#ifdef cisonew
       ocetra(i,j,1,isco213)=ocetra(i,j,1,isco213)+(flux13d-flux13u)/pddpo(i,j,1)
       ocetra(i,j,1,isco214)=ocetra(i,j,1,isco214)+(flux14d-flux14u)/pddpo(i,j,1)
#endif

! Surface flux of oxygen
       oxflux=kwo2*dtbgc*(ocetra(i,j,1,ioxygen)-satoxy(i,j,1)*(ato2/196800)*rpp0)
       ocetra(i,j,1,ioxygen)=ocetra(i,j,1,ioxygen)-oxflux/pddpo(i,j,1)
! Surface flux of gaseous nitrogen (same piston velocity as for O2)
       niflux=kwn2*dtbgc*(ocetra(i,j,1,igasnit)-anisa*(atn2/802000)*rpp0) 
       ocetra(i,j,1,igasnit)=ocetra(i,j,1,igasnit)-niflux/pddpo(i,j,1)
! Surface flux of laughing gas (same piston velocity as for O2 and N2)
       n2oflux=kwn2o*dtbgc*(ocetra(i,j,1,ian2o)-satn2o(i,j)*atn2o*rpp0) 
       ocetra(i,j,1,ian2o)=ocetra(i,j,1,ian2o)-n2oflux/pddpo(i,j,1)
#ifdef CFC
! Surface fluxes for CFC: eqn. (1a) in ocmip2 howto doc(hyc)
!     flux of CFC: downward direction (mol/m**2/s)
!      flx11=kw_11*(a_11*cfc11_atm(i,j)*ppair/p0-trc(i,j,1,1))
!      flx12=kw_12*(a_12*cfc12_atm(i,j)*ppair/p0-trc(i,j,1,2))
!      unit should be in [kmol cfc m-2]
!      unit of [cfc11_atm(i,j)*ppair/p0] should be in [pptv]
!      unit of [flx11-12] is in [kmol / m2]

      IF (pglat(i,j).GE.10) THEN
       atm_cfc11=atm_cfc11_nh
       atm_cfc12=atm_cfc12_nh
       atm_sf6=atm_sf6_nh
      ELSE IF (pglat(i,j).LE.-10) THEN
       atm_cfc11=atm_cfc11_sh
       atm_cfc12=atm_cfc12_sh
       atm_sf6=atm_sf6_sh
      ELSE
       fact=(pglat(i,j)-(-10))/20.
       atm_cfc11=fact*atm_cfc11_nh+(1-fact)*atm_cfc11_sh
       atm_cfc12=fact*atm_cfc12_nh+(1-fact)*atm_cfc12_sh
       atm_sf6=fact*atm_sf6_nh+(1-fact)*atm_sf6_sh
      ENDIF

! Surface flux of cfc11
      flx11=kw_11*dtbgc*                                               &
     & (a_11*atm_cfc11*ppao(i,j)*9.86923*1e-6-ocetra(i,j,1,icfc11))
      ocetra(i,j,1,icfc11)=ocetra(i,j,1,icfc11)+flx11/pddpo(i,j,1)
! Surface flux of cfc12
      flx12=kw_12*dtbgc*                                               &
     & (a_12*atm_cfc12*ppao(i,j)*9.86923*1e-6-ocetra(i,j,1,icfc12))
      ocetra(i,j,1,icfc12)=ocetra(i,j,1,icfc12)+flx12/pddpo(i,j,1)
! Surface flux of sf6
      flxsf=kw_sf*dtbgc*                                               &
     & (a_sf*atm_sf6*ppao(i,j)*9.86923*1e-6-ocetra(i,j,1,isf6))
      ocetra(i,j,1,isf6)=ocetra(i,j,1,isf6)+flxsf/pddpo(i,j,1)
#endif

! Surface flux of dms
       dmsflux = kwdms*dtbgc*ocetra(i,j,1,idms)  
       ocetra(i,j,1,idms)=ocetra(i,j,1,idms)-dmsflux/pddpo(i,j,1)


! Save surface fluxes 
       atmflx(i,j,iatmco2)=fluxu-fluxd
       atmflx(i,j,iatmo2)=oxflux
       atmflx(i,j,iatmn2)=niflux
       atmflx(i,j,iatmn2o)=n2oflux
       atmflx(i,j,iatmdms)=dmsflux ! [kmol dms m-2 timestep-1]
#ifdef cisonew
       atmflx(i,j,iatmc13)=flux13u-flux13d
       atmflx(i,j,iatmc14)=flux14u-flux14d
#endif
#ifdef CFC
       atmflx(i,j,iatmf11)=flx11
       atmflx(i,j,iatmf12)=flx12
       atmflx(i,j,iatmsf6)=flxsf
#endif
#ifdef natDIC
       atmflx(i,j,iatmnco2)=natfluxu-natfluxd
#endif

! Save up- and downward components of carbon fluxes for output
       co2fxd(i,j)  = fluxd 
       co2fxu(i,j)  = fluxu 
#ifdef cisonew
       co213fxd(i,j)= flux13d 
       co213fxu(i,j)= flux13u 
       co214fxd(i,j)= flux14d 
       co214fxu(i,j)= flux14u 
#endif

! Save pco2 w.r.t. dry air for output
       pco2d(i,j) = cu * 1.e6 / Khd
#ifdef natDIC
       natpco2d(i,j) = natcu * 1.e6 / Khd
#endif

! Save product of piston velocity and solubility for output
       kwco2sol(i,j) = kwco2*Kh*1e-6


      endif ! k==1

! -----------------------------------------------------------------
! Deep ocean processes

! Determine Omega Calcite/Aragonite and dissolution of caco3 based on OmegaC:
!   omegaC=([CO3]*[Ca])/([CO3]sat*[Ca]sat)
!   Following Sarmiento and Gruber book, assumed that [Ca]=[Ca]sat
!   Thus, [CO3]sat=[CO3]/OmegaC. 
      omega = ( calcon * s / 35. ) * cc
      OmegaA(i,j,k) = omega / Kspa
      OmegaC(i,j,k) = omega / Kspc
      supsat=co3(i,j,k)-co3(i,j,k)/OmegaC(i,j,k)
      undsa=MAX(0.,-supsat)
      dissol=MIN(undsa,0.05*ocetra(i,j,k,icalc))
#ifdef natDIC
      natomega = ( calcon * s / 35. ) * natcc
      natOmegaA(i,j,k) = natomega / Kspa
      natOmegaC(i,j,k) = natomega / Kspc
      natsupsat=natco3(i,j,k)-natco3(i,j,k)/natOmegaC(i,j,k)
      natundsa=MAX(0.,-natsupsat)
      natdissol=MIN(natundsa,0.05*ocetra(i,j,k,inatcalc))
#endif
#ifdef cisonew
      dissol13=dissol*ocetra(i,j,k,icalc13)/(ocetra(i,j,k,icalc)+safediv)
      dissol14=dissol*ocetra(i,j,k,icalc14)/(ocetra(i,j,k,icalc)+safediv)
#endif
      ocetra(i,j,k,icalc)=ocetra(i,j,k,icalc)-dissol
      ocetra(i,j,k,ialkali)=ocetra(i,j,k,ialkali)+2.*dissol
      ocetra(i,j,k,isco212)=ocetra(i,j,k,isco212)+dissol
#ifdef natDIC
      ocetra(i,j,k,inatcalc)=ocetra(i,j,k,inatcalc)-natdissol
      ocetra(i,j,k,inatalkali)=ocetra(i,j,k,inatalkali)+2.*natdissol
      ocetra(i,j,k,inatsco212)=ocetra(i,j,k,inatsco212)+natdissol
#endif
#ifdef cisonew
      ocetra(i,j,k,icalc13)=ocetra(i,j,k,icalc13)-dissol13
      ocetra(i,j,k,isco213)=ocetra(i,j,k,isco213)+dissol13
      ocetra(i,j,k,icalc14)=ocetra(i,j,k,icalc14)-dissol14
      ocetra(i,j,k,isco214)=ocetra(i,j,k,isco214)+dissol14
#endif


#ifdef cisonew
! Decay of the ocean tracers that contain radioactive carbon 14C
      ocetra(i,j,k,isco214) = ocetra(i,j,k,isco214)*c14dec
      ocetra(i,j,k,idet14)  = ocetra(i,j,k,idet14) *c14dec
      ocetra(i,j,k,icalc14) = ocetra(i,j,k,icalc14)*c14dec
      ocetra(i,j,k,idoc14)  = ocetra(i,j,k,idoc14)*c14dec
      ocetra(i,j,k,iphy14)  = ocetra(i,j,k,iphy14)*c14dec
      ocetra(i,j,k,izoo14)  = ocetra(i,j,k,izoo14)*c14dec
#endif

      ! Save bottom level dissociation konstants for use in sediment module
      if( k==kbo(i,j) ) then

        keqb( 1,i,j)  = K1
        keqb( 2,i,j)  = K2
        keqb( 3,i,j)  = Kb
        keqb( 4,i,j)  = Kw
        keqb( 5,i,j)  = Ks1
        keqb( 6,i,j)  = Kf
        keqb( 7,i,j)  = Ksi
        keqb( 8,i,j)  = K1p
        keqb( 9,i,j)  = K2p
        keqb(10,i,j)  = K3p
        keqb(11,i,j)  = Kspc

      end if

      ENDIF ! omask>0.5
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO



! C14 decay in the sediment (could be moved to sediment part)
#ifdef cisonew
#ifndef sedbypass
        do k=1,ks
!$OMP PARALLEL DO  
        do j=1,kpje
        do i=1,kpie
        if(omask(i,j).gt.0.5) then
        sedlay(i,j,k,issso14)=sedlay(i,j,k,issso14)*c14dec
        sedlay(i,j,k,isssc14)=sedlay(i,j,k,isssc14)*c14dec
        powtra(i,j,k,ipowc14)=powtra(i,j,k,ipowc14)*c14dec
        endif
        enddo
        enddo
!$OMP END PARALLEL DO
        enddo
#endif
#endif


      RETURN
      END

