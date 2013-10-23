      SUBROUTINE CARCHM(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,psao,ppao,   &
     &                  ptho,prho,psicomo,pfu10,ptiestu,omask)

!**********************************************************************
!
!**** *CARCHM* - .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!     - rename: ssso12(i,j,k)=sedlay(i,j,k,issso12 ) etc.; no equivalence statements
!     - rename: powasi(i,j,k )=powtra(i,j,1,ipowasi) etc.; no equivalence statements
!     - interfacing with ocean model
!
!     J.Tjiputra, BCCR 09.18.08
!     - modified all carbon chemistry formulations following the OCMIP protocols
!
!     J.Schwinger       *GFI, UiB*       2013-04-22    
!     - Use density prho consistent with MICOM for conversion to mol/kg
!     - Calculate solubility of O2 and N2 every timestep, consistent with
!       what is done for carbon chemistry. Array chemcm not used any more.
!     - Added J.Tjiputras code for cfc- and sf6-fluxes
!     - Cautious code clean-up
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
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st dimension) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (2nd dimension) [m].
!     *REAL*    *psao*    - salinity [psu].
!     *REAL*    *ppao*    - sea level presure [Pascal].
!     *REAL*    *ptho*    - potential temperature.
!     *REAL*    *prho*    - density [g/cm^3].
!     *REAL*    *psicomo* - sea ice.
!     *REAL*    *pfu10*   - forcing field wind speed.
!     *REAL*    *ptiestu* - depth of layer centres [m].
!     *REAL*    *omask*   - ocean mask.
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
!      USE mo_timeser_bgc
      USE mo_control_bgc
      USE mo_bgcmean
      USE mo_param1_bgc 

      implicit none

      INTEGER :: kpie,kpje,kpke      
      REAL    :: pddpo(kpie,kpje,kpke)
      REAL    :: pdlxp(kpie,kpje)
      REAL    :: pdlyp(kpie,kpje)
      REAL    :: psao(kpie,kpje,kpke)
      REAL    :: ppao(kpie,kpje)
      REAL    :: ptho(kpie,kpje,kpke)
      REAL    :: prho(kpie,kpje,kpke)
      REAL    :: psicomo(kpie,kpje)
      REAL    :: pfu10(kpie,kpje)
      REAL    :: ptiestu(kpie,kpje,kpke)
      REAL    :: omask(kpie,kpje)

      REAL    :: supsat, undsa, dissol
      REAL    :: fluxd,fluxu
      REAL    :: kwco2,kwo2,kwdms
      REAL    :: scco2,sco2,scdms
      REAL    :: contppm, Xconvxa
      REAL    :: oxflux,niflux,dmsflux,n2oflux
      REAL    :: ato2,atn2,atco2,pco2
      REAL    :: oxy,ani,anisa
!Tjiputra update=for list of new variables==============================
      REAL    :: t,tk,tk100,rs,prb,s,rrho,invtk,dlogtk,is,is2,sqrtis,s2,sqrts,s15,scl
      REAL    :: tmp,nkhwe74,Kh,K1,K2,Kb,K1p,K2p,K3p,Ksi,Kw,Ks1,Kf,Kspc,Kspa
      REAL    :: zprb,zprb2,deltav,deltak,borat,sti,ft,tc,ta,sit,pt,ah1
      REAL    :: hso4,hsi,hf,hpo4,ab,aw,ac,ah2o,ah2,erel,cu,cb,cc
      REAL    :: omega


      REAL, DIMENSION(11) :: a0, a1, a2, b0, b1, b2, lnkpok0
      DATA a0 /-25.5, -15.82, -29.48, -25.60, -18.03, -9.78, -48.76, &
     &       -46., -14.51, -23.12, -26.57/
      DATA a1 /0.1271, -0.0219, 0.1622, 0.2324, 0.0466, -0.0090, &
     &       0.5304, 0.5304, 0.1211, 0.1758, 0.2020/
      DATA a2 /0.0, 0.0, 2.608e-3, -3.6246e-3, 0.316e-3, &
     &       -0.942e-3, 0.0, 0.0, -0.321e-3, -2.647e-3, -3.042e-3/
      DATA b0 /-3.08e-3, 1.13e-3, -2.84e-3, -5.13e-3, -4.53e-3,  &
     &       -3.91e-3, -11.76e-3, -11.76e-3, -2.67e-3, -5.15e-3,& 
     &       -4.08e-3/
      DATA b1 /0.0877e-3, -0.1475e-3, 0.0, 0.0794e-3, 0.09e-3, &
     &       0.054e-3, 0.3692e-3, 0.3692e-3, 0.0427e-3, &
     &       0.09e-3, 0.0714e-3/
      DATA b2 /0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      INTEGER            :: i,j,k,l,js,jit,iter
      INTEGER, parameter :: niter=20
      REAL,    parameter :: eps=5.e-5
      REAL,    parameter :: rgas = 83.131  ! Gas constant, value as 
                                           ! used by Millero 1995

!End of Tjiputra update=for list of variables===========================
#ifdef __c_isotopes
      REAL :: r13,r14,rat13,rat14
      REAL :: flux14d,flux14u,flux13d,flux13u
      REAL :: atc13,atc14
#endif

#ifdef ANTC14
      REAL :: fantc14d,fantc14u
#endif
 
      INTEGER, DIMENSION(kpie,kpje) :: ind1,ind2
      REAL, DIMENSION(kpie,kpje,ddm) :: wghts
      REAL, DIMENSION(kpie,kpje) :: aux2d_co2fxd
      REAL, DIMENSION(kpie,kpje) :: aux2d_co2fxu
      REAL, DIMENSION(kpie,kpje) :: aux2d_pco2
      REAL, DIMENSION(kpie,kpje) :: aux2d_kwco2
      REAL, DIMENSION(kpie,kpje) :: aux2d_oxflux
      REAL, DIMENSION(kpie,kpje) :: aux2d_niflux
      REAL, DIMENSION(kpie,kpje) :: aux2d_dmsflux
      REAL, DIMENSION(kpie,kpje) :: aux2d_dms
#ifdef CFC
      REAL, DIMENSION(kpie,kpje) :: aux2d_cfc11
      REAL, DIMENSION(kpie,kpje) :: aux2d_cfc12
      REAL, DIMENSION(kpie,kpje) :: aux2d_sf6
      REAL :: sch_11,sch_12,sch_sf,kw_11,kw_12,kw_sf
      REAL :: flx11,flx12,flxsf,a_11,a_12,a_sf
#endif
      aux2d_co2fxd=0
      aux2d_co2fxu=0
      aux2d_pco2=0
      aux2d_kwco2=0
      aux2d_oxflux=0
      aux2d_niflux=0
      aux2d_dmsflux=0
      aux2d_dms=0      
#ifdef CFC
      aux2d_cfc11=0.
      aux2d_cfc12=0.
      aux2d_sf6=0.
#endif


      contppm=1/0.35e-3

!Tjiputra update===========================================================================

!$OMP PARALLEL DO                                                               &
#ifdef CFC
!$OMP&PRIVATE(t,tk,tk100s,rs,prb,invtk,dlogtk,is,is2,sqrtis,s2,sqrts,s15,scl,   &
!$OMP&        tmp,nKhwe74,Kh,K1,K2,Kb,K1p,K2p,K3p,Ksi,Kw,Ks1,Kf,Kspc,           &
!$OMP&        Kspa,deltav,deltak,zprb,zprb2,lnkpok0,borat,sti,ft,tc,ta,         &
!$OMP&        sit,pt,ah1,iter,jit,hso4,hf,hsi,hpo4,ab,aw,ac,ah2o,ah2,           &
!$OMP&        erel,cu,cb,cc,pco2,scco2,scdms,sco2,oxy,ani,anisa,Xconvxa,        &
!$OMP&        kwco2,kwdms,kwo2,atco2,ato2,atn2,fluxd,fluxu,oxflux,              &
!$OMP&        niflux,n2oflux,dmsflux,omega,sch_11,sch_12,sch_sf,kw_11,          &
!$OMP&        kw_12,kw_sf,a_11,a_12,a_sf,flx11,flx12,flxsf)
#else
!$OMP&PRIVATE(t,tk,tk100s,rs,prb,invtk,dlogtk,is,is2,sqrtis,s2,sqrts,s15,scl,   &
!$OMP&        tmp,nKhwe74,Kh,K1,K2,Kb,K1p,K2p,K3p,Ksi,Kw,Ks1,Kf,Kspc,           &
!$OMP&        Kspa,deltav,deltak,zprb,zprb2,lnkpok0,borat,sti,ft,tc,ta,         &
!$OMP&        sit,pt,ah1,iter,jit,hso4,hf,hsi,hpo4,ab,aw,ac,ah2o,ah2,           &
!$OMP&        erel,cu,cb,cc,pco2,scco2,scdms,sco2,oxy,ani,anisa,Xconvxa,        &
!$OMP&        kwco2,kwdms,kwo2,atco2,ato2,atn2,fluxd,fluxu,oxflux,              &
!$OMP&        niflux,n2oflux,dmsflux,omega)
#endif
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
      IF(omask(i,j).gt.0.5.and.pddpo(i,j,k).GT.dp_min) THEN
       t = ptho(i,j,k)
       tk = t + tzero
       tk100 = tk/100.0
       s = MAX(25.,psao(i,j,k))
       rrho = prho(i,j,k)
       prb = ptiestu(i,j,k)*98060*1.027e-6 ! pressure in unit bars, 98060 = onem
       invtk = 1.0 / tk
       dlogtk = log(tk)
       is = 19.924 * s / ( 1000. - 1.005 * s )
       is2 = is * is
       sqrtis = SQRT(is)
       s2     = s * s
       sqrts  = SQRT(s)
       s15    = s**1.5
       scl    = s * salchl

! Kh = [CO2]/ p CO2
! Weiss (1974)   [mol/kg/atm]
       tmp = 9345.17 * invtk - 60.2409 + 23.3585 * log( tk/100. )
       nKhwe74 = tmp + s * ( 0.023517 - 0.00023656 * tk + 0.0047036e-4 * tk * tk )
       Kh      = exp( nKhwe74 )
! K1 = [H][HCO3]/[H2CO3]   ; K2 = [H][CO3]/[HCO3]
! Millero p.664 (1995) using Mehrbach et al. data on seawater scale
       K1 = 10**( -1.0 * ( 3670.7 * invtk - 62.008 + 9.7944 * dlogtk - 0.0118 * s + 0.000116 * s2 ) )
       K2 = 10**( -1.0 * ( 1394.7 * invtk + 4.777 - 0.0184 * s + 0.000118 * s2 ) )
! Kb = [H][BO2]/[HBO2] !
! Millero p.669 (1995) using DATA from Dickson (1990)
       Kb = exp( ( -8966.90 - 2890.53 * sqrts - 77.942 * s + 1.728 * s15 - 0.0996 * s2 ) * invtk +     &
     &    ( 148.0248 + 137.1942 * sqrts + 1.62142 * s ) + ( -24.4344 - 25.085 * sqrts - 0.2474 * s ) * &
     &    dlogtk + 0.053105 * sqrts * tk )
! K1p = [H][H2PO4]/[H3PO4] ; K2p = [H][HPO4]/[H2PO4] ; K3p = [H][PO4]/[HPO4]
! DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
       K1p = exp( -4576.752 * invtk + 115.525 - 18.453 * dlogtk + ( -106.736 * invtk + 0.69171 ) *     &
     &    sqrts + ( -0.65643 * invtk - 0.01844 ) * s )
       K2p = exp( -8814.715 * invtk + 172.0883 - 27.927 * dlogtk + ( -160.340 * invtk + 1.3566 ) *     &
     &    sqrts + ( 0.37335 * invtk - 0.05778 ) *s );
       K3p = exp( -3070.75 * invtk - 18.141 + ( 17.27039 * invtk + 2.81197 ) * sqrts + ( -44.99486 *   &
     &    invtk - 0.09984 ) * s );
! Ksi = [H][SiO(OH)3]/[Si(OH)4]
! Millero p.671 (1995) using data from Yao and Millero (1995)
       Ksi = exp( -8904.2 * invtk + 117.385 - 19.334 * dlogtk + ( -458.79 * invtk + 3.5913 ) * sqrtis  &
     &    + ( 188.74 * invtk - 1.5998) * is + ( -12.1652 * invtk + 0.07871) * is2 + log(1.0-0.001005*s))
! Kw = [H][OH] 
! Millero p.670 (1995) using composite data
       Kw = exp( -13847.26 * invtk + 148.9652 - 23.6521 * dlogtk + ( 118.67 * invtk - 5.977 + 1.0495 * &
     &    dlogtk ) * sqrts - 0.01615 * s)
! Ks = [H][SO4]/[HSO4]
! Dickson (1990, J. chem. Thermodynamics 22, 113)
       Ks1 = exp( -4276.1 * invtk + 141.328 - 23.093 * dlogtk + ( -13856. * invtk + 324.57 - 47.986 *   &
     &    dlogtk ) * sqrtis + ( 35474. * invtk - 771.54 + 114.723 * dlogtk ) * is - 2698. * invtk *    &
     &    is**1.5 + 1776. * invtk * is2 + log(1.0 - 0.001005 * s ) )
! Kf = [H][F]/[HF]
! Dickson and Riley (1979) -- change pH scale to total
       Kf = exp( 1590.2 * invtk - 12.641 + 1.525 * sqrtis + log( 1.0 - 0.001005 * s ) + log( 1.0 + (   &
     &    0.1400 / 96.062 ) * scl / Ks1 ) )
! Kspc (calcite)
! apparent solubility product of calcite : Kspc = [Ca2+]T [CO32-]T
! where $[]_T$ refers to the equilibrium total (free + complexed) ion concentration.
!          Mucci 1983 mol/kg-soln
       Kspc = 10**( -171.9065 - 0.077993 * tk + 2839.319 / tk + 71.595 * log10( tk ) + ( - 0.77712 +    &
     &    0.0028426 * tk + 178.34 / tk ) * sqrts - 0.07711 * s + 0.0041249 * s15 );
! Kspa (aragonite)
! apparent solubility product of aragonite : Kspa = [Ca2+]T [CO32-]T
! where $[]_T$ refers to the equilibrium total (free + complexed) ion concentration.
!          Mucci 1983 mol/kg-soln
!jt    Kspa = 10^( -171.945 - 0.077993 * tk + 2903.293 / tk  + 71.595 * log10( tk ) ...
!jt           +( -0.068393 + 0.0017276 * tk + 88.135 / tk ) * sqrts ...
!jt           - 0.10018 * s + 0.0059415 * s15 );
!---------------------- Pressure effect on Ks (Millero, 95) ----------
! index: K1 1, K2 2, Kb 3, Kw 4, Ks 5, Kf 6, 
!        Kspc 7, Kspa 8, K1p 9, K2p 10, K3p 11
       DO js = 1,11
        deltav      = a0(js) + a1(js) * t + a2(js) * t * t
        deltak      = b0(js) + b1(js) * t + b2(js) * t * t
        zprb        = prb / ( rgas * tk )
        zprb2       = prb * zprb
        lnkpok0(js) = - ( deltav * zprb + 0.5 * deltak * zprb2 )
       ENDDO
       K1   = K1   * exp( lnkpok0(1)  )
       K2   = K2   * exp( lnkpok0(2)  )
       Kb   = Kb   * exp( lnkpok0(3)  )
       Kw   = Kw   * exp( lnkpok0(4)  )
       Ks1  = Ks1  * exp( lnkpok0(5)  )
       Kf   = Kf   * exp( lnkpok0(6)  )
       Kspc = Kspc * exp( lnkpok0(7)  )
!JT    Kspa = Kspa * exp( lnkpok0(8)  )
       K1p  = K1p  * exp( lnkpok0(9)  )
       K2p  = K2p  * exp( lnkpok0(10) )
       K3p  = K3p  * exp( lnkpok0(11) )
!
! Calculate concentrations for borate, sulfate, and fluoride
! Uppstrom (1974)
      borat = bor1 * scl * bor2
! Morris & Riley (1966)
      sti = 0.14 * scl / 96.062
! Riley (1965)
      ft = 0.000067 * scl / 18.9984
!         ========================
!         ==                    ==
!         ==        MAIN        ==
!         ==        ====        ==
!         ==                    ==
!         ========================
      tc   = ocetra(i,j,k,isco212) / rrho  ! convert to mol/kg
      ta   = ocetra(i,j,k,ialkali) / rrho
      sit  = ocetra(i,j,k,isilica) / rrho
      pt   = ocetra(i,j,k,iphosph) / rrho
      ah1  = hi(i,j,k)
      iter = 0
      DO jit = 1,niter
       hso4 = sti / ( 1. + Ks1 / ( ah1 / ( 1. + sti / Ks1 ) ) )
       hf   = 1. / ( 1. + Kf / ah1 )
       hsi  = 1./ ( 1. + ah1 / Ksi )
       hpo4 = ( K1p * K2p * ( ah1 + 2. * K3p ) - ah1**3 ) /    & 
     &        ( ah1**3 + K1p * ah1**2 + K1p * K2p * ah1 + K1p * K2p * K3p )
       ab   = borat / ( 1. + ah1 / Kb )
       aw   = Kw / ah1 - ah1 / ( 1. + sti / Ks1 )
       ac   = ta + hso4 - sit * hsi - ab - aw + ft * hf - pt * hpo4
       ah2o = SQRT( ( tc - ac )**2 + 4. * ( ac * K2 / K1 ) * ( 2. * tc - ac ) )
       ah2  = 0.5 * K1 / ac *( ( tc - ac ) + ah2o )
       erel = ( ah2 - ah1 ) / ah2
       if (abs( erel ).ge.eps) then
        ah1 = ah2
        iter = iter + 1
       else
        ah1 = ah2
        exit
       endif
      ENDDO

      if(ah1.gt.0.) then 
        hi(i,j,k)=max(1.e-20,ah1)
      endif
! Determine CO2*, HCO3- and CO3-- concentrations (in mol/kg soln)
      cu = ( 2. * tc - ac ) / ( 2. + K1 / ah1 )
      cb = K1 * cu / ah1
      cc = K2 * cb / ah1
! Carbonate ion concentration, convert from mol/kg to kmol/m^3 
      co3(i,j,k)  = cc * rrho 

! Determine CO2 pressure and fugacity (in micoatm)
! NOTE: equation below for pCO2 needs requires CO2 in mol/kg
      pco2 = cu * 1.e6 / Kh

! solubility of O2 (Weiss, R. F. (1970), Deep-Sea Res., 17, 721-735)
      oxy=ox0+ox1/tk100+ox2*alog(tk100)+ox3*tk100+s*(ox4+ox5*tk100+ox6*tk100**2)
      satoxy(i,j,k)=exp(oxy)*oxyco

      if (k.eq.1) then
       scco2 = 2073.1 - 125.62*ptho(i,j,k) + 3.6276*ptho(i,j,k)**2  &
     &       - 0.043219*ptho(i,j,1)**3   ! Schmidt number
       scdms = 2674.0-147.12*ptho(i,j,1)+3.726*ptho(i,j,1)**2       &
     &       - 0.038*ptho(i,j,1)**3
       sco2  = 1638.0 - 81.83*ptho(i,j,1) + 1.483*ptho(i,j,1)**2    & 
     &       - 0.008004*ptho(i,j,1)**3  
#ifdef CFC
! Schmidt number for CFC11,12 based on Zheng et al. 98 (JGR), SF6 taken from Wanninkhof 1992
      sch_11 = 3501.8 - 210.31*ptho(i,j,1) + 6.1851*ptho(i,j,1)**2  &
     &         -0.07513*ptho(i,j,1)**3
      sch_12 = 3845.4 - 228.95*ptho(i,j,1) + 6.1908*ptho(i,j,1)**2  &
     &         -0.06743*ptho(i,j,1)**3
      sch_sf = 3531.6 - 231.40*ptho(i,j,1) + 7.2168*ptho(i,j,1)**2  &
     &         -0.090558*ptho(i,j,1)**3 
#endif

! solubility of N2 (Weiss, R. F. (1970), Deep-Sea Res., 17, 721-735)
       ani=an0+an1/tk100+an2*alog(tk100)+an3*tk100+s*(an4+an5*tk100+an6*tk100**2)
       anisa=exp(ani)*oxyco

! solubility of laughing gas  (Weiss, 1974)
       rs=al1+al2*(100./tk)+al3*log(tk100)+s*(bl1+bl2*(tk100)+bl3*(tk100)**2)
       satn2o(i,j)=atn2o*exp(rs)

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


! Compute the transfer (piston) velocity kw in m/s according to Wanninkhof(1992,eq.3) 
! see ocmip2 guide 'howto' chapter 2.2 eqn. (2)
       Xconvxa = 9.3611e-07
       kwco2 = (1.-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660./scco2)**0.5
       kwdms = (1.-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660./scdms)**0.5 
       kwo2  = (1.-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660./sco2)**0.5 
#ifdef CFC
       kw_11 = (1.-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660./sch_11)**0.5
       kw_12 = (1.-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660./sch_12)**0.5
       kw_sf = (1.-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660./sch_sf)**0.5
#endif


#ifdef DIFFAT      
         atco2 = atm(i,j,iatmco2)
         ato2  = atm(i,j,iatmo2)
         atn2  = atm(i,j,iatmn2)
#ifdef __c_isotopes
         atc13 = atm(i,j,iatmc13)
         atc14 = atm(i,j,iatmc14)
#endif
#elif CCSMCOUPLED
         atco2 = atm(i,j,iatmco2)
         ato2  = atm_o2
         atn2  = atm_n2
#ifdef __c_isotopes
         atc13 = atco2*atm_c13/atm_co2
         atc14 = atco2*atm_c14/atm_co2
#endif
#else
         atco2 = atm_co2
         ato2  = atm_o2
         atn2  = atm_n2
#ifdef __c_isotopes
         atc13=  atm_c13
         atc14=  atm_c14
#endif /*__c_isotopes*/
#endif


       fluxd=atco2*kwco2*dtbgc*Kh*1e-6 ! JT replaced ak0 with Kh*1e-6
       fluxu=pco2 *kwco2*dtbgc*Kh*1e-6 ! JT replaced ak0 with Kh*1e-6


!proxies d13C stuff
#ifdef __c_isotopes
	 Roc13=ocetra(i,j,1,isco213)/(ocetra(i,j,1,isco212)+1e-25)
	 Roc14=ocetra(i,j,1,isco214)/(ocetra(i,j,1,isco212)+1e-25)
#ifdef DIFFAT
         rat13=atm(i,j,iatmc13)/atm(i,j,iatmco2)
         rat14=atm(i,j,iatmc14)/atm(i,j,iatmco2)
#else
         rat13=atc13/atco2
         rat14=atc14/atco2
#endif

	 flux13d=atc13*kwco2*dtbgc*Kh*1e-6              ! *ppao(i,j)/101300. (atm to ocean)
         flux13u=pco2 *kwco2*dtbgc*Kh*1e-6*Roc13*0.9935 ! *ppao(i,j)/101300. (ocean to atm,
                                                                        !twofold fractionation through evaporation)

	 flux14d=atc14*kwco2*dtbgc*Kh*1e-6              ! *ppao(i,j)/101300.
         flux14u=pco2 *kwco2*dtbgc*Kh*1e-6*Roc14*0.987  ! *ppao(i,j)/101300. (twofold fractionation through evaporation)
#endif /*__c_isotopes*/
#if defined(DIFFAT) || defined(CCSMCOUPLED) 	   
!         atm(i,j,iatmco2)=atm(i,j,iatmco2)+(fluxu-fluxd)*contppm
!jt         atmflx(i,j,iatmco2)=(fluxu-fluxd)*contppm
         atmflx(i,j,iatmco2)=(fluxu-fluxd)
#ifdef __c_isotopes
!         atm(i,j,iatmc13)=atm(i,j,iatmc13)+(flux13u-flux13d)*contppm
!         atm(i,j,iatmc14)=atm(i,j,iatmc14)+(flux14u-flux14d)*contppm
         atmflx(i,j,iatmc13)=(flux13u-flux13d)*contppm
         atmflx(i,j,iatmc14)=(flux14u-flux14d)*contppm
#endif
#endif /*DIFFAT*/	 

       ocetra(i,j,1,isco212)=ocetra(i,j,1,isco212)+(fluxd-fluxu)/pddpo(i,j,1)
#ifdef __c_isotopes
         ocetra(i,j,1,isco213)=                                     &
     &   ocetra(i,j,1,isco213)+(flux13d-flux13u)/pddpo(i,j,1)
         ocetra(i,j,1,isco214)=                                     &
     &   ocetra(i,j,1,isco214)+(flux14d-flux14u)/pddpo(i,j,1)
#endif

! Surface flux of oxygen
       oxflux=kwo2*dtbgc*(ocetra(i,j,1,ioxygen)-satoxy(i,j,1)*(ato2/196800)) ! *ppao(i,j)/101300.
       ocetra(i,j,1,ioxygen)=ocetra(i,j,1,ioxygen)-oxflux/pddpo(i,j,1)
! Surface flux of gaseous nitrogen (same piston velocity as for O2)
       niflux=kwo2*dtbgc*(ocetra(i,j,1,igasnit)-anisa*(atn2/802000)) ! *ppao(i,j)/101300.
       ocetra(i,j,1,igasnit)=ocetra(i,j,1,igasnit)-niflux/pddpo(i,j,1)
! Surface flux of laughing gas (same piston velocity as for O2 and N2)
       n2oflux=kwo2*dtbgc*(ocetra(i,j,1,ian2o)-satn2o(i,j)) ! *ppao(i,j)/101300.
       ocetra(i,j,1,ian2o)=ocetra(i,j,1,ian2o)-n2oflux/pddpo(i,j,1)
#ifdef CFC
! Surface fluxes for CFC: eqn. (1a) in ocmip2 howto doc(hyc)
!     flux of CFC: downward direction (mol/m**2/s)
!      flx11=kw_11*(a_11*cfc11_atm(i,j)*ppair/p0-trc(i,j,1,1))
!      flx12=kw_12*(a_12*cfc12_atm(i,j)*ppair/p0-trc(i,j,1,2))
!      unit should be in [kmol cfc m-2]
!      unit of [cfc11_atm(i,j)*ppair/p0] should be in [pptv]
!      unit of [flx11-12] is in [kmol / m2]

! Surface flux of cfc11
      flx11=kw_11*dtbgc*                                                &
     & (a_11*atm_cfc11*ppao(i,j)*9.86923*1e-6-ocetra(i,j,1,icfc11))
      ocetra(i,j,1,icfc11)=ocetra(i,j,1,icfc11)+flx11/pddpo(i,j,1)
! Surface flux of cfc12
      flx12=kw_12*dtbgc*                                                &
     & (a_12*atm_cfc12*ppao(i,j)*9.86923*1e-6-ocetra(i,j,1,icfc12))
      ocetra(i,j,1,icfc12)=ocetra(i,j,1,icfc12)+flx12/pddpo(i,j,1)
! Surface flux of sf6
      flxsf=kw_sf*dtbgc*                                                &
     & (a_sf*atm_sf6*ppao(i,j)*9.86923*1e-6-ocetra(i,j,1,isf6))
      ocetra(i,j,1,isf6)=ocetra(i,j,1,isf6)+flxsf/pddpo(i,j,1)
#endif
#ifdef DIFFAT	 	 
!         atm(i,j,iatmo2)=atm(i,j,iatmo2) + (oxflux + 0.5*n2oflux)*contppm
!         atm(i,j,iatmn2)=atm(i,j,iatmn2) + (niflux + n2oflux)*contppm
       atmflx(i,j,iatmo2)=(oxflux + 0.5*n2oflux)*contppm
       atmflx(i,j,iatmn2)=(niflux + n2oflux) *contppm
#endif	 
       atmflx(i,j,iatmn2o)=n2oflux

! Surface flux of dms
       dmsflux = kwdms*dtbgc*ocetra(i,j,1,idms)  
       ocetra(i,j,1,idms)=ocetra(i,j,1,idms)-dmsflux/pddpo(i,j,1)

       aux2d_co2fxd(i,j)  = aux2d_co2fxd(i,j)  + fluxd 
       aux2d_co2fxu(i,j)  = aux2d_co2fxu(i,j)  + fluxu 
       aux2d_pco2(i,j)    = aux2d_pco2(i,j)    + pco2 
       aux2d_kwco2(i,j)   = aux2d_kwco2(i,j)   + kwco2*Kh*1e-6 ! JT replaced ak0 with Kh*1e-6
       aux2d_oxflux(i,j)  = aux2d_oxflux(i,j)  + oxflux 
       aux2d_niflux(i,j)  = aux2d_niflux(i,j)  + niflux
       aux2d_dmsflux(i,j) = aux2d_dmsflux(i,j) + dmsflux
       aux2d_dms(i,j)     = aux2d_dms(i,j)     + ocetra(i,j,1,idms) 
#ifdef CFC
       aux2d_cfc11(i,j)   = aux2d_cfc11(i,j)   + flx11
       aux2d_cfc12(i,j)   = aux2d_cfc12(i,j)   + flx12
       aux2d_sf6(i,j)     = aux2d_sf6(i,j)     + flxsf
#endif

      endif ! k==1

! Determine Omega Calcite et Aragonite
      omega = ( calcon * s / 35. ) * cc
!      OmegaA(i,j,k) = omega / Kspa
      OmegaC(i,j,k) = omega / Kspc

#ifdef __c_isotopes
       ocetra(i,j,k,isco214)=ocetra(i,j,k,isco214)*c14ret
       ocetra(i,j,k,idet14) =ocetra(i,j,k,idet14) *c14ret
       ocetra(i,j,k,icalc14)=ocetra(i,j,k,icalc14)*c14ret
#endif

      ! Save bottom level dissociation konstants for use in sediment module
      if( k==kbo(i,j) ) then

        k1b(i,j)  = K1
        k2b(i,j)  = K2
        kbb(i,j)  = Kb
        kwb(i,j)  = Kw
        kspb(i,j) = Kspc

      end if

      ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

! Accumulate diagnostic 2d variables   
      call accsrf(jco2fxd,aux2d_co2fxd,omask,0)
      call accsrf(jco2fxu,aux2d_co2fxu,omask,0)
      call accsrf(jpco2,aux2d_pco2,omask,0)
      call accsrf(jkwco2,aux2d_kwco2,omask,0)
      call accsrf(joxflux,aux2d_oxflux,omask,0)
      call accsrf(jniflux,aux2d_niflux,omask,0)
      call accsrf(jdmsflux,aux2d_dmsflux,omask,0)
      call accsrf(jdms,aux2d_dms,omask,0)
      call accsrf(jn2ofx,atmflx(1,1,iatmn2o),omask,0)
#ifdef CFC
      call accsrf(jcfc11fx,aux2d_cfc11,omask,0)
      call accsrf(jcfc12fx,aux2d_cfc12,omask,0)
      call accsrf(jsf6fx,aux2d_sf6,omask,0)
#endif

! Accumulate diagnostic layer variables 
      call acclyr(jomegac,OmegaC,pddpo,1)
      call acclyr(jn2o,ocetra(1,1,1,ian2o),pddpo,1) 
#ifdef CFC
      call acclyr(jcfc11,ocetra(1,1,1,icfc11),pddpo,1)
      call acclyr(jcfc12,ocetra(1,1,1,icfc12),pddpo,1)
      call acclyr(jsf6,ocetra(1,1,1,isf6),pddpo,1)
#endif

! Accumulate diagnostic level variables 
      IF (sum(jlvlomegac).NE.0) THEN
        DO k=1,kpke
          call bgczlv(pddpo,k,ind1,ind2,wghts)
          call acclvl(jlvlomegac,OmegaC,k,ind1,ind2,wghts)
          call acclvl(jlvln2o,ocetra(1,1,1,ian2o),k,ind1,ind2,wghts)          
        ENDDO 
      ENDIF
#ifdef CFC
      IF (SUM(jlvlcfc11+jlvlcfc12+jlvlsf6).NE.0) THEN
        DO k=1,kpke
          call bgczlv(pddpo,k,ind1,ind2,wghts)
          call acclvl(jlvlcfc11,ocetra(1,1,1,icfc11),k,ind1,ind2,wghts)
          call acclvl(jlvlcfc12,ocetra(1,1,1,icfc12),k,ind1,ind2,wghts)
          call acclvl(jlvlsf6,ocetra(1,1,1,isf6),k,ind1,ind2,wghts)
        ENDDO
      ENDIF
#endif


! -----------------------------------------------------------------
! Deep ocean processes

! C14 decay in the sediment
#ifdef __c_isotopes
        do k=1,ks
!$OMP PARALLEL DO  
        do j=1,kpje
        do i=1,kpie
        if(bolay(i,j).gt.0.) then
        sedlay(i,j,k,issso14)=sedlay(i,j,k,issso14)*c14ret
        sedlay(i,j,k,isssc14)=sedlay(i,j,k,isssc14)*c14ret
        powtra(i,j,k,ipowc14)=powtra(i,j,k,ipowc14)*c14ret
        endif
        enddo
        enddo
!$OMP END PARALLEL DO
        enddo
!        atm_c14=atm_c14*c14ret+c14prod*c14ret
#endif

!
! Dissolution of calcium
! Note : mixed layer (k=1) is assumed to be always supersaturated
!        (saturation really depends on temperature/DIC/alkalinity)
!
        DO 111 k=2,kpke
#ifdef __c_isotopes
!$OMP PARALLEL DO PRIVATE(supsat,undsa,dissol,r13,r14)
#else
!$OMP PARALLEL DO PRIVATE(supsat,undsa,dissol)
#endif
        DO 11 j=1,kpje
        DO 11 i=1,kpie
         IF(omask(i,j).GT.0.5.and.pddpo(i,j,k).GT.dp_min) THEN
! Dissolution of caco3 based on OmegaC:
!   omegaC=([CO3]*[Ca])/([CO3]sat*[Ca]sat)
!   Following Sarmiento and Gruber book, assumed that [Ca]=[Ca]sat
!   Thus, [CO3]sat=[CO3]/OmegaC. 
           supsat=co3(i,j,k)-co3(i,j,k)/OmegaC(i,j,k)
           undsa=MAX(0.,-supsat)
           dissol=MIN(undsa,0.05*ocetra(i,j,k,icalc))
#ifdef __c_isotopes
           r13=dissol*ocetra(i,j,k,icalc13)                        &
      &             /(ocetra(i,j,k,icalc)+1.e-25)
           r14=dissol*ocetra(i,j,k,icalc14)                        &
      &             /(ocetra(i,j,k,icalc)+1.e-25)
#endif
           ocetra(i,j,k,icalc)=ocetra(i,j,k,icalc)-dissol
           ocetra(i,j,k,ialkali)=ocetra(i,j,k,ialkali)+2.*dissol
           ocetra(i,j,k,isco212)=ocetra(i,j,k,isco212)+dissol
#ifdef __c_isotopes
           ocetra(i,j,k,icalc13)=ocetra(i,j,k,icalc13)-r13
           ocetra(i,j,k,isco213)=ocetra(i,j,k,isco213)+r13               ! remineralized calcite shells
           ocetra(i,j,k,icalc14)=ocetra(i,j,k,icalc14)-r14
           ocetra(i,j,k,isco214)=ocetra(i,j,k,isco214)+r14
#endif
         ENDIF
11    CONTINUE
!$OMP END PARALLEL DO
111   CONTINUE

!
      IF( kchck .EQ. 1) THEN
!         CALL CHCK_BGC(io_stdo_bgc,icyclibgc,                          &
!     &       'Check values of ocean tracer at exit from SBR CARCHM :', &
!     &       kpie,kpje,kpke,pddpo)

         DO  k=1,kpke
            DO  j=1,kpje
               DO  i=1,kpie
                  IF( hi(i,j,k) .LT. 0.0 ) THEN
                     WRITE(io_stdo_bgc,*)                              &
     &                   'CARCHM: invalid values of hi at i,j,k=',     &
     &                   i,j,k,hi(i,j,k)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      RETURN
      END

