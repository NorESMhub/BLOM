      SUBROUTINE CARCHM                                               &
     &           (kpie,kpje,kpke,pddpo,pdlxp,pdlyp,psao,ptho,psicomo, &
     &            pfu10,kplyear,kplmon,kplday,ndtdaybgc,ldtrunbgc,    &
     &            ptiestu,kmonlen,omask)

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
!     Tjiputra, BCCR 09.18.08
!     - modified all carbon chemistry formulations following the OCMIP protocols
!
!     Purpose
!     -------
!     Inorganic carbon cycle.
!
!     Method
!     -------
!     Surface fluxes of CO2 / N2O / dms
!     Dissolution of calcium
!     Note: O2 solubility in seawater is calculated in chemin
!     kchck=1 can be used to check max/min of bgc arrays on wet/dry cells.
!
!     To do :
!     -------
!
!     Surface fluxes should not go through sea ice !
!SL" the transfer velocity (piston velocity) is frequently expressed as
!SL  3 m/day or 12.5 cm/hour for typical moderate wind conditions
!SL
!SL  We set Vst = 3m/86400 s (to be modified in the coupled mode with
!SL  the actual scalar windspeed according to Wanninkhof).
!SL  Vst/Dzw(1) then is the time constant for the surface layer"
!SL
!SL  da ist dann klar, dasz noch mit dt multipliziert werden musz.
!
!
!     *CALL* *CARCHM(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,psao,ptho,psicomo,
!                                 pfu10,kplyear,kplmon,kplday,kmonlen)*
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st REAL of model grid.
!     *INTEGER* *kpje*    - 2nd REAL of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) REAL of model grid.
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd REAL) [m].
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st REAL) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (2nd REAL) [m].
!     *REAL*    *psao*    - salinity [psu.].
!     *REAL*    *ptho*    - potential temperature (3rd REAL).
!     *REAL*    *ppao*    - sea level pressure [Pascal].
!     *REAL*    *psicomo* - sea ice (2nd REAL).
!     *REAL*    *pfu10*   - forcing field wind speed (2nd REAL).
!
!     Externals
!     ---------
!     .
!
!**********************************************************************

      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
!      USE mo_timeser_bgc
      USE mo_control_bgc
      USE mo_bgcmean
      use mo_param1_bgc 

      implicit none
      INTEGER :: i,j,k,l,kpie,kpje,kpke
      INTEGER :: kplyear,kplmon,kplday,kmonlen,laumo1,ldtrunbgc,ndtdaybgc

      character rungen*3
      
      REAL psao(kpie,kpje,kpke)
      REAL pddpo(kpie,kpje,kpke)
      REAL pdlxp(kpie,kpje)
      REAL pdlyp(kpie,kpje)
      REAL omask(kpie,kpje)
      REAL psicomo(kpie,kpje)
      REAL pfu10(kpie,kpje)
      REAL ptho(kpie,kpje,kpke)

      REAL :: supsat, undsa, dissol
      REAL :: fluxd,fluxu
      REAL :: dddhhh,dadh,a,h,c,alk,t1,t2
      REAL :: akbi,ak2i,ak1i,dtja
      REAL :: kwco2,kwo2,kwdms
      REAL :: scco2,sco2,scdms
      REAL :: contppm, Xconvxa
      REAL :: oxflux,niflux,dmsflux,n2oflux
      REAL :: ato2, atn2, atco2,pco2
      REAL :: AHI, ABE,RMONLEN,RPLDAY
      REAL :: AK0,AK1,AK2,AKB,AKW,BT,oxysa,anisa
!Tjiputra update=for list of new variables==============================
      REAL :: t,tk,prb,s,invtk,dlogtk,is,is2,sqrtis,s2,sqrts,s15,scl
      REAL :: tmp,nkhwe74,Kh,K1,K2,Kb,K1p,K2p,K3p,Ksi,Kw,Ks1,Kf,Kspc,Kspa
      REAL :: zprb,zprb2,deltav,deltak,borat,sti,ft,tc,ta,sit,pt,ah1
      REAL :: hso4,hsi,hf,hpo4,ab,aw,ac,ah2o,ah2,erel,RHO2,cu,cb,cc
      REAL :: co2(kpie,kpje,kpke),hco3(kpie,kpje,kpke)
      REAL :: R,B,fco2,omega,OmegaA(kpie,kpje,kpke)
!      REAL :: R,B,fco2,omega,OmegaA(kpie,kpje,kpke),OmegaC(kpie,kpje,kpke)
      REAL :: &
        rgas = 83.131, &
        miss_val = 1.e+20, &
        eps = 5.e-5
      REAL ptiestu(kpie,kpje,kpke)
      REAL, DIMENSION(11) :: a0, a1, a2, b0, b1, b2, lnkpok0
      INTEGER :: js,niter,jit
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
!End of Tjiputra update=for list of variables===========================
#ifdef __c_isotopes
      REAL :: r13,r14,rat13,rat14
      REAL :: flux14d,flux14u,flux13d,flux13u
      REAL :: atc13,atc14
#endif

      REAL :: thickness
      INTEGER :: iter

#ifdef ANTC14
      REAL :: fantc14d,fantc14u
#endif

!      WRITE(*,*) 'CARCHM called with REDUCED :',                     &
!     &           kpie,kpje,kpke,pddpo(40,40,1),psao(40,40,1),        &
!     &           ptho(40,40,1),psicomo(40,40),pfu10(40,40),          &
!     &           kplyear,kplmon,kplday,kmonlen,ldtrunbgc
      laumo1=kplmon+1
      IF(laumo1.GT.12) laumo1=1
      
      RMONLEN=kmonlen
      RPLDAY=kplday
      AHI=RPLDAY/RMONLEN
      ABE=1.-AHI

      contppm=1/0.35e-3

!      if(ldtrunbgc.eq.1)then
!       if (kplyear.eq.1850.and.kplmon.eq.1.and.kplday.eq.1) then 
!        write(*,*)'jttesting annCO2Flx is set to 0.',kplyear,kplmon,kplday,ldtrunbgc
!        annCO2Flx=0.
!       else
!        open(79,file='../tmpOceCO2Flx'//rungen,form='FORMATTED')
!        read(79,'(E18.11)')annCO2flx
!        close(79)
!       endif
!      endif
!      if (kplmon.eq.1.and.kplday.eq.1.and.MOD(ldtrunbgc,18).eq.0) then 
!       write(*,*)'jttesting writing OCECO2FLX=',kplyear,annCO2flx*1e-15
!       open(79,file='../AnnOCECO2FLX'//rungen,form='FORMATTED',position='append')
!! Write annual CO2 flux from ocean in [Peta mol of CO2] unit
!       write(79,'(I4,E18.9)')kplyear-1,annCO2flx*1e-12
!       close(79)
!       annCO2flx=0.
!      endif
!
!     -----------------------------------------------------------------
!*         1. SET HALF PRECISION CONSTANTS
!             --- ---- --------- ---------
!
!
!Tjiputra update===========================================================================
      niter=20
!$OMP PARALLEL DO                                                       &
!$OMP&PRIVATE(t,tk,s,R,prb,invtk,dlogtk,is,is2,sqrtis,s2,sqrts,s15,scl, &
!$OMP&        tmp,nKhwe74,Kh,K1,K2,Kb,K1p,K2p,K3p,Ksi,Kw,Ks1,Kf,Kspc,   &
!$OMP&        Kspa,deltav,deltak,zprb,zprb2,lnkpok0,borat,sti,ft,tc,ta, &
!$OMP&        sit,pt,ah1,iter,jit,hso4,hf,hsi,hpo4,ab,aw,ac,ah2o,ah2,   &
!$OMP&        erel,cu,cb,cc,pco2,scco2,scdms,sco2,oxysa,anisa,Xconvxa,  &
!$OMP&        kwco2,kwdms,kwo2,atco2,ato2,atn2,fluxd,fluxu,oxflux,      &
!$OMP&        niflux,n2oflux,dmsflux,omega)
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
      IF(omask(i,j).gt.0.5.and.pddpo(i,j,k).GT.1.e-12) THEN
       t = ptho(i,j,k)
       tk = 273.15 + t
       s = psao(i,j,k)
       R = RHO2(MAX(25.,s),ptho(i,j,k),0.)
       prb = 1.025E-1*ptiestu(i,j,k) ! pressure in unit bars
       invtk = 1.0 / tk
       dlogtk = log(tk)
       is = 19.924 * s / ( 1000. - 1.005 * s )
       is2 = is * is
       sqrtis = SQRT(is)
       s2     = s * s
       sqrts  = SQRT(s)
       s15    = s**1.5
       scl    = s / 1.80655
!---------------------- Kh (K Henry) ---------------------------------
!              CO2(g) <-> CO2(aq.)
!              Kh      = [CO2]/ p CO2
!              Weiss (1974)   [mol/kg/atm]
!
       tmp = 9345.17 * invtk - 60.2409 + 23.3585 * log( tk/100. )
       nKhwe74 = tmp + s * ( 0.023517 - 0.00023656 * tk + 0.0047036e-4 * tk * tk )
       Kh      = exp( nKhwe74 )
!  K1 = [H][HCO3]/[H2CO3]   ; K2 = [H][CO3]/[HCO3]
! Millero p.664 (1995) using Mehrbach et al. data on seawater scale
       K1 = 10**( -1.0 * ( 3670.7 * invtk - 62.008 + 9.7944 * dlogtk - 0.0118 * s + 0.000116 * s2 ) )
       K2 = 10**( -1.0 * ( 1394.7 * invtk + 4.777 - 0.0184 * s + 0.000118 * s2 ) )
!  Kb = [H][BO2]/[HBO2] !
!  Millero p.669 (1995) using DATA from Dickson (1990)
       Kb = exp( ( -8966.90 - 2890.53 * sqrts - 77.942 * s + 1.728 * s15 - 0.0996 * s2 ) * invtk +     &
     &    ( 148.0248 + 137.1942 * sqrts + 1.62142 * s ) + ( -24.4344 - 25.085 * sqrts - 0.2474 * s ) * &
     &    dlogtk + 0.053105 * sqrts * tk )
!  K1p = [H][H2PO4]/[H3PO4] ; K2p = [H][HPO4]/[H2PO4] ; K3p = [H][PO4]/[HPO4]
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
!       if(i.eq.100.and.j.eq.40)write(*,*)'jtcarchm6=',prb,R,ptiestu(i,j,k)
!       if(i.eq.100.and.j.eq.40)write(*,*)'jtcarchm5=',Kw,invtk,dlogtk,sqrts,s,tk,prb
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
!  Calculate concentrations for borate, sulfate, and fluoride
! Uppstrom (1974)
      borat = 0.000232 * scl / 10.811
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
      tc  = ocetra(i,j,k,isco212) / R *1000.
      ta  = ocetra(i,j,k,ialkali) / R *1000.
      sit = ocetra(i,j,k,isilica) / R *1000.
      pt  = ocetra(i,j,k,iphosph) / R *1000.
      ah1  = 1.e-8
      iter = 0
      DO jit = 1,niter
       hso4 = sti / ( 1. + Ks1 / ( ah1 / ( 1. + sti / Ks1 ) ) )
       hf   = 1. / ( 1. + Kf / ah1 )
       hsi  = 1./ ( 1. + ah1 / Ksi )
       hpo4 = ( K1p * K2p * ( ah1 + 2. * K3p ) - ah1**3 ) /    & 
     &        ( ah1**3 + K1p * ah1**2 + K1p * K2p * ah1 + K1p * K2p * K3p )
       ab   = borat / ( 1. + ah1 / Kb )
       aw   = Kw / ah1 - ah1 / ( 1. + sti / Ks1 )
!       if(i.eq.100.and.j.eq.40)write(*,*)'jtcarchm4=',aw,Kw,ah1,sti,Ks1,scl 
       ac   = ta + hso4 - sit * hsi - ab - aw + ft * hf - pt * hpo4
       ah2o  = SQRT( ( tc - ac )**2 + 4. * ( ac * K2 / K1 ) * ( 2. * tc - ac ) )
       ah2  = 0.5 * K1 / ac *( ( tc - ac ) + ah2o )
       erel = ( ah2 - ah1 ) / ah2
       if (abs( erel ).ge.eps) then
        ah1 = ah2
        iter = iter + 1
       endif
      ENDDO
      if  (ah1.gt.0.) then 
!       ocetra(i,j,k,iph)=-log10(ah1)
       hi(i,j,k)=max(1.e-20,ah1)
      endif
! Determine CO2*, HCO3- and CO3-- concentrations (in mol/kg soln)
      cu = ( 2. * tc - ac ) / ( 2. + K1 / ah1 )
      cb = K1 * cu / ah1
      cc = K2 * cb / ah1
! Convert from mol/kg to mol/m^3 (xfac = 1.028 kg/L  x  1000 L/m^3)
      co2 (i,j,k)  = cu * R
      hco3(i,j,k)  = cb * R
      co3 (i,j,k)  = cc * R
!JT Convert from mol/m^3 to kmol/m^3
      co3 (i,j,k)  = co3(i,j,k)/1000.
!      ocetra(i,j,k,ico3)=co3(i,j,k)
! Determine CO2 pressure and fugacity (in micoatm)
! NOTE: equation below for pCO2 needs requires CO2 in mol/kg
      pco2 = cu * 1.e6 / Kh
!       if(i.eq.100.and.j.eq.40)write(*,*)'jtcarchm2=',pco2,cu,Kh,tc,ac,K1,ah1 

      if (k.eq.1) then
       scco2 = 2073.1 - 125.62*ptho(i,j,k) + 3.6276*ptho(i,j,k)**2  &
     &       - 0.043219*ptho(i,j,1)**3   ! Schmidt number
       scdms = 2674.0-147.12*ptho(i,j,1)+3.726*ptho(i,j,1)**2       &
     &       - 0.038*ptho(i,j,1)**3
       sco2  = 1638.0 - 81.83*ptho(i,j,1) + 1.483*ptho(i,j,1)**2    & 
     &       - 0.008004*ptho(i,j,1)**3  
      oxysa=AHI*CHEMCM(i,j,7,LAUMO1)+ABE*CHEMCM(i,j,7,kplmon)
      anisa=AHI*CHEMCM(i,j,8,LAUMO1)+ABE*CHEMCM(i,j,8,kplmon)
!  Compute the transfer (piston) velocity in m/s
       Xconvxa = 9.3611e-07
       kwco2 = (1-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660/scco2)**0.5
       kwdms = (1-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660/scdms)**0.5 
       kwo2  = (1-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660/sco2)**0.5 

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
!      if (Roc13.lt.0.8.or.Roc13.gt.1.2) then
!        write (io_stdo_bgc,*) 'Roc',i,j,Roc13,ocetra(i,j,1,isco213),ocetra(i,j,1,isco212)
!      endif
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
         atmflx(i,j,iatmco2)=(fluxu-fluxd)*contppm
#ifdef __c_isotopes
!         atm(i,j,iatmc13)=atm(i,j,iatmc13)+(flux13u-flux13d)*contppm
!         atm(i,j,iatmc14)=atm(i,j,iatmc14)+(flux14u-flux14d)*contppm
         atmflx(i,j,iatmc13)=(flux13u-flux13d)*contppm
         atmflx(i,j,iatmc14)=(flux14u-flux14d)*contppm
#endif
#endif /*DIFFAT*/	 

!       if(i.eq.100.and.j.eq.40)write(*,*)'jtcarchm1=',fluxu,fluxd,atco2,pco2,kwco2,Kh,dtbgc 
!       if(i.eq.100.and.j.eq.40)write(*,*)'jtcarchm0=',CHEMCM(i,j,5,kplmon) ,Kh
!       annCO2flx=annCO2flx+(fluxu-fluxd)*pdlyp(i,j)*pdlxp(i,j)
       ocetra(i,j,1,isco212)=ocetra(i,j,1,isco212)+(fluxd-fluxu)/pddpo(i,j,1)
#ifdef __c_isotopes
         ocetra(i,j,1,isco213)=                                     &
     &   ocetra(i,j,1,isco213)+(flux13d-flux13u)/pddpo(i,j,1)
         ocetra(i,j,1,isco214)=                                     &
     &   ocetra(i,j,1,isco214)+(flux14d-flux14u)/pddpo(i,j,1)
!      if (ocetra(i,j,1,isco213).le.0.0) then
!        write (io_stdo_bgc,*) 'isco213',ocetra(i,j,1,isco213),Roc13,pddpo(i,j,k)
!      endif
!      if (ocetra(i,j,1,isco214).le.0.0) then
!        write (io_stdo_bgc,*) 'isco214',ocetra(i,j,1,isco214),Roc14,pddpo(i,j,k)
!      endif
#endif

! Surface flux of oxygen
       oxflux=kwo2*dtbgc*(ocetra(i,j,1,ioxygen)-oxysa*(ato2/196800)) ! *ppao(i,j)/101300.
       ocetra(i,j,1,ioxygen)=ocetra(i,j,1,ioxygen)-oxflux/pddpo(i,j,1)
! Surface flux of gaseous nitrogen (same piston velocity as for O2)
       niflux=kwo2*dtbgc*(ocetra(i,j,1,igasnit)-anisa*(atn2/802000)) ! *ppao(i,j)/101300.
       ocetra(i,j,1,igasnit)=ocetra(i,j,1,igasnit)-niflux/pddpo(i,j,1)
! Surface flux of laughing gas (same piston velocity as for O2 and N2)
       n2oflux=kwo2*dtbgc*(ocetra(i,j,1,ian2o)-satn2o(i,j)) ! *ppao(i,j)/101300.
       ocetra(i,j,1,ian2o)=ocetra(i,j,1,ian2o)-n2oflux/pddpo(i,j,1)
#ifdef DIFFAT	 	 
!         atm(i,j,iatmo2)=atm(i,j,iatmo2) + (oxflux + 0.5*n2oflux)*contppm
!         atm(i,j,iatmn2)=atm(i,j,iatmn2) + (niflux + n2oflux)*contppm
       atmflx(i,j,iatmo2)=(oxflux + 0.5*n2oflux)*contppm
       atmflx(i,j,iatmn2)=(niflux + n2oflux) *contppm
#endif	 

! Surface flux of dms
       dmsflux = kwdms*dtbgc*ocetra(i,j,1,idms)  
       ocetra(i,j,1,idms)=ocetra(i,j,1,idms)-dmsflux/pddpo(i,j,1)

       bgcm2d(i,j,jco2fxd)=bgcm2d(i,j,jco2fxd)+fluxd
       bgcm2d(i,j,jco2fxu)=bgcm2d(i,j,jco2fxu)+fluxu
       bgcm2d(i,j,jpco2)  =bgcm2d(i,j,jpco2)  +pco2
       bgcm2d(i,j,jkwco2) =bgcm2d(i,j,jkwco2) +kwco2*Kh*1e-6 ! JT replaced ak0 with Kh*1e-6

       bgcm2d(i,j,joxflux)=bgcm2d(i,j,joxflux)+oxflux+0.5*n2oflux
       bgcm2d(i,j,jniflux)=bgcm2d(i,j,jniflux)+niflux+n2oflux
       bgcm2d(i,j,jdmsflux)=bgcm2d(i,j,jdmsflux) + dmsflux
       bgcm2d(i,j,jdms)=bgcm2d(i,j,jdms) + ocetra(i,j,1,idms)*pddpo(i,j,1)
      endif

!JT      B = ( -1636.75 + 12.0408 * tk - 0.0327957 * ( tk * tk ) + 0.0000316528 * ( tk * tk * tk ) ) * 1.e-6
!JT      fCO2(i,j,k) = pCO2(i,j,k) * exp( ( ( prb + 1. ) * 100000. ) * ( B + 2. * ( 57.7 - 0.118 * tk ) * 1.e-6 ) / ( 8.314 * tk ) )
! Determine Omega Calcite et Aragonite
      omega = ( 0.01028 * s / 35. ) * cc
!JT      OmegaA(i,j,k) = omega / Kspa
      OmegaC(i,j,k) = omega / Kspc
!JT Convert from mol/m^3 to kmol/m^3
      OmegaC(i,j,k) = OmegaC(i,j,k)/1000.
            bgcm3d(i,j,k,jomegac) = 			&
     &      bgcm3d(i,j,k,jomegac) + OmegaC(i,j,k)       *pddpo(i,j,k) 
!      ocetra(i,j,k,iomegac)=OmegaC(i,j,k)
!  Determine BetaD (from Seacarb software)
!JT      ahkb = ah1 + Kb;
!JT      ahk2 = ah1 + 2. * K2; 
!JT      cccb = 2. * cc + cb;
!JT      DD = -( ( -Kb * borat ) / ( ahkb * ahkb ) ) - ( -Kw / ( ah1 * ah1 ) ) + 1.;
!JT      A  = ( 2. * K2 * cccb + ah1 * ahk2 * DD ) / ( ahk2 * ahk2 );
!JT      B  = ( ( cccb * ah1 ) / ( ahk2 * K1 ) + ( ah1 / K1 ) * A );
!JT      C  = ( -K2 * cccb + K2 * ahk2 * DD ) / ( ahk2 * ahk2 );
!JT      PhiD (jk,jj,ji) = -1. / ( ah1 * log( 10.) * ( B + A + C ) );
!JT      BetaD(jk,jj,ji) = -ah1 * log( 10. ) * tc / cu * B * PhiD(jk,jj,ji);
#ifdef __c_isotopes
       ocetra(i,j,k,isco214)=ocetra(i,j,k,isco214)*c14ret
       ocetra(i,j,k,idet14) =ocetra(i,j,k,idet14) *c14ret
       ocetra(i,j,k,icalc14)=ocetra(i,j,k,icalc14)*c14ret
#endif

      ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
!End of Tjiputra update====================================================================

!      if(ldtrunbgc.eq.540)then
!       WRITE(*,*)'JT MONTH, CO2FLX (kmol/month)',kplmon,annCO2flx
!       open(79,file='../tmpOceCO2Flx'//rungen,form='FORMATTED')
!       write(79,'(E18.11)')annCO2flx
!       close(79)
!      endif


!     -----------------------------------------------------------------
!*        22. CHEMICAL CONSTANTS - DEEP OCEAN


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
!$OMP PARALLEL DO PRIVATE(supsat,undsa,dissol,r13,r14)
        DO 11 j=1,kpje
        DO 11 i=1,kpie
         IF(omask(i,j).GT.0.5.and.pddpo(i,j,k).GT.1.e-12) THEN
!update on dissolution of caco3 based on OmegaC concentration
           supsat=co3(i,j,k)-co3(i,j,k)/OmegaC(i,j,k)
! OLD           supsat=co3(i,j,k)-97.*aksp(i,j,k)
           undsa=MAX(0.,-supsat)
           dissol=MIN(undsa,0.02*ocetra(i,j,k,icalc))
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
!      if (ocetra(i,j,1,isco213).le.0.0) then
!      if (r13.lt.0.8.or.r13.gt.1.2) then
!        write (io_stdo_bgc,*) 'isco2132',ocetra(i,j,k,isco213),r13,pddpo(i,j,k)
!      endif
!      if (r14.lt.0.8.or.r14.gt.1.2) then
!        write (io_stdo_bgc,*) 'isco2142',ocetra(i,j,k,isco214),r14,pddpo(i,j,k)
!      endif
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

