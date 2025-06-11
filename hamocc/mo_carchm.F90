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

module mo_carchm

  use mo_kind, only: rp

  implicit none
  private

  public  :: carchm
  public  :: carchm_solve

  ! Maximum number of iterations for carbon chemistry solver
  integer, parameter :: niter=20

  ! Accuracy for test of convergence in carbon chemistry solver
  real(rp), parameter :: eps=5.e-5_rp

  ! Avoid division by zero
  real(rp), parameter :: eps_safe = epsilon(1._rp)

  ! Minimum and maximum SST/SSS set for carbon chemistry and gas exchange calculations
  real(rp), parameter :: temp_min = -1.0_rp
  real(rp), parameter :: temp_max = 40.0_rp
  real(rp), parameter :: saln_min =  5.0_rp
  real(rp), parameter :: saln_max = 40.0_rp

  ! Minimum and maximum H+ concentrations allowed in the iterative carbon chemistry solver
  ! The values set below (corresponding to pH 5 and 11) should never be reached in any resonable 
  ! oceanographic context
  real(rp), parameter :: ah_min = 1.0e-11_rp
  real(rp), parameter :: ah_max = 1.0e-5_rp

contains

  subroutine carchm(kpie,kpje,kpke,kbnd,pdlxp,pdlyp,pddpo,prho,pglat,omask,psicomo,ppao,pfu10,     &
       &            ptho,psao)

    !***********************************************************************************************
    ! Inorganic carbon cycle including dissolution of CaCO3 and air-sea gass exchange
    !
    !
    !  Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
    !
    !  Modified:
    !  S.Legutke,             *MPI-MaD, HH*    10.04.01
    !  - rename: ssso12(i,j,k)=sedlay(i,j,k,issso12 ) etc.; no equivalence statements
    !  - rename: powasi(i,j,k )=powtra(i,j,1,ipowasi) etc.; no equivalence statements
    !  - interfacing with ocean model
    !  J.Tjiputra,            *BCCR*           09.18.08
    !  - modified all carbon chemistry formulations following the OCMIP protocols
    !  J.Schwinger,           *GFI, UiB*       2013-04-22
    !  - Use density prho consistent with MICOM for conversion to mol/kg
    !  - Calculate solubility of O2 and N2 every timestep, consistent with
    !    what is done for carbon chemistry. Array chemcm not used any more.
    !  - Added J.Tjiputras code for cfc- and sf6-fluxes
    !  - Cautious code clean-up
    !  J.Schwinger,           *UNI-RESEARCH*   2017-08-30
    !   - Moved the accumulation of global fields for output to routine
    !     hamocc4bgc.
    !  A.Moree,          *GFI, Bergen*   2018-04-12
    !  - new version of carbon isotope code
    !  J.Tjiputra,       *Uni Research, Bergen*   2018-04-12
    !  - added preformed and saturated DIC tracers
    !  J.Schwinger,      *Uni Research, Bergen*   2018-04-12
    !  - moved accumulation of all output fields to seperate subroutine,
    !    related code-restructuring
    !  - dissolution of CaCO3 moved into main loop
    !  - added sediment bypass preprocessor option
    !***********************************************************************************************

    use mo_carbch,      only: atm,atmflx,co2fxd,co2fxu,co2star,co3,hi,keqb,kwco2sol,               &
                              ocetra,omegaa,omegac,fco2,pco2,xco2,pco2_gex,satn2o,satoxy,          &
                              kwco2a,co2sol,pn2om
    use mo_chemcon,     only: al1,al2,al3,al4,an0,an1,an2,an3,an4,an5,an6,                         &
                              bl1,bl2,bl3,calcon,ox0,ox1,ox2,ox3,ox4,ox5,ox6,                      &
                              oxyco,tzero,                                                         &
                              SV0_air,SV1_air,SV2_air,SV3_air,SV4_air,SD0_air,SD1_air,SD2_air,     &
                              SD3_air,Vb_nh3,M_nh3,kappa
    use mo_control_bgc, only: dtbgc,use_cisonew,use_natDIC,use_CFC,use_BROMO,                      &
                              use_cisonew,use_sedbypass,use_extNcycle
    use mo_param1_bgc,  only: ialkali,iatmo2,iatmco2,iatmdms,iatmn2,iatmn2o,ian2o,icalc,           &
                              idicsat,idms,igasnit,ioxygen,iphosph,isco212,isilica,                &
                              iatmf11,iatmf12,iatmsf6,icfc11,icfc12,isf6,                          &
                              iatmc13,iatmc14,icalc13,icalc14,idet14,idoc14,iphy14,                &
                              isco213,isco214,izoo14,safediv,                                      &
                              iatmnco2,inatalkali,inatcalc,inatsco212,                             &
                              ks,issso14,isssc14,ipowc14,                                          &
                              iatmbromo,ibromo,iatmnh3,ianh4
    use mo_param_bgc,   only: dremcalc,srfdic_min,c14dec,atm_co2_nat,atm_n2o
    use mo_vgrid,       only: dp_min,kmle,kbo,ptiestu
    use mo_carbch,      only: atm_cfc11_nh,atm_cfc11_sh,atm_cfc12_nh,atm_cfc12_sh,                 &
                              atm_sf6_nh,atm_sf6_sh,                                               &
                              co213fxd,co213fxu,co214fxd,co214fxu,                                 &
                              nathi,natco3,natpco2,natomegaa,natomegac,pnh3
    use mo_sedmnt,      only: sedlay,powtra,burial

    ! Arguments
    integer, intent(in) :: kpie                                              ! 1st dimension of model grid.
    integer, intent(in) :: kpje                                              ! 2nd dimension of model grid.
    integer, intent(in) :: kpke                                              ! 3rd (vertical) dimension of model grid.
    integer, intent(in) :: kbnd                                              ! nb of halo grid points
    real(rp),intent(in) :: pdlxp(kpie,kpje)                                  ! size of grid cell (1st dimension) [m].
    real(rp),intent(in) :: pdlyp(kpie,kpje)                                  ! size of grid cell (2nd dimension) [m].
    real(rp),intent(in) :: pddpo(kpie,kpje,kpke)                             ! size of grid cell (3rd dimension) [m].
    real(rp),intent(in) :: prho(kpie,kpje,kpke)                              ! density [g/cm^3].
    real(rp),intent(in) :: pglat(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)          ! latitude of grid cells [deg north].
    real(rp),intent(in) :: omask(kpie,kpje)                                  ! ocean mask.
    real(rp),intent(in) :: psicomo(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)        ! sea ice.
    real(rp),intent(in) :: ppao(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)           ! sea level pressure [pascal].
    real(rp),intent(in) :: pfu10(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)          ! forcing field wind speed.
    real(rp),intent(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)      ! potential temperature.
    real(rp),intent(in) :: psao(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)      ! salinity [psu].

    ! Local variables
    integer     :: i,j,k,l,js
    real(rp)    :: supsat, undsa, dissol
    real(rp)    :: rpp0,fluxd,fluxu
    real(rp)    :: kwco2,kwo2,kwn2,kwdms,kwn2o
    real(rp)    :: scco2,sco2,scn2,scdms,scn2o
    real(rp)    :: xconvxa
    real(rp)    :: oxflux,niflux,dmsflux,n2oflux
    real(rp)    :: ato2,atn2,atco2,atn2o
    real(rp)    :: oxy,ani,anisa
    real(rp)    :: rrho,t,t2,t3,t4,tk,tk100,prb,s,rs
    real(rp)    :: Kh0,K1,K2,Kb,K1p,K2p,K3p,Ksi,Kw,Ks1,Kf,Kspc,Kspa
    real(rp)    :: tc,ta,sit,pt,ah1,ac,cu,cu_sat,cb,cc,tc_sat
    real(rp)    :: Bvir,delta,fc,pH2O,omega
    real(rp)    :: atm_cfc11,atm_cfc12,atm_sf6,fact                    ! CFC
    real(rp)    :: sch_11,sch_12,sch_sf,kw_11,kw_12,kw_sf              ! CFC
    real(rp)    :: flx11,flx12,flxsf,a_11,a_12,a_sf                    ! CFC
    real(rp)    :: natcu,natcu_sat,natcb,natcc                         ! natDIC
    real(rp)    :: natfluxd,natfluxu,natomega                          ! natDIC
    real(rp)    :: natsupsat,natundsa,natdissol                        ! natDIC
    real(rp)    :: rco213,rco214                                       ! cisonew
    real(rp)    :: dissol13,dissol14,cu13,cu14,cu_sat13,cu_sat14       ! cisonew
    real(rp)    :: flux14d,flux14u,flux13d,flux13u                     ! cisonew
    real(rp)    :: atco213,atco214                                     ! cisonew
    real(rp)    :: frac_k,frac_aqg,frac_dicg                           ! cisonew
    real(rp)    :: flx_bromo,sch_bromo,kw_bromo,a_bromo,atbrf,Kb1,lsub ! BROMO
    ! extNcycle
    real(rp)    :: flx_nh3,sch_nh3_a,sch_nh3_w,kw_nh3,ka_nh3,atnh3,diff_nh3_a,diff_nh3_w,mu_air,mu_w,p_dbar,rho_air
    real(rp)    :: h_nh3,hstar_nh3,pKa_nh3,Kh_nh3,cD_wind,u_star


    ! set variables for diagnostic output to zero
    atmflx (:,:,:)=0._rp
    co2fxd   (:,:)=0._rp
    co2fxu   (:,:)=0._rp
    fco2     (:,:)=0._rp
    pco2     (:,:)=0._rp
    xco2     (:,:)=0._rp
    pco2_gex (:,:)=0._rp
    kwco2a   (:,:)=0._rp
    co2sol   (:,:)=0._rp
    kwco2sol (:,:)=0._rp
    co2star(:,:,:)=0._rp
    co3    (:,:,:)=0._rp
    satoxy (:,:,:)=0._rp
    omegaA (:,:,:)=0._rp
    omegaC (:,:,:)=0._rp
    pn2om    (:,:)=0._rp
    if (use_cisonew) then
      co213fxd (:,:)=0._rp
      co213fxu (:,:)=0._rp
      co214fxd (:,:)=0._rp
      co214fxu (:,:)=0._rp
    endif
    if (use_natDIC) then
      natpco2    (:,:)=0._rp
      natco3   (:,:,:)=0._rp
      natomegaA(:,:,:)=0._rp
      natomegaC(:,:,:)=0._rp
    endif
    if (use_extNcycle) then
      pnh3       (:,:)=0._rp
    endif

    ! -----------------------------------------------------------------
    ! Solve CO2 system
    ! -----------------------------------------------------------------
    !$OMP PARALLEL DO PRIVATE(t,t2,t3,t4,tk,tk100,s,rs,prb,Kh0,K1,K2                &
    !$OMP ,Kb,K1p,K2p,K3p,Ksi,Kw,Ks1,Kf,Kspc,Kspa,tc,ta,sit,pt,ah1,ac               &
    !$OMP ,cu,cb,cc,cu_sat,Bvir,delta,fc,pH2O,rpp0,scco2,scdms,sco2,oxy,ani,anisa   &
    !$OMP ,Xconvxa ,kwco2,kwdms,kwo2,atco2,ato2,atn2,atn2o,fluxd,fluxu,oxflux       &
    !$OMP ,tc_sat,niflux,n2oflux,dmsflux,omega,supsat,undsa,dissol                  &
    !$OMP ,sch_11,sch_12,sch_sf,kw_11,kw_12,kw_sf,a_11,a_12,a_sf,flx11              &
    !$OMP ,flx12,flxsf,atm_cfc11,atm_cfc12,atm_sf6,fact                             &
    !$OMP ,natcu,natcu_sat,natcb,natcc,natfluxd,natfluxu,natomega                   &
    !$OMP ,natsupsat,natundsa,natdissol,atco213,atco214,rco213,rco214,frac_aqg      &
    !$OMP ,frac_dicg,flux13d,flux13u,flux14d,flux14u,dissol13,dissol14,cu13,cu14    &
    !$OMP ,cu_sat13,cu_sat14,flx_bromo,sch_bromo,kw_bromo,a_bromo,atbrf,Kb1,lsub    &
    !$OMP ,flx_nh3,sch_nh3_a,sch_nh3_w,kw_nh3,ka_nh3,atnh3                          &
    !$OMP ,diff_nh3_a,diff_nh3_w,mu_air,mu_w,p_dbar,rho_air,h_nh3                   &
    !$OMP ,hstar_nh3,pKa_nh3,Kh_nh3,cD_wind,u_star,k,j,i,rrho,scn2,scn2o,kwn2,kwn2o)
    do k=1,kpke
      do j=1,kpje
        do i=1,kpie

          if (omask(i,j) > 0.5_rp .and. pddpo(i,j,k) > dp_min) then

            ! Carbon chemistry: Calculate equilibrium constants and solve for [H+] and
            ! carbonate alkalinity (ac)
            t    = min(temp_max,max(temp_min,ptho(i,j,k)))
            s    = min(saln_max,max(saln_min,psao(i,j,k)))
            tk   = t + tzero
            tk100= tk/100.0_rp

            rrho = prho(i,j,k)                   ! seawater density [g/cm3]
            prb  = ptiestu(i,j,k)*98060.0_rp*1.027e-6_rp ! pressure in unit bars, 98060 = onem

            tc   = ocetra(i,j,k,isco212) / rrho  ! convert to mol/kg
            ta   = ocetra(i,j,k,ialkali) / rrho
            sit  = ocetra(i,j,k,isilica) / rrho
            pt   = ocetra(i,j,k,iphosph) / rrho
            ah1  = hi(i,j,k)

            call CARCHM_KEQUI(t,s,prb,Kh0,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,Kspc,Kspa)

            call CARCHM_SOLVE(s,tc,ta,sit,pt,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,ah1,ac)

            hi(i,j,k)=ah1

            ! Determine CO2*, HCO3- and CO3-- concentrations (in mol/kg soln)
            cu = ( 2._rp * tc - ac ) / ( 2._rp + K1 / ah1 )
            cb = K1 * cu / ah1
            cc = K2 * cb / ah1
            co2star(i,j,k)=cu

            ! Carbonate ion concentration, convert from mol/kg to kmol/m^3
            co3(i,j,k)  = cc * rrho

            if (use_natDIC) then
              tc   = ocetra(i,j,k,inatsco212) / rrho  ! convert to mol/kg
              ta   = ocetra(i,j,k,inatalkali) / rrho
              ah1  = nathi(i,j,k)

              call CARCHM_SOLVE(s,tc,ta,sit,pt,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,ah1,ac)

              nathi(i,j,k)=ah1

              ! Determine natural CO2*, HCO3- and CO3-- concentrations (in mol/kg soln)
              natcu = ( 2._rp * tc - ac ) / ( 2._rp + K1 / ah1 )
              natcb = K1 * natcu / ah1
              natcc = K2 * natcb / ah1
              ! Natural carbonate ion concentration, convert from mol/kg to kmol/m^3
              natco3(i,j,k) = natcc * rrho
            endif

            ! solubility of O2 (Weiss, R.F. 1970, Deep-Sea Res., 17, 721-735) for moist air at
            ! 1 atm; multiplication with oxyco converts to kmol/m^3/atm, temp=[-1,40], saln=[0,40]
            oxy=ox0+ox1/tk100+ox2*alog(tk100)+ox3*tk100+s*(ox4+ox5*tk100+ox6*tk100**2)
            satoxy(i,j,k)=exp(oxy)*oxyco

            if (k==1) then

              ! -----------------------------------------------------------------
              ! Calculate Schmidt numbers, solubilities, and piston velocities
              ! -----------------------------------------------------------------
              t2 = t**2
              t3 = t**3
              t4 = t**4

              ! Schmidt numbers according to Wanninkhof (2014), Table 1, for temp=[-2,40]
              scco2 = 2116.8_rp - 136.25_rp*t + 4.7353_rp*t2 - 0.092307_rp*t3 + 0.0007555_rp *t4
              sco2  = 1920.4_rp - 135.6_rp *t + 5.2122_rp*t2 - 0.10939_rp *t3 + 0.00093777_rp*t4
              scn2  = 2304.8_rp - 162.75_rp*t + 6.2557_rp*t2 - 0.13129_rp *t3 + 0.0011255_rp *t4
              scdms = 2855.7_rp - 177.63_rp*t + 6.0438_rp*t2 - 0.11645_rp *t3 + 0.00094743_rp*t4
              scn2o = 2356.2_rp - 166.38_rp*t + 6.3952_rp*t2 - 0.13422_rp *t3 + 0.0011506_rp *t4
              if (use_CFC) then
                sch_11 = 3579.2_rp - 222.63_rp*t + 7.5749_rp*t2 - 0.14595_rp *t3 + 0.0011874_rp *t4
                sch_12 = 3828.1_rp - 249.86_rp*t + 8.7603_rp*t2 - 0.1716_rp  *t3 + 0.001408_rp  *t4
                sch_sf = 3177.5_rp - 200.57_rp*t + 6.8865_rp*t2 - 0.13335_rp *t3 + 0.0010877_rp *t4
              endif
              if (use_BROMO) then
                ! Stemmler et al. (2015; Biogeosciences) Eq. (9); Quack and Wallace
                ! (2003; GBC), temp=[-1,30]
                sch_bromo = 4662.8_rp - 319.45_rp*t + 9.9012_rp*t2 - 0.1159_rp*t3
              endif
              if (use_extNcycle) then
                ! Tsilingiris 2008 Eq.(45) for moist air (kg/m s)
                mu_air   = SV0_air + SV1_air*t + SV2_air*t2 + SV3_air*t3 + SV4_air*t4

                ! Tsinlingiris(44) moist air density (kg/m3)
                rho_air  = SD0_air + SD1_air*t + SD2_air*t2 + SD3_air*t3

                ! molecular viscosity of sea water
                ! (Matthaeus 1972, Richards 1998,assuming salinity s in per mille = ~PSU)
                p_dbar =  ppao(i,j)*1.0e-4_rp  ! sea level pressure (Pa *1e-5 -> bar *10-> dbar
                mu_w   = 1.79e-2_rp - 6.1299e-4_rp * t + 1.4467e-5_rp * t2 - 1.6826e-7_rp * t3     &
                       & - 1.8266e-7_rp * p_dbar + 9.8972e-12_rp * p_dbar*p_dbar + 2.4727e-5_rp * s&
                       & + s * (4.8429e-7_rp * t - 4.7172e-8_rp * t2 + 7.5986e-10_rp * t3)         &
                       & + s * (1.3817e-8_rp * t - 2.6363e-10_rp * t2)                             &
                       & - p_dbar*p_dbar * (6.3255e-13_rp * t - 1.2116e-14_rp * t2)
                mu_w   = mu_w * 0.1_rp ! conversion from g/(cm s) to kg/(m s)
 
                ! diffusion coeff in air (m2/s) Fuller 1966 / Johnson 2010
                ! division by pressure: ppao [Pa]; in Fuller, p is a factor for denominator [atm]
                diff_nh3_a =  1.0e-7_rp * tk**1.75_rp * M_nh3 / (ppao(i,j)/101325.0_rp)

                ! Johnson 2010 - (34) cm2/s -> m2/s (1e-8*1e-4=1e-12)
                ! closer to fit for Li & Gregory of: 9.874e-6*exp(2.644e-2*t)
                ! mu_w*1000: kg/(m s) -> cPoise as in Eq.(34) of Johnson 2010
                diff_nh3_w = 1.25e-12_rp*tk**1.52_rp *(mu_w*1000._rp)**(9.58_rp/Vb_nh3 -1.12_rp)   &
                           & *(Vb_nh3**(-0.19_rp) - 0.292_rp)

                ! Schmidt number air phase
                sch_nh3_a  = mu_air /(diff_nh3_a * rho_air)
                ! Schmidt number water phase
                sch_nh3_w  = mu_w   /(diff_nh3_w * rrho * 1000._rp)
              endif

              ! solubility of N2 (Weiss, R.F. 1970, Deep-Sea Res., 17, 721-735) for moist air
              ! at 1 atm; multiplication with oxyco converts to kmol/m^3/atm, temp=[-1,40], saln=[0,40]
              ani=an0+an1/tk100+an2*alog(tk100)+an3*tk100+s*(an4+an5*tk100+an6*tk100**2)
              anisa=exp(ani)*oxyco

              ! solubility of laughing gas  (Weiss and Price 1980, Marine Chemistry, 8, 347-359)
              ! for moist air at 1 atm in kmol/m^3/atm, temp=[-1,40], saln=[0,40]
              rs=al1+al2/tk100+al3*log(tk100)+al4*tk100**2+s*(bl1+bl2*tk100+bl3*tk100**2)
              satn2o(i,j)=exp(rs)

              if (use_CFC) then
                ! solubility of cfc11,12 (mol/(l*atm)) (Warner and Weiss 1985) and
                ! sf6 from eq. 6 of Bullister et al. (2002), temp=[-1,40], saln=[0,40]
                ! These are the alpha in (1b) of the ocmpic2 howto
                a_11 = exp(-229.9261_rp + 319.6552_rp*(100._rp/tk) + 119.4471_rp*log(tk100)        &
                     &     -1.39165_rp*(tk100)**2 + s*(-0.142382_rp + 0.091459_rp*(tk100)          &
                     &     -0.0157274_rp*(tk100)**2))
                a_12 = exp(-218.0971_rp + 298.9702_rp*(100._rp/tk) + 113.8049_rp*log(tk100)        &
                     &     -1.39165_rp*(tk100)**2 + s*(-0.143566_rp + 0.091015_rp*(tk100)          &
                     &     -0.0153924_rp*(tk100)**2))
                a_sf = exp(-80.0343_rp  + 117.232_rp *(100._rp/tk) +  29.5817_rp*log(tk100)        &
                     &     +s*(0.033518_rp-0.0373942_rp*(tk100)+0.00774862_rp*(tk100)**2))
                ! conversion from mol/(l * atm) to kmol/(m3 * pptv)
                a_11 = 1.0e-12_rp * a_11
                a_12 = 1.0e-12_rp * a_12
                a_sf = 1.0e-12_rp * a_sf
              endif
              if (use_BROMO) then
                !Henry's law constant [dimensionless] for Bromoform from Quack and Wallace (2003; GBC)
                ! temp=[0,20]
                a_bromo = exp(13.16_rp - 4973.0_rp*(1.0_rp/tk))
              endif
              if (use_extNcycle) then
                !Henry number for NH3 (Paulot et al. 2015, )
                h_nh3 =  (17.93_rp*tk/tzero * exp(4092._rp/tk - 9.7_rp))**(-1)
                ! Dissociation constant (Paulot et al. 2015, Bell 2007/2008)
                pKa_nh3 = 10.0423_rp - 3.15536e-2_rp*t + 3.071e-3_rp*s
                ! effective gas-over-liquid Henry constant (Paulot et al. 2015)
                hstar_nh3   = h_nh3/(1._rp + 10._rp**(log10(hi(i,j,k))+pKa_nh3))
              endif

              ! Transfer (piston) velocity kw according to Wanninkhof (2014), in units of ms-1
              Xconvxa = 6.97e-07_rp   ! Wanninkhof's a=0.251 converted from [cm hr-1]/[m s-1]^2 to [ms-1]/[m s-1]^2
              kwco2 = (1._rp-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660._rp/scco2)**0.5_rp
              kwo2  = (1._rp-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660._rp/sco2)**0.5_rp
              kwn2  = (1._rp-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660._rp/scn2)**0.5_rp
              kwdms = (1._rp-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660._rp/scdms)**0.5_rp
              kwn2o = (1._rp-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660._rp/scn2o)**0.5_rp
              if (use_CFC) then
                kw_11 = (1._rp-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660._rp/sch_11)**0.5_rp
                kw_12 = (1._rp-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660._rp/sch_12)**0.5_rp
                kw_sf = (1._rp-psicomo(i,j)) * Xconvxa * pfu10(i,j)**2*(660._rp/sch_sf)**0.5_rp
              endif
              if (use_BROMO) then
                ! Stemmler et al. (2015; Biogeosciences) Eq. (8)
                !  1.e-2/3600 = conversion from [cm hr-1]/[m s-1]^2 to [ms-1]/[m s-1]^2
                kw_bromo=(1._rp-psicomo(i,j)) * 1.e-2_rp/3600._rp *                       &
                     &   (0.222_rp*pfu10(i,j)**2+0.33_rp*pfu10(i,j))*(660._rp/sch_bromo)**0.5_rp
              endif
              if (use_extNcycle) then
                ! Paulot et al. 2015 / Johnson 2010
                ! friction velocity of wind (m/s)
                u_star  = pfu10(i,j)*sqrt(6.1e-4_rp + 6.3e-5_rp*pfu10(i,j))
                ! wind drag coeff (-)
                cD_wind = (u_star / (pfu10(i,j) + eps_safe))**2
                ! gas transfer velocity on gas phase side (m/s)
                ka_nh3  = 1.0e-3_rp + u_star/ (13.3_rp*sch_nh3_a + (eps_safe + cD_wind)**(-0.5_rp) &
                        &  - 5._rp + log(sch_nh3_a)/(2._rp*kappa))
                ! gas transfer velocity on liquid phase side (m/s) Nightingale 2000b - 3600*100: cm/h -> m/s
                kw_nh3  =  (0.24_rp*pfu10(i,j)**2 + 0.061_rp*pfu10(i,j))*sqrt(600._rp/sch_nh3_w)/360000._rp

                ! total effective gas transfer velocity (m/s)
                Kh_nh3  = (1._rp/(ka_nh3 + eps_safe) + hstar_nh3/(kw_nh3 + eps_safe))**(-1._rp)
                ! account for ice
                Kh_nh3  = (1._rp-psicomo(i,j)) * Kh_nh3
              endif


              ! -----------------------------------------------------------------
              ! Calculate and apply surface fluxes
              ! -----------------------------------------------------------------

              atco2 = atm(i,j,iatmco2)
              ato2  = atm(i,j,iatmo2)
              atn2  = atm(i,j,iatmn2)
              atn2o = atm(i,j,iatmn2o)
              if (use_cisonew) then
                atco213 = atm(i,j,iatmc13)
                atco214 = atm(i,j,iatmc14)
              endif
              if (use_BROMO) then
                atbrf = atm(i,j,iatmbromo)
              endif
              if (use_extNcycle) then
                atnh3  = atm(i,j,iatmnh3)
              endif

              ! Sea level pressure in atm. This is used in all surface flux calculations where
              ! atmospheric concentration is given as a mixing ratio (i.e. partial pressure = mixing ratio*SLP/P_0 [atm])
              rpp0  = ppao(i,j)/101325.0_rp

              ! calculate correction for non-ideality of CO2 (fugacity coefficient, Weiss and Price 1980)
              Bvir  = -1636.75_rp + 12.0408_rp*tk -0.0327957_rp*tk**2 + 0.0000316528_rp*tk**3   ! temp=[-12,47]
              delta = 57.7_rp-0.118_rp*tk                                                 ! temp=[0,40]
              fc    = exp( rpp0*(Bvir + 2.0_rp*delta)/(82.057_rp*tk) )

              ! calculate water vapor pressure [atm] as function of temperature and salinity (Weiss and Price 1980)
              ! temp=[0,40], saln=[0,40]
              pH2O  = exp(24.4543_rp - 67.4509_rp*(100.0_rp/tk) - 4.8489_rp*log(tk/100.0_rp) - 0.000544_rp*s)

              ! Calculate the CO2 concentration in equilibrium with atmospheric x' (atco2=mole fraction of CO2 in dry air [ppm])
              cu_sat = Kh0*atco2*1.0e-6_rp*(rpp0-pH2O)*fc

              fluxd = cu_sat*kwco2*dtbgc*rrho ! cu_sat and cu are in mol/kg. Multiply by rrho (g/cm^3) 
              fluxu = cu    *kwco2*dtbgc*rrho ! to get fluxes in kmol/m^2

              ! Set limit for CO2 outgassing to avoid negative DIC concentration
              fluxu = min(fluxu,fluxd-(srfdic_min - ocetra(i,j,k,isco212))*pddpo(i,j,1))

              if (use_natDIC) then
                natcu_sat = Kh0*atm_co2_nat*1.0e-6_rp*(rpp0-pH2O)*fc
                natfluxd  = natcu_sat*kwco2*dtbgc*rrho
                natfluxu  = natcu    *kwco2*dtbgc*rrho
                natfluxu  = min(natfluxu,natfluxd-(srfdic_min - ocetra(i,j,k,inatsco212))*pddpo(i,j,1))
              endif

              ! Calculate saturation DIC concentration in mixed layer
              call carchm_solve_DICsat(s,cu_sat,ta,sit,pt,Kh0,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,tc_sat)
              ocetra(i,j,1:kmle(i,j),idicsat) = tc_sat * rrho ! convert mol/kg to kmlo/m^3

              if (use_cisonew ) then
                ! Ocean-Atmosphere fluxes for carbon isotopes
                rco213=ocetra(i,j,1,isco213)/(ocetra(i,j,1,isco212)+safediv) ! Fraction DIC13 over total DIC
                rco214=ocetra(i,j,1,isco214)/(ocetra(i,j,1,isco212)+safediv) ! Fraction DIC14 over total DIC

                cu13     = cu * rco213 ! concentration of dissolved CO213
                cu14     = cu * rco214 ! concentration of dissolved CO214
                cu_sat13 = Kh0*atco213*1.0e-6_rp*(rpp0-pH2O)*fc
                cu_sat14 = Kh0*atco214*1.0e-6_rp*(rpp0-pH2O)*fc

                ! fractionation factors for 13C during air-sea gas exchange (Zhang et al. 1995, Orr et al. 2017)
                frac_k    = 0.99912_rp                                 ! Constant kinetic fractionation
                frac_aqg  = (0.0049_rp*t - 1.31_rp)/1000._rp + 1._rp            ! Gas dissolution fractionation
                frac_dicg = (0.0144_rp*t*(cc/(cc+cu+cb)) - 0.107_rp*t + 10.53_rp)/1000._rp + 1._rp !DIC to CO2 frac
                flux13d   = cu_sat13*kwco2*dtbgc*rrho*frac_aqg*frac_k
                flux13u   = cu13    *kwco2*dtbgc*rrho*frac_aqg*frac_k/frac_dicg
                flux14d   = cu_sat14*kwco2*dtbgc*rrho*(frac_aqg**2)*(frac_k**2)
                flux14u   = cu14    *kwco2*dtbgc*rrho*(frac_aqg**2)*(frac_k**2)/(frac_dicg**2)
              endif

              ! Update DIC
              ocetra(i,j,1,isco212)=ocetra(i,j,1,isco212)+(fluxd-fluxu)/pddpo(i,j,1)
              if (use_natDIC) then
                ocetra(i,j,1,inatsco212)=ocetra(i,j,1,inatsco212)+(natfluxd-natfluxu)/pddpo(i,j,1)
              endif
              if (use_cisonew) then
                ocetra(i,j,1,isco213)=ocetra(i,j,1,isco213)+(flux13d-flux13u)/pddpo(i,j,1)
                ocetra(i,j,1,isco214)=ocetra(i,j,1,isco214)+(flux14d-flux14u)/pddpo(i,j,1)
              endif

              ! Surface flux of oxygen
              oxflux=kwo2*dtbgc*(ocetra(i,j,1,ioxygen)-satoxy(i,j,1)*(ato2/196800.0_rp)*rpp0)
              ocetra(i,j,1,ioxygen)=ocetra(i,j,1,ioxygen)-oxflux/pddpo(i,j,1)
              ! Surface flux of gaseous nitrogen (same piston velocity as for O2)
              niflux=kwn2*dtbgc*(ocetra(i,j,1,igasnit)-anisa*(atn2/802000.0_rp)*rpp0)
              ocetra(i,j,1,igasnit)=ocetra(i,j,1,igasnit)-niflux/pddpo(i,j,1)
              ! Surface flux of laughing gas (same piston velocity as for O2 and N2)
              n2oflux=kwn2o*dtbgc*(ocetra(i,j,1,ian2o)-satn2o(i,j)*atn2o*1.0e-12_rp*rpp0)
              ! pN2O under moist air assumption at normal pressure
              pn2om(i,j) = 1.0e9_rp * ocetra(i,j,1,ian2o)/satn2o(i,j)
              ocetra(i,j,1,ian2o)=ocetra(i,j,1,ian2o)-n2oflux/pddpo(i,j,1)
              if (use_CFC) then
                ! Surface fluxes for CFC: eqn. (1a) in ocmip2 howto doc(hyc)
                !     flux of CFC: downward direction (mol/m**2/s)
                !      flx11=kw_11*(a_11*cfc11_atm(i,j)*ppair/p0-trc(i,j,1,1))
                !      flx12=kw_12*(a_12*cfc12_atm(i,j)*ppair/p0-trc(i,j,1,2))
                !      unit should be in [kmol cfc m-2]
                !      unit of [cfc11_atm(i,j)*ppair/p0] should be in [pptv]
                !      unit of [flx11-12] is in [kmol / m2]

                if (pglat(i,j) >= 10.0_rp) then
                  atm_cfc11=atm_cfc11_nh
                  atm_cfc12=atm_cfc12_nh
                  atm_sf6=atm_sf6_nh
                else if (pglat(i,j) <= -10.0_rp) then
                  atm_cfc11=atm_cfc11_sh
                  atm_cfc12=atm_cfc12_sh
                  atm_sf6=atm_sf6_sh
                else
                  fact=(pglat(i,j)-(-10.0_rp))/20.0_rp
                  atm_cfc11=fact*atm_cfc11_nh+(1.0_rp-fact)*atm_cfc11_sh
                  atm_cfc12=fact*atm_cfc12_nh+(1.0_rp-fact)*atm_cfc12_sh
                  atm_sf6  =fact*atm_sf6_nh  +(1.0_rp-fact)*atm_sf6_sh
                endif

                ! Surface flux of cfc11, cfc12, and sf6
                flx11=kw_11*dtbgc*(a_11*atm_cfc11*rpp0-ocetra(i,j,1,icfc11))
                flx12=kw_12*dtbgc*(a_12*atm_cfc12*rpp0-ocetra(i,j,1,icfc12))
                flxsf=kw_sf*dtbgc*(a_sf*atm_sf6  *rpp0-ocetra(i,j,1,isf6))
                ocetra(i,j,1,icfc11)=ocetra(i,j,1,icfc11)+flx11/pddpo(i,j,1)
                ocetra(i,j,1,icfc12)=ocetra(i,j,1,icfc12)+flx12/pddpo(i,j,1)
                ocetra(i,j,1,isf6)  =ocetra(i,j,1,isf6)  +flxsf/pddpo(i,j,1)
              endif

              ! Surface flux of dms
              ! Note that kwdms already has the open ocean fraction in the term
              dmsflux = kwdms*dtbgc*ocetra(i,j,1,idms)
              ocetra(i,j,1,idms) = ocetra(i,j,1,idms) - dmsflux/pddpo(i,j,1)

              if (use_BROMO) then
                ! Quack and Wallace (2003) eq. 1
                ! flux = kw*(Cw - Ca/H) ; kw[m s-1]; Cw[kmol m-3];
                ! Convert Ca(atbrf) from
                !  [pptv]    to [ppp]      by multiplying with 1e-12 (ppp = parts per part, dimensionless)
                !  [ppp]     to [mol L-1]  by multiplying with pressure[bar]/(SST[K]*R[L bar K-1 mol-1]); R=0,083
                !  [mol L-1] to [kmol m-3] by multiplying with 1

                flx_bromo = kw_bromo*dtbgc*(atbrf/a_bromo*1.0e-12_rp*ppao(i,j)*1.0e-5_rp/(tk*0.083_rp) - ocetra(i,j,1,ibromo))
                ocetra(i,j,1,ibromo) = ocetra(i,j,1,ibromo) + flx_bromo/pddpo(i,j,1)
              endif
              if (use_extNcycle) then
                ! surface flux NH3 - currently assumed atNH3 in pptv
                flx_nh3 =  Kh_nh3*dtbgc*(atnh3*1.0e-12_rp*ppao(i,j)*1.0e-5_rp/(tk*0.08314510_rp) - hstar_nh3*ocetra(i,j,1,ianh4))
                ocetra(i,j,1,ianh4) = ocetra(i,j,1,ianh4) + flx_nh3/pddpo(i,j,1)

                ! pNH3 in natm
                pnh3(i,j) =  hstar_nh3*ocetra(i,j,1,ianh4)  * 8.20573660809596e-5_rp * tk * 1.0e12_rp
              endif

              ! -----------------------------------------------------------------
              ! Save variables for output
              ! -----------------------------------------------------------------

              ! Save surface fluxes
              atmflx(i,j,iatmco2)=fluxu-fluxd
              atmflx(i,j,iatmo2)=oxflux
              atmflx(i,j,iatmn2)=niflux
              atmflx(i,j,iatmn2o)=n2oflux   ! positive to atmosphere [kmol N2O m-2 timestep-1]
              atmflx(i,j,iatmdms)=dmsflux   ! positive to atmosphere [kmol dms m-2 timestep-1]
              if (use_cisonew) then
                atmflx(i,j,iatmc13)=flux13u-flux13d
                atmflx(i,j,iatmc14)=flux14u-flux14d
              endif
              if (use_CFC) then
                atmflx(i,j,iatmf11)=flx11
                atmflx(i,j,iatmf12)=flx12
                atmflx(i,j,iatmsf6)=flxsf
              endif
              if (use_natDIC) then
                atmflx(i,j,iatmnco2)=natfluxu-natfluxd
              endif
              if (use_BROMO) then
                atmflx(i,j,iatmbromo)=-flx_bromo
              endif
              if (use_extNcycle) then
                atmflx(i,j,iatmnh3)=-flx_nh3 ! positive to atmosphere [kmol NH3 m-2 timestep-1]
              endif

              ! Save up- and downward components of carbon fluxes for output
              co2fxd(i,j)  = fluxd
              co2fxu(i,j)  = fluxu
              if (use_cisonew) then
                co213fxd(i,j)= flux13d
                co213fxu(i,j)= flux13u
                co214fxd(i,j)= flux14d
                co214fxu(i,j)= flux14u
              endif

              ! Save pCO2-related diagnostics for output
              fco2(i,j)     = cu * 1.0e6_rp / Kh0       ! Equilibrium CO2 fugacity at the air sea interface [micro atm]
              pco2(i,j)     = fco2(i,j)/fc           ! Equilibrium CO2 partial pressure at the air sea interface [micro atm]
              xco2(i,j)     = pco2(i,j)/(rpp0-pH2O)  ! Equilibrium CO2 dry air mixing raio [ppm]
              pco2_gex(i,j) = atco2*(rpp0-pH2O)      ! Actual CO2 partial pressure at the air sea interface [micro atm]
              if (use_natDIC) then
                natpco2(i,j) = natcu * 1.0e6_rp / Kh0 / fc
              endif

              ! Save product of piston velocity and solubility for output
              kwco2sol(i,j) = kwco2*Kh0*1.0e-6_rp ! m/s mol/kg/muatm
              kwco2a(i,j)   = kwco2            ! m/s (incl. ice fraction!)
              co2sol(i,j)   = Kh0*1.0e-6_rp       ! mol/kg/uatm

            endif ! k==1


            ! -----------------------------------------------------------------
            ! Deep ocean processes
            ! -----------------------------------------------------------------

            if (use_BROMO) then
              ! Degradation to hydrolysis (Eq. 2-4 of Stemmler et al., 2015)
              ! A1=1.23e17 mol min-1 => 2.05e12 kmol sec-1
              Kb1=2.05e12_rp*exp(-1.073e5_rp/(8.314_rp*tk))*dtbgc
              ocetra(i,j,k,ibromo)=ocetra(i,j,k,ibromo)*(1.0_rp-(Kb1*Kw/ah1))
              ! Degradation to halogen substitution (Eq. 5-6 of Stemmler et al., 2015)
              lsub=7.33e-10_rp*exp(1.250713e4_rp*(1.0_rp/298._rp-1.0_rp/tk))*dtbgc
              ocetra(i,j,k,ibromo)=ocetra(i,j,k,ibromo)*(1.0_rp-lsub)
            endif

            ! Determine Omega Calcite/Aragonite and dissolution of caco3 based on OmegaC:
            !   omegaC=([CO3]*[Ca])/([CO3]sat*[Ca]sat)
            !   Following Sarmiento and Gruber book, assumed that [Ca]=[Ca]sat
            !   Thus, [CO3]sat=[CO3]/OmegaC.
            omega = ( calcon * s / 35.0_rp ) * cc
            OmegaA(i,j,k) = omega / Kspa
            OmegaC(i,j,k) = omega / Kspc
            supsat=co3(i,j,k)-co3(i,j,k)/OmegaC(i,j,k)
            undsa=max(0.0_rp,-supsat)
            dissol=min(undsa,dremcalc*ocetra(i,j,k,icalc))
            if (use_natDIC) then
              natomega = ( calcon * s / 35.0_rp ) * natcc
              natOmegaA(i,j,k) = natomega / Kspa
              natOmegaC(i,j,k) = natomega / Kspc
              natsupsat=natco3(i,j,k)-natco3(i,j,k)/natOmegaC(i,j,k)
              natundsa=max(0.0_rp,-natsupsat)
              natdissol=min(natundsa,dremcalc*ocetra(i,j,k,inatcalc))
            endif
            if (use_cisonew) then
              dissol13=dissol*ocetra(i,j,k,icalc13)/(ocetra(i,j,k,icalc)+safediv)
              dissol14=dissol*ocetra(i,j,k,icalc14)/(ocetra(i,j,k,icalc)+safediv)
            endif
            ocetra(i,j,k,icalc)=ocetra(i,j,k,icalc)-dissol
            ocetra(i,j,k,ialkali)=ocetra(i,j,k,ialkali)+2._rp*dissol
            ocetra(i,j,k,isco212)=ocetra(i,j,k,isco212)+dissol
            if (use_natDIC) then
              ocetra(i,j,k,inatcalc)=ocetra(i,j,k,inatcalc)-natdissol
              ocetra(i,j,k,inatalkali)=ocetra(i,j,k,inatalkali)+2._rp*natdissol
              ocetra(i,j,k,inatsco212)=ocetra(i,j,k,inatsco212)+natdissol
            endif
            if (use_cisonew) then
              ocetra(i,j,k,icalc13)=ocetra(i,j,k,icalc13)-dissol13
              ocetra(i,j,k,isco213)=ocetra(i,j,k,isco213)+dissol13
              ocetra(i,j,k,icalc14)=ocetra(i,j,k,icalc14)-dissol14
              ocetra(i,j,k,isco214)=ocetra(i,j,k,isco214)+dissol14
            endif


            if (use_cisonew) then
              ! Decay of the ocean tracers that contain radioactive carbon 14C
              ocetra(i,j,k,isco214) = ocetra(i,j,k,isco214)*c14dec
              ocetra(i,j,k,idet14)  = ocetra(i,j,k,idet14) *c14dec
              ocetra(i,j,k,icalc14) = ocetra(i,j,k,icalc14)*c14dec
              ocetra(i,j,k,idoc14)  = ocetra(i,j,k,idoc14)*c14dec
              ocetra(i,j,k,iphy14)  = ocetra(i,j,k,iphy14)*c14dec
              ocetra(i,j,k,izoo14)  = ocetra(i,j,k,izoo14)*c14dec
            endif

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

            endif

          endif ! omask>0.5
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    ! C14 decay in the sediment (could be moved to sediment part)
    if (use_cisonew .and. .not. use_sedbypass) then
      do k=1,ks
        !$OMP PARALLEL DO PRIVATE(i,j)
        do j=1,kpje
          do i=1,kpie
            if(omask(i,j) > 0.5_rp) then
              sedlay(i,j,k,issso14)=sedlay(i,j,k,issso14)*c14dec
              sedlay(i,j,k,isssc14)=sedlay(i,j,k,isssc14)*c14dec
              powtra(i,j,k,ipowc14)=powtra(i,j,k,ipowc14)*c14dec
            endif
          enddo
        enddo
        !$OMP END PARALLEL DO
      enddo

      !$OMP PARALLEL DO PRIVATE(i)
      do j=1,kpje
        do i=1,kpie
          if(omask(i,j) > 0.5_rp) then
            burial(i,j,issso14) = burial(i,j,issso14)*c14dec
            burial(i,j,isssc14) = burial(i,j,isssc14)*c14dec
          endif
        enddo
      enddo
      !$OMP END PARALLEL DO
    endif ! end of use_cisonew and not use_sedbypass

  end subroutine carchm
  

  subroutine carchm_kequi(temp,saln,prb,Kh0,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,Kspc,Kspa)

    !***********************************************************************************************
    ! Calculate equilibrium constants for the carbonate system
    !***********************************************************************************************

    use mo_chemcon, only: tzero,rgas,bor1,bor2,salchl,ac1,ac2,ac3,ac4,bc1,bc2,bc3,ad1,ad2,ad3,     &
                          bd1,bd2,bd3,a0,a1,a2,b0,b1,b2

    ! Arguments
    real(rp),    intent(in)    :: temp   ! potential temperature [degr C].
    real(rp),    intent(in)    :: saln   ! salinity [psu].
    real(rp),    intent(in)    :: prb    ! pressure [bar].
    real(rp),    intent(out)   :: Kh0    ! equilibrium constant Kh0 = [CO2]/fCO2, moist air.
    real(rp),    intent(out)   :: K1     ! equilibrium constant K1  = [H][HCO3]/[H2CO3].
    real(rp),    intent(out)   :: K2     ! equilibrium constant K2  = [H][CO3]/[HCO3].
    real(rp),    intent(out)   :: Kb     ! equilibrium constant Kb  = [H][BO2]/[HBO2].
    real(rp),    intent(out)   :: Kw     ! equilibrium constant Kw  = [H][OH].
    real(rp),    intent(out)   :: Ks1    ! equilibrium constant Ks1 = [H][SO4]/[HSO4].
    real(rp),    intent(out)   :: Kf     ! equilibrium constant Kf  = [H][F]/[HF].
    real(rp),    intent(out)   :: Ksi    ! equilibrium constant Ksi = [H][SiO(OH)3]/[Si(OH)4].
    real(rp),    intent(out)   :: K1p    ! equilibrium constant K1p = [H][H2PO4]/[H3PO4].
    real(rp),    intent(out)   :: K2p    ! equilibrium constant K2p = [H][HPO4]/[H2PO4].
    real(rp),    intent(out)   :: K3p    ! equilibrium constant K3p = [H][PO4]/[HPO4].
    real(rp),    intent(out)   :: Kspc   ! equilibrium constant Kspc= [Ca2+]T [CO3]T.
    real(rp),    intent(out)   :: Kspa   ! equilibrium constant Kspa= [Ca2+]T [CO3]T.

    ! Local varibles
    integer                :: js
    real(rp)                   :: tk,tk100,invtk,dlogtk
    real(rp)                   :: t,s,is,is2,sqrtis,s15,s2,sqrts,scl
    real(rp)                   :: pK01,pK02
    real(rp)                   :: nKhwe74,deltav,deltak,zprb,zprb2
    real(rp)                   :: lnkpok0(11)

    t = min(temp_max,max(temp_min,temp))
    s = min(saln_max,max(saln_min,saln))
    tk = t + tzero
    tk100 = tk/100.0_rp
    invtk = 1.0_rp / tk
    dlogtk = log(tk)
    is = 19.924_rp * s / ( 1000._rp - 1.005_rp * s )
    is2 = is * is
    sqrtis = sqrt(is)
    s15    = s**1.5_rp
    s2     = s * s
    sqrts  = sqrt(s)
    scl    = s * salchl

    ! Note: temperature and salinity ranges given in comments indicate the valid range of the
    !       fit. This is usually taken to be the range for which measurements have been available.


    ! Kh0 = [CO2]/ fCO2 (fCO2 = fugacity of CO2 in air)
    ! Weiss (1974), note this does not include a correction for the effect of moist air at 
    ! air-sea interface [mol/kg/atm], temp=[-1,45], saln=[0,45]
    nKhwe74 = ad1+ad2/tk100+ad3*log(tk100)+s*(bd1+bd2*tk100+bd3*tk100**2)
    Kh0     = exp( nKhwe74 )
    ! K1 = [H][HCO3]/[H2CO3]   ; K2 = [H][CO3]/[HCO3]
    ! Waters et al. (2014) on total pH scale, temp=[0,45], saln=[0,45]
    pK01 = -126.34048_rp + 6320.813_rp*invtk + 19.568224_rp*dlogtk
    pK02 =  -90.18333_rp + 5143.692_rp*invtk + 14.613358_rp*dlogtk
    K1 = 10**( -1.0_rp*(pK01 + 13.568513_rp*sqrts + 0.031645_rp  *s - 5.3834e-5_rp*s2 - 539.2304_rp*sqrts*invtk - &
         &      5.635_rp  *s*invtk - 2.0901396_rp*sqrts*dlogtk) )
    K2 = 10**( -1.0_rp*(pK02 + 21.389248_rp*sqrts + 0.12452358_rp*s - 3.7447e-4_rp*s2 - 787.3736_rp*sqrts*invtk - &
         &     19.84233_rp*s*invtk - 3.3773006_rp*sqrts*dlogtk) )
    ! Kb = [H][BO2]/[HBO2] !
    ! Millero p.669 (1995) using data from Dickson (1990), temp=[0,45], saln=[5,45]
    Kb = exp( ( -8966.90_rp - 2890.53_rp  * sqrts - 77.942_rp  * s + 1.728_rp * s15 - 0.0996_rp * s2 ) * invtk +  &
         &    ( 148.0248_rp + 137.1942_rp * sqrts + 1.62142_rp * s ) +                                      &
         &    ( -24.4344_rp - 25.085_rp   * sqrts - 0.2474_rp  * s ) * dlogtk + 0.053105_rp * sqrts * tk )
    ! K1p = [H][H2PO4]/[H3PO4] ; K2p = [H][HPO4]/[H2PO4] ; K3p = [H][PO4]/[HPO4]
    ! DOE(1994) eq 7.2.20 with footnote using data from Millero (1974), temp=[5,30], saln=[5,40]
    K1p = exp( -4576.752_rp * invtk + 115.525_rp - 18.453_rp * dlogtk + ( -106.736_rp * invtk + 0.69171_rp ) *    &
         &     sqrts + ( -0.65643_rp * invtk - 0.01844_rp ) * s )
    K2p = exp( -8814.715_rp * invtk + 172.0883_rp - 27.927_rp * dlogtk + ( -160.340_rp * invtk + 1.3566_rp ) *    &
         &     sqrts + ( 0.37335_rp * invtk - 0.05778_rp ) *s );
    K3p = exp( -3070.75_rp * invtk - 18.141_rp + ( 17.27039_rp * invtk + 2.81197_rp ) * sqrts + ( -44.99486_rp *  &
         &     invtk - 0.09984_rp ) * s );
    ! Ksi = [H][SiO(OH)3]/[Si(OH)4]
    ! Millero p.671 (1995) using data from Yao and Millero (1995)
    Ksi = exp( -8904.2_rp * invtk + 117.385_rp - 19.334_rp * dlogtk + ( -458.79_rp * invtk + 3.5913_rp ) * sqrtis &
         & + ( 188.74_rp * invtk - 1.5998_rp) * is + ( -12.1652_rp * invtk + 0.07871_rp) * is2 +               &
         &     log(1.0_rp-0.001005_rp*s))
    ! Kw = [H][OH]
    ! Millero p.670 (1995) using composite data, temp=[0,45], saln=[0,45]
    Kw = exp( -13847.26_rp * invtk + 148.9652_rp - 23.6521_rp * dlogtk + ( 118.67_rp * invtk - 5.977_rp +         &
         &     1.0495_rp * dlogtk ) * sqrts - 0.01615_rp * s)
    ! Ks = [H][SO4]/[HSO4]
    ! Dickson (1990, J. chem. Thermodynamics 22, 113), temp=[0,45], saln=[5,45]
    Ks1 = exp( -4276.1_rp * invtk + 141.328_rp - 23.093_rp * dlogtk + ( -13856._rp * invtk + 324.57_rp - 47.986_rp * &
         &     dlogtk ) * sqrtis + ( 35474._rp * invtk - 771.54_rp + 114.723_rp * dlogtk ) * is - 2698._rp *   &
         &     invtk * is**1.5_rp + 1776._rp * invtk * is2 + log(1.0_rp - 0.001005_rp * s ) )
    ! Kf = [H][F]/[HF]
    ! Dickson and Riley (1979) -- change pH scale to total, temp=[0,45], saln=[0,45]
    Kf = exp( 1590.2_rp * invtk - 12.641_rp + 1.525_rp * sqrtis + log( 1.0_rp - 0.001005_rp * s ) + log( 1.0_rp + (  &
         &    0.1400_rp / 96.062_rp ) * scl / Ks1 ) )
    ! Kspc (calcite)
    ! apparent solubility product of calcite : Kspc = [Ca2+]T [CO32-]T
    ! where $[]_T$ refers to the equilibrium total (free + complexed) ion concentration.
    !          Mucci 1983 mol/kg-soln, temp=[0,40], saln=[5,35]
    Kspc = 10**( -171.9065_rp - 0.077993_rp * tk + 2839.319_rp / tk + 71.595_rp * log10( tk ) + ( - 0.77712_rp +  &
         &       0.0028426_rp * tk + 178.34_rp / tk ) * sqrts - 0.07711_rp * s + 0.0041249_rp * s15 );
    ! Kspa (aragonite)
    ! apparent solubility product of aragonite : Kspa = [Ca2+]T [CO32-]T
    ! where $[]_T$ refers to the equilibrium total (free + complexed) ion concentration.
    !          Mucci 1983 mol/kg-soln, temp=[0,40], saln=[5,35]
    Kspa = 10**( -171.945_rp - 0.077993_rp * tk + 2903.293_rp / tk  + 71.595_rp * log10( tk ) + ( -0.068393_rp +  &
         &       0.0017276_rp * tk + 88.135_rp / tk ) * sqrts - 0.10018_rp * s + 0.0059415_rp * s15 );


    !---------------------- Pressure effect on Ks (Millero, 95) --------------------
    ! index: K1 1, K2 2, Kb 3, Kw 4, Ks 5, Kf 6, Kspc 7, Kspa 8, K1p 9, K2p 10, K3p 11
    ! The fit below is valid for saln=35
    do js = 1,11
      deltav      = a0(js) + a1(js) * t + a2(js) * t * t
      deltak      = b0(js) + b1(js) * t + b2(js) * t * t
      zprb        = prb / ( rgas * tk )
      zprb2       = prb * zprb
      lnkpok0(js) = - ( deltav * zprb + 0.5_rp * deltak * zprb2 )
    enddo

    K1   = K1   * exp( lnkpok0(1)  )
    K2   = K2   * exp( lnkpok0(2)  )
    Kb   = Kb   * exp( lnkpok0(3)  )
    Kw   = Kw   * exp( lnkpok0(4)  )
    Ks1  = Ks1  * exp( lnkpok0(5)  )
    Kf   = Kf   * exp( lnkpok0(6)  )
    Kspc = Kspc * exp( lnkpok0(7)  )
    Kspa = Kspa * exp( lnkpok0(8)  )
    K1p  = K1p  * exp( lnkpok0(9)  )
    K2p  = K2p  * exp( lnkpok0(10) )
    K3p  = K3p  * exp( lnkpok0(11) )

  end subroutine carchm_kequi


  subroutine carchm_solve(saln,tc,ta,sit,pt,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,ah1,ac)

    !***********************************************************************************************
    ! Solve carbon chemistry.
    !***********************************************************************************************

    use mo_chemcon, only: bor1,bor2,salchl

    ! Arguments
    real(rp),    intent(in)    :: saln  ! salinity [psu].
    real(rp),    intent(in)    :: tc    ! total DIC concentraion [mol/kg].
    real(rp),    intent(in)    :: ta    ! total alkalinity [eq/kg].
    real(rp),    intent(in)    :: sit   ! silicate concentration [mol/kg].
    real(rp),    intent(in)    :: pt    ! phosphate concentration [mol/kg].
    real(rp),    intent(in)    :: K1    ! equilibrium constant K1  = [H][HCO3]/[H2CO3].
    real(rp),    intent(in)    :: K2    ! equilibrium constant K2  = [H][CO3]/[HCO3].
    real(rp),    intent(in)    :: Kb    ! equilibrium constant Kb  = [H][BO2]/[HBO2].
    real(rp),    intent(in)    :: Kw    ! equilibrium constant Kw  = [H][OH].
    real(rp),    intent(in)    :: Ks1   ! equilibrium constant Ks1 = [H][SO4]/[HSO4].
    real(rp),    intent(in)    :: Kf    ! equilibrium constant Kf  = [H][F]/[HF].
    real(rp),    intent(in)    :: Ksi   ! equilibrium constant Ksi = [H][SiO(OH)3]/[Si(OH)4].
    real(rp),    intent(in)    :: K1p   ! equilibrium constant K1p = [H][H2PO4]/[H3PO4].
    real(rp),    intent(in)    :: K2p   ! equilibrium constant K2p = [H][HPO4]/[H2PO4].
    real(rp),    intent(in)    :: K3p   ! equilibrium constant K3p = [H][PO4]/[HPO4].
    real(rp),    intent(inout) :: ah1   ! hydrogen ion concentration.
    real(rp),    intent(out)   :: ac    ! carbonate alkalinity.

    ! Local variables
    integer :: jit
    real(rp)    :: s,scl,borat,sti,ft
    real(rp)    :: hso4,hf,hsi,hpo4,ab,aw,ah2o,ah2,erel

    ! Calculate concentrations for borate, sulfate, and fluoride; see Dickson, A.G.,
    ! Sabine, C.L. and Christian, J.R. (Eds.) 2007. Guide to best practices
    ! for ocean CO2 measurements. PICES Special Publication 3, chapter 5 p. 10

    s = min(saln_max,max(saln_min,saln))
    scl = s * salchl
    borat = bor1 * scl * bor2      ! Uppstrom (1974)
    sti = 0.14_rp * scl / 96.062_rp      ! Morris & Riley (1966)
    ft = 0.000067_rp * scl / 18.9984_rp  ! Riley (1965)

    iflag: do jit = 1,niter
      hso4 = sti / ( 1._rp + Ks1 / ( ah1 / ( 1._rp + sti / Ks1 ) ) )
      hf   = 1._rp / ( 1._rp + Kf / ah1 )
      hsi  = 1._rp/ ( 1._rp + ah1 / Ksi )
      hpo4 = ( K1p * K2p * ( ah1 + 2._rp * K3p ) - ah1**3 ) /    &
           ( ah1**3 + K1p * ah1**2 + K1p * K2p * ah1 + K1p * K2p * K3p )
      ab   = borat / ( 1._rp + ah1 / Kb )
      aw   = Kw / ah1 - ah1 / ( 1._rp + sti / Ks1 )
      ac   = ta + hso4 - sit * hsi - ab - aw + ft * hf - pt * hpo4
      ah2o = sqrt( ( tc - ac )**2 + 4._rp * ( ac * K2 / K1 ) * ( 2._rp * tc - ac ) )
      ah2  = 0.5_rp * K1 / ac *( ( tc - ac ) + ah2o )
      erel = ( ah2 - ah1 ) / ah2
      if (abs( erel ) >= eps) then
        ! make sure [H+] stays between ah_min and ah_max
        ah1 = max(ah_min,min(ah2,ah_max))
      else
        ah1 = max(ah_min,min(ah2,ah_max))
        exit iflag
      endif
    enddo iflag

  end subroutine carchm_solve


  subroutine carchm_solve_dicsat(saln,co2_sat,ta,sit,pt,Kh0,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,tc_sat)

    !***********************************************************************************************
    ! Solve DICsat from TALK and h2co3_sat (saturation concentration of dissolved CO2).
    !
    ! J. Tjiputra,    *BCCR, Bergen*    25.01.17
    !***********************************************************************************************

    use mo_chemcon, only: bor1,bor2,salchl

    real(rp),    intent(in)  :: saln    ! salinity [psu].
    real(rp),    intent(in)  :: co2_sat ! saturation concentraion of dissoveld CO2 [mol/kg].
    real(rp),    intent(in)  :: ta      ! total alkalinity [eq/kg].
    real(rp),    intent(in)  :: sit     ! silicate concentration [mol/kg].
    real(rp),    intent(in)  :: pt      ! phosphate concentration [mol/kg].
    real(rp),    intent(in)  :: Kh0     ! equilibrium constant Kh0 = [H2CO3]/fCO2.
    real(rp),    intent(in)  :: K1      ! equilibrium constant K1  = [H][HCO3]/[H2CO3].
    real(rp),    intent(in)  :: K2      ! equilibrium constant K2  = [H][CO3]/[HCO3].
    real(rp),    intent(in)  :: Kb      ! equilibrium constant Kb  = [H][BO2]/[HBO2].
    real(rp),    intent(in)  :: Kw      ! equilibrium constant Kw  = [H][OH].
    real(rp),    intent(in)  :: Ks1     ! equilibrium constant Ks1 = [H][SO4]/[HSO4].
    real(rp),    intent(in)  :: Kf      ! equilibrium constant Kf  = [H][F]/[HF].
    real(rp),    intent(in)  :: Ksi     ! equilibrium constant Ksi = [H][SiO(OH)3]/[Si(OH)4].
    real(rp),    intent(in)  :: K1p     ! equilibrium constant K1p = [H][H2PO4]/[H3PO4].
    real(rp),    intent(in)  :: K2p     ! equilibrium constant K2p = [H][HPO4]/[H2PO4].
    real(rp),    intent(in)  :: K3p     ! equilibrium constant K3p = [H][PO4]/[HPO4].
    real(rp),    intent(out) :: tc_sat  ! saturated total DIC concentration [mol/kg].

    ! Local varibles
    integer :: jit
    real(rp)    :: s,scl,borat,sti,ft
    real(rp)    :: hso4,hf,hsi,hpo4,ab,aw,ah2o,ah2,erel
    real(rp)    :: hco3_sat,co3_sat,ah1,ac

    ! Calculate concentrations for borate, sulfate, and fluoride; see Dickson, A.G.,
    ! Sabine, C.L. and Christian, J.R. (Eds.) 2007. Guide to best practices
    ! for ocean CO2 measurements. PICES Special Publication 3, chapter 5 p. 10

    s = min(saln_max,max(saln_min,saln))
    scl = s * salchl
    borat = bor1 * scl * bor2      ! Uppstrom (1974)
    sti = 0.14_rp * scl / 96.062_rp      ! Morris & Riley (1966)
    ft = 0.000067_rp * scl / 18.9984_rp  ! Riley (1965)
    ah1=1.e-8_rp

    iflag: do jit = 1,niter
      hso4 = sti / ( 1._rp + Ks1 / ( ah1 / ( 1._rp + sti / Ks1 ) ) )
      hf   = 1._rp / ( 1._rp + Kf / ah1 )
      hsi  = 1._rp/ ( 1._rp + ah1 / Ksi )
      hpo4 = ( K1p * K2p * ( ah1 + 2._rp * K3p ) - ah1**3 ) /    &
           ( ah1**3 + K1p * ah1**2 + K1p * K2p * ah1 + K1p * K2p * K3p )
      ab   = borat / ( 1._rp + ah1 / Kb )
      aw   = Kw / ah1 - ah1 / ( 1._rp + sti / Ks1 )
      ac   = ta + hso4 - sit * hsi - ab - aw + ft * hf - pt * hpo4
      ah2o = sqrt((K1*co2_sat)**2 + 4._rp*ac*2._rp*K1*k2*co2_sat)
      ah2  = (K1*co2_sat + ah2o)/(2._rp*ac)
      erel = ( ah2 - ah1 ) / ah2
      if (abs( erel ) >= eps) then
        ! make sure [H+] stays between ah_min and ah_max
        ah1 = max(ah_min,min(ah2,ah_max))
      else
        ah1 = max(ah_min,min(ah2,ah_max))
        exit iflag
      endif
    enddo iflag

    hco3_sat  = K1 *      co2_sat / ah1
    co3_sat   = K1 * K2 * co2_sat / ah1**2
    tc_sat    = co2_sat + hco3_sat + co3_sat

  end subroutine carchm_solve_dicsat

end module mo_carchm
