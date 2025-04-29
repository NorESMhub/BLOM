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

module mo_ocprod

  implicit none
  private

  public :: ocprod

contains

  subroutine ocprod(kpie,kpje,kpke,kbnd,pdlxp,pdlyp,pddpo,omask,ptho,pi_ph,psao,ppao,prho)

    !***********************************************************************************************
    !  Biological production, remineralization and particle sinking.
    !
    !  Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
    !
    !  Modified
    !  S.Legutke,             *MPI-MaD, HH*    2010-04-01
    !  J.Schwinger,           *GFI, UiB*       2013-04-22
    !   - Corrected bug in light penetration formulation
    !   - Cautious code clean-up
    !  J.Tjiputra,            *UNI-RESEARCH*   2015-11-25
    !   - Implemented natural DIC/ALK/CALC
    !  I.Kriest,              *GEOMAR*         2016-08-11
    !   - Modified stoichiometry for denitrification (affects NO3, N2, Alk)
    !  J.Schwinger,           *UNI-RESEARCH*   2017-08-30
    !   - Removed split of the layer that only partly falls into the euphotic zone. Loops are
    !     now calculated over
    !      (1) layers that are completely or partly in the euphotoc zone
    !      (2) layers that do not lie within the euphotic zone.
    !   - Moved the accumulation of global fields for output to routine hamocc4bgc.
    !     The accumulation of local fields has been moved to the end of this routine.
    !  A.Moree,          *GFI, Bergen*   2018-04-12
    !   - new version of carbon isotope code
    !  J.Schwinger,      *Uni Research, Bergen*   2018-04-12
    !   - moved accumulation of all output fields to seperate subroutine,
    !     related code-restructuring
    !   - added sediment bypass preprocessor option and related code
    !  J.Schwinger,      *NORCE Climate, Bergen*   2020-05-29
    !   - Cleaned up parameter list
    !   - Dust deposition field now passed as an argument
    !  T. Bourgeois,     *NORCE climate, Bergen*   2025-04-14
    !  - implement R2OMIP protocol
    !***********************************************************************************************

    use mod_xc,           only: mnproc
    use mo_carbch,        only: ocetra,satoxy,hi,co2star
    use mo_sedmnt,        only: prcaca,produs,prorca,silpro,pror13,pror14,prca13,prca14
    use mo_param_bgc,     only: drempoc,drempoc_anaerob,bkox_drempoc,dremn2o,dremopal,dremsul,     &
                                dyphy,ecan,epsher,fesoly,                                          &
                                gammap,gammaz,grami,grazra,pi_alpha,phytomi,                       &
                                rcalc,rcar,rdn2o1,rdn2o2,rdnit0,rdnit1,rdnit2,                     &
                                relaxfe,remido,riron,rnit,rnoi,ro2ut,ropal,                        &
                                spemor,wcal_const,wdust_const,wopal_const,wpoc_const,              &
                                zinges,alar1,alar2,alar3,                                          &
                                alow1,alow2,alow3,calmax,cellmass,                                 &
                                cellsink,dustd1,dustd2,dustd3,dustsink,fractdim,                   &
                                fse,fsh,nmldmin,plower,pupper,sinkexp,stick,tmfac,                 &
                                tsfac,vsmall,zdis,wmin,wmax,wlin,rbro,                             &
                                dmsp1,dmsp2,dmsp3,dmsp4,dmsp5,dmsp6,dms_gamma,                     &
                                fbro1,fbro2,atten_f,atten_c,atten_uv,atten_w,bkopal,bkphy,bkzoo,   &
                                POM_remin_q10,POM_remin_Tref,opal_remin_q10,opal_remin_Tref,       &
                                bkphyanh4,bkphyano3,bkphosph,bkiron,ro2utammo,max_limiter,         &
                                O2thresh_aerob,O2thresh_hypoxic,NO3thresh_sulf,                    &
                                rcar_tdoclc,rcar_tdochc,rnit_tdoclc,rnit_tdochc,ro2ut_tdoclc,      &
                                ro2ut_tdochc,ro2utammo_tdoclc,ro2utammo_tdochc,rem_tdoclc,         &
                                rem_tdochc
    use mo_biomod,        only: bsiflx0100,bsiflx0500,bsiflx1000,bsiflx2000,bsiflx4000,bsiflx_bot, &
                                calflx0100,calflx0500,calflx1000,calflx2000,calflx4000,calflx_bot, &
                                carflx0100,carflx0500,carflx1000,carflx2000,carflx4000,carflx_bot, &
                                dustflx0100,dustflx0500,dustflx1000,dustflx2000,dustflx4000,       &
                                dustflx_bot,                                                       &
                                expoor,exposi,expoca,intdnit,intdms_bac,intdmsprod,intdms_uv,      &
                                intphosy,int_chbr3_prod,int_chbr3_uv,                              &
                                phosy3d,abs_oce,strahl,asize3d,wmass,wnumb,eps3d,phosy_NH4,        &
                                phosy_NO3,remin_aerob,remin_sulf
    use mo_param1_bgc,    only: ialkali,ian2o,iano3,icalc,idet,idms,idoc,ifdust,                   &
                                itdoc_lc,itdoc_hc,itdoc_lc13,itdoc_hc13,itdoc_lc14,itdoc_hc14,     &
                                igasnit,iiron,iopal,ioxygen,iphosph,iphy,isco212,                  &
                                isilica,izoo,iadust,inos,ibromo,                                   &
                                icalc13,icalc14,idet13,idet14,idoc13,idoc14,                       &
                                iphy13,iphy14,isco213,isco214,izoo13,izoo14,safediv,               &
                                inatalkali,inatcalc,inatsco212,ianh4
    use mo_control_bgc,   only: dtb,io_stdo_bgc,with_dmsph,                                        &
                                use_BROMO,use_AGG,use_PBGC_OCNP_TIMESTEP,use_FB_BGC_OCE,           &
                                use_AGG,use_cisonew,use_natDIC, use_WLIN,use_sedbypass,use_M4AGO,  &
                                use_extNcycle,lkwrbioz_off,lTO2depremin,use_river2omip
    use mo_vgrid,         only: dp_min,dp_min_sink,k0100,k0500,k1000,k2000,k4000,kwrbioz,ptiestu
    use mo_vgrid,         only: kmle
    use mo_clim_swa,      only: swa_clim
    use mo_inventory_bgc, only: inventory_bgc
    use mo_ihamocc4m4ago, only: ihamocc_mean_aggregate_sinking_speed,ws_agg
    use mo_extNwatercol,  only: nitrification,denit_NO3_to_NO2,anammox,denit_dnra,extN_inv_check


    ! Arguments
    integer, intent(in) :: kpie                                         ! 1st dimension of model grid.
    integer, intent(in) :: kpje                                         ! 2nd dimension of model grid.
    integer, intent(in) :: kpke                                         ! 3rd (vertical) dimension of model grid.
    integer, intent(in) :: kbnd                                         ! nb of halo grid points
    real,    intent(in) :: pdlxp(kpie,kpje)                             ! size of grid cell (1st dimension) [m].
    real,    intent(in) :: pdlyp(kpie,kpje)                             ! size of grid cell (2nd dimension) [m].
    real,    intent(in) :: pddpo(kpie,kpje,kpke)                        ! size of grid cell (3rd dimension) [m].
    real,    intent(in) :: omask(kpie,kpje)                             ! land/ocean mask (1=ocean)
    real,    intent(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) ! potential temperature [deg C].
    real,    intent(in) :: pi_ph(kpie,kpje)
    real,    intent(in) :: psao(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) ! salinity [psu].
    real,    intent(in) :: ppao(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)      ! sea level pressure [Pascal].
    real,    intent(in) :: prho(kpie,kpje,kpke)                         ! density [g/cm^3].

    ! Local variables
    integer, parameter :: nsinkmax = 12
    integer :: i,j,k,l
    integer :: is,kdonor
    real :: abs_bgc(kpie,kpje,kpke)
    real :: tco(nsinkmax),tcn(nsinkmax),q(nsinkmax)
    real :: atten,avphy,avanut,avanfe,pho,xa,xn,ya,yn,phosy
    real :: avgra,grazing,avsil,avdic,graton
    real :: gratpoc,grawa,bacfra,phymor,zoomor,excdoc,exud
    real :: export, delsil, delcar, sterph, sterzo, remin
    real :: docrem,opalrem,remin2o,aou,refra,pocrem,phyrem,tdoclc_rem,tdochc_rem
    real :: zoothresh,phythresh
    real :: temp,temfa,phofa                  ! temperature and irradiation factor for photosynthesis
    real :: absorption,absorption_uv
    real :: dmsprod,dms_bac,dms_uv,dms_ph
    real :: dtr,dz
    real :: wpocd,wcald,wopald,wdustd,dagg
    real :: wcal,wdust,wopal,wpoc
    real :: o2lim ! O2 limitation of ammonification (POC remin)
    ! sedbypass
    real :: florca,flcaca,flsil
    ! cisonew
    real :: phygrowth
    real :: phosy13,phosy14
    real :: growth_co2
    real :: bifr13,bifr14,bifr13_perm
    real :: grazing13,grazing14
    real :: graton13,graton14
    real :: gratpoc13,gratpoc14
    real :: bacfra13,bacfra14
    real :: phymor13,phymor14
    real :: grawa13,grawa14
    real :: zoomor13,zoomor14
    real :: excdoc13,excdoc14
    real :: exud13,exud14
    real :: export13,export14
    real :: delcar13,delcar14
    real :: dtr13,dtr14
    real :: sterph13,sterph14
    real :: sterzo13,sterzo14
    real :: pocrem13,pocrem14
    real :: docrem13,docrem14
    real :: tdoclc_rem13,tdochc_rem13,tdoclc_rem14,tdochc_rem14
    real :: phyrem13,phyrem14
    real :: rem13,rem14
    real :: rco213,rco214,rdoc13,rdoc14,rdet13,rdet14
    real :: rtdoclc13,rtdochc13,rtdoclc14,rtdochc14
    real :: rphy13,rphy14,rzoo13,rzoo14
    ! sedbypass
    real :: flor13,flor14,flca13,flca14
    ! AGG
    real :: aggregate(kpie,kpje,kpke)
    real :: dustagg(kpie,kpje,kpke)
    real :: avmass, avnos, anosloss
    real :: zmornos, eps, e1,e2,e3,e4,es1,es3
    real :: TopM,TopF, snow,fshear,sagg1,sagg2,sagg4
    real :: sett_agg,shear_agg,effsti,dfirst,dshagg,dsett
    real :: wnos,wnosd
    ! BROMO
    real :: bro_beta,bro_uv
    real :: abs_uv(kpie,kpje,kpke)
    ! extNcycle
    character(len=:), allocatable :: inv_message
    real :: ano3up_inh,nutlim,anh4lim,nlim,grlim,nh4uptfrac

    ! set variables for diagnostic output to zero
    expoor     (:,:) = 0.
    expoca     (:,:) = 0.
    exposi     (:,:) = 0.
    carflx0100 (:,:) = 0.
    carflx0500 (:,:) = 0.
    carflx1000 (:,:) = 0.
    carflx2000 (:,:) = 0.
    carflx4000 (:,:) = 0.
    carflx_bot (:,:) = 0.
    bsiflx0100 (:,:) = 0.
    bsiflx0500 (:,:) = 0.
    bsiflx1000 (:,:) = 0.
    bsiflx2000 (:,:) = 0.
    bsiflx4000 (:,:) = 0.
    bsiflx_bot (:,:) = 0.
    calflx0100 (:,:) = 0.
    calflx0500 (:,:) = 0.
    calflx1000 (:,:) = 0.
    calflx2000 (:,:) = 0.
    calflx4000 (:,:) = 0.
    calflx_bot (:,:) = 0.
    dustflx0100(:,:) = 0.
    dustflx0500(:,:) = 0.
    dustflx1000(:,:) = 0.
    dustflx2000(:,:) = 0.
    dustflx4000(:,:) = 0.
    dustflx_bot(:,:) = 0.
    intdnit    (:,:) = 0.
    intphosy   (:,:) = 0.
    intdmsprod (:,:) = 0.
    intdms_bac (:,:) = 0.
    intdms_uv  (:,:) = 0.
    phosy3d  (:,:,:) = 0.

    if (use_BROMO) then
      int_chbr3_uv  (:,:) = 0.
      int_chbr3_prod(:,:) = 0.
    end if
    if (use_AGG) then
      eps3d(:,:,:)    = 0.
      asize3d(:,:,:)  = 0.
    endif
    if (use_extNcycle) then
      phosy_NH4(:,:,:)   = 0.
      phosy_NO3(:,:,:)   = 0.
      remin_aerob(:,:,:) = 0.
      remin_sulf(:,:,:)  = 0.
    endif

    if (use_PBGC_OCNP_TIMESTEP) then
      if (mnproc == 1) then
        write(io_stdo_bgc,*)' '
        write(io_stdo_bgc,*)'beginning of OCRPOD '
      endif
      call inventory_bgc(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
    endif

    ! Calculate swr absorption by water and phytoplankton

    abs_bgc(:,:,:) = 0.
    if (use_BROMO) then
      abs_uv(:,:,:) = 0.
    endif
    if (use_FB_BGC_OCE) then
      abs_oce(:,:,:) = 0.
      abs_oce(:,:,1) = 1.
    endif

    !$OMP PARALLEL DO PRIVATE(i,k,absorption,absorption_uv,atten,dz)
    do j = 1,kpje
      do i = 1,kpie

        if(omask(i,j) > 0.5) then

          absorption    = 1.
          absorption_uv = 1.

          vloop: do k = 1,merge(kpke,kwrbioz(i,j),lkwrbioz_off)

            if(pddpo(i,j,k) > 0.0) then

              dz = pddpo(i,j,k)

              ! Average light intensity in layer k
              atten = atten_w + atten_c * max(0.,ocetra(i,j,k,iphy))
              abs_bgc(i,j,k) = ((absorption/atten)*      (1.-exp(-atten*dz)))/dz
              if (use_BROMO) then
                abs_uv(i,j,k)  = ((absorption_uv/atten_uv)*(1.-exp(-atten_uv*dz)))/dz
              endif
              if (use_FB_BGC_OCE) then
                abs_oce(i,j,k) = abs_oce(i,j,k) * absorption
                if (k == 2) then
                  abs_oce(i,j,2) = atten_f * absorption
                endif
              endif

              ! Radiation intensity I_0 at the top of next layer
              absorption    = absorption    * exp(-atten*dz)
              absorption_uv = absorption_uv * exp(-atten_uv*dz)

            endif
          enddo vloop

        endif ! omask > 0.5

      enddo
    enddo
    !$OMP END PARALLEL DO

    if (use_M4AGO) then
      ! even though we loose detritus, etc. we call the calculation for settling velocity by M4AGO here
      ! to enable further future development... - assuming that the operator splitting decently functions
      call ihamocc_mean_aggregate_sinking_speed(kpie, kpje, kpke, kbnd, pddpo, omask, ptho, psao, ppao, prho)
    endif

    !$OMP PARALLEL DO PRIVATE(avphy,avgra,avsil,avanut,avanfe,pho,xa,xn   &
    !$OMP  ,phosy,ya,yn,grazing,graton,gratpoc,grawa,bacfra,phymor        &
    !$OMP  ,zoomor,excdoc,exud,export,delsil,delcar,dmsprod               &
    !$OMP  ,dms_bac,dms_uv,dtr,phofa,temfa,zoothresh,dms_ph,dz,opalrem    &
    !$OMP  ,avmass,avnos,zmornos,tdoclc_rem,tdochc_rem                    &
    !$OMP  ,tdoclc_rem13,tdochc_rem13,tdoclc_rem14,tdochc_rem14           &
    !$OMP  ,rco213,rco214,rphy13,rphy14,rzoo13,rzoo14,grazing13,grazing14 &
    !$OMP  ,graton13,graton14,gratpoc13,gratpoc14,grawa13,grawa14         &
    !$OMP  ,phosy13,phosy14,bacfra13,bacfra14,phymor13,phymor14,zoomor13  &
    !$OMP  ,zoomor14,excdoc13,excdoc14,exud13,exud14,export13,export14    &
    !$OMP  ,delcar13,delcar14,dtr13,dtr14,bifr13,bifr14,bifr13_perm       &
    !$OMP  ,growth_co2,phygrowth                                          &
    !$OMP  ,bro_beta,bro_uv                                               &
    !$OMP  ,ano3up_inh,nutlim,anh4lim,nlim,grlim,nh4uptfrac               &
    !$OMP  ,i,k)

    loop1: do j = 1,kpje
      do i = 1,kpie
        do k = 1,merge(kpke,kwrbioz(i,j),lkwrbioz_off)

          if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then


            if (use_AGG) then
              avmass = ocetra(i,j,k,iphy) + ocetra(i,j,k,idet)
            endif

            temp = min(40.,max(-3.,ptho(i,j,k)))
            phofa = pi_alpha * strahl(i,j) * abs_bgc(i,j,k)
            temfa = 0.6 * 1.066**temp
            !taylor:     temfa= 0.6*(1. + 0.0639*ptho(i,j,k) *               &
            !    &              (1. + 0.0639*ptho(i,j,k)/2. * (1. + 0.0639*ptho(i,j,k)/3.)))
            pho = dtb * phofa * temfa / sqrt(phofa**2 + temfa**2)

            avphy = max(phytomi,ocetra(i,j,k,iphy))                   ! 'available' phytoplankton
            avgra = max(grami,ocetra(i,j,k,izoo))                     ! 'available' zooplankton
            avsil = max(0.,ocetra(i,j,k,isilica))
            avdic = max(0.,ocetra(i,j,k,isco212))
            if (use_extNcycle)then
              ano3up_inh = bkphyanh4/(bkphyanh4 + ocetra(i,j,k,ianh4)) ! inhibition of NO3 uptake
              nutlim     = min(ocetra(i,j,k,iphosph)/(ocetra(i,j,k,iphosph)+bkphosph),                              &
                               ocetra(i,j,k,iiron)/(ocetra(i,j,k,iiron)+bkiron))
              anh4lim    = ocetra(i,j,k,ianh4)/(ocetra(i,j,k,ianh4) + bkphyanh4)
              nlim       = ano3up_inh*ocetra(i,j,k,iano3)/(ocetra(i,j,k,iano3) +  bkphyano3) + anh4lim
              grlim      = min(nutlim,nlim) ! growth limitation

              nh4uptfrac = anh4lim/(nlim+epsilon(1.))
              ! re-check avnut - can sum N avail exceed indiv. contrib?
              avanut     = max(0.,min(ocetra(i,j,k,iphosph), ocetra(i,j,k,iiron)/riron,                              &
                         &        rnoi*((1.-nh4uptfrac)*ocetra(i,j,k,iano3) + nh4uptfrac*ocetra(i,j,k,ianh4))))

              xn         = avphy/(1. - pho*grlim)       ! phytoplankton growth
              phosy      = max(0.,min(xn-avphy,max_limiter*avanut)) ! limit PP growth to available nutr.
            else
              avanut = max(0.,min(ocetra(i,j,k,iphosph),rnoi*ocetra(i,j,k,iano3)))
              avanfe = max(0.,min(avanut,ocetra(i,j,k,iiron)/riron))
              xa = avanfe
              xn = xa/(1.+pho*avphy/(xa+bkphy))
              phosy = max(0.,xa-xn)
            endif
            phosy = merge(avdic/rcar, phosy, avdic <= rcar*phosy)     ! limit phosy by available DIC
            ya = avphy+phosy
            yn = (ya+grazra*avgra*phytomi/(avphy+bkzoo))/(1.+grazra*avgra/(avphy+bkzoo))
            grazing = max(0.,ya-yn)
            graton = epsher*(1.-zinges)*grazing
            gratpoc = (1.-epsher)*grazing
            grawa = epsher*zinges*grazing
            phythresh = max(0.,(ocetra(i,j,k,iphy)-2.*phytomi))
            phymor = dyphy*phythresh
            zoothresh = max(0.,(ocetra(i,j,k,izoo)-2.*grami))
            if (lkwrbioz_off) then
              bacfra = 0.
              if (use_river2omip) then
                tdoclc_rem = 0.
                tdochc_rem = 0.
              endif
            else
              bacfra = remido*ocetra(i,j,k,idoc)
              if (use_river2omip) then
                tdoclc_rem = rem_tdoclc*ocetra(i,j,k,itdoc_lc)
                tdochc_rem = rem_tdochc*ocetra(i,j,k,itdoc_hc)
              endif
            endif
            exud = gammap*phythresh
            zoomor = spemor*zoothresh*zoothresh           ! *10 compared to linear in tropics (tinka)
            excdoc = gammaz*zoothresh                     ! excretion of doc by zooplankton
            export = zoomor*(1.-ecan) + phymor + gratpoc  ! ecan=.95, gratpoc= .2*grazing

            if (use_cisonew) then
              ! calculation of isotope fractionation during photosynthesis (Laws 1997)
              if(ocetra(i,j,k,iphy) < phytomi) then
                bifr13 = 1.
              else
                phygrowth   = ((ocetra(i,j,k,iphy)+phosy)/ocetra(i,j,k,iphy))/dtb ! Growth rate phytoplankton [1/d]
                growth_co2  = phygrowth/(co2star(i,j,k)*1.e6+safediv)             ! CO2* in [mol/kg]
                bifr13_perm = (6.03 + 5.5*growth_co2)/(0.225 + growth_co2)        ! Permil (~20)
                bifr13_perm = max(5.,min(26.,bifr13_perm))                        ! Limit the range to [5,26]
                bifr13      = (1000. - bifr13_perm) / 1000.                       ! Fractionation factor 13c (~0.98)
              endif

              bifr14 = bifr13**2

              ! calculation of 13C and 14C equivalent of biology
              rco213 = ocetra(i,j,k,isco213)/(ocetra(i,j,k,isco212)+safediv)
              rco214 = ocetra(i,j,k,isco214)/(ocetra(i,j,k,isco212)+safediv)
              rphy13 = ocetra(i,j,k,iphy13)/(ocetra(i,j,k,iphy)+safediv)
              rphy14 = ocetra(i,j,k,iphy14)/(ocetra(i,j,k,iphy)+safediv)
              rzoo13 = ocetra(i,j,k,izoo13)/(ocetra(i,j,k,izoo)+safediv)
              rzoo14 = ocetra(i,j,k,izoo14)/(ocetra(i,j,k,izoo)+safediv)

              phosy13 = phosy*bifr13*rco213
              phosy14 = phosy*bifr14*rco214

              grazing13 = grazing*rphy13
              grazing14 = grazing*rphy14

              graton13 = epsher*(1.-zinges)*grazing13
              graton14 = epsher*(1.-zinges)*grazing14

              gratpoc13 = (1.-epsher)*grazing13
              gratpoc14 = (1.-epsher)*grazing14

              grawa13 = epsher*zinges*grazing13
              grawa14 = epsher*zinges*grazing14

              if (lkwrbioz_off) then
                bacfra13 = 0.
                bacfra14 = 0.
                if (use_river2omip) then
                  tdoclc_rem13 = 0.
                  tdochc_rem13 = 0.
                  tdoclc_rem14 = 0.
                  tdochc_rem14 = 0.
                endif
              else
                bacfra13 = remido*ocetra(i,j,k,idoc13)
                bacfra14 = remido*ocetra(i,j,k,idoc14)
                if (use_river2omip) then
                  tdoclc_rem13 = rem_tdoclc*ocetra(i,j,k,itdoc_lc13)
                  tdochc_rem13 = rem_tdochc*ocetra(i,j,k,itdoc_hc13)
                  tdoclc_rem14 = rem_tdoclc*ocetra(i,j,k,itdoc_lc14)
                  tdochc_rem14 = rem_tdochc*ocetra(i,j,k,itdoc_hc14)
                endif
              endif

              phymor13 = phymor*rphy13
              phymor14 = phymor*rphy14

              zoomor13 = zoomor*rzoo13
              zoomor14 = zoomor*rzoo14

              excdoc13 = excdoc*rzoo13
              excdoc14 = excdoc*rzoo14

              exud13 = exud*rphy13
              exud14 = exud*rphy14

              export13 = zoomor13*(1.-ecan) + phymor13 + gratpoc13
              export14 = zoomor14*(1.-ecan) + phymor14 + gratpoc14
            endif

            if (use_AGG) then
              delsil = min(ropal*phosy*avsil/(avsil+bkopal),0.5*avsil)
              delcar = rcalc*min(calmax*phosy,(phosy-delsil/ropal))
              ! definition of delcar13/14 for the AGG scheme currently missing
            else
              delsil = min(ropal*export*avsil/(avsil+bkopal),0.5*avsil)
              delcar = rcalc * export * bkopal/(avsil+bkopal)
              if (use_cisonew) then
                delcar13 = rcalc * export13 * bkopal/(avsil+bkopal)
                delcar14 = rcalc * export14 * bkopal/(avsil+bkopal)
              endif
            endif

            if(with_dmsph) then
              dms_ph  = 1. + (-log10(hi(i,j,1)) - pi_ph(i,j))*dms_gamma
            else
              dms_ph  = 1.
            endif
            dmsprod = (dmsp5*delsil+dmsp4*delcar)*(1.+1./(temp+dmsp1)**2)*dms_ph
            if (lkwrbioz_off) then
               dms_bac = 0.
            else
               dms_bac = dmsp3*abs(temp+3.)*ocetra(i,j,k,idms)                 &
                       &             *(ocetra(i,j,k,idms)/(dmsp6+ocetra(i,j,k,idms)))
            endif
            dms_uv  = dmsp2*phofa/pi_alpha*ocetra(i,j,k,idms)

            dtr = bacfra-phosy+graton+ecan*zoomor

            ocetra(i,j,k,iphosph) = ocetra(i,j,k,iphosph)+dtr
            if (.not. use_extNcycle) then
              ocetra(i,j,k,iano3)   = ocetra(i,j,k,iano3)+dtr*rnit
              ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali)-2.*delcar-(rnit+1)*dtr
              ocetra(i,j,k,ioxygen) = ocetra(i,j,k,ioxygen)-dtr*ro2ut
            else
              ocetra(i,j,k,iano3)   = ocetra(i,j,k,iano3) - (1.-nh4uptfrac)*phosy*rnit
              ocetra(i,j,k,ianh4)   = ocetra(i,j,k,ianh4) - nh4uptfrac*phosy*rnit + (dtr+phosy)*rnit
              ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali) - nh4uptfrac*phosy*(rnit-1.)             & ! NH4 + PO4 Uptake
                                    &                       + (1.-nh4uptfrac)*phosy*(rnit+1.)        & ! NO3 + PO4 Uptake
                                    &                       + (dtr+phosy)*(rnit-1.)  - 2.*delcar       ! Remin to (NH4 + PO4) and CaCO3 formation
              ocetra(i,j,k,ioxygen) = ocetra(i,j,k,ioxygen) + nh4uptfrac*phosy*ro2utammo             & ! NH4 uptake
                                    &                       + (1.-nh4uptfrac)*phosy*ro2ut            & ! NO3 uptake
                                    &                       - (dtr+phosy)*ro2utammo                    ! Remin to NH4
              ! Output
              phosy_NH4(i,j,k)   =  nh4uptfrac*phosy*rnit      ! kmol N/m3/dtb - NH4 uptake during PP growth
              phosy_NO3(i,j,k)   = (1.-nh4uptfrac)*phosy*rnit  ! kmol N/m3/dtb - NO3 uptake during PP growth
              remin_aerob(i,j,k) = (dtr+phosy)*rnit            ! kmol N/m3/dtb - Aerob remin to ammonium  (var. sources)
            endif
            if (use_river2omip) then
              ocetra(i,j,k,iphosph) = ocetra(i,j,k,iphosph) + tdoclc_rem+tdochc_rem
              if (.not. use_extNcycle) then
                ocetra(i,j,k,iano3)   = ocetra(i,j,k,iano3)   + tdoclc_rem*rnit_tdoclc             &
                                      &                       + tdochc_rem*rnit_tdochc
                ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali) - (rnit_tdoclc+1.)*tdoclc_rem        &
                                      &                       - (rnit_tdochc+1.)*tdochc_rem
                ocetra(i,j,k,ioxygen) = ocetra(i,j,k,ioxygen) - tdoclc_rem*ro2ut_tdoclc            &
                                      &                       - tdochc_rem*ro2ut_tdochc
              else
                ocetra(i,j,k,ianh4)   = ocetra(i,j,k,ianh4)   + tdoclc_rem*rnit_tdoclc             &
                                      &                       + tdochc_rem*rnit_tdochc
                ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali) + tdoclc_rem*(rnit_tdoclc-1.)        &
                                      &                       + tdochc_rem*(rnit_tdochc-1.)
                ocetra(i,j,k,ioxygen) = ocetra(i,j,k,ioxygen) - tdoclc_rem*ro2utammo_tdoclc        &
                                                              - tdochc_rem*ro2utammo_tdochc
                remin_aerob(i,j,k)    = remin_aerob(i,j,k)    + tdoclc_rem*rnit_tdoclc             &
                                      &                       + tdochc_rem*rnit_tdochc
              endif
              ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212) + tdoclc_rem*rcar_tdoclc               &
                                    &                       + tdochc_rem*rcar_tdochc
              ocetra(i,j,k,iiron)   = ocetra(i,j,k,iiron)   + (tdoclc_rem+tdochc_rem)*riron
            endif
            ocetra(i,j,k,idet) = ocetra(i,j,k,idet)+export
            ocetra(i,j,k,idms) = ocetra(i,j,k,idms)+dmsprod-dms_bac-dms_uv
            ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212)-delcar+rcar*dtr
            ocetra(i,j,k,iphy) = ocetra(i,j,k,iphy)+phosy-grazing-phymor-exud
            ocetra(i,j,k,izoo) = ocetra(i,j,k,izoo)+grawa-excdoc-zoomor
            ocetra(i,j,k,idoc) = ocetra(i,j,k,idoc)-bacfra+excdoc+exud
            if (use_river2omip) then
              ocetra(i,j,k,itdoc_lc) = ocetra(i,j,k,itdoc_lc)-tdoclc_rem
              ocetra(i,j,k,itdoc_hc) = ocetra(i,j,k,itdoc_hc)-tdochc_rem
            endif
            ocetra(i,j,k,icalc) = ocetra(i,j,k,icalc)+delcar
            if (use_cisonew) then
              dtr13 = bacfra13-phosy13+graton13+ecan*zoomor13
              dtr14 = bacfra14-phosy14+graton14+ecan*zoomor14

              ocetra(i,j,k,idet13) = ocetra(i,j,k,idet13)+export13
              ocetra(i,j,k,idet14) = ocetra(i,j,k,idet14)+export14
              ocetra(i,j,k,isco213) = ocetra(i,j,k,isco213)-delcar13+rcar*dtr13
              ocetra(i,j,k,isco214) = ocetra(i,j,k,isco214)-delcar14+rcar*dtr14
              ocetra(i,j,k,iphy13) = ocetra(i,j,k,iphy13)+phosy13-grazing13-phymor13-exud13
              ocetra(i,j,k,iphy14) = ocetra(i,j,k,iphy14)+phosy14-grazing14-phymor14-exud14
              ocetra(i,j,k,izoo13) = ocetra(i,j,k,izoo13)+grawa13-excdoc13-zoomor13
              ocetra(i,j,k,izoo14) = ocetra(i,j,k,izoo14)+grawa14-excdoc14-zoomor14
              ocetra(i,j,k,idoc13) = ocetra(i,j,k,idoc13)-bacfra13+excdoc13+exud13
              ocetra(i,j,k,idoc14) = ocetra(i,j,k,idoc14)-bacfra14+excdoc14+exud14
              ocetra(i,j,k,icalc13) = ocetra(i,j,k,icalc13)+delcar13
              ocetra(i,j,k,icalc14) = ocetra(i,j,k,icalc14)+delcar14
              if (use_river2omip) then
                ocetra(i,j,k,isco213) = ocetra(i,j,k,isco213) + tdoclc_rem13*rcar_tdoclc           &
                                      &                       + tdochc_rem13*rcar_tdochc
                ocetra(i,j,k,isco214) = ocetra(i,j,k,isco214) + tdoclc_rem14*rcar_tdoclc           &
                                      &                       + tdochc_rem14*rcar_tdochc    
                ocetra(i,j,k,itdoc_lc13) = ocetra(i,j,k,itdoc_lc13)-tdoclc_rem13
                ocetra(i,j,k,itdoc_hc13) = ocetra(i,j,k,itdoc_hc13)-tdochc_rem13
                ocetra(i,j,k,itdoc_lc14) = ocetra(i,j,k,itdoc_lc14)-tdoclc_rem14
                ocetra(i,j,k,itdoc_hc14) = ocetra(i,j,k,itdoc_hc14)-tdochc_rem14
              endif
            endif
            if (use_natDIC) then
              ocetra(i,j,k,inatsco212) = ocetra(i,j,k,inatsco212)-delcar+rcar*dtr
              ocetra(i,j,k,inatalkali) = ocetra(i,j,k,inatalkali)-2.*delcar-(rnit+1)*dtr
              ocetra(i,j,k,inatcalc) = ocetra(i,j,k,inatcalc)+delcar
            endif
            if (lkwrbioz_off) then
                  opalrem = 0.
            else
               if (use_M4AGO) then
                  opalrem = dremopal*opal_remin_q10**((ptho(i,j,k)-opal_remin_Tref)/10.)*ocetra(i,j,k,iopal)
               else
                  opalrem = dremopal*ocetra(i,j,k,iopal)
               endif
            endif
            ocetra(i,j,k,isilica) = ocetra(i,j,k,isilica)-delsil+opalrem
            ocetra(i,j,k,iopal) = ocetra(i,j,k,iopal)+delsil-opalrem
            ocetra(i,j,k,iiron) = ocetra(i,j,k,iiron)+dtr*riron                     &
                 &                - relaxfe*max(ocetra(i,j,k,iiron)-fesoly,0.)

            if (use_BROMO) then
              ! Bromo source from phytoplankton production and sink to photolysis
              ! Hense and Quack (200) Pg537 Decay time scale is 30days =0.0333/day
              ! sinks owing to degradation by nitrifiers (Pg 538 of Hense and Quack,
              ! 2009) is omitted because the magnitude is more than 2 order smaller
              ! than sink through halide substitution & hydrolysis (Fig. 3)
              ! Assume that only 30% of incoming radiation are UV (i.e. 50% of non-PAR
              ! radiation; PAR radiationis assume to be 40% of incoming radiation)
              bro_beta = rbro*(fbro1*avsil/(avsil+bkopal)+fbro2*bkopal/(avsil+bkopal))
              if (swa_clim(i,j,1) > 0.) then
                bro_uv = 0.0333*dtb*0.3*(strahl(i,j)/swa_clim(i,j,1))*abs_uv(i,j,k)*ocetra(i,j,k,ibromo)
              else
                bro_uv = 0.0
              endif
              ocetra(i,j,k,ibromo) = ocetra(i,j,k,ibromo)+bro_beta*phosy-bro_uv
            endif

            if (use_AGG) then

              !***********************************************************************
              ! effects of biological processes on number of particles:
              ! photosynthesis creates POM
              ! exudation deletes POM
              ! grazing deletes POM; but only the fraction that is not egested as
              ! fecal pellets again (grawa remains in zoo, graton goes to po4)
              ! none of the processes at the current time is assumed to change
              ! the size distribution (subject to change)
              ! NOTE that phosy, exud etc. are in kmol/m3!
              ! Thus divide by avmass (kmol/m3)
              !**********************************************************************

              if(avmass > 0.) then
                avnos = ocetra(i,j,k,inos)
                anosloss = (phosy-exud-graton-grawa)*avnos/avmass
                ocetra(i,j,k,inos) = ocetra(i,j,k,inos)+anosloss
              endif

              !***********************************************************************
              ! dead zooplankton corpses come with their own, flat distribution
              ! this flow even takes place if there is neither nos nor mass
              ! NOTE: zoomor is in kmol/m3!! Thus multiply flow by 1.e+6
              !***********************************************************************

              zmornos = zoomor * (1.-ecan) * zdis * 1.e+6
              ocetra(i,j,k,inos) = ocetra(i,j,k,inos)+zmornos
            endif

            ! add up for total inventory and output
            dz = pddpo(i,j,k)

            expoor(i,j)     = expoor(i,j)    +export*rcar*dz
            expoca(i,j)     = expoca(i,j)    +delcar*dz
            exposi(i,j)     = exposi(i,j)    +delsil*dz
            intdmsprod(i,j) = intdmsprod(i,j)+dmsprod*dz
            intdms_bac(i,j) = intdms_bac(i,j)+dms_bac*dz
            intdms_uv(i,j)  = intdms_uv (i,j)+dms_uv*dz

            if (use_BROMO) then
              int_chbr3_uv(i,j)  = int_chbr3_uv (i,j) + bro_uv*dz
              int_chbr3_prod(i,j)  = int_chbr3_prod (i,j) + bro_beta*phosy*dz
            endif

            intphosy(i,j)   = intphosy(i,j)  +phosy*rcar*dz  ! primary production in kmol C m-2
            phosy3d(i,j,k)  = phosy*rcar                     ! primary production in kmol C m-3


          endif         ! pddpo(i,j,k) > dp_min
        enddo         ! kwrbioz
      enddo         ! kpie
    enddo loop1   ! kpje

    !$OMP END PARALLEL DO

    if (use_PBGC_OCNP_TIMESTEP) then
      if (mnproc == 1) then
        write(io_stdo_bgc,*)' '
        write(io_stdo_bgc,*)'in OCRPOD after 1st bio prod'
      endif
      call inventory_bgc(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
    endif

    !$OMP PARALLEL DO PRIVATE(phythresh,zoothresh,sterph,sterzo,remin     &
    !$OMP  ,opalrem,aou,refra,dms_bac,pocrem,docrem,phyrem,dz,o2lim       &
    !$OMP  ,avmass,avnos,zmornos,tdoclc_rem,tdochc_rem                    &
    !$OMP  ,tdoclc_rem13,tdochc_rem13,tdoclc_rem14,tdochc_rem14           &
    !$OMP  ,rtdoclc13,rtdochc13,rtdoclc14,rtdochc14                       &
    !$OMP  ,rphy13,rphy14,rzoo13,rzoo14,rdet13,rdet14,rdoc13,rdoc14       &
    !$OMP  ,sterph13,sterph14,sterzo13,sterzo14,pocrem13,pocrem14         &
    !$OMP  ,docrem13,docrem14,phyrem13,phyrem14                           &
    !$OMP  ,i,k)

    loop2: do j = 1,kpje
      do i = 1,kpie
        do k = merge(1,kwrbioz(i,j)+1,lkwrbioz_off),kpke
          if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then

            if (use_AGG) then
              avmass = ocetra(i,j,k,iphy)+ocetra(i,j,k,idet)
            endif
            temp = min(40.,max(-3.,ptho(i,j,k)))
            phythresh = max(0.,(ocetra(i,j,k,iphy)-2.*phytomi))
            zoothresh = max(0.,(ocetra(i,j,k,izoo)-2.*grami))
            sterph = 0.5*dyphy*phythresh                                ! phytoplankton to detritus
            sterzo = spemor*zoothresh*zoothresh                         ! quadratic mortality
            if (use_cisonew) then
              rphy13 = ocetra(i,j,k,iphy13)/(ocetra(i,j,k,iphy)+safediv)
              rphy14 = ocetra(i,j,k,iphy14)/(ocetra(i,j,k,iphy)+safediv)
              rzoo13 = ocetra(i,j,k,izoo13)/(ocetra(i,j,k,izoo)+safediv)
              rzoo14 = ocetra(i,j,k,izoo14)/(ocetra(i,j,k,izoo)+safediv)
              rdet13 = ocetra(i,j,k,idet13)/(ocetra(i,j,k,idet)+safediv)
              rdet14 = ocetra(i,j,k,idet14)/(ocetra(i,j,k,idet)+safediv)
              rdoc13 = ocetra(i,j,k,idoc13)/(ocetra(i,j,k,idoc)+safediv)
              rdoc14 = ocetra(i,j,k,idoc14)/(ocetra(i,j,k,idoc)+safediv)
              if (use_river2omip) then
                rtdoclc13 = ocetra(i,j,k,itdoc_lc13)/(ocetra(i,j,k,itdoc_lc)+safediv)
                rtdochc13 = ocetra(i,j,k,itdoc_hc13)/(ocetra(i,j,k,itdoc_hc)+safediv)
                rtdoclc14 = ocetra(i,j,k,itdoc_lc14)/(ocetra(i,j,k,itdoc_lc)+safediv)
                rtdochc14 = ocetra(i,j,k,itdoc_hc14)/(ocetra(i,j,k,itdoc_hc)+safediv)
              endif 
              sterph13 = sterph*rphy13
              sterph14 = sterph*rphy14
              sterzo13 = sterzo*rzoo13
              sterzo14 = sterzo*rzoo14
            endif

            if (lkwrbioz_off) then ! dying before in PP loop
              sterph   = 0.
              sterzo   = 0.
              sterph13 = 0.
              sterph14 = 0.
              sterzo13 = 0.
              sterzo14 = 0.
            endif

            ocetra(i,j,k,iphy) = ocetra(i,j,k,iphy)-sterph
            ocetra(i,j,k,izoo) = ocetra(i,j,k,izoo)-sterzo
            if (use_cisonew) then
              ocetra(i,j,k,iphy13) = ocetra(i,j,k,iphy13)-sterph13
              ocetra(i,j,k,iphy14) = ocetra(i,j,k,iphy14)-sterph14
              ocetra(i,j,k,izoo13) = ocetra(i,j,k,izoo13)-sterzo13
              ocetra(i,j,k,izoo14) = ocetra(i,j,k,izoo14)-sterzo14
            endif

            if(ocetra(i,j,k,ioxygen) > O2thresh_aerob) then
              if (lTO2depremin) then
                ! Both, use_M4AGO and use_extNcycle switch lTO2depremin to true!
                o2lim  = ocetra(i,j,k,ioxygen)/(ocetra(i,j,k,ioxygen) + bkox_drempoc)
                pocrem = drempoc*o2lim*POM_remin_q10**((ptho(i,j,k)-POM_remin_Tref)/10.)*ocetra(i,j,k,idet)
              else
                pocrem = drempoc*ocetra(i,j,k,idet)
              endif

              if (.not. use_extNcycle) then
                pocrem = min(pocrem,                    0.33*ocetra(i,j,k,ioxygen)/ro2ut)
                docrem = min(remido*ocetra(i,j,k,idoc), 0.33*ocetra(i,j,k,ioxygen)/ro2ut)
                phyrem = min(0.5*dyphy*phythresh,       0.33*ocetra(i,j,k,ioxygen)/ro2ut)
              else
                pocrem = min(pocrem,                    0.33*ocetra(i,j,k,ioxygen)/ro2utammo)
                docrem = min(remido*ocetra(i,j,k,idoc), 0.33*ocetra(i,j,k,ioxygen)/ro2utammo)
                phyrem = min(0.5*dyphy*phythresh,       0.33*ocetra(i,j,k,ioxygen)/ro2utammo)
              endif

              if (use_river2omip) then
                tdoclc_rem = rem_tdoclc*ocetra(i,j,k,itdoc_lc)
                tdochc_rem = rem_tdochc*ocetra(i,j,k,itdoc_hc)
                tdoclc_rem = min(rem_tdoclc*ocetra(i,j,k,idoc), 0.33*ocetra(i,j,k,ioxygen)/ro2ut_tdoclc)
                tdochc_rem = min(rem_tdochc*ocetra(i,j,k,idoc), 0.33*ocetra(i,j,k,ioxygen)/ro2ut_tdochc)
                ocetra(i,j,k,itdoc_lc) = ocetra(i,j,k,itdoc_lc) - tdoclc_rem
                ocetra(i,j,k,itdoc_hc) = ocetra(i,j,k,itdoc_hc) - tdochc_rem
              endif

              if (lkwrbioz_off) then ! dying before in PP loop
                phyrem = 0.
              endif

              if (use_cisonew) then
                pocrem13 = pocrem*rdet13
                pocrem14 = pocrem*rdet14
                docrem13 = docrem*rdoc13
                docrem14 = docrem*rdoc14
                phyrem13 = phyrem*rphy13
                phyrem14 = phyrem*rphy14
                if (use_river2omip) then
                  tdoclc_rem13 = tdoclc_rem*rtdoclc13
                  tdochc_rem13 = tdochc_rem*rtdoclc13
                  tdoclc_rem14 = tdoclc_rem*rtdoclc14
                  tdochc_rem14 = tdochc_rem*rtdoclc14
                endif
              endif
            else
              pocrem = 0.
              docrem = 0.
              phyrem = 0.
              if (use_river2omip) then
                tdoclc_rem = 0.
                tdochc_rem = 0.
              endif
              if (use_cisonew) then
                pocrem13 = 0.
                docrem13 = 0.
                phyrem13 = 0.
                pocrem14 = 0.
                docrem14 = 0.
                phyrem14 = 0.
                if (use_river2omip) then
                  tdoclc_rem13 = 0.
                  tdochc_rem13 = 0.
                  tdoclc_rem14 = 0.
                  tdochc_rem14 = 0.
                endif
              endif
            endif

            ocetra(i,j,k,idet) = ocetra(i,j,k,idet) - pocrem + sterph + sterzo
            ocetra(i,j,k,idoc) = ocetra(i,j,k,idoc) - docrem
            ocetra(i,j,k,iphy) = ocetra(i,j,k,iphy) - phyrem

            remin = pocrem + docrem + phyrem

            ocetra(i,j,k,iphosph) = ocetra(i,j,k,iphosph)+remin

            if (.not. use_extNcycle) then
              ocetra(i,j,k,iano3) = ocetra(i,j,k,iano3)+remin*rnit
              ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali)-(rnit+1)*remin
              ocetra(i,j,k,ioxygen) = ocetra(i,j,k,ioxygen)-ro2ut*remin
            else
              ocetra(i,j,k,ianh4) = ocetra(i,j,k,ianh4) + remin*rnit
              ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali) + (rnit-1.)*remin
              ocetra(i,j,k,ioxygen) = ocetra(i,j,k,ioxygen) - ro2utammo*remin
              remin_aerob(i,j,k)  = remin_aerob(i,j,k)+remin*rnit ! kmol/NH4/dtb - remin to NH4 from various sources
            endif
            
            if (use_river2omip) then
              ocetra(i,j,k,iphosph) = ocetra(i,j,k,iphosph) + tdoclc_rem + tdochc_rem
              if (.not. use_extNcycle) then
                ocetra(i,j,k,iano3)   = ocetra(i,j,k,iano3)   + tdoclc_rem*rnit_tdoclc             &
                                      &                       + tdochc_rem*rnit_tdochc
                ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali) - (rnit_tdoclc+1.)*tdoclc_rem        &
                                      &                       - (rnit_tdochc+1.)*tdochc_rem
                ocetra(i,j,k,ioxygen) = ocetra(i,j,k,ioxygen) - tdoclc_rem*ro2ut_tdoclc            &
                                      &                       - tdochc_rem*ro2ut_tdochc
              else
                ocetra(i,j,k,ianh4)   = ocetra(i,j,k,ianh4)   + tdoclc_rem*rnit_tdoclc             &
                                      &                       + tdochc_rem*rnit_tdochc
                ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali) + tdoclc_rem*(rnit_tdoclc-1.)        &
                                      &                       + tdochc_rem*(rnit_tdochc-1.)
                ocetra(i,j,k,ioxygen) = ocetra(i,j,k,ioxygen) - tdoclc_rem*ro2utammo_tdoclc        &
                                                              - tdochc_rem*ro2utammo_tdochc
                remin_aerob(i,j,k)    = remin_aerob(i,j,k)    + tdoclc_rem*rnit_tdoclc             &
                                      &                       + tdochc_rem*rnit_tdochc
              endif
              ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212) + tdoclc_rem*rcar_tdoclc               &
                                    &                       + tdochc_rem*rcar_tdochc
              ocetra(i,j,k,iiron) = ocetra(i,j,k,iiron)     + (tdoclc_rem+tdochc_rem)*riron
            endif
            
            ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212)+rcar*remin
            ocetra(i,j,k,iiron) = ocetra(i,j,k,iiron)+remin*riron           &
                 &             -relaxfe*max(ocetra(i,j,k,iiron)-fesoly,0.)
            if (use_natDIC) then
              ocetra(i,j,k,inatsco212) = ocetra(i,j,k,inatsco212)+rcar*remin
              ocetra(i,j,k,inatalkali) = ocetra(i,j,k,inatalkali)-(rnit+1)*remin
            endif
            if (use_cisonew) then
              ocetra(i,j,k,idet13) = ocetra(i,j,k,idet13)-pocrem13+sterph13+sterzo13
              ocetra(i,j,k,idet14) = ocetra(i,j,k,idet14)-pocrem14+sterph14+sterzo14
              ocetra(i,j,k,idoc13) = ocetra(i,j,k,idoc13)-docrem13
              ocetra(i,j,k,idoc14) = ocetra(i,j,k,idoc14)-docrem14
              ocetra(i,j,k,iphy13) = ocetra(i,j,k,iphy13)-phyrem13
              ocetra(i,j,k,iphy14) = ocetra(i,j,k,iphy14)-phyrem14

              ocetra(i,j,k,isco213) = ocetra(i,j,k,isco213)+rcar*(pocrem13+docrem13+phyrem13)
              ocetra(i,j,k,isco214) = ocetra(i,j,k,isco214)+rcar*(pocrem14+docrem14+phyrem14)
              if (use_river2omip) then
                ocetra(i,j,k,itdoc_lc13) = ocetra(i,j,k,itdoc_lc13)-tdoclc_rem13
                ocetra(i,j,k,itdoc_hc13) = ocetra(i,j,k,itdoc_hc13)-tdochc_rem13
                ocetra(i,j,k,itdoc_lc14) = ocetra(i,j,k,itdoc_lc14)-tdoclc_rem14
                ocetra(i,j,k,itdoc_hc14) = ocetra(i,j,k,itdoc_hc14)-tdochc_rem14
              endif
            endif
            !***********************************************************************
            ! as ragueneau (2000) notes, Si(OH)4sat is about 1000 umol, but
            ! Si(OH)4 varies only between 0-100 umol
            ! so the expression dremopal*(Si(OH)4sat-Si(OH)4) would change the
            ! rate only from 0 to 100%
            !***********************************************************************
            if (use_M4AGO) then
              opalrem = dremopal*opal_remin_q10**((ptho(i,j,k)-opal_remin_Tref)/10.)*ocetra(i,j,k,iopal)
            else
              opalrem = dremopal*0.1*(temp+3.)*ocetra(i,j,k,iopal)
            endif
            ocetra(i,j,k,iopal) = ocetra(i,j,k,iopal)-opalrem
            ocetra(i,j,k,isilica) = ocetra(i,j,k,isilica)+opalrem

            if (.not. use_extNcycle) then
              !***********************************************************************
              !           There is about 1.e4 O2 on 1 N2O molecule (Broeker&Peng)
              !           refra : Tim Rixton, private communication
              !***********************************************************************
              aou = satoxy(i,j,k)-ocetra(i,j,k,ioxygen)
              refra = 1.+3.*(0.5+sign(0.5,aou-1.97e-4))
              ocetra(i,j,k,ian2o) = ocetra(i,j,k,ian2o)+remin*1.e-4*ro2ut*refra
              ocetra(i,j,k,igasnit) = ocetra(i,j,k,igasnit)-remin*1.e-4*ro2ut*refra
              ocetra(i,j,k,ioxygen) = ocetra(i,j,k,ioxygen)-remin*1.e-4*ro2ut*refra*0.5
            endif

            dms_bac = dmsp3 * abs(temp+3.) * ocetra(i,j,k,idms)                     &
                    &    * (ocetra(i,j,k,idms) / (dmsp6+ocetra(i,j,k,idms)))
            ocetra(i,j,k,idms) = ocetra(i,j,k,idms)-dms_bac

            dz = pddpo(i,j,k)
            intdms_bac(i,j) =  intdms_bac(i,j)+dms_bac*dz

            if (use_AGG) then
              !***********************************************************************
              ! loss of snow numbers due to remineralization of poc
              ! gain of snow numbers due to zooplankton mortality
              ! NOTE that remin is in kmol/m3. Thus divide by avmass (kmol/m3)
              !***********************************************************************
              if(avmass > 0.) then
                avnos = ocetra(i,j,k,inos)
                ocetra(i,j,k,inos) = ocetra(i,j,k,inos)-remin*avnos/avmass
              endif
              !***********************************************************************
              ! dead zooplankton corpses come with their own, flat distribution
              ! this flow even takes place if there is neither nos nor mass
              ! NOTE: zoomor is in kmol/m3!! Thus multiply flow by 1.e+6
              !***********************************************************************
              zmornos = sterzo * zdis * 1.e+6
              ocetra(i,j,k,inos) = ocetra(i,j,k,inos) + zmornos
            endif/*AGG*/

          endif
        enddo
      enddo
    enddo loop2
    !$OMP END PARALLEL DO

    if (use_PBGC_OCNP_TIMESTEP) then
      if (mnproc == 1) then
        write(io_stdo_bgc,*)' '
        write(io_stdo_bgc,*)'in OCRPOD after poc remin'
      endif
      call inventory_bgc(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
    endif

    if (.not. use_extNcycle) then
      ! =====>>>> Regular CMIP6 iHAMOCC version for denitrification wo extended nitrogen cycle =====>>>>
      !$OMP PARALLEL DO PRIVATE(remin,remin2o,dz,avmass,avnos,rem13,rem14,i,k)
      loop3: do j = 1,kpje
        do i = 1,kpie
          do k = merge(1,kwrbioz(i,j)+1,lkwrbioz_off),kpke
            if(omask(i,j) > 0.5) then
              if(ocetra(i,j,k,ioxygen) < O2thresh_hypoxic .and. pddpo(i,j,k) > dp_min) then
                if (use_AGG) then
                  avmass = ocetra(i,j,k,iphy) + ocetra(i,j,k,idet)
                endif

                remin   = drempoc_anaerob*min(ocetra(i,j,k,idet),0.5   *ocetra(i,j,k,iano3)/rdnit1)
                remin2o =      dremn2o*min(ocetra(i,j,k,idet),0.003 *ocetra(i,j,k,ian2o)/rdn2o1)

                if (use_cisonew) then
                  rem13 = (remin+remin2o)*ocetra(i,j,k,idet13)/(ocetra(i,j,k,idet)+safediv)
                  rem14 = (remin+remin2o)*ocetra(i,j,k,idet14)/(ocetra(i,j,k,idet)+safediv)
                endif
                ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali)+(rdnit1-1)*remin-remin2o
                ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212)+rcar*(remin+remin2o)
                ocetra(i,j,k,idet) = ocetra(i,j,k,idet)-(remin+remin2o)
                ocetra(i,j,k,iphosph) = ocetra(i,j,k,iphosph)+(remin+remin2o)
                ocetra(i,j,k,iano3) = ocetra(i,j,k,iano3)-rdnit1*remin
                ocetra(i,j,k,igasnit) = ocetra(i,j,k,igasnit)+rdnit2*remin+rdn2o2*remin2o
                ocetra(i,j,k,iiron) = ocetra(i,j,k,iiron)+riron*(remin+remin2o)
                ocetra(i,j,k,ian2o) = ocetra(i,j,k,ian2o)-rdn2o1*remin2o
                if (use_natDIC) then
                  ocetra(i,j,k,inatalkali) = ocetra(i,j,k,inatalkali)+(rdnit1-1)*remin-remin2o
                  ocetra(i,j,k,inatsco212) = ocetra(i,j,k,inatsco212)+rcar*(remin+remin2o)
                endif
                if (use_cisonew) then
                  ocetra(i,j,k,isco213) = ocetra(i,j,k,isco213)+rcar*rem13
                  ocetra(i,j,k,isco214) = ocetra(i,j,k,isco214)+rcar*rem14
                  ocetra(i,j,k,idet13) = ocetra(i,j,k,idet13)-rem13
                  ocetra(i,j,k,idet14) = ocetra(i,j,k,idet14)-rem14
                endif

                ! nitrate loss through denitrification in kmol N m-2
                dz = pddpo(i,j,k)
                intdnit(i,j) = intdnit(i,j) + rdnit0*remin*dz

                if (use_AGG) then
                  !***********************************************************************
                  ! loss of snow numbers due to remineralization of poc
                  ! NOTE that remin is in kmol/m3. Thus divide by avmass (kmol/m3)
                  !***********************************************************************
                  if(avmass > 0.) then
                    avnos = ocetra(i,j,k,inos)
                    ocetra(i,j,k,inos) = ocetra(i,j,k,inos)-(remin+remin2o)*avnos/avmass
                  endif
                endif/*AGG*/

              endif
            endif
          enddo
        enddo
      enddo loop3
      !$OMP END PARALLEL DO

      if (use_PBGC_OCNP_TIMESTEP) then
        if (mnproc == 1) then
          write(io_stdo_bgc,*)' '
          write(io_stdo_bgc,*)'in OCRPOD after remin n2o'
        endif
        call inventory_bgc(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
      endif
      ! <<<<===== end of CMIP6 version denitrification processes without extended nitrogen cycle <<<<=====
    else
      !======>>>> extended nitrogen cycle processes (aerobic and anaerobic) that follow ammonification
      inv_message = 'in OCPROD after extNcycle nitrification'
      call nitrification(kpie,kpje,kpke,kbnd,pddpo,omask,ptho)
      call extN_inv_check(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,inv_message)

      inv_message = 'in OCPROD after extNcycle denitrification NO3 -> NO2'
      call denit_NO3_to_NO2(kpie,kpje,kpke,kbnd,pddpo,omask,ptho)
      call extN_inv_check(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,inv_message)

      inv_message = 'in OCPROD after extNcycle anammox'
      call anammox(kpie,kpje,kpke,kbnd,pddpo,omask,ptho)
      call extN_inv_check(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,inv_message)

      inv_message = 'in OCPROD after extNcycle denitrification / DNRA'
      call denit_dnra(kpie,kpje,kpke,kbnd,pddpo,omask,ptho)
      call extN_inv_check(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,inv_message)
    endif

    !sulphate reduction   ! introduced 11.5.2007 to improve poc-remineralisation in the
    !                       oxygen minimum zone in the subsurface equatorial Pacific
    !                       assumption of endless pool of SO4 (typical concentration are on the order of mmol/l)
    !      js 02072007:    for other runs than current millenium (cosmos-setup) experiments this seems
    !                      to cause trouble as phosphate concentrations are too high at the depth of the oxygen
    !                      minimum in the equatorial pacific/atlantic
    !                      does it make sense to check for oxygen and nitrate deficit?

    !$OMP PARALLEL DO PRIVATE(remin,avmass,avnos,rem13,rem14,i,k)
    loop4: do j = 1,kpje
      do i = 1,kpie
        do k = merge(1,kwrbioz(i,j)+1,lkwrbioz_off),kpke
          if(omask(i,j) > 0.5 .and. pddpo(i,j,k) > dp_min) then
            if(ocetra(i,j,k,ioxygen) < O2thresh_hypoxic .and. ocetra(i,j,k,iano3) < NO3thresh_sulf ) then

              if (use_AGG) then
                avmass = ocetra(i,j,k,iphy)+ocetra(i,j,k,idet)
              endif
              remin = dremsul*ocetra(i,j,k,idet)
              if (use_cisonew) then
                rem13 = remin*ocetra(i,j,k,idet13)/(ocetra(i,j,k,idet)+safediv)
                rem14 = remin*ocetra(i,j,k,idet14)/(ocetra(i,j,k,idet)+safediv)
              endif
              ocetra(i,j,k,idet) = ocetra(i,j,k,idet)-remin
              ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali)-(rnit+1)*remin
              ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212)+rcar*remin
              ocetra(i,j,k,iphosph) = ocetra(i,j,k,iphosph)+remin
              ocetra(i,j,k,iano3) = ocetra(i,j,k,iano3)+rnit*remin
              ocetra(i,j,k,iiron) = ocetra(i,j,k,iiron)+riron*remin
              if (use_natDIC) then
                ocetra(i,j,k,inatalkali) = ocetra(i,j,k,inatalkali)-(rnit+1)*remin
                ocetra(i,j,k,inatsco212) = ocetra(i,j,k,inatsco212)+rcar*remin
              endif
              if (use_cisonew) then
                ocetra(i,j,k,idet13) = ocetra(i,j,k,idet13)-rem13
                ocetra(i,j,k,idet14) = ocetra(i,j,k,idet14)-rem14
                ocetra(i,j,k,isco213) = ocetra(i,j,k,isco213)+rcar*rem13
                ocetra(i,j,k,isco214) = ocetra(i,j,k,isco214)+rcar*rem14
              endif
              if (use_extNcycle) then
                ! Output
                remin_sulf(i,j,k) = remin ! kmol P/m3/dtb
              endif
              if (use_AGG) then
                !***********************************************************************
                ! loss of snow numbers due to remineralization of poc
                ! NOTE that remin is in kmol/m3. Thus divide by avmass (kmol/m3)
                !***********************************************************************
                if(avmass > 0.) then
                  avnos = ocetra(i,j,k,inos)
                  ocetra(i,j,k,inos) = ocetra(i,j,k,inos)-remin*avnos/avmass
                endif
              endif

            endif
          endif
        enddo
      enddo
    enddo loop4
    !$OMP END PARALLEL DO
    !    end sulphate reduction

    if (use_PBGC_OCNP_TIMESTEP) then
      if (mnproc == 1) then
        write(io_stdo_bgc,*)' '
        write(io_stdo_bgc,*)'in OCRPOD after sulphate reduction '
      endif
      call inventory_bgc(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
    endif


    if (use_AGG) then

      !**********************AGGREGATION***************************************
      ! General:
      ! Sinking speed, size distribution and aggregation are calculated
      ! as in Kriest and Evans, 2000. I assume that opal and calcium carbonate
      ! sink at the same speed as P (mass).
      !
      ! Sinking speed and aggregation: I assume that if there is no phosphorous mass,
      ! the sinking speed is the minimum sinking speed of aggregates. I further
      ! assume that then there are no particles, and that the rate of aggregation
      ! is 0. This scheme removes no P in the absence of P, but still opal and/or
      ! calcium carbonate.
      ! This could or should be changed, because silica as well as carbonate
      ! shell will add to the aggregate mass, and should be considered.
      ! Puh. Does anyone know functional relationships between
      ! size and Si or CaCO3? Perhaps on a later version, I have to
      ! take the relationship bewteen weight and size?
      !
      ! Size distribution and resulting loss of marine snow aggregates due to
      ! aggregation (aggregate(i,j,k)) and sinking speed of mass and numbers
      ! (wmass(i,j,k) and wnumb(i,j,k) are calculated in a loop over 2-kpke.
      !
      !************************************************************************

      wmass(:,:,:)     = 0.0
      wnumb(:,:,:)     = 0.0
      aggregate(:,:,:) = 0.0
      dustagg(:,:,:)   = 0.0

      do k = 1,kpke
        do j = 1,kpje
          do i = 1,kpie

            if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then

              !***********************************************************************
              !  Have a special resetting for numbers, that fixes their conc. to one
              !  depending on mass of marine snow:
              !  Compartments have already been set to 0 in
              !  ADVECTION_BGC.h and OCTDIFF_BGC.h.
              !  Ensure that if there is no mass, there are no particles, and
              !  that the number of particles is in the right range (this is crude, but
              !  is supposed to happen only due to numerical errors such as truncation or
              !  overshoots during advection)
              ! (1) avnos<<avmass, such that eps = FractDim + 1: increase numbers
              !     such that eps = FractDim + 1 + safe (currently set to 1.e-6 in BELEG_PARM)
              ! (2) avnos>>avmass, such that  Nbar (=Mass/Nos/cellmass) <=1: decrease numbers
              !     such that Nbar=1.1 (i.e. 1.1 cells per aggregate, set in BELEG_PARM)
              !************************************************************************
              avmass = ocetra(i,j,k,iphy)+ocetra(i,j,k,idet)
              snow  = avmass*1.e+6

              if(avmass > 0.) then
                ! Set minimum particle number to nmldmin in the mixed layer. This is to prevent
                ! very small values of nos (and asscociated high sinking speed if there is mass)
                ! in high latitudes during winter
                if ( k <= kmle(i,j) ) then
                  ocetra(i,j,k,inos) = max(nmldmin,ocetra(i,j,k,inos))
                endif

                ocetra(i,j,k,inos) = max(snow*pupper,ocetra(i,j,k,inos))
                ocetra(i,j,k,inos) = min(snow*plower,ocetra(i,j,k,inos))

                avnos = ocetra(i,j,k,inos)
                eps   = ((1.+ FractDim)*snow-avnos*cellmass)/(snow-avnos*cellmass)

                ! prevent epsilon from becoming exactly one of the values which are
                ! needed for the division (guide from??js)
                if (abs(eps-3.) < 1.e-15) eps = 3.+ vsmall
                if (abs(eps-4.) < 1.e-15) eps = 4.+ vsmall
                if (abs(eps-3.-SinkExp) < 1.e-15)          eps = 3.+SinkExp+vsmall
                if (abs(eps-1.-SinkExp-FractDim) < 1.e-15) eps = 1.+SinkExp+FractDim+vsmall

                e1 = 1. - eps
                e2 = 2. - eps
                e3 = 3. - eps
                e4 = 4. - eps
                es1 = e1 + SinkExp
                es3 = e3 + SinkExp
                TopF = (alar1/alow1)**e1
                TopM = TopF * TMFac

                ! SINKING SPEED FOR THIS LAYER
                wmass(i,j,k) = cellsink * ( (FractDim+e1)/ (FractDim+es1)    &
                     &         + TopM * TSFac * SinkExp / (FractDim+es1))
                wnumb(i,j,k) = cellsink * (e1/es1 + TopF*TSFac*SinkExp/es1)

                ! AGGREGATION

                ! As a first step, assume that shear in the mixed layer is high and
                ! zero below.
                if ( k <= kmle(i,j) ) then
                  fshear = fsh
                else
                  fshear = 0.
                endif


                ! shear kernel:
                sagg1 = (TopF-1.) * (TopF*alar3-alow3) * e1 / e4                     &
                     &  + 3. * (TopF*alar1-alow1)                                    &
                     &  * (TopF*alar2-alow2) * e1 * e1 / (e2*e3)
                sagg2 = TopF*((alar3 + 3.                                            &
                     &  * (alar2*alow1*e1/e2 + alar1*alow2*e1/e3) + alow3*e1/e4)     &
                     &  - TopF*alar3*(1.+3*(       e1/e2+       e1/e3)+     e1/e4))
                sagg4 = TopF * TopF * 4. * alar3
                shear_agg = (sagg1+sagg2+sagg4) * fshear

                ! settlement kernel:
                sagg1 = (TopF * TopF * alar2 * TSFac - alow2)                        &
                     &   * SinkExp / (es3 * e3 * (es3 + e1))                         &
                     &   + alow2 * ((1. - TopF * TSFac) / (e3 * es1)                 &
                     &   - (1. - TopF) / (es3*e1))
                sagg2 = TopF * e1 * (TSFac * ( alow2 - TopF * alar2) / e3            &
                     &   - (alow2 - TopF * alar2 * TSFac) / es3)
                sett_agg =  (e1*e1*sagg1+sagg2) * fse

                effsti = Stick * (ocetra(i,j,k,iopal)*1.e+6/ropal)/                  &
                     &  ((ocetra(i,j,k,iopal) * 1.e+6 / ropal) + snow)

                aggregate(i,j,k) = (shear_agg+sett_agg) * effsti * avnos * avnos

                ! dust aggregation:
                ! shear kernel:
                dfirst = dustd3 + 3. * dustd2 * alar1 + 3. * dustd1 * alar2 + alar3
                dshagg = e1 * fsh * (dfirst * TopF / e1 - (                          &
                     &   (TopF-1.)/e1*dustd3 + 3.*(TopF*alar1-alow1)/e2*dustd2       &
                     &   + 3.*(TopF*alar2-alow2)/e3*dustd1 + (TopF*alar3-alow3)/e4))

                ! settlement kernel:
                dsett = fse * dustd2 * ((e1+SinkExp*TopF*TSFac)/es1-dustsink/cellsink)

                dustagg(i,j,k) = effsti * avnos * ocetra(i,j,k,ifdust)               &
                     &           * (dshagg+dsett)

                eps3d(i,j,k)   = eps
                asize3d(i,j,k) = snow / avnos / cellmass

              else

                wmass(i,j,k) = cellsink
                wnumb(i,j,k) = 0.
                aggregate(i,j,k) = 0.
                dustagg(i,j,k) = 0.
                ocetra(i,j,k,inos) = 0.

                eps3d(i,j,k)   = 1.
                asize3d(i,j,k) = 0.

              endif ! avmass > 0

            endif ! pddpo > dp_min .and. omask > 0.5
          enddo ! i=1,kpie
        enddo ! j=1,kpje
      enddo ! k=1,kpke

    endif ! use_AGG

    !
    ! implicit method for sinking of particles:
    ! C(k,T+dt)=C(k,T) + (w*dt/ddpo(k))*(C(k-1,T+1)-C(k,T+1))
    ! -->
    ! C(k,T+dt)=(ddpo(k)*C(k,T)+w*dt*C(k-1,T+dt))/(ddpo(k)+w*dt)
    ! sedimentation=w*dt*C(ks,T+dt)
    !
    !$OMP PARALLEL DO PRIVATE(kdonor,wpoc,wpocd,wcal,wcald,wopal,wopald,wdust,wdustd,tco,tcn,q,wnos,wnosd,dagg,i,k) ORDERED
    do j = 1,kpje
      do i = 1,kpie

        tco(:) = 0.0
        tcn(:) = 0.0

        if(omask(i,j) > 0.5) then

          kdonor = 1
          do k = 1,kpke
            !$OMP ORDERED
            ! Sum up total column inventory before sinking scheme
            if( pddpo(i,j,k) > dp_min ) then
              tco( 1) = tco( 1) + ocetra(i,j,k,idet  )*pddpo(i,j,k)
              tco( 2) = tco( 2) + ocetra(i,j,k,icalc )*pddpo(i,j,k)
              if (use_natDIC) then
                tco( 3) = tco( 3) + ocetra(i,j,k,inatcalc)*pddpo(i,j,k)
              endif
              tco( 4) = tco( 4) + ocetra(i,j,k,iopal )*pddpo(i,j,k)
              tco( 5) = tco( 5) + ocetra(i,j,k,ifdust)*pddpo(i,j,k)
              if (use_AGG) then
                tco( 6) = tco( 6) + ocetra(i,j,k,iphy  )*pddpo(i,j,k)
                tco( 7) = tco( 7) + ocetra(i,j,k,inos  )*pddpo(i,j,k)
                tco( 8) = tco( 8) + ocetra(i,j,k,iadust)*pddpo(i,j,k)
              endif
              if (use_cisonew) then
                tco( 9) = tco( 9) + ocetra(i,j,k,idet13 )*pddpo(i,j,k)
                tco(10) = tco(10) + ocetra(i,j,k,idet14 )*pddpo(i,j,k)
                tco(11) = tco(11) + ocetra(i,j,k,icalc13)*pddpo(i,j,k)
                tco(12) = tco(12) + ocetra(i,j,k,icalc14)*pddpo(i,j,k)
              endif
            endif

            if(pddpo(i,j,k) > dp_min_sink) then

              if (use_AGG) then
                wpoc   = wmass(i,j,k)
                wpocd  = wmass(i,j,kdonor)
                wcal   = wmass(i,j,k)
                wcald  = wmass(i,j,kdonor)
                wopal  = wmass(i,j,k)
                wopald = wmass(i,j,kdonor)
                wnos   = wnumb(i,j,k)
                wnosd  = wnumb(i,j,kdonor)
                wdust  = dustsink
                wdustd = dustsink
                dagg   = dustagg(i,j,k)
              else if (use_WLIN) then
                wpoc   = min(wmin+wlin*ptiestu(i,j,k),     wmax)
                wpocd  = min(wmin+wlin*ptiestu(i,j,kdonor),wmax)
                wcal   = wcal_const
                wcald  = wcal_const
                wopal  = wopal_const
                wopald = wopal_const
                wdust  = wdust_const
                wdustd = wdust_const
                dagg   = 0.0
              else if (use_M4AGO) then
                wpoc   = ws_agg(i,j,k)
                wpocd  = ws_agg(i,j,kdonor)
                wcal   = ws_agg(i,j,k)
                wcald  = ws_agg(i,j,kdonor)
                wopal  = ws_agg(i,j,k)
                wopald = ws_agg(i,j,kdonor)
                wdust  = ws_agg(i,j,k)
                wdustd = ws_agg(i,j,kdonor)
                dagg   = 0.0
              else
                wpoc   = wpoc_const
                wpocd  = wpoc_const
                wcal   = wcal_const
                wcald  = wcal_const
                wopal  = wopal_const
                wopald = wopal_const
                wdust  = wdust_const
                wdustd = wdust_const
                dagg   = 0.0
              endif

              if( k == 1 ) then
                wpocd  = 0.0
                wcald  = 0.0
                wopald = 0.0
                wdustd = 0.0
                if (use_AGG) then
                  wnosd  = 0.0
                else if (use_WLIN) then
                  wpoc = wmin
                else if (use_M4AGO) then
                  wpoc = ws_agg(i,j,k)
                endif
              endif

              ocetra(i,j,k,idet)   = (ocetra(i,j,k,     idet)  *pddpo(i,j,k)                       &
                                   +  ocetra(i,j,kdonor,idet)  *wpocd) / (pddpo(i,j,k)+wpoc)
              ocetra(i,j,k,icalc)  = (ocetra(i,j,k,     icalc) *pddpo(i,j,k)                       &
                   &               +  ocetra(i,j,kdonor,icalc) *wcald) / (pddpo(i,j,k)+wcal)
              ocetra(i,j,k,iopal)  = (ocetra(i,j,k,     iopal) *pddpo(i,j,k)                       &
                   &               +  ocetra(i,j,kdonor,iopal) *wopald)/ (pddpo(i,j,k)+wopal)
              ocetra(i,j,k,ifdust) = (ocetra(i,j,k,     ifdust)*pddpo(i,j,k)                       &
                   &               +  ocetra(i,j,kdonor,ifdust)*wdustd)/ (pddpo(i,j,k)+wdust)      &
                                   -  dagg
              if (use_cisonew) then
                ocetra(i,j,k,idet13)  = (ocetra(i,j,k,     idet13) *pddpo(i,j,k)                   &
                     &                +  ocetra(i,j,kdonor,idet13) *wpocd) / (pddpo(i,j,k)+wpoc)
                ocetra(i,j,k,idet14)  = (ocetra(i,j,k,     idet14) *pddpo(i,j,k)                   &
                     &                +  ocetra(i,j,kdonor,idet14) *wpocd) / (pddpo(i,j,k)+wpoc)
                ocetra(i,j,k,icalc13) = (ocetra(i,j,k,     icalc13)*pddpo(i,j,k)                   &
                     &                +  ocetra(i,j,kdonor,icalc13)*wcald) / (pddpo(i,j,k)+wcal)
                ocetra(i,j,k,icalc14) = (ocetra(i,j,k,     icalc14)*pddpo(i,j,k)                   &
                     &                +  ocetra(i,j,kdonor,icalc14)*wcald) / (pddpo(i,j,k)+wcal)
              endif
              if (use_natDIC) then
                ocetra(i,j,k,inatcalc)= (ocetra(i,j,k,     inatcalc)*pddpo(i,j,k)                  &
                     &                +  ocetra(i,j,kdonor,inatcalc)*wcald) / (pddpo(i,j,k)+wcal)
              endif
              if (use_AGG) then
                ocetra(i,j,k,iphy)    = (ocetra(i,j,k,     iphy)*pddpo(i,j,k)                      &
                     &                +  ocetra(i,j,kdonor,iphy)*wpocd) / (pddpo(i,j,k)+wpoc)
                ocetra(i,j,k,inos)    = (ocetra(i,j,k,     inos)*pddpo(i,j,k)                      &
                     &                +  ocetra(i,j,kdonor,inos)*wnosd) / (pddpo(i,j,k)+wnos)      &
                                      -  aggregate(i,j,k)
                ocetra(i,j,k,iadust)  = (ocetra(i,j,k,    iadust)*pddpo(i,j,k)                     &
                     &                +  ocetra(i,j,kdonor,iadust)*wpocd) / (pddpo(i,j,k)+wpoc)    &
                                      +  dagg
              endif
              kdonor = k

            else if( pddpo(i,j,k) > dp_min ) then

              ocetra(i,j,k,idet)   = ocetra(i,j,kdonor,idet)
              ocetra(i,j,k,icalc)  = ocetra(i,j,kdonor,icalc)
              if (use_cisonew) then
                ocetra(i,j,k,idet13) = ocetra(i,j,kdonor,idet13)
                ocetra(i,j,k,idet14) = ocetra(i,j,kdonor,idet14)
                ocetra(i,j,k,icalc13) = ocetra(i,j,kdonor,icalc13)
                ocetra(i,j,k,icalc14) = ocetra(i,j,kdonor,icalc14)
              endif
              if (use_natDIC) then
                ocetra(i,j,k,inatcalc) = ocetra(i,j,kdonor,inatcalc)
              endif
              ocetra(i,j,k,iopal)  = ocetra(i,j,kdonor,iopal)
              ocetra(i,j,k,ifdust) = ocetra(i,j,kdonor,ifdust)
              if (use_AGG) then
                ocetra(i,j,k,iphy)   = ocetra(i,j,kdonor,iphy)
                ocetra(i,j,k,inos)   = ocetra(i,j,kdonor,inos)
                ocetra(i,j,k,iadust) = ocetra(i,j,kdonor,iadust)
              endif

            endif  ! pddpo > dp_min_sink

            ! Sum up total column inventory after sinking scheme
            ! flux to sediment added after kpke-loop
            if( pddpo(i,j,k) > dp_min ) then
              tcn( 1) = tcn( 1) + ocetra(i,j,k,idet  )*pddpo(i,j,k)
              tcn( 2) = tcn( 2) + ocetra(i,j,k,icalc )*pddpo(i,j,k)
              if (use_natDIC) then
                tcn( 3) = tcn( 3) + ocetra(i,j,k,inatcalc)*pddpo(i,j,k)
              endif
              tcn( 4) = tcn( 4) + ocetra(i,j,k,iopal )*pddpo(i,j,k)
              tcn( 5) = tcn( 5) + ocetra(i,j,k,ifdust)*pddpo(i,j,k)
              if (use_AGG) then
                tcn( 6) = tcn( 6) + ocetra(i,j,k,iphy  )*pddpo(i,j,k)
                tcn( 7) = tcn( 7) + ocetra(i,j,k,inos  )*pddpo(i,j,k)
                tcn( 8) = tcn( 8) + ocetra(i,j,k,iadust)*pddpo(i,j,k)
              endif
              if (use_cisonew) then
                tcn( 9) = tcn( 9) + ocetra(i,j,k,idet13 )*pddpo(i,j,k)
                tcn(10) = tcn(10) + ocetra(i,j,k,idet14 )*pddpo(i,j,k)
                tcn(11) = tcn(11) + ocetra(i,j,k,icalc13)*pddpo(i,j,k)
                tcn(12) = tcn(12) + ocetra(i,j,k,icalc14)*pddpo(i,j,k)
              endif
            endif
            !$OMP END ORDERED
          enddo  ! loop k=1,kpke


          ! Add fluxes to sediment to new total column inventory
          tcn( 1) = tcn( 1) + ocetra(i,j,kdonor,idet  )*wpoc
          tcn( 2) = tcn( 2) + ocetra(i,j,kdonor,icalc )*wcal
          if (use_natDIC) then
            tcn( 3) = tcn( 3) + ocetra(i,j,kdonor,inatcalc)*wcal
          endif
          tcn( 4) = tcn( 4) + ocetra(i,j,kdonor,iopal )*wopal
          tcn( 5) = tcn( 5) + ocetra(i,j,kdonor,ifdust)*wdust
          if (use_AGG) then
            tcn( 6) = tcn( 6) + ocetra(i,j,kdonor,iphy  )*wpoc
            tcn( 7) = tcn( 7) + ocetra(i,j,kdonor,inos  )*wnos
            tcn( 8) = tcn( 8) + ocetra(i,j,kdonor,iadust)*wpoc
          endif
          if (use_cisonew) then
            tcn( 9) = tcn( 9) + ocetra(i,j,kdonor,idet13 )*wpoc
            tcn(10) = tcn(10) + ocetra(i,j,kdonor,idet14 )*wpoc
            tcn(11) = tcn(11) + ocetra(i,j,kdonor,icalc13)*wcal
            tcn(12) = tcn(12) + ocetra(i,j,kdonor,icalc14)*wcal
          endif

          ! Do columnwise multiplicative mass conservation correction
          q(:) = 1.0
          do is = 1,nsinkmax
            if( tco(is) > 1.e-12 .and. tcn(is) > 1.e-12 ) q(is) = tco(is)/tcn(is)
          enddo
          do k = 1,kpke
            if( pddpo(i,j,k) > dp_min ) then
              ocetra(i,j,k,idet  ) = ocetra(i,j,k,idet  )*q(1)
              ocetra(i,j,k,icalc ) = ocetra(i,j,k,icalc )*q(2)
              if (use_natDIC) then
                ocetra(i,j,k,inatcalc) = ocetra(i,j,k,inatcalc)*q(3)
              endif
              ocetra(i,j,k,iopal ) = ocetra(i,j,k,iopal )*q(4)
              ocetra(i,j,k,ifdust) = ocetra(i,j,k,ifdust)*q(5)
              if (use_AGG) then
                ocetra(i,j,k,iphy  ) = ocetra(i,j,k,iphy  )*q(6)
                ocetra(i,j,k,inos  ) = ocetra(i,j,k,inos  )*q(7)
                ocetra(i,j,k,iadust) = ocetra(i,j,k,iadust)*q(8)
              endif
              if (use_cisonew) then
                ocetra(i,j,k,idet13 ) = ocetra(i,j,k,idet13 )*q(9)
                ocetra(i,j,k,idet14 ) = ocetra(i,j,k,idet14 )*q(10)
                ocetra(i,j,k,icalc13) = ocetra(i,j,k,icalc13)*q(11)
                ocetra(i,j,k,icalc14) = ocetra(i,j,k,icalc14)*q(12)
              endif
            endif
          enddo

          ! Fluxes to sediment, layers thinner than dp_min_sink are ignored.
          ! Note that kdonor=kbo(i,j) by definition since kbo is the lowermost
          ! layer thicker than dp_min_sink.
          if (use_AGG) then
            prorca(i,j) = ocetra(i,j,kdonor,iphy  )*wpoc  + ocetra(i,j,kdonor,idet  )*wpoc
            prcaca(i,j) = ocetra(i,j,kdonor,icalc )*wcal
            silpro(i,j) = ocetra(i,j,kdonor,iopal )*wopal
            produs(i,j) = ocetra(i,j,kdonor,ifdust)*wdust + ocetra(i,j,kdonor,iadust)*wpoc

            if (use_cisonew) then
              pror13(i,j) = ocetra(i,j,kdonor,iphy13)*wpoc + ocetra(i,j,kdonor,idet13)*wpoc
              pror14(i,j) = ocetra(i,j,kdonor,iphy14)*wpoc + ocetra(i,j,kdonor,idet14)*wpoc
              prca13(i,j) = ocetra(i,j,kdonor,icalc13)*wcal
              prca14(i,j) = ocetra(i,j,kdonor,icalc14)*wcal
            endif
          else
            prorca(i,j) = ocetra(i,j,kdonor,idet  )*wpoc
            prcaca(i,j) = ocetra(i,j,kdonor,icalc )*wcal
            silpro(i,j) = ocetra(i,j,kdonor,iopal )*wopal
            produs(i,j) = ocetra(i,j,kdonor,ifdust)*wdust
            if (use_cisonew) then
              pror13(i,j) = ocetra(i,j,kdonor,idet13 )*wpoc
              prca13(i,j) = ocetra(i,j,kdonor,icalc13)*wcal
              pror14(i,j) = ocetra(i,j,kdonor,idet14 )*wpoc
              prca14(i,j) = ocetra(i,j,kdonor,icalc14)*wcal
            endif
          endif

        endif  ! omask > 0.5
      enddo    ! loop i=1,kpie
    enddo    ! loop j=1,kpje
    !$OMP END PARALLEL DO


    ! Calculate mass sinking flux for carbon, opal and calcium carbonate
    ! through the 100 m, 500 m, 1000 m, 2000 m, and 4000 m depth surfaces. These
    ! fluxes are intentionally calculated using values at the NEW timelevel
    ! to be fully consistent with the implicit sinking scheme

    !$OMP PARALLEL DO PRIVATE(i,k,wpoc,wcal,wopal)
    do j = 1,kpje
      do i = 1,kpie
        if(omask(i,j) > 0.5) then

          ! 100 m
          k = k0100(i,j)
          if(k > 0) then
            if (use_AGG) then
              wpoc  = wmass(i,j,k)
              wcal  = wmass(i,j,k)
              wopal = wmass(i,j,k)
              wdust = dustsink
            else if (use_WLIN) then
              wpoc  = min(wmin+wlin*ptiestu(i,j,k), wmax)
              wdust = wdust_const
            else if (use_M4AGO) then
              wpoc   = ws_agg(i,j,k)
              wcal   = ws_agg(i,j,k)
              wopal  = ws_agg(i,j,k)
              wdust  = ws_agg(i,j,k)
            else
              wpoc   = wpoc_const
              wcal   = wcal_const
              wopal  = wopal_const
              wdust  = wdust_const
            endif

            if (use_AGG) then
              carflx0100(i,j) = (ocetra(i,j,k,idet)+ocetra(i,j,k,iphy))*rcar*wpoc
              dustflx0100(i,j)= ocetra(i,j,k,ifdust)*wdust + ocetra(i,j,k,iadust)*wpoc
            else
              carflx0100(i,j) = ocetra(i,j,k,idet)*rcar*wpoc
              dustflx0100(i,j)= ocetra(i,j,k,ifdust)*wdust
            endif
            bsiflx0100(i,j) = ocetra(i,j,k,iopal)*wopal
            calflx0100(i,j) = ocetra(i,j,k,icalc)*wcal
          endif

          ! 500 m
          k = k0500(i,j)
          if(k > 0) then
            if (use_AGG) then
              wpoc  = wmass(i,j,k)
              wcal  = wmass(i,j,k)
              wopal = wmass(i,j,k)
              wdust = dustsink
            else if (use_WLIN) then
              wpoc  = min(wmin+wlin*ptiestu(i,j,k), wmax)
              wdust = wdust_const
            else if (use_M4AGO) then
              wpoc   = ws_agg(i,j,k)
              wcal   = ws_agg(i,j,k)
              wopal  = ws_agg(i,j,k)
              wdust  = ws_agg(i,j,k)
            else
              wpoc   = wpoc_const
              wcal   = wcal_const
              wopal  = wopal_const
              wdust  = wdust_const
            endif

            if (use_AGG) then
              carflx0500(i,j) = (ocetra(i,j,k,idet)+ocetra(i,j,k,iphy))*rcar*wpoc
              dustflx0500(i,j)= ocetra(i,j,k,ifdust)*wdust + ocetra(i,j,k,iadust)*wpoc
            else
              carflx0500(i,j) = ocetra(i,j,k,idet)*rcar*wpoc
              dustflx0500(i,j)= ocetra(i,j,k,ifdust)*wdust
            endif
            bsiflx0500(i,j) = ocetra(i,j,k,iopal)*wopal
            calflx0500(i,j) = ocetra(i,j,k,icalc)*wcal
          endif

          ! 1000 m
          k = k1000(i,j)
          if(k > 0) then
            if (use_AGG) then
              wpoc  = wmass(i,j,k)
              wcal  = wmass(i,j,k)
              wopal = wmass(i,j,k)
              wdust = dustsink
            else if (use_WLIN) then
              wpoc  = min(wmin+wlin*ptiestu(i,j,k), wmax)
              wdust = wdust_const
            else if (use_M4AGO) then
              wpoc   = ws_agg(i,j,k)
              wcal   = ws_agg(i,j,k)
              wopal  = ws_agg(i,j,k)
              wdust  = ws_agg(i,j,k)
            else
              wpoc   = wpoc_const
              wcal   = wcal_const
              wopal  = wopal_const
              wdust  = wdust_const
            endif

            if (use_AGG) then
              carflx1000(i,j) = (ocetra(i,j,k,idet)+ocetra(i,j,k,iphy))*rcar*wpoc
              dustflx1000(i,j)= ocetra(i,j,k,ifdust)*wdust + ocetra(i,j,k,iadust)*wpoc
            else
              carflx1000(i,j) = ocetra(i,j,k,idet)*rcar*wpoc
              dustflx1000(i,j)= ocetra(i,j,k,ifdust)*wdust
            endif
            bsiflx1000(i,j) = ocetra(i,j,k,iopal)*wopal
            calflx1000(i,j) = ocetra(i,j,k,icalc)*wcal
          endif

          ! 2000 m
          k = k2000(i,j)
          if(k > 0) then
            if (use_AGG) then
              wpoc  = wmass(i,j,k)
              wcal  = wmass(i,j,k)
              wopal = wmass(i,j,k)
              wdust = dustsink
            else if (use_WLIN) then
              wpoc  = min(wmin+wlin*ptiestu(i,j,k), wmax)
              wdust = wdust_const
            else if (use_M4AGO) then
              wpoc   = ws_agg(i,j,k)
              wcal   = ws_agg(i,j,k)
              wopal  = ws_agg(i,j,k)
              wdust  = ws_agg(i,j,k)
            else
              wpoc   = wpoc_const
              wcal   = wcal_const
              wopal  = wopal_const
              wdust  = wdust_const
            endif

            if (use_AGG) then
              carflx2000(i,j) = (ocetra(i,j,k,idet)+ocetra(i,j,k,iphy))*rcar*wpoc
              dustflx2000(i,j)= ocetra(i,j,k,ifdust)*wdust + ocetra(i,j,k,iadust)*wpoc
            else
              carflx2000(i,j) = ocetra(i,j,k,idet)*rcar*wpoc
              dustflx2000(i,j)= ocetra(i,j,k,ifdust)*wdust
            endif
            bsiflx2000(i,j) = ocetra(i,j,k,iopal)*wopal
            calflx2000(i,j) = ocetra(i,j,k,icalc)*wcal
          endif

          ! 4000 m
          k = k4000(i,j)
          if(k > 0) then
            if (use_AGG) then
              wpoc  = wmass(i,j,k)
              wcal  = wmass(i,j,k)
              wopal = wmass(i,j,k)
              wdust = dustsink
            else if (use_WLIN) then
              wpoc  = min(wmin+wlin*ptiestu(i,j,k), wmax)
              wdust = wdust_const
            else if (use_M4AGO) then
              wpoc   = ws_agg(i,j,k)
              wcal   = ws_agg(i,j,k)
              wopal  = ws_agg(i,j,k)
              wdust  = ws_agg(i,j,k)
            else
              wpoc   = wpoc_const
              wcal   = wcal_const
              wopal  = wopal_const
              wdust  = wdust_const
            endif

            if (use_AGG) then
              carflx4000(i,j) = (ocetra(i,j,k,idet)+ocetra(i,j,k,iphy))*rcar*wpoc
              dustflx4000(i,j)= ocetra(i,j,k,ifdust)*wdust + ocetra(i,j,k,iadust)*wpoc
            else
              carflx4000(i,j) = ocetra(i,j,k,idet)*rcar*wpoc
              dustflx4000(i,j)= ocetra(i,j,k,ifdust)*wdust
            endif
            bsiflx4000(i,j) = ocetra(i,j,k,iopal)*wopal
            calflx4000(i,j) = ocetra(i,j,k,icalc)*wcal
          endif

          ! bottom fluxes
          carflx_bot(i,j) = prorca(i,j)*rcar
          bsiflx_bot(i,j) = silpro(i,j)
          calflx_bot(i,j) = prcaca(i,j)
          dustflx_bot(i,j)= produs(i,j)

        endif ! omask > 0.5
      enddo
    enddo
    !$OMP END PARALLEL DO

    if (use_sedbypass) then

      ! If sediment bypass is activated, fluxes to the sediment are distributed
      ! over the water column. Detritus is kept as detritus, while opal and CaCO3
      ! are remineralised instantanously

      !$OMP PARALLEL DO PRIVATE(dz,florca,flcaca,flsil,flor13,flor14,flca13,flca14,i,k) ORDERED
      do j=1,kpje
        do i = 1,kpie
          if(omask(i,j) > 0.5) then

            ! calculate depth of water column
            dz = 0.0
            do k = 1,kpke
              !$OMP ORDERED
              if( pddpo(i,j,k) > dp_min ) dz = dz+pddpo(i,j,k)
              !$OMP END ORDERED
            enddo

            florca = prorca(i,j)/dz
            flcaca = prcaca(i,j)/dz
            flsil = silpro(i,j)/dz
            prorca(i,j) = 0.
            prcaca(i,j) = 0.
            silpro(i,j) = 0.
            if (use_cisonew) then
              flor13 = pror13(i,j)/dz
              flor14 = pror13(i,j)/dz
              flca13 = prca13(i,j)/dz
              flca14 = prca14(i,j)/dz
              pror13(i,j) = 0.
              pror14(i,j) = 0.
              prca13(i,j) = 0.
              prca14(i,j) = 0.
            endif

            do k = 1,kpke
              if( pddpo(i,j,k) <= dp_min ) cycle

              ocetra(i,j,k,idet)    = ocetra(i,j,k,idet)+florca
              ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali)+2.*flcaca
              ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212)+flcaca
              ocetra(i,j,k,isilica) = ocetra(i,j,k,isilica)+flsil
              if (use_cisonew) then
                ocetra(i,j,k,idet13)  = ocetra(i,j,k,idet13)+flor13
                ocetra(i,j,k,idet14)  = ocetra(i,j,k,idet14)+flor14
                ocetra(i,j,k,isco213) = ocetra(i,j,k,isco213)+flca13
                ocetra(i,j,k,isco214) = ocetra(i,j,k,isco214)+flca14
              endif
            enddo ! k=1,kpke

          endif ! omask > 0.5
        enddo
      enddo
      !$OMP END PARALLEL DO
    endif ! use_sedbypass

    if (use_PBGC_OCNP_TIMESTEP) then
      if (mnproc == 1) then
        write(io_stdo_bgc,*)' '
        write(io_stdo_bgc,*)'in OCRPOD after sinking poc '
      endif
      call inventory_bgc(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
    endif

  end subroutine ocprod

end module mo_ocprod
