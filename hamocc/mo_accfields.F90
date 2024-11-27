! Copyright (C) 2020  J. Schwinger, A. Moree
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

module mo_accfields

  implicit none
  private

  public :: accfields

contains

  subroutine accfields(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask)

    !***********************************************************************************************
    !
    !     J.Schwinger,    *UNI-RESEARCH*    2018-03-22
    !
    !     Modified
    !     --------
    !
    !     Purpose
    !     -------
    !     Accumulate fields for time-averaged output and write output
    !
    !***********************************************************************************************

    use mod_xc,           only: mnproc
    use mod_dia,          only: ddm
    use mo_carbch,        only: atm,atmflx,co2fxd,co2fxu,co3,hi,kwco2sol,                          &
                                ndepnoyflx,rivinflx,oalkflx,ocetra,omegaa,omegac,fco2,pco2,xco2,   &
                                pco2_gex,satoxy,sedfluxo,sedfluxb,kwco2a,co2sol,pn2om,             &
                                co213fxd,co213fxu,co214fxd,co214fxu,                               &
                                natco3,nathi,natomegaa,natomegac,natpco2,pnh3,ndepnhxflx
    use mo_biomod,        only: bsiflx_bot,bsiflx0100,bsiflx0500,bsiflx1000,                       &
                                bsiflx2000,bsiflx4000,calflx_bot,calflx0100,calflx0500,            &
                                calflx1000,calflx2000,calflx4000,carflx_bot,carflx0100,            &
                                carflx0500,carflx1000,carflx2000,carflx4000,                       &
                                dustflx_bot,dustflx0100,                                           &
                                dustflx0500,dustflx1000,dustflx2000,dustflx4000,                   &
                                expoca,expoor,exposi,intdms_bac,intdms_uv,intdmsprod,              &
                                intdnit,intnfix,intphosy,phosy3d,                                  &
                                int_chbr3_prod,int_chbr3_uv,asize3d,eps3d,wnumb,wmass,             &
                                nitr_NH4,nitr_NO2,nitr_N2O_prod,nitr_NH4_OM,nitr_NO2_OM,denit_NO3, &
                                denit_NO2,denit_N2O,DNRA_NO2,anmx_N2_prod,anmx_OM_prod,phosy_NH4,  &
                                phosy_NO3,remin_aerob,remin_sulf
    use mo_param_bgc,     only: c14fac,re1312,re14to
    use mo_bgcmean,       only: domassfluxes,jalkali,jano3,jasize,jatmco2,                         &
                                jbsiflx0100,jbsiflx0500,jbsiflx1000,jbsiflx2000,                   &
                                jbsiflx4000,jbsiflx_bot,jcalc,jcalflx0100,jcalflx0500,             &
                                jcalflx1000,jcalflx2000,jcalflx4000,                               &
                                jcalflx_bot,jcarflx0100,jcarflx0500,                               &
                                jcarflx1000,jcarflx2000,jcarflx4000,jcarflx_bot,                   &
                                jdustflx0100,jdustflx0500,                                         &
                                jdustflx1000,jdustflx2000,jdustflx4000,jdustflx_bot,               &
                                jsediffic,jsediffal,jsediffph,jsediffox,                           &
                                jburflxsso12,jburflxsssc12,jburflxssssil,jburflxssster,            &
                                jsediffn2,jsediffno3,jsediffsi,jco2flux,                           &
                                jco2fxd,jco2fxu,jco3,jdic,jdicsat,                                 &
                                jdms,jdms_bac,jdms_uv,jdmsflux,                                    &
                                jdmsprod,jdoc,jdp,jeps,jexpoca,                                    &
                                jexport,jexposi,jgrazer,jintdnit,jintnfix,jintphosy,               &
                                jiralk,jirdet,jirdin,jirdip,jirdoc,jiriron,                        &
                                jiron,jirsi,jkwco2,jlvlalkali,jlvlano3,jlvlasize,                  &
                                jlvlbigd14c,jlvlbromo,jlvlcalc,jlvlcalc13,                         &
                                jlvlcfc11,jlvlcfc12,jlvlco3,jlvld13c,jlvld14c,                     &
                                jlvldic,jlvldic13,jlvldic14,jlvldicsat,jlvldoc,                    &
                                jlvldoc13,jlvleps,jlvlgrazer,jlvlgrazer13,jlvliron,                &
                                jlvln2o,jlvlnatalkali,jlvlnatcalc,jlvlnatco3,                      &
                                jlvlnatdic,jlvlnatomegaa,jlvlnatomegac,jlvlnos,                    &
                                jlvlo2sat,jlvlomegaa,jlvlomegac,jlvlopal,jlvloxygen,               &
                                jlvlph,jlvlphosph,jlvlphosy,jlvlphyto,jlvlphyto13,                 &
                                jlvlpoc,jlvlpoc13,jlvlprefalk,jlvlprefdic,                         &
                                jlvlprefo2,jlvlprefpo4,jlvlsf6,jlvlsilica,                         &
                                jlvlwnos,jlvlwphy,jn2flux,jn2o,jn2oflux,jn2ofx,                    &
                                jprorca,jprcaca,jsilpro,jpodiic,jpodial,jpodiph,                   &
                                jpodiox,jpodin2,jpodino3,jpodisi,jndepnoy,jndepnhx,joalk,          &
                                jniflux,jnos,jo2flux,jo2sat,jomegaa,jomegac,jopal,                 &
                                joxflux,joxygen,jfco2,jpco2,jxco2,jpco2_gex,jkwco2sol,jco2sol,     &
                                jph,jphosph,jphosy,jphyto,jpoc,jprefalk,                           &
                                jprefdic,jprefo2,jprefpo4,jsilica,jsrfalkali,                      &
                                jsrfano3,jsrfdic,jsrfiron,jsrfoxygen,jsrfphosph,                   &
                                jsrfphyto,jsrfsilica,jsrfph,jwnos,jwphy,jndepnoyfx,                &
                                joalkfx,nbgc,nacc_bgc,bgcwrt,glb_inventory,                        &
                                bgct2d,acclvl,acclyr,accsrf,bgczlv,                                &
                                jatmbromo,jbromo,jbromo_prod,jbromo_uv,jbromofx,jsrfbromo,         &
                                jcfc11,jcfc11fx,jcfc12,jcfc12fx,jsf6,jsf6fx,                       &
                                jatmc13,jatmc14,jbigd14c,jcalc13,jco213fxd,jco213fxu,              &
                                jco214fxd,jco214fxu,jd13c,jd14c,jdic13,jdic14,                     &
                                jdoc13,jgrazer13,jphyto13,jpoc13,                                  &
                                jlvlnatph,jnatalkali,jnatcalc,jnatco2fx,jnatco3,                   &
                                jnatdic,jnatomegaa,jnatomegac,jnatpco2,jnatph,                     &
                                jsrfnatalk,jsrfnatdic,jsrfnatph,                                   &
                                jbursssc12,jburssso12,jburssssil,jburssster,                       &
                                jpowaal,jpowaic,jpowaox,jpowaph,jpowaph,jpowasi,jpown2,            &
                                jpowno3,jsssc12,jssso12,jssssil,jssster,accbur,accsdm,             &
                                jatmco2,jatmn2,jatmo2,                                             &
                                jlvlanh4,jlvlano2,                                                 &
                                jlvl_nitr_NH4, jsrfpn2om,                                          &
                                jlvl_nitr_NO2,jlvl_nitr_N2O_prod,jlvl_nitr_NH4_OM,jlvl_nitr_NO2_OM,&
                                jlvl_denit_NO3,jlvl_denit_NO2,jlvl_denit_N2O,jlvl_DNRA_NO2,        &
                                jlvl_anmx_N2_prod,jlvl_anmx_OM_prod,jlvl_phosy_NH4,jlvl_phosy_NO3, &
                                jlvl_remin_aerob,jlvl_remin_sulf,                                  &
                                jagg_ws,jdynvis,jagg_stick,jagg_stickf,jagg_dmax,jagg_avdp,        &
                                jagg_avrhop,jagg_avdC,jagg_df,jagg_b,jagg_Vrhof,jagg_Vpor,         &
                                jlvl_agg_ws,jlvl_dynvis,jlvl_agg_stick,jlvl_agg_stickf,            &
                                jlvl_agg_dmax,jlvl_agg_avdp,jlvl_agg_avrhop,jlvl_agg_avdC,         &
                                jlvl_agg_df,jlvl_agg_b,jlvl_agg_Vrhof,jlvl_agg_Vpor,               &
                                jprefsilica,jlvlprefsilica,                                        &
                                jnh3flux,janh3fx,janh4,jano2,jsrfanh4,jsrfano2,jsrfpnh3,           &
                                jnitr_NH4,jnitr_NO2,jnitr_N2O_prod,jnitr_NH4_OM,jnitr_NO2_OM,      &
                                jdenit_NO3,jdenit_NO2,jdenit_N2O,jDNRA_NO2,janmx_N2_prod,          &
                                janmx_OM_prod,jphosy_NH4,jphosy_NO3,jremin_aerob,jremin_sulf,      &
                                jpownh4,jpown2o,jpowno2,jsdm_nitr_NH4,jsdm_nitr_NO2,               &
                                jsdm_nitr_N2O_prod,jsdm_nitr_NH4_OM,jsdm_nitr_NO2_OM,              &
                                jsdm_denit_NO3,jsdm_denit_NO2,jsdm_denit_N2O,jsdm_DNRA_NO2,        &
                                jsdm_anmx_N2_prod,jsdm_anmx_OM_prod,jsdm_remin_aerob,              &
                                jsdm_remin_sulf,jsediffnh4,jsediffn2o,jsediffno2,jatmn2o,jatmnh3,  &
                                jndepnhxfx,jshelfage,jlvlshelfage
    use mo_control_bgc,   only: io_stdo_bgc,dtb,use_BROMO,use_AGG,use_WLIN,use_natDIC,             &
                                use_CFC,use_sedbypass,use_cisonew,use_BOXATM,use_M4AGO,            &
                                use_extNcycle,use_pref_tracers,use_shelfsea_res_time
    use mo_param1_bgc,    only: ialkali,ian2o,iano3,iatmco2,iatmdms,iatmn2,iatmn2o,iatmo2,         &
                                icalc,idet,idms,idicsat,idoc,iiron,iopal,                          &
                                ioxygen,iphosph,iphy,iprefalk,iprefdic,                            &
                                iprefpo4,iprefo2,isco212,isilica,izoo,                             &
                                irdin,irdip,irsi,iralk,iriron,irdoc,irdet,inos,iatmbromo,ibromo,   &
                                iatmf11,iatmf12,iatmsf6,icfc11,icfc12,isf6,                        &
                                iatmc13,iatmc14,icalc13,idet13,idoc13,iphy13,isco213,isco214,      &
                                izoo13,safediv,                                                    &
                                iatmnco2,inatalkali,inatcalc,inatsco212,                           &
                                ipowaal,ipowaic,ipowaox,ipowaph,ipowasi,                           &
                                ipown2,ipowno3,isssc12,issso12,issssil,issster,                    &
                                issso12,isssc12,issssil,issster,iprefsilica,iatmnh3,ianh4,iano2,   &
                                ipownh4,ipown2o,ipowno2,ishelfage
    use mo_sedmnt,        only: powtra,sedlay,burial
    use mo_vgrid,         only: dp_min
    use mo_inventory_bgc, only: inventory_bgc
    use mo_ncwrt_bgc    , only: ncwrt_bgc
    use mo_ihamocc4m4ago, only: aggregate_diagnostics,kav_dp,kav_rho_p,kav_d_C,kws_agg,kdf_agg,    &
                                kstickiness_agg,kb_agg,kstickiness_frustule,kLmax_agg,kdynvis,     &
                                kav_rhof_V,kav_por_V
    use mo_extNsediment,  only: extNsed_diagnostics,ised_nitr_NH4,ised_nitr_NO2,ised_nitr_N2O_prod,&
                                ised_nitr_NH4_OM,ised_nitr_NO2_OM,ised_denit_NO3,ised_denit_NO2,   &
                                ised_denit_N2O,ised_DNRA_NO2,ised_anmx_N2_prod,ised_anmx_OM_prod,  &
                                ised_remin_aerob,ised_remin_sulf

    ! Arguments
    integer , intent(in) :: kpie                  ! 1st dimension of model grid.
    integer , intent(in) :: kpje                  ! 2nd dimension of model grid.
    integer , intent(in) :: kpke                  ! 3rd (vertical) dimension of model grid.
    real    , intent(in) :: pdlxp(kpie,kpje)      ! size of grid cell (1st dimension) [m].
    real    , intent(in) :: pdlyp(kpie,kpje)      ! size of grid cell (2nd dimension) [m].
    real    , intent(in) :: pddpo(kpie,kpje,kpke) ! size of grid cell (3rd dimension) [m].
    real    , intent(in) :: omask(kpie,kpje)      ! land/ocean mask

    ! Local variables
    integer :: i,j,k,l
    integer :: ind1(kpie,kpje),ind2(kpie,kpje)
    real    :: wghts(kpie,kpje,ddm)
    real    :: di12C                   ! cisonew
    real    :: d13C(kpie,kpje,kpke)    ! cisonew
    real    :: d14C(kpie,kpje,kpke)    ! cisonew
    real    :: bigd14C(kpie,kpje,kpke) ! cisonew

    if (use_cisonew) then
      ! Calculation d13C, d14C and Dd14C: Delta notation for output
      d13C(:,:,:)=0.
      d14C(:,:,:)=0.
      bigd14C(:,:,:)=0.
      do k=1,kpke
        do j=1,kpje
          do i=1,kpie
            if(omask(i,j).gt.0.5.and.pddpo(i,j,k).gt.dp_min) then

              di12C=max(ocetra(i,j,k,isco212)-ocetra(i,j,k,isco213),0.)
              d13C(i,j,k)=(ocetra(i,j,k,isco213)/(di12C+safediv)/re1312-1.)*1000.
              d14C(i,j,k)=(ocetra(i,j,k,isco214)*c14fac/                                           &
                          (ocetra(i,j,k,isco212)+safediv)/re14to-1.)*1000.
              bigd14C(i,j,k)=d14C(i,j,k)-2.*(d13C(i,j,k)+25.)*(1.+d14C(i,j,k)/1000.)

            endif
          enddo
        enddo
      enddo
    endif


    ! Accumulated fluxes for inventory.F90. Note that these are currently not written to restart!
    ! Division by 2 is to account for leap-frog timestepping (but this is not exact)
    do j=1,kpje
      do i=1,kpie
        if(omask(i,j).gt.0.5) then

          ! Atmosphere-ocean fluxes
          bgct2d(i,j,jco2flux) = bgct2d(i,j,jco2flux) + atmflx(i,j,iatmco2)/2.0
          bgct2d(i,j,jo2flux)  = bgct2d(i,j,jo2flux)  + atmflx(i,j,iatmo2)/2.0
          bgct2d(i,j,jn2flux)  = bgct2d(i,j,jn2flux)  + atmflx(i,j,iatmn2)/2.0
          bgct2d(i,j,jn2oflux) = bgct2d(i,j,jn2oflux) + atmflx(i,j,iatmn2o)/2.0
          if (use_extNcycle) then
            bgct2d(i,j,jnh3flux) = bgct2d(i,j,jnh3flux) + atmflx(i,j,iatmnh3)/2.0
            bgct2d(i,j,jndepnhx) = bgct2d(i,j,jndepnhx) + ndepnhxflx(i,j)/2.0
          endif
          ! Particle fluxes between water-column and sediment
          bgct2d(i,j,jprorca)  = bgct2d(i,j,jprorca)  + carflx_bot(i,j)/2.0
          bgct2d(i,j,jprcaca)  = bgct2d(i,j,jprcaca)  + calflx_bot(i,j)/2.0
          bgct2d(i,j,jsilpro)  = bgct2d(i,j,jsilpro)  + bsiflx_bot(i,j)/2.0
          if (.not. use_sedbypass) then
            ! Diffusive fluxes between water-column and sediment
            bgct2d(i,j,jpodiic)  = bgct2d(i,j,jpodiic)  + sedfluxo(i,j,ipowaic)/2.0
            bgct2d(i,j,jpodial)  = bgct2d(i,j,jpodial)  + sedfluxo(i,j,ipowaal)/2.0
            bgct2d(i,j,jpodiph)  = bgct2d(i,j,jpodiph)  + sedfluxo(i,j,ipowaph)/2.0
            bgct2d(i,j,jpodiox)  = bgct2d(i,j,jpodiox)  + sedfluxo(i,j,ipowaox)/2.0
            bgct2d(i,j,jpodin2)  = bgct2d(i,j,jpodin2)  + sedfluxo(i,j,ipown2)/2.0
            bgct2d(i,j,jpodino3) = bgct2d(i,j,jpodino3) + sedfluxo(i,j,ipowno3)/2.0
            bgct2d(i,j,jpodisi)  = bgct2d(i,j,jpodisi)  + sedfluxo(i,j,ipowasi)/2.0
          endif
          ! N-deposition, ocean alkalinization, and riverine input fluxes
          bgct2d(i,j,jndepnoy) = bgct2d(i,j,jndepnoy) + ndepnoyflx(i,j)/2.0
          bgct2d(i,j,joalk)    = bgct2d(i,j,joalk)    + oalkflx(i,j)/2.0
          bgct2d(i,j,jirdin)   = bgct2d(i,j,jirdin)   + rivinflx(i,j,irdin)/2.0
          bgct2d(i,j,jirdip)   = bgct2d(i,j,jirdip)   + rivinflx(i,j,irdip)/2.0
          bgct2d(i,j,jirsi)    = bgct2d(i,j,jirsi)    + rivinflx(i,j,irsi)/2.0
          bgct2d(i,j,jiralk)   = bgct2d(i,j,jiralk)   + rivinflx(i,j,iralk)/2.0
          bgct2d(i,j,jiriron)  = bgct2d(i,j,jiriron)  + rivinflx(i,j,iriron)/2.0
          bgct2d(i,j,jirdoc)   = bgct2d(i,j,jirdoc)   + rivinflx(i,j,irdoc)/2.0
          bgct2d(i,j,jirdet)   = bgct2d(i,j,jirdet)   + rivinflx(i,j,irdet)/2.0

        endif
      enddo
    enddo

    ! Accumulate atmosphere fields and fluxes
    call accsrf(jatmco2,atm(1,1,iatmco2),omask,0)
    if (use_BOXATM) then
      call accsrf(jatmo2 ,atm(1,1,iatmo2),omask,0)
      call accsrf(jatmn2 ,atm(1,1,iatmn2),omask,0)
    endif
    call accsrf(joxflux,atmflx(1,1,iatmo2),omask,0)
    call accsrf(jniflux,atmflx(1,1,iatmn2),omask,0)
    call accsrf(jn2ofx,atmflx(1,1,iatmn2o),omask,0)
    call accsrf(jdmsflux,atmflx(1,1,iatmdms),omask,0)
    if (use_CFC) then
      call accsrf(jcfc11fx,atmflx(1,1,iatmf11),omask,0)
      call accsrf(jcfc12fx,atmflx(1,1,iatmf12),omask,0)
      call accsrf(jsf6fx,atmflx(1,1,iatmsf6),omask,0)
    endif
    if (use_natDIC) then
      call accsrf(jnatco2fx,atmflx(1,1,iatmnco2),omask,0)
    endif
    if (use_BROMO) then
      call accsrf(jatmbromo,atm(1,1,iatmbromo),omask,0)
      call accsrf(jbromofx,atmflx(1,1,iatmbromo),omask,0)
    endif
    if (use_cisonew) then
      call accsrf(jatmc13,atm(1,1,iatmc13),omask,0)
      call accsrf(jatmc14,atm(1,1,iatmc14),omask,0)
    endif
    if (use_extNcycle) then
      call accsrf(janh3fx,atmflx(1,1,iatmnh3),omask,0)
      call accsrf(jatmnh3,atm(1,1,iatmnh3),omask,0)
      call accsrf(jatmn2o,atm(1,1,iatmn2o),omask,0)
    endif
    ! Save up and downward fluxes for CO2 seperately
    call accsrf(jco2fxd,co2fxd,omask,0)
    call accsrf(jco2fxu,co2fxu,omask,0)
    if (use_cisonew) then
      call accsrf(jco213fxd,co213fxd,omask,0)
      call accsrf(jco213fxu,co213fxu,omask,0)
      call accsrf(jco214fxd,co214fxd,omask,0)
      call accsrf(jco214fxu,co214fxu,omask,0)
    endif

    ! Accumulate 2d diagnostics
    call accsrf(jfco2,fco2,omask,0)
    call accsrf(jpco2,pco2,omask,0)
    call accsrf(jxco2,xco2,omask,0)
    call accsrf(jpco2_gex,pco2_gex,omask,0)
    call accsrf(jkwco2sol,kwco2sol,omask,0)
    call accsrf(jkwco2,kwco2a,omask,0)
    call accsrf(jco2sol,co2sol,omask,0)
    call accsrf(jsrfphosph,ocetra(1,1,1,iphosph),omask,0)
    call accsrf(jsrfoxygen,ocetra(1,1,1,ioxygen),omask,0)
    call accsrf(jsrfiron,ocetra(1,1,1,iiron),omask,0)
    call accsrf(jsrfano3,ocetra(1,1,1,iano3),omask,0)
    call accsrf(jsrfalkali,ocetra(1,1,1,ialkali),omask,0)
    call accsrf(jsrfsilica,ocetra(1,1,1,isilica),omask,0)
    call accsrf(jsrfdic,ocetra(1,1,1,isco212),omask,0)
    call accsrf(jsrfphyto,ocetra(1,1,1,iphy),omask,0)
    call accsrf(jsrfph,hi(1,1,1),omask,0)
    call accsrf(jdms,ocetra(1,1,1,idms),omask,0)
    call accsrf(jsrfpn2om,pn2om,omask,0)
    call accsrf(jexport,expoor,omask,0)
    call accsrf(jexpoca,expoca,omask,0)
    call accsrf(jexposi,exposi,omask,0)
    call accsrf(jdmsprod,intdmsprod,omask,0)
    call accsrf(jdms_uv,intdms_uv,omask,0)
    call accsrf(jdms_bac,intdms_bac,omask,0)
    call accsrf(jintphosy,intphosy,omask,0)
    call accsrf(jintdnit,intdnit,omask,0)
    call accsrf(jintnfix,intnfix,omask,0)
    if (use_natDIC) then
      call accsrf(jsrfnatdic,ocetra(1,1,1,inatsco212),omask,0)
      call accsrf(jsrfnatalk,ocetra(1,1,1,inatalkali),omask,0)
      call accsrf(jnatpco2,natpco2,omask,0)
      call accsrf(jsrfnatph,nathi(1,1,1),omask,0)
    endif
    if (use_BROMO) then
      call accsrf(jsrfbromo,ocetra(1,1,1,ibromo),omask,0)
      call accsrf(jbromo_prod,int_chbr3_prod,omask,0)
      call accsrf(jbromo_uv,int_chbr3_uv,omask,0)
    endif

    ! Accumulate fluxes due to N-deposition, ocean alkalinization
    call accsrf(jndepnoyfx,ndepnoyflx,omask,0)
    call accsrf(joalkfx,oalkflx,omask,0)

    if (use_extNcycle) then
      call accsrf(jsrfanh4,ocetra(1,1,1,ianh4),omask,0)
      call accsrf(jsrfpnh3,pnh3,omask,0)
      call accsrf(jsrfano2,ocetra(1,1,1,iano2),omask,0)
      call accsrf(jndepnhxfx,ndepnhxflx,omask,0)
    endif

    ! Accumulate the diagnostic mass sinking field
    if( domassfluxes ) then
      call accsrf(jcarflx0100,carflx0100,omask,0)
      call accsrf(jbsiflx0100,bsiflx0100,omask,0)
      call accsrf(jcalflx0100,calflx0100,omask,0)
      call accsrf(jdustflx0100,dustflx0100,omask,0)
      call accsrf(jcarflx0500,carflx0500,omask,0)
      call accsrf(jbsiflx0500,bsiflx0500,omask,0)
      call accsrf(jcalflx0500,calflx0500,omask,0)
      call accsrf(jdustflx0500,dustflx0500,omask,0)
      call accsrf(jcarflx1000,carflx1000,omask,0)
      call accsrf(jbsiflx1000,bsiflx1000,omask,0)
      call accsrf(jcalflx1000,calflx1000,omask,0)
      call accsrf(jdustflx1000,dustflx1000,omask,0)
      call accsrf(jcarflx2000,carflx2000,omask,0)
      call accsrf(jbsiflx2000,bsiflx2000,omask,0)
      call accsrf(jcalflx2000,calflx2000,omask,0)
      call accsrf(jdustflx2000,dustflx2000,omask,0)
      call accsrf(jcarflx4000,carflx4000,omask,0)
      call accsrf(jbsiflx4000,bsiflx4000,omask,0)
      call accsrf(jcalflx4000,calflx4000,omask,0)
      call accsrf(jdustflx4000,dustflx4000,omask,0)
      call accsrf(jcarflx_bot,carflx_bot,omask,0)
      call accsrf(jbsiflx_bot,bsiflx_bot,omask,0)
      call accsrf(jcalflx_bot,calflx_bot,omask,0)
      call accsrf(jdustflx_bot,dustflx_bot,omask,0)
    endif

    if (.not. use_sedbypass) then
      ! Accumulate diffusive fluxes between water column and sediment
      call accsrf(jsediffic,sedfluxo(1,1,ipowaic),omask,0)
      call accsrf(jsediffal,sedfluxo(1,1,ipowaal),omask,0)
      call accsrf(jsediffph,sedfluxo(1,1,ipowaph),omask,0)
      call accsrf(jsediffox,sedfluxo(1,1,ipowaox),omask,0)
      call accsrf(jsediffn2,sedfluxo(1,1,ipown2),omask,0)
      call accsrf(jsediffno3,sedfluxo(1,1,ipowno3),omask,0)
      call accsrf(jsediffsi,sedfluxo(1,1,ipowasi),omask,0)
      call accsrf(jburflxsso12,sedfluxb(1,1,issso12),omask,0)
      call accsrf(jburflxsssc12,sedfluxb(1,1,isssc12),omask,0)
      call accsrf(jburflxssssil,sedfluxb(1,1,issssil),omask,0)
      call accsrf(jburflxssster,sedfluxb(1,1,issster),omask,0)
      if (use_extNcycle) then
        call accsrf(jsediffnh4,sedfluxo(1,1,ipownh4),omask,0)
        call accsrf(jsediffn2o,sedfluxo(1,1,ipown2o),omask,0)
        call accsrf(jsediffno2,sedfluxo(1,1,ipowno2),omask,0)
      endif
    endif
    ! Accumulate layer diagnostics
    call acclyr(jdp,pddpo,pddpo,0)
    call acclyr(jphyto,ocetra(1,1,1,iphy),pddpo,1)
    call acclyr(jgrazer,ocetra(1,1,1,izoo),pddpo,1)
    call acclyr(jphosph,ocetra(1,1,1,iphosph),pddpo,1)
    call acclyr(joxygen,ocetra(1,1,1,ioxygen),pddpo,1)
    call acclyr(jiron,ocetra(1,1,1,iiron),pddpo,1)
    call acclyr(jano3,ocetra(1,1,1,iano3),pddpo,1)
    call acclyr(jalkali,ocetra(1,1,1,ialkali),pddpo,1)
    call acclyr(jsilica,ocetra(1,1,1,isilica),pddpo,1)
    call acclyr(jdic,ocetra(1,1,1,isco212),pddpo,1)
    call acclyr(jdoc,ocetra(1,1,1,idoc),pddpo,1)
    call acclyr(jpoc,ocetra(1,1,1,idet),pddpo,1)
    call acclyr(jcalc,ocetra(1,1,1,icalc),pddpo,1)
    call acclyr(jopal,ocetra(1,1,1,iopal),pddpo,1)
    call acclyr(jn2o,ocetra(1,1,1,ian2o),pddpo,1)
    call acclyr(jco3,co3,pddpo,1)
    call acclyr(jph,hi,pddpo,1)
    call acclyr(jomegaa,OmegaA,pddpo,1)
    call acclyr(jomegac,OmegaC,pddpo,1)
    call acclyr(jphosy,phosy3d,pddpo,1)
    call acclyr(jo2sat,satoxy,pddpo,1)
    call acclyr(jdicsat,ocetra(1,1,1,idicsat),pddpo,1)
    if (use_pref_tracers) then
      call acclyr(jprefo2,ocetra(1,1,1,iprefo2),pddpo,1)
      call acclyr(jprefpo4,ocetra(1,1,1,iprefpo4),pddpo,1)
      call acclyr(jprefsilica,ocetra(1,1,1,iprefsilica),pddpo,1)
      call acclyr(jprefalk,ocetra(1,1,1,iprefalk),pddpo,1)
      call acclyr(jprefdic,ocetra(1,1,1,iprefdic),pddpo,1)
    endif
    if (use_shelfsea_res_time) then
      call acclyr(jshelfage,ocetra(1,1,1,ishelfage),pddpo,1)
    endif
    if (use_natDIC) then
      call acclyr(jnatalkali,ocetra(1,1,1,inatalkali),pddpo,1)
      call acclyr(jnatdic,ocetra(1,1,1,inatsco212),pddpo,1)
      call acclyr(jnatcalc,ocetra(1,1,1,inatcalc),pddpo,1)
      call acclyr(jnatco3,natco3,pddpo,1)
      call acclyr(jnatph,nathi,pddpo,1)
      call acclyr(jnatomegaa,natOmegaA,pddpo,1)
      call acclyr(jnatomegac,natOmegaC,pddpo,1)
    endif
    if (use_cisonew) then
      call acclyr(jdic13,ocetra(1,1,1,isco213),pddpo,1)
      call acclyr(jdic14,ocetra(1,1,1,isco214),pddpo,1)
      call acclyr(jd13c,d13C,pddpo,1)
      call acclyr(jd14c,d14C,pddpo,1)
      call acclyr(jbigd14c,bigd14C,pddpo,1)
      call acclyr(jpoc13,ocetra(1,1,1,idet13),pddpo,1)
      call acclyr(jdoc13,ocetra(1,1,1,idoc13),pddpo,1)
      call acclyr(jcalc13,ocetra(1,1,1,icalc13),pddpo,1)
      call acclyr(jphyto13,ocetra(1,1,1,iphy13),pddpo,1)
      call acclyr(jgrazer13,ocetra(1,1,1,izoo13),pddpo,1)
    endif
    if (use_AGG) then
      call acclyr(jnos,ocetra(1,1,1,inos),pddpo,1)
      call acclyr(jwphy, wmass/dtb,pddpo,1)
      call acclyr(jwnos, wnumb/dtb,pddpo,1)
      call acclyr(jeps,  eps3d,    pddpo,1)
      call acclyr(jasize,asize3d,  pddpo,1)
    endif
    if (use_CFC) then
      call acclyr(jcfc11,ocetra(1,1,1,icfc11),pddpo,1)
      call acclyr(jcfc12,ocetra(1,1,1,icfc12),pddpo,1)
      call acclyr(jsf6,ocetra(1,1,1,isf6),pddpo,1)
    endif
    if (use_BROMO) then
      call acclyr(jbromo,ocetra(1,1,1,ibromo),pddpo,1)
    endif
    if (use_extNcycle) then
      call acclyr(janh4,ocetra(1,1,1,ianh4),pddpo,1)
      call acclyr(jano2,ocetra(1,1,1,iano2),pddpo,1)
      call acclyr(jnitr_NH4,nitr_NH4,pddpo,1)
      call acclyr(jnitr_NO2,nitr_NO2,pddpo,1)
      call acclyr(jnitr_N2O_prod,nitr_N2O_prod,pddpo,1)
      call acclyr(jnitr_NH4_OM,nitr_NH4_OM,pddpo,1)
      call acclyr(jnitr_NO2_OM,nitr_NO2_OM,pddpo,1)
      call acclyr(jdenit_NO3,denit_NO3,pddpo,1)
      call acclyr(jdenit_NO2,denit_NO2,pddpo,1)
      call acclyr(jdenit_N2O,denit_N2O,pddpo,1)
      call acclyr(jDNRA_NO2,DNRA_NO2,pddpo,1)
      call acclyr(janmx_N2_prod,anmx_N2_prod,pddpo,1)
      call acclyr(janmx_OM_prod,anmx_OM_prod,pddpo,1)
      call acclyr(jphosy_NH4,phosy_NH4,pddpo,1)
      call acclyr(jphosy_NO3,phosy_NO3,pddpo,1)
      call acclyr(jremin_aerob,remin_aerob,pddpo,1)
      call acclyr(jremin_sulf,remin_sulf,pddpo,1)
    endif
    if (use_M4AGO) then
      ! M4AGO
      call acclyr(jagg_ws,aggregate_diagnostics(1,1,1,kws_agg),pddpo,1)
      call acclyr(jdynvis,aggregate_diagnostics(1,1,1,kdynvis),pddpo,1)
      call acclyr(jagg_stick,aggregate_diagnostics(1,1,1,kstickiness_agg),pddpo,1)
      call acclyr(jagg_stickf,aggregate_diagnostics(1,1,1,kstickiness_frustule),pddpo,1)
      call acclyr(jagg_dmax,aggregate_diagnostics(1,1,1,kLmax_agg),pddpo,1)
      call acclyr(jagg_avdp,aggregate_diagnostics(1,1,1,kav_dp),pddpo,1)
      call acclyr(jagg_avrhop,aggregate_diagnostics(1,1,1,kav_rho_p),pddpo,1)
      call acclyr(jagg_avdC,aggregate_diagnostics(1,1,1,kav_d_C),pddpo,1)
      call acclyr(jagg_df,aggregate_diagnostics(1,1,1,kdf_agg),pddpo,1)
      call acclyr(jagg_b,aggregate_diagnostics(1,1,1,kb_agg),pddpo,1)
      call acclyr(jagg_Vrhof,aggregate_diagnostics(1,1,1,kav_rhof_V),pddpo,1)
      call acclyr(jagg_Vpor,aggregate_diagnostics(1,1,1,kav_por_V),pddpo,1)
    endif
    ! Accumulate level diagnostics
    if (SUM(jlvlphyto+jlvlgrazer+jlvlphosph+jlvloxygen+jlvliron+             &
         &  jlvlano3+jlvlalkali+jlvlsilica+jlvldic+jlvldoc+jlvlpoc+jlvlcalc+ &
         &  jlvlopal+jlvln2o+jlvlco3+jlvlph+jlvlomegaa+jlvlomegac+jlvlphosy+ &
         &  jlvlo2sat+jlvlprefo2+jlvlprefpo4+jlvlprefalk+jlvlprefdic+        &
         &  jlvlprefsilica+jlvlshelfage+                                     &
         &  jlvldicsat+jlvlnatdic+jlvlnatalkali+jlvlnatcalc+jlvlnatco3+      &
         &  jlvlnatomegaa+jlvlnatomegac+jlvldic13+jlvldic14+jlvld13c+        &
         &  jlvld14c+jlvlbigd14c+jlvlpoc13+jlvldoc13+jlvlcalc13+jlvlphyto13+ &
         &  jlvlgrazer13+jlvlnos+jlvlwphy+jlvlwnos+jlvleps+jlvlasize+        &
         &  jlvlcfc11+jlvlcfc12+jlvlsf6+jlvlbromo+jlvlanh4+jlvlano2+        &
         &  jlvl_nitr_NH4+jlvl_nitr_NO2+jlvl_nitr_N2O_prod+jlvl_nitr_NH4_OM+&
         &  jlvl_nitr_NO2_OM+jlvl_denit_NO3+jlvl_denit_NO2+jlvl_denit_N2O+  &
         &  jlvl_DNRA_NO2+jlvl_anmx_N2_prod+jlvl_anmx_OM_prod+              &
         &  jlvl_phosy_NH4+jlvl_phosy_NO3+jlvl_remin_aerob+jlvl_remin_sulf+ &
         &  jlvl_agg_ws+jlvl_dynvis+jlvl_agg_stick+jlvl_agg_stickf+         &
         &  jlvl_agg_dmax+jlvl_agg_avdp+jlvl_agg_avrhop+jlvl_agg_avdC+      &
         &  jlvl_agg_df+jlvl_agg_b+jlvl_agg_Vrhof+jlvl_agg_Vpor             &
         &  ) /= 0) then
      do k=1,kpke
        call bgczlv(pddpo,k,ind1,ind2,wghts)
        call acclvl(jlvlphyto,ocetra(1,1,1,iphy),k,ind1,ind2,wghts)
        call acclvl(jlvlgrazer,ocetra(1,1,1,izoo),k,ind1,ind2,wghts)
        call acclvl(jlvlphosph,ocetra(1,1,1,iphosph),k,ind1,ind2,wghts)
        call acclvl(jlvloxygen,ocetra(1,1,1,ioxygen),k,ind1,ind2,wghts)
        call acclvl(jlvliron,ocetra(1,1,1,iiron),k,ind1,ind2,wghts)
        call acclvl(jlvlano3,ocetra(1,1,1,iano3),k,ind1,ind2,wghts)
        call acclvl(jlvlalkali,ocetra(1,1,1,ialkali),k,ind1,ind2,wghts)
        call acclvl(jlvlsilica,ocetra(1,1,1,isilica),k,ind1,ind2,wghts)
        call acclvl(jlvldic,ocetra(1,1,1,isco212),k,ind1,ind2,wghts)
        call acclvl(jlvldoc,ocetra(1,1,1,idoc),k,ind1,ind2,wghts)
        call acclvl(jlvlpoc,ocetra(1,1,1,idet),k,ind1,ind2,wghts)
        call acclvl(jlvlcalc,ocetra(1,1,1,icalc),k,ind1,ind2,wghts)
        call acclvl(jlvlopal,ocetra(1,1,1,iopal),k,ind1,ind2,wghts)
        call acclvl(jlvln2o,ocetra(1,1,1,ian2o),k,ind1,ind2,wghts)
        call acclvl(jlvlco3,co3,k,ind1,ind2,wghts)
        call acclvl(jlvlph,hi,k,ind1,ind2,wghts)
        call acclvl(jlvlomegaa,OmegaA,k,ind1,ind2,wghts)
        call acclvl(jlvlomegac,OmegaC,k,ind1,ind2,wghts)
        call acclvl(jlvlphosy,phosy3d,k,ind1,ind2,wghts)
        call acclvl(jlvlo2sat,satoxy,k,ind1,ind2,wghts)
        call acclvl(jlvldicsat,ocetra(1,1,1,idicsat),k,ind1,ind2,wghts)
        if (use_pref_tracers) then
          call acclvl(jlvlprefo2,ocetra(1,1,1,iprefo2),k,ind1,ind2,wghts)
          call acclvl(jlvlprefpo4,ocetra(1,1,1,iprefpo4),k,ind1,ind2,wghts)
          call acclvl(jlvlprefsilica,ocetra(1,1,1,iprefsilica),k,ind1,ind2,wghts)
          call acclvl(jlvlprefalk,ocetra(1,1,1,iprefalk),k,ind1,ind2,wghts)
          call acclvl(jlvlprefdic,ocetra(1,1,1,iprefdic),k,ind1,ind2,wghts)
        endif
        if (use_shelfsea_res_time) then
          call acclvl(jlvlshelfage,ocetra(1,1,1,ishelfage),k,ind1,ind2,wghts)
        endif
        if (use_natDIC) then
          call acclvl(jlvlnatdic,ocetra(1,1,1,inatsco212),k,ind1,ind2,wghts)
          call acclvl(jlvlnatalkali,ocetra(1,1,1,inatalkali),k,ind1,ind2,wghts)
          call acclvl(jlvlnatcalc,ocetra(1,1,1,inatcalc),k,ind1,ind2,wghts)
          call acclvl(jlvlnatco3,natco3,k,ind1,ind2,wghts)
          call acclvl(jlvlnatph,nathi,k,ind1,ind2,wghts)
          call acclvl(jlvlnatomegaa,natOmegaA,k,ind1,ind2,wghts)
          call acclvl(jlvlnatomegac,natOmegaC,k,ind1,ind2,wghts)
        endif
        if (use_cisonew) then
          call acclvl(jlvld13c,d13C,k,ind1,ind2,wghts)
          call acclvl(jlvld14c,d14C,k,ind1,ind2,wghts)
          call acclvl(jlvlbigd14c,bigd14C,k,ind1,ind2,wghts)
          call acclvl(jlvldic13,ocetra(1,1,1,isco213),k,ind1,ind2,wghts)
          call acclvl(jlvldic14,ocetra(1,1,1,isco214),k,ind1,ind2,wghts)
          call acclvl(jlvlpoc13,ocetra(1,1,1,idet13),k,ind1,ind2,wghts)
          call acclvl(jlvldoc13,ocetra(1,1,1,idoc13),k,ind1,ind2,wghts)
          call acclvl(jlvlcalc13,ocetra(1,1,1,icalc13),k,ind1,ind2,wghts)
          call acclvl(jlvlphyto13,ocetra(1,1,1,iphy13),k,ind1,ind2,wghts)
          call acclvl(jlvlgrazer13,ocetra(1,1,1,izoo13),k,ind1,ind2,wghts)
        endif
        if (use_AGG) then
          call acclvl(jlvlnos,ocetra(1,1,1,inos),k,ind1,ind2,wghts)
          call acclvl(jlvlwphy, wmass/dtb,k,ind1,ind2,wghts)
          call acclvl(jlvlwnos, wnumb/dtb,k,ind1,ind2,wghts)
          call acclvl(jlvleps,  eps3d,    k,ind1,ind2,wghts)
          call acclvl(jlvlasize,asize3d,  k,ind1,ind2,wghts)
        endif
        if (use_CFC) then
          call acclvl(jlvlcfc11,ocetra(1,1,1,icfc11),k,ind1,ind2,wghts)
          call acclvl(jlvlcfc12,ocetra(1,1,1,icfc12),k,ind1,ind2,wghts)
          call acclvl(jlvlsf6,ocetra(1,1,1,isf6),k,ind1,ind2,wghts)
        endif
        if (use_BROMO) then
          call acclvl(jlvlbromo,ocetra(1,1,1,ibromo),k,ind1,ind2,wghts)
        endif
        if (use_extNcycle) then
          call acclvl(jlvlanh4,ocetra(1,1,1,ianh4),k,ind1,ind2,wghts)
          call acclvl(jlvlano2,ocetra(1,1,1,iano2),k,ind1,ind2,wghts)
          call acclvl(jlvl_nitr_NH4,nitr_NH4,k,ind1,ind2,wghts)
          call acclvl(jlvl_nitr_NO2,nitr_NO2,k,ind1,ind2,wghts)
          call acclvl(jlvl_nitr_N2O_prod,nitr_N2O_prod,k,ind1,ind2,wghts)
          call acclvl(jlvl_nitr_NH4_OM,nitr_NH4_OM,k,ind1,ind2,wghts)
          call acclvl(jlvl_nitr_NO2_OM,nitr_NO2_OM,k,ind1,ind2,wghts)
          call acclvl(jlvl_denit_NO3,denit_NO3,k,ind1,ind2,wghts)
          call acclvl(jlvl_denit_NO2,denit_NO2,k,ind1,ind2,wghts)
          call acclvl(jlvl_denit_N2O,denit_N2O,k,ind1,ind2,wghts)
          call acclvl(jlvl_DNRA_NO2,DNRA_NO2,k,ind1,ind2,wghts)
          call acclvl(jlvl_anmx_N2_prod,anmx_N2_prod,k,ind1,ind2,wghts)
          call acclvl(jlvl_anmx_OM_prod,anmx_OM_prod,k,ind1,ind2,wghts)
          call acclvl(jlvl_phosy_NH4,phosy_NH4,k,ind1,ind2,wghts)
          call acclvl(jlvl_phosy_NO3,phosy_NO3,k,ind1,ind2,wghts)
          call acclvl(jlvl_remin_aerob,remin_aerob,k,ind1,ind2,wghts)
          call acclvl(jlvl_remin_sulf,remin_sulf,k,ind1,ind2,wghts)
        endif
        if (use_M4AGO) then
          !M4AGO
          call acclvl(jlvl_agg_ws,aggregate_diagnostics(1,1,1,kws_agg),k,ind1,ind2,wghts)
          call acclvl(jlvl_dynvis,aggregate_diagnostics(1,1,1,kdynvis),k,ind1,ind2,wghts)
          call acclvl(jlvl_agg_stick,aggregate_diagnostics(1,1,1,kstickiness_agg),k,ind1,ind2,wghts)
          call acclvl(jlvl_agg_stickf,aggregate_diagnostics(1,1,1,kstickiness_frustule),k,ind1,ind2,wghts)
          call acclvl(jlvl_agg_dmax,aggregate_diagnostics(1,1,1,kLmax_agg),k,ind1,ind2,wghts)
          call acclvl(jlvl_agg_avdp,aggregate_diagnostics(1,1,1,kav_dp),k,ind1,ind2,wghts)
          call acclvl(jlvl_agg_avrhop,aggregate_diagnostics(1,1,1,kav_rho_p),k,ind1,ind2,wghts)
          call acclvl(jlvl_agg_avdC,aggregate_diagnostics(1,1,1,kav_d_C),k,ind1,ind2,wghts)
          call acclvl(jlvl_agg_df,aggregate_diagnostics(1,1,1,kdf_agg),k,ind1,ind2,wghts)
          call acclvl(jlvl_agg_b,aggregate_diagnostics(1,1,1,kb_agg),k,ind1,ind2,wghts)
          call acclvl(jlvl_agg_Vrhof,aggregate_diagnostics(1,1,1,kav_rhof_V),k,ind1,ind2,wghts)
          call acclvl(jlvl_agg_Vpor,aggregate_diagnostics(1,1,1,kav_por_V),k,ind1,ind2,wghts)
        endif
      enddo
    endif

    if (.not. use_sedbypass) then
      ! Accumulate sediments
      call accsdm(jpowaic,powtra(1,1,1,ipowaic))
      call accsdm(jpowaal,powtra(1,1,1,ipowaal))
      call accsdm(jpowaph,powtra(1,1,1,ipowaph))
      call accsdm(jpowaox,powtra(1,1,1,ipowaox))
      call accsdm(jpown2 ,powtra(1,1,1,ipown2) )
      call accsdm(jpowno3,powtra(1,1,1,ipowno3))
      call accsdm(jpowasi,powtra(1,1,1,ipowasi))
      call accsdm(jssso12,sedlay(1,1,1,issso12))
      call accsdm(jssssil,sedlay(1,1,1,issssil))
      call accsdm(jsssc12,sedlay(1,1,1,isssc12))
      call accsdm(jssster,sedlay(1,1,1,issster))

      ! Accumulate sediment burial
      call accbur(jburssso12,burial(1,1,issso12))
      call accbur(jburssssil,burial(1,1,issssil))
      call accbur(jbursssc12,burial(1,1,isssc12))
      call accbur(jburssster,burial(1,1,issster))
      if (use_extNcycle) then
        call accsdm(jpownh4,powtra(1,1,1,ipownh4))
        call accsdm(jpown2o,powtra(1,1,1,ipown2o))
        call accsdm(jpowno2,powtra(1,1,1,ipowno2))
        call accsdm(jsdm_nitr_NH4      ,extNsed_diagnostics(1,1,1,ised_nitr_NH4))
        call accsdm(jsdm_nitr_NO2      ,extNsed_diagnostics(1,1,1,ised_nitr_NO2))
        call accsdm(jsdm_nitr_N2O_prod ,extNsed_diagnostics(1,1,1,ised_nitr_N2O_prod))
        call accsdm(jsdm_nitr_NH4_OM   ,extNsed_diagnostics(1,1,1,ised_nitr_NH4_OM))
        call accsdm(jsdm_nitr_NO2_OM   ,extNsed_diagnostics(1,1,1,ised_nitr_NO2_OM))
        call accsdm(jsdm_denit_NO3     ,extNsed_diagnostics(1,1,1,ised_denit_NO3))
        call accsdm(jsdm_denit_NO2     ,extNsed_diagnostics(1,1,1,ised_denit_NO2))
        call accsdm(jsdm_denit_N2O     ,extNsed_diagnostics(1,1,1,ised_denit_N2O))
        call accsdm(jsdm_DNRA_NO2      ,extNsed_diagnostics(1,1,1,ised_DNRA_NO2))
        call accsdm(jsdm_anmx_N2_prod  ,extNsed_diagnostics(1,1,1,ised_anmx_N2_prod))
        call accsdm(jsdm_anmx_OM_prod  ,extNsed_diagnostics(1,1,1,ised_anmx_OM_prod))
        call accsdm(jsdm_remin_aerob   ,extNsed_diagnostics(1,1,1,ised_remin_aerob))
        call accsdm(jsdm_remin_sulf    ,extNsed_diagnostics(1,1,1,ised_remin_sulf))
      endif
    endif

    ! Write output if requested
    do l=1,nbgc
      nacc_bgc(l)=nacc_bgc(l)+1
      if (bgcwrt(l)) then
        if (GLB_INVENTORY(l).ne.0) then
          call INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,l)
        endif
        call ncwrt_bgc(l)
        nacc_bgc(l)=0
      endif
    enddo

    atmflx=0. ! nullifying atm flux here to have zero fluxes for stepwise inventory fluxes
    ndepnoyflx=0.
    oalkflx=0.
    rivinflx=0.
    if (use_extNcycle) then
      ndepnhxflx=0.
    endif

  end subroutine accfields

end module mo_accfields
