! Copyright (C) 2020  I Bethke, J. Tjiputra, J. Schwinger, A. Moree, M.
!                     Bentsen
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


subroutine ncwrt_bgc(iogrp)
!
! --- -------------------------------------------
! --- output routine for HAMOCC diagnostic fields
! --- -------------------------------------------
!
  use mod_time,       only: date0,date,calendar,nstep,nstep_in_day,             &
       &                    nday_of_year,time0,time
  use mod_xc,         only: kdm,mnproc,itdm,jtdm,lp
  use mod_grid,       only: depths
  use mod_dia,        only: diafnm,sigmar1,iotype,ddm,depthslev,                &
       &                    depthslev_bnds
  use mo_control_bgc, only: dtbgc
  use mo_vgrid,       only: k0100,k0500,k1000,k2000,k4000
  use mo_param1_bgc,  only: ks
  use mod_nctools,    only: ncwrt1,ncdims,nctime,ncfcls,ncfopn,                 &
       &                    ncdimc
  use mo_bgcmean,     only: domassfluxes,                                       &
       &                    flx_cal0100,flx_cal0500,flx_cal1000,                &
       &                    flx_cal2000,flx_cal4000,flx_cal_bot,                &
       &                    flx_car0100,flx_car0500,flx_car1000,                &
       &                    flx_car2000,flx_car4000,flx_car_bot,                &
       &                    flx_bsi0100,flx_bsi0500,flx_bsi1000,                &
       &                    flx_bsi2000,flx_bsi4000,flx_bsi_bot,                &
       &                    flx_sediffic,flx_sediffal,flx_sediffph,             &
       &                    flx_sediffox,flx_sediffn2,flx_sediffno3,            &
       &                    flx_sediffsi,                                       &
       &                    jsediffic,jsediffal,jsediffph,jsediffox,            &
       &                    jsediffn2,jsediffno3,jsediffsi,                     &
       &                    jalkali,jano3,jasize,jatmco2,                       &
       &                    jbsiflx0100,jbsiflx0500,jbsiflx1000,                &
       &                    jbsiflx2000,jbsiflx4000,jbsiflx_bot,                &
       &                    jcalc,jcalflx0100,jcalflx0500,jcalflx1000,          &
       &                    jcalflx2000,jcalflx4000,jcalflx_bot,                &
       &                    jcarflx0100,jcarflx0500,jcarflx1000,                &
       &                    jcarflx2000,jcarflx4000,jcarflx_bot,                &
       &                    jco2flux,jco2fxd,jco2fxu,jco3,jdic,jdicsat,         &
       &                    jdms,jdms_bac,jdms_uv,jdmsflux,jdmsprod,            &
       &                    jdoc,jdp,jeps,jexpoca,jexport,jexposi,              &
       &                    jgrazer,                                            &
       &                    jintdnit,jintnfix,jintphosy,jiralk,jirdet,          &
       &                    jirdin,jirdip,jirdoc,jiriron,jiron,jirsi,           &
       &                    jkwco2,jlvlalkali,jlvlano3,jlvlasize,               &
       &                    jlvlbigd14c,jlvlbromo,jlvlcalc,jlvlcalc13,          &
       &                    jlvlcfc11,jlvlcfc12,jlvlco3,jlvld13c,               &
       &                    jlvld14c,jlvldic,jlvldic13,jlvldic14,               &
       &                    jlvldicsat,jlvldoc,jlvldoc13,jlvleps,               &
       &                    jlvlgrazer,jlvlgrazer13,jlvliron,jlvln2o,           &
       &                    jlvlnatalkali,jlvlnatcalc,jlvlnatco3,               &
       &                    jlvlnatdic,jlvlnatomegaa,jlvlnatomegac,             &
       &                    jlvlnos,jlvlo2sat,jlvlomegaa,jlvlomegac,            &
       &                    jlvlopal,jlvloxygen,jlvlph,jlvlphosph,              &
       &                    jlvlphosy,jlvlphyto,jlvlphyto13,jlvlpoc,            &
       &                    jlvlpoc13,jlvlprefalk,jlvlprefdic,                  &
       &                    jlvlprefo2,jlvlprefpo4,jlvlsf6,jlvlsilica,          &
       &                    jlvlwnos,jlvlwphy,jn2flux,jn2o,jn2oflux,            &
       &                    jn2ofx,jndep,jniflux,jnos,jo2flux,jo2sat,           &
       &                    jomegaa,jomegac,jopal,joxflux,joxygen,jpco2,        &
       &                    jph,jphosph,jphosy,jphyto,jpoc,jprefalk,            &
       &                    jprefdic,jprefo2,jprefpo4,jsilica,                  &
       &                    jsrfalkali,jsrfano3,jsrfdic,jsrfiron,               &
       &                    jsrfoxygen,jsrfphosph,jsrfphyto,jsrfsilica,         &
       &                    jwnos,jwphy,                                        &
       &                    lyr_dp,lyr_dic,lyr_alkali,lyr_phosph,               &
       &                    lyr_oxygen,lyr_ano3,lyr_silica,lyr_doc,             &
       &                    lyr_phyto,lyr_grazer,lyr_poc,lyr_calc,              &
       &                    lyr_opal,lyr_iron,lyr_phosy,lyr_co3,lyr_ph,         &
       &                    lyr_omegaa,lyr_omegac,lyr_n2o,lyr_prefo2,           &
       &                    lyr_o2sat,lyr_prefpo4,lyr_prefalk,                  &
       &                    lyr_prefdic,lyr_dicsat,                             &
       &                    lvl_dic,lvl_alkali,                                 &
       &                    lvl_phosph,lvl_oxygen,lvl_ano3,lvl_silica,          &
       &                    lvl_doc,lvl_phyto,lvl_grazer,lvl_poc,               &
       &                    lvl_calc,lvl_opal,lvl_iron,lvl_phosy,               &
       &                    lvl_co3,lvl_ph,lvl_omegaa,lvl_omegac,               &
       &                    lvl_n2o,lvl_prefo2,lvl_o2sat,lvl_prefpo4,           &
       &                    lvl_prefalk,lvl_prefdic,lvl_dicsat,                 &
       &                    lvl_o2sat,srf_n2ofx,srf_atmco2,srf_kwco2,           &
       &                    srf_pco2,srf_dmsflux,srf_co2fxd,                    &
       &                    srf_co2fxu,srf_oxflux,srf_niflux,srf_dms,           &
       &                    srf_dmsprod,srf_dms_bac,srf_dms_uv,                 &
       &                    srf_export,srf_exposi,srf_expoca,srf_dic,           &
       &                    srf_alkali,srf_phosph,srf_oxygen,srf_ano3,          &
       &                    srf_silica,srf_iron,srf_phyto,                      &
       &                    int_phosy,int_nfix,int_dnit,                        &
       &                    nbgc,nacc_bgc,bgcwrt,glb_inventory,bgct2d,          &
       &                    nbgcmax,glb_ncformat,glb_compflag,                  &
       &                    glb_fnametag,filefq_bgc,diagfq_bgc,                 &
       &                    filemon_bgc,fileann_bgc,ip,wrtlyr,wrtlvl,           &
       &                    loglyr,inilvl,inilyr,inisrf,loglvl,                 &
       &                    msklvl,wrtsrf,msksrf,finlyr,                        &
       &                    lyr_agg_ws,lyr_dynvis,lyr_agg_stick,                &
       &                    lyr_agg_stickf,lyr_agg_dmax,lyr_agg_avdp,           &
       &                    lyr_agg_avrhop,lyr_agg_avdC,lyr_agg_df,             &
       &                    lyr_agg_b,lyr_agg_Vrhof,lyr_agg_Vpor,               &
       &                    lvl_agg_ws,lvl_dynvis,lvl_agg_stick,                &
       &                    lvl_agg_stickf,lvl_agg_dmax,lvl_agg_avdp,           &
       &                    lvl_agg_avrhop,lvl_agg_avdC,lvl_agg_df,             &
       &                    lvl_agg_b,lvl_agg_Vrhof,lvl_agg_Vpor,               &
       &                    jagg_ws,jdynvis,jagg_stick,                         &
       &                    jagg_stickf,jagg_dmax,jagg_avdp,                    &
       &                    jagg_avrhop,jagg_avdC,jagg_df,                      &
       &                    jagg_b,jagg_Vrhof,jagg_Vpor,                        &
       &                    jlvl_agg_ws,jlvl_dynvis,jlvl_agg_stick,             &
       &                    jlvl_agg_stickf,jlvl_agg_dmax,jlvl_agg_avdp,        &
       &                    jlvl_agg_avrhop,jlvl_agg_avdC,jlvl_agg_df,          &
       &                    jlvl_agg_b,jlvl_agg_Vrhof,jlvl_agg_Vpor 
#ifdef AGG
  use mo_bgcmean, only: lyr_nos,lyr_wphy, lyr_wnos,lyr_eps,                     &
       &                lyr_asize,lvl_nos,lvl_wphy,lvl_wnos,lvl_eps,            &
       &                lvl_asize
#endif
#ifdef BROMO
  use mo_bgcmean, only: jbromo,jbromofx,jsrfbromo,jbromo_prod,                  &
       &                jbromo_uv,jatmbromo,lvl_bromo,srf_bromofx,              &
       &                srf_bromo,int_bromopro,int_bromouv,                     &
       &                srf_atmbromo,lyr_bromo
#endif
#ifdef CFC
  use mo_bgcmean,only: jcfc11,jcfc12,jsf6,jcfc11fx,jcfc12fx,jsf6fx,             &
       &               lvl_cfc11,lvl_cfc12,lvl_sf6,srf_cfc11,                   &
       &               srf_cfc12,srf_sf6,lyr_cfc11,lyr_cfc12,                   &
       &               lyr_sf6
#endif
#ifdef cisonew
  use mo_biomod, only: c14fac
  use mo_bgcmean, only: jdic13,jdic14,jd13c,jd14c,jbigd14c,jpoc13,              &
       &                jdoc13,jcalc13,jphyto13,jgrazer13,jco213fxd,            &
       &                jco213fxu,jco214fxd,jco214fxu,jatmc13,                  &
       &                jatmc14,jdic13,jdic14,jd13c,jd14c,jbigd14c,             &
       &                srf_co213fxd,srf_co213fxu,srf_co214fxd,                 &
       &                srf_co214fxu,srf_atmc13,srf_atmc14,lyr_dic13,           &
       &                lyr_dic14,lyr_d13c,lyr_d14c,lyr_bigd14c,                &
       &                lyr_poc13,lyr_doc13,lyr_calc13,lyr_phyto13,             &
       &                lyr_grazer13,lvl_dic13,lvl_dic14,lvl_d13c,              &
       &                lvl_d14c,lvl_bigd14c,lvl_poc13,lvl_doc13,               &
       &                lvl_calc13,lvl_phyto13,lvl_grazer13
#endif
#ifdef natDIC
  use mo_bgcmean, only: jnatalkali,jnatdic,jnatcalc,jnatco3,jnatph,             &
       &                jnatomegaa,jnatomegac,lyr_natph,jlvlnatph,              &
       &                lvl_natph,jsrfnatdic,                                   &
       &                jsrfnatalk,jnatpco2,jnatco2fx,lyr_natco3,               &
       &                lyr_natalkali,lyr_natdic,lyr_natcalc,                   &
       &                lyr_natomegaa,lyr_natomegac,lvl_natco3,                 &
       &                lvl_natalkali,lvl_natdic,lvl_natcalc,                   &
       &                lvl_natomegaa,lvl_natomegac,srf_natdic,                 &
       &                srf_natalkali,srf_natpco2,srf_natco2fx
#endif
#ifndef sedbypass
  use mo_bgcmean, only: jpowaic,jpowaal,jpowaph,jpowaox,jpown2,                 &
       &                jpowno3,jpowasi,jssso12,jssssil,jssster,                &
       &                jsssc12,jbursssc12,jburssssil,jburssster,               &
       &                sdm_powaic,sdm_powaal,sdm_powaph,sdm_powaox,            &
       &                sdm_pown2,sdm_powno3,sdm_powasi,sdm_ssso12,             &
       &                sdm_ssssil,sdm_sssc12,sdm_ssster,jburssso12,            &
       &                bur_sssc12,bur_ssssil,bur_ssster,bur_ssso12,            &
       &                inisdm,inibur,wrtsdm,accbur,accsdm,wrtbur
#endif
#ifdef extNcycle
  use mo_bgcmean, only: janh4,jano2,jlvlanh4,jlvlano2,jsrfanh4,                 &
       &                jsrfano2,janh3fx,srf_anh4,srf_ano2,                     &
       &                srf_anh3fx,lyr_anh4,lyr_ano2,lvl_anh4,                  &
       &                lvl_ano2,                                               &
       &                LYR_nitr_NH4,LYR_nitr_NO2,LYR_nitr_N2O_prod,            &
       &                LYR_nitr_NH4_OM,LYR_nitr_NO2_OM,                        &
       &                LYR_denit_NO3,LYR_denit_NO2,LYR_denit_N2O,              &
       &                LYR_DNRA_NO2,LYR_anmx_N2_prod,                          &
       &                LYR_anmx_OM_prod,LYR_phosy_NH4,                         &
       &                LYR_phosy_NO3,LYR_remin_aerob,LYR_remin_sulf,           &
       &                LVL_nitr_NH4,LVL_nitr_NO2,LVL_nitr_N2O_prod,            &
       &                LVL_nitr_NH4_OM,LVL_nitr_NO2_OM,                        &
       &                LVL_denit_NO3,LVL_denit_NO2,LVL_denit_N2O,              &
       &                LVL_DNRA_NO2,LVL_anmx_N2_prod,                          &
       &                LVL_anmx_OM_prod,LVL_phosy_NH4,                         &
       &                LVL_phosy_NO3,LVL_remin_aerob,LVL_remin_sulf,           &
       &                jnitr_NH4,jnitr_NO2,jnitr_N2O_prod,                     &
       &                jnitr_NH4_OM,jnitr_NO2_OM,jdenit_NO3,                   &
       &                jdenit_NO2,jdenit_N2O,jDNRA_NO2,                        &
       &                janmx_N2_prod,janmx_OM_prod,jphosy_NH4,                 &
       &                jphosy_NO3,jremin_aerob,jremin_sulf,                    &
       &                jlvl_nitr_NH4,jlvl_nitr_NO2,                            &
       &                jlvl_nitr_N2O_prod,jlvl_nitr_NH4_OM,                    &
       &                jlvl_nitr_NO2_OM,jlvl_denit_NO3,                        &
       &                jlvl_denit_NO2,jlvl_denit_N2O,jlvl_DNRA_NO2,            &
       &                jlvl_anmx_N2_prod,jlvl_anmx_OM_prod,                    &
       &                jlvl_phosy_NH4,jlvl_phosy_NO3,                          &
       &                jlvl_remin_aerob,jlvl_remin_sulf  
#endif
#if defined(extNcycle) && ! defined(sedbypass)
  use mo_bgcmean, only: jpownh4,jpown2o,jpowno2,jsdm_nitr_NH4,jsdm_nitr_NO2,    &
       &                jsdm_nitr_N2O_prod,jsdm_nitr_NH4_OM,jsdm_nitr_NO2_OM,   &
       &                jsdm_denit_NO3,jsdm_denit_NO2,jsdm_denit_N2O,           &
       &                jsdm_DNRA_NO2,jsdm_anmx_N2_prod,jsdm_anmx_OM_prod,      &
       &                jsdm_remin_aerob,jsdm_remin_sulf, SDM_POWNH4,SDM_POWN2O,&
       &                SDM_POWNO2,SDM_nitr_NH4,SDM_nitr_NO2,SDM_nitr_N2O_prod, &
       &                SDM_nitr_NH4_OM,SDM_nitr_NO2_OM,SDM_denit_NO3,          &
       &                SDM_denit_NO2,SDM_denit_N2O,SDM_DNRA_NO2,               &
       &                SDM_anmx_N2_prod,SDM_anmx_OM_prod,SDM_remin_aerob,      &
       &                SDM_remin_sulf,jsediffnh4,jsediffn2o,jsediffno2,        &
       &                FLX_SEDIFFNH4,FLX_SEDIFFN2O,FLX_SEDIFFNO2 
#endif
  implicit none

  integer iogrp

  integer i,j,k,l,nt
  integer ny,nm,nd,dayfrac,irec(nbgcmax),cmpflg
  character*256 fname(nbgcmax)
  character startdate*20,timeunits*30
  real datenum,rnacc
  logical append2file(nbgcmax)
  data append2file /nbgcmax*.false./
  save fname,irec,append2file

  ! --- set time information
  timeunits=' '
  startdate=' '
  write(timeunits,'(a11,i4.4,a1,i2.2,a1,i2.2,a6)')                              &
       &   'days since ',min(1800,date0%year),'-',1,'-',1,' 00:00'
  write(startdate,'(i4.4,a1,i2.2,a1,i2.2,a6)')                                  &
       &   date0%year,'-',date0%month,'-',date0%day,' 00:00'
  datenum=time-time0-0.5*diagfq_bgc(iogrp)/nstep_in_day

  ! --- get file name
  if (.not.append2file(iogrp)) then
     call diafnm(GLB_FNAMETAG(iogrp),                                           &
          &   filefq_bgc(iogrp)/real(nstep_in_day),                             &
          &   filemon_bgc(iogrp),fileann_bgc(iogrp),fname(iogrp))
     append2file(iogrp)=.true.
     irec(iogrp)=1
  else
     irec(iogrp)=irec(iogrp)+1
  endif
  if (((fileann_bgc(iogrp).and.nday_of_year.eq.1.or.                            &
       &   filemon_bgc(iogrp).and.date%day.eq.1).and.                           &
       &   mod(nstep,nstep_in_day).eq.0).or.                                    &
       &   .not.(fileann_bgc(iogrp).or.filemon_bgc(iogrp)).and.                 &
       &   mod(nstep+.5,filefq_bgc(iogrp)).lt.1.) then
     append2file(iogrp)=.false.
  endif

  ! --- prepare output fields
  if (mnproc.eq.1) then
     write (lp,'(a,f6.2,a)') ' ncwrt_bgc: fields averaged over ',               &
          &   real(nacc_bgc(iogrp)),' steps'
     write(lp,*) 'irec(iogrp)',irec(iogrp)
  endif
  rnacc=1./real(nacc_bgc(iogrp))
  cmpflg=GLB_COMPFLAG(iogrp)

  ! --- create output file
  if (GLB_NCFORMAT(iogrp).eq.1) then
     call ncfopn(fname(iogrp),'w','6',irec(iogrp),iotype)
  elseif (GLB_NCFORMAT(iogrp).eq.2) then
     call ncfopn(fname(iogrp),'w','h',irec(iogrp),iotype)
  else
     call ncfopn(fname(iogrp),'w','c',irec(iogrp),iotype)
  endif

  ! --- define spatial and time dimensions
  if (cmpflg.ne.0) then
     call ncdimc('pcomp',ip,0)
  else
     call ncdims('x',itdm)
     call ncdims('y',jtdm)
  endif
  call ncdims('sigma',kdm)
  call ncdims('depth',ddm)
  call ncdims('ks',ks)
  call ncdims('bounds',2)
  call ncdims('time',0)
  call hamoccvardef(iogrp,timeunits,calendar,cmpflg)
  call nctime(datenum,calendar,timeunits,startdate)

  ! --- write auxillary dimension information
  call ncwrt1('sigma','sigma',sigmar1)
  call ncwrt1('depth','depth',depthslev)
  call ncwrt1('depth_bnds','bounds depth',depthslev_bnds)

  ! --- finalize accumulation
  call finlyr(jphyto(iogrp),jdp(iogrp))
  call finlyr(jgrazer(iogrp),jdp(iogrp))
  call finlyr(jdoc(iogrp),jdp(iogrp))
  call finlyr(jphosy(iogrp),jdp(iogrp))
  call finlyr(jphosph(iogrp),jdp(iogrp))
  call finlyr(joxygen(iogrp),jdp(iogrp))
  call finlyr(jiron(iogrp),jdp(iogrp))
  call finlyr(jano3(iogrp),jdp(iogrp))
  call finlyr(jalkali(iogrp),jdp(iogrp))
  call finlyr(jsilica(iogrp),jdp(iogrp))
  call finlyr(jdic(iogrp),jdp(iogrp))
  call finlyr(jpoc(iogrp),jdp(iogrp))
  call finlyr(jcalc(iogrp),jdp(iogrp))
  call finlyr(jopal(iogrp),jdp(iogrp))
  call finlyr(jco3(iogrp),jdp(iogrp))
  call finlyr(jph(iogrp),jdp(iogrp))
  call finlyr(jomegaa(iogrp),jdp(iogrp))
  call finlyr(jomegac(iogrp),jdp(iogrp))
  call finlyr(jn2o(iogrp),jdp(iogrp))
  call finlyr(jprefo2(iogrp),jdp(iogrp))
  call finlyr(jo2sat(iogrp),jdp(iogrp))
  call finlyr(jprefpo4(iogrp),jdp(iogrp))
  call finlyr(jprefalk(iogrp),jdp(iogrp))
  call finlyr(jprefdic(iogrp),jdp(iogrp))
  call finlyr(jdicsat(iogrp),jdp(iogrp))
#ifdef cisonew
  call finlyr(jdic13(iogrp),jdp(iogrp))
  call finlyr(jdic14(iogrp),jdp(iogrp))
  call finlyr(jd13c(iogrp),jdp(iogrp))
  call finlyr(jd14c(iogrp),jdp(iogrp))
  call finlyr(jbigd14c(iogrp),jdp(iogrp))
  call finlyr(jpoc13(iogrp),jdp(iogrp))
  call finlyr(jdoc13(iogrp),jdp(iogrp))
  call finlyr(jcalc13(iogrp),jdp(iogrp))
  call finlyr(jphyto13(iogrp),jdp(iogrp))
  call finlyr(jgrazer13(iogrp),jdp(iogrp))
#endif
#ifdef AGG
  call finlyr(jnos(iogrp),jdp(iogrp))
  call finlyr(jwphy(iogrp),jdp(iogrp))
  call finlyr(jwnos(iogrp),jdp(iogrp))
  call finlyr(jeps(iogrp),jdp(iogrp))
  call finlyr(jasize(iogrp),jdp(iogrp))
#endif
#ifdef CFC
  call finlyr(jcfc11(iogrp),jdp(iogrp))
  call finlyr(jcfc12(iogrp),jdp(iogrp))
  call finlyr(jsf6(iogrp),jdp(iogrp))
#endif
#ifdef natDIC
  call finlyr(jnatalkali(iogrp),jdp(iogrp))
  call finlyr(jnatdic(iogrp),jdp(iogrp))
  call finlyr(jnatcalc(iogrp),jdp(iogrp))
  call finlyr(jnatco3(iogrp),jdp(iogrp))
  call finlyr(jnatph(iogrp),jdp(iogrp))
  call finlyr(jnatomegaa(iogrp),jdp(iogrp))
  call finlyr(jnatomegac(iogrp),jdp(iogrp))
#endif
#ifdef BROMO
  call finlyr(jbromo(iogrp),jdp(iogrp))
#endif
#ifdef extNcycle
  call finlyr(janh4(iogrp),jdp(iogrp))
  call finlyr(jano2(iogrp),jdp(iogrp))
  call finlyr(jnitr_NH4(iogrp),jdp(iogrp))
  call finlyr(jnitr_NO2(iogrp),jdp(iogrp))
  call finlyr(jnitr_N2O_prod(iogrp),jdp(iogrp))
  call finlyr(jnitr_NH4_OM(iogrp),jdp(iogrp))
  call finlyr(jnitr_NO2_OM(iogrp),jdp(iogrp))
  call finlyr(jdenit_NO3(iogrp),jdp(iogrp))
  call finlyr(jdenit_NO2(iogrp),jdp(iogrp))
  call finlyr(jdenit_N2O(iogrp),jdp(iogrp))
  call finlyr(jDNRA_NO2(iogrp),jdp(iogrp))
  call finlyr(janmx_N2_prod(iogrp),jdp(iogrp))
  call finlyr(janmx_OM_prod(iogrp),jdp(iogrp))
  call finlyr(jphosy_NH4(iogrp),jdp(iogrp))
  call finlyr(jphosy_NO3(iogrp),jdp(iogrp))
  call finlyr(jremin_aerob(iogrp),jdp(iogrp))
  call finlyr(jremin_sulf(iogrp),jdp(iogrp))
#endif
  !  M4AGO
  call finlyr(jagg_ws(iogrp),jdp(iogrp))
  call finlyr(jdynvis(iogrp),jdp(iogrp))
  call finlyr(jagg_stick(iogrp),jdp(iogrp))
  call finlyr(jagg_stickf(iogrp),jdp(iogrp))
  call finlyr(jagg_dmax(iogrp),jdp(iogrp))
  call finlyr(jagg_avdp(iogrp),jdp(iogrp))
  call finlyr(jagg_avrhop(iogrp),jdp(iogrp))
  call finlyr(jagg_avdC(iogrp),jdp(iogrp))
  call finlyr(jagg_df(iogrp),jdp(iogrp))
  call finlyr(jagg_b(iogrp),jdp(iogrp))
  call finlyr(jagg_Vrhof(iogrp),jdp(iogrp))
  call finlyr(jagg_Vpor(iogrp),jdp(iogrp))

  ! --- Mask sea floor in mass fluxes
  call msksrf(jcarflx0100(iogrp),k0100)
  call msksrf(jcarflx0500(iogrp),k0500)
  call msksrf(jcarflx1000(iogrp),k1000)
  call msksrf(jcarflx2000(iogrp),k2000)
  call msksrf(jcarflx4000(iogrp),k4000)
  call msksrf(jbsiflx0100(iogrp),k0100)
  call msksrf(jbsiflx0500(iogrp),k0500)
  call msksrf(jbsiflx1000(iogrp),k1000)
  call msksrf(jbsiflx2000(iogrp),k2000)
  call msksrf(jbsiflx4000(iogrp),k4000)
  call msksrf(jcalflx0100(iogrp),k0100)
  call msksrf(jcalflx0500(iogrp),k0500)
  call msksrf(jcalflx1000(iogrp),k1000)
  call msksrf(jcalflx2000(iogrp),k2000)
  call msksrf(jcalflx4000(iogrp),k4000)

  ! --- Mask sea floor in level data
  call msklvl(jlvlphyto(iogrp),depths)
  call msklvl(jlvlgrazer(iogrp),depths)
  call msklvl(jlvldoc(iogrp),depths)
  call msklvl(jlvlphosy(iogrp),depths)
  call msklvl(jlvlphosph(iogrp),depths)
  call msklvl(jlvloxygen(iogrp),depths)
  call msklvl(jlvliron(iogrp),depths)
  call msklvl(jlvlano3(iogrp),depths)
  call msklvl(jlvlalkali(iogrp),depths)
  call msklvl(jlvlsilica(iogrp),depths)
  call msklvl(jlvldic(iogrp),depths)
  call msklvl(jlvlpoc(iogrp),depths)
  call msklvl(jlvlcalc(iogrp),depths)
  call msklvl(jlvlopal(iogrp),depths)
  call msklvl(jlvlco3(iogrp),depths)
  call msklvl(jlvlph(iogrp),depths)
  call msklvl(jlvlomegaa(iogrp),depths)
  call msklvl(jlvlomegac(iogrp),depths)
  call msklvl(jlvln2o(iogrp),depths)
  call msklvl(jlvlprefo2(iogrp),depths)
  call msklvl(jlvlo2sat(iogrp),depths)
  call msklvl(jlvlprefpo4(iogrp),depths)
  call msklvl(jlvlprefalk(iogrp),depths)
  call msklvl(jlvlprefdic(iogrp),depths)
  call msklvl(jlvldicsat(iogrp),depths)
#ifdef cisonew
  call msklvl(jlvldic13(iogrp),depths)
  call msklvl(jlvldic14(iogrp),depths)
  call msklvl(jlvld13c(iogrp),depths)
  call msklvl(jlvld14c(iogrp),depths)
  call msklvl(jlvlbigd14c(iogrp),depths)
  call msklvl(jlvlpoc13(iogrp),depths)
  call msklvl(jlvldoc13(iogrp),depths)
  call msklvl(jlvlcalc13(iogrp),depths)
  call msklvl(jlvlphyto13(iogrp),depths)
  call msklvl(jlvlgrazer13(iogrp),depths)
#endif
#ifdef AGG
  call msklvl(jlvlnos(iogrp),depths)
  call msklvl(jlvlwphy(iogrp),depths)
  call msklvl(jlvlwnos(iogrp),depths)
  call msklvl(jlvleps(iogrp),depths)
  call msklvl(jlvlasize(iogrp),depths)
#endif
#ifdef CFC
  call msklvl(jlvlcfc11(iogrp),depths)
  call msklvl(jlvlcfc12(iogrp),depths)
  call msklvl(jlvlsf6(iogrp),depths)
#endif
#ifdef natDIC
  call msklvl(jlvlnatalkali(iogrp),depths)
  call msklvl(jlvlnatdic(iogrp),depths)
  call msklvl(jlvlnatcalc(iogrp),depths)
  call msklvl(jlvlnatco3(iogrp),depths)
  call msklvl(jlvlnatph(iogrp),depths)
  call msklvl(jlvlnatomegaa(iogrp),depths)
  call msklvl(jlvlnatomegac(iogrp),depths)
#endif
#ifdef BROMO
  call msklvl(jlvlbromo(iogrp),depths)
#endif
#ifdef extNcycle
  call msklvl(jlvlanh4(iogrp),depths)
  call msklvl(jlvlano2(iogrp),depths)
  call msklvl(jlvl_nitr_NH4(iogrp),depths)
  call msklvl(jlvl_nitr_NO2(iogrp),depths)
  call msklvl(jlvl_nitr_N2O_prod(iogrp),depths)
  call msklvl(jlvl_nitr_NH4_OM(iogrp),depths)
  call msklvl(jlvl_nitr_NO2_OM(iogrp),depths)
  call msklvl(jlvl_denit_NO3(iogrp),depths)
  call msklvl(jlvl_denit_NO2(iogrp),depths)
  call msklvl(jlvl_denit_N2O(iogrp),depths)
  call msklvl(jlvl_DNRA_NO2(iogrp),depths)
  call msklvl(jlvl_anmx_N2_prod(iogrp),depths)
  call msklvl(jlvl_anmx_OM_prod(iogrp),depths)
  call msklvl(jlvl_phosy_NH4(iogrp),depths)
  call msklvl(jlvl_phosy_NO3(iogrp),depths)
  call msklvl(jlvl_remin_aerob(iogrp),depths)
  call msklvl(jlvl_remin_sulf(iogrp),depths)
#endif
  !   M4AGO
  call msklvl(jlvl_agg_ws(iogrp),depths)
  call msklvl(jlvl_dynvis(iogrp),depths)
  call msklvl(jlvl_agg_stick(iogrp),depths)
  call msklvl(jlvl_agg_stickf(iogrp),depths)
  call msklvl(jlvl_agg_dmax(iogrp),depths)
  call msklvl(jlvl_agg_avdp(iogrp),depths)
  call msklvl(jlvl_agg_avrhop(iogrp),depths)
  call msklvl(jlvl_agg_avdC(iogrp),depths)
  call msklvl(jlvl_agg_df(iogrp),depths)
  call msklvl(jlvl_agg_b(iogrp),depths)
  call msklvl(jlvl_agg_Vrhof(iogrp),depths)
  call msklvl(jlvl_agg_Vpor(iogrp),depths)

  ! --- Compute log10 of pH
  if (LYR_PH(iogrp).ne.0) call loglyr(jph(iogrp),1.,0.)
  if (LVL_PH(iogrp).ne.0) call loglvl(jlvlph(iogrp),rnacc,0.)
#ifdef natDIC
  if (LYR_NATPH(iogrp).ne.0) call loglyr(jnatph(iogrp),1.,0.)
  if (LVL_NATPH(iogrp).ne.0) call loglvl(jlvlnatph(iogrp),rnacc,0.)
#endif

  ! --- Store 2d fields
  call wrtsrf(jkwco2(iogrp),SRF_KWCO2(iogrp),rnacc,0.,cmpflg,                   &
       &   'kwco2',' ',' ',' ')
  call wrtsrf(jpco2(iogrp),SRF_PCO2(iogrp),rnacc,0.,cmpflg,                     &
       &   'pco2','Surface PCO2',' ','uatm')
  call wrtsrf(jdmsflux(iogrp),SRF_DMSFLUX(iogrp),rnacc*1e3/dtbgc,0.,            &
       &   cmpflg,'dmsflux','DMS flux',' ','mol DMS m-2 s-1')
  call wrtsrf(jco2fxd(iogrp),SRF_CO2FXD(iogrp),rnacc*12./dtbgc,0.,              &
       &   cmpflg,'co2fxd','Downward CO2 flux',' ','kg C m-2 s-1')
  call wrtsrf(jco2fxu(iogrp),SRF_CO2FXU(iogrp),rnacc*12./dtbgc,0.,              &
       &   cmpflg,'co2fxu','Upward CO2 flux',' ','kg C m-2 s-1')
  call wrtsrf(joxflux(iogrp),SRF_OXFLUX(iogrp),rnacc*1e3/dtbgc,0.,              &
       &   cmpflg,'fgo2','Oxygen flux',' ','mol O2 m-2 s-1')
  call wrtsrf(jniflux(iogrp),SRF_NIFLUX(iogrp),rnacc*1e3/dtbgc,0.,              &
       &   cmpflg,'fgn2','Nitrogen flux',' ','mol N2 m-2 s-1')
  call wrtsrf(jdms(iogrp),SRF_DMS(iogrp),rnacc,0.,cmpflg,                       &
       &   'dms','DMS',' ','kmol DMS m-3')
  call wrtsrf(jdmsprod(iogrp),SRF_DMSPROD(iogrp),rnacc*1e3/dtbgc,0.,            &
       &   cmpflg,'dmsprod','DMS production from phytoplankton production',     &
       &   ' ','mol DMS m-2 s-1')
  call wrtsrf(jdms_bac(iogrp),SRF_DMS_BAC(iogrp),rnacc*1e3/dtbgc,0.,            &
       &   cmpflg,'dms_bac','DMS bacterial consumption',' ',                    &
       &   'mol DMS m-2 s-1')
  call wrtsrf(jdms_uv(iogrp),SRF_DMS_UV(iogrp),rnacc*1e3/dtbgc,0.,              &
       &   cmpflg,'dms_uv','DMS photolysis reduction',' ','mol DMS m-2 s-1')
  call wrtsrf(jexport(iogrp),SRF_EXPORT(iogrp),rnacc*1e3/dtbgc,0.,              &
       &   cmpflg,'epc100','Export production',' ','mol C m-2 s-1')
  call wrtsrf(jexposi(iogrp),SRF_EXPOSI(iogrp),rnacc*1e3/dtbgc,0.,              &
       &   cmpflg,'epsi100','Si export production',' ','mol Si m-2 s-1')
  call wrtsrf(jexpoca(iogrp),SRF_EXPOCA(iogrp),rnacc*1e3/dtbgc,0.,              &
       &   cmpflg,'epcalc100','Ca export production',' ','mol Ca m-2 s-1')
  call wrtsrf(jsrfdic(iogrp),SRF_DIC(iogrp),                                    &
       &   rnacc*1e3,0.,cmpflg,'srfdissic',                                     &
       &   'Surface dissolved inorganic carbon',' ','mol C m-3')
  call wrtsrf(jsrfalkali(iogrp),SRF_ALKALI(iogrp),                              &
       &   rnacc*1e3,0.,cmpflg,'srftalk',                                       &
       &   'Surface alkalinity',' ','eq m-3')
  call wrtsrf(jsrfphosph(iogrp),SRF_PHOSPH(iogrp),                              &
       &   rnacc*1e3,0.,cmpflg,'srfpo4',                                        &
       &   'Surface phosphorus',' ','mol P m-3')
  call wrtsrf(jsrfoxygen(iogrp),SRF_OXYGEN(iogrp),                              &
       &   rnacc*1e3,0.,cmpflg,'srfo2',                                         &
       &   'Surface oxygen',' ','mol O2 m-3')
  call wrtsrf(jsrfano3(iogrp),SRF_ANO3(iogrp),                                  &
       &   rnacc*1e3,0.,cmpflg,'srfno3',                                        &
       &   'Surface nitrate',' ','mol N m-3')
  call wrtsrf(jsrfsilica(iogrp),SRF_SILICA(iogrp),                              &
       &   rnacc*1e3,0.,cmpflg,'srfsi',                                         &
       &   'Surface silicate',' ','mol Si m-3')
  call wrtsrf(jsrfiron(iogrp),SRF_IRON(iogrp),                                  &
       &   rnacc*1e3,0.,cmpflg,'srfdfe',                                        &
       &   'Surface dissolved iron',' ','mol Fe m-3')
  call wrtsrf(jsrfphyto(iogrp),SRF_PHYTO(iogrp),                                &
       &   rnacc*1e3,0.,cmpflg,'srfphyc',                                       &
       &   'Surface phytoplankton',' ','mol P m-3')
  call wrtsrf(jintphosy(iogrp),INT_PHOSY(iogrp),                                &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'ppint',                                   &
       &   'Integrated primary production',' ','mol C m-2 s-1')
  call wrtsrf(jintnfix(iogrp),INT_NFIX(iogrp),                                  &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'nfixint',                                 &
       &   'Integrated nitrogen fixation',' ','mol N m-2 s-1')
  call wrtsrf(jintdnit(iogrp),INT_DNIT(iogrp),                                  &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'dnitint',                                 &
       &   'Integrated denitrification',' ','mol N m-2 s-1')
  call wrtsrf(jcarflx0100(iogrp),FLX_CAR0100(iogrp),                            &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'carflx0100',                              &
       &   'C flux at 100m',' ','mol C m-2 s-1')
  call wrtsrf(jcarflx0500(iogrp),FLX_CAR0500(iogrp),                            &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'carflx0500',                              &
       &   'C flux at 500m',' ','mol C m-2 s-1')
  call wrtsrf(jcarflx1000(iogrp),FLX_CAR1000(iogrp),                            &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'carflx1000',                              &
       &   'C flux at 1000m',' ','mol C m-2 s-1')
  call wrtsrf(jcarflx2000(iogrp),FLX_CAR2000(iogrp),                            &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'carflx2000',                              &
       &   'C flux at 2000m',' ','mol C m-2 s-1')
  call wrtsrf(jcarflx4000(iogrp),FLX_CAR4000(iogrp),                            &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'carflx4000',                              &
       &   'C flux at 4000m',' ','mol C m-2 s-1')
  call wrtsrf(jcarflx_bot(iogrp),FLX_CAR_BOT(iogrp),                            &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'carflx_bot',                              &
       &   'C flux to sediment',' ','mol C m-2 s-1')
  call wrtsrf(jbsiflx0100(iogrp),FLX_BSI0100(iogrp),                            &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'bsiflx0100',                              &
       &   'Opal flux at 100m',' ','mol Si m-2 s-1')
  call wrtsrf(jbsiflx0500(iogrp),FLX_BSI0500(iogrp),                            &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'bsiflx0500',                              &
       &   'Opal flux at 500m',' ','mol Si m-2 s-1')
  call wrtsrf(jbsiflx1000(iogrp),FLX_BSI1000(iogrp),                            &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'bsiflx1000',                              &
       &   'Opal flux at 1000m',' ','mol Si m-2 s-1')
  call wrtsrf(jbsiflx2000(iogrp),FLX_BSI2000(iogrp),                            &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'bsiflx2000',                              &
       &   'Opal flux at 2000m',' ','mol Si m-2 s-1')
  call wrtsrf(jbsiflx4000(iogrp),FLX_BSI4000(iogrp),                            &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'bsiflx4000',                              &
       &   'Opal flux at 4000m',' ','mol Si m-2 s-1')
  call wrtsrf(jbsiflx_bot(iogrp),FLX_BSI_BOT(iogrp),                            &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'bsiflx_bot',                              &
       &   'Opal flux to sediment',' ','mol Si m-2 s-1')
  call wrtsrf(jcalflx0100(iogrp),FLX_CAL0100(iogrp),                            &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'calflx0100',                              &
       &   'CaCO3 flux at 100m',' ','mol Ca m-2 s-1')
  call wrtsrf(jcalflx0500(iogrp),FLX_CAL0500(iogrp),                            &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'calflx0500',                              &
       &   'CaCO3 flux at 500m',' ','mol Ca m-2 s-1')
  call wrtsrf(jcalflx1000(iogrp),FLX_CAL1000(iogrp),                            &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'calflx1000',                              &
       &   'CaCO3 flux at 1000m',' ','mol Ca m-2 s-1')
  call wrtsrf(jcalflx2000(iogrp),FLX_CAL2000(iogrp),                            &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'calflx2000',                              &
       &   'CaCO3 flux at 2000m',' ','mol Ca m-2 s-1')
  call wrtsrf(jcalflx4000(iogrp),FLX_CAL4000(iogrp),                            &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'calflx4000',                              &
       &   'CaCO3 flux at 4000m',' ','mol Ca m-2 s-1')
  call wrtsrf(jcalflx_bot(iogrp),FLX_CAL_BOT(iogrp),                            &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'calflx_bot',                              &
       &   'CaCO3 flux to sediment',' ','mol Ca m-2 s-1')
#ifndef sedbypass
  call wrtsrf(jsediffic(iogrp),FLX_SEDIFFIC(iogrp),                             &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'sedfdic',' ',' ',' ')
  call wrtsrf(jsediffal(iogrp),FLX_SEDIFFAL(iogrp),                             &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'sedfalk',' ',' ',' ')
  call wrtsrf(jsediffph(iogrp),FLX_SEDIFFPH(iogrp),                             &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'sedfpho',' ',' ',' ')
  call wrtsrf(jsediffox(iogrp),FLX_SEDIFFOX(iogrp),                             &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'sedfox',' ',' ',' ')
  call wrtsrf(jsediffn2(iogrp),FLX_SEDIFFN2(iogrp),                             &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'sedfn2',' ',' ',' ')
  call wrtsrf(jsediffno3(iogrp),FLX_SEDIFFNO3(iogrp),                           &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'sedfno3',' ',' ',' ')
  call wrtsrf(jsediffsi(iogrp),FLX_SEDIFFSI(iogrp),                             &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'sedfsi',' ',' ',' ')
#endif
#if defined(extNcycle) && ! defined(sedbypass)
  call wrtsrf(jsediffnh4(iogrp),FLX_SEDIFFNH4(iogrp),                             &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'sedfnh4',' ',' ',' ')
  call wrtsrf(jsediffn2o(iogrp),FLX_SEDIFFN2O(iogrp),                             &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'sedfn2o',' ',' ',' ')
  call wrtsrf(jsediffno2(iogrp),FLX_SEDIFFNO2(iogrp),                             &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'sedfno2',' ',' ',' ')
#endif
  call wrtsrf(jn2ofx(iogrp),SRF_N2OFX(iogrp),rnacc*1e3/dtbgc,0.,                &
       &   cmpflg,'n2oflux','N2O flux',' ','mol N2O m-2 s-1')
#ifdef cisonew
  call wrtsrf(jco213fxd(iogrp),SRF_CO213FXD(iogrp),                             &
       &   rnacc*12./dtbgc,0.,cmpflg,'co213fxd',                                &
       &   'Downward 13CO2 flux',' ','kg C m-2 s-1')
  call wrtsrf(jco213fxu(iogrp),SRF_CO213FXU(iogrp),                             &
       &   rnacc*12./dtbgc,0.,cmpflg,'co213fxu',                                &
       &   'Upward 13CO2 flux',' ','kg C m-2 s-1')
  call wrtsrf(jco214fxd(iogrp),SRF_CO214FXD(iogrp),                             &
       &   rnacc*12.*c14fac/dtbgc,0.,cmpflg,'co214fxd',                         &
       &   'Downward 14CO2 flux',' ','kg C m-2 s-1')
  call wrtsrf(jco214fxu(iogrp),SRF_CO214FXU(iogrp),                             &
       &   rnacc*12.*c14fac/dtbgc,0.,cmpflg,'co214fxu',                         &
       &   'Upward 14CO2 flux',' ','kg C m-2 s-1')
#endif
#ifdef CFC
  call wrtsrf(jcfc11fx(iogrp),SRF_CFC11(iogrp),rnacc*1e3/dtbgc,0.,              &
       &   cmpflg,'cfc11flux','CFC-11 flux',' ','mol CFC12 m-2 s-1')
  call wrtsrf(jcfc12fx(iogrp),SRF_CFC12(iogrp),rnacc*1e3/dtbgc,0.,              &
       &   cmpflg,'cfc12flux','CFC-12 flux',' ','mol CFC12 m-2 s-1')
  call wrtsrf(jsf6fx(iogrp),SRF_SF6(iogrp),rnacc*1e3/dtbgc,0.,                  &
       &   cmpflg,'sf6flux','SF-6 flux',' ','mol SF6 m-2 s-1')
#endif
#ifdef natDIC
  call wrtsrf(jsrfnatdic(iogrp),SRF_NATDIC(iogrp),                              &
       &   rnacc*1e3,0.,cmpflg,'srfnatdissic',                                  &
       &   'Surface natural dissolved inorganic carbon',' ','mol C m-3')
  call wrtsrf(jsrfnatalk(iogrp),SRF_NATALKALI(iogrp),                           &
       &   rnacc*1e3,0.,cmpflg,'srfnattalk',                                    &
       &   'Surface natural alkalinity',' ','eq m-3')
  call wrtsrf(jnatpco2(iogrp),SRF_NATPCO2(iogrp),rnacc,0.,cmpflg,               &
       &   'natpco2','Surface natural PCO2',' ','uatm')
  call wrtsrf(jnatco2fx(iogrp),SRF_NATCO2FX(iogrp),                             &
       &   rnacc*12./dtbgc,0.,cmpflg,'natco2fx',                                &
       &   'Natural CO2 flux',' ','kg C m-2 s-1')
#endif
#ifdef BROMO
  call wrtsrf(jbromofx(iogrp),SRF_BROMOFX(iogrp),rnacc*1e3/dtbgc,               &
       &   0.,cmpflg,'bromofx','Bromoform flux',' ','mol CHBr3 m-2 s-1')
  call wrtsrf(jsrfbromo(iogrp),SRF_BROMO(iogrp),rnacc*1e3,0.,                   &
       &   cmpflg,'srfbromo','Surface bromoform',' ','mol CHBr3 m-3')
  call wrtsrf(jbromo_prod(iogrp),INT_BROMOPRO(iogrp),                           &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'intbromoprod',                            &
       &   'Integrated bromoform production',' ','mol CHBr3 m-2 s-1')
  call wrtsrf(jbromo_uv(iogrp),INT_BROMOUV(iogrp),                              &
       &   rnacc*1e3/dtbgc,0.,cmpflg,'intbromouv',                              &
       &   'Integrated bromoform loss to photolysis',' ',                       &
       &   'mol CHBr3 m-2 s-1')
  call wrtsrf(jatmbromo(iogrp),SRF_ATMBROMO(iogrp),rnacc,0.,                    &
       &   cmpflg,'atmbromo','Atmospheric bromoform',' ','ppt')
#endif


  call wrtsrf(jatmco2(iogrp),SRF_ATMCO2(iogrp),rnacc,0.,cmpflg,                 &
       &   'atmco2','Atmospheric CO2',' ','ppm')
#if defined(BOXATM)
  call wrtsrf(jatmo2(iogrp),SRF_ATMO2(iogrp),rnacc,0.,cmpflg,                   &
       &   'atmo2','Atmospheric O2',' ','ppm')
  call wrtsrf(jatmn2(iogrp),SRF_ATMN2(iogrp),rnacc,0.,cmpflg,                   &
       &   'atmn2','Atmospheric N2',' ','ppm')
#endif
#ifdef cisonew
  call wrtsrf(jatmc13(iogrp),SRF_ATMC13(iogrp),rnacc,0.,cmpflg,                 &
       &   'atmc13','Atmospheric 13CO2',' ','ppm')
  call wrtsrf(jatmc14(iogrp),SRF_ATMC14(iogrp),rnacc,0.,cmpflg,                 &
       &   'atmc14','Atmospheric 14CO2',' ','ppm')
#endif
#ifdef extNcycle
  call wrtsrf(jsrfanh4(iogrp),SRF_ANH4(iogrp),                                  &
       &  rnacc*1e3,0.,cmpflg,'srfnh4',                                         &
       &  'Surface ammonium',' ','mol N m-3')
  call wrtsrf(jsrfano2(iogrp),SRF_ANO2(iogrp),                                  &
       &  rnacc*1e3,0.,cmpflg,'srfno2',                                         &
       &  'Surface nitrite',' ','mol N m-3')
      call wrtsrf(janh3fx(iogrp),SRF_ANH3FX(iogrp),rnacc*1e3/dtbgc,0.,          &
       &  cmpflg,'nh3flux','NH3 flux',' ','mol NH3 m-2 s-1')
#endif

  ! --- Store 3d layer fields
  call wrtlyr(jdp(iogrp),LYR_DP(iogrp),rnacc,0.,cmpflg,                         &
       &   'pddpo','Layer thickness',' ','m')
  call wrtlyr(jdic(iogrp),LYR_DIC(iogrp),1e3,0.,cmpflg,                         &
       &   'dissic','Dissolved inorganic carbon',' ','mol C m-3')
  call wrtlyr(jalkali(iogrp),LYR_ALKALI(iogrp),1e3,0.,cmpflg,                   &
       &   'talk','Alkalinity',' ','eq m-3')
  call wrtlyr(jphosph(iogrp),LYR_PHOSPH(iogrp),1e3,0.,cmpflg,                   &
       &   'po4','Phosphorus',' ','mol P m-3')
  call wrtlyr(joxygen(iogrp),LYR_OXYGEN(iogrp),1e3,0.,cmpflg,                   &
       &   'o2','Oxygen',' ','mol O2 m-3')
  call wrtlyr(jano3(iogrp),LYR_ANO3(iogrp),1e3,0.,cmpflg,                       &
       &   'no3','Nitrate',' ','mol N m-3')
  call wrtlyr(jsilica(iogrp),LYR_SILICA(iogrp),1e3,0.,cmpflg,                   &
       &   'si','Silicate',' ','mol Si m-3')
  call wrtlyr(jdoc(iogrp),LYR_DOC(iogrp),1e3,0.,cmpflg,                         &
       &   'dissoc','Dissolved organic carbon',' ','mol P m-3')
  call wrtlyr(jphyto(iogrp),LYR_PHYTO(iogrp),1e3,0.,cmpflg,                     &
       &   'phyc','Phytoplankton',' ','mol P m-3')
  call wrtlyr(jgrazer(iogrp),LYR_GRAZER(iogrp),1e3,0.,cmpflg,                   &
       &   'zooc','Zooplankton',' ','mol P m-3')
  call wrtlyr(jpoc(iogrp),LYR_POC(iogrp),1e3,0.,cmpflg,                         &
       &   'detoc','Detritus',' ','mol P m-3')
  call wrtlyr(jcalc(iogrp),LYR_CALC(iogrp),1e3,0.,cmpflg,                       &
       &   'calc','CaCO3 shells',' ','mol C m-3')
  call wrtlyr(jopal(iogrp),LYR_OPAL(iogrp),1e3,0.,cmpflg,                       &
       &   'opal','Opal shells',' ','mol Si m-3')
  call wrtlyr(jiron(iogrp),LYR_IRON(iogrp),1e3,0.,cmpflg,                       &
       &   'dfe','Dissolved iron',' ','mol Fe m-3')
  call wrtlyr(jphosy(iogrp),LYR_PHOSY(iogrp),1e3/dtbgc,0.,cmpflg,               &
       &   'pp','Primary production',' ','mol C m-3 s-1')
  call wrtlyr(jco3(iogrp),LYR_CO3(iogrp),1e3,0.,cmpflg,                         &
       &   'co3','Carbonate ions',' ','mol C m-3')
  call wrtlyr(jph(iogrp),LYR_PH(iogrp),-1.,0.,cmpflg,                           &
       &   'ph','pH',' ','-log10([h+])')
  call wrtlyr(jomegaa(iogrp),LYR_OMEGAA(iogrp),1.,0.,cmpflg,                    &
       &   'omegaa','OmegaA',' ','-')
  call wrtlyr(jomegac(iogrp),LYR_OMEGAC(iogrp),1.,0.,cmpflg,                    &
       &   'omegac','OmegaC',' ','-')
  call wrtlyr(jn2o(iogrp),LYR_N2O(iogrp),1e3,0.,cmpflg,                         &
       &   'n2o','N2O',' ','mol N2O m-3')
  call wrtlyr(jprefo2(iogrp),LYR_PREFO2(iogrp),1e3,0.,cmpflg,                   &
       &   'p_o2','Preformed oxygen',' ','mol O2 m-3')
  call wrtlyr(jo2sat(iogrp),LYR_O2SAT(iogrp),1e3,0.,cmpflg,                     &
       &   'satoxy','Saturated oxygen',' ','mol O2 m-3')
  call wrtlyr(jprefpo4(iogrp),LYR_PREFPO4(iogrp),1e3,0.,cmpflg,                 &
       &   'p_po4','Preformed phosphorus',' ','mol P m-3')
  call wrtlyr(jprefalk(iogrp),LYR_PREFALK(iogrp),1e3,0.,cmpflg,                 &
       &   'p_talk','Preformed alkalinity',' ','eq m-3')
  call wrtlyr(jprefdic(iogrp),LYR_PREFDIC(iogrp),1e3,0.,cmpflg,                 &
       &   'p_dic','Preformed DIC',' ','mol C m-3')
  call wrtlyr(jdicsat(iogrp),LYR_DICSAT(iogrp),1e3,0.,cmpflg,                   &
       &   'sat_dic','Saturated DIC',' ','mol C m-3')
#ifdef cisonew
  call wrtlyr(jdic13(iogrp),LYR_DIC13(iogrp),1.e3,0.,cmpflg,                    &
       &   'dissic13','Dissolved C13',' ','mol 13C m-3')
  call wrtlyr(jdic14(iogrp),LYR_DIC14(iogrp),1.e3*c14fac,0.,cmpflg,             &
       &   'dissic14','Dissolved C14',' ','mol 14C m-3')
  call wrtlyr(jd13c(iogrp),LYR_D13C(iogrp),1.,0.,cmpflg,                        &
       &   'delta13c','delta13C of DIC',' ','permil')
  call wrtlyr(jd14c(iogrp),LYR_D14C(iogrp),1.,0.,cmpflg,                        &
       &   'delta14c','delta14C of DIC',' ','permil')
  call wrtlyr(jbigd14c(iogrp),LYR_BIGD14C(iogrp),1.,0.,cmpflg,                  &
       &   'bigdelta14c','big delta14C of DIC',' ','permil')
  call wrtlyr(jpoc13(iogrp),LYR_POC13(iogrp),1e3,0.,cmpflg,                     &
       &   'detoc13','Detritus13',' ','mol P m-3')
  call wrtlyr(jdoc13(iogrp),LYR_DOC13(iogrp),1e3,0.,cmpflg,                     &
       &   'dissoc13','Dissolved organic carbon13',' ','mol P m-3')
  call wrtlyr(jcalc13(iogrp),LYR_CALC13(iogrp),1e3,0.,cmpflg,                   &
       &   'calc13','Ca13CO3 shells',' ','mol 13C m-3')
  call wrtlyr(jphyto13(iogrp),LYR_PHYTO13(iogrp),1e3,0.,cmpflg,                 &
       &   'phyc13','Phytoplankton13',' ','mol P m-3')
  call wrtlyr(jgrazer13(iogrp),LYR_GRAZER13(iogrp),1e3,0.,cmpflg,               &
       &   'zooc13','Zooplankton13',' ','mol P m-3')
#endif
#ifdef AGG
  call wrtlyr(jnos(iogrp),LYR_NOS(iogrp),1.,0.,cmpflg,                          &
       &   'nos','Marine snow aggregates per cm^3 sea water',' ','1/cm^3')
  call wrtlyr(jwphy(iogrp),LYR_WPHY(iogrp),1.,0.,cmpflg,                        &
       &   'wphy','Av. mass sinking speed of marine snow',' ','m/day')
  call wrtlyr(jwnos(iogrp),LYR_WNOS(iogrp),1.,0.,cmpflg,                        &
       &   'wnos','Av. number sinking speed of marine snow',' ','m/day')
  call wrtlyr(jeps(iogrp),LYR_EPS(iogrp),1.,0.,cmpflg,                          &
       &   'eps','Av. size distribution exponent',' ','-')
  call wrtlyr(jasize(iogrp),LYR_ASIZE(iogrp),1.,0.,cmpflg,                      &
       &   'asize','Av. size of marine snow aggregates',' ','nb. of cells')
#endif
#ifdef CFC
  call wrtlyr(jcfc11(iogrp),LYR_CFC11(iogrp),1e3,0.,cmpflg,                     &
       &   'cfc11','CFC-11',' ','mol cfc11 m-3')
  call wrtlyr(jcfc12(iogrp),LYR_CFC12(iogrp),1e3,0.,cmpflg,                     &
       &   'cfc12','CFC-12',' ','mol cfc12 m-3')
  call wrtlyr(jsf6(iogrp),LYR_SF6(iogrp),1e3,0.,cmpflg,                         &
       &   'sf6','SF-6',' ','mol sf6 m-3')
#endif
#ifdef natDIC
  call wrtlyr(jnatco3(iogrp),LYR_NATCO3(iogrp),1e3,0.,cmpflg,                   &
       &   'natco3','Natural Carbonate ions',' ','mol C m-3')
  call wrtlyr(jnatalkali(iogrp),LYR_NATALKALI(iogrp),1e3,0.,cmpflg,             &
       &   'nattalk','Natural alkalinity',' ','eq m-3')
  call wrtlyr(jnatdic(iogrp),LYR_NATDIC(iogrp),1e3,0.,cmpflg,                   &
       &   'natdissic','Natural dissolved inorganic carbon',' ',                &
       &   'mol C m-3')
  call wrtlyr(jnatcalc(iogrp),LYR_NATCALC(iogrp),1e3,0.,cmpflg,                 &
       &   'natcalc','Natural CaCO3 shells',' ','mol C m-3')
  call wrtlyr(jnatph(iogrp),LYR_NATPH(iogrp),-1.,0.,cmpflg,                     &
       &   'natph','Natural pH',' ','-log10([h+])')
  call wrtlyr(jnatomegaa(iogrp),LYR_NATOMEGAA(iogrp),1.,0.,cmpflg,              &
       &   'natomegaa','Natural OmegaA',' ','-')
  call wrtlyr(jnatomegac(iogrp),LYR_NATOMEGAC(iogrp),1.,0.,cmpflg,              &
       &   'natomegac','Natural OmegaC',' ','-')
#endif
#ifdef BROMO
  call wrtlyr(jbromo(iogrp),LYR_BROMO(iogrp),1e3,0.,cmpflg,                     &
       &   'bromo','Bromoform',' ','mol CHBr3 m-3')
#endif
#ifdef extNcycle
  call wrtlyr(janh4(iogrp),LYR_ANH4(iogrp),1e3,0.,cmpflg,                       &         
     &  'nh4','Ammonium',' ','mol N m-3')
  call wrtlyr(jano2(iogrp),LYR_ANO2(iogrp),1e3,0.,cmpflg,                       &
     &  'no2','Nitrite',' ','mol N m-3')
  call wrtlyr(jnitr_NH4(iogrp),LYR_nitr_NH4(iogrp),1e3/dtbgc,0.,                &
     &   cmpflg,                                                                &
     &  'nh4nitr','NH4 nitrification rate',' ','mol N m-3 s-1')
  call wrtlyr(jnitr_NO2(iogrp),LYR_nitr_NO2(iogrp),1e3/dtbgc,0.,                &
     &   cmpflg,                                                                &
     &  'no2nitr','NO2 nitrification rate',' ','mol N m-3 s-1')
  call wrtlyr(jnitr_N2O_prod(iogrp),LYR_nitr_N2O_prod(iogrp),                   &
     &  1e3/dtbgc,0.,cmpflg,                                                    &
     &  'nitr_n2o','N2O prod during NH4 nitrification',' ',                     &
     &  'mol N2O m-3 s-1')
  call wrtlyr(jnitr_NH4_OM(iogrp),LYR_nitr_NH4_OM(iogrp),1e3/dtbgc,             &
     &  0.,cmpflg,                                                              &
     &  'nh4nitr_om','OM production during NH4 nitrification',' ',              &
     &  'mol P m-3 s-1')
  call wrtlyr(jnitr_NO2_OM(iogrp),LYR_nitr_NO2_OM(iogrp),1e3/dtbgc,             &
     &  0.,cmpflg,                                                              &
     &  'no2nitr_om','OM production during NO2 nitrification',' ',              &
     &  'mol P m-3 s-1')
  call wrtlyr(jdenit_NO3(iogrp),LYR_denit_NO3(iogrp),1e3/dtbgc,0.,              &
     &  cmpflg,                                                                 &
     &  'no3denit','NO3 denitrification rate',' ','mol N m-3 s-1') 
  call wrtlyr(jdenit_NO2(iogrp),LYR_denit_NO2(iogrp),1e3/dtbgc,0.,              &
     &  cmpflg,                                                                 &
     &  'no2denit','NO2 denitrification rate',' ','mol N m-3 s-1')
  call wrtlyr(jdenit_N2O(iogrp),LYR_denit_N2O(iogrp),1e3/dtbgc,0.,              &
     &  cmpflg,                                                                 &
     &  'n2odenit','N2O denitrification rate',' ','mol N2O m-3 s-1')
  call wrtlyr(jDNRA_NO2(iogrp),LYR_DNRA_NO2(iogrp),1e3/dtbgc,0.,                &
     &  cmpflg,                                                                 &
     &  'no2dnra','NO2 DNRA rate',' ','mol N m-3 s-1')
  call wrtlyr(janmx_N2_prod(iogrp),LYR_anmx_N2_prod(iogrp),                     &
     &  1e3/dtbgc,0.,cmpflg,                                                    &
     &  'anmx_n2','Anammox N2 production rate',' ','mol N2 m-3 s-1')
  call wrtlyr(janmx_OM_prod(iogrp),LYR_anmx_OM_prod(iogrp),                     &
     &  1e3/dtbgc,0.,cmpflg,                                                    &
     &  'anmx_om','Anammox OM production rate',' ','mol P m-3 s-1')
  call wrtlyr(jphosy_NH4(iogrp),LYR_phosy_NH4(iogrp),1e3/dtbgc,0.,              &
     &  cmpflg,                                                                 &
     &  'phosy_nh4','PP consumption rate of NH4',' ','mol N m-3 s-1')
  call wrtlyr(jphosy_NO3(iogrp),LYR_phosy_NO3(iogrp),1e3/dtbgc,0.,              &
     &  cmpflg,                                                                 &
     &  'phosy_no3','PP consumption rate of NO3',' ','mol N m-3 s-1')
  call wrtlyr(jremin_aerob(iogrp),LYR_remin_aerob(iogrp),1e3/dtbgc,             &
     &  0.,cmpflg,                                                              &
     &  'remina','Aerob remineralization rate',' ','mol N m-3 s-1')
  call wrtlyr(jremin_sulf(iogrp),LYR_remin_sulf(iogrp),1e3/dtbgc,               &
     &  0.,cmpflg,                                                              &
     &  'remins','Sulfate remineralization rate',' ','mol P m-3 s-1')
#endif
!      M4AGO
  call wrtlyr(jagg_ws(iogrp),LYR_agg_ws(iogrp),1.,0.,cmpflg,                    &         
     &  'agg_ws','aggregate mean settling velocity',' ','m d-1')
  call wrtlyr(jdynvis(iogrp),LYR_dynvis(iogrp),1.,0.,cmpflg,                    &
     &  'dynvis','dynamic viscosity of sea water',' ','kg m-1 s-1')
  call wrtlyr(jagg_stick(iogrp),LYR_agg_stick(iogrp),1.,0.,cmpflg,              &
     &  'agg_stick','aggregate mean stickiness',' ','-')
  call wrtlyr(jagg_stickf(iogrp),LYR_agg_stickf(iogrp),1.,0.,cmpflg,            &
     &  'agg_stickf','opal frustule stickiness',' ','-')
  call wrtlyr(jagg_dmax(iogrp),LYR_agg_dmax(iogrp),1.,0.,cmpflg,                &
     &  'agg_dmax','aggregate maximum diameter',' ','m')
  call wrtlyr(jagg_avdp(iogrp),LYR_agg_avdp(iogrp),1.,0.,cmpflg,                &
     &  'agg_avdp','mean primary particle diameter',' ','m')
  call wrtlyr(jagg_avrhop(iogrp),LYR_agg_avrhop(iogrp),1.,0.,cmpflg,            &
     &  'agg_avrhop','mean primary particle density',' ','kg m-3')
  call wrtlyr(jagg_avdC(iogrp),LYR_agg_avdC(iogrp),1.,0.,cmpflg,                &
     &  'agg_avdC','Conc.-weighted mean aggregate diameter',' ','m')
  call wrtlyr(jagg_df(iogrp),LYR_agg_df(iogrp),1.,0.,cmpflg,                    &
     &  'agg_df','aggregate fractal dimension',' ','-')
  call wrtlyr(jagg_b(iogrp),LYR_agg_b(iogrp),1.,0.,cmpflg,                      &
     &  'agg_b','aggregate number distribution slope',' ','-')
  call wrtlyr(jagg_Vrhof(iogrp),LYR_agg_Vrhof(iogrp),1.,0.,cmpflg,              &
     &  'agg_Vrhof','V-weighted aggregate mean density',' ','kg m-3')
  call wrtlyr(jagg_Vpor(iogrp),LYR_agg_Vpor(iogrp),1.,0.,cmpflg,                &
     &  'agg_Vpor','V-weighted aggregate mean porosity',' ','-')

  ! --- Store 3d level fields
  call wrtlvl(jlvldic(iogrp),LVL_DIC(iogrp),rnacc*1e3,0.,cmpflg,                &
       &   'dissiclvl','Dissolved inorganic carbon',' ','mol C m-3')
  call wrtlvl(jlvlalkali(iogrp),LVL_ALKALI(iogrp),rnacc*1e3,0.,                 &
       &   cmpflg, 'talklvl','Alkalinity',' ','eq m-3')
  call wrtlvl(jlvlphosph(iogrp),LVL_PHOSPH(iogrp),rnacc*1e3,0.,                 &
       &   cmpflg,'po4lvl','Phosphorus',' ','mol P m-3')
  call wrtlvl(jlvloxygen(iogrp),LVL_OXYGEN(iogrp),rnacc*1e3,0.,                 &
       &   cmpflg,'o2lvl','Oxygen',' ','mol O2 m-3')
  call wrtlvl(jlvlano3(iogrp),LVL_ANO3(iogrp),rnacc*1e3,0.,cmpflg,              &
       &   'no3lvl','Nitrate',' ','mol N m-3')
  call wrtlvl(jlvlsilica(iogrp),LVL_SILICA(iogrp),rnacc*1e3,0.,                 &
       &   cmpflg, 'silvl','Silicate',' ','mol Si m-3')
  call wrtlvl(jlvldoc(iogrp),LVL_DOC(iogrp),rnacc*1e3,0.,cmpflg,                &
       &   'dissoclvl','Dissolved organic carbon',' ','mol P m-3')
  call wrtlvl(jlvlphyto(iogrp),LVL_PHYTO(iogrp),rnacc*1e3,0.,cmpflg,            &
       &   'phyclvl','Phytoplankton',' ','mol P m-3')
  call wrtlvl(jlvlgrazer(iogrp),LVL_GRAZER(iogrp),rnacc*1e3,0.,                 &
       &   cmpflg,'zooclvl','Zooplankton',' ','mol P m-3')
  call wrtlvl(jlvlpoc(iogrp),LVL_POC(iogrp),rnacc*1e3,0.,cmpflg,                &
       &   'detoclvl','Detritus',' ','mol P m-3')
  call wrtlvl(jlvlcalc(iogrp),LVL_CALC(iogrp),rnacc*1e3,0.,cmpflg,              &
       &   'calclvl','CaCO3 shells',' ','mol C m-3')
  call wrtlvl(jlvlopal(iogrp),LVL_OPAL(iogrp),rnacc*1e3,0.,cmpflg,              &
       &   'opallvl','Opal shells',' ','mol Si m-3')
  call wrtlvl(jlvliron(iogrp),LVL_IRON(iogrp),rnacc*1e3,0.,cmpflg,              &
       &   'dfelvl','Dissolved iron',' ','mol Fe m-3')
  call wrtlvl(jlvlphosy(iogrp),LVL_PHOSY(iogrp),rnacc*1e3/dtbgc,0.,             &
       &   cmpflg,'pplvl','Primary production',' ','mol C m-3 s-1')
  call wrtlvl(jlvlco3(iogrp),LVL_CO3(iogrp),rnacc*1e3,0.,cmpflg,                &
       &   'co3lvl','Carbonate ions',' ','mol C m-3')
  call wrtlvl(jlvlph(iogrp),LVL_PH(iogrp),-1.,0.,cmpflg,                        &
       &   'phlvl','pH',' ','-log10([h+])')
  call wrtlvl(jlvlomegaa(iogrp),LVL_OMEGAA(iogrp),rnacc,0.,cmpflg,              &
       &   'omegaalvl','OmegaA',' ','-')
  call wrtlvl(jlvlomegac(iogrp),LVL_OMEGAC(iogrp),rnacc,0.,cmpflg,              &
       &   'omegaclvl','OmegaC',' ','-')
  call wrtlvl(jlvln2o(iogrp),LVL_N2O(iogrp),rnacc*1e3,0.,cmpflg,                &
       &   'n2olvl','N2O',' ','mol N2O m-3')
  call wrtlvl(jlvlprefo2(iogrp),LVL_PREFO2(iogrp),rnacc*1e3,0.,                 &
       &   cmpflg,'p_o2lvl','Preformed oxygen',' ','mol O2 m-3')
  call wrtlvl(jlvlo2sat(iogrp),LVL_O2SAT(iogrp),rnacc*1e3,0.,                   &
       &   cmpflg,'satoxylvl','Saturated oxygen',' ','mol O2 m-3')
  call wrtlvl(jlvlprefpo4(iogrp),LVL_PREFPO4(iogrp),rnacc*1e3,0.,               &
       &   cmpflg,'p_po4lvl','Preformed phosphorus',' ','mol P m-3')
  call wrtlvl(jlvlprefalk(iogrp),LVL_PREFALK(iogrp),rnacc*1e3,0.,               &
       &   cmpflg, 'p_talklvl','Preformed alkalinity',' ','eq m-3')
  call wrtlvl(jlvlprefdic(iogrp),LVL_PREFDIC(iogrp),rnacc*1e3,0.,               &
       &   cmpflg, 'p_diclvl','Preformed DIC',' ','mol C m-3')
  call wrtlvl(jlvldicsat(iogrp),LVL_DICSAT(iogrp),rnacc*1e3,0.,                 &
       &   cmpflg, 'sat_diclvl','Saturated DIC',' ','mol C m-3')
#ifdef cisonew
  call wrtlvl(jlvldic13(iogrp),LVL_DIC13(iogrp),rnacc*1.e3,                     &
       &   0.,cmpflg,'dissic13lvl','Dissolved C13',' ','mol 13C m-3')
  call wrtlvl(jlvldic14(iogrp),LVL_DIC14(iogrp),rnacc*1.e3*c14fac,              &
       &   0.,cmpflg,'dissic14lvl','Dissolved C14',' ','mol 14C m-3')
  call wrtlvl(jlvld13c(iogrp),LVL_D13C(iogrp),rnacc,                            &
       &   0.,cmpflg,'delta13clvl','delta13C of DIC',' ','permil')
  call wrtlvl(jlvld14c(iogrp),LVL_D14C(iogrp),rnacc,                            &
       &   0.,cmpflg,'delta14clvl','delta14C of DIC',' ','permil')
  call wrtlvl(jlvlbigd14c(iogrp),LVL_BIGD14C(iogrp),rnacc,                      &
       &   0.,cmpflg,'bigdelta14clvl','big delta14C of DIC',' ','permil')
  call wrtlvl(jlvlpoc13(iogrp),LVL_POC13(iogrp),rnacc*1e3,                      &
       &   0.,cmpflg,'detoc13lvl','Detritus13',' ','mol P m-3')
  call wrtlvl(jlvldoc13(iogrp),LVL_DOC13(iogrp),rnacc*1e3,                      &
       &   0.,cmpflg,'dissoc13lvl','Dissolved organic carbon13',' ',            &
       &   'mol P m-3')
  call wrtlvl(jlvlcalc13(iogrp),LVL_CALC13(iogrp),rnacc*1e3,                    &
       &   0.,cmpflg,'calc13lvl','Ca13CO3 shells',' ','mol 13C m-3')
  call wrtlvl(jlvlphyto13(iogrp),LVL_PHYTO13(iogrp),rnacc*1e3,                  &
       &   0.,cmpflg,'phyc13lvl','Phytoplankton13',' ','mol P m-3')
  call wrtlvl(jlvlgrazer13(iogrp),LVL_GRAZER13(iogrp),rnacc*1e3,                &
       &   0.,cmpflg,'zooc13lvl','Zooplankton13',' ','mol P m-3')
#endif
#ifdef AGG
  call wrtlvl(jlvlnos(iogrp),LVL_NOS(iogrp),                                    &
       &   rnacc,0.,cmpflg,'noslvl',                                            &
       &   'Marine snow aggregates per cm^3 sea water',' ','1/cm^3')
  call wrtlvl(jlvlwphy(iogrp),LVL_WPHY(iogrp),                                  &
       &   rnacc,0.,cmpflg,'wphylvl',                                           &
       &   'Av. mass sinking speed of marine snow',' ','m/day')
  call wrtlvl(jlvlwnos(iogrp),LVL_WNOS(iogrp),                                  &
       &   rnacc,0.,cmpflg,'wnoslvl',                                           &
       &   'Av. number sinking speed of marine snow',' ','m/day')
  call wrtlvl(jlvleps(iogrp),LVL_EPS(iogrp),                                    &
       &   rnacc,0.,cmpflg,'epslvl',                                            &
       &   'Av. size distribution exponent',' ','-')
  call wrtlvl(jlvlasize(iogrp),LVL_ASIZE(iogrp),                                &
       &   rnacc,0.,cmpflg,'asizelvl',                                          &
       &   'Av. size of marine snow aggregates',' ','nb. of cells')
#endif
#ifdef CFC
  call wrtlvl(jlvlcfc11(iogrp),LVL_CFC11(iogrp),rnacc*1e3,0.,cmpflg,            &
       &   'cfc11lvl','CFC-11',' ','mol cfc11 m-3')
  call wrtlvl(jlvlcfc12(iogrp),LVL_CFC12(iogrp),rnacc*1e3,0.,cmpflg,            &
       &   'cfc12lvl','CFC-12',' ','mol cfc12 m-3')
  call wrtlvl(jlvlsf6(iogrp),LVL_SF6(iogrp),rnacc*1e3,0.,cmpflg,                &
       &   'sf6lvl','SF-6',' ','mol sf6 m-3')
#endif
#ifdef natDIC
  call wrtlvl(jlvlnatco3(iogrp),LVL_NATCO3(iogrp),                              &
       &   rnacc*1e3,0.,cmpflg,'natco3lvl',                                     &
       &   'Natural carbonate ions',' ','mol C m-3')
  call wrtlvl(jlvlnatalkali(iogrp),LVL_NATALKALI(iogrp),                        &
       &   rnacc*1e3,0.,cmpflg,'nattalklvl',                                    &
       &   'Natural alkalinity',' ','eq m-3')
  call wrtlvl(jlvlnatdic(iogrp),LVL_NATDIC(iogrp),                              &
       &   rnacc*1e3,0.,cmpflg,'natdissiclvl',                                  &
       &   'Natural dissolved inorganic carbon',' ','mol C m-3')
  call wrtlvl(jlvlnatcalc(iogrp),LVL_NATCALC(iogrp),                            &
       &   rnacc*1e3,0.,cmpflg,'natcalclvl',                                    &
       &   'Natural CaCO3 shells',' ','mol C m-3')
  call wrtlvl(jlvlnatph(iogrp),LVL_NATPH(iogrp),-1.,0.,cmpflg,                  &
       &   'natphlvl','Natural pH',' ','-log10([h+])')
  call wrtlvl(jlvlnatomegaa(iogrp),LVL_NATOMEGAA(iogrp),                        &
       &   rnacc,0.,cmpflg,'natomegaalvl',                                      &
       &   'Natural OmegaA',' ','-')
  call wrtlvl(jlvlnatomegac(iogrp),LVL_NATOMEGAC(iogrp),                        &
       &   rnacc,0.,cmpflg,'natomegaclvl',                                      &
       &   'Natural OmegaC',' ','-')
#endif
#ifdef BROMO
  call wrtlvl(jlvlbromo(iogrp),LVL_BROMO(iogrp),rnacc*1e3,0.,                   &
       &   cmpflg,'bromolvl','Bromoform',' ','mol CHBr3 m-3')
#endif
#ifdef extNcycle
  call wrtlvl(jlvlanh4(iogrp),LVL_ANH4(iogrp),rnacc*1e3,0.,cmpflg,              &
     &  'nh4lvl','Ammonium',' ','mol N m-3')
  call wrtlvl(jlvlano2(iogrp),LVL_ANO2(iogrp),rnacc*1e3,0.,cmpflg,              &
     &  'no2lvl','Nitrite',' ','mol N m-3')
  call wrtlvl(jlvl_nitr_NH4(iogrp),LVL_nitr_NH4(iogrp),                         &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'nh4nitrlvl','NH4 nitrification rate',' ','mol N m-3 s-1')
  call wrtlvl(jlvl_nitr_NO2(iogrp),LVL_nitr_NO2(iogrp),                         &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'no2nitrlvl','NO2 nitrification rate',' ','mol N m-3 s-1')
  call wrtlvl(jlvl_nitr_N2O_prod(iogrp),LVL_nitr_N2O_prod(iogrp),               &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'nitr_n2olvl','N2O prod during NH4 nitrification',' ',                  &
     &  'mol N2O m-3 s-1')
  call wrtlvl(jlvl_nitr_NH4_OM(iogrp),LVL_nitr_NH4_OM(iogrp),                   &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'nh4nitr_omlvl','OM production during NH4 nitrification',' ',           &
     &  'mol P m-3 s-1')
  call wrtlvl(jlvl_nitr_NO2_OM(iogrp),LVL_nitr_NO2_OM(iogrp),                   &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'no2nitr_omlvl','OM production during NO2 nitrification',' ',           &
     &  'mol P m-3 s-1')
  call wrtlvl(jlvl_denit_NO3(iogrp),LVL_denit_NO3(iogrp),                       &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'no3denitlvl','NO3 denitrification rate',' ','mol N m-3 s-1')
  call wrtlvl(jlvl_denit_NO2(iogrp),LVL_denit_NO2(iogrp),                       &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'no2denitlvl','NO2 denitrification rate',' ','mol N m-3 s-1')
  call wrtlvl(jlvl_denit_N2O(iogrp),LVL_denit_N2O(iogrp),                       &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'n2odenitlvl','N2O denitrification rate',' ','mol N2O m-3 s-1')
  call wrtlvl(jlvl_DNRA_NO2(iogrp),LVL_DNRA_NO2(iogrp),                         &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'no2dnralvl','NO2 DNRA rate',' ','mol N m-3 s-1')
  call wrtlvl(jlvl_anmx_N2_prod(iogrp),LVL_anmx_N2_prod(iogrp),                 &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'anmx_n2lvl','Anammox N2 production rate',' ','mol N2 m-3 s-1')
  call wrtlvl(jlvl_anmx_OM_prod(iogrp),LVL_anmx_OM_prod(iogrp),                 &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'anmx_omlvl','Anammox OM production rate',' ','mol P m-3 s-1')
  call wrtlvl(jlvl_phosy_NH4(iogrp),LVL_phosy_NH4(iogrp),                       &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'phosy_nh4lvl','PP consumption rate of NH4',' ',                        &
     &  'mol N m-3 s-1')
  call wrtlvl(jlvl_phosy_NO3(iogrp),LVL_phosy_NO3(iogrp),                       &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'phosy_no3lvl','PP consumption rate of NO3',' ',                        &
     &  'mol N m-3 s-1')
  call wrtlvl(jlvl_remin_aerob(iogrp),LVL_remin_aerob(iogrp),                   &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'reminalvl','Aerob remineralization rate',' ','mol N m-3 s-1')
  call wrtlvl(jlvl_remin_sulf(iogrp),LVL_remin_sulf(iogrp),                     &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'reminslvl','Sulfate remineralization rate',' ','mol P m-3 s-1')
#endif
!      M4AGO
  call wrtlvl(jlvl_agg_ws(iogrp),LVL_agg_ws(iogrp),rnacc,0.,cmpflg,             &
     &  'agg_wslvl','aggregate mean settling velocity',' ','m d-1')
  call wrtlvl(jlvl_dynvis(iogrp),LVL_dynvis(iogrp),rnacc,0.,cmpflg,             &
     &  'dynvislvl','dynamic viscosity of sea water',' ','kg m-1 s-1')
  call wrtlvl(jlvl_agg_stick(iogrp),LVL_agg_stick(iogrp),rnacc,0.,              &
     &   cmpflg,                                                                &
     &  'agg_sticklvl','aggregate mean stickiness',' ','-')
  call wrtlvl(jlvl_agg_stickf(iogrp),LVL_agg_stickf(iogrp),rnacc,0.,            &
     &   cmpflg,                                                                &
     &  'agg_stickflvl','opal frustule stickiness',' ','-')
  call wrtlvl(jlvl_agg_dmax(iogrp),LVL_agg_dmax(iogrp),rnacc,0.,                &
     &   cmpflg,                                                                &
     &  'agg_dmaxlvl','aggregate maximum diameter',' ','m')
  call wrtlvl(jlvl_agg_avdp(iogrp),LVL_agg_avdp(iogrp),rnacc,0.,                &
     &   cmpflg,                                                                &
     &  'agg_avdplvl','mean primary particle diameter',' ','m')
  call wrtlvl(jlvl_agg_avrhop(iogrp),LVL_agg_avrhop(iogrp),rnacc,0.,            &
     &   cmpflg,                                                                &
     &  'agg_avrhoplvl','mean primary particle density',' ','kg m-3')
  call wrtlvl(jlvl_agg_avdC(iogrp),LVL_agg_avdC(iogrp),rnacc,0.,                &
     &   cmpflg,                                                                &
     &  'agg_avdClvl','Conc.-weighted mean aggregate diameter',' ','m')
  call wrtlvl(jlvl_agg_df(iogrp),LVL_agg_df(iogrp),rnacc,0.,cmpflg,             &
     &  'agg_dflvl','aggregate fractal dimension',' ','-')
  call wrtlvl(jlvl_agg_b(iogrp),LVL_agg_b(iogrp),rnacc,0.,cmpflg,               &
     &  'agg_blvl','aggregate number distribution slope',' ','-')
  call wrtlvl(jlvl_agg_Vrhof(iogrp),LVL_agg_Vrhof(iogrp),rnacc,0.,              &
     &   cmpflg,                                                                &
     &  'agg_Vrhoflvl','V-weighted aggregate mean density',' ','kg m-3')
  call wrtlvl(jlvl_agg_Vpor(iogrp),LVL_agg_Vpor(iogrp),rnacc,0.,                &
     &   cmpflg,                                                                &
     &  'agg_Vporlvl','V-weighted aggregate mean porosity',' ','-')

  ! --- Store sediment fields
#ifndef sedbypass
  call wrtsdm(jpowaic(iogrp),SDM_POWAIC(iogrp),rnacc*1e3,0.,cmpflg,             &
       &   'powdic','PoWa DIC',' ','mol C m-3')
  call wrtsdm(jpowaal(iogrp),SDM_POWAAL(iogrp),rnacc*1e3,0.,cmpflg,             &
       &   'powalk','PoWa alkalinity',' ','eq m-3')
  call wrtsdm(jpowaph(iogrp),SDM_POWAPH(iogrp),rnacc*1e3,0.,cmpflg,             &
       &   'powpho','PoWa phosphorus',' ','mol P m-3')
  call wrtsdm(jpowaox(iogrp),SDM_POWAOX(iogrp),rnacc*1e3,0.,cmpflg,             &
       &   'powox','PoWa oxygen',' ','mol O2 m-3')
  call wrtsdm(jpown2(iogrp),SDM_POWN2(iogrp),  rnacc*1e3,0.,cmpflg,             &
       &   'pown2','PoWa N2',' ','mol N2 m-3')
  call wrtsdm(jpowno3(iogrp),SDM_POWNO3(iogrp),rnacc*1e3,0.,cmpflg,             &
       &   'powno3','PoWa nitrate',' ','mol N m-3')
  call wrtsdm(jpowasi(iogrp),SDM_POWASI(iogrp),rnacc*1e3,0.,cmpflg,             &
       &   'powsi','PoWa silicate',' ','mol Si m-3')
  call wrtsdm(jssso12(iogrp),SDM_SSSO12(iogrp),rnacc*1e3,0.,cmpflg,             &
       &   'ssso12','Sediment detritus',' ','mol P m-3')
  call wrtsdm(jssssil(iogrp),SDM_SSSSIL(iogrp),rnacc*1e3,0.,cmpflg,             &
       &   'ssssil','Sediment silicate',' ','mol Si m-3')
  call wrtsdm(jsssc12(iogrp),SDM_SSSC12(iogrp),rnacc*1e3,0.,cmpflg,             &
       &   'sssc12','Sediment CaCO3',' ','mol C m-3')
  call wrtsdm(jssster(iogrp),SDM_SSSTER(iogrp),rnacc*1e3,0.,cmpflg,             &
       &   'ssster','Sediment clay',' ','mol m-3')

  ! --- Store sediment burial fields
  call wrtbur(jburssso12(iogrp),BUR_SSSO12(iogrp),rnacc*1e3,0.,                 &
       &   cmpflg,'buro12','Burial org carbon',' ','mol P m-2')
  call wrtbur(jbursssc12(iogrp),BUR_SSSC12(iogrp),rnacc*1e3,0.,                 &
       &   cmpflg,'burc12','Burial calcium ',' ','mol C m-2')
  call wrtbur(jburssssil(iogrp),BUR_SSSSIL(iogrp),rnacc*1e3,0.,                 &
       &   cmpflg,'bursil','Burial silicate',' ','mol Si m-2')
  call wrtbur(jburssster(iogrp),BUR_SSSTER(iogrp),rnacc*1e3,0.,                 &
       &   cmpflg,'burter','Burial clay',' ','mol  m-2')
#endif
#if defined(extNcycle) && ! defined(sedbypass)
  call wrtsdm(jpownh4(iogrp),SDM_POWNH4(iogrp),rnacc*1e3,0.,cmpflg,             &
       &   'pownh4','PoWa ammonium',' ','mol N m-3')
  call wrtsdm(jpown2o(iogrp),SDM_POWN2O(iogrp),rnacc*1e3,0.,cmpflg,             &
       &   'pown2o','PoWa nitrous oxide',' ','mol N2O m-3')
  call wrtsdm(jpowno2(iogrp),SDM_POWNO2(iogrp),rnacc*1e3,0.,cmpflg,             &
       &   'powno2','PoWa nitrite',' ','mol N m-3')
  call wrtsdm(jsdm_nitr_NH4(iogrp),sdm_nitr_NH4(iogrp),                         &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'nh4nitrsdm','NH4 nitrification rate sediment',' ','mol N m-3 s-1')
  call wrtsdm(jsdm_nitr_NO2(iogrp),sdm_nitr_NO2(iogrp),                         &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'no2nitrsdm','NO2 nitrification rate sediment',' ','mol N m-3 s-1')
  call wrtsdm(jsdm_nitr_N2O_prod(iogrp),sdm_nitr_N2O_prod(iogrp),               &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'nitr_n2osdm','N2O prod during NH4 nitrification sediment',' ',         &
     &  'mol N2O m-3 s-1')
  call wrtsdm(jsdm_nitr_NH4_OM(iogrp),sdm_nitr_NH4_OM(iogrp),                   &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'nh4nitr_omsdm','OM production during NH4 nitrification sediment',' ',  &
     &  'mol P m-3 s-1')
  call wrtsdm(jsdm_nitr_NO2_OM(iogrp),sdm_nitr_NO2_OM(iogrp),                   &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'no2nitr_omsdm','OM production during NO2 nitrification sediment',' ',  &
     &  'mol P m-3 s-1')
  call wrtsdm(jsdm_denit_NO3(iogrp),sdm_denit_NO3(iogrp),                       &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'no3denitsdm','NO3 denitrification rate sediment',' ','mol N m-3 s-1')
  call wrtsdm(jsdm_denit_NO2(iogrp),sdm_denit_NO2(iogrp),                       &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'no2denitsdm','NO2 denitrification rate sediment',' ','mol N m-3 s-1')
  call wrtsdm(jsdm_denit_N2O(iogrp),sdm_denit_N2O(iogrp),                       &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'n2odenitsdm','N2O denitrification rate sediment',' ','mol N2O m-3 s-1')
  call wrtsdm(jsdm_DNRA_NO2(iogrp),sdm_DNRA_NO2(iogrp),                         &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'no2dnrasdm','NO2 DNRA rate sediment',' ','mol N m-3 s-1')
  call wrtsdm(jsdm_anmx_N2_prod(iogrp),sdm_anmx_N2_prod(iogrp),                 &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'anmx_n2sdm','Anammox N2 production rate sediment',' ','mol N2 m-3 s-1')
  call wrtsdm(jsdm_anmx_OM_prod(iogrp),sdm_anmx_OM_prod(iogrp),                 &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'anmx_omsdm','Anammox OM production rate sediment',' ','mol P m-3 s-1')
  call wrtsdm(jsdm_remin_aerob(iogrp),sdm_remin_aerob(iogrp),                   &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'reminasdm','Aerob remineralization rate sediment',' ','mol N m-3 s-1')
  call wrtsdm(jsdm_remin_sulf(iogrp),sdm_remin_sulf(iogrp),                     &
     &  rnacc*1e3/dtbgc,0.,cmpflg,                                              &
     &  'reminssdm','Sulfate remineralization rate sediment',' ','mol P m-3 s-1') 
#endif

  ! --- close netcdf file
  call ncfcls

  ! --- Initialise fields
  call inisrf(jkwco2(iogrp),0.)
  call inisrf(jpco2(iogrp),0.)
  call inisrf(jdmsflux(iogrp),0.)
  call inisrf(jco2fxd(iogrp),0.)
  call inisrf(jco2fxu(iogrp),0.)
  call inisrf(joxflux(iogrp),0.)
  call inisrf(jniflux(iogrp),0.)
  call inisrf(jn2ofx(iogrp),0.)
  call inisrf(jdms(iogrp),0.)
  call inisrf(jdmsprod(iogrp),0.)
  call inisrf(jdms_bac(iogrp),0.)
  call inisrf(jdms_uv(iogrp),0.)
  call inisrf(jexport(iogrp),0.)
  call inisrf(jexposi(iogrp),0.)
  call inisrf(jexpoca(iogrp),0.)
  call inisrf(jsrfdic(iogrp),0.)
  call inisrf(jsrfalkali(iogrp),0.)
  call inisrf(jsrfphosph(iogrp),0.)
  call inisrf(jsrfoxygen(iogrp),0.)
  call inisrf(jsrfano3(iogrp),0.)
  call inisrf(jsrfsilica(iogrp),0.)
  call inisrf(jsrfiron(iogrp),0.)
  call inisrf(jsrfphyto(iogrp),0.)
  call inisrf(jintphosy(iogrp),0.)
  call inisrf(jintnfix(iogrp),0.)
  call inisrf(jintdnit(iogrp),0.)
  call inisrf(jcarflx0100(iogrp),0.)
  call inisrf(jcarflx0500(iogrp),0.)
  call inisrf(jcarflx1000(iogrp),0.)
  call inisrf(jcarflx2000(iogrp),0.)
  call inisrf(jcarflx4000(iogrp),0.)
  call inisrf(jcarflx_bot(iogrp),0.)
  call inisrf(jbsiflx0100(iogrp),0.)
  call inisrf(jbsiflx0500(iogrp),0.)
  call inisrf(jbsiflx1000(iogrp),0.)
  call inisrf(jbsiflx2000(iogrp),0.)
  call inisrf(jbsiflx4000(iogrp),0.)
  call inisrf(jbsiflx_bot(iogrp),0.)
  call inisrf(jcalflx0100(iogrp),0.)
  call inisrf(jcalflx0500(iogrp),0.)
  call inisrf(jcalflx1000(iogrp),0.)
  call inisrf(jcalflx2000(iogrp),0.)
  call inisrf(jcalflx4000(iogrp),0.)
  call inisrf(jcalflx_bot(iogrp),0.)
#ifndef sedbypass
  call inisrf(jsediffic(iogrp),0.)
  call inisrf(jsediffal(iogrp),0.)
  call inisrf(jsediffph(iogrp),0.)
  call inisrf(jsediffox(iogrp),0.)
  call inisrf(jsediffn2(iogrp),0.)
  call inisrf(jsediffno3(iogrp),0.)
  call inisrf(jsediffsi(iogrp),0.)
#endif
#ifdef cisonew
  call inisrf(jco213fxd(iogrp),0.)
  call inisrf(jco213fxu(iogrp),0.)
  call inisrf(jco214fxd(iogrp),0.)
  call inisrf(jco214fxu(iogrp),0.)
#endif
#ifdef CFC
  call inisrf(jcfc11fx(iogrp),0.)
  call inisrf(jcfc12fx(iogrp),0.)
  call inisrf(jsf6fx(iogrp),0.)
#endif
#ifdef natDIC
  call inisrf(jsrfnatdic(iogrp),0.)
  call inisrf(jsrfnatalk(iogrp),0.)
  call inisrf(jnatpco2(iogrp),0.)
  call inisrf(jnatco2fx(iogrp),0.)
#endif
#ifdef BROMO
  call inisrf(jsrfbromo(iogrp),0.)
  call inisrf(jbromofx(iogrp),0.)
  call inisrf(jbromo_prod(iogrp),0.)
  call inisrf(jbromo_uv(iogrp),0.)
  call inisrf(jatmbromo(iogrp),0.)
#endif


  call inisrf(jatmco2(iogrp),0.)
#if defined(BOXATM)
  call inisrf(jatmo2(iogrp),0.)
  call inisrf(jatmn2(iogrp),0.)
#endif
#ifdef cisonew
  call inisrf(jatmc13(iogrp),0.)
  call inisrf(jatmc14(iogrp),0.)
#endif
#ifdef extNcycle
  call inisrf(jsrfanh4(iogrp),0.)
  call inisrf(jsrfano2(iogrp),0.)
  call inisrf(janh3fx(iogrp),0.)
#endif
#if defined(extNcycle) && ! defined(sedbypass)
  call inisrf(jsediffnh4(iogrp),0.)
  call inisrf(jsediffn2o(iogrp),0.)
  call inisrf(jsediffno2(iogrp),0.)
#endif
  call inilyr(jdp(iogrp),0.)
  call inilyr(jdic(iogrp),0.)
  call inilyr(jalkali(iogrp),0.)
  call inilyr(jphosy(iogrp),0.)
  call inilyr(jphosph(iogrp),0.)
  call inilyr(joxygen(iogrp),0.)
  call inilyr(jano3(iogrp),0.)
  call inilyr(jsilica(iogrp),0.)
  call inilyr(jdoc(iogrp),0.)
  call inilyr(jphyto(iogrp),0.)
  call inilyr(jgrazer(iogrp),0.)
  call inilyr(jpoc(iogrp),0.)
  call inilyr(jcalc(iogrp),0.)
  call inilyr(jopal(iogrp),0.)
  call inilyr(jiron(iogrp),0.)
  call inilyr(jco3(iogrp),0.)
  call inilyr(jph(iogrp),0.)
  call inilyr(jomegaa(iogrp),0.)
  call inilyr(jomegac(iogrp),0.)
  call inilyr(jn2o(iogrp),0.)
  call inilyr(jprefo2(iogrp),0.)
  call inilyr(jo2sat(iogrp),0.)
  call inilyr(jprefpo4(iogrp),0.)
  call inilyr(jprefalk(iogrp),0.)
  call inilyr(jprefdic(iogrp),0.)
  call inilyr(jdicsat(iogrp),0.)
#ifdef cisonew
  call inilyr(jdic13(iogrp),0.)
  call inilyr(jdic14(iogrp),0.)
  call inilyr(jd13c(iogrp),0.)
  call inilyr(jd14c(iogrp),0.)
  call inilyr(jbigd14c(iogrp),0.)
  call inilyr(jpoc13(iogrp),0.)
  call inilyr(jdoc13(iogrp),0.)
  call inilyr(jcalc13(iogrp),0.)
  call inilyr(jphyto13(iogrp),0.)
  call inilyr(jgrazer13(iogrp),0.)
#endif
#ifdef AGG
  call inilyr(jnos(iogrp),0.)
  call inilyr(jwphy(iogrp),0.)
  call inilyr(jwnos(iogrp),0.)
  call inilyr(jeps(iogrp),0.)
  call inilyr(jasize(iogrp),0.)
#endif
#ifdef CFC
  call inilyr(jcfc11(iogrp),0.)
  call inilyr(jcfc12(iogrp),0.)
  call inilyr(jsf6(iogrp),0.)
#endif
#ifdef natDIC
  call inilyr(jnatco3(iogrp),0.)
  call inilyr(jnatalkali(iogrp),0.)
  call inilyr(jnatdic(iogrp),0.)
  call inilyr(jnatcalc(iogrp),0.)
  call inilyr(jnatph(iogrp),0.)
  call inilyr(jnatomegaa(iogrp),0.)
  call inilyr(jnatomegac(iogrp),0.)
#endif
#ifdef BROMO
  call inilyr(jbromo(iogrp),0.)
#endif
#ifdef extNcycle
  call inilyr(janh4(iogrp),0.)       
  call inilyr(jano2(iogrp),0.)      
  call inilyr(jnitr_NH4(iogrp),0.)
  call inilyr(jnitr_NO2(iogrp),0.)
  call inilyr(jnitr_N2O_prod(iogrp),0.)
  call inilyr(jnitr_NH4_OM(iogrp),0.)
  call inilyr(jnitr_NO2_OM(iogrp),0.)
  call inilyr(jdenit_NO3(iogrp),0.)
  call inilyr(jdenit_NO2(iogrp),0.)
  call inilyr(jdenit_N2O(iogrp),0.)
  call inilyr(jDNRA_NO2(iogrp),0.)
  call inilyr(janmx_N2_prod(iogrp),0.)
  call inilyr(janmx_OM_prod(iogrp),0.)
  call inilyr(jphosy_NH4(iogrp),0.)
  call inilyr(jphosy_NO3(iogrp),0.)
  call inilyr(jremin_aerob(iogrp),0.)
  call inilyr(jremin_sulf(iogrp),0.) 
#endif
  !   M4AGO
  call inilyr(jagg_ws(iogrp),0.)
  call inilyr(jdynvis(iogrp),0.)
  call inilyr(jagg_stick(iogrp),0.)
  call inilyr(jagg_stickf(iogrp),0.)
  call inilyr(jagg_dmax(iogrp),0.)
  call inilyr(jagg_avdp(iogrp),0.)
  call inilyr(jagg_avrhop(iogrp),0.)
  call inilyr(jagg_avdC(iogrp),0.)
  call inilyr(jagg_df(iogrp),0.)
  call inilyr(jagg_b(iogrp),0.)
  call inilyr(jagg_Vrhof(iogrp),0.)
  call inilyr(jagg_Vpor(iogrp),0.)

  call inilvl(jlvldic(iogrp),0.)
  call inilvl(jlvlalkali(iogrp),0.)
  call inilvl(jlvlphosy(iogrp),0.)
  call inilvl(jlvlphosph(iogrp),0.)
  call inilvl(jlvloxygen(iogrp),0.)
  call inilvl(jlvlano3(iogrp),0.)
  call inilvl(jlvlsilica(iogrp),0.)
  call inilvl(jlvldoc(iogrp),0.)
  call inilvl(jlvlphyto(iogrp),0.)
  call inilvl(jlvlgrazer(iogrp),0.)
  call inilvl(jlvlpoc(iogrp),0.)
  call inilvl(jlvlcalc(iogrp),0.)
  call inilvl(jlvlopal(iogrp),0.)
  call inilvl(jlvliron(iogrp),0.)
  call inilvl(jlvlco3(iogrp),0.)
  call inilvl(jlvlph(iogrp),0.)
  call inilvl(jlvlomegaa(iogrp),0.)
  call inilvl(jlvlomegac(iogrp),0.)
  call inilvl(jlvln2o(iogrp),0.)
  call inilvl(jlvlprefo2(iogrp),0.)
  call inilvl(jlvlo2sat(iogrp),0.)
  call inilvl(jlvlprefpo4(iogrp),0.)
  call inilvl(jlvlprefalk(iogrp),0.)
  call inilvl(jlvlprefdic(iogrp),0.)
  call inilvl(jlvldicsat(iogrp),0.)
#ifdef cisonew
  call inilvl(jlvldic13(iogrp),0.)
  call inilvl(jlvldic14(iogrp),0.)
  call inilvl(jlvld13c(iogrp),0.)
  call inilvl(jlvld14c(iogrp),0.)
  call inilvl(jlvlbigd14c(iogrp),0.)
  call inilvl(jlvlpoc13(iogrp),0.)
  call inilvl(jlvldoc13(iogrp),0.)
  call inilvl(jlvlcalc13(iogrp),0.)
  call inilvl(jlvlphyto13(iogrp),0.)
  call inilvl(jlvlgrazer13(iogrp),0.)
#endif
#ifdef AGG
  call inilvl(jlvlnos(iogrp),0.)
  call inilvl(jlvlwphy(iogrp),0.)
  call inilvl(jlvlwnos(iogrp),0.)
  call inilvl(jlvleps(iogrp),0.)
  call inilvl(jlvlasize(iogrp),0.)
#endif
#ifdef CFC
  call inilvl(jlvlcfc11(iogrp),0.)
  call inilvl(jlvlcfc12(iogrp),0.)
  call inilvl(jlvlsf6(iogrp),0.)
#endif
#ifdef natDIC
  call inilvl(jlvlnatco3(iogrp),0.)
  call inilvl(jlvlnatalkali(iogrp),0.)
  call inilvl(jlvlnatdic(iogrp),0.)
  call inilvl(jlvlnatcalc(iogrp),0.)
  call inilvl(jlvlnatph(iogrp),0.)
  call inilvl(jlvlnatomegaa(iogrp),0.)
  call inilvl(jlvlnatomegac(iogrp),0.)
#endif
#ifdef BROMO
  call inilvl(jlvlbromo(iogrp),0.)
#endif
#ifdef extNcycle
  call inilvl(jlvlanh4(iogrp),0.)  
  call inilvl(jlvlano2(iogrp),0.)  
  call inilvl(jlvl_nitr_NH4(iogrp),0.)
  call inilvl(jlvl_nitr_NO2(iogrp),0.)
  call inilvl(jlvl_nitr_N2O_prod(iogrp),0.)
  call inilvl(jlvl_nitr_NH4_OM(iogrp),0.)
  call inilvl(jlvl_nitr_NO2_OM(iogrp),0.)
  call inilvl(jlvl_denit_NO3(iogrp),0.)
  call inilvl(jlvl_denit_NO2(iogrp),0.)
  call inilvl(jlvl_denit_N2O(iogrp),0.)
  call inilvl(jlvl_DNRA_NO2(iogrp),0.)
  call inilvl(jlvl_anmx_N2_prod(iogrp),0.)
  call inilvl(jlvl_anmx_OM_prod(iogrp),0.)
  call inilvl(jlvl_phosy_NH4(iogrp),0.)
  call inilvl(jlvl_phosy_NO3(iogrp),0.)
  call inilvl(jlvl_remin_aerob(iogrp),0.)
  call inilvl(jlvl_remin_sulf(iogrp),0.) 
#endif
  !  M4AGO
  call inilvl(jlvl_agg_ws(iogrp),0.)
  call inilvl(jlvl_dynvis(iogrp),0.)
  call inilvl(jlvl_agg_stick(iogrp),0.)
  call inilvl(jlvl_agg_stickf(iogrp),0.)
  call inilvl(jlvl_agg_dmax(iogrp),0.)
  call inilvl(jlvl_agg_avdp(iogrp),0.)
  call inilvl(jlvl_agg_avrhop(iogrp),0.)
  call inilvl(jlvl_agg_avdC(iogrp),0.)
  call inilvl(jlvl_agg_df(iogrp),0.)
  call inilvl(jlvl_agg_b(iogrp),0.)
  call inilvl(jlvl_agg_Vrhof(iogrp),0.)
  call inilvl(jlvl_agg_Vpor(iogrp),0.)

#ifndef sedbypass
  call inisdm(jpowaic(iogrp),0.)
  call inisdm(jpowaal(iogrp),0.)
  call inisdm(jpowaph(iogrp),0.)
  call inisdm(jpowaox(iogrp),0.)
  call inisdm(jpown2(iogrp),0.)
  call inisdm(jpowno3(iogrp),0.)
  call inisdm(jpowasi(iogrp),0.)
  call inisdm(jssso12(iogrp),0.)
  call inisdm(jssssil(iogrp),0.)
  call inisdm(jsssc12(iogrp),0.)
  call inisdm(jssster(iogrp),0.)

  call inibur(jburssso12(iogrp),0.)
  call inibur(jbursssc12(iogrp),0.)
  call inibur(jburssssil(iogrp),0.)
  call inibur(jburssster(iogrp),0.)
#endif
#if defined(extNcycle) && ! defined(sedbypass)
  call inisdm(jpownh4(iogrp),0.)
  call inisdm(jpown2o(iogrp),0.)
  call inisdm(jpowno2(iogrp),0.)
  call inisdm(jsdm_nitr_NH4(iogrp),0.)
  call inisdm(jsdm_nitr_NO2(iogrp),0.)
  call inisdm(jsdm_nitr_N2O_prod(iogrp),0.) 
  call inisdm(jsdm_nitr_NH4_OM(iogrp),0.) 
  call inisdm(jsdm_nitr_NO2_OM(iogrp),0.) 
  call inisdm(jsdm_denit_NO3(iogrp),0.) 
  call inisdm(jsdm_denit_NO2(iogrp),0.) 
  call inisdm(jsdm_denit_N2O(iogrp),0.)
  call inisdm(jsdm_DNRA_NO2(iogrp),0.)
  call inisdm(jsdm_anmx_N2_prod(iogrp),0.)      
  call inisdm(jsdm_anmx_OM_prod(iogrp),0.)
  call inisdm(jsdm_remin_aerob(iogrp),0.)
  call inisdm(jsdm_remin_sulf(iogrp),0.)
#endif

  nacc_bgc(iogrp)=0

end subroutine ncwrt_bgc


subroutine hamoccvardef(iogrp,timeunits,calendar,cmpflg)
  use mod_nctools, only:ncdefvar,ncattr,ncfopn,ncdimc,ncdims,                   &
       &   nctime,ncfcls,ncedef,ncdefvar3d,ndouble

  use mo_bgcmean, only: srf_kwco2,srf_pco2,srf_dmsflux,srf_co2fxd,              &
       &   srf_co2fxu,srf_oxflux,srf_niflux,srf_dms,srf_dmsprod,                &
       &   srf_dms_bac,srf_dms_uv,srf_export,srf_exposi,srf_expoca,             &
       &   srf_dic,srf_alkali,srf_phosph,srf_oxygen,srf_ano3,srf_silica,        &
       &   srf_iron,srf_phyto,int_phosy,int_nfix,int_dnit,flx_car0100,          &
       &   flx_car0500,flx_car1000,flx_car2000,flx_car4000,flx_car_bot,         &
       &   flx_bsi0100,flx_bsi0500,flx_bsi1000,flx_bsi2000,flx_bsi4000,         &
       &   flx_bsi_bot,flx_cal0100,flx_cal0500,flx_cal1000,flx_cal2000,         &
       &   flx_cal4000,flx_cal_bot,flx_sediffic,flx_sediffal,                   &
       &   flx_sediffph,flx_sediffox,flx_sediffn2,flx_sediffno3,                &
       &   flx_sediffsi,srf_n2ofx,srf_atmco2,lyr_dp,lyr_dic,                    &
       &   lyr_alkali,lyr_phosph,lyr_oxygen,lyr_ano3,lyr_silica,lyr_doc,        &
       &   lyr_phyto,lyr_grazer,lyr_poc,lyr_calc,lyr_opal,lyr_iron,             &
       &   lyr_phosy,lyr_co3,lyr_ph,lyr_omegaa,lyr_omegac,lyr_n2o,              &
       &   lyr_prefo2,lyr_o2sat,lyr_prefpo4,lyr_prefalk,lyr_prefdic,            &
       &   lyr_dicsat,lvl_dic,lvl_alkali,lvl_phosph,lvl_oxygen,lvl_ano3,        &
       &   lvl_silica,lvl_doc,lvl_phyto,lvl_grazer,lvl_poc,lvl_calc,            &
       &   lvl_opal,lvl_iron,lvl_phosy,lvl_co3,lvl_ph,lvl_omegaa,               &
       &   lvl_omegac,lvl_n2o,lvl_prefo2,lvl_o2sat,lvl_prefpo4,                 &
       &   lvl_prefalk,lvl_prefdic,lvl_dicsat,                                  &
       &   lyr_agg_ws,lyr_dynvis,lyr_agg_stick,                                 &
       &   lyr_agg_stickf,lyr_agg_dmax,lyr_agg_avdp,                            &
       &   lyr_agg_avrhop,lyr_agg_avdC,lyr_agg_df,                              &
       &   lyr_agg_b,lyr_agg_Vrhof,lyr_agg_Vpor,                                &
       &   lvl_agg_ws,lvl_dynvis,lvl_agg_stick,                                 &
       &   lvl_agg_stickf,lvl_agg_dmax,lvl_agg_avdp,                            &
       &   lvl_agg_avrhop,lvl_agg_avdC,lvl_agg_df,                              &
       &   lvl_agg_b,lvl_agg_Vrhof,lvl_agg_Vpor
#ifdef AGG
  use mo_bgcmean, only: lyr_nos,lyr_wphy,lyr_wnos,lyr_eps,                      &
       &   lyr_asize,lvl_nos,lvl_wphy,lvl_wnos,lvl_eps,lvl_asize
#endif
#if defined(BOXATM)
  use mo_bgcmean, only: srf_atmo2,srf_atmn2
#endif

#ifdef BROMO
  use mo_bgcmean, only:srf_bromo,srf_bromofx,int_bromopro,                      &
       &   int_bromouv,srf_atmbromo,lyr_bromo,lvl_bromo
#endif
#ifdef CFC
  use mo_bgcmean, only: srf_cfc11,srf_cfc12,srf_sf6,lyr_cfc11,                  &
       &   lyr_cfc12,lyr_sf6,lvl_cfc11,lvl_cfc12,lvl_sf6
#endif
#ifdef cisonew
  use mo_bgcmean, only: srf_co213fxd,srf_co213fxu,srf_co214fxd,                 &
       &   srf_co214fxu,srf_atmc13,srf_atmc14,lyr_dic13,lyr_dic14,              &
       &   lyr_d13c,lyr_d14c,lyr_bigd14c,lyr_poc13,lyr_doc13,                   &
       &   lyr_calc13,lyr_phyto13,lyr_grazer13,lvl_dic13,lvl_dic14,             &
       &   lvl_d13c,lvl_d14c,lvl_bigd14c,lvl_poc13,lvl_doc13,                   &
       &   lvl_calc13,lvl_phyto13,lvl_grazer13
#endif
#ifdef natDIC
  use mo_bgcmean, only: srf_natdic,srf_natalkali,srf_natpco2,                   &
       &   srf_natco2fx,lyr_natco3,lyr_natalkali,lyr_natdic,                    &
       &   lyr_natcalc,lyr_natph,lyr_natomegaa,lyr_natomegac,                   &
       &   lvl_natalkali,lvl_natdic,lvl_natcalc,lvl_natph,                      &
       &   lvl_natomegaa,lvl_natomegac,lvl_natco3
#endif
#ifndef sedbypass
  use mo_bgcmean, only: sdm_powaic,sdm_powaal,sdm_powaph,sdm_powaox,            &
       &   sdm_pown2,sdm_powno3,sdm_powasi,sdm_ssso12,sdm_ssssil,               &
       &   sdm_sssc12,sdm_ssster,bur_ssso12,bur_sssc12,bur_ssssil,              &
       &   bur_ssster
#endif
#ifdef extNcycle
  use mo_bgcmean, only: janh4,jano2,jlvlanh4,jlvlano2,jsrfanh4,                 &
       &                     jsrfano2,janh3fx,srf_anh4,srf_ano2,                &
       &                     srf_anh3fx,lyr_anh4,lyr_ano2,lvl_anh4,             &
       &                     lvl_ano2,                                          &
       &                     LYR_nitr_NH4,LYR_nitr_NO2,LYR_nitr_N2O_prod,       &
       &                     LYR_nitr_NH4_OM,LYR_nitr_NO2_OM,                   &
       &                     LYR_denit_NO3,LYR_denit_NO2,LYR_denit_N2O,         &
       &                     LYR_DNRA_NO2,LYR_anmx_N2_prod,                     &
       &                     LYR_anmx_OM_prod,LYR_phosy_NH4,                    &
       &                     LYR_phosy_NO3,LYR_remin_aerob,LYR_remin_sulf,      &
       &                     LVL_nitr_NH4,LVL_nitr_NO2,LVL_nitr_N2O_prod,       &
       &                     LVL_nitr_NH4_OM,LVL_nitr_NO2_OM,                   &
       &                     LVL_denit_NO3,LVL_denit_NO2,LVL_denit_N2O,         &
       &                     LVL_DNRA_NO2,LVL_anmx_N2_prod,                     &
       &                     LVL_anmx_OM_prod,LVL_phosy_NH4,                    &
       &                     LVL_phosy_NO3,LVL_remin_aerob,LVL_remin_sulf,      &
       &                     jnitr_NH4,jnitr_NO2,jnitr_N2O_prod,                &
       &                     jnitr_NH4_OM,jnitr_NO2_OM,jdenit_NO3,              &
       &                     jdenit_NO2,jdenit_N2O,jDNRA_NO2,                   &
       &                     janmx_N2_prod,janmx_OM_prod,jphosy_NH4,            &
       &                     jphosy_NO3,jremin_aerob,jremin_sulf,               &
       &                     jlvl_nitr_NH4,jlvl_nitr_NO2,                       &
       &                     jlvl_nitr_N2O_prod,jlvl_nitr_NH4_OM,               &
       &                     jlvl_nitr_NO2_OM,jlvl_denit_NO3,                   &
       &                     jlvl_denit_NO2,jlvl_denit_N2O,jlvl_DNRA_NO2,       &
       &                     jlvl_anmx_N2_prod,jlvl_anmx_OM_prod,               &
       &                     jlvl_phosy_NH4,jlvl_phosy_NO3,                     &
       &                     jlvl_remin_aerob,jlvl_remin_sulf      
#endif
#if defined(extNcycle) && ! defined(sedbypass)
  use mo_bgcmean, only: jpownh4,jpown2o,jpowno2,jsdm_nitr_NH4,jsdm_nitr_NO2,    &
       &                jsdm_nitr_N2O_prod,jsdm_nitr_NH4_OM,jsdm_nitr_NO2_OM,   &
       &                jsdm_denit_NO3,jsdm_denit_NO2,jsdm_denit_N2O,           &
       &                jsdm_DNRA_NO2,jsdm_anmx_N2_prod,jsdm_anmx_OM_prod,      &
       &                jsdm_remin_aerob,jsdm_remin_sulf, SDM_POWNH4,SDM_POWN2O,&
       &                SDM_POWNO2,SDM_nitr_NH4,SDM_nitr_NO2,SDM_nitr_N2O_prod, &
       &                SDM_nitr_NH4_OM,SDM_nitr_NO2_OM,SDM_denit_NO3,          &
       &                SDM_denit_NO2,SDM_denit_N2O,SDM_DNRA_NO2,               &
       &                SDM_anmx_N2_prod,SDM_anmx_OM_prod,SDM_remin_aerob,      &
       &                SDM_remin_sulf,jsediffnh4,jsediffn2o,jsediffno2,        &
       &                FLX_SEDIFFNH4,FLX_SEDIFFN2O,FLX_SEDIFFNO2 
#endif



  implicit none

  integer iogrp,cmpflg
  character timeunits*30,calendar*19
  call ncdefvar('time','time',ndouble,0)
  call ncattr('long_name','time')
  call ncattr('units',timeunits)
  call ncattr('calendar',calendar)
  call ncdefvar('sigma','sigma',ndouble,8)
  call ncattr('long_name','Potential density')
  call ncattr('standard_name','sea_water_sigma_theta')
  call ncattr('units','kg m-3')
  call ncattr('positive','down')
  call ncdefvar('depth','depth',ndouble,8)
  call ncattr('long_name','z level')
  call ncattr('units','m')
  call ncattr('positive','down')
  call ncattr('bounds','depth_bnds')
  call ncdefvar('depth_bnds','bounds depth',ndouble,8)
  call ncdefvar3d(SRF_KWCO2(iogrp),cmpflg,'p',                                  &
       &   'kwco2',' ',' ',' ',0)
  call ncdefvar3d(SRF_PCO2(iogrp),cmpflg,'p',                                   &
       &   'pco2','Surface PCO2',' ','uatm',0)
  call ncdefvar3d(SRF_DMSFLUX(iogrp),                                           &
       &   cmpflg,'p','dmsflux','DMS flux',' ','mol DMS m-2 s-1',0)
  call ncdefvar3d(SRF_CO2FXD(iogrp),                                            &
       &   cmpflg,'p','co2fxd','Downward CO2 flux',' ','kg C m-2 s-1',0)
  call ncdefvar3d(SRF_CO2FXU(iogrp),                                            &
       &   cmpflg,'p','co2fxu','Upward CO2 flux',' ','kg C m-2 s-1',0)
  call ncdefvar3d(SRF_OXFLUX(iogrp),                                            &
       &   cmpflg,'p','fgo2','Oxygen flux',' ','mol O2 m-2 s-1',0)
  call ncdefvar3d(SRF_NIFLUX(iogrp),                                            &
       &   cmpflg,'p','fgn2','Nitrogen flux',' ','mol N2 m-2 s-1',0)
  call ncdefvar3d(SRF_DMS(iogrp),cmpflg,'p',                                    &
       &   'dms','DMS',' ','kmol DMS m-3',0)
  call ncdefvar3d(SRF_DMSPROD(iogrp),cmpflg,'p',                                &
       &   'dmsprod','DMS production from phytoplankton production',' ',        &
       &   'mol DMS m-2 s-1',0)
  call ncdefvar3d(SRF_DMS_BAC(iogrp),cmpflg,'p',                                &
       &   'dms_bac','DMS bacterial consumption',' ','mol DMS m-2 s-1',0)
  call ncdefvar3d(SRF_DMS_UV(iogrp),cmpflg,'p',                                 &
       &   'dms_uv','DMS photolysis reduction',' ','mol DMS m-2 s-1',0)
  call ncdefvar3d(SRF_EXPORT(iogrp),                                            &
       &   cmpflg,'p','epc100','Export production',' ','mol C m-2 s-1',0)
  call ncdefvar3d(SRF_EXPOSI(iogrp),cmpflg,'p',                                 &
       &   'epsi100','Si export production',' ','mol Si m-2 s-1',0)
  call ncdefvar3d(SRF_EXPOCA(iogrp),cmpflg,'p',                                 &
       &   'epcalc100','Ca export production',' ','mol Ca m-2 s-1',0)
  call ncdefvar3d(SRF_DIC(iogrp),cmpflg,'p','srfdissic',                        &
       &   'Surface dissolved inorganic carbon',' ','mol C m-3',0)
  call ncdefvar3d(SRF_ALKALI(iogrp),cmpflg,'p','srftalk',                       &
       &   'Surface alkalinity',' ','eq m-3',0)
  call ncdefvar3d(SRF_PHOSPH(iogrp),cmpflg,'p','srfpo4',                        &
       &   'Surface phosphorus',' ','mol P m-3',0)
  call ncdefvar3d(SRF_OXYGEN(iogrp),cmpflg,'p','srfo2',                         &
       &   'Surface oxygen',' ','mol O2 m-3',0)
  call ncdefvar3d(SRF_ANO3(iogrp),cmpflg,'p','srfno3',                          &
       &   'Surface nitrate',' ','mol N m-3',0)
  call ncdefvar3d(SRF_SILICA(iogrp),cmpflg,'p','srfsi',                         &
       &   'Surface silicate',' ','mol Si m-3',0)
  call ncdefvar3d(SRF_IRON(iogrp),cmpflg,'p','srfdfe',                          &
       &   'Surface dissolved iron',' ','mol Fe m-3',0)
  call ncdefvar3d(SRF_PHYTO(iogrp),cmpflg,'p','srfphyc',                        &
       &   'Surface phytoplankton',' ','mol P m-3',0)
  call ncdefvar3d(INT_PHOSY(iogrp),cmpflg,'p','ppint',                          &
       &   'Integrated primary production',' ','mol C m-2 s-1',0)
  call ncdefvar3d(INT_NFIX(iogrp),cmpflg,'p','nfixint',                         &
       &   'Integrated nitrogen fixation',' ','mol N m-2 s-1',0)
  call ncdefvar3d(INT_DNIT(iogrp),cmpflg,'p','dnitint',                         &
       &   'Integrated denitrification',' ','mol N m-2 s-1',0)
  call ncdefvar3d(FLX_CAR0100(iogrp),cmpflg,'p','carflx0100',                   &
       &   'C flux at 100m',' ','mol C m-2 s-1',0)
  call ncdefvar3d(FLX_CAR0500(iogrp),cmpflg,'p','carflx0500',                   &
       &   'C flux at 500m',' ','mol C m-2 s-1',0)
  call ncdefvar3d(FLX_CAR1000(iogrp),cmpflg,'p','carflx1000',                   &
       &   'C flux at 1000m',' ','mol C m-2 s-1',0)
  call ncdefvar3d(FLX_CAR2000(iogrp),cmpflg,'p','carflx2000',                   &
       &   'C flux at 2000m',' ','mol C m-2 s-1',0)
  call ncdefvar3d(FLX_CAR4000(iogrp),cmpflg,'p','carflx4000',                   &
       &   'C flux at 4000m',' ','mol C m-2 s-1',0)
  call ncdefvar3d(FLX_CAR_BOT(iogrp),cmpflg,'p','carflx_bot',                   &
       &   'C flux to sediment',' ','mol C m-2 s-1',0)
  call ncdefvar3d(FLX_BSI0100(iogrp),cmpflg,'p','bsiflx0100',                   &
       &   'Opal flux at 100m',' ','mol Si m-2 s-1',0)
  call ncdefvar3d(FLX_BSI0500(iogrp),cmpflg,'p','bsiflx0500',                   &
       &   'Opal flux at 500m',' ','mol Si m-2 s-1',0)
  call ncdefvar3d(FLX_BSI1000(iogrp),cmpflg,'p','bsiflx1000',                   &
       &   'Opal flux at 1000m',' ','mol Si m-2 s-1',0)
  call ncdefvar3d(FLX_BSI2000(iogrp),cmpflg,'p','bsiflx2000',                   &
       &   'Opal flux at 2000m',' ','mol Si m-2 s-1',0)
  call ncdefvar3d(FLX_BSI4000(iogrp),cmpflg,'p','bsiflx4000',                   &
       &   'Opal flux at 4000m',' ','mol Si m-2 s-1',0)
  call ncdefvar3d(FLX_BSI_BOT(iogrp),cmpflg,'p','bsiflx_bot',                   &
       &   'Opal flux to sediment',' ','mol Si m-2 s-1',0)
  call ncdefvar3d(FLX_CAL0100(iogrp),cmpflg,'p','calflx0100',                   &
       &   'CaCO3 flux at 100m',' ','mol Ca m-2 s-1',0)
  call ncdefvar3d(FLX_CAL0500(iogrp),cmpflg,'p','calflx0500',                   &
       &   'CaCO3 flux at 500m',' ','mol Ca m-2 s-1',0)
  call ncdefvar3d(FLX_CAL1000(iogrp),cmpflg,'p','calflx1000',                   &
       &   'CaCO3 flux at 1000m',' ','mol Ca m-2 s-1',0)
  call ncdefvar3d(FLX_CAL2000(iogrp),cmpflg,'p','calflx2000',                   &
       &   'CaCO3 flux at 2000m',' ','mol Ca m-2 s-1',0)
  call ncdefvar3d(FLX_CAL4000(iogrp),cmpflg,'p','calflx4000',                   &
       &   'CaCO3 flux at 4000m',' ','mol Ca m-2 s-1',0)
  call ncdefvar3d(FLX_CAL_BOT(iogrp),cmpflg,'p','calflx_bot',                   &
       &   'CaCO3 flux to sediment',' ','mol Ca m-2 s-1',0)
  call ncdefvar3d(SRF_N2OFX(iogrp),cmpflg,'p','n2oflux',                        &
       &   'N2O flux',' ','mol N2O m-2 s-1',0)
#ifndef sedbypass
  call ncdefvar3d(FLX_SEDIFFIC(iogrp),cmpflg,'p','sedfdic',                     &
       &   'diffusive DIC flux to sediment (positive downwards)',               &
       &   ' ','mol C m-2 s-1',0)
  call ncdefvar3d(FLX_SEDIFFAL(iogrp),cmpflg,'p','sedfalk',                     &
       &   'diffusive alkalinity flux to sediment (positive downwards)',        &
       &   ' ','mol m-2 s-1',0)
  call ncdefvar3d(FLX_SEDIFFPH(iogrp),cmpflg,'p','sedfpho',                     &
       &   'diffusive phosphate flux to sediment (positive downwards)',         &
       &   ' ','mol m-2 s-1',0)
  call ncdefvar3d(FLX_SEDIFFOX(iogrp),cmpflg,'p','sedfox',                      &
       &   'diffusive oxygen flux to sediment (positive downwards)',            &
       &   ' ','mol O2 m-2 s-1',0)
  call ncdefvar3d(FLX_SEDIFFN2(iogrp),cmpflg,'p','sedfn2',                      &
       &   'diffusive N2 flux to sediment (positive downwards)',                &
       &   ' ','mol N2 m-2 s-1',0)
  call ncdefvar3d(FLX_SEDIFFNO3(iogrp),cmpflg,'p','sedfno3',                    &
       &   'diffusive nitrate flux to sediment (positive downwards)',           &
       &   ' ','mol NO3 m-2 s-1',0)
  call ncdefvar3d(FLX_SEDIFFSI(iogrp),cmpflg,'p','sedfsi',                      &
       &   'diffusive silica flux to sediment (positive downwards)',            &
       &   ' ','mol Si m-2 s-1',0)
#endif
#if defined(extNcycle) && ! defined(sedbypass)
  call ncdefvar3d(FLX_SEDIFFNH4(iogrp),cmpflg,'p','sedfnh4',                    &
       &   'diffusive ammonium flux to sediment (positive downwards)',          &
       &   ' ','mol NH4 m-2 s-1',0)
  call ncdefvar3d(FLX_SEDIFFN2O(iogrp),cmpflg,'p','sedfn2o',                    &
       &   'diffusive nitrous oxide flux to sediment (positive downwards)',     &
       &   ' ','mol N2O m-2 s-1',0)
  call ncdefvar3d(FLX_SEDIFFNO2(iogrp),cmpflg,'p','sedfno2',                    &
       &   'diffusive nitrite flux to sediment (positive downwards)',           &
       &   ' ','mol NO2 m-2 s-1',0)
#endif
#ifdef cisonew
  call ncdefvar3d(SRF_CO213FXD(iogrp),cmpflg,'p','co213fxd',                    &
       &   'Downward 13CO2 flux',' ','kg C m-2 s-1',0)
  call ncdefvar3d(SRF_CO213FXU(iogrp),cmpflg,'p','co213fxu',                    &
       &   'Upward 13CO2 flux',' ','kg C m-2 s-1',0)
  call ncdefvar3d(SRF_CO214FXD(iogrp),cmpflg,'p','co214fxd',                    &
       &   'Downward 14CO2 flux',' ','kg C m-2 s-1',0)
  call ncdefvar3d(SRF_CO214FXU(iogrp),cmpflg,'p','co214fxu',                    &
       &   'Upward 14CO2 flux',' ','kg C m-2 s-1',0)
#endif
#ifdef CFC
  call ncdefvar3d(SRF_CFC11(iogrp),cmpflg,'p','cfc11flux',                      &
       &   'CFC-11 flux',' ','mol CFC12 m-2 s-1',0)
  call ncdefvar3d(SRF_CFC12(iogrp),                                             &
       &   cmpflg,'p','cfc12flux','CFC-12 flux',' ','mol CFC12 m-2 s-1',0)
  call ncdefvar3d(SRF_SF6(iogrp),                                               &
       &   cmpflg,'p','sf6flux','SF-6 flux',' ','mol SF6 m-2 s-1',0)
#endif
#ifdef natDIC
  call ncdefvar3d(SRF_NATDIC(iogrp),cmpflg,'p','srfnatdissic',                  &
       &   'Surface natural dissolved inorganic carbon',' ','mol C m-3',0)
  call ncdefvar3d(SRF_NATALKALI(iogrp),cmpflg,'p','srfnattalk',                 &
       &   'Surface natural alkalinity',' ','eq m-3',0)
  call ncdefvar3d(SRF_NATPCO2(iogrp),cmpflg,'p',                                &
       &   'natpco2','Surface natural PCO2',' ','uatm',0)
  call ncdefvar3d(SRF_NATCO2FX(iogrp),                                          &
       &   cmpflg,'p','natco2fx','Natural CO2 flux',' ','kg C m-2 s-1',0)
#endif
#ifdef BROMO
  call ncdefvar3d(SRF_BROMO(iogrp),cmpflg,'p','srfbromo',                       &
       &   'Surface bromoform',' ','mol CHBr3 m-3',0)
  call ncdefvar3d(SRF_BROMOfx(iogrp),cmpflg,'p','bromofx',                      &
       &   'Surface bromoform flux',' ','mol CHBr3 m-2 s-1',0)
  call ncdefvar3d(INT_BROMOPRO(iogrp),cmpflg,'p','intbromoprod',                &
       &   'Integrated bromoform production',' ','mol CHBr3 m-2 s-1',0)
  call ncdefvar3d(INT_BROMOUV(iogrp),cmpflg,'p','intbromouv',                   &
       &   'Integrated bromoform loss to photolysis',' ',                       &
       &   'mol CHBr3 m-2 s-1',0)
  call ncdefvar3d(SRF_ATMBROMO(iogrp),cmpflg,'p',                               &
       &   'atmbromo','Atmospheric bromoform',' ','ppt',0)
#endif

  call ncdefvar3d(SRF_ATMCO2(iogrp),cmpflg,'p',                                 &
       &   'atmco2','Atmospheric CO2',' ','ppm',0)
#if defined(BOXATM)
  call ncdefvar3d(SRF_ATMO2(iogrp),cmpflg,'p',                                  &
       &   'atmo2','Atmospheric O2',' ','ppm',0)
  call ncdefvar3d(SRF_ATMN2(iogrp),cmpflg,'p',                                  &
       &   'atmn2','Atmospheric N2',' ','ppm',0)
#endif
#ifdef cisonew
  call ncdefvar3d(SRF_ATMC13(iogrp),cmpflg,'p',                                 &
       &   'atmc13','Atmospheric 13CO2',' ','ppm',0)
  call ncdefvar3d(SRF_ATMC14(iogrp),cmpflg,'p',                                 &
       &   'atmc14','Atmospheric 14CO2',' ','ppm',0)
#endif
#ifdef extNcycle
  call ncdefvar3d(SRF_ANH4(iogrp),cmpflg,'p','srfnh4',                          &
     &  'Surface ammonium',' ','mol N m-3',0)
  call ncdefvar3d(SRF_ANO2(iogrp),cmpflg,'p','srfno2',                          &
     &  'Surface nitrite',' ','mol N m-3',0)
  call ncdefvar3d(SRF_ANH3FX(iogrp),cmpflg,'p','nh3flux',                       &
     &  'NH3 flux',' ','mol NH3 m-2 s-1',0)
#endif
  ! --- define 3d layer fields
  call ncdefvar3d(LYR_DP(iogrp),cmpflg,'p',                                     &
       &   'pddpo','Layer thickness',' ','m',1)
  call ncdefvar3d(LYR_DIC(iogrp),cmpflg,'p',                                    &
       &   'dissic','Dissolved inorganic carbon',' ','mol C m-3',1)
  call ncdefvar3d(LYR_ALKALI(iogrp),cmpflg,'p',                                 &
       &   'talk','Alkalinity',' ','eq m-3',1)
  call ncdefvar3d(LYR_PHOSPH(iogrp),cmpflg,'p',                                 &
       &   'po4','Phosphorus',' ','mol P m-3',1)
  call ncdefvar3d(LYR_OXYGEN(iogrp),cmpflg,'p',                                 &
       &   'o2','Oxygen',' ','mol O2 m-3',1)
  call ncdefvar3d(LYR_ANO3(iogrp),cmpflg,'p',                                   &
       &   'no3','Nitrate',' ','mol N m-3',1)
  call ncdefvar3d(LYR_SILICA(iogrp),cmpflg,'p',                                 &
       &   'si','Silicate',' ','mol Si m-3',1)
  call ncdefvar3d(LYR_DOC(iogrp),cmpflg,'p',                                    &
       &   'dissoc','Dissolved organic carbon',' ','mol P m-3',1)
  call ncdefvar3d(LYR_PHYTO(iogrp),cmpflg,'p',                                  &
       &   'phyc','Phytoplankton',' ','mol P m-3',1)
  call ncdefvar3d(LYR_GRAZER(iogrp),cmpflg,'p',                                 &
       &   'zooc','Zooplankton',' ','mol P m-3',1)
  call ncdefvar3d(LYR_POC(iogrp),cmpflg,'p',                                    &
       &   'detoc','Detritus',' ','mol P m-3',1)
  call ncdefvar3d(LYR_CALC(iogrp),cmpflg,'p',                                   &
       &   'calc','CaCO3 shells',' ','mol C m-3',1)
  call ncdefvar3d(LYR_OPAL(iogrp),cmpflg,'p',                                   &
       &   'opal','Opal shells',' ','mol Si m-3',1)
  call ncdefvar3d(LYR_IRON(iogrp),cmpflg,'p',                                   &
       &   'dfe','Dissolved iron',' ','mol Fe m-3',1)
  call ncdefvar3d(LYR_PHOSY(iogrp),cmpflg,'p',                                  &
       &   'pp','Primary production',' ','mol C m-3 s-1',1)
  call ncdefvar3d(LYR_CO3(iogrp),cmpflg,'p',                                    &
       &   'co3','Carbonate ions',' ','mol C m-3',1)
  call ncdefvar3d(LYR_PH(iogrp),cmpflg,'p',                                     &
       &   'ph','pH',' ','-log10([h+])',1)
  call ncdefvar3d(LYR_OMEGAA(iogrp),cmpflg,'p',                                 &
       &   'omegaa','OmegaA',' ','1',1)
  call ncdefvar3d(LYR_OMEGAC(iogrp),cmpflg,'p',                                 &
       &   'omegac','OmegaC',' ','1',1)
  call ncdefvar3d(LYR_N2O(iogrp),cmpflg,'p',                                    &
       &   'n2o','N2O',' ','mol N2O m-3',1)
  call ncdefvar3d(LYR_PREFO2(iogrp),cmpflg,'p',                                 &
       &   'p_o2','Preformed oxygen',' ','mol O2 m-3',1)
  call ncdefvar3d(LYR_O2SAT(iogrp),cmpflg,'p',                                  &
       &   'satoxy','Saturated oxygen',' ','mol O2 m-3',1)
  call ncdefvar3d(LYR_PREFPO4(iogrp),cmpflg,'p',                                &
       &   'p_po4','Preformed phosphorus',' ','mol P m-3',1)
  call ncdefvar3d(LYR_PREFALK(iogrp),cmpflg,'p',                                &
       &   'p_talk','Preformed alkalinity',' ','eq m-3',1)
  call ncdefvar3d(LYR_PREFDIC(iogrp),cmpflg,'p',                                &
       &   'p_dic','Preformed DIC',' ','mol C m-3',1)
  call ncdefvar3d(LYR_DICSAT(iogrp),cmpflg,'p',                                 &
       &   'sat_dic','Saturated DIC',' ','mol C m-3',1)
#ifdef cisonew
  call ncdefvar3d(LYR_DIC13(iogrp),cmpflg,'p',                                  &
       &   'dissic13','Dissolved C13',' ','mol 13C m-3',1)
  call ncdefvar3d(LYR_DIC14(iogrp),cmpflg,'p',                                  &
       &   'dissic14','Dissolved C14',' ','mol 14C m-3',1)
  call ncdefvar3d(LYR_D13C(iogrp),cmpflg,'p',                                   &
       &   'delta13c','delta13C of DIC',' ','permil',1)
  call ncdefvar3d(LYR_D14C(iogrp),cmpflg,'p',                                   &
       &   'delta14c','delta14C of DIC',' ','permil',1)
  call ncdefvar3d(LYR_BIGD14C(iogrp),cmpflg,'p',                                &
       &   'bigdelta14c','big delta14C of DIC',' ','permil',1)
  call ncdefvar3d(LYR_POC13(iogrp),cmpflg,'p',                                  &
       &   'detoc13','Detritus13',' ','mol P m-3',1)
  call ncdefvar3d(LYR_DOC13(iogrp),cmpflg,'p',                                  &
       &   'dissoc13','Dissolved organic carbon13',' ','mol P m-3',1)
  call ncdefvar3d(LYR_CALC13(iogrp),cmpflg,'p',                                 &
       &   'calc13','Ca13CO3 shells',' ','mol 13C m-3',1)
  call ncdefvar3d(LYR_PHYTO13(iogrp),cmpflg,'p',                                &
       &   'phyc13','Phytoplankton13',' ','mol P m-3',1)
  call ncdefvar3d(LYR_GRAZER13(iogrp),cmpflg,'p',                               &
       &   'zooc13','Zooplankton13',' ','mol P m-3',1)
#endif
#ifdef AGG
  call ncdefvar3d(LYR_NOS(iogrp),cmpflg,'p',                                    &
       &   'nos','Marine snow aggregates per cm^3 sea water',' ','1/cm^3',1)
  call ncdefvar3d(LYR_WPHY(iogrp),cmpflg,'p',                                   &
       &   'wphy','Av. mass sinking speed of marine snow',' ','m/day',1)
  call ncdefvar3d(LYR_WNOS(iogrp),cmpflg,'p',                                   &
       &   'wnos','Av. number sinking speed of marine snow',' ','m/day',1)
  call ncdefvar3d(LYR_EPS(iogrp),cmpflg,'p',                                    &
       &   'eps','Av. size distribution exponent',' ','-',1)
  call ncdefvar3d(LYR_ASIZE(iogrp),cmpflg,'p',                                  &
       &   'asize','Av. size of marine snow aggregates',' ','nb. of cells',1)
#endif
#ifdef CFC
  call ncdefvar3d(LYR_CFC11(iogrp),cmpflg,'p',                                  &
       &   'cfc11','CFC-11',' ','mol cfc11 m-3',1)
  call ncdefvar3d(LYR_CFC12(iogrp),cmpflg,'p',                                  &
       &   'cfc12','CFC-12',' ','mol cfc12 m-3',1)
  call ncdefvar3d(LYR_SF6(iogrp),cmpflg,'p',                                    &
       &   'sf6','SF-6',' ','mol sf6 m-3',1)
#endif
#ifdef natDIC
  call ncdefvar3d(LYR_NATCO3(iogrp),cmpflg,'p',                                 &
       &   'natco3','Natural Carbonate ions',' ','mol C m-3',1)
  call ncdefvar3d(LYR_NATALKALI(iogrp),cmpflg,'p','nattalk',                    &
       &   'Natural alkalinity',' ','eq m-3',1)
  call ncdefvar3d(LYR_NATDIC(iogrp),cmpflg,'p','natdissic',                     &
       &   'Natural dissolved inorganic carbon',' ','mol C m-3',1)
  call ncdefvar3d(LYR_NATCALC(iogrp),cmpflg,'p','natcalc',                      &
       &   'Natural CaCO3',' ','mol C m-3',1)
  call ncdefvar3d(LYR_NATPH(iogrp),cmpflg,'p',                                  &
       &   'natph','Natural pH',' ','-log10([h+])',1)
  call ncdefvar3d(LYR_NATOMEGAA(iogrp),cmpflg,'p','natomegaa',                  &
       &   'Natural OmegaA',' ','1',1)
  call ncdefvar3d(LYR_NATOMEGAC(iogrp),cmpflg,'p','natomegac',                  &
       &   'Natural OmegaC',' ','1',1)
#endif
#ifdef BROMO
  call ncdefvar3d(LYR_BROMO(iogrp),cmpflg,'p',                                  &
       &   'bromo','Bromoform',' ','mol CHBr3 m-3',1)
#endif
#ifdef extNcycle
  call ncdefvar3d(LYR_ANH4(iogrp),cmpflg,'p',                                   &
       &  'nh4','Ammonium',' ','mol N m-3',1)
  call ncdefvar3d(LYR_ANO2(iogrp),cmpflg,'p',                                   &
       &  'no2','Nitrite',' ','mol N m-3',1)
  call ncdefvar3d(LYR_nitr_NH4(iogrp),cmpflg,'p',                               &
       &  'nh4nitr','NH4 nitrification rate',' ','mol N m-3 s-1',1)
  call ncdefvar3d(LYR_nitr_NO2(iogrp),cmpflg,'p',                               &
       &  'no2nitr','NO2 nitrification rate',' ','mol N m-3 s-1',1)
  call ncdefvar3d(LYR_nitr_N2O_prod(iogrp),cmpflg,'p',                          &
       &  'nitr_n2o','N2O prod during NH4 nitrification',' ',                   &
       &  'mol N2O m-3 s-1',1)
  call ncdefvar3d(LYR_nitr_NH4_OM(iogrp),cmpflg,'p',                            &
       &  'nh4nitr_om','OM production during NH4 nitrification',' ',            &
       &  'mol P m-3 s-1',1)
  call ncdefvar3d(LYR_nitr_NO2_OM(iogrp),cmpflg,'p',                            &
       &  'no2nitr_om','OM production during NO2 nitrification',' ',            &
       &  'mol P m-3 s-1',1)
  call ncdefvar3d(LYR_denit_NO3(iogrp),cmpflg,'p',                              &
       &  'no3denit','NO3 denitrification rate',' ','mol N m-3 s-1',1)
  call ncdefvar3d(LYR_denit_NO2(iogrp),cmpflg,'p',                              &
       &  'no2denit','NO2 denitrification rate',' ','mol N m-3 s-1',1)
  call ncdefvar3d(LYR_denit_N2O(iogrp),cmpflg,'p',                              &
       &  'n2odenit','N2O denitrification rate',' ','mol N2O m-3 s-1',1)
  call ncdefvar3d(LYR_DNRA_NO2(iogrp),cmpflg,'p',                               &
       &  'no2dnra','NO2 DNRA rate',' ','mol N m-3 s-1',1)
  call ncdefvar3d(LYR_anmx_N2_prod(iogrp),cmpflg,'p',                           &
       &  'anmx_n2','Anammox N2 production rate',' ','mol N2 m-3 s-1',1)
  call ncdefvar3d(LYR_anmx_OM_prod(iogrp),cmpflg,'p',                           &
       &  'anmx_om','Anammox OM production rate',' ','mol P m-3 s-1',1)
  call ncdefvar3d(LYR_phosy_NH4(iogrp),cmpflg,'p',                              &
       &  'phosy_nh4','PP consumption rate of NH4',' ','mol N m-3 s-1',1)
  call ncdefvar3d(LYR_phosy_NO3(iogrp),cmpflg,'p',                              &
       &  'phosy_no3','PP consumption rate of NO3',' ','mol N m-3 s-1',1)
  call ncdefvar3d(LYR_remin_aerob(iogrp),cmpflg,'p',                            &
       &  'remina','Aerob remineralization rate',' ','mol N m-3 s-1',1)
  call ncdefvar3d(LYR_remin_sulf(iogrp),cmpflg,'p',                             &
       &  'remins','Sulfate remineralization rate',' ','mol P m-3 s-1',1)
#endif
  !      M4AGO
  call ncdefvar3d(LYR_agg_ws(iogrp),cmpflg,'p',                                 &
       &  'agg_ws','aggregate mean settling velocity',' ','m d-1',1)
  call ncdefvar3d(LYR_dynvis(iogrp),cmpflg,'p',                                 &
       &  'dynvis','dynamic viscosity of sea water',' ','kg m-1 s-1',1)
  call ncdefvar3d(LYR_agg_stick(iogrp),cmpflg,'p',                              &
       &  'agg_stick','aggregate mean stickiness',' ','-',1)
  call ncdefvar3d(LYR_agg_stickf(iogrp),cmpflg,'p',                             &
       &  'agg_stickf','opal frustule stickiness',' ','-',1)
  call ncdefvar3d(LYR_agg_dmax(iogrp),cmpflg,'p',                               &
       &  'agg_dmax','aggregate maximum diameter',' ','m',1)
  call ncdefvar3d(LYR_agg_avdp(iogrp),cmpflg,'p',                               &
       &  'agg_avdp','mean primary particle diameter',' ','m',1)
  call ncdefvar3d(LYR_agg_avrhop(iogrp),cmpflg,'p',                             &
       &  'agg_avrhop','mean primary particle density',' ','kg m-3',1)
  call ncdefvar3d(LYR_agg_avdC(iogrp),cmpflg,'p',                               &
       &  'agg_avdC','Conc.-weighted mean aggregate diameter',' ','m',1)
  call ncdefvar3d(LYR_agg_df(iogrp),cmpflg,'p',                                 &
       &  'agg_df','aggregate fractal dimension',' ','-',1)
  call ncdefvar3d(LYR_agg_b(iogrp),cmpflg,'p',                                  &
       &  'agg_b','aggregate number distribution slope',' ','-',1)
  call ncdefvar3d(LYR_agg_Vrhof(iogrp),cmpflg,'p',                              &
       &  'agg_Vrhof','V-weighted aggregate mean density',' ','kg m-3',1)
  call ncdefvar3d(LYR_agg_Vpor(iogrp),cmpflg,'p',                               &
       &  'agg_Vpor','V-weighted aggregate mean porosity',' ','-',1)

  ! --- define 3d level fields
  call ncdefvar3d(LVL_DIC(iogrp),cmpflg,'p',                                    &
       &   'dissiclvl','Dissolved inorganic carbon',' ','mol C m-3',2)
  call ncdefvar3d(LVL_ALKALI(iogrp),cmpflg,'p',                                 &
       &   'talklvl','Alkalinity',' ','eq m-3',2)
  call ncdefvar3d(LVL_PHOSPH(iogrp),cmpflg,'p',                                 &
       &   'po4lvl','Phosphorus',' ','mol P m-3',2)
  call ncdefvar3d(LVL_OXYGEN(iogrp),cmpflg,'p',                                 &
       &   'o2lvl','Oxygen',' ','mol O2 m-3',2)
  call ncdefvar3d(LVL_ANO3(iogrp),cmpflg,'p',                                   &
       &   'no3lvl','Nitrate',' ','mol N m-3',2)
  call ncdefvar3d(LVL_SILICA(iogrp),cmpflg,'p',                                 &
       &   'silvl','Silicate',' ','mol Si m-3',2)
  call ncdefvar3d(LVL_DOC(iogrp),cmpflg,'p',                                    &
       &   'dissoclvl','Dissolved organic carbon',' ','mol P m-3',2)
  call ncdefvar3d(LVL_PHYTO(iogrp),cmpflg,'p',                                  &
       &   'phyclvl','Phytoplankton',' ','mol P m-3',2)
  call ncdefvar3d(LVL_GRAZER(iogrp),cmpflg,'p',                                 &
       &   'zooclvl','Zooplankton',' ','mol P m-3',2)
  call ncdefvar3d(LVL_POC(iogrp),cmpflg,'p',                                    &
       &   'detoclvl','Detritus',' ','mol P m-3',2)
  call ncdefvar3d(LVL_CALC(iogrp),cmpflg,'p',                                   &
       &   'calclvl','CaCO3 shells',' ','mol C m-3',2)
  call ncdefvar3d(LVL_OPAL(iogrp),cmpflg,'p',                                   &
       &   'opallvl','Opal shells',' ','mol Si m-3',2)
  call ncdefvar3d(LVL_IRON(iogrp),cmpflg,'p',                                   &
       &   'dfelvl','Dissolved iron',' ','mol Fe m-3',2)
  call ncdefvar3d(LVL_PHOSY(iogrp),cmpflg,'p',                                  &
       &   'pplvl','Primary production',' ','mol C m-3 s-1',2)
  call ncdefvar3d(LVL_CO3(iogrp),cmpflg,'p',                                    &
       &   'co3lvl','Carbonate ions',' ','mol C m-3',2)
  call ncdefvar3d(LVL_PH(iogrp),cmpflg,'p',                                     &
       &   'phlvl','pH',' ','-log10([h+])',2)
  call ncdefvar3d(LVL_OMEGAA(iogrp),cmpflg,'p',                                 &
       &   'omegaalvl','OmegaA',' ','1',2)
  call ncdefvar3d(LVL_OMEGAC(iogrp),cmpflg,'p',                                 &
       &   'omegaclvl','OmegaC',' ','1',2)
  call ncdefvar3d(LVL_N2O(iogrp),cmpflg,'p',                                    &
       &   'n2olvl','N2O',' ','mol N2O m-3',2)
  call ncdefvar3d(LVL_PREFO2(iogrp),cmpflg,'p',                                 &
       &   'p_o2lvl','Preformed oxygen',' ','mol O2 m-3',2)
  call ncdefvar3d(LVL_O2SAT(iogrp),cmpflg,'p',                                  &
       &   'satoxylvl','Saturated oxygen',' ','mol O2 m-3',2)
  call ncdefvar3d(LVL_PREFPO4(iogrp),cmpflg,'p',                                &
       &   'p_po4lvl','Preformed phosphorus',' ','mol P m-3',2)
  call ncdefvar3d(LVL_PREFALK(iogrp),cmpflg,'p',                                &
       &   'p_talklvl','Preformed alkalinity',' ','eq m-3',2)
  call ncdefvar3d(LVL_PREFDIC(iogrp),cmpflg,'p',                                &
       &   'p_diclvl','Preformed DIC',' ','mol C m-3',2)
  call ncdefvar3d(LVL_DICSAT(iogrp),cmpflg,'p',                                 &
       &   'sat_diclvl','Saturated DIC',' ','mol C m-3',2)
#ifdef cisonew
  call ncdefvar3d(LVL_DIC13(iogrp),cmpflg,'p',                                  &
       &   'dissic13lvl','Dissolved C13',' ','mol 13C m-3',2)
  call ncdefvar3d(LVL_DIC14(iogrp),cmpflg,'p',                                  &
       &   'dissic14lvl','Dissolved C14',' ','mol 14C m-3',2)
  call ncdefvar3d(LVL_D13C(iogrp),cmpflg,'p',                                   &
       &   'delta13clvl','delta13C of DIC',' ','permil',2)
  call ncdefvar3d(LVL_D14C(iogrp),cmpflg,'p',                                   &
       &   'delta14clvl','delta14C of DIC',' ','permil',2)
  call ncdefvar3d(LVL_BIGD14C(iogrp),cmpflg,'p',                                &
       &   'bigdelta14clvl','big delta14C of DIC',' ','permil',2)
  call ncdefvar3d(LVL_POC13(iogrp),cmpflg,'p',                                  &
       &   'detoc13lvl','Detritus13',' ','mol P m-3',2)
  call ncdefvar3d(LVL_DOC13(iogrp),cmpflg,'p',                                  &
       &   'dissoc13lvl','Dissolved organic carbon13',' ','mol P m-3',2)
  call ncdefvar3d(LVL_CALC13(iogrp),cmpflg,'p',                                 &
       &   'calc13lvl','Ca13CO3 shells',' ','mol 13C m-3',2)
  call ncdefvar3d(LVL_PHYTO13(iogrp),cmpflg,'p',                                &
       &   'phyc13lvl','Phytoplankton13',' ','mol P m-3',2)
  call ncdefvar3d(LVL_GRAZER13(iogrp),cmpflg,'p',                               &
       &   'zooc13lvl','Zooplankton13',' ','mol P m-3',2)
#endif
#ifdef AGG
  call ncdefvar3d(LVL_NOS(iogrp),cmpflg,'p','noslvl',                           &
       &   'Marine snow aggregates per cm^3 sea water',' ','1/cm^3',2)
  call ncdefvar3d(LVL_WPHY(iogrp),cmpflg,'p','wphylvl',                         &
       &   'Av. mass sinking speed of marine snow',' ','m/day',2)
  call ncdefvar3d(LVL_WNOS(iogrp),cmpflg,'p','wnoslvl',                         &
       &   'Av. number sinking speed of marine snow',' ','m/day',2)
  call ncdefvar3d(LVL_EPS(iogrp),cmpflg,'p','epslvl',                           &
       &   'Av. size distribution exponent',' ','-',2)
  call ncdefvar3d(LVL_ASIZE(iogrp),cmpflg,'p','asizelvl',                       &
       &   'Av. size of marine snow aggregates',' ','nb. of cells',2)
#endif
#ifdef CFC
  call ncdefvar3d(LVL_CFC11(iogrp),cmpflg,'p',                                  &
       &   'cfc11lvl','CFC-11',' ','mol cfc11 m-3',2)
  call ncdefvar3d(LVL_CFC12(iogrp),cmpflg,'p',                                  &
       &   'cfc12lvl','CFC-12',' ','mol cfc12 m-3',2)
  call ncdefvar3d(LVL_SF6(iogrp),cmpflg,'p',                                    &
       &   'sf6lvl','SF-6',' ','mol sf6 m-3',2)
#endif
#ifdef natDIC
  call ncdefvar3d(LVL_NATCO3(iogrp),cmpflg,'p',                                 &
       &   'natco3lvl','Natural Carbonate ions',' ','mol C m-3',2)
  call ncdefvar3d(LVL_NATALKALI(iogrp),cmpflg,'p','nattalklvl',                 &
       &   'Natural alkalinity',' ','eq m-3',2)
  call ncdefvar3d(LVL_NATDIC(iogrp),cmpflg,'p','natdissiclvl',                  &
       &   'Natual dissolved inorganic carbon',' ','mol C m-3',2)
  call ncdefvar3d(LVL_NATCALC(iogrp),cmpflg,'p',                                &
       &   'natcalclvl','Natural CaCO3 shells',' ','mol C m-3',2)
  call ncdefvar3d(LVL_NATPH(iogrp),cmpflg,'p',                                  &
       &   'natphlvl','Natural pH',' ','-log10([h+])',2)
  call ncdefvar3d(LVL_NATOMEGAA(iogrp),cmpflg,'p',                              &
       &   'natomegaalvl','Natural OmegaA',' ','1',2)
  call ncdefvar3d(LVL_NATOMEGAC(iogrp),cmpflg,'p',                              &
       &   'natomegaclvl','Natural OmegaC',' ','1',2)
#endif
#ifdef BROMO
  call ncdefvar3d(LVL_BROMO(iogrp),cmpflg,'p',                                  &
       &   'bromolvl','Bromoform',' ','mol CHBr3 m-3',2)
#endif
#ifdef extNcycle
  call ncdefvar3d(LVL_ANH4(iogrp),cmpflg,'p',                                   &
       &  'nh4lvl','Ammonium',' ','mol N m-3',2)
  call ncdefvar3d(LVL_ANO2(iogrp),cmpflg,'p',                                   &
       &  'no2lvl','Nitrite',' ','mol N m-3',2)
  call ncdefvar3d(LVL_nitr_NH4(iogrp),cmpflg,'p',                               &
       &  'nh4nitrlvl','NH4 nitrification rate',' ','mol N m-3 s-1',2)
  call ncdefvar3d(LVL_nitr_NO2(iogrp),cmpflg,'p',                               &
       &  'no2nitrlvl','NO2 nitrification rate',' ','mol N m-3 s-1',2)
  call ncdefvar3d(LVL_nitr_N2O_prod(iogrp),cmpflg,'p',                          &
       &  'nitr_n2olvl','N2O prod during NH4 nitrification',' ',                &
       &  'mol N2O m-3 s-1',2)
  call ncdefvar3d(LVL_nitr_NH4_OM(iogrp),cmpflg,'p',                            &
       &  'nh4nitr_omlvl','OM production during NH4 nitrification',' ',         &
       &  'mol P m-3 s-1',2)
  call ncdefvar3d(LVL_nitr_NO2_OM(iogrp),cmpflg,'p',                            &
       &  'no2nitr_omlvl','OM production during NO2 nitrification',' ',         &
       &  'mol P m-3 s-1',2)
  call ncdefvar3d(LVL_denit_NO3(iogrp),cmpflg,'p',                              &
       &  'no3denitlvl','NO3 denitrification rate',' ','mol N m-3 s-1',2)
  call ncdefvar3d(LVL_denit_NO2(iogrp),cmpflg,'p',                              &
       &  'no2denitlvl','NO2 denitrification rate',' ','mol N m-3 s-1',2)
  call ncdefvar3d(LVL_denit_N2O(iogrp),cmpflg,'p',                              &
       &  'n2odenitlvl','N2O denitrification rate',' ',                         &
       &  'mol N2O m-3 s-1',2)
  call ncdefvar3d(LVL_DNRA_NO2(iogrp),cmpflg,'p',                               &
       &  'no2dnralvl','NO2 DNRA rate',' ','mol N m-3 s-1',2)
  call ncdefvar3d(LVL_anmx_N2_prod(iogrp),cmpflg,'p',                           &
       &  'anmx_n2lvl','Anammox N2 production rate',' ',                        &
       &  'mol N2 m-3 s-1',2)
  call ncdefvar3d(LVL_anmx_OM_prod(iogrp),cmpflg,'p',                           &
       &  'anmx_omlvl','Anammox OM production rate',' ','mol P m-3 s-1',2)
  call ncdefvar3d(LVL_phosy_NH4(iogrp),cmpflg,'p',                              &
       &  'phosy_nh4lvl','PP consumption rate of NH4',' ',                      &
       &  'mol N m-3 s-1',2)
  call ncdefvar3d(LVL_phosy_NO3(iogrp),cmpflg,'p',                              &
       &  'phosy_no3lvl','PP consumption rate of NO3',' ',                      &
       &  'mol N m-3 s-1',2)
  call ncdefvar3d(LVL_remin_aerob(iogrp),cmpflg,'p',                            &
       &  'reminalvl','Aerob remineralization rate',' ',                        &
       &  'mol N m-3 s-1',2)
  call ncdefvar3d(LVL_remin_sulf(iogrp),cmpflg,'p',                             &
       &  'reminslvl','Sulfate remineralization rate',' ',                      &
       &  'mol P m-3 s-1',2)
#endif
  !      M4AGO
  call ncdefvar3d(LVL_agg_ws(iogrp),cmpflg,'p',                                 &
       &  'agg_wslvl','aggregate mean settling velocity',' ','m d-1',2)
  call ncdefvar3d(LVL_dynvis(iogrp),cmpflg,'p',                                 &
       &  'dynvislvl','dynamic viscosity of sea water',' ','kg m-1 s-1',        &
       &   2)
  call ncdefvar3d(LVL_agg_stick(iogrp),cmpflg,'p',                              &
       &  'agg_sticklvl','aggregate mean stickiness',' ','-',2)
  call ncdefvar3d(LVL_agg_stickf(iogrp),cmpflg,'p',                             &
       &  'agg_stickflvl','opal frustule stickiness',' ','-',2)
  call ncdefvar3d(LVL_agg_dmax(iogrp),cmpflg,'p',                               &
       &  'agg_dmaxlvl','aggregate maximum diameter',' ','m',2)
  call ncdefvar3d(LVL_agg_avdp(iogrp),cmpflg,'p',                               &
       &  'agg_avdplvl','mean primary particle diameter',' ','m',2)
  call ncdefvar3d(LVL_agg_avrhop(iogrp),cmpflg,'p',                             &
       &  'agg_avrhoplvl','mean primary particle density',' ','kg m-3',2)
  call ncdefvar3d(LVL_agg_avdC(iogrp),cmpflg,'p',                               &
       &  'agg_avdClvl','Conc.-weighted mean aggregate diameter',' ',           &
       &  'm',2)
  call ncdefvar3d(LVL_agg_df(iogrp),cmpflg,'p',                                 &
       &  'agg_dflvl','aggregate fractal dimension',' ','-',2)
  call ncdefvar3d(LVL_agg_b(iogrp),cmpflg,'p',                                  &
       &  'agg_blvl','aggregate number distribution slope',' ','-',2)
  call ncdefvar3d(LVL_agg_Vrhof(iogrp),cmpflg,'p',                              &
       &  'agg_Vrhoflvl','V-weighted aggregate mean density',' ',               &
       &  'kg m-3',2)
  call ncdefvar3d(LVL_agg_Vpor(iogrp),cmpflg,'p',                               &
       &  'agg_Vporlvl','V-weighted aggregate mean porosity',' ','-',2)

  ! --- define sediment fields
#ifndef sedbypass
  call ncdefvar3d(SDM_POWAIC(iogrp),cmpflg,'p',                                 &
       &   'powdic','PoWa DIC',' ','mol C m-3',3)
  call ncdefvar3d(SDM_POWAAL(iogrp),cmpflg,'p',                                 &
       &   'powalk','PoWa alkalinity',' ','eq m-3',3)
  call ncdefvar3d(SDM_POWAPH(iogrp),cmpflg,'p',                                 &
       &   'powpho','PoWa phosphorus',' ','mol P m-3',3)
  call ncdefvar3d(SDM_POWAOX(iogrp),cmpflg,'p',                                 &
       &   'powox','PoWa oxygen',' ','mol O2 m-3',3)
  call ncdefvar3d(SDM_POWN2(iogrp), cmpflg,'p',                                 &
       &   'pown2','PoWa N2',' ','mol N2 m-3',3)
  call ncdefvar3d(SDM_POWNO3(iogrp),cmpflg,'p',                                 &
       &   'powno3','PoWa nitrate',' ','mol N m-3',3)
  call ncdefvar3d(SDM_POWASI(iogrp),cmpflg,'p',                                 &
       &   'powsi','PoWa silicate',' ','mol Si m-3',3)
  call ncdefvar3d(SDM_SSSO12(iogrp),cmpflg,'p',                                 &
       &   'ssso12','Sediment detritus',' ','mol P m-3',3)
  call ncdefvar3d(SDM_SSSSIL(iogrp),cmpflg,'p',                                 &
       &   'ssssil','Sediment silicate',' ','mol Si m-3',3)
  call ncdefvar3d(SDM_SSSC12(iogrp),cmpflg,'p',                                 &
       &   'sssc12','Sediment CaCO3',' ','mol C m-3',3)
  call ncdefvar3d(SDM_SSSTER(iogrp),cmpflg,'p',                                 &
       &   'ssster','Sediment clay',' ','mol m-3',3)

  ! --- define sediment burial fields
  call ncdefvar3d(BUR_SSSO12(iogrp),                                            &
       &   cmpflg,'p','buro12','Burial org carbon',' ','mol P m-2',4)
  call ncdefvar3d(BUR_SSSC12(iogrp),                                            &
       &   cmpflg,'p','burc12','Burial calcium ',' ','mol C m-2',4)
  call ncdefvar3d(BUR_SSSSIL(iogrp),                                            &
       &   cmpflg,'p','bursil','Burial silicate',' ','mol Si m-2',4)
  call ncdefvar3d(BUR_SSSTER(iogrp),                                            &
       &   cmpflg,'p','burter','Burial clay',' ','mol  m-2',4)
#endif
#if defined(extNcycle) && ! defined(sedbypass)
  call ncdefvar3d(SDM_POWNH4(iogrp),cmpflg,'p',                                 &
       &   'pownh4','PoWa ammonium',' ','mol N m-3',3)
  call ncdefvar3d(SDM_POWN2O(iogrp),cmpflg,'p',                                 &
       &   'pown2o','PoWa nitrous oxide',' ','mol N m-3',3)
  call ncdefvar3d(SDM_POWNO2(iogrp),cmpflg,'p',                                 &
       &   'powno2','PoWa nitrite',' ','mol N m-3',3)

  call ncdefvar3d(sdm_nitr_NH4(iogrp),cmpflg,'p',                               &
     &  'nh4nitrsdm','NH4 nitrification rate sediment',' ','mol N m-3 s-1',3)
  call ncdefvar3d(sdm_nitr_NO2(iogrp),cmpflg,'p',                               &
     &  'no2nitrsdm','NO2 nitrification rate sediment',' ','mol N m-3 s-1',3)
  call ncdefvar3d(sdm_nitr_N2O_prod(iogrp),cmpflg,'p',                          &
     &  'nitr_n2osdm','N2O prod during NH4 nitrification sediment',' ',         &
     &  'mol N2O m-3 s-1',3)
  call ncdefvar3d(sdm_nitr_NH4_OM(iogrp),cmpflg,'p',                            &
     &  'nh4nitr_omsdm','OM production during NH4 nitrification sediment',' ',  &
     &  'mol P m-3 s-1',3)
  call ncdefvar3d(sdm_nitr_NO2_OM(iogrp),cmpflg,'p',                            &
     &  'no2nitr_omsdm','OM production during NO2 nitrification sediment',' ',  &
     &  'mol P m-3 s-1',3)
  call ncdefvar3d(sdm_denit_NO3(iogrp),cmpflg,'p',                              &
     &  'no3denitsdm','NO3 denitrification rate sediment',' ','mol N m-3 s-1',3)
  call ncdefvar3d(sdm_denit_NO2(iogrp),cmpflg,'p',                              &
     &  'no2denitsdm','NO2 denitrification rate sediment',' ','mol N m-3 s-1',3)
  call ncdefvar3d(sdm_denit_N2O(iogrp),cmpflg,'p',                              &
     &  'n2odenitsdm','N2O denitrification rate sediment',' ',                  &
     &  'mol N2O m-3 s-1',3)
  call ncdefvar3d(sdm_DNRA_NO2(iogrp),cmpflg,'p',                               &
     &  'no2dnrasdm','NO2 DNRA rate sediment',' ','mol N m-3 s-1',3)
  call ncdefvar3d(sdm_anmx_N2_prod(iogrp),cmpflg,'p',                           &
     &  'anmx_n2sdm','Anammox N2 production rate sediment',' ',                 &
     &  'mol N2 m-3 s-1',3)
  call ncdefvar3d(sdm_anmx_OM_prod(iogrp),cmpflg,'p',                           &
     &  'anmx_omsdm','Anammox OM production rate sediment',' ',                 &
     &  'mol P m-3 s-1',3)
  call ncdefvar3d(sdm_remin_aerob(iogrp),cmpflg,'p',                            &
     &  'reminasdm','Aerob remineralization rate sediment',' ',                 &
     &  'mol N m-3 s-1',3)
  call ncdefvar3d(sdm_remin_sulf(iogrp),cmpflg,'p',                             &
     &  'reminssdm','Sulfate remineralization rate sediment',' ',               &
     &  'mol P m-3 s-1',3) 
#endif 
  ! --- enddef netcdf file
  call ncedef
end subroutine hamoccvardef
