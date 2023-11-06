! Copyright (C) 2002  P. Wetzel
! Copyright (C) 2020  I. Bethke, J. Tjiputra, J. Schwinger, A. Moree,
!                     P.-G. Chiu, M. Bentsen
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


      MODULE mo_bgcmean
!***********************************************************************
!
!**** *MODULE mo_bgcmean* - Variables for bgcmean.
!
!
!     Patrick Wetzel    *MPI-Met, HH*    09.12.02
!     Ingo Bethke       *Bjer.NE. C.*    05.11.09 
!     J. Schwinger      *GFI, UiB        10.02.12
!      - added variables and functions for sediment burial
!      - added variables for CFC output
!      - added initialisation of namelist variables and 
!        index arrays
!
!     Tjiputra          *UNI-RESEARCH    25.11.15
!      - added natural DIC/ALK/CALC/OMEGAC variables 
!  
!     A.Moree,          *GFI, Bergen*   2018-04-12
!     - new version of carbon isotope code
!
!     J.Tjiputra,       *Uni Research, Bergen*   2018-04-12
!     - added preformed and saturated DIC tracers
!
!     J.Schwinger,      *Uni Research, Bergen*   2018-04-12
!     - changed naming of particle fluxes
!     - removed output of AOU and added O2_sat instead
!     - added output of omegaA
!     - added sediment bypass preprocessor option
!
!     Purpose
!     -------
!     - declaration and memory allocation
!     - declaration of auxiliary functions  
!
!**********************************************************************
      use mod_xc,         only: ii,jj,kk,idm,jdm,kdm,nbdy,ifp,isp,ilp,mnproc,ip
      use mod_dia,        only: ddm,depthslev,depthslev_bnds,nstepinday,pbath
      use mod_nctools,    only: ncpack,nccomp,nccopa,ncwrtr
      use netcdf,         only: nf90_fill_double
      use mo_param1_bgc,  only: ks
      use mo_control_bgc, only: use_sedbypass,use_cisonew,use_CFC,use_natDIC,use_BROMO,use_BOXATM,use_AGG

      implicit none

      PRIVATE :: ii,jj,kk,idm,jdm,kdm,nbdy,ifp,isp,ilp                
      PUBLIC  :: ks,ddm,depthslev,depthslev_bnds

! --- Averaging and writing frequencies for diagnostic output     
      integer, save :: nbgc
      integer, parameter :: nbgcmax=10
      real, dimension(nbgcmax), save ::  diagfq_bgc,filefq_bgc
      integer, dimension(nbgcmax), save :: nacc_bgc
      logical, dimension(nbgcmax), save :: diagmon_bgc,diagann_bgc,     &
     &  filemon_bgc,fileann_bgc,bgcwrt
 
! --- Namelist for diagnostic output 
      integer, dimension(nbgcmax), save ::                              &
     & SRF_KWCO2     =0    ,SRF_PCO2      =0    ,SRF_DMSFLUX   =0    ,  &
     & SRF_KWCO2KHM  =0    ,SRF_CO2KHM    =0    ,SRF_CO2KH     =0    ,  &
     & SRF_PCO2M     =0    ,                                            &
     & SRF_CO2FXD    =0    ,SRF_CO2FXU    =0    ,SRF_CO213FXD  =0    ,  &
     & SRF_CO213FXU  =0    ,SRF_CO214FXD  =0    ,SRF_CO214FXU  =0    ,  &
     & SRF_OXFLUX    =0    ,SRF_NIFLUX    =0    ,SRF_DMS       =0    ,  &
     & SRF_DMSPROD   =0    ,SRF_DMS_BAC   =0    ,SRF_DMS_UV    =0    ,  &
     & SRF_EXPORT    =0    ,SRF_EXPOSI    =0    ,SRF_EXPOCA    =0    ,  &
     & SRF_ATMCO2    =0    ,SRF_ATMO2     =0    ,SRF_ATMN2     =0    ,  &
     & SRF_ATMC13    =0    ,SRF_ATMC14    =0    ,                       &
     & SRF_N2OFX     =0    ,SRF_CFC11     =0    ,SRF_CFC12     =0    ,  &
     & SRF_SF6       =0    ,SRF_PHOSPH    =0    ,SRF_OXYGEN    =0    ,  &
     & SRF_IRON      =0    ,SRF_ANO3      =0    ,SRF_ALKALI    =0    ,  &
     & SRF_SILICA    =0    ,SRF_DIC       =0    ,SRF_PHYTO     =0    ,  &
     & SRF_PH        =0    ,                                            &
     & SRF_NATDIC    =0    ,SRF_NATALKALI =0    ,SRF_NATPCO2   =0    ,  &
     & SRF_NATCO2FX  =0    ,SRF_NATPH     =0    ,                       &
     & SRF_ATMBROMO  =0    ,SRF_BROMO     =0    ,SRF_BROMOFX   =0    ,  &
     & INT_BROMOPRO  =0    ,INT_BROMOUV   =0    ,                       &
     & INT_PHOSY     =0    ,INT_NFIX      =0    ,INT_DNIT      =0    ,  &
     & FLX_NDEP      =0    ,FLX_OALK      =0    ,                       &
     & FLX_CAR0100   =0    ,FLX_CAR0500   =0    ,FLX_CAR1000   =0    ,  &
     & FLX_CAR2000   =0    ,FLX_CAR4000   =0    ,FLX_CAR_BOT   =0    ,  &
     & FLX_BSI0100   =0    ,FLX_BSI0500   =0    ,FLX_BSI1000   =0    ,  &
     & FLX_BSI2000   =0    ,FLX_BSI4000   =0    ,FLX_BSI_BOT   =0    ,  &
     & FLX_CAL0100   =0    ,FLX_CAL0500   =0    ,FLX_CAL1000   =0    ,  &
     & FLX_CAL2000   =0    ,FLX_CAL4000   =0    ,FLX_CAL_BOT   =0    ,  &
     & FLX_SEDIFFIC  =0    ,FLX_SEDIFFAL  =0    ,FLX_SEDIFFPH  =0    ,  &
     & FLX_SEDIFFOX  =0    ,FLX_SEDIFFN2  =0    ,FLX_SEDIFFNO3 =0    ,  &
     & FLX_SEDIFFSI  =0    ,                                            &   
     & LYR_PHYTO     =0    ,LYR_GRAZER    =0    ,LYR_DOC       =0    ,  &
     & LYR_PHOSY     =0    ,LYR_PHOSPH    =0    ,LYR_OXYGEN    =0    ,  &
     & LYR_IRON      =0    ,LYR_ANO3      =0    ,LYR_ALKALI    =0    ,  &
     & LYR_SILICA    =0    ,LYR_DIC       =0    ,LYR_POC       =0    ,  &
     & LYR_CALC      =0    ,LYR_OPAL      =0    ,LYR_CO3       =0    ,  &
     & LYR_PH        =0    ,LYR_OMEGAA    =0    ,LYR_OMEGAC    =0    ,  &
     & LYR_DIC13     =0    ,LYR_DIC14     =0    ,LYR_DP        =0    ,  &
     & LYR_NOS       =0    ,LYR_WPHY      =0    ,LYR_WNOS      =0    ,  &
     & LYR_EPS       =0    ,LYR_ASIZE     =0    ,LYR_N2O       =0    ,  &
     & LYR_PREFO2    =0    ,LYR_O2SAT     =0    ,LYR_PREFPO4   =0    ,  &
     & LYR_PREFALK   =0    ,LYR_PREFDIC   =0    ,LYR_DICSAT    =0    ,  &
     & LYR_CFC11     =0    ,LYR_CFC12     =0    ,LYR_SF6       =0    ,  &
     & LYR_NATDIC    =0    ,LYR_NATALKALI =0    ,LYR_NATCALC   =0    ,  &
     & LYR_NATPH     =0    ,LYR_NATOMEGAA =0    ,LYR_NATOMEGAC =0    ,  &
     & LYR_NATCO3    =0    ,                                            &
     & LYR_BROMO     =0    ,                                            &
     & LYR_D13C      =0    ,LYR_D14C      =0    ,LYR_BIGD14C   =0    ,  &
     & LYR_POC13     =0    ,LYR_DOC13     =0    ,LYR_CALC13    =0    ,  &
     & LYR_PHYTO13   =0    ,LYR_GRAZER13  =0    ,                       &
     & LVL_PHYTO     =0    ,LVL_GRAZER    =0    ,LVL_DOC       =0    ,  &
     & LVL_PHOSY     =0    ,LVL_PHOSPH    =0    ,LVL_OXYGEN    =0    ,  &
     & LVL_IRON      =0    ,LVL_ANO3      =0    ,LVL_ALKALI    =0    ,  &
     & LVL_SILICA    =0    ,LVL_DIC       =0    ,LVL_POC       =0    ,  &
     & LVL_CALC      =0    ,LVL_OPAL      =0    ,LVL_CO3       =0    ,  &
     & LVL_PH        =0    ,LVL_OMEGAA    =0    ,LVL_OMEGAC    =0    ,  &
     & LVL_DIC13     =0    ,LVL_DIC14     =0    ,LVL_NOS       =0    ,  &
     & LVL_WPHY      =0    ,LVL_WNOS      =0    ,LVL_EPS       =0    ,  &
     & LVL_ASIZE     =0    ,LVL_N2O       =0    ,LVL_PREFO2    =0    ,  &
     & LVL_O2SAT     =0    ,LVL_PREFPO4   =0    ,LVL_PREFALK   =0    ,  &
     & LVL_PREFDIC   =0    ,LVL_DICSAT    =0    ,                       &
     & LVL_CFC11     =0    ,LVL_CFC12     =0    ,LVL_SF6       =0    ,  &
     & LVL_NATDIC    =0    ,LVL_NATALKALI =0    ,LVL_NATCALC   =0    ,  &
     & LVL_NATPH     =0    ,LVL_NATOMEGAA =0    ,LVL_NATOMEGAC =0    ,  &
     & LVL_NATCO3    =0    ,                                            &
     & LVL_BROMO     =0    ,                                            &
     & LVL_D13C      =0    ,LVL_D14C      =0    ,LVL_BIGD14C   =0    ,  &
     & LVL_POC13     =0    ,LVL_DOC13     =0    ,LVL_CALC13    =0    ,  &
     & LVL_PHYTO13   =0    ,LVL_GRAZER13  =0    ,                       &
     & SDM_POWAIC    =0    ,SDM_POWAAL    =0    ,SDM_POWAPH    =0    ,  &
     & SDM_POWAOX    =0    ,SDM_POWN2     =0    ,SDM_POWNO3    =0    ,  &
     & SDM_POWASI    =0    ,SDM_SSSO12    =0    ,SDM_SSSSIL    =0    ,  &
     & SDM_SSSC12    =0    ,SDM_SSSTER    =0                         ,  &
     & BUR_SSSO12    =0    ,BUR_SSSC12    =0    ,BUR_SSSSIL    =0    ,  &
     & BUR_SSSTER    =0                                              ,  &
     & GLB_AVEPERIO  =0    ,GLB_FILEFREQ  =0    ,GLB_COMPFLAG  =0    ,  &
     & GLB_NCFORMAT  =0    ,GLB_INVENTORY =0 
      character(LEN=10), dimension(nbgcmax), save :: GLB_FNAMETAG
      namelist /DIABGC/                                                 &
     & SRF_KWCO2         ,SRF_PCO2          ,SRF_DMSFLUX       ,        &
     & SRF_KWCO2KHM      ,SRF_CO2KHM        ,SRF_CO2KH         ,        &
     & SRF_PCO2M         ,                                              &
     & SRF_CO2FXD        ,SRF_CO2FXU        ,SRF_CO213FXD      ,        &
     & SRF_CO213FXU      ,SRF_CO214FXD      ,SRF_CO214FXU      ,        &
     & SRF_OXFLUX        ,SRF_NIFLUX        ,SRF_DMS           ,        &
     & SRF_DMSPROD       ,SRF_DMS_BAC       ,SRF_DMS_UV        ,        &
     & SRF_EXPORT        ,SRF_EXPOSI        ,SRF_EXPOCA        ,        &
     & SRF_ATMCO2        ,SRF_ATMO2         ,SRF_ATMN2         ,        &
     & SRF_ATMC13        ,SRF_ATMC14        ,                           &
     & SRF_N2OFX         ,SRF_CFC11         ,SRF_CFC12         ,        &
     & SRF_SF6           ,SRF_PHOSPH        ,SRF_OXYGEN        ,        &
     & SRF_IRON          ,SRF_ANO3          ,SRF_ALKALI        ,        &
     & SRF_SILICA        ,SRF_DIC           ,SRF_PHYTO         ,        &
     & SRF_PH            ,                                              &
     & SRF_NATDIC        ,SRF_NATALKALI     ,SRF_NATPCO2       ,        &
     & SRF_NATCO2FX      ,SRF_NATPH         ,                           &
     & SRF_ATMBROMO      ,SRF_BROMO         ,SRF_BROMOFX       ,        &
     & INT_BROMOPRO      ,INT_BROMOUV       ,                           &
     & INT_PHOSY         ,INT_NFIX          ,INT_DNIT          ,        &
     & FLX_NDEP          ,FLX_OALK          ,                           &
     & FLX_CAR0100       ,FLX_CAR0500       ,FLX_CAR1000       ,        &
     & FLX_CAR2000       ,FLX_CAR4000       ,FLX_CAR_BOT       ,        &
     & FLX_BSI0100       ,FLX_BSI0500       ,FLX_BSI1000       ,        &
     & FLX_BSI2000       ,FLX_BSI4000       ,FLX_BSI_BOT       ,        &
     & FLX_CAL0100       ,FLX_CAL0500       ,FLX_CAL1000       ,        &
     & FLX_CAL2000       ,FLX_CAL4000       ,FLX_CAL_BOT       ,        &
     & FLX_SEDIFFIC      ,FLX_SEDIFFAL      ,FLX_SEDIFFPH      ,        &
     & FLX_SEDIFFOX      ,FLX_SEDIFFN2      ,FLX_SEDIFFNO3     ,        &
     & FLX_SEDIFFSI      ,                                              &   
     & LYR_PHYTO         ,LYR_GRAZER        ,LYR_DOC           ,        &
     & LYR_PHOSY         ,LYR_PHOSPH        ,LYR_OXYGEN        ,        &
     & LYR_IRON          ,LYR_ANO3          ,LYR_ALKALI        ,        &
     & LYR_SILICA        ,LYR_DIC           ,LYR_POC           ,        &
     & LYR_CALC          ,LYR_OPAL          ,LYR_CO3           ,        &
     & LYR_PH            ,LYR_OMEGAA        ,LYR_OMEGAC        ,        &
     & LYR_DIC13         ,LYR_DIC14         ,LYR_DP            ,        &
     & LYR_NOS           ,LYR_WPHY          ,LYR_WNOS          ,        &
     & LYR_EPS           ,LYR_ASIZE         ,LYR_N2O           ,        &
     & LYR_PREFO2        ,LYR_O2SAT         ,LYR_PREFPO4       ,        &
     & LYR_PREFALK       ,LYR_PREFDIC       ,LYR_DICSAT        ,        &
     & LYR_CFC11         ,LYR_CFC12         ,LYR_SF6           ,        &
     & LYR_NATDIC        ,LYR_NATALKALI     ,LYR_NATCALC       ,        &
     & LYR_NATPH         ,LYR_NATOMEGAA     ,LYR_NATOMEGAC     ,        &
     & LYR_NATCO3        ,                                              &
     & LYR_BROMO         ,                                              &
     & LYR_D13C          ,LYR_D14C          ,LYR_BIGD14C       ,        &
     & LYR_PHYTO13       ,LYR_GRAZER13      ,LYR_POC13         ,        &
     & LYR_DOC13         ,LYR_CALC13        ,                           &
     & LVL_PHYTO         ,LVL_GRAZER        ,LVL_DOC           ,        &
     & LVL_PHOSY         ,LVL_PHOSPH        ,LVL_OXYGEN        ,        &
     & LVL_IRON          ,LVL_ANO3          ,LVL_ALKALI        ,        &
     & LVL_SILICA        ,LVL_DIC           ,LVL_POC           ,        &
     & LVL_CALC          ,LVL_OPAL          ,LVL_CO3           ,        &
     & LVL_PH            ,LVL_OMEGAA        ,LVL_OMEGAC        ,        &
     & LVL_DIC13         ,LVL_DIC14         ,LVL_NOS           ,        &
     & LVL_WPHY          ,LVL_WNOS          ,LVL_EPS           ,        &
     & LVL_ASIZE         ,LVL_N2O           ,LVL_PREFO2        ,        &
     & LVL_O2SAT         ,LVL_PREFPO4       ,LVL_PREFALK       ,        &
     & LVL_PREFDIC       ,LVL_DICSAT        ,                           &
     & LVL_CFC11         ,LVL_CFC12         ,LVL_SF6           ,        &
     & LVL_NATDIC        ,LVL_NATALKALI     ,LVL_NATCALC       ,        &
     & LVL_NATPH         ,LVL_NATOMEGAA     ,LVL_NATOMEGAC     ,        &
     & LVL_NATCO3        ,                                              &
     & LVL_BROMO         ,                                              &
     & LVL_D13C          ,LVL_D14C          ,LVL_BIGD14C       ,        &
     & LVL_PHYTO13       ,LVL_GRAZER13      ,LVL_POC13         ,        &
     & LVL_DOC13         ,LVL_CALC13        ,                           &
     & SDM_POWAIC        ,SDM_POWAAL        ,SDM_POWAPH        ,        &
     & SDM_POWAOX        ,SDM_POWN2         ,SDM_POWNO3        ,        &
     & SDM_POWASI        ,SDM_SSSO12        ,SDM_SSSSIL        ,        &
     & SDM_SSSC12        ,SDM_SSSTER                           ,        &
     & BUR_SSSO12        ,BUR_SSSC12        ,BUR_SSSSIL        ,        &
     & BUR_SSSTER                                              ,        &
     & GLB_AVEPERIO      ,GLB_FILEFREQ      ,GLB_COMPFLAG      ,        &
     & GLB_NCFORMAT      ,GLB_FNAMETAG      ,GLB_INVENTORY

!----------------------------------------------------------------      
! declarations for inventory_bgc.F90
! order and increments of river (jir...) indices require to be the same 
! as in mo_riverinpt 
      integer, parameter ::                                             &
     &          jco2flux  =1,                                           &
     &          jo2flux   =2,                                           &
     &          jn2flux   =3,                                           &
     &          jn2oflux  =4,                                           &
     &          jprorca   =5,                                           &
     &          jprcaca   =6,                                           &
     &          jsilpro   =7,                                           &
     &          jpodiic   =8,                                           &
     &          jpodial   =9,                                           &
     &          jpodiph   =10,                                          &
     &          jpodiox   =11,                                          &
     &          jpodin2   =12,                                          &
     &          jpodino3  =13,                                          &
     &          jpodisi   =14,                                          &
     &          jndep     =15,                                          &
     &          joalk     =16,                                          &
     &          jirdin    =17,                                          &
     &          jirdip    =18,                                          &
     &          jirsi     =19,                                          &
     &          jiralk    =20,                                          &
     &          jiriron   =21,                                          &
     &          jirdoc    =22,                                          &
     &          jirdet    =23,                                          &
     &          nbgct2d   =23
      
!----------------------------------------------------------------      
      integer, save :: i_bsc_m2d 
      integer, dimension(nbgcmax), save ::                              &
     &          jkwco2     = 0 ,                                        &
     &          jkwco2khm  = 0 ,                                        &
     &          jco2kh     = 0 ,                                        &
     &          jco2khm    = 0 ,                                        &
     &          jpco2      = 0 ,                                        &
     &          jpco2m     = 0 ,                                        &
     &          jdmsflux   = 0 ,                                        &
     &          jco2fxd    = 0 ,                                        &
     &          jco2fxu    = 0 ,                                        &
     &          jco213fxd  = 0 ,                                        &
     &          jco213fxu  = 0 ,                                        &
     &          jco214fxd  = 0 ,                                        &
     &          jco214fxu  = 0 ,                                        &
     &          joxflux    = 0 ,                                        &
     &          jniflux    = 0 ,                                        &
     &          jn2ofx     = 0 ,                                        &
     &          jdms       = 0 ,                                        &
     &          jdmsprod   = 0 ,                                        &
     &          jdms_bac   = 0 ,                                        &
     &          jdms_uv    = 0 ,                                        &
     &          jexport    = 0 ,                                        &
     &          jexpoca    = 0 ,                                        &
     &          jexposi    = 0 ,                                        &
     &          jcfc11fx   = 0 ,                                        &
     &          jcfc12fx   = 0 ,                                        &
     &          jsf6fx     = 0 ,                                        &
     &          jsrfphosph = 0 ,                                        &
     &          jsrfoxygen = 0 ,                                        &
     &          jsrfiron   = 0 ,                                        &
     &          jsrfano3   = 0 ,                                        &
     &          jsrfalkali = 0 ,                                        &
     &          jsrfsilica = 0 ,                                        &
     &          jsrfdic    = 0 ,                                        &
     &          jsrfphyto  = 0 ,                                        &
     &          jsrfph     = 0 ,                                        &
     &          jintphosy  = 0 ,                                        &
     &          jintnfix   = 0 ,                                        &
     &          jintdnit   = 0 ,                                        &
     &          jndepfx    = 0 ,                                        &
     &          joalkfx    = 0 ,                                        &
     &          jcarflx0100= 0 ,                                        &
     &          jcarflx0500= 0 ,                                        &
     &          jcarflx1000= 0 ,                                        &
     &          jcarflx2000= 0 ,                                        &
     &          jcarflx4000= 0 ,                                        &
     &          jcarflx_bot= 0 ,                                        &
     &          jbsiflx0100= 0 ,                                        &
     &          jbsiflx0500= 0 ,                                        &
     &          jbsiflx1000= 0 ,                                        &
     &          jbsiflx2000= 0 ,                                        &
     &          jbsiflx4000= 0 ,                                        &
     &          jbsiflx_bot= 0 ,                                        &
     &          jcalflx0100= 0 ,                                        &
     &          jcalflx0500= 0 ,                                        &
     &          jcalflx1000= 0 ,                                        &
     &          jcalflx2000= 0 ,                                        &
     &          jcalflx4000= 0 ,                                        &
     &          jcalflx_bot= 0

      integer, dimension(nbgcmax), save ::                              &
     &          jsediffic  = 0 ,                                        &
     &          jsediffal  = 0 ,                                        &
     &          jsediffph  = 0 ,                                        &
     &          jsediffox  = 0 ,                                        &
     &          jsediffn2  = 0 ,                                        &
     &          jsediffno3 = 0 ,                                        &
                jsediffsi  = 0

      integer, dimension(nbgcmax), save ::                              &
     &          jsrfnatdic = 0 ,                                        &
     &          jsrfnatalk = 0 ,                                        &
     &          jnatpco2   = 0 ,                                        &
     &          jnatco2fx  = 0 ,                                        &
     &          jsrfnatph  = 0

      integer, dimension(nbgcmax), save ::                              &
     &          jbromofx   = 0 ,                                        &
     &          jsrfbromo  = 0 ,                                        &
     &          jbromo_prod= 0 ,                                        &
     &          jbromo_uv  = 0

      integer, save :: i_atm_m2d  
      integer, dimension(nbgcmax), save ::                              &
     &          jatmco2  = 0 ,                                          &
     &          jatmo2   = 0 ,                                          &
     &          jatmn2   = 0 ,                                          &
     &          jatmc13  = 0 ,                                          &
     &          jatmc14  = 0 ,                                          &
     &          jatmbromo= 0  

      integer, save :: nbgcm2d 

      logical, save :: domassfluxes = .false.

!----------------------------------------------------------------
      integer, save :: i_bsc_m3d,ilvl_bsc_m3d 
      integer, dimension(nbgcmax), save ::                              &
     &          jdp        = 0 ,                                        &
     &          jphyto     = 0 ,                                        &
     &          jgrazer    = 0 ,                                        &
     &          jdoc       = 0 ,                                        &
     &          jphosy     = 0 ,                                        &
     &          jphosph    = 0 ,                                        &
     &          joxygen    = 0 ,                                        &
     &          jiron      = 0 ,                                        &
     &          jano3      = 0 ,                                        &
     &          jalkali    = 0 ,                                        &
     &          jsilica    = 0 ,                                        &
     &          jdic       = 0 ,                                        &
     &          jpoc       = 0 ,                                        &
     &          jcalc      = 0 ,                                        &
     &          jopal      = 0 ,                                        &
     &          jco3       = 0 ,                                        &
     &          jph        = 0 ,                                        &
     &          jomegaa    = 0 ,                                        &
     &          jomegac    = 0 ,                                        &
     &          jn2o       = 0 ,                                        &
     &          jprefo2    = 0 ,                                        &
     &          jo2sat     = 0 ,                                        &
     &          jprefpo4   = 0 ,                                        &
     &          jprefalk   = 0 ,                                        &
     &          jprefdic   = 0 ,                                        &
     &          jdicsat    = 0 ,                                        &
     &          jcfc11     = 0 ,                                        &
     &          jcfc12     = 0 ,                                        &
     &          jsf6       = 0 ,                                        &
     &          jlvlphyto  = 0 ,                                        &
     &          jlvlgrazer = 0 ,                                        &
     &          jlvldoc    = 0 ,                                        &
     &          jlvlphosy  = 0 ,                                        &
     &          jlvlphosph = 0 ,                                        &
     &          jlvloxygen = 0 ,                                        &
     &          jlvliron   = 0 ,                                        &
     &          jlvlano3   = 0 ,                                        &
     &          jlvlalkali = 0 ,                                        &
     &          jlvlsilica = 0 ,                                        &
     &          jlvldic    = 0 ,                                        &
     &          jlvlpoc    = 0 ,                                        &
     &          jlvlcalc   = 0 ,                                        &
     &          jlvlopal   = 0 ,                                        &
     &          jlvlco3    = 0 ,                                        &
     &          jlvlph     = 0 ,                                        &
     &          jlvlomegaa = 0 ,                                        &
     &          jlvlomegac = 0 ,                                        &
     &          jlvln2o    = 0 ,                                        &
     &          jlvlprefo2 = 0 ,                                        &
     &          jlvlo2sat  = 0 ,                                        &
     &          jlvlprefpo4= 0 ,                                        &
     &          jlvlprefalk= 0 ,                                        &
     &          jlvlprefdic= 0 ,                                        &
     &          jlvldicsat = 0 ,                                        &
     &          jlvlcfc11  = 0 ,                                        &
     &          jlvlcfc12  = 0 ,                                        &
     &          jlvlsf6    = 0
    
      integer, dimension(nbgcmax), save ::                              &
     &          jdic13     = 0 ,                                        &
     &          jdic14     = 0 ,                                        &
     &          jd13c      = 0 ,                                        &
     &          jd14c      = 0 ,                                        &
     &          jbigd14c   = 0 ,                                        &
     &          jpoc13     = 0 ,                                        &
     &          jdoc13     = 0 ,                                        &
     &          jcalc13    = 0 ,                                        &
     &          jphyto13   = 0 ,                                        &
     &          jgrazer13  = 0 ,                                        &
     &          jlvldic13  = 0 ,                                        &
     &          jlvldic14  = 0 ,                                        &
     &          jlvld13c   = 0 ,                                        &
     &          jlvld14c   = 0 ,                                        &
     &          jlvlbigd14c= 0 ,                                        &
     &          jlvlpoc13  = 0 ,                                        &
     &          jlvldoc13  = 0 ,                                        &
     &          jlvlcalc13 = 0 ,                                        &
     &          jlvlphyto13 = 0,                                        &
     &          jlvlgrazer13= 0                                                                                        

      integer, dimension(nbgcmax), save ::                              &
     &          jnos       = 0 ,                                        &
     &          jwphy      = 0 ,                                        &
     &          jwnos      = 0 ,                                        &
     &          jeps       = 0 ,                                        &
     &          jasize     = 0 ,                                        &
     &          jlvlnos    = 0 ,                                        &
     &          jlvlwphy   = 0 ,                                        &
     &          jlvlwnos   = 0 ,                                        &
     &          jlvleps    = 0 ,                                        &
     &          jlvlasize  = 0

      integer, dimension(nbgcmax), save ::                              &
     &          jnatco3       = 0 ,                                     &
     &          jnatalkali    = 0 ,                                     &
     &          jnatdic       = 0 ,                                     &
     &          jnatcalc      = 0 ,                                     &
     &          jnatph        = 0 ,                                     &
     &          jnatomegaa    = 0 ,                                     &
     &          jnatomegac    = 0 ,                                     &
     &          jlvlnatco3    = 0 ,                                     &
     &          jlvlnatalkali = 0 ,                                     &
     &          jlvlnatdic    = 0 ,                                     &
     &          jlvlnatcalc   = 0 ,                                     &
     &          jlvlnatph     = 0 ,                                     &
     &          jlvlnatomegaa = 0 ,                                     &
     &          jlvlnatomegac = 0 

      integer, dimension(nbgcmax), save ::                              &
     &          jbromo     = 0 ,                                        &
     &          jlvlbromo  = 0               

      integer, save :: nbgcm3d,nbgcm3dlvl 

!----------------------------------------------------------------
! sediment
      integer, save :: i_bsc_sed 
      integer, dimension(nbgcmax), save ::                              &
     &          jpowaic = 0 ,                                           &
     &          jpowaal = 0 ,                                           &
     &          jpowaph = 0 ,                                           &
     &          jpowaox = 0 ,                                           &
     &          jpown2  = 0 ,                                           &
     &          jpowno3 = 0 ,                                           &
     &          jpowasi = 0 ,                                           &
     &          jssso12 = 0 ,                                           &
     &          jssssil = 0 ,                                           &
     &          jsssc12 = 0 ,                                           &
     &          jssster = 0              


      integer, save :: nbgct_sed    

!----------------------------------------------------------------
!  burial
      integer, save :: i_bsc_bur 
      integer, dimension(nbgcmax), save ::                              &
     &          jburssso12 = 0 ,                                        &
     &          jbursssc12 = 0 ,                                        &
     &          jburssssil = 0 ,                                        &
     &          jburssster = 0


      integer, save :: nbgct_bur

!----------------------------------------------------------------

      real, dimension (:,:,:),   allocatable :: bgct2d
      real, dimension (:,:,:),   allocatable :: bgcm2d
      real, dimension (:,:,:,:), allocatable :: bgcm3d,bgcm3dlvl
      real, dimension (:,:,:,:), allocatable :: bgct_sed 
      real, dimension (:,:,:),   allocatable :: bgct_bur 
     
 
      CONTAINS


      SUBROUTINE ALLOC_MEM_BGCMEAN(kpie,kpje,kpke)

      use mo_control_bgc, only: io_stdo_bgc,bgc_namelist,get_bgc_namelist

      implicit none 
     
      integer, intent(in) :: kpie,kpje,kpke

      integer             :: m,n,errstat,iounit,checkdp(nbgcmax)

!     Read namelist for diagnostic output
      GLB_AVEPERIO=0
      if(.not. allocated(bgc_namelist)) call get_bgc_namelist
      OPEN (newunit=iounit, file=bgc_namelist,                          &
           status='old', action='read', recl=80)
      READ (iounit,nml=diabgc)
      CLOSE (iounit)

!     Determine number of output groups 
      nbgc=0 
      do n=1,nbgcmax
        if (GLB_AVEPERIO(n).NE.0) then 
          nbgc=nbgc+1
          nacc_bgc(n)=0 
        endif
      enddo

      do n=1,nbgc
        GLB_FILEFREQ(n)=max(GLB_AVEPERIO(n),GLB_FILEFREQ(n))
        if (GLB_AVEPERIO(n).LT.0) then
          diagfq_bgc(n)=-real(nstepinday)/GLB_AVEPERIO(n)
        else
          diagfq_bgc(n)=nstepinday*max(1,GLB_AVEPERIO(n))
        endif
        diagmon_bgc(n)=.false.
        diagann_bgc(n)=.false.
        if (GLB_AVEPERIO(n).EQ.30) then
          diagmon_bgc(n)=.true.
        elseIF (GLB_AVEPERIO(n).EQ.365) then
          diagann_bgc(n)=.true.
        endif
        if (GLB_FILEFREQ(n).LT.0) then
          filefq_bgc(n)=-real(nstepinday)/GLB_FILEFREQ(n)
        else
          filefq_bgc(n)=nstepinday*max(1,GLB_FILEFREQ(n))
        endif
        filemon_bgc(n)=.false.
        fileann_bgc(n)=.false.
        if (GLB_FILEFREQ(n).EQ.30) then
          filemon_bgc(n)=.true.
        elseIF (GLB_FILEFREQ(n).EQ.365) then
          fileann_bgc(n)=.true.
        endif
      enddo

!     Re-define index variables according to namelist 
      i_bsc_m2d=0 
      do n=1,nbgc  
        if (SRF_KWCO2(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jkwco2(n)=i_bsc_m2d*min(1,SRF_KWCO2(n))
        if (SRF_KWCO2KHM(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jkwco2khm(n)=i_bsc_m2d*min(1,SRF_KWCO2KHM(n))
        if (SRF_CO2KH(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jco2kh(n)=i_bsc_m2d*min(1,SRF_CO2KH(n))
        if (SRF_CO2KHM(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jco2khm(n)=i_bsc_m2d*min(1,SRF_CO2KHM(n))
        if (SRF_PCO2(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jpco2(n)=i_bsc_m2d*min(1,SRF_PCO2(n))
        if (SRF_PCO2M(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jpco2m(n)=i_bsc_m2d*min(1,SRF_PCO2M(n))
        if (SRF_DMSFLUX(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jdmsflux(n)=i_bsc_m2d*min(1,SRF_DMSFLUX(n))
        if (SRF_CO2FXD(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jco2fxd(n)=i_bsc_m2d*min(1,SRF_CO2FXD(n))
        if (SRF_CO2FXU(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jco2fxu(n)=i_bsc_m2d*min(1,SRF_CO2FXU(n))
        if (SRF_OXFLUX(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        joxflux(n)=i_bsc_m2d*min(1,SRF_OXFLUX(n))
        if (SRF_NIFLUX(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jniflux(n)=i_bsc_m2d*min(1,SRF_NIFLUX(n))
        if (SRF_DMS(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jdms(n)=i_bsc_m2d*min(1,SRF_DMS(n))
        if (SRF_DMSPROD(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jdmsprod(n)=i_bsc_m2d*min(1,SRF_DMSPROD(n))
        if (SRF_DMS_BAC(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jdms_bac(n)=i_bsc_m2d*min(1,SRF_DMS_BAC(n))
        if (SRF_DMS_UV(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jdms_uv(n)=i_bsc_m2d*min(1,SRF_DMS_UV(n))
        if (SRF_EXPORT(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jexport(n)=i_bsc_m2d*min(1,SRF_EXPORT(n))
        if (SRF_EXPOCA(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jexpoca(n)=i_bsc_m2d*min(1,SRF_EXPOCA(n))
        if (SRF_EXPOSI(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jexposi(n)=i_bsc_m2d*min(1,SRF_EXPOSI(n))
        if (SRF_N2OFX(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jn2ofx(n)=i_bsc_m2d*min(1,SRF_N2OFX(n))
        if (SRF_PHOSPH(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jsrfphosph(n)=i_bsc_m2d*min(1,SRF_PHOSPH(n))
        if (SRF_OXYGEN(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jsrfoxygen(n)=i_bsc_m2d*min(1,SRF_OXYGEN(n))
        if (SRF_IRON(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jsrfiron(n)=i_bsc_m2d*min(1,SRF_IRON(n))
        if (SRF_ANO3(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jsrfano3(n)=i_bsc_m2d*min(1,SRF_ANO3(n))
        if (SRF_ALKALI(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jsrfalkali(n)=i_bsc_m2d*min(1,SRF_ALKALI(n))
        if (SRF_SILICA(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jsrfsilica(n)=i_bsc_m2d*min(1,SRF_SILICA(n))
        if (SRF_DIC(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jsrfdic(n)=i_bsc_m2d*min(1,SRF_DIC(n))
        if (SRF_PHYTO(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jsrfphyto(n)=i_bsc_m2d*min(1,SRF_PHYTO(n))
        if (SRF_PH(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jsrfph(n)=i_bsc_m2d*min(1,SRF_PH(n))
        if (INT_PHOSY(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jintphosy(n)=i_bsc_m2d*min(1,INT_PHOSY(n))
        if (INT_NFIX(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jintnfix(n)=i_bsc_m2d*min(1,INT_NFIX(n))
        if (INT_DNIT(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jintdnit(n)=i_bsc_m2d*min(1,INT_DNIT(n))
        if (FLX_NDEP(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jndepfx(n)=i_bsc_m2d*min(1,FLX_NDEP(n))
        if (FLX_OALK(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        joalkfx(n)=i_bsc_m2d*min(1,FLX_OALK(n))
        if (FLX_CAR0100(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcarflx0100(n)=i_bsc_m2d*min(1,FLX_CAR0100(n))
        if (FLX_CAR0500(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcarflx0500(n)=i_bsc_m2d*min(1,FLX_CAR0500(n))
        if (FLX_CAR1000(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcarflx1000(n)=i_bsc_m2d*min(1,FLX_CAR1000(n))
        if (FLX_CAR2000(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcarflx2000(n)=i_bsc_m2d*min(1,FLX_CAR2000(n))
        if (FLX_CAR4000(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcarflx4000(n)=i_bsc_m2d*min(1,FLX_CAR4000(n))
        if (FLX_CAR_BOT(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcarflx_bot(n)=i_bsc_m2d*min(1,FLX_CAR_BOT(n))
        if (FLX_BSI0100(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jbsiflx0100(n)=i_bsc_m2d*min(1,FLX_BSI0100(n))
        if (FLX_BSI0500(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jbsiflx0500(n)=i_bsc_m2d*min(1,FLX_BSI0500(n))
        if (FLX_BSI1000(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jbsiflx1000(n)=i_bsc_m2d*min(1,FLX_BSI1000(n))
        if (FLX_BSI2000(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jbsiflx2000(n)=i_bsc_m2d*min(1,FLX_BSI2000(n))
        if (FLX_BSI4000(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jbsiflx4000(n)=i_bsc_m2d*min(1,FLX_BSI4000(n))
        if (FLX_BSI_BOT(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jbsiflx_bot(n)=i_bsc_m2d*min(1,FLX_BSI_BOT(n))
        if (FLX_CAL0100(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcalflx0100(n)=i_bsc_m2d*min(1,FLX_CAL0100(n))
        if (FLX_CAL0500(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcalflx0500(n)=i_bsc_m2d*min(1,FLX_CAL0500(n))
        if (FLX_CAL1000(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcalflx1000(n)=i_bsc_m2d*min(1,FLX_CAL1000(n))
        if (FLX_CAL2000(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcalflx2000(n)=i_bsc_m2d*min(1,FLX_CAL2000(n))
        if (FLX_CAL4000(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcalflx4000(n)=i_bsc_m2d*min(1,FLX_CAL4000(n))
        if (FLX_CAL_BOT(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcalflx_bot(n)=i_bsc_m2d*min(1,FLX_CAL_BOT(n))
        if (.not. use_sedbypass) then
           if (FLX_SEDIFFIC(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
           jsediffic(n)=i_bsc_m2d*min(1,FLX_SEDIFFIC(n))
           if (FLX_SEDIFFAL(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
           jsediffal(n)=i_bsc_m2d*min(1,FLX_SEDIFFAL(n))
           if (FLX_SEDIFFPH(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
           jsediffph(n)=i_bsc_m2d*min(1,FLX_SEDIFFph(n))
           if (FLX_SEDIFFOX(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
           jsediffox(n)=i_bsc_m2d*min(1,FLX_SEDIFFOX(n))
           if (FLX_SEDIFFN2(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
           jsediffn2(n)=i_bsc_m2d*min(1,FLX_SEDIFFN2(n))
           if (FLX_SEDIFFNO3(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
           jsediffno3(n)=i_bsc_m2d*min(1,FLX_SEDIFFNO3(n))
           if (FLX_SEDIFFSI(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
           jsediffsi(n)=i_bsc_m2d*min(1,FLX_SEDIFFSI(n))
        endif
        if (use_cisonew) then
           if (SRF_CO213FXD(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
           jco213fxd(n)=i_bsc_m2d*min(1,SRF_CO213FXD(n))
           if (SRF_CO213FXU(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
           jco213fxu(n)=i_bsc_m2d*min(1,SRF_CO213FXU(n))
           if (SRF_CO214FXD(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
           jco214fxd(n)=i_bsc_m2d*min(1,SRF_CO214FXD(n))
           if (SRF_CO214FXU(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
           jco214fxu(n)=i_bsc_m2d*min(1,SRF_CO214FXU(n))
        endif
        if (use_CFC) then
           if (SRF_CFC11(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
           jcfc11fx(n)=i_bsc_m2d*min(1,SRF_CFC11(n))
           if (SRF_CFC12(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
           jcfc12fx(n)=i_bsc_m2d*min(1,SRF_CFC12(n))
           if (SRF_SF6(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
           jsf6fx(n)=i_bsc_m2d*min(1,SRF_SF6(n))
        endif
        if (use_natDIC) then
           if (SRF_NATDIC(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
           jsrfnatdic(n)=i_bsc_m2d*min(1,SRF_NATDIC(n))
           if (SRF_NATALKALI(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
           jsrfnatalk(n)=i_bsc_m2d*min(1,SRF_NATALKALI(n))
           if (SRF_NATPCO2(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
           jnatpco2(n)=i_bsc_m2d*min(1,SRF_NATPCO2(n))
           if (SRF_NATCO2FX(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
           jnatco2fx(n)=i_bsc_m2d*min(1,SRF_NATCO2FX(n))
           if (SRF_NATPH(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
           jsrfnatph(n)=i_bsc_m2d*min(1,SRF_NATPH(n))
        endif
        if (use_BROMO ) then
           if (SRF_BROMO(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
           jsrfbromo(n)=i_bsc_m2d*min(1,SRF_BROMO(n))
           if (SRF_BROMOFX(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
           jbromofx(n)=i_bsc_m2d*min(1,SRF_BROMOFX(n))
           if (INT_BROMOPRO(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
           jbromo_prod(n)=i_bsc_m2d*min(1,INT_BROMOPRO(n))
           if (INT_BROMOUV(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
           jbromo_uv(n)=i_bsc_m2d*min(1,INT_BROMOUV(n))
        endif
      enddo 

      domassfluxes = any(                                    &
        jcarflx0100+jcarflx0500+jcarflx1000+                 &
        jcarflx2000+jcarflx4000+jcarflx_bot+                 &
        jbsiflx0100+jbsiflx0500+jbsiflx1000+                 &
        jbsiflx2000+jbsiflx4000+jbsiflx_bot+                 &
        jcalflx0100+jcalflx0500+jcalflx1000+                 &
        jcalflx2000+jcalflx4000+jcalflx_bot  > 0)

      i_atm_m2d=i_bsc_m2d
      do n=1,nbgc
        if (SRF_ATMCO2(n).GT.0) i_atm_m2d=i_atm_m2d+1
        jatmco2(n)=i_atm_m2d*min(1,SRF_ATMCO2(n))
        if (use_BOXATM) then
           if (SRF_ATMO2(n).GT.0) i_atm_m2d=i_atm_m2d+1
           jatmo2(n)=i_atm_m2d*min(1,SRF_ATMO2(n))
           if (SRF_ATMN2(n).GT.0) i_atm_m2d=i_atm_m2d+1
           jatmn2(n)=i_atm_m2d*min(1,SRF_ATMN2(n))
        endif
        if (use_cisonew) then
           if (SRF_ATMC13(n).GT.0) i_atm_m2d=i_atm_m2d+1
           jatmc13(n)=i_atm_m2d*min(1,SRF_ATMC13(n))
           if (SRF_ATMC14(n).GT.0) i_atm_m2d=i_atm_m2d+1
           jatmc14(n)=i_atm_m2d*min(1,SRF_ATMC14(n))
        endif
        if (use_BROMO ) then
           if (SRF_ATMBROMO(n).GT.0) i_atm_m2d=i_atm_m2d+1
           jatmbromo(n)=i_atm_m2d*min(1,SRF_ATMBROMO(n))
        endif
      enddo 
      i_atm_m2d=i_atm_m2d-i_bsc_m2d

      i_bsc_m3d=0 
      ilvl_bsc_m3d=0 
      do n=1,nbgc 
        checkdp(n)=0

        if (LYR_PHYTO(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jphyto(n)=i_bsc_m3d*min(1,LYR_PHYTO(n))
        if (LYR_GRAZER(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jgrazer(n)=i_bsc_m3d*min(1,LYR_GRAZER(n))
        if (LYR_DOC(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jdoc(n)=i_bsc_m3d*min(1,LYR_DOC(n))
        if (LYR_PHOSY(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jphosy(n)=i_bsc_m3d*min(1,LYR_PHOSY(n))
        if (LYR_PHOSPH(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jphosph(n)=i_bsc_m3d*min(1,LYR_PHOSPH(n))
        if (LYR_OXYGEN(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        joxygen(n)=i_bsc_m3d*min(1,LYR_OXYGEN(n))
        if (LYR_IRON(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jiron(n)=i_bsc_m3d*min(1,LYR_IRON(n))
        if (LYR_ANO3(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jano3(n)=i_bsc_m3d*min(1,LYR_ANO3(n))
        if (LYR_ALKALI(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jalkali(n)=i_bsc_m3d*min(1,LYR_ALKALI(n))
        if (LYR_SILICA(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jsilica(n)=i_bsc_m3d*min(1,LYR_SILICA(n))
        if (LYR_DIC(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jdic(n)=i_bsc_m3d*min(1,LYR_DIC(n))
        if (LYR_POC(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jpoc(n)=i_bsc_m3d*min(1,LYR_POC(n))
        if (LYR_CALC(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jcalc(n)=i_bsc_m3d*min(1,LYR_CALC(n))
        if (LYR_OPAL(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jopal(n)=i_bsc_m3d*min(1,LYR_OPAL(n))
        if (LYR_CO3(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jco3(n)=i_bsc_m3d*min(1,LYR_CO3(n))
        if (LYR_PH(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jph(n)=i_bsc_m3d*min(1,LYR_PH(n))
        if (LYR_OMEGAA(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jomegaa(n)=i_bsc_m3d*min(1,LYR_OMEGAA(n))
        if (LYR_OMEGAC(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jomegac(n)=i_bsc_m3d*min(1,LYR_OMEGAC(n))
        if (LYR_N2O(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jn2o(n)=i_bsc_m3d*min(1,LYR_N2O(n))
        if (LYR_PREFO2(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jprefo2(n)=i_bsc_m3d*min(1,LYR_PREFO2(n))
        if (LYR_O2SAT(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jo2sat(n)=i_bsc_m3d*min(1,LYR_O2SAT(n))
        if (LYR_PREFPO4(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jprefpo4(n)=i_bsc_m3d*min(1,LYR_PREFPO4(n))
        if (LYR_PREFALK(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jprefalk(n)=i_bsc_m3d*min(1,LYR_PREFALK(n))
        if (LYR_PREFDIC(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jprefdic(n)=i_bsc_m3d*min(1,LYR_PREFDIC(n))
        if (LYR_DICSAT(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jdicsat(n)=i_bsc_m3d*min(1,LYR_DICSAT(n))
        if (LYR_DP(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jdp(n)=i_bsc_m3d*min(1,LYR_DP(n))
        if (use_CFC) then
           if (LYR_CFC11(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jcfc11(n)=i_bsc_m3d*min(1,LYR_CFC11(n))
           if (LYR_CFC12(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jcfc12(n)=i_bsc_m3d*min(1,LYR_CFC12(n))
           if (LYR_SF6(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jsf6(n)=i_bsc_m3d*min(1,LYR_SF6(n))
        endif
        if (use_cisonew) then
           if (LYR_DIC13(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jdic13(n)=i_bsc_m3d*min(1,LYR_DIC13(n))
           if (LYR_DIC14(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jdic14(n)=i_bsc_m3d*min(1,LYR_DIC14(n))
           if (LYR_D13C(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jd13c(n)=i_bsc_m3d*min(1,LYR_D13C(n))
           if (LYR_D14C(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jd14c(n)=i_bsc_m3d*min(1,LYR_D14C(n))
           if (LYR_BIGD14C(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jbigd14c(n)=i_bsc_m3d*min(1,LYR_BIGD14C(n))
           if (LYR_POC13(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jpoc13(n)=i_bsc_m3d*min(1,LYR_POC13(n))
           if (LYR_DOC13(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jdoc13(n)=i_bsc_m3d*min(1,LYR_DOC13(n))
           if (LYR_CALC13(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jcalc13(n)=i_bsc_m3d*min(1,LYR_CALC13(n))
           if (LYR_PHYTO13(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jphyto13(n)=i_bsc_m3d*min(1,LYR_PHYTO13(n))
           if (LYR_GRAZER13(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jgrazer13(n)=i_bsc_m3d*min(1,LYR_GRAZER13(n))
        endif
        if (use_AGG) then
           if (LYR_NOS(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jnos(n)=i_bsc_m3d*min(1,LYR_NOS(n))
           if (LYR_WPHY(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jwphy(n)=i_bsc_m3d*min(1,LYR_WPHY(n))
           if (LYR_WNOS(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jwnos(n)=i_bsc_m3d*min(1,LYR_WNOS(n))
           if (LYR_EPS(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jeps(n)=i_bsc_m3d*min(1,LYR_EPS(n))
           if (LYR_ASIZE(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jasize(n)=i_bsc_m3d*min(1,LYR_ASIZE(n))
        endif
        if (use_natDIC) then
           if (LYR_NATCO3(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jnatco3(n)=i_bsc_m3d*min(1,LYR_NATCO3(n))
           if (LYR_NATALKALI(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jnatalkali(n)=i_bsc_m3d*min(1,LYR_NATALKALI(n))
           if (LYR_NATDIC(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jnatdic(n)=i_bsc_m3d*min(1,LYR_NATDIC(n))
           if (LYR_NATCALC(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jnatcalc(n)=i_bsc_m3d*min(1,LYR_NATCALC(n))
           if (LYR_NATPH(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jnatph(n)=i_bsc_m3d*min(1,LYR_NATPH(n))
           if (LYR_NATOMEGAA(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jnatomegaa(n)=i_bsc_m3d*min(1,LYR_NATOMEGAA(n))
           if (LYR_NATOMEGAC(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jnatomegac(n)=i_bsc_m3d*min(1,LYR_NATOMEGAC(n))
        endif
        if (use_BROMO) then
           if (LYR_BROMO(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
           jbromo(n)=i_bsc_m3d*min(1,LYR_BROMO(n))
        endif

        if (LVL_PHYTO(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlphyto(n)=ilvl_bsc_m3d*min(1,LVL_PHYTO(n))
        if (LVL_GRAZER(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlgrazer(n)=ilvl_bsc_m3d*min(1,LVL_GRAZER(n))
        if (LVL_DOC(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvldoc(n)=ilvl_bsc_m3d*min(1,LVL_DOC(n))
        if (LVL_PHOSY(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlphosy(n)=ilvl_bsc_m3d*min(1,LVL_PHOSY(n))
        if (LVL_PHOSPH(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlphosph(n)=ilvl_bsc_m3d*min(1,LVL_PHOSPH(n))
        if (LVL_OXYGEN(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvloxygen(n)=ilvl_bsc_m3d*min(1,LVL_OXYGEN(n))
        if (LVL_IRON(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvliron(n)=ilvl_bsc_m3d*min(1,LVL_IRON(n))
        if (LVL_ANO3(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlano3(n)=ilvl_bsc_m3d*min(1,LVL_ANO3(n))
        if (LVL_ALKALI(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlalkali(n)=ilvl_bsc_m3d*min(1,LVL_ALKALI(n))
        if (LVL_SILICA(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlsilica(n)=ilvl_bsc_m3d*min(1,LVL_SILICA(n))
        if (LVL_DIC(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvldic(n)=ilvl_bsc_m3d*min(1,LVL_DIC(n))
        if (LVL_POC(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlpoc(n)=ilvl_bsc_m3d*min(1,LVL_POC(n))
        if (LVL_CALC(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlcalc(n)=ilvl_bsc_m3d*min(1,LVL_CALC(n))
        if (LVL_OPAL(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlopal(n)=ilvl_bsc_m3d*min(1,LVL_OPAL(n))
        if (LVL_CO3(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlco3(n)=ilvl_bsc_m3d*min(1,LVL_CO3(n))
        if (LVL_PH(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlph(n)=ilvl_bsc_m3d*min(1,LVL_PH(n))
        if (LVL_OMEGAA(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlomegaa(n)=ilvl_bsc_m3d*min(1,LVL_OMEGAA(n))
        if (LVL_OMEGAC(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlomegac(n)=ilvl_bsc_m3d*min(1,LVL_OMEGAC(n))
        if (LVL_N2O(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvln2o(n)=ilvl_bsc_m3d*min(1,LVL_N2O(n))
        if (LVL_PREFO2(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlprefo2(n)=ilvl_bsc_m3d*min(1,LVL_PREFO2(n))
        if (LVL_O2SAT(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlo2sat(n)=ilvl_bsc_m3d*min(1,LVL_O2SAT(n))
        if (LVL_PREFPO4(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlprefpo4(n)=ilvl_bsc_m3d*min(1,LVL_PREFPO4(n))
        if (LVL_PREFALK(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlprefalk(n)=ilvl_bsc_m3d*min(1,LVL_PREFALK(n))
        if (LVL_PREFDIC(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlprefdic(n)=ilvl_bsc_m3d*min(1,LVL_PREFDIC(n))
        if (LVL_DICSAT(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvldicsat(n)=ilvl_bsc_m3d*min(1,LVL_DICSAT(n))
        if (use_CFC) then
           if (LVL_CFC11(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlcfc11(n)=ilvl_bsc_m3d*min(1,LVL_CFC11(n))
           if (LVL_CFC12(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlcfc12(n)=ilvl_bsc_m3d*min(1,LVL_CFC12(n))
           if (LVL_SF6(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlsf6(n)=ilvl_bsc_m3d*min(1,LVL_SF6(n))
        endif
        if (use_cisonew) then
           if (LVL_DIC13(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvldic13(n)=ilvl_bsc_m3d*min(1,LVL_DIC13(n))
           if (LVL_DIC14(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvldic14(n)=ilvl_bsc_m3d*min(1,LVL_DIC14(n))
           if (LVL_D13C(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvld13c(n)=ilvl_bsc_m3d*min(1,LVL_D13C(n))
           if (LVL_D14C(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvld14c(n)=ilvl_bsc_m3d*min(1,LVL_D14C(n))
           if (LVL_BIGD14C(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlbigd14c(n)=ilvl_bsc_m3d*min(1,LVL_BIGD14C(n))
           if (LVL_POC13(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlpoc13(n)=ilvl_bsc_m3d*min(1,LVL_POC13(n))
           if (LVL_DOC13(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvldoc13(n)=ilvl_bsc_m3d*min(1,LVL_DOC13(n))
           if (LVL_CALC13(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlcalc13(n)=ilvl_bsc_m3d*min(1,LVL_CALC13(n))
           if (LVL_PHYTO13(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlphyto13(n)=ilvl_bsc_m3d*min(1,LVL_PHYTO13(n))
           if (LVL_GRAZER13(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlgrazer13(n)=ilvl_bsc_m3d*min(1,LVL_GRAZER13(n))
        endif
        if (use_AGG) then
           if (LVL_NOS(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlnos(n)=ilvl_bsc_m3d*min(1,LVL_NOS(n))
           if (LVL_WPHY(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlwphy(n)=ilvl_bsc_m3d*min(1,LVL_WPHY(n))
           if (LVL_WNOS(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlwnos(n)=ilvl_bsc_m3d*min(1,LVL_WNOS(n))
           if (LVL_EPS(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvleps(n)=ilvl_bsc_m3d*min(1,LVL_EPS(n))
           if (LVL_ASIZE(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlasize(n)=ilvl_bsc_m3d*min(1,LVL_ASIZE(n))
        endif
        if (use_natDIC) then
           if (LVL_NATCO3(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlnatco3(n)=ilvl_bsc_m3d*min(1,LVL_NATCO3(n))
           if (LVL_NATALKALI(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlnatalkali(n)=ilvl_bsc_m3d*min(1,LVL_NATALKALI(n))
           if (LVL_NATDIC(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlnatdic(n)=ilvl_bsc_m3d*min(1,LVL_NATDIC(n))
           if (LVL_NATCALC(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlnatcalc(n)=ilvl_bsc_m3d*min(1,LVL_NATCALC(n))
           if (LVL_NATPH(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlnatph(n)=ilvl_bsc_m3d*min(1,LVL_NATPH(n))
           if (LVL_NATOMEGAA(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlnatomegaa(n)=ilvl_bsc_m3d*min(1,LVL_NATOMEGAA(n))
           if (LVL_NATOMEGAC(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlnatomegac(n)=ilvl_bsc_m3d*min(1,LVL_NATOMEGAC(n))
        endif
        if (use_BROMO) then
           if (LVL_BROMO(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
           jlvlbromo(n)=ilvl_bsc_m3d*min(1,LVL_BROMO(n))
        endif

        if (i_bsc_m3d.NE.0) checkdp(n)=1
      enddo 
      

!     add dp required 
      do n=1,nbgc
        if (checkdp(n).NE.0.AND.LYR_DP(n).EQ.0) then 
          i_bsc_m3d=i_bsc_m3d+1
          jdp(n)=i_bsc_m3d
        endif
      enddo 
  
      i_bsc_sed=0
      i_bsc_bur=0
      if (.not. use_sedbypass) then
         do n=1,nbgc
            if (SDM_POWAIC(n).GT.0) i_bsc_sed=i_bsc_sed+1
            jpowaic(n)=i_bsc_sed*min(1,SDM_POWAIC(n))
            if (SDM_POWAAL(n).GT.0) i_bsc_sed=i_bsc_sed+1
            jpowaal(n)=i_bsc_sed*min(1,SDM_POWAAL(n))
            if (SDM_POWAPH(n).GT.0) i_bsc_sed=i_bsc_sed+1
            jpowaph(n)=i_bsc_sed*min(1,SDM_POWAPH(n))
            if (SDM_POWAOX(n).GT.0) i_bsc_sed=i_bsc_sed+1
            jpowaox(n)=i_bsc_sed*min(1,SDM_POWAOX(n))
            if (SDM_POWN2(n) .GT.0) i_bsc_sed=i_bsc_sed+1
            jpown2(n) =i_bsc_sed*min(1,SDM_POWN2(n))
            if (SDM_POWNO3(n).GT.0) i_bsc_sed=i_bsc_sed+1
            jpowno3(n)=i_bsc_sed*min(1,SDM_POWNO3(n))
            if (SDM_POWASI(n).GT.0) i_bsc_sed=i_bsc_sed+1
            jpowasi(n)=i_bsc_sed*min(1,SDM_POWASI(n))
            if (SDM_SSSO12(n).GT.0) i_bsc_sed=i_bsc_sed+1
            jssso12(n)=i_bsc_sed*min(1,SDM_SSSO12(n))
            if (SDM_SSSSIL(n).GT.0) i_bsc_sed=i_bsc_sed+1
            jssssil(n)=i_bsc_sed*min(1,SDM_SSSSIL(n))
            if (SDM_SSSC12(n).GT.0) i_bsc_sed=i_bsc_sed+1
            jsssc12(n)=i_bsc_sed*min(1,SDM_SSSC12(n))
            if (SDM_SSSTER(n).GT.0) i_bsc_sed=i_bsc_sed+1
            jssster(n)=i_bsc_sed*min(1,SDM_SSSTER(n))
         enddo

         do n=1,nbgc
            if (BUR_SSSO12(n).GT.0) i_bsc_bur=i_bsc_bur+1
            jburssso12(n)=i_bsc_bur*min(1,BUR_SSSO12(n))
            if (BUR_SSSC12(n).GT.0) i_bsc_bur=i_bsc_bur+1
            jbursssc12(n)=i_bsc_bur*min(1,BUR_SSSC12(n))
            if (BUR_SSSSIL(n).GT.0) i_bsc_bur=i_bsc_bur+1
            jburssssil(n)=i_bsc_bur*min(1,BUR_SSSSIL(n))
            if (BUR_SSSTER(n).GT.0) i_bsc_bur=i_bsc_bur+1
            jburssster(n)=i_bsc_bur*min(1,BUR_SSSTER(n))
         enddo
      endif
         
      nbgcm2d    = i_bsc_m2d+i_atm_m2d
      nbgcm3d    = i_bsc_m3d
      nbgcm3dlvl = ilvl_bsc_m3d
      nbgct_sed  = i_bsc_sed
      nbgct_bur  = i_bsc_bur

!     Allocate buffers 

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)' '
        write(io_stdo_bgc,*)'***************************************************'
        write(io_stdo_bgc,*)'Memory allocation for averaging model output :'
        write(io_stdo_bgc,*)' '
      endif


      if (mnproc.EQ.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable bgct2d ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',nbgct2d
      endif

      ALLOCATE (bgct2d(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,nbgct2d),      &
     &  stat=errstat)
      if (errstat.NE.0) STOP 'not enough memory bgct2d'
      if (nbgct2d.NE.0) bgct2d=0.     
 
      if (mnproc.EQ.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable bgcm2d ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',nbgcm2d
      endif

      ALLOCATE (bgcm2d(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,nbgcm2d),      &
     &  stat=errstat)
      if (errstat.NE.0) STOP 'not enough memory bgcm2d'
      if (nbgcm2d.NE.0) bgcm2d=0.

      if (mnproc.EQ.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable bgcm3d ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',kpke
        write(io_stdo_bgc,*)'Forth dimension    : ',nbgcm3d
      endif

      ALLOCATE (bgcm3d(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,kpke,nbgcm3d), &
     &  stat=errstat)
      if (errstat.NE.0) STOP 'not enough memory bgcm3d'
      if (nbgcm3d.NE.0) bgcm3d=0.

      if (mnproc.EQ.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable bgcm3dlvl '
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',kpke
        write(io_stdo_bgc,*)'Forth dimension    : ',nbgcm3dlvl
      endif

      ALLOCATE (bgcm3dlvl(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,ddm,        &
     &  nbgcm3dlvl),stat=errstat)
      if (errstat.NE.0) STOP 'not enough memory bgcm3dlvl'
      if (nbgcm3dlvl.NE.0) bgcm3dlvl=0.

      if (.not. use_sedbypass) then
         if (mnproc.EQ.1) then
            write(io_stdo_bgc,*)'Memory allocation for variable bgctsed ...'
            write(io_stdo_bgc,*)'First dimension    : ',kpie
            write(io_stdo_bgc,*)'Second dimension   : ',kpje
            write(io_stdo_bgc,*)'Third dimension    : ',ks
            write(io_stdo_bgc,*)'Forth dimension    : ',nbgct_sed
         endif

         ALLOCATE (bgct_sed(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,ks,          &
              &  nbgct_sed),stat=errstat)
         if (errstat.NE.0) STOP 'not enough memory bgct_sed'
         if (nbgct_sed.NE.0) bgct_sed=0. 

         if (mnproc.EQ.1) then
            write(io_stdo_bgc,*)'Memory allocation for variable bgctbur ...'
            write(io_stdo_bgc,*)'First dimension    : ',kpie
            write(io_stdo_bgc,*)'Second dimension   : ',kpje
            write(io_stdo_bgc,*)'Third dimension    : ',nbgct_bur
         endif

         ALLOCATE (bgct_bur(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,             &
              &  nbgct_bur),stat=errstat)
         if (errstat.NE.0) STOP 'not enough memory bgct_sed'
         if (nbgct_bur.NE.0) bgct_bur=0. 
      endif

      END SUBROUTINE ALLOC_MEM_BGCMEAN
 


      SUBROUTINE inisrf(pos,inival)
!
! --- ------------------------------------------------------------------
! --- Description: initialise 2d diagnostic field
! ---   
! --- Arguments:
! ---   int  pos      (in)     : position in common buffer  
! ---   real inival   (in)     : value used for initalisation
! --- ------------------------------------------------------------------
!
      implicit none
! 
      integer :: pos
      real ::inival
! 
      integer :: i,j,l
!
! --- Check whether field should be initialised
      if (pos.EQ.0) return
!
!$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            bgcm2d(i,j,pos)=inival
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
!
      END SUBROUTINE inisrf



      SUBROUTINE inilyr(pos,inival)
!
! --- ------------------------------------------------------------------
! --- Description: initialise layer diagnostic field
! ---   
! --- Arguments:
! ---   int  pos      (in)     : position in common buffer  
! ---   real inival   (in)     : value used for initalisation
! --- ------------------------------------------------------------------
!
      implicit none
! 
      integer :: pos
      real ::inival
! 
      integer :: i,j,k,l
!
! --- Check whether field should be initialised
      if (pos.EQ.0) return
!
      do k=1,kdm
!$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              bgcm3d(i,j,k,pos)=inival
            enddo
          enddo
        enddo
!$OMP END PARALLEL DO
      enddo
!
      END SUBROUTINE inilyr



      SUBROUTINE inilvl(pos,inival)
!
! --- ------------------------------------------------------------------
! --- Description: initialise level diagnostic field
! ---   
! --- Arguments:
! ---   int  pos      (in)     : position in common buffer  
! ---   real inival   (in)     : value used for initalisation
! --- ------------------------------------------------------------------
!
      implicit none
! 
      integer :: pos
      real ::inival
! 
      integer :: i,j,k,l
!
! --- Check whether field should be initialised
      if (pos.EQ.0) return
!
      do k=1,ddm
!$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              bgcm3dlvl(i,j,k,pos)=inival
            enddo
          enddo
        enddo
!$OMP END PARALLEL DO
      enddo
!
      END SUBROUTINE inilvl



      SUBROUTINE inisdm(pos,inival)
!
! --- ------------------------------------------------------------------
! --- Description: initialise sediment diagnostic field
! ---   
! --- Arguments:
! ---   int  pos      (in)     : position in common buffer  
! ---   real inival   (in)     : value used for initalisation
! --- ------------------------------------------------------------------
!
      implicit none
! 
      integer :: pos
      real ::inival
! 
      integer :: i,j,k,l
!
! --- Check whether field should be initialised
      if (pos.EQ.0) return
!
      do k=1,ks
!$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              bgct_sed(i,j,k,pos)=inival
            enddo
          enddo
        enddo
!$OMP END PARALLEL DO
      enddo
!
      END SUBROUTINE inisdm



      SUBROUTINE inibur(pos,inival)
!
! --- ------------------------------------------------------------------
! --- Description: initialise sediment burial diagnostic field
! ---   
! --- Arguments:
! ---   int  pos      (in)     : position in common buffer  
! ---   real inival   (in)     : value used for initalisation
! --- ------------------------------------------------------------------
!
      implicit none
! 
      integer :: pos
      real ::inival
! 
      integer :: i,j,l
!
! --- Check whether field should be initialised
      if (pos.EQ.0) return
!
!$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            bgct_bur(i,j,pos)=inival
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
!
      END SUBROUTINE inibur



      SUBROUTINE accsrf(pos,fld,wghts,wghtsflg)
!
! --- ------------------------------------------------------------------
! --- Description: accumulate 2d fields 
! ---  
! --- Arguments: 
! ---   int  pos      (in)     : position in 2d buffer  
! ---   real fld      (in)     : input data used for accumulation
! ---   real wghts    (in)     : weights used for accumulation
! --- ------------------------------------------------------------------
!
      implicit none
! 
      integer :: pos(nbgcmax),wghtsflg
      real, dimension(idm,jdm) :: fld,wghts
! 
      integer :: i,j,l,o
!
! --- Check whether field should be accumulated
      do o=1,nbgc
        if (pos(o).EQ.0) cycle
!
          if (wghtsflg.eq.0) then 
!$OMP PARALLEL DO PRIVATE(l,i)
            do j=1,jj
              do l=1,isp(j)
                do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                  bgcm2d(i,j,pos(o))=bgcm2d(i,j,pos(o))+fld(i,j)
                enddo
              enddo
            enddo
!$OMP END PARALLEL DO
          else
!$OMP PARALLEL DO PRIVATE(l,i)
            do j=1,jj
              do l=1,isp(j)
                do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                  bgcm2d(i,j,pos(o))=bgcm2d(i,j,pos(o))+fld(i,j)*       &
     &              wghts(i,j)
                enddo
              enddo
            enddo
!$OMP END PARALLEL DO
          endif 
!
      enddo
!   
      END SUBROUTINE accsrf



      SUBROUTINE acclyr(pos,fld,wghts,wghtsflg)
!
! --- ------------------------------------------------------------------
! --- Description: accumulate layer fields 
! ---  
! --- Arguments: 
! ---   int  pos      (in)     : position in 3d layer buffer  
! ---   real fld      (in)     : input data used for accumulation
! ---   real wghts    (in)     : weights used for accumulation
! --- ------------------------------------------------------------------
!
      implicit none
! 
      integer :: pos(nbgcmax),wghtsflg
      real, dimension(idm,jdm,kdm) :: fld,wghts
! 
      integer :: i,j,k,l,o
!
! --- Check whether field should be accumulated
      do o=1,nbgc
        if (pos(o).EQ.0) cycle
!
          if (wghtsflg.eq.0) then 
            do k=1,kk
!$OMP PARALLEL DO PRIVATE(l,i)
              do j=1,jj
                do l=1,isp(j)
                  do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                    bgcm3d(i,j,k,pos(o))=bgcm3d(i,j,k,pos(o))+          &
     &                fld(i,j,k)
                  enddo
                enddo
              enddo
!$OMP END PARALLEL DO
            enddo
          else
            do k=1,kk
!$OMP PARALLEL DO PRIVATE(l,i)
              do j=1,jj
                do l=1,isp(j)
                  do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                    bgcm3d(i,j,k,pos(o))=bgcm3d(i,j,k,pos(o))+          &
     &                fld(i,j,k)*wghts(i,j,k)
                  enddo
                enddo
              enddo
!$OMP END PARALLEL DO
            enddo
          endif
!
      enddo
      END SUBROUTINE acclyr



      SUBROUTINE acclvl(pos,fld,k,ind1,ind2,wghts)
!
! --- ------------------------------------------------------------------
! --- Description: accumulate 3d level fields
! ---  
! --- Arguments: 
! ---   int  pos      (in)     : position in buffer  
! ---   real fld      (in)     : input data used for accumulation
! ---   int  k        (in)     : layer index of fld  
! ---   int  ind1     (in)     : index field for first accumulated level 
! ---   int  ind2     (in)     : index field for last accumulated level 
! ---   real wghts    (in)     : weights used for accumulation
! --- ------------------------------------------------------------------
!
      implicit none
! 
      integer :: pos(nbgcmax),k
      integer, dimension(idm,jdm) :: ind1,ind2
      real, dimension(idm,jdm,ddm) :: wghts
      real, dimension(idm,jdm,kdm) :: fld
! 
      integer :: d,i,j,l,o
!
! --- Check whether field should be accumulated
      do o=1,nbgc
        if (pos(o).EQ.0) cycle
!
!$OMP PARALLEL DO PRIVATE(l,i,d)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              do d=ind1(i,j),ind2(i,j)
                bgcm3dlvl(i,j,d,pos(o))=bgcm3dlvl(i,j,d,pos(o))+        &
     &            fld(i,j,k)*wghts(i,j,d)
              enddo
            enddo
          enddo
        enddo
!$OMP END PARALLEL DO
      enddo
!     
      END SUBROUTINE acclvl



      SUBROUTINE accsdm(pos,fld)
!
! --- ------------------------------------------------------------------
! --- Description: accumulate sediment fields 
! ---  
! --- Arguments: 
! ---   int  pos      (in)     : position in 3d layer buffer  
! ---   real fld      (in)     : input data used for accumulation
! --- ------------------------------------------------------------------
!
      implicit none
! 
      integer :: pos(nbgcmax)
      real, dimension(idm,jdm,ks) :: fld
! 
      integer :: i,j,k,l,o
!
! --- Check whether field should be accumulated
      do o=1,nbgc
        if (pos(o).EQ.0) cycle
!
        do k=1,ks
!$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                bgct_sed(i,j,k,pos(o))=bgct_sed(i,j,k,pos(o))+fld(i,j,k)
              enddo
            enddo
          enddo
!$OMP END PARALLEL DO
        enddo
      enddo
!   
      END SUBROUTINE accsdm



      SUBROUTINE accbur(pos,fld)
!
! --- ------------------------------------------------------------------
! --- Description: accumulate sediment burial fields 
! ---  
! --- Arguments: 
! ---   int  pos      (in)     : position in 3d layer buffer  
! ---   real fld      (in)     : input data used for accumulation
! --- ------------------------------------------------------------------
!
      implicit none
! 
      integer :: pos(nbgcmax)
      real, dimension(idm,jdm) :: fld
! 
      integer :: i,j,l,o
!
! --- Check whether field should be accumulated
      do o=1,nbgc
        if (pos(o).EQ.0) cycle
!
!$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              bgct_bur(i,j,pos(o))=bgct_bur(i,j,pos(o))+fld(i,j)
            enddo
          enddo
        enddo
!$OMP END PARALLEL DO
      enddo
!   
      END SUBROUTINE accbur



      SUBROUTINE finsrf(posacc,poswgt)
!
! --- ------------------------------------------------------------------
! --- Description: finalise accumulation of weighted 2d fields 
! ---   
! --- Arguments:
! ---   real posacc   (in)     : position of accumulated field in buffer
! ---   real poswgt   (in)     : position of accumulated weights 
! --- ------------------------------------------------------------------
!
      implicit none
! 
      integer :: posacc,poswgt
! 
      integer :: i,j,l
      real, parameter :: epsil=1e-11
!
! --- Check whether field should be initialised
      if (posacc.EQ.0) return
!
!$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            bgcm2d(i,j,posacc)=bgcm2d(i,j,posacc)/                      &
     &        max(epsil,bgcm2d(i,j,poswgt))
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
!     
      END SUBROUTINE finsrf



      SUBROUTINE finlyr(posacc,poswgt)
!
! --- ------------------------------------------------------------------
! --- Description: finalise accumulation of weighted 3d layer fields 
! ---   
! --- Arguments:
! ---   real posacc   (in)     : position of accumulated field in buffer
! ---   real poswgt)  (in)     : position of accumulated weights 
! --- ------------------------------------------------------------------
!
      implicit none
! 
      integer :: posacc,poswgt
! 
      integer :: i,j,k,l
      real, parameter :: epsil=1e-11
!
! --- Check whether field should be initialised
      if (posacc.EQ.0) return
!
      do k=1,kk
!$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (bgcm3d(i,j,k,poswgt).GT.epsil) then
                bgcm3d(i,j,k,posacc)=bgcm3d(i,j,k,posacc)/              &
     &            bgcm3d(i,j,k,poswgt) 
              else 
                bgcm3d(i,j,k,posacc)=nf90_fill_double
              endif
            enddo
          enddo
        enddo
!$OMP END PARALLEL DO
      enddo
!     
      END SUBROUTINE finlyr



      SUBROUTINE wrtsrf(pos,frmt,sfac,offs,cmpflg,vnm)
!
! --- ------------------------------------------------------------------
! --- Description: writes diagnostic 2d field to file  
! ---   
! --- Arguments:
! ---   int  pos      (in)     : variable position in common buffer
! ---   int  frmt     (in)     : format/precision of output 
! ---                            0=field is not written  
! ---                            2=field is written as int2 with scale 
! ---                              factor and offset 
! ---                            4=field is written as real4
! ---                            8=field is written as real8
! ---   real sfac     (in)     : user def.NE. scale factor to be applied   
! ---   real offs     (in)     : user def.NE. offset to be added 
! ---   int  cmpflg   (in)     : compression flag; only wet points are 
! ---                            written if flag is set to 1 
! ---   char vnm      (in)     : variable name used in nc-file 
! --- ------------------------------------------------------------------
!
      implicit none
! 
      real,            intent(in) :: sfac,offs
      integer,         intent(in) :: frmt,cmpflg,pos
      character(LEN=*),intent(in) :: vnm
!
      integer                     :: n
      character(LEN=100)          :: dims
!
! --- Check whether field should be written
      if (pos.EQ.0 .OR. frmt.EQ.0) return
!
! --- Create dimension string 
      if (cmpflg.EQ.1) then
        dims='pcomp time'
      else
        dims='x y time'
      endif
!
! --- Check output format
      if (frmt.EQ.2) then
        if (cmpflg.EQ.1) then
          call nccopa(vnm,dims,bgcm2d(1-nbdy,1-nbdy,pos),ip,sfac,       &
     &      offs)
        else
          call ncpack(vnm,dims,bgcm2d(1-nbdy,1-nbdy,pos),ip,2,          &
     &      sfac,offs)
        endif
      elseIF (frmt.EQ.4) then
        if (cmpflg.EQ.1) then
          call nccomp(vnm,dims,bgcm2d(1-nbdy,1-nbdy,pos),ip,sfac,       &
     &      offs,4)
        else
          call ncwrtr(vnm,dims,bgcm2d(1-nbdy,1-nbdy,pos),ip,1,          &
     &      sfac,offs,4)
        endif
      elseIF (frmt.EQ.8) then
        if (cmpflg.EQ.1) then
          call nccomp(vnm,dims,bgcm2d(1-nbdy,1-nbdy,pos),ip,sfac,       &
     &      offs,8)
        else
          call ncwrtr(vnm,dims,bgcm2d(1-nbdy,1-nbdy,pos),ip,1,          &
     &      sfac,offs,8)
        endif
      else
        STOP 'unknown output format '
      endif
!
      END SUBROUTINE wrtsrf



      SUBROUTINE wrtlyr(pos,frmt,sfac,offs,cmpflg,vnm)
!
! --- ------------------------------------------------------------------
! --- Description: writes diagnostic layer field to file  
! ---   
! --- Arguments:
! ---   int  pos      (in)     : variable position in common buffer
! ---   int  frmt     (in)     : format/precision of output 
! ---                            0=field is not written  
! ---                            2=field is written as int2 with scale 
! ---                              factor and offset 
! ---                            4=field is written as real4
! ---                            8=field is written as real8
! ---   real sfac     (in)     : user def.NE. scale factor to be applied   
! ---   real offs     (in)     : user def.NE. offset to be added 
! ---   int  cmpflg   (in)     : compression flag; only wet points are 
! ---                            written if flag is set to 1 
! ---   char vnm      (in)     : variable name used in nc-file 
! --- ------------------------------------------------------------------
!
      implicit none
! 
      real,            intent(in) :: sfac,offs
      integer,         intent(in) :: frmt,cmpflg,pos
      character(LEN=*),intent(in) :: vnm
!
      integer                     :: n
      character(LEN=100)          :: dims
!
! --- Check whether field should be written
      if (pos.EQ.0 .OR. frmt.EQ.0) return
!
! --- Create dimension string 
      if (cmpflg.EQ.1) then
        dims='pcomp sigma time'
      else
        dims='x y sigma time'
      endif
!
! --- Check output format
      if (frmt.EQ.2) then
        if (cmpflg.EQ.1) then
          call nccopa(vnm,dims,bgcm3d(1-nbdy,1-nbdy,1,pos),ip,sfac,     &
     &      offs)
        else
          call ncpack(vnm,dims,bgcm3d(1-nbdy,1-nbdy,1,pos),ip,2,        &
     &      sfac,offs)
        endif
      elseIF (frmt.EQ.4) then
        if (cmpflg.EQ.1) then
          call nccomp(vnm,dims,bgcm3d(1-nbdy,1-nbdy,1,pos),ip,sfac,     &
     &      offs,4)
        else
          call ncwrtr(vnm,dims,bgcm3d(1-nbdy,1-nbdy,1,pos),ip,2,        &
     &      sfac,offs,4)
        endif
      elseIF (frmt.EQ.8) then
        if (cmpflg.EQ.1) then
          call nccomp(vnm,dims,bgcm3d(1-nbdy,1-nbdy,1,pos),ip,sfac,     &
     &      offs,8)
        else
          call ncwrtr(vnm,dims,bgcm3d(1-nbdy,1-nbdy,1,pos),ip,2,        &
     &      sfac,offs,8)
        endif
      else
        STOP 'unknown output format '
      endif
!
      END SUBROUTINE wrtlyr



      SUBROUTINE wrtlvl(pos,frmt,sfac,offs,cmpflg,vnm)
!
! --- ------------------------------------------------------------------
! --- Description: writes diagnostic level field to file  
! ---   
! --- Arguments:
! ---   int  pos      (in)     : variable position in common buffer
! ---   int  frmt     (in)     : format/precision of output 
! ---                            0=field is not written  
! ---                            2=field is written as int2 with scale 
! ---                              factor and offset 
! ---                            4=field is written as real4
! ---                            8=field is written as real8
! ---   real sfac     (in)     : user def.NE. scale factor to be applied   
! ---   real offs     (in)     : user def.NE. offset to be added 
! ---   int  cmpflg   (in)     : compression flag; only wet points are 
! ---                            written if flag is set to 1 
! ---   char vnm      (in)     : variable name used in nc-file 
! --- ------------------------------------------------------------------
!
      implicit none
! 
      real,            intent(in) :: sfac,offs
      integer,         intent(in) :: frmt,cmpflg,pos
      character(LEN=*),intent(in) :: vnm
!
      integer                     :: n
      character(LEN=100)          :: dims
!
! --- Check whether field should be written
      if (pos.EQ.0 .OR. frmt.EQ.0) return
!
! --- Create dimension string 
      if (cmpflg.EQ.1) then
        dims='pcomp depth time'
      else
        dims='x y depth time'
      endif
!
! --- Check output format
      if (frmt.EQ.2) then
        if (cmpflg.EQ.1) then
          call nccopa(vnm,dims,bgcm3dlvl(1-nbdy,1-nbdy,1,pos),ip,sfac,  &
     &      offs)
        else
          call ncpack(vnm,dims,bgcm3dlvl(1-nbdy,1-nbdy,1,pos),ip,2,     &
     &      sfac,offs)
        endif
      elseIF (frmt.EQ.4) then
        if (cmpflg.EQ.1) then
          call nccomp(vnm,dims,bgcm3dlvl(1-nbdy,1-nbdy,1,pos),ip,sfac,  &
     &      offs,4)
        else
          call ncwrtr(vnm,dims,bgcm3dlvl(1-nbdy,1-nbdy,1,pos),ip,2,     &
     &      sfac,offs,4)
        endif
      elseIF (frmt.EQ.8) then
        if (cmpflg.EQ.1) then
          call nccomp(vnm,dims,bgcm3dlvl(1-nbdy,1-nbdy,1,pos),ip,sfac,  &
     &      offs,8)
        else
          call ncwrtr(vnm,dims,bgcm3dlvl(1-nbdy,1-nbdy,1,pos),ip,2,     &
     &      sfac,offs,8)
        endif
      else
        STOP 'unknown output format '
      endif
!
      END SUBROUTINE wrtlvl



      SUBROUTINE wrtsdm(pos,frmt,sfac,offs,cmpflg,vnm)
!
! --- ------------------------------------------------------------------
! --- Description: writes diagnostic sediment field to file  
! ---   
! --- Arguments:
! ---   int  pos      (in)     : variable position in common buffer
! ---   int  frmt     (in)     : format/precision of output 
! ---                            0=field is not written  
! ---                            2=field is written as int2 with scale 
! ---                              factor and offset 
! ---                            4=field is written as real4
! ---                            8=field is written as real8
! ---   real sfac     (in)     : user def.NE. scale factor to be applied   
! ---   real offs     (in)     : user def.NE. offset to be added 
! ---   int  cmpflg   (in)     : compression flag; only wet points are 
! ---                            written if flag is set to 1 
! ---   char vnm      (in)     : variable name used in nc-file 
! --- ------------------------------------------------------------------
!
      implicit none
! 
      real,            intent(in) :: sfac,offs
      integer,         intent(in) :: frmt,cmpflg,pos
      character(LEN=*),intent(in) :: vnm
!
      integer                     :: n
      character(LEN=100)          :: dims
!
! --- Check whether field should be written
      if (pos.EQ.0 .OR. frmt.EQ.0) return
!
! --- Create dimension string 
      if (cmpflg.EQ.1) then
        dims='pcomp ks time'
      else
        dims='x y ks time'
      endif
!
! --- Check output format
      if (frmt.EQ.2) then
        if (cmpflg.EQ.1) then
          call nccopa(vnm,dims,bgct_sed(1-nbdy,1-nbdy,1,pos),ip,sfac,   &
     &      offs)
        else
          call ncpack(vnm,dims,bgct_sed(1-nbdy,1-nbdy,1,pos),ip,1,      &
     &      sfac,offs)
        endif
      elseIF (frmt.EQ.4) then
        if (cmpflg.EQ.1) then
          call nccomp(vnm,dims,bgct_sed(1-nbdy,1-nbdy,1,pos),ip,sfac,   &
     &      offs,4)
        else
          call ncwrtr(vnm,dims,bgct_sed(1-nbdy,1-nbdy,1,pos),ip,1,      &
     &      sfac,offs,4)
        endif
      elseIF (frmt.EQ.8) then
        if (cmpflg.EQ.1) then
          call nccomp(vnm,dims,bgct_sed(1-nbdy,1-nbdy,1,pos),ip,sfac,   &
     &      offs,8)
        else
          call ncwrtr(vnm,dims,bgct_sed(1-nbdy,1-nbdy,1,pos),ip,1,      &
     &      sfac,offs,8)
        endif
      else
        STOP 'unknown output format '
      endif
!
      END SUBROUTINE wrtsdm



      SUBROUTINE wrtbur(pos,frmt,sfac,offs,cmpflg,vnm)
!
! --- ------------------------------------------------------------------
! --- Description: writes diagnostic sediment burial field to file  
! ---   
! --- Arguments:
! ---   int  pos      (in)     : variable position in common buffer
! ---   int  frmt     (in)     : format/precision of output 
! ---                            0=field is not written  
! ---                            2=field is written as int2 with scale 
! ---                              factor and offset 
! ---                            4=field is written as real4
! ---                            8=field is written as real8
! ---   real sfac     (in)     : user def.NE. scale factor to be applied   
! ---   real offs     (in)     : user def.NE. offset to be added 
! ---   int  cmpflg   (in)     : compression flag; only wet points are 
! ---                            written if flag is set to 1 
! ---   char vnm      (in)     : variable name used in nc-file 
! --- ------------------------------------------------------------------
!
      implicit none
! 
      real,            intent(in) :: sfac,offs
      integer,         intent(in) :: frmt,cmpflg,pos
      character(LEN=*),intent(in) :: vnm
!
      integer                     :: n
      character(LEN=100)          :: dims
!
! --- Check whether field should be written
      if (pos.EQ.0 .OR. frmt.EQ.0) return
!
! --- Create dimension string 
      if (cmpflg.EQ.1) then
        dims='pcomp time'
      else
        dims='x y time'
      endif
!
! --- Check output format
      if (frmt.EQ.2) then
        if (cmpflg.EQ.1) then
          call nccopa(vnm,dims,bgct_bur(1-nbdy,1-nbdy,pos),ip,sfac,     &
     &      offs)
        else
          call ncpack(vnm,dims,bgct_bur(1-nbdy,1-nbdy,pos),ip,1,        &
     &      sfac,offs)
        endif
      elseIF (frmt.EQ.4) then
        if (cmpflg.EQ.1) then
          call nccomp(vnm,dims,bgct_bur(1-nbdy,1-nbdy,pos),ip,sfac,     &
     &      offs,4)
        else
          call ncwrtr(vnm,dims,bgct_bur(1-nbdy,1-nbdy,pos),ip,1,        &
     &      sfac,offs,4)
        endif
      elseIF (frmt.EQ.8) then
        if (cmpflg.EQ.1) then
          call nccomp(vnm,dims,bgct_bur(1-nbdy,1-nbdy,pos),ip,sfac,     &
     &      offs,8)
        else
          call ncwrtr(vnm,dims,bgct_bur(1-nbdy,1-nbdy,pos),ip,1,        &
     &      sfac,offs,8)
        endif
      else
        STOP 'unknown output format '
      endif
!
      END SUBROUTINE wrtbur



      SUBROUTINE logsrf(pos,sfac,offs)
!
! --- ------------------------------------------------------------------
! --- Description: replace 2d field with log10(field) 
! ---   
! --- Arguments:
! ---   int  pos      (in)     : field position in layer buffer 
! ---   real sfac     (in)     : scale factor to be applied before log10   
! ---   real offs     (in)     : offset to be added before log10   
! --- ------------------------------------------------------------------
!
      implicit none
! 
      real ::sfac,offs
      integer :: pos
! 
      integer :: i,j,l
      real ::epsil=1e-11
!
! --- Check whether field should be processed
      if (pos.EQ.0) return
!
!$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if (bgcm2d(i,j,pos).LT.epsil) then
              bgcm2d(i,j,pos)=0.
            else
              bgcm2d(i,j,pos)=log10(bgcm2d(i,j,pos)*sfac+offs)
            endif
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
!
      END SUBROUTINE logsrf



      SUBROUTINE loglyr(pos,sfac,offs)
!
! --- ------------------------------------------------------------------
! --- Description: replace layer field with log10(field) 
! ---   
! --- Arguments:
! ---   int  pos      (in)     : field position in layer buffer 
! ---   real sfac     (in)     : scale factor to be applied before log10   
! ---   real offs     (in)     : offset to be added before log10   
! --- ------------------------------------------------------------------
!
      implicit none
! 
      real ::sfac,offs
      integer :: pos
! 
      integer :: i,j,k,l
      real ::epsil=1e-11
!
! --- Check whether field should be processed
      if (pos.EQ.0) return
!
      do k=1,kdm
!$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (bgcm3d(i,j,k,pos).LT.epsil) then
                bgcm3d(i,j,k,pos)=0.
              elseIF (bgcm3d(i,j,k,pos).NE.nf90_fill_double) then
                bgcm3d(i,j,k,pos)=log10(bgcm3d(i,j,k,pos)*sfac+offs)
              endif
            enddo
          enddo
        enddo
!$OMP END PARALLEL DO
      enddo 
!
      END SUBROUTINE loglyr



      SUBROUTINE loglvl(pos,sfac,offs)
!
! --- ------------------------------------------------------------------
! --- Description: replace level field with log10(field) 
! ---   
! --- Arguments:
! ---   int  pos      (in)     : field position in layer buffer 
! ---   real sfac     (in)     : scale factor to be applied before log10   
! ---   real offs     (in)     : offset to be added before log10   
! --- ------------------------------------------------------------------
!
      implicit none
! 
      real ::sfac,offs
      integer :: pos
! 
      integer :: i,j,k,l
      real ::epsil=1e-11
!
! --- Check whether field should be processed
      if (pos.EQ.0) return
!
      do k=1,ddm
!$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (bgcm3dlvl(i,j,k,pos).LT.epsil) then
                bgcm3dlvl(i,j,k,pos)=0.
              elseIF (bgcm3dlvl(i,j,k,pos).NE.nf90_fill_double) then
                bgcm3dlvl(i,j,k,pos)=log10(bgcm3dlvl(i,j,k,pos)*sfac+   &
     &            offs)
              endif
            enddo
          enddo
        enddo
!$OMP END PARALLEL DO
      enddo
!
      END SUBROUTINE loglvl



      SUBROUTINE logsdm(pos,sfac,offs)
!
! --- ------------------------------------------------------------------
! --- Description: replace sediment field with log10(field) 
! ---   
! --- Arguments:
! ---   int  pos      (in)     : field position in layer buffer 
! ---   real sfac     (in)     : scale factor to be applied before log10   
! ---   real offs     (in)     : offset to be added before log10   
! --- ------------------------------------------------------------------
!
      implicit none
! 
      real ::sfac,offs
      integer :: pos
! 
      integer :: i,j,k,l
      real ::epsil=1e-11
!
! --- Check whether field should be processed
      if (pos.EQ.0) return
!
      do k=1,ks
!$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (bgct_sed(i,j,k,pos).LT.epsil) then
                bgct_sed(i,j,k,pos)=0.
              else
                bgct_sed(i,j,k,pos)=log10(bgct_sed(i,j,k,pos)*sfac+offs)
              endif
            enddo
          enddo
        enddo
!$OMP END PARALLEL DO
      enddo
!
      END SUBROUTINE logsdm


      SUBROUTINE msksrf(pos,idepth)
!
! --- ------------------------------------------------------------------
! --- Description: set sea floor points to NaN in mass flux fields
! ---   
! --- Arguments:
! ---   int  pos      (in)     : field position in level buffer 
! ---   int  idepth   (in)     : k-index field used to define the
! ---                            depth surface
! --- ------------------------------------------------------------------
!
      implicit none
!
      integer :: pos
      integer, dimension(idm,jdm) :: idepth 
! 
      integer :: i,j,l
      real, parameter :: mskval=nf90_fill_double
!
! --- Check whether field should be initia
      if (pos.EQ.0) return
!
!$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if( idepth(i,j) <= 0 ) bgcm2d(i,j,pos)=mskval
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
!
      END SUBROUTINE msksrf


      SUBROUTINE msklvl(pos,depths)
!
! --- ------------------------------------------------------------------
! --- Description: set sea floor points to NaN in level fields 
! ---   
! --- Arguments:
! ---   int  pos      (in)     : field position in level buffer 
! ---   int  depths   (in)     : bathymetry field 
! --- ------------------------------------------------------------------
!
      implicit none
!
      integer :: pos
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: depths 
! 
      integer :: i,j,k,l
      logical :: iniflg=.true.
      integer, dimension(idm,jdm) :: kmax
      real, parameter :: mskval=nf90_fill_double
!
      save iniflg,kmax
!     
! --- Check whether field should be processed
      if (pos.EQ.0) return
!
! --- Prepare index fields for masking

      if (iniflg) then
!$OMP PARALLEL DO PRIVATE(i)
        do j=1,jj
          do i=1,ii
            kmax(i,j)=0
          enddo
        enddo
!$OMP END PARALLEL DO
        do k=1,ddm
!$OMP PARALLEL DO PRIVATE(i)
          do j=1,jj
            do i=1,ii
              if (depths(i,j).GT.depthslev_bnds(1,k)) kmax(i,j)=k
            enddo
          enddo
!$OMP END PARALLEL DO
        enddo
        iniflg=.false.
      endif
!
!$OMP PARALLEL DO PRIVATE(i,k)
      do j=1,jj
        do i=1,ii
          do k=kmax(i,j)+1,ddm
            bgcm3dlvl(i,j,k,pos)=mskval
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
!
      END SUBROUTINE msklvl



      SUBROUTINE bgczlv(pddpo,kin,ind1,ind2,weights)
!-----------------------------------------------------------------------
!
!
      implicit none
!
      integer :: d,i,j,k,l,kin
      integer, dimension(idm,jdm) :: ind1,ind2
! 
      real, parameter :: eps=1e-10
      real, dimension(idm,jdm,kdm) :: pddpo,ztop,zbot
      real, dimension(idm,jdm,ddm) :: weights,dlev
!
      logical :: iniflg=.true.
!
      save ztop,zbot,dlev,iniflg
!
! --- Adjust bounds of levitus levels according to model bathymetry
      if (iniflg) then
        do d=1,ddm
!$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                 dlev(i,j,d)=max(eps,min(pbath(i,j),                    &
     &             depthslev_bnds(2,d))-depthslev_bnds(1,d))
              enddo
            enddo
          enddo
!$OMP END PARALLEL DO
        enddo
        iniflg=.false.
      endif
!
! --- Compute top and bottom depths of density layers 
      if (kin.EQ.1) then       
!$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              zbot(i,j,1)=pddpo(i,j,1)
            enddo
          enddo
        enddo
!$OMP END PARALLEL DO
        do k=2,kk
!$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                zbot(i,j,k)=zbot(i,j,k-1)+pddpo(i,j,k)
              enddo
            enddo
          enddo
!$OMP END PARALLEL DO
        enddo
!$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              zbot(i,j,1)=zbot(i,j,1)*pbath(i,j)/zbot(i,j,kk)
              ztop(i,j,1)=0.
              ind1(i,j)=1
            enddo
          enddo
        enddo
!$OMP END PARALLEL DO
        do k=2,kk
!$OMP PARALLEL DO PRIVATE(l,i)
          do j=1,jj
            do l=1,isp(j)
              do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                zbot(i,j,k)=zbot(i,j,k)*pbath(i,j)/zbot(i,j,kk)
                ztop(i,j,k)=zbot(i,j,k-1)
              enddo
            enddo
          enddo
!$OMP END PARALLEL DO
        enddo
      endif
!
! --- Compute interpolation weights 
!$OMP PARALLEL DO PRIVATE(l,i,d)
      do j=1,jj
        do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            ind2(i,j)=0
            if (pddpo(i,j,kin).GT.eps) then
              do d=ind1(i,j),ddm
                if (depthslev_bnds(2,d).LE.ztop(i,j,kin)) then
                  ind1(i,j)=d+1
                  CYCLE
                elseIF (depthslev_bnds(1,d).GE.zbot(i,j,kin)) then
                  EXIT
                endif
                ind2(i,j)=d
                weights(i,j,d)=(min(zbot(i,j,kin),                      &
     &            depthslev_bnds(2,d))-max(ztop(i,j,kin),               &
     &            depthslev_bnds(1,d)))/dlev(i,j,d)
              enddo
            endif
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
!
      END SUBROUTINE bgczlv

      END MODULE mo_bgcmean
