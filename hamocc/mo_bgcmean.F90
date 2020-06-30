! Copyright (C) 2002  P. Wetzel
! Copyright (C) 2020  I. Bethke, J. Tjiputra, J. Schwinger, A. Moree
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
      USE mod_xc, only: ii,jj,kk,idm,jdm,kdm,nbdy,ifp,isp,ilp
      USE mod_dia, only: ddm,depthslev,depthslev_bnds,nstepinday,pbath
      USE mod_nctools
      USE mo_param1_bgc, only: ks 

      IMPLICIT NONE

      PRIVATE :: ii,jj,kk,idm,jdm,kdm,nbdy,ifp,isp,ilp                
      PUBLIC :: ks,ddm,depthslev,depthslev_bnds

! --- Averaging and writing frequencies for diagnostic output     
      INTEGER, SAVE :: nbgc
      INTEGER, PARAMETER :: nbgcmax=10
      REAL, DIMENSION(nbgcmax), SAVE ::  diagfq_bgc,filefq_bgc,bgcwrt
      INTEGER, DIMENSION(nbgcmax), SAVE :: nacc_bgc
      LOGICAL, DIMENSION(nbgcmax), SAVE :: diagmon_bgc,diagann_bgc,     &
     &  filemon_bgc,fileann_bgc
 
! --- Namelist for diagnostic output 
      INTEGER, DIMENSION(nbgcmax), SAVE ::                              &
     & SRF_KWCO2     =0    ,SRF_PCO2      =0    ,SRF_DMSFLUX   =0    ,  &
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
     & SRF_NATDIC    =0    ,SRF_NATALKALI =0    ,SRF_NATPCO2   =0    ,  &
     & SRF_NATCO2FX  =0    ,                                            &
     & INT_PHOSY     =0    ,INT_NFIX      =0    ,INT_DNIT      =0    ,  &
     & FLX_CAR0100   =0    ,FLX_CAR0500   =0    ,FLX_CAR1000   =0    ,  &
     & FLX_CAR2000   =0    ,FLX_CAR4000   =0    ,FLX_CAR_BOT   =0    ,  &
     & FLX_BSI0100   =0    ,FLX_BSI0500   =0    ,FLX_BSI1000   =0    ,  &
     & FLX_BSI2000   =0    ,FLX_BSI4000   =0    ,FLX_BSI_BOT   =0    ,  &
     & FLX_CAL0100   =0    ,FLX_CAL0500   =0    ,FLX_CAL1000   =0    ,  &
     & FLX_CAL2000   =0    ,FLX_CAL4000   =0    ,FLX_CAL_BOT   =0    ,  &
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
      CHARACTER(LEN=10), DIMENSION(nbgcmax), SAVE :: GLB_FNAMETAG
      namelist /DIABGC/                                                 &
     & SRF_KWCO2         ,SRF_PCO2          ,SRF_DMSFLUX       ,        &
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
     & SRF_NATDIC        ,SRF_NATALKALI     ,SRF_NATPCO2       ,        &
     & SRF_NATCO2FX      ,                                              &
     & INT_PHOSY         ,INT_NFIX          ,INT_DNIT          ,        &
     & FLX_CAR0100       ,FLX_CAR0500       ,FLX_CAR1000       ,        &
     & FLX_CAR2000       ,FLX_CAR4000       ,FLX_CAR_BOT       ,        &
     & FLX_BSI0100       ,FLX_BSI0500       ,FLX_BSI1000       ,        &
     & FLX_BSI2000       ,FLX_BSI4000       ,FLX_BSI_BOT       ,        &
     & FLX_CAL0100       ,FLX_CAL0500       ,FLX_CAL1000       ,        &
     & FLX_CAL2000       ,FLX_CAL4000       ,FLX_CAL_BOT       ,        &
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
      INTEGER, parameter ::                                             &
     &          jco2flux  =1,                                           &
     &          jco214f   =2,                                           &
     &          jo2flux   =3,                                           &
     &          jn2flux   =4,                                           &
     &          jn2oflux  =5,                                           &
     &          jprorca   =6,                                           &
     &          jprcaca   =7,                                           &
     &          jsilpro   =8,                                           &
     &          jprodus   =9,                                           &
     &          nbgct2d   =9
      
!----------------------------------------------------------------      
      INTEGER, SAVE :: i_bsc_m2d 
      INTEGER, DIMENSION(nbgcmax), SAVE ::                              &
     &          jkwco2     = 0 ,                                        &
     &          jpco2      = 0 ,                                        &
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
     &          jintphosy  = 0 ,                                        &
     &          jintnfix   = 0 ,                                        &
     &          jintdnit   = 0 ,                                        &
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

      INTEGER, DIMENSION(nbgcmax), SAVE ::                              &
     &          jsrfnatdic = 0 ,                                        &
     &          jsrfnatalk = 0 ,                                        &
     &          jnatpco2   = 0 ,                                        &
     &          jnatco2fx  = 0


      INTEGER, SAVE :: i_atm_m2d  
      INTEGER, DIMENSION(nbgcmax), SAVE ::                              &
     &          jatmco2  = 0 ,                                          &
     &          jatmo2   = 0 ,                                          &
     &          jatmn2   = 0 ,                                          &
     &          jatmc13  = 0 ,                                          &
     &          jatmc14  = 0  

      INTEGER, SAVE :: nbgcm2d 

      LOGICAL, SAVE :: domassfluxes = .false.

!----------------------------------------------------------------
      INTEGER, SAVE :: i_bsc_m3d,ilvl_bsc_m3d 
      INTEGER, DIMENSION(nbgcmax), SAVE ::                              &
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
    
      INTEGER, DIMENSION(nbgcmax), SAVE ::                              &
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

      INTEGER, DIMENSION(nbgcmax), SAVE ::                              &
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

      INTEGER, DIMENSION(nbgcmax), SAVE ::                              &
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

      INTEGER, SAVE :: nbgcm3d,nbgcm3dlvl 

!----------------------------------------------------------------
! sediment
      INTEGER, SAVE :: i_bsc_sed 
      INTEGER, DIMENSION(nbgcmax), SAVE ::                              &
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


      INTEGER, SAVE :: nbgct_sed    

!----------------------------------------------------------------
!  burial
      INTEGER, SAVE :: i_bsc_bur 
      INTEGER, DIMENSION(nbgcmax), SAVE ::                              &
     &          jburssso12 = 0 ,                                        &
     &          jbursssc12 = 0 ,                                        &
     &          jburssssil = 0 ,                                        &
     &          jburssster = 0


      INTEGER, SAVE :: nbgct_bur

!----------------------------------------------------------------

      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: bgct2d
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: bgcm2d
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: bgcm3d,bgcm3dlvl
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: bgct_sed 
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: bgct_bur 
     
 
      CONTAINS


      SUBROUTINE ALLOC_MEM_BGCMEAN(kpie,kpje,kpke)

      USE mod_xc
      USE mo_control_bgc
      USE mo_param1_bgc 

      IMPLICIT NONE 
     
      INTEGER, intent(in) :: kpie,kpje,kpke

      INTEGER             :: m,n,errstat,iounit,checkdp(nbgcmax)
      LOGICAL             :: isopen,exists      

!     Read namelist for diagnostic output 
      GLB_AVEPERIO=0 
      DO iounit=10,99
        INQUIRE(unit=iounit,opened=isopen)
        IF (.NOT.isopen) EXIT
      ENDDO
      INQUIRE(file='ocn_in',exist=exists)
      IF (exists) THEN
        OPEN (iounit,file='ocn_in',status='old',action='read',recl=80)
      ELSE   
        INQUIRE(file='limits',exist=exists)
        IF (exists) THEN
          OPEN (iounit,file='limits',status='old',action='read',      &
     &      recl=80)
        ELSE
          STOP 'cannot find limits file' 
        ENDIF 
      ENDIF  
      READ (iounit,nml=diabgc)
      CLOSE (iounit)

!     Determine number of output groups 
      nbgc=0 
      DO n=1,nbgcmax
        IF (GLB_AVEPERIO(n).NE.0) THEN 
          nbgc=nbgc+1
          nacc_bgc(n)=0 
        ENDIF
      ENDDO

      DO n=1,nbgc
        GLB_FILEFREQ(n)=max(GLB_AVEPERIO(n),GLB_FILEFREQ(n))
        IF (GLB_AVEPERIO(n).LT.0) THEN
          diagfq_bgc(n)=-real(nstepinday)/GLB_AVEPERIO(n)
        ELSE
          diagfq_bgc(n)=nstepinday*max(1,GLB_AVEPERIO(n))
        ENDIF
        diagmon_bgc(n)=.false.
        diagann_bgc(n)=.false.
        IF (GLB_AVEPERIO(n).EQ.30) THEN
          diagmon_bgc(n)=.true.
        ELSEIF (GLB_AVEPERIO(n).EQ.365) THEN
          diagann_bgc(n)=.true.
        ENDIF
        IF (GLB_FILEFREQ(n).LT.0) THEN
          filefq_bgc(n)=-real(nstepinday)/GLB_FILEFREQ(n)
        ELSE
          filefq_bgc(n)=nstepinday*max(1,GLB_FILEFREQ(n))
        ENDIF
        filemon_bgc(n)=.false.
        fileann_bgc(n)=.false.
        IF (GLB_FILEFREQ(n).EQ.30) THEN
          filemon_bgc(n)=.true.
        ELSEIF (GLB_FILEFREQ(n).EQ.365) THEN
          fileann_bgc(n)=.true.
        ENDIF
      ENDDO

!     Re-define index variables according to namelist 
      i_bsc_m2d=0 
      DO n=1,nbgc  
        IF (SRF_KWCO2(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jkwco2(n)=i_bsc_m2d*min(1,SRF_KWCO2(n))
        IF (SRF_PCO2(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jpco2(n)=i_bsc_m2d*min(1,SRF_PCO2(n))
        IF (SRF_DMSFLUX(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jdmsflux(n)=i_bsc_m2d*min(1,SRF_DMSFLUX(n))
        IF (SRF_CO2FXD(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jco2fxd(n)=i_bsc_m2d*min(1,SRF_CO2FXD(n))
        IF (SRF_CO2FXU(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jco2fxu(n)=i_bsc_m2d*min(1,SRF_CO2FXU(n))
        IF (SRF_OXFLUX(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        joxflux(n)=i_bsc_m2d*min(1,SRF_OXFLUX(n))
        IF (SRF_NIFLUX(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jniflux(n)=i_bsc_m2d*min(1,SRF_NIFLUX(n))
        IF (SRF_DMS(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jdms(n)=i_bsc_m2d*min(1,SRF_DMS(n))
        IF (SRF_DMSPROD(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jdmsprod(n)=i_bsc_m2d*min(1,SRF_DMSPROD(n))
        IF (SRF_DMS_BAC(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jdms_bac(n)=i_bsc_m2d*min(1,SRF_DMS_BAC(n))
        IF (SRF_DMS_UV(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jdms_uv(n)=i_bsc_m2d*min(1,SRF_DMS_UV(n))
        IF (SRF_EXPORT(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jexport(n)=i_bsc_m2d*min(1,SRF_EXPORT(n))
        IF (SRF_EXPOCA(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jexpoca(n)=i_bsc_m2d*min(1,SRF_EXPOCA(n))
        IF (SRF_EXPOSI(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jexposi(n)=i_bsc_m2d*min(1,SRF_EXPOSI(n))
        IF (SRF_N2OFX(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jn2ofx(n)=i_bsc_m2d*min(1,SRF_N2OFX(n))
        IF (SRF_PHOSPH(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jsrfphosph(n)=i_bsc_m2d*min(1,SRF_PHOSPH(n))
        IF (SRF_OXYGEN(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jsrfoxygen(n)=i_bsc_m2d*min(1,SRF_OXYGEN(n))
        IF (SRF_IRON(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jsrfiron(n)=i_bsc_m2d*min(1,SRF_IRON(n))
        IF (SRF_ANO3(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jsrfano3(n)=i_bsc_m2d*min(1,SRF_ANO3(n))
        IF (SRF_ALKALI(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jsrfalkali(n)=i_bsc_m2d*min(1,SRF_ALKALI(n))
        IF (SRF_SILICA(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jsrfsilica(n)=i_bsc_m2d*min(1,SRF_SILICA(n))
        IF (SRF_DIC(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jsrfdic(n)=i_bsc_m2d*min(1,SRF_DIC(n))
        IF (SRF_PHYTO(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jsrfphyto(n)=i_bsc_m2d*min(1,SRF_PHYTO(n))
        IF (INT_PHOSY(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jintphosy(n)=i_bsc_m2d*min(1,INT_PHOSY(n))
        IF (INT_NFIX(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jintnfix(n)=i_bsc_m2d*min(1,INT_NFIX(n))
        IF (INT_DNIT(n).GT.0) i_bsc_m2d=i_bsc_m2d+1
        jintdnit(n)=i_bsc_m2d*min(1,INT_DNIT(n))
        IF (FLX_CAR0100(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcarflx0100(n)=i_bsc_m2d*min(1,FLX_CAR0100(n))
        IF (FLX_CAR0500(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcarflx0500(n)=i_bsc_m2d*min(1,FLX_CAR0500(n))
        IF (FLX_CAR1000(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcarflx1000(n)=i_bsc_m2d*min(1,FLX_CAR1000(n))
        IF (FLX_CAR2000(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcarflx2000(n)=i_bsc_m2d*min(1,FLX_CAR2000(n))
        IF (FLX_CAR4000(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcarflx4000(n)=i_bsc_m2d*min(1,FLX_CAR4000(n))
        IF (FLX_CAR_BOT(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcarflx_bot(n)=i_bsc_m2d*min(1,FLX_CAR_BOT(n))
        IF (FLX_BSI0100(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jbsiflx0100(n)=i_bsc_m2d*min(1,FLX_BSI0100(n))
        IF (FLX_BSI0500(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jbsiflx0500(n)=i_bsc_m2d*min(1,FLX_BSI0500(n))
        IF (FLX_BSI1000(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jbsiflx1000(n)=i_bsc_m2d*min(1,FLX_BSI1000(n))
        IF (FLX_BSI2000(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jbsiflx2000(n)=i_bsc_m2d*min(1,FLX_BSI2000(n))
        IF (FLX_BSI4000(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jbsiflx4000(n)=i_bsc_m2d*min(1,FLX_BSI4000(n))
        IF (FLX_BSI_BOT(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jbsiflx_bot(n)=i_bsc_m2d*min(1,FLX_BSI_BOT(n))
        IF (FLX_CAL0100(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcalflx0100(n)=i_bsc_m2d*min(1,FLX_CAL0100(n))
        IF (FLX_CAL0500(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcalflx0500(n)=i_bsc_m2d*min(1,FLX_CAL0500(n))
        IF (FLX_CAL1000(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcalflx1000(n)=i_bsc_m2d*min(1,FLX_CAL1000(n))
        IF (FLX_CAL2000(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcalflx2000(n)=i_bsc_m2d*min(1,FLX_CAL2000(n))
        IF (FLX_CAL4000(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcalflx4000(n)=i_bsc_m2d*min(1,FLX_CAL4000(n))
        IF (FLX_CAL_BOT(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcalflx_bot(n)=i_bsc_m2d*min(1,FLX_CAL_BOT(n))
#ifdef cisonew
        IF (SRF_CO213FXD(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jco213fxd(n)=i_bsc_m2d*min(1,SRF_CO213FXD(n))
        IF (SRF_CO213FXU(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jco213fxu(n)=i_bsc_m2d*min(1,SRF_CO213FXU(n))
        IF (SRF_CO214FXD(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jco214fxd(n)=i_bsc_m2d*min(1,SRF_CO214FXD(n))
        IF (SRF_CO214FXU(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jco214fxu(n)=i_bsc_m2d*min(1,SRF_CO214FXU(n))
#endif
#ifdef CFC
        IF (SRF_CFC11(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcfc11fx(n)=i_bsc_m2d*min(1,SRF_CFC11(n))
        IF (SRF_CFC12(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jcfc12fx(n)=i_bsc_m2d*min(1,SRF_CFC12(n))
        IF (SRF_SF6(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jsf6fx(n)=i_bsc_m2d*min(1,SRF_SF6(n))
#endif
#ifdef natDIC
        IF (SRF_NATDIC(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jsrfnatdic(n)=i_bsc_m2d*min(1,SRF_NATDIC(n))
        IF (SRF_NATALKALI(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jsrfnatalk(n)=i_bsc_m2d*min(1,SRF_NATALKALI(n))
        IF (SRF_NATPCO2(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jnatpco2(n)=i_bsc_m2d*min(1,SRF_NATPCO2(n))
        IF (SRF_NATCO2FX(n).GT.0) i_bsc_m2d=i_bsc_m2d+1 
        jnatco2fx(n)=i_bsc_m2d*min(1,SRF_NATCO2FX(n))
#endif
      ENDDO 

      domassfluxes = any(                                    &
        jcarflx0100+jcarflx0500+jcarflx1000+                 &
        jcarflx2000+jcarflx4000+jcarflx_bot+                 &
        jbsiflx0100+jbsiflx0500+jbsiflx1000+                 &
        jbsiflx2000+jbsiflx4000+jbsiflx_bot+                 &
        jcalflx0100+jcalflx0500+jcalflx1000+                 &
        jcalflx2000+jcalflx4000+jcalflx_bot  > 0)

      i_atm_m2d=i_bsc_m2d
      DO n=1,nbgc
        IF (SRF_ATMCO2(n).GT.0) i_atm_m2d=i_atm_m2d+1
        jatmco2(n)=i_atm_m2d*min(1,SRF_ATMCO2(n))
#if defined(BOXATM)
        IF (SRF_ATMO2(n).GT.0) i_atm_m2d=i_atm_m2d+1
        jatmo2(n)=i_atm_m2d*min(1,SRF_ATMO2(n))
        IF (SRF_ATMN2(n).GT.0) i_atm_m2d=i_atm_m2d+1
        jatmn2(n)=i_atm_m2d*min(1,SRF_ATMN2(n))
#endif 
#ifdef cisonew
        IF (SRF_ATMC13(n).GT.0) i_atm_m2d=i_atm_m2d+1
        jatmc13(n)=i_atm_m2d*min(1,SRF_ATMC13(n))
        IF (SRF_ATMC14(n).GT.0) i_atm_m2d=i_atm_m2d+1
        jatmc14(n)=i_atm_m2d*min(1,SRF_ATMC14(n))
#endif
      ENDDO 
      i_atm_m2d=i_atm_m2d-i_bsc_m2d

      i_bsc_m3d=0 
      ilvl_bsc_m3d=0 
      DO n=1,nbgc 
        checkdp(n)=0

        IF (LYR_PHYTO(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jphyto(n)=i_bsc_m3d*min(1,LYR_PHYTO(n))
        IF (LYR_GRAZER(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jgrazer(n)=i_bsc_m3d*min(1,LYR_GRAZER(n))
        IF (LYR_DOC(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jdoc(n)=i_bsc_m3d*min(1,LYR_DOC(n))
        IF (LYR_PHOSY(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jphosy(n)=i_bsc_m3d*min(1,LYR_PHOSY(n))
        IF (LYR_PHOSPH(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jphosph(n)=i_bsc_m3d*min(1,LYR_PHOSPH(n))
        IF (LYR_OXYGEN(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        joxygen(n)=i_bsc_m3d*min(1,LYR_OXYGEN(n))
        IF (LYR_IRON(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jiron(n)=i_bsc_m3d*min(1,LYR_IRON(n))
        IF (LYR_ANO3(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jano3(n)=i_bsc_m3d*min(1,LYR_ANO3(n))
        IF (LYR_ALKALI(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jalkali(n)=i_bsc_m3d*min(1,LYR_ALKALI(n))
        IF (LYR_SILICA(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jsilica(n)=i_bsc_m3d*min(1,LYR_SILICA(n))
        IF (LYR_DIC(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jdic(n)=i_bsc_m3d*min(1,LYR_DIC(n))
        IF (LYR_POC(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jpoc(n)=i_bsc_m3d*min(1,LYR_POC(n))
        IF (LYR_CALC(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jcalc(n)=i_bsc_m3d*min(1,LYR_CALC(n))
        IF (LYR_OPAL(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jopal(n)=i_bsc_m3d*min(1,LYR_OPAL(n))
        IF (LYR_CO3(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jco3(n)=i_bsc_m3d*min(1,LYR_CO3(n))
        IF (LYR_PH(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jph(n)=i_bsc_m3d*min(1,LYR_PH(n))
        IF (LYR_OMEGAA(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jomegaa(n)=i_bsc_m3d*min(1,LYR_OMEGAA(n))
        IF (LYR_OMEGAC(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jomegac(n)=i_bsc_m3d*min(1,LYR_OMEGAC(n))
        IF (LYR_N2O(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jn2o(n)=i_bsc_m3d*min(1,LYR_N2O(n))
        IF (LYR_PREFO2(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jprefo2(n)=i_bsc_m3d*min(1,LYR_PREFO2(n))
        IF (LYR_O2SAT(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jo2sat(n)=i_bsc_m3d*min(1,LYR_O2SAT(n))
        IF (LYR_PREFPO4(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jprefpo4(n)=i_bsc_m3d*min(1,LYR_PREFPO4(n))
        IF (LYR_PREFALK(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jprefalk(n)=i_bsc_m3d*min(1,LYR_PREFALK(n))
        IF (LYR_PREFDIC(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jprefdic(n)=i_bsc_m3d*min(1,LYR_PREFDIC(n))
        IF (LYR_DICSAT(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jdicsat(n)=i_bsc_m3d*min(1,LYR_DICSAT(n))
        IF (LYR_DP(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jdp(n)=i_bsc_m3d*min(1,LYR_DP(n))
#ifdef CFC
        IF (LYR_CFC11(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jcfc11(n)=i_bsc_m3d*min(1,LYR_CFC11(n))
        IF (LYR_CFC12(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jcfc12(n)=i_bsc_m3d*min(1,LYR_CFC12(n))
        IF (LYR_SF6(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jsf6(n)=i_bsc_m3d*min(1,LYR_SF6(n))
#endif
#ifdef cisonew
        IF (LYR_DIC13(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jdic13(n)=i_bsc_m3d*min(1,LYR_DIC13(n))
        IF (LYR_DIC14(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jdic14(n)=i_bsc_m3d*min(1,LYR_DIC14(n))
        IF (LYR_D13C(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jd13c(n)=i_bsc_m3d*min(1,LYR_D13C(n))
        IF (LYR_D14C(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jd14c(n)=i_bsc_m3d*min(1,LYR_D14C(n))
        IF (LYR_BIGD14C(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jbigd14c(n)=i_bsc_m3d*min(1,LYR_BIGD14C(n))
        IF (LYR_POC13(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jpoc13(n)=i_bsc_m3d*min(1,LYR_POC13(n))
        IF (LYR_DOC13(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jdoc13(n)=i_bsc_m3d*min(1,LYR_DOC13(n))
        IF (LYR_CALC13(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jcalc13(n)=i_bsc_m3d*min(1,LYR_CALC13(n))
        IF (LYR_PHYTO13(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jphyto13(n)=i_bsc_m3d*min(1,LYR_PHYTO13(n))
        IF (LYR_GRAZER13(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jgrazer13(n)=i_bsc_m3d*min(1,LYR_GRAZER13(n))
#endif 
#ifdef AGG
        IF (LYR_NOS(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jnos(n)=i_bsc_m3d*min(1,LYR_NOS(n))
        IF (LYR_WPHY(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jwphy(n)=i_bsc_m3d*min(1,LYR_WPHY(n))
        IF (LYR_WNOS(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jwnos(n)=i_bsc_m3d*min(1,LYR_WNOS(n))
        IF (LYR_EPS(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jeps(n)=i_bsc_m3d*min(1,LYR_EPS(n))
        IF (LYR_ASIZE(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jasize(n)=i_bsc_m3d*min(1,LYR_ASIZE(n))
#endif
#ifdef natDIC
        IF (LYR_NATCO3(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jnatco3(n)=i_bsc_m3d*min(1,LYR_NATCO3(n))
        IF (LYR_NATALKALI(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jnatalkali(n)=i_bsc_m3d*min(1,LYR_NATALKALI(n))
        IF (LYR_NATDIC(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jnatdic(n)=i_bsc_m3d*min(1,LYR_NATDIC(n))
        IF (LYR_NATCALC(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jnatcalc(n)=i_bsc_m3d*min(1,LYR_NATCALC(n))
        IF (LYR_NATPH(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jnatph(n)=i_bsc_m3d*min(1,LYR_NATPH(n))
        IF (LYR_NATOMEGAA(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jnatomegaa(n)=i_bsc_m3d*min(1,LYR_NATOMEGAA(n))
        IF (LYR_NATOMEGAC(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jnatomegac(n)=i_bsc_m3d*min(1,LYR_NATOMEGAC(n))
#endif

        IF (LVL_PHYTO(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlphyto(n)=ilvl_bsc_m3d*min(1,LVL_PHYTO(n))
        IF (LVL_GRAZER(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlgrazer(n)=ilvl_bsc_m3d*min(1,LVL_GRAZER(n))
        IF (LVL_DOC(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvldoc(n)=ilvl_bsc_m3d*min(1,LVL_DOC(n))
        IF (LVL_PHOSY(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlphosy(n)=ilvl_bsc_m3d*min(1,LVL_PHOSY(n))
        IF (LVL_PHOSPH(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlphosph(n)=ilvl_bsc_m3d*min(1,LVL_PHOSPH(n))
        IF (LVL_OXYGEN(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvloxygen(n)=ilvl_bsc_m3d*min(1,LVL_OXYGEN(n))
        IF (LVL_IRON(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvliron(n)=ilvl_bsc_m3d*min(1,LVL_IRON(n))
        IF (LVL_ANO3(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlano3(n)=ilvl_bsc_m3d*min(1,LVL_ANO3(n))
        IF (LVL_ALKALI(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlalkali(n)=ilvl_bsc_m3d*min(1,LVL_ALKALI(n))
        IF (LVL_SILICA(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlsilica(n)=ilvl_bsc_m3d*min(1,LVL_SILICA(n))
        IF (LVL_DIC(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvldic(n)=ilvl_bsc_m3d*min(1,LVL_DIC(n))
        IF (LVL_POC(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlpoc(n)=ilvl_bsc_m3d*min(1,LVL_POC(n))
        IF (LVL_CALC(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlcalc(n)=ilvl_bsc_m3d*min(1,LVL_CALC(n))
        IF (LVL_OPAL(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlopal(n)=ilvl_bsc_m3d*min(1,LVL_OPAL(n))
        IF (LVL_CO3(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlco3(n)=ilvl_bsc_m3d*min(1,LVL_CO3(n))
        IF (LVL_PH(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlph(n)=ilvl_bsc_m3d*min(1,LVL_PH(n))
        IF (LVL_OMEGAA(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlomegaa(n)=ilvl_bsc_m3d*min(1,LVL_OMEGAA(n))
        IF (LVL_OMEGAC(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlomegac(n)=ilvl_bsc_m3d*min(1,LVL_OMEGAC(n))
        IF (LVL_N2O(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvln2o(n)=ilvl_bsc_m3d*min(1,LVL_N2O(n))
        IF (LVL_PREFO2(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlprefo2(n)=ilvl_bsc_m3d*min(1,LVL_PREFO2(n))
        IF (LVL_O2SAT(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlo2sat(n)=ilvl_bsc_m3d*min(1,LVL_O2SAT(n))
        IF (LVL_PREFPO4(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlprefpo4(n)=ilvl_bsc_m3d*min(1,LVL_PREFPO4(n))
        IF (LVL_PREFALK(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlprefalk(n)=ilvl_bsc_m3d*min(1,LVL_PREFALK(n))
        IF (LVL_PREFDIC(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlprefdic(n)=ilvl_bsc_m3d*min(1,LVL_PREFDIC(n))
        IF (LVL_DICSAT(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvldicsat(n)=ilvl_bsc_m3d*min(1,LVL_DICSAT(n))
#ifdef CFC
        IF (LVL_CFC11(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlcfc11(n)=ilvl_bsc_m3d*min(1,LVL_CFC11(n))
        IF (LVL_CFC12(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlcfc12(n)=ilvl_bsc_m3d*min(1,LVL_CFC12(n))
        IF (LVL_SF6(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlsf6(n)=ilvl_bsc_m3d*min(1,LVL_SF6(n))
#endif
#ifdef cisonew
        IF (LVL_DIC13(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvldic13(n)=ilvl_bsc_m3d*min(1,LVL_DIC13(n))
        IF (LVL_DIC14(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvldic14(n)=ilvl_bsc_m3d*min(1,LVL_DIC14(n))
        IF (LVL_D13C(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvld13c(n)=ilvl_bsc_m3d*min(1,LVL_D13C(n))
        IF (LVL_D14C(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvld14c(n)=ilvl_bsc_m3d*min(1,LVL_D14C(n))
        IF (LVL_BIGD14C(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlbigd14c(n)=ilvl_bsc_m3d*min(1,LVL_BIGD14C(n))
        IF (LVL_POC13(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlpoc13(n)=ilvl_bsc_m3d*min(1,LVL_POC13(n))
        IF (LVL_DOC13(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvldoc13(n)=ilvl_bsc_m3d*min(1,LVL_DOC13(n))
        IF (LVL_CALC13(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlcalc13(n)=ilvl_bsc_m3d*min(1,LVL_CALC13(n))
        IF (LVL_PHYTO13(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlphyto13(n)=ilvl_bsc_m3d*min(1,LVL_PHYTO13(n))
        IF (LVL_GRAZER13(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlgrazer13(n)=ilvl_bsc_m3d*min(1,LVL_GRAZER13(n))
#endif
#ifdef AGG
        IF (LVL_NOS(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlnos(n)=ilvl_bsc_m3d*min(1,LVL_NOS(n))
        IF (LVL_WPHY(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlwphy(n)=ilvl_bsc_m3d*min(1,LVL_WPHY(n))
        IF (LVL_WNOS(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlwnos(n)=ilvl_bsc_m3d*min(1,LVL_WNOS(n))
        IF (LVL_EPS(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvleps(n)=ilvl_bsc_m3d*min(1,LVL_EPS(n))
        IF (LVL_ASIZE(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlasize(n)=ilvl_bsc_m3d*min(1,LVL_ASIZE(n))
#endif
#ifdef natDIC
        IF (LVL_NATCO3(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlnatco3(n)=ilvl_bsc_m3d*min(1,LVL_NATCO3(n))
        IF (LVL_NATALKALI(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlnatalkali(n)=ilvl_bsc_m3d*min(1,LVL_NATALKALI(n))
        IF (LVL_NATDIC(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlnatdic(n)=ilvl_bsc_m3d*min(1,LVL_NATDIC(n))
        IF (LVL_NATCALC(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlnatcalc(n)=ilvl_bsc_m3d*min(1,LVL_NATCALC(n))
        IF (LVL_NATPH(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlnatph(n)=ilvl_bsc_m3d*min(1,LVL_NATPH(n))
        IF (LVL_NATOMEGAA(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlnatomegaa(n)=ilvl_bsc_m3d*min(1,LVL_NATOMEGAA(n))
        IF (LVL_NATOMEGAC(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlnatomegac(n)=ilvl_bsc_m3d*min(1,LVL_NATOMEGAC(n))
#endif

        IF (i_bsc_m3d.NE.0) checkdp(n)=1
      ENDDO 
      

!     add dp required 
      DO n=1,nbgc
        IF (checkdp(n).NE.0.AND.LYR_DP(n).EQ.0) THEN 
          i_bsc_m3d=i_bsc_m3d+1
          jdp(n)=i_bsc_m3d
        ENDIF
      ENDDO 
  
      i_bsc_sed=0
      i_bsc_bur=0
#ifndef sedbypass
      DO n=1,nbgc
        IF (SDM_POWAIC(n).GT.0) i_bsc_sed=i_bsc_sed+1
        jpowaic(n)=i_bsc_sed*min(1,SDM_POWAIC(n))
        IF (SDM_POWAAL(n).GT.0) i_bsc_sed=i_bsc_sed+1
        jpowaal(n)=i_bsc_sed*min(1,SDM_POWAAL(n))
        IF (SDM_POWAPH(n).GT.0) i_bsc_sed=i_bsc_sed+1
        jpowaph(n)=i_bsc_sed*min(1,SDM_POWAPH(n))
        IF (SDM_POWAOX(n).GT.0) i_bsc_sed=i_bsc_sed+1
        jpowaox(n)=i_bsc_sed*min(1,SDM_POWAOX(n))
        IF (SDM_POWN2(n) .GT.0) i_bsc_sed=i_bsc_sed+1
        jpown2(n) =i_bsc_sed*min(1,SDM_POWN2(n))
        IF (SDM_POWNO3(n).GT.0) i_bsc_sed=i_bsc_sed+1
        jpowno3(n)=i_bsc_sed*min(1,SDM_POWNO3(n))
        IF (SDM_POWASI(n).GT.0) i_bsc_sed=i_bsc_sed+1
        jpowasi(n)=i_bsc_sed*min(1,SDM_POWASI(n))
        IF (SDM_SSSO12(n).GT.0) i_bsc_sed=i_bsc_sed+1
        jssso12(n)=i_bsc_sed*min(1,SDM_SSSO12(n))
        IF (SDM_SSSSIL(n).GT.0) i_bsc_sed=i_bsc_sed+1
        jssssil(n)=i_bsc_sed*min(1,SDM_SSSSIL(n))
        IF (SDM_SSSC12(n).GT.0) i_bsc_sed=i_bsc_sed+1
        jsssc12(n)=i_bsc_sed*min(1,SDM_SSSC12(n))
        IF (SDM_SSSTER(n).GT.0) i_bsc_sed=i_bsc_sed+1
        jssster(n)=i_bsc_sed*min(1,SDM_SSSTER(n))
      ENDDO

      DO n=1,nbgc
        IF (BUR_SSSO12(n).GT.0) i_bsc_bur=i_bsc_bur+1
        jburssso12(n)=i_bsc_bur*min(1,BUR_SSSO12(n))
        IF (BUR_SSSC12(n).GT.0) i_bsc_bur=i_bsc_bur+1
        jbursssc12(n)=i_bsc_bur*min(1,BUR_SSSC12(n))
        IF (BUR_SSSSIL(n).GT.0) i_bsc_bur=i_bsc_bur+1
        jburssssil(n)=i_bsc_bur*min(1,BUR_SSSSIL(n))
        IF (BUR_SSSTER(n).GT.0) i_bsc_bur=i_bsc_bur+1
        jburssster(n)=i_bsc_bur*min(1,BUR_SSSTER(n))
      ENDDO
#endif
         
      nbgcm2d    = i_bsc_m2d+i_atm_m2d
      nbgcm3d    = i_bsc_m3d
      nbgcm3dlvl = ilvl_bsc_m3d
      nbgct_sed  = i_bsc_sed
      nbgct_bur  = i_bsc_bur

!     Allocate buffers 

      IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)' '
        WRITE(io_stdo_bgc,*)'***************************************************'
        WRITE(io_stdo_bgc,*)'Memory allocation for averaging model output :'
        WRITE(io_stdo_bgc,*)' '
      ENDIF


      IF (mnproc.EQ.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgct2d ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',nbgct2d
      ENDIF

      ALLOCATE (bgct2d(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,nbgct2d),      &
     &  stat=errstat)
      IF (errstat.NE.0) STOP 'not enough memory bgct2d'
      IF (nbgct2d.NE.0) bgct2d=0.     
 
      IF (mnproc.EQ.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgcm2d ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',nbgcm2d
      ENDIF

      ALLOCATE (bgcm2d(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,nbgcm2d),      &
     &  stat=errstat)
      IF (errstat.NE.0) STOP 'not enough memory bgcm2d'
      IF (nbgcm2d.NE.0) bgcm2d=0.

      IF (mnproc.EQ.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgcm3d ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        WRITE(io_stdo_bgc,*)'Forth dimension    : ',nbgcm3d
      ENDIF

      ALLOCATE (bgcm3d(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,kpke,nbgcm3d), &
     &  stat=errstat)
      IF (errstat.NE.0) STOP 'not enough memory bgcm3d'
      IF (nbgcm3d.NE.0) bgcm3d=0.

      IF (mnproc.EQ.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgcm3dlvl '
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        WRITE(io_stdo_bgc,*)'Forth dimension    : ',nbgcm3dlvl
      ENDIF

      ALLOCATE (bgcm3dlvl(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,ddm,        &
     &  nbgcm3dlvl),stat=errstat)
      IF (errstat.NE.0) STOP 'not enough memory bgcm3dlvl'
      IF (nbgcm3dlvl.NE.0) bgcm3dlvl=0.

#ifndef sedbypass
      IF (mnproc.EQ.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgctsed ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',ks
        WRITE(io_stdo_bgc,*)'Forth dimension    : ',nbgct_sed
      ENDIF

      ALLOCATE (bgct_sed(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,ks,          &
     &  nbgct_sed),stat=errstat)
      IF (errstat.NE.0) STOP 'not enough memory bgct_sed'
      IF (nbgct_sed.NE.0) bgct_sed=0. 

      IF (mnproc.EQ.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgctbur ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',nbgct_bur
      ENDIF

      ALLOCATE (bgct_bur(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,             &
     &  nbgct_bur),stat=errstat)
      IF (errstat.NE.0) STOP 'not enough memory bgct_sed'
      IF (nbgct_bur.NE.0) bgct_bur=0. 
#endif

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
      IMPLICIT NONE
! 
      INTEGER :: pos
      REAL ::inival
! 
      INTEGER :: i,j,l
!
! --- Check whether field should be initialised
      IF (pos.EQ.0) RETURN
!
!$OMP PARALLEL DO
      DO j=1,jj
        DO l=1,isp(j)
          DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            bgcm2d(i,j,pos)=inival
          ENDDO
        ENDDO
      ENDDO
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
      IMPLICIT NONE
! 
      INTEGER :: pos
      REAL ::inival
! 
      INTEGER :: i,j,k,l
!
! --- Check whether field should be initialised
      IF (pos.EQ.0) RETURN
!
      DO k=1,kdm
!$OMP PARALLEL DO
        DO j=1,jj
          DO l=1,isp(j)
            DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              bgcm3d(i,j,k,pos)=inival
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDDO
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
      IMPLICIT NONE
! 
      INTEGER :: pos
      REAL ::inival
! 
      INTEGER :: i,j,k,l
!
! --- Check whether field should be initialised
      IF (pos.EQ.0) RETURN
!
      DO k=1,ddm
!$OMP PARALLEL DO
        DO j=1,jj
          DO l=1,isp(j)
            DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              bgcm3dlvl(i,j,k,pos)=inival
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDDO
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
      IMPLICIT NONE
! 
      INTEGER :: pos
      REAL ::inival
! 
      INTEGER :: i,j,k,l
!
! --- Check whether field should be initialised
      IF (pos.EQ.0) RETURN
!
      DO k=1,ks
!$OMP PARALLEL DO
        DO j=1,jj
          DO l=1,isp(j)
            DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              bgct_sed(i,j,k,pos)=inival
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDDO
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
      IMPLICIT NONE
! 
      INTEGER :: pos
      REAL ::inival
! 
      INTEGER :: i,j,l
!
! --- Check whether field should be initialised
      IF (pos.EQ.0) RETURN
!
!$OMP PARALLEL DO
      DO j=1,jj
        DO l=1,isp(j)
          DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            bgct_bur(i,j,pos)=inival
          ENDDO
        ENDDO
      ENDDO
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
      IMPLICIT NONE
! 
      INTEGER :: pos(nbgcmax),wghtsflg
      REAL, DIMENSION(idm,jdm) :: fld,wghts
! 
      INTEGER :: i,j,l,o
!
! --- Check whether field should be accumulated
      DO o=1,nbgc
        IF (pos(o).EQ.0) cycle
!
          IF (wghtsflg.eq.0) then 
!$OMP PARALLEL DO 
            DO j=1,jj
              DO l=1,isp(j)
                DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                  bgcm2d(i,j,pos(o))=bgcm2d(i,j,pos(o))+fld(i,j)
                ENDDO
              ENDDO
            ENDDO
!$OMP END PARALLEL DO
          ELSE
!$OMP PARALLEL DO 
            DO j=1,jj
              DO l=1,isp(j)
                DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                  bgcm2d(i,j,pos(o))=bgcm2d(i,j,pos(o))+fld(i,j)*       &
     &              wghts(i,j)
                ENDDO
              ENDDO
            ENDDO
!$OMP END PARALLEL DO
          ENDIF 
!
      ENDDO
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
      IMPLICIT NONE
! 
      INTEGER :: pos(nbgcmax),wghtsflg
      REAL, DIMENSION(idm,jdm,kdm) :: fld,wghts
! 
      INTEGER :: i,j,k,l,o
!
! --- Check whether field should be accumulated
      DO o=1,nbgc
        IF (pos(o).EQ.0) cycle
!
          IF (wghtsflg.eq.0) then 
            DO k=1,kk
!$OMP PARALLEL DO 
              DO j=1,jj
                DO l=1,isp(j)
                  DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                    bgcm3d(i,j,k,pos(o))=bgcm3d(i,j,k,pos(o))+          &
     &                fld(i,j,k)
                  ENDDO
                ENDDO
              ENDDO
!$OMP END PARALLEL DO
            ENDDO
          ELSE
            DO k=1,kk
!$OMP PARALLEL DO 
              DO j=1,jj
                DO l=1,isp(j)
                  DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                    bgcm3d(i,j,k,pos(o))=bgcm3d(i,j,k,pos(o))+          &
     &                fld(i,j,k)*wghts(i,j,k)
                  ENDDO
                ENDDO
              ENDDO
!$OMP END PARALLEL DO
            ENDDO
          ENDIF
!
      ENDDO
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
      IMPLICIT NONE
! 
      INTEGER :: pos(nbgcmax),k
      INTEGER, DIMENSION(idm,jdm) :: ind1,ind2
      REAL, DIMENSION(idm,jdm,ddm) :: wghts
      REAL, DIMENSION(idm,jdm,kdm) :: fld
! 
      INTEGER :: d,i,j,l,o
!
! --- Check whether field should be accumulated
      DO o=1,nbgc
        IF (pos(o).EQ.0) cycle
!
!$OMP PARALLEL DO
        DO j=1,jj
          DO l=1,isp(j)
            DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              DO d=ind1(i,j),ind2(i,j)
                bgcm3dlvl(i,j,d,pos(o))=bgcm3dlvl(i,j,d,pos(o))+        &
     &            fld(i,j,k)*wghts(i,j,d)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDDO
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
      IMPLICIT NONE
! 
      INTEGER :: pos(nbgcmax)
      REAL, DIMENSION(idm,jdm,ks) :: fld
! 
      INTEGER :: i,j,k,l,o
!
! --- Check whether field should be accumulated
      DO o=1,nbgc
        IF (pos(o).EQ.0) cycle
!
        DO k=1,ks
!$OMP PARALLEL DO 
          DO j=1,jj
            DO l=1,isp(j)
              DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                bgct_sed(i,j,k,pos(o))=bgct_sed(i,j,k,pos(o))+fld(i,j,k)
              ENDDO
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ENDDO
      ENDDO
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
      IMPLICIT NONE
! 
      INTEGER :: pos(nbgcmax)
      REAL, DIMENSION(idm,jdm) :: fld
! 
      INTEGER :: i,j,l,o
!
! --- Check whether field should be accumulated
      DO o=1,nbgc
        IF (pos(o).EQ.0) cycle
!
!$OMP PARALLEL DO 
        DO j=1,jj
          DO l=1,isp(j)
            DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              bgct_bur(i,j,pos(o))=bgct_bur(i,j,pos(o))+fld(i,j)
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDDO
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
      IMPLICIT NONE
! 
      INTEGER :: posacc,poswgt
! 
      INTEGER :: i,j,l
      REAL, parameter :: epsil=1e-11
!
! --- Check whether field should be initialised
      IF (posacc.EQ.0) RETURN
!
!$OMP PARALLEL DO
      DO j=1,jj
        DO l=1,isp(j)
          DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            bgcm2d(i,j,posacc)=bgcm2d(i,j,posacc)/                      &
     &        max(epsil,bgcm2d(i,j,poswgt))
          ENDDO
        ENDDO
      ENDDO
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
      IMPLICIT NONE
! 
      INTEGER :: posacc,poswgt
! 
      INTEGER :: i,j,k,l
      REAL, parameter :: epsil=1e-11
!
! --- Check whether field should be initialised
      IF (posacc.EQ.0) RETURN
!
      DO k=1,kk
!$OMP PARALLEL DO
        DO j=1,jj
          DO l=1,isp(j)
            DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              IF (bgcm3d(i,j,k,poswgt).GT.epsil) THEN
                bgcm3d(i,j,k,posacc)=bgcm3d(i,j,k,posacc)/              &
     &            bgcm3d(i,j,k,poswgt) 
              ELSE 
                bgcm3d(i,j,k,posacc)=nf90_fill_double
              ENDIF
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDDO
!     
      END SUBROUTINE finlyr



      SUBROUTINE wrtsrf(pos,frmt,sfac,offs,cmpflg,vnm,vlngnm,vstdnm,    &
     &  vunits)
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
! ---                            written IF flag is set to 1 
! ---   char vnm      (in)     : variable name used in nc-file 
! ---   char vlngnm   (in)     : variable long name (skipped IF ' ') 
! ---   char vstdnm   (in)     : variable standard name (skipped IF ' ') 
! ---   char vunits   (in)     : variable units (skipped IF ' ') 
! --- ------------------------------------------------------------------
!
      IMPLICIT NONE
! 
      REAL ::sfac,offs
      INTEGER :: frmt,cmpflg,pos,n
      CHARACTER(LEN=*) :: vnm,vlngnm,vstdnm,vunits
!
      CHARACTER(LEN=100) :: dims
!
! --- Check whether field should be written
      IF (pos.EQ.0 .OR. frmt.EQ.0) RETURN
!
! --- Create dimension string 
      IF (cmpflg.EQ.1) THEN
        dims='pcomp time'
      ELSE
        dims='x y time'
      ENDIF
!
! --- Check output format
      IF (frmt.EQ.2) THEN
        IF (cmpflg.EQ.1) THEN
          CALL nccopa(vnm,dims,bgcm2d(1-nbdy,1-nbdy,pos),ip,sfac,       &
     &      offs)
        ELSE
          CALL ncpack(vnm,dims,bgcm2d(1-nbdy,1-nbdy,pos),ip,2,          &
     &      sfac,offs)
        ENDIF
      ELSEIF (frmt.EQ.4) THEN
        IF (cmpflg.EQ.1) THEN
          CALL nccomp(vnm,dims,bgcm2d(1-nbdy,1-nbdy,pos),ip,sfac,       &
     &      offs,4)
        ELSE
          CALL ncwrtr(vnm,dims,bgcm2d(1-nbdy,1-nbdy,pos),ip,1,          &
     &      sfac,offs,4)
        ENDIF
      ELSEIF (frmt.EQ.8) THEN
        IF (cmpflg.EQ.1) THEN
          CALL nccomp(vnm,dims,bgcm2d(1-nbdy,1-nbdy,pos),ip,sfac,       &
     &      offs,8)
        ELSE
          CALL ncwrtr(vnm,dims,bgcm2d(1-nbdy,1-nbdy,pos),ip,1,          &
     &      sfac,offs,8)
        ENDIF
      ELSE
        STOP 'unknown output format '
      ENDIF
!
! --- Def.NE.attributes
!      IF (len(trim(vunits)).NE.0) CALL ncattr('units',vunits)
!      IF (len(trim(vlngnm)).NE.0) CALL ncattr('long_name',vlngnm)
!      IF (len(trim(vstdnm)).NE.0) CALL ncattr('standard_name',vstdnm)
!      CALL ncattr('coordinates','plon plat')
!      CALL ncattr('cell_measures','area: parea')
!
      END SUBROUTINE wrtsrf



      SUBROUTINE wrtlyr(pos,frmt,sfac,offs,cmpflg,vnm,vlngnm,vstdnm,    &
     &  vunits)
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
! ---                            written IF flag is set to 1 
! ---   char vnm      (in)     : variable name used in nc-file 
! ---   char vlngnm   (in)     : variable long name (skipped IF ' ') 
! ---   char vstdnm   (in)     : variable standard name (skipped IF ' ') 
! ---   char vunits   (in)     : variable units (skipped IF ' ') 
! --- ------------------------------------------------------------------
!
      IMPLICIT NONE
! 
      REAL ::sfac,offs
      INTEGER :: frmt,cmpflg,pos,n
      CHARACTER(LEN=*) :: vnm,vlngnm,vstdnm,vunits
!
      CHARACTER(LEN=100) :: dims
!
! --- Check whether field should be written
      IF (pos.EQ.0 .OR. frmt.EQ.0) RETURN
!
! --- Create dimension string 
      IF (cmpflg.EQ.1) THEN
        dims='pcomp sigma time'
      ELSE
        dims='x y sigma time'
      ENDIF
!
! --- Check output format
      IF (frmt.EQ.2) THEN
        IF (cmpflg.EQ.1) THEN
          CALL nccopa(vnm,dims,bgcm3d(1-nbdy,1-nbdy,1,pos),ip,sfac,     &
     &      offs)
        ELSE
          CALL ncpack(vnm,dims,bgcm3d(1-nbdy,1-nbdy,1,pos),ip,2,        &
     &      sfac,offs)
        ENDIF
      ELSEIF (frmt.EQ.4) THEN
        IF (cmpflg.EQ.1) THEN
          CALL nccomp(vnm,dims,bgcm3d(1-nbdy,1-nbdy,1,pos),ip,sfac,     &
     &      offs,4)
        ELSE
          CALL ncwrtr(vnm,dims,bgcm3d(1-nbdy,1-nbdy,1,pos),ip,2,        &
     &      sfac,offs,4)
        ENDIF
      ELSEIF (frmt.EQ.8) THEN
        IF (cmpflg.EQ.1) THEN
          CALL nccomp(vnm,dims,bgcm3d(1-nbdy,1-nbdy,1,pos),ip,sfac,     &
     &      offs,8)
        ELSE
          CALL ncwrtr(vnm,dims,bgcm3d(1-nbdy,1-nbdy,1,pos),ip,2,        &
     &      sfac,offs,8)
        ENDIF
      ELSE
        STOP 'unknown output format '
      ENDIF
!
! --- Def.NE.attributes
!      IF (len(trim(vunits)).NE.0) CALL ncattr('units',vunits)
!      IF (len(trim(vlngnm)).NE.0) CALL ncattr('long_name',vlngnm)
!      IF (len(trim(vstdnm)).NE.0) CALL ncattr('standard_name',vstdnm)
!      CALL ncattr('coordinates','plon plat')
!      CALL ncattr('cell_measures','area: parea')
!
      END SUBROUTINE wrtlyr



      SUBROUTINE wrtlvl(pos,frmt,sfac,offs,cmpflg,vnm,vlngnm,vstdnm,    &
     &  vunits)
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
! ---                            written IF flag is set to 1 
! ---   char vnm      (in)     : variable name used in nc-file 
! ---   char vlngnm   (in)     : variable long name (skipped IF ' ') 
! ---   char vstdnm   (in)     : variable standard name (skipped IF ' ') 
! ---   char vunits   (in)     : variable units (skipped IF ' ') 
! --- ------------------------------------------------------------------
!
      IMPLICIT NONE
! 
      REAL ::sfac,offs
      INTEGER :: frmt,cmpflg,pos,n
      CHARACTER(LEN=*) :: vnm,vlngnm,vstdnm,vunits
!
      CHARACTER(LEN=100) :: dims
!
! --- Check whether field should be written
      IF (pos.EQ.0 .OR. frmt.EQ.0) RETURN
!
! --- Create dimension string 
      IF (cmpflg.EQ.1) THEN
        dims='pcomp depth time'
      ELSE
        dims='x y depth time'
      ENDIF
!
! --- Check output format
      IF (frmt.EQ.2) THEN
        IF (cmpflg.EQ.1) THEN
          CALL nccopa(vnm,dims,bgcm3dlvl(1-nbdy,1-nbdy,1,pos),ip,sfac,  &
     &      offs)
        ELSE
          CALL ncpack(vnm,dims,bgcm3dlvl(1-nbdy,1-nbdy,1,pos),ip,2,     &
     &      sfac,offs)
        ENDIF
      ELSEIF (frmt.EQ.4) THEN
        IF (cmpflg.EQ.1) THEN
          CALL nccomp(vnm,dims,bgcm3dlvl(1-nbdy,1-nbdy,1,pos),ip,sfac,  &
     &      offs,4)
        ELSE
          CALL ncwrtr(vnm,dims,bgcm3dlvl(1-nbdy,1-nbdy,1,pos),ip,2,     &
     &      sfac,offs,4)
        ENDIF
      ELSEIF (frmt.EQ.8) THEN
        IF (cmpflg.EQ.1) THEN
          CALL nccomp(vnm,dims,bgcm3dlvl(1-nbdy,1-nbdy,1,pos),ip,sfac,  &
     &      offs,8)
        ELSE
          CALL ncwrtr(vnm,dims,bgcm3dlvl(1-nbdy,1-nbdy,1,pos),ip,2,     &
     &      sfac,offs,8)
        ENDIF
      ELSE
        STOP 'unknown output format '
      ENDIF
!
! --- Def.NE.attributes
!      IF (len(trim(vunits)).NE.0) CALL ncattr('units',vunits)
!      IF (len(trim(vlngnm)).NE.0) CALL ncattr('long_name',vlngnm)
!      IF (len(trim(vstdnm)).NE.0) CALL ncattr('standard_name',vstdnm)
!      CALL ncattr('coordinates','plon plat')
!      CALL ncattr('cell_measures','area: parea')
!
      END SUBROUTINE wrtlvl



      SUBROUTINE wrtsdm(pos,frmt,sfac,offs,cmpflg,vnm,vlngnm,vstdnm,    &
     &  vunits)
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
! ---                            written IF flag is set to 1 
! ---   char vnm      (in)     : variable name used in nc-file 
! ---   char vlngnm   (in)     : variable long name (skipped IF ' ') 
! ---   char vstdnm   (in)     : variable standard name (skipped IF ' ') 
! ---   char vunits   (in)     : variable units (skipped IF ' ') 
! --- ------------------------------------------------------------------
!
      IMPLICIT NONE
! 
      REAL ::sfac,offs
      INTEGER :: frmt,cmpflg,pos,n
      CHARACTER(LEN=*) :: vnm,vlngnm,vstdnm,vunits
!
      CHARACTER(LEN=100) :: dims
!
! --- Check whether field should be written
      IF (pos.EQ.0 .OR. frmt.EQ.0) RETURN
!
! --- Create dimension string 
      IF (cmpflg.EQ.1) THEN
        dims='pcomp ks time'
      ELSE
        dims='x y ks time'
      ENDIF
!
! --- Check output format
      IF (frmt.EQ.2) THEN
        IF (cmpflg.EQ.1) THEN
          CALL nccopa(vnm,dims,bgct_sed(1-nbdy,1-nbdy,1,pos),ip,sfac,   &
     &      offs)
        ELSE
          CALL ncpack(vnm,dims,bgct_sed(1-nbdy,1-nbdy,1,pos),ip,1,      &
     &      sfac,offs)
        ENDIF
      ELSEIF (frmt.EQ.4) THEN
        IF (cmpflg.EQ.1) THEN
          CALL nccomp(vnm,dims,bgct_sed(1-nbdy,1-nbdy,1,pos),ip,sfac,   &
     &      offs,4)
        ELSE
          CALL ncwrtr(vnm,dims,bgct_sed(1-nbdy,1-nbdy,1,pos),ip,1,      &
     &      sfac,offs,4)
        ENDIF
      ELSEIF (frmt.EQ.8) THEN
        IF (cmpflg.EQ.1) THEN
          CALL nccomp(vnm,dims,bgct_sed(1-nbdy,1-nbdy,1,pos),ip,sfac,   &
     &      offs,8)
        ELSE
          CALL ncwrtr(vnm,dims,bgct_sed(1-nbdy,1-nbdy,1,pos),ip,1,      &
     &      sfac,offs,8)
        ENDIF
      ELSE
        STOP 'unknown output format '
      ENDIF
!
! --- Def.NE.attributes
!      IF (len(trim(vunits)).NE.0) CALL ncattr('units',vunits)
!      IF (len(trim(vlngnm)).NE.0) CALL ncattr('long_name',vlngnm)
!      IF (len(trim(vstdnm)).NE.0) CALL ncattr('standard_name',vstdnm)
!      CALL ncattr('coordinates','plon plat')
!      CALL ncattr('cell_measures','area: parea')
!
      END SUBROUTINE wrtsdm



      SUBROUTINE wrtbur(pos,frmt,sfac,offs,cmpflg,vnm,vlngnm,vstdnm,    &
     &  vunits)
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
! ---                            written IF flag is set to 1 
! ---   char vnm      (in)     : variable name used in nc-file 
! ---   char vlngnm   (in)     : variable long name (skipped IF ' ') 
! ---   char vstdnm   (in)     : variable standard name (skipped IF ' ') 
! ---   char vunits   (in)     : variable units (skipped IF ' ') 
! --- ------------------------------------------------------------------
!
      IMPLICIT NONE
! 
      REAL ::sfac,offs
      INTEGER :: frmt,cmpflg,pos,n
      CHARACTER(LEN=*) :: vnm,vlngnm,vstdnm,vunits
!
      CHARACTER(LEN=100) :: dims
!
! --- Check whether field should be written
      IF (pos.EQ.0 .OR. frmt.EQ.0) RETURN
!
! --- Create dimension string 
      IF (cmpflg.EQ.1) THEN
        dims='pcomp time'
      ELSE
        dims='x y time'
      ENDIF
!
! --- Check output format
      IF (frmt.EQ.2) THEN
        IF (cmpflg.EQ.1) THEN
          CALL nccopa(vnm,dims,bgct_bur(1-nbdy,1-nbdy,pos),ip,sfac,     &
     &      offs)
        ELSE
          CALL ncpack(vnm,dims,bgct_bur(1-nbdy,1-nbdy,pos),ip,1,        &
     &      sfac,offs)
        ENDIF
      ELSEIF (frmt.EQ.4) THEN
        IF (cmpflg.EQ.1) THEN
          CALL nccomp(vnm,dims,bgct_bur(1-nbdy,1-nbdy,pos),ip,sfac,     &
     &      offs,4)
        ELSE
          CALL ncwrtr(vnm,dims,bgct_bur(1-nbdy,1-nbdy,pos),ip,1,        &
     &      sfac,offs,4)
        ENDIF
      ELSEIF (frmt.EQ.8) THEN
        IF (cmpflg.EQ.1) THEN
          CALL nccomp(vnm,dims,bgct_bur(1-nbdy,1-nbdy,pos),ip,sfac,     &
     &      offs,8)
        ELSE
          CALL ncwrtr(vnm,dims,bgct_bur(1-nbdy,1-nbdy,pos),ip,1,        &
     &      sfac,offs,8)
        ENDIF
      ELSE
        STOP 'unknown output format '
      ENDIF
!
! --- Def.NE.attributes
!      IF (len(trim(vunits)).NE.0) CALL ncattr('units',vunits)
!      IF (len(trim(vlngnm)).NE.0) CALL ncattr('long_name',vlngnm)
!      IF (len(trim(vstdnm)).NE.0) CALL ncattr('standard_name',vstdnm)
!      CALL ncattr('coordinates','plon plat')
!      CALL ncattr('cell_measures','area: parea')
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
      IMPLICIT NONE
! 
      REAL ::sfac,offs
      INTEGER :: pos
! 
      INTEGER :: i,j,l
      REAL ::epsil=1e-11
!
! --- Check whether field should be processed
      IF (pos.EQ.0) RETURN
!
!$OMP PARALLEL DO 
      DO j=1,jj
        DO l=1,isp(j)
          DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            IF (bgcm2d(i,j,pos).LT.epsil) THEN
              bgcm2d(i,j,pos)=0.
            ELSE
              bgcm2d(i,j,pos)=log10(bgcm2d(i,j,pos)*sfac+offs)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
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
      IMPLICIT NONE
! 
      REAL ::sfac,offs
      INTEGER :: pos
! 
      INTEGER :: i,j,k,l
      REAL ::epsil=1e-11
!
! --- Check whether field should be processed
      IF (pos.EQ.0) RETURN
!
      DO k=1,kdm
!$OMP PARALLEL DO 
        DO j=1,jj
          DO l=1,isp(j)
            DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              IF (bgcm3d(i,j,k,pos).LT.epsil) THEN
                bgcm3d(i,j,k,pos)=0.
              ELSEIF (bgcm3d(i,j,k,pos).NE.nf90_fill_double) THEN
                bgcm3d(i,j,k,pos)=log10(bgcm3d(i,j,k,pos)*sfac+offs)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDDO 
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
      IMPLICIT NONE
! 
      REAL ::sfac,offs
      INTEGER :: pos
! 
      INTEGER :: i,j,k,l
      REAL ::epsil=1e-11
!
! --- Check whether field should be processed
      IF (pos.EQ.0) RETURN
!
      DO k=1,ddm
!$OMP PARALLEL DO 
        DO j=1,jj
          DO l=1,isp(j)
            DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              IF (bgcm3dlvl(i,j,k,pos).LT.epsil) THEN
                bgcm3dlvl(i,j,k,pos)=0.
              ELSEIF (bgcm3dlvl(i,j,k,pos).NE.nf90_fill_double) THEN
                bgcm3dlvl(i,j,k,pos)=log10(bgcm3dlvl(i,j,k,pos)*sfac+   &
     &            offs)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDDO
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
      IMPLICIT NONE
! 
      REAL ::sfac,offs
      INTEGER :: pos
! 
      INTEGER :: i,j,k,l
      REAL ::epsil=1e-11
!
! --- Check whether field should be processed
      IF (pos.EQ.0) RETURN
!
      DO k=1,ks
!$OMP PARALLEL DO 
        DO j=1,jj
          DO l=1,isp(j)
            DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              IF (bgct_sed(i,j,k,pos).LT.epsil) THEN
                bgct_sed(i,j,k,pos)=0.
              ELSE
                bgct_sed(i,j,k,pos)=log10(bgct_sed(i,j,k,pos)*sfac+offs)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDDO
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
      IMPLICIT NONE
!
      INTEGER :: pos
      INTEGER, DIMENSION(idm,jdm) :: idepth 
! 
      INTEGER :: i,j,l
      REAL, parameter :: mskval=nf90_fill_double
!
! --- Check whether field should be initia
      IF (pos.EQ.0) RETURN
!
!$OMP PARALLEL DO
      DO j=1,jj
        DO l=1,isp(j)
          DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            if( idepth(i,j) <= 0 ) bgcm2d(i,j,pos)=mskval
          ENDDO
        ENDDO
      ENDDO
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
      IMPLICIT NONE
!
      INTEGER :: pos
      REAL, DIMENSION(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: depths 
! 
      INTEGER :: i,j,k,l
      LOGICAL :: iniflg=.true.
      INTEGER, DIMENSION(idm,jdm) :: kmax
      REAL, parameter :: mskval=nf90_fill_double
!
      SAVE iniflg,kmax
!     
! --- Check whether field should be processed
      IF (pos.EQ.0) RETURN
!
! --- Prepare index fields for masking

      IF (iniflg) THEN
!$OMP PARALLEL DO
        DO j=1,jj
          DO i=1,ii
            kmax(i,j)=0
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
        DO k=1,ddm
!$OMP PARALLEL DO
          DO j=1,jj
            DO i=1,ii
              IF (depths(i,j).GT.depthslev_bnds(1,k)) kmax(i,j)=k
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ENDDO
        iniflg=.false.
      ENDIF
!
!$OMP PARALLEL DO
      DO j=1,jj
        DO i=1,ii
          DO k=kmax(i,j)+1,ddm
            bgcm3dlvl(i,j,k,pos)=mskval
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
!
      END SUBROUTINE msklvl



      SUBROUTINE bgczlv(pddpo,kin,ind1,ind2,weights)
!-----------------------------------------------------------------------
!
      USE mod_xc
!
      IMPLICIT NONE
!
      INTEGER :: d,i,j,k,l,kin
      INTEGER, DIMENSION(idm,jdm) :: ind1,ind2
! 
      REAL, PARAMETER :: eps=1e-10
      REAL, DIMENSION(idm,jdm,kdm) :: pddpo,ztop,zbot
      REAL, DIMENSION(idm,jdm,ddm) :: weights,dlev
!
      LOGICAL :: iniflg=.true.
!
      SAVE ztop,zbot,dlev,iniflg
!
! --- Adjust bounds of levitus levels according to model bathymetry
      IF (iniflg) THEN
        DO d=1,ddm
!$OMP PARALLEL DO  
          DO j=1,jj
            DO l=1,isp(j)
              DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                 dlev(i,j,d)=max(eps,min(pbath(i,j),                    &
     &             depthslev_bnds(2,d))-depthslev_bnds(1,d))
              ENDDO
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ENDDO
        iniflg=.false.
      ENDIF
!
! --- Compute top and bottom depths of density layers 
      IF (kin.EQ.1) THEN       
!$OMP PARALLEL DO  
        DO j=1,jj
          DO l=1,isp(j)
            DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              zbot(i,j,1)=pddpo(i,j,1)
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
        DO k=2,kk
!$OMP PARALLEL DO  
          DO j=1,jj
            DO l=1,isp(j)
              DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                zbot(i,j,k)=zbot(i,j,k-1)+pddpo(i,j,k)
              ENDDO
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ENDDO
!$OMP PARALLEL DO  
        DO j=1,jj
          DO l=1,isp(j)
            DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              zbot(i,j,1)=zbot(i,j,1)*pbath(i,j)/zbot(i,j,kk)
              ztop(i,j,1)=0.
              ind1(i,j)=1
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
        DO k=2,kk
!$OMP PARALLEL DO  
          DO j=1,jj
            DO l=1,isp(j)
              DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
                zbot(i,j,k)=zbot(i,j,k)*pbath(i,j)/zbot(i,j,kk)
                ztop(i,j,k)=zbot(i,j,k-1)
              ENDDO
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ENDDO
      ENDIF
!
! --- Compute interpolation weights 
!$OMP PARALLEL DO 
      DO j=1,jj
        DO l=1,isp(j)
          DO i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            ind2(i,j)=0
            IF (pddpo(i,j,kin).GT.eps) THEN
              DO d=ind1(i,j),ddm
                IF (depthslev_bnds(2,d).LE.ztop(i,j,kin)) THEN
                  ind1(i,j)=d+1
                  CYCLE
                ELSEIF (depthslev_bnds(1,d).GE.zbot(i,j,kin)) THEN
                  EXIT
                ENDIF
                ind2(i,j)=d
                weights(i,j,d)=(min(zbot(i,j,kin),                      &
     &            depthslev_bnds(2,d))-max(ztop(i,j,kin),               &
     &            depthslev_bnds(1,d)))/dlev(i,j,d)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
!
      END SUBROUTINE bgczlv

      END MODULE mo_bgcmean
