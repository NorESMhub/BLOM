      MODULE mo_bgcmean
!***********************************************************************
!
!**** *MODULE mo_bgcmean* - Variables for bgcmean.
!
!     Ingo Bethke       *Bjer.NE. C.*    05.11.09 
!     Patrick Wetzel    *MPI-Met, HH*    09.12.02
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
     & SRF_KWCO2         ,SRF_PCO2          ,SRF_DMSFLUX       ,        &
     & SRF_CO2FXD        ,SRF_CO2FXU        ,SRF_OXFLUX        ,        &
     & SRF_NIFLUX        ,SRF_DMS           ,SRF_DMSPROD       ,        &
     & SRF_DMS_BAC       ,SRF_DMS_UV        ,SRF_EXPORT        ,        &
     & SRF_EXPOSI        ,SRF_EXPOCA        ,SRF_ATMCO2        ,        &
     & SRF_ATMO2         ,SRF_ATMN2         ,                           &
     & LYR_PHYTO         ,LYR_GRAZER        ,LYR_DOC           ,        &
     & LYR_PHOSY         ,LYR_PHOSPH        ,LYR_OXYGEN        ,        &
     & LYR_IRON          ,LYR_ANO3          ,LYR_ALKALI        ,        &
     & LYR_SILICA        ,LYR_DIC           ,LYR_POC           ,        &
     & LYR_CALC          ,LYR_OPAL          ,LYR_CO3           ,        &
     & LYR_PH            ,LYR_OMEGAC        ,LYR_DIC13         ,        &
     & LYR_DIC14         ,LYR_DP            ,LYR_NOS           ,        &
     & LVL_PHYTO         ,LVL_GRAZER        ,LVL_DOC           ,        &
     & LVL_PHOSY         ,LVL_PHOSPH        ,LVL_OXYGEN        ,        &
     & LVL_IRON          ,LVL_ANO3          ,LVL_ALKALI        ,        &
     & LVL_SILICA        ,LVL_DIC           ,LVL_POC           ,        &
     & LVL_CALC          ,LVL_OPAL          ,LVL_CO3           ,        &
     & LVL_PH            ,LVL_OMEGAC        ,LVL_DIC13         ,        &
     & LVL_DIC14         ,LVL_DP            ,LVL_NOS           ,        &
     & SDM_POWAIC        ,SDM_POWAAL        ,SDM_POWAPH        ,        &
     & SDM_POWAOX        ,SDM_POWN2         ,SDM_POWNO3        ,        &
     & SDM_POWASI        ,SDM_SSSO12        ,SDM_SSSSIL        ,        &
     & SDM_SSSC12        ,SDM_SSSTER        ,                           &
     & GLB_AVEPERIO      ,GLB_FILEFREQ      ,GLB_COMPFLAG      ,        &
     & GLB_NCFORMAT      ,GLB_INVENTORY 
      CHARACTER(LEN=10), DIMENSION(nbgcmax), SAVE :: GLB_FNAMETAG
      namelist /DIABGC/                                                 &
     & SRF_KWCO2         ,SRF_PCO2          ,SRF_DMSFLUX       ,        &
     & SRF_CO2FXD        ,SRF_CO2FXU        ,SRF_OXFLUX        ,        &
     & SRF_NIFLUX        ,SRF_DMS           ,SRF_DMSPROD       ,        &
     & SRF_DMS_BAC       ,SRF_DMS_UV        ,SRF_EXPORT        ,        &
     & SRF_EXPOSI        ,SRF_EXPOCA        ,SRF_ATMCO2        ,        &
     & SRF_ATMO2         ,SRF_ATMN2         ,                           &
     & LYR_PHYTO         ,LYR_GRAZER        ,LYR_DOC           ,        &
     & LYR_PHOSY         ,LYR_PHOSPH        ,LYR_OXYGEN        ,        &
     & LYR_IRON          ,LYR_ANO3          ,LYR_ALKALI        ,        &
     & LYR_SILICA        ,LYR_DIC           ,LYR_POC           ,        &
     & LYR_CALC          ,LYR_OPAL          ,LYR_CO3           ,        &
     & LYR_PH            ,LYR_OMEGAC        ,LYR_DIC13         ,        &
     & LYR_DIC14         ,LYR_DP            ,LYR_NOS           ,        &
     & LVL_PHYTO         ,LVL_GRAZER        ,LVL_DOC           ,        &
     & LVL_PHOSY         ,LVL_PHOSPH        ,LVL_OXYGEN        ,        &
     & LVL_IRON          ,LVL_ANO3          ,LVL_ALKALI        ,        &
     & LVL_SILICA        ,LVL_DIC           ,LVL_POC           ,        &
     & LVL_CALC          ,LVL_OPAL          ,LVL_CO3           ,        &
     & LVL_PH            ,LVL_OMEGAC        ,LVL_DIC13         ,        &
     & LVL_DIC14         ,LVL_NOS           ,                           &
     & SDM_POWAIC        ,SDM_POWAAL        ,SDM_POWAPH        ,        &
     & SDM_POWAOX        ,SDM_POWN2         ,SDM_POWNO3        ,        &
     & SDM_POWASI        ,SDM_SSSO12        ,SDM_SSSSIL        ,        &
     & SDM_SSSC12        ,SDM_SSSTER        ,                           &
     & GLB_AVEPERIO      ,GLB_FILEFREQ      ,GLB_COMPFLAG      ,        &
     & GLB_NCFORMAT      ,GLB_FNAMETAG      ,GLB_INVENTORY

!----------------------------------------------------------------      

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
     &          jkwco2     ,                                            &
     &          jpco2      ,                                            &
     &          jdmsflux   ,                                            &
     &          jco2fxd    ,                                            &
     &          jco2fxu    ,                                            &
     &          joxflux    ,                                            &
     &          jniflux    ,                                            &
     &          jdms       ,                                            &
     &          jdmsprod   ,                                            &
     &          jdms_bac   ,                                            &
     &          jdms_uv    ,                                            &
     &          jexport    ,                                            &
     &          jexpoca    ,                                            &
     &          jexposi     

      INTEGER, SAVE :: i_atm_m2d  
      INTEGER, DIMENSION(nbgcmax), SAVE ::                              &
     &          jatmco2  ,                                              &
     &          jatmo2   ,                                              &
     &          jatmn2           

      INTEGER, SAVE :: nbgcm2d 

!----------------------------------------------------------------
! 2d: jdms,jexport, jexpoca, jexposi
! out : jpoc, jatten, jcalc, jopal
      INTEGER, SAVE :: i_bsc_m3d,ilvl_bsc_m3d 
      INTEGER, DIMENSION(nbgcmax), SAVE ::                              &
     &          jdp     ,                                               &
     &          jphyto  ,                                               &
     &          jgrazer ,                                               &
     &          jdoc    ,                                               &
     &          jphosy  ,                                               &
     &          jphosph ,                                               &
     &          joxygen ,                                               &
     &          jiron   ,                                               &
     &          jano3   ,                                               &
     &          jalkali ,                                               &
     &          jsilica ,                                               &
     &          jdic    ,                                               &
     &          jpoc    ,                                               &
     &          jcalc   ,                                               &
     &          jopal   ,                                               &
     &          jco3    ,                                               &
     &          jph     ,                                               &
     &          jomegac ,                                               &
     &          jlvlphyto  ,                                            &
     &          jlvlgrazer ,                                            &
     &          jlvldoc    ,                                            &
     &          jlvlphosy  ,                                            &
     &          jlvlphosph ,                                            &
     &          jlvloxygen ,                                            &
     &          jlvliron   ,                                            &
     &          jlvlano3   ,                                            &
     &          jlvlalkali ,                                            &
     &          jlvlsilica ,                                            &
     &          jlvldic    ,                                            &
     &          jlvlpoc    ,                                            &
     &          jlvlcalc   ,                                            &
     &          jlvlopal   ,                                            &
     &          jlvlco3    ,                                            &
     &          jlvlph     ,                                            &
     &          jlvlomegac
    

      INTEGER, SAVE :: i_iso_m3d,ilvl_iso_m3d
      INTEGER, DIMENSION(nbgcmax), SAVE ::                              &
     &          jdic13  ,                                               &
     &          jdic14  ,                                               &
     &          jlvldic13  ,                                            &
     &          jlvldic14   

      INTEGER, SAVE :: i_ant_m3d,ilvl_ant_m3d 

      INTEGER, SAVE :: i_agg_m3d,ilvl_agg_m3d 
      INTEGER, DIMENSION(nbgcmax), SAVE ::                              &
     &          jnos ,                                                  &
     &          jlvlnos     

      INTEGER, SAVE :: i_cfc_t3d,ilvl_cfc_t3d
      INTEGER, DIMENSION(nbgcmax), SAVE ::                              &
     &          jcfc11_t ,                                              &
     &          jcfc12_t ,                                              &
     &          jac14_t  ,                                              &
     &          jlvlcfc11_t ,                                           &
     &          jlvlcfc12_t ,                                           &
     &          jlvlac14_t 

      INTEGER, SAVE :: nbgcm3d,nbgcm3dlvl 

!----------------------------------------------------------------
! js: sediment
!
      INTEGER, SAVE :: i_bsc_sed 
      INTEGER, DIMENSION(nbgcmax), SAVE ::                              &
     &          jpowaic  ,                                              &
     &          jpowaal  ,                                              &
     &          jpowaph  ,                                              &
     &          jpowaox  ,                                              &
     &          jpown2   ,                                              &
     &          jpowno3  ,                                              &
     &          jpowasi  ,                                              &
     &          jssso12  ,                                              &
     &          jssssil  ,                                              &
     &          jsssc12  ,                                              &
     &          jssster              


      INTEGER, SAVE :: nbgct_sed    

!----------------------------------------------------------------

      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: bgct2d
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: bgcm2d
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: bgcm3d,bgcm3dlvl
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: bgct_sed 
     
 
      CONTAINS


      SUBROUTINE ALLOC_MEM_BGCMEAN(kpie,kpje,kpke)

      USE mod_xc
      USE mo_control_bgc
      USE mo_param1_bgc 

      IMPLICIT NONE 
     
      INTEGER :: kpie,kpje,kpke

      INTEGER :: m,n,errstat,iounit,checkdp(nbgcmax)
      LOGICAL :: isopen,exists      

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

!     Determ.NE.number of output groups 
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

!     Re-def.NE.index variables according to namelist 
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
      ENDDO 

      i_atm_m2d=i_bsc_m2d
      DO n=1,nbgc
#if defined(DIFFAT) || defined(CCSMCOUPLED)
        IF (SRF_ATMCO2(n).GT.0) i_atm_m2d=i_atm_m2d+1
        jatmco2(n)=i_atm_m2d*min(1,SRF_ATMCO2(n))
#endif
#ifdef DIFFAT
        IF (SRF_ATMO2(n).GT.0) i_atm_m2d=i_atm_m2d+1
        jatmo2(n)=i_atm_m2d*min(1,SRF_ATMO2(n))
        IF (SRF_ATMN2(n).GT.0) i_atm_m2d=i_atm_m2d+1
        jatmn2(n)=i_atm_m2d*min(1,SRF_ATMN2(n))
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
        IF (LYR_OMEGAC(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jomegac(n)=i_bsc_m3d*min(1,LYR_OMEGAC(n))
        IF (LYR_NOS(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jnos(n)=i_bsc_m3d*min(1,LYR_NOS(n))
        IF (LYR_DP(n).GT.0) i_bsc_m3d=i_bsc_m3d+1
        jdp(n)=i_bsc_m3d*min(1,LYR_DP(n))

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
        IF (LVL_OMEGAC(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlomegac(n)=ilvl_bsc_m3d*min(1,LVL_OMEGAC(n))
        IF (LVL_NOS(n).GT.0) ilvl_bsc_m3d=ilvl_bsc_m3d+1
        jlvlnos(n)=ilvl_bsc_m3d*min(1,LVL_NOS(n))

        IF (i_bsc_m3d.NE.0) checkdp(n)=1
      ENDDO 
      
      i_iso_m3d=i_bsc_m3d
      ilvl_iso_m3d=ilvl_bsc_m3d
      DO n=1,nbgc
        m=i_iso_m3d

        IF (LYR_DIC13(n).GT.0) i_iso_m3d=i_iso_m3d+1
        jdic13(n)=i_iso_m3d*min(1,LYR_DIC13(n))
        IF (LYR_DIC14(n).GT.0) i_iso_m3d=i_iso_m3d+1
        jdic14(n)=i_iso_m3d*min(1,LYR_DIC14(n))

        IF (LVL_DIC13(n).GT.0) ilvl_iso_m3d=ilvl_iso_m3d+1
        jlvldic13(n)=ilvl_iso_m3d*min(1,LVL_DIC13(n))
        IF (LVL_DIC14(n).GT.0) ilvl_iso_m3d=ilvl_iso_m3d+1
        jlvldic14(n)=ilvl_iso_m3d*min(1,LVL_DIC14(n))

        IF (i_iso_m3d.NE.m) checkdp(n)=1
      ENDDO 
      i_iso_m3d=i_iso_m3d-i_bsc_m3d
      ilvl_iso_m3d=ilvl_iso_m3d-ilvl_bsc_m3d

      i_ant_m3d=0
      ilvl_ant_m3d=0

      i_agg_m3d=i_bsc_m3d+i_iso_m3d+i_ant_m3d
      ilvl_agg_m3d=ilvl_bsc_m3d+ilvl_iso_m3d+ilvl_ant_m3d
      DO n=1,nbgc
        m=i_agg_m3d
 
        IF (LYR_NOS(n).GT.0) i_agg_m3d=i_agg_m3d+1
        jnos(n)=i_agg_m3d*min(1,LYR_NOS(n))

        IF (LVL_NOS(n).GT.0) ilvl_agg_m3d=ilvl_agg_m3d+1
        jlvlnos(n)=ilvl_agg_m3d*min(1,LVL_NOS(n))

        IF (i_agg_m3d.NE.m) checkdp(n)=1
      ENDDO
      i_agg_m3d=i_agg_m3d-i_bsc_m3d-i_iso_m3d-i_ant_m3d
      ilvl_agg_m3d=ilvl_agg_m3d-ilvl_bsc_m3d-ilvl_iso_m3d-ilvl_ant_m3d
 
      i_cfc_t3d=0 
      ilvl_cfc_t3d=0 

!     add dp required 
      DO n=1,nbgc
        IF (checkdp(n).NE.0.AND.LYR_DP(n).EQ.0) THEN 
          i_bsc_m3d=i_bsc_m3d+1
          jdp(n)=i_bsc_m3d+i_agg_m3d+i_iso_m3d+i_ant_m3d+i_cfc_t3d
        ENDIF
      ENDDO 
  
      i_bsc_sed=0
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
         
      nbgcm2d=i_bsc_m2d+i_atm_m2d
      nbgcm3d=i_bsc_m3d+i_iso_m3d+i_agg_m3d
      nbgcm3dlvl=ilvl_bsc_m3d+ilvl_iso_m3d+ilvl_agg_m3d
      nbgct_sed = i_bsc_sed

!     Allocate buffers 

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
      IF (pos.EQ.0) RETURN
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
          CALL ncwrtr(vnm,dims,bgcm2d(1-nbdy,1-nbdy,pos),ip,2,          &
     &      sfac,offs,4)
        ENDIF
      ELSEIF (frmt.EQ.8) THEN
        IF (cmpflg.EQ.1) THEN
          CALL nccomp(vnm,dims,bgcm2d(1-nbdy,1-nbdy,pos),ip,sfac,       &
     &      offs,8)
        ELSE
          CALL ncwrtr(vnm,dims,bgcm2d(1-nbdy,1-nbdy,pos),ip,2,          &
     &      sfac,offs,8)
        ENDIF
      ELSE
        STOP 'unknown output format '
      ENDIF
!
! --- Def.NE.attributes
      IF (len(trim(vunits)).NE.0) CALL ncattr('units',vunits)
      IF (len(trim(vlngnm)).NE.0) CALL ncattr('long_name',vlngnm)
      IF (len(trim(vstdnm)).NE.0) CALL ncattr('standard_name',vstdnm)
      CALL ncattr('coordinates','plon plat')
      CALL ncattr('cell_measures','area: parea')
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
      IF (frmt.EQ.0) RETURN
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
      IF (len(trim(vunits)).NE.0) CALL ncattr('units',vunits)
      IF (len(trim(vlngnm)).NE.0) CALL ncattr('long_name',vlngnm)
      IF (len(trim(vstdnm)).NE.0) CALL ncattr('standard_name',vstdnm)
      CALL ncattr('coordinates','plon plat')
      CALL ncattr('cell_measures','area: parea')
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
      IF (frmt.EQ.0) RETURN
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
      IF (len(trim(vunits)).NE.0) CALL ncattr('units',vunits)
      IF (len(trim(vlngnm)).NE.0) CALL ncattr('long_name',vlngnm)
      IF (len(trim(vstdnm)).NE.0) CALL ncattr('standard_name',vstdnm)
      CALL ncattr('coordinates','plon plat')
      CALL ncattr('cell_measures','area: parea')
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
      IF (frmt.EQ.0) RETURN
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
          CALL ncpack(vnm,dims,bgct_sed(1-nbdy,1-nbdy,1,pos),ip,2,      &
     &      sfac,offs)
        ENDIF
      ELSEIF (frmt.EQ.4) THEN
        IF (cmpflg.EQ.1) THEN
          CALL nccomp(vnm,dims,bgct_sed(1-nbdy,1-nbdy,1,pos),ip,sfac,   &
     &      offs,4)
        ELSE
          CALL ncwrtr(vnm,dims,bgct_sed(1-nbdy,1-nbdy,1,pos),ip,2,      &
     &      sfac,offs,4)
        ENDIF
      ELSEIF (frmt.EQ.8) THEN
        IF (cmpflg.EQ.1) THEN
          CALL nccomp(vnm,dims,bgct_sed(1-nbdy,1-nbdy,1,pos),ip,sfac,   &
     &      offs,8)
        ELSE
          CALL ncwrtr(vnm,dims,bgct_sed(1-nbdy,1-nbdy,1,pos),ip,2,      &
     &      sfac,offs,8)
        ENDIF
      ELSE
        STOP 'unknown output format '
      ENDIF
!
! --- Def.NE.attributes
      IF (len(trim(vunits)).NE.0) CALL ncattr('units',vunits)
      IF (len(trim(vlngnm)).NE.0) CALL ncattr('long_name',vlngnm)
      IF (len(trim(vstdnm)).NE.0) CALL ncattr('standard_name',vstdnm)
      CALL ncattr('coordinates','plon plat')
      CALL ncattr('cell_measures','area: parea')
!
      END SUBROUTINE wrtsdm



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
