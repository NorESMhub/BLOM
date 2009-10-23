      MODULE mo_bgcmean


!***********************************************************************
!
!**** *MODULE mo_bgcmean* - Variables for bgcmean.
!
!     Patrick Wetzel    *MPI-Met, HH*    09.12.02
!  
!     Purpose
!     -------
!     - declaration and memory allocation
!
!**********************************************************************
      implicit none

      INTEGER, PARAMETER ::                                   &
     &          jco2flux  =1,                                 &
     &          jco214f   =2,                                 &         
     &          jo2flux   =3,                                 &
     &          jn2flux   =4,                                 &
     &          jn2oflux  =5,                                 &
     &          jprorca   =6,                                 &
     &          jprcaca   =7,                                 &
     &          jsilpro   =8,                                 &
     &          jprodus   =9,                                 &
     &          nbgct2d   =9
      
!----------------------------------------------------------------      
      
      INTEGER, PARAMETER ::                                   &
     &          i_bsc_m2d=14,                                 &
     &          jkwco2   =1,                                  &
     &          jpco2    =2,                                  &
     &          jdmsflux =3,                                  &
     &          jco2fxd  =4,                                  &
     &          jco2fxu  =5,                                  &
     &          joxflux  =6,                                  &
     &          jniflux  =7,                                  &
     &          jdms     =8,                                  &
     &          jdmsprod =9,                                  &
     &          jdms_bac =10,                                 &
     &          jdms_uv  =11,                                 &
     &          jexport  =12,                                 &
     &          jexpoca  =13,                                 &
     &          jexposi  =14

#ifdef DIFFAT
      INTEGER, PARAMETER ::                                   &
     &          i_atm_m2d=3,                                  &
     &          jatmco2  =i_bsc_m2d+1,                        &
     &          jatmo2   =i_bsc_m2d+2,                        &
     &          jatmn2   =i_bsc_m2d+3
#elif CCSMCOUPLED
      INTEGER, PARAMETER ::                                   &
     &          i_atm_m2d=1,                                  &
     &          jatmco2  =i_bsc_m2d+1
#else
      INTEGER, PARAMETER :: i_atm_m2d = 0
#endif

      INTEGER, PARAMETER :: i_ant_m2d = 0

#ifdef PCFC
      INTEGER, PARAMETER ::                                   &
     &          i_cfc_m2d=5,                                  &     
     &          jac14fx  =i_bsc_m2d+i_atm_m2d+i_ant_m2d+1,    & 
     &          jcfc11fx =i_bsc_m2d+i_atm_m2d+i_ant_m2d+2,    & 
     &          jcfc12fx =i_bsc_m2d+i_atm_m2d+i_ant_m2d+3,    & 
     &          jpcfc11  =i_bsc_m2d+i_atm_m2d+i_ant_m2d+4,    & 
     &          jpcfc12  =i_bsc_m2d+i_atm_m2d+i_ant_m2d+5
#else
      INTEGER, PARAMETER :: i_cfc_m2d = 0
#endif         

#ifdef PCOMPONENT_ANALYSIS
      INTEGER, PARAMETER ::                                                 &
     &          i_com_m2d=6,                                                &
     &          jdpco2_dalk =i_bsc_m2d+i_atm_m2d+i_ant_m2d+i_cfc_m2d+1,     &
     &          jdpco2_ddic =i_bsc_m2d+i_atm_m2d+i_ant_m2d+i_cfc_m2d+2,     & 
     &          jdpco2_dsst =i_bsc_m2d+i_atm_m2d+i_ant_m2d+i_cfc_m2d+3,     & 
     &          jdpco2_dsss =i_bsc_m2d+i_atm_m2d+i_ant_m2d+i_cfc_m2d+4,     &  
     &          jdpco2_dalk_ant =i_bsc_m2d+i_atm_m2d+i_ant_m2d+i_cfc_m2d+5, & 
     &          jdpco2_ddic_ant =i_bsc_m2d+i_atm_m2d+i_ant_m2d+i_cfc_m2d+6
#else
      INTEGER, PARAMETER :: i_com_m2d = 0
#endif       

      INTEGER, PARAMETER :: nbgcm2d=i_bsc_m2d+i_atm_m2d+i_ant_m2d+i_cfc_m2d+i_com_m2d

!----------------------------------------------------------------
! 2d: jdms,jexport, jexpoca, jexposi
! out : jpoc, jatten, jcalc, jopal

      INTEGER, PARAMETER ::                                   &
     &          i_bsc_m3d=17,                                 &
     &          jphyto  =1,                                   &
     &          jgrazer =2,                                   &  
     &          jdoc    =3,                                   & 
     &          jphosy  =4,                                   &  
     &          jphosph =5,                                   &  
     &          joxygen =6,                                   &  
     &          jiron   =7,                                   &
     &          jano3   =8,                                   &
     &          jalkali =9,                                   & 
     &          jsilica =10,                                  &   
     &          jdic    =11,                                  &  
     &          jpoc    =12,                                  & 
     &          jcalc   =13,                                  & 
     &          jopal   =14,                                  &
     &          jco3    =15,                                  &    
     &          jph     =16,                                  &                                        
     &          jomegac =17         
#ifdef __c_isotopes
      INTEGER, PARAMETER ::                                   &
     &          i_iso_m3d=2,                                  &         
     &          jdic13  =i_bsc_m3d+1,                         &
     &          jdic14  =i_bsc_m3d+2 
#else
      INTEGER, PARAMETER :: i_iso_m3d = 0
#endif

      INTEGER, PARAMETER :: i_ant_m3d = 0

#ifdef AGG
      INTEGER, PARAMETER ::                                   &
     &          i_agg_m3d=1,                                  &
     &          jnos     =i_bsc_m3d+i_ant_m3d+1
#else 
      INTEGER, PARAMETER :: i_agg_m3d = 0
#endif 

#ifdef PCFC
      INTEGER, PARAMETER ::                                   &
     &          i_cfc_t3d=3,                                  &         
     &          jcfc11_t=i_bsc_m3d+i_ant_m3d+1,               & 
     &          jcfc12_t=i_bsc_m3d+i_ant_m3d+2,               & 
     &          jac14_t =i_bsc_m3d+i_ant_m3d+3
#else
      INTEGER, PARAMETER :: i_cfc_t3d = 0
#endif

      INTEGER, PARAMETER :: nbgcm3d =i_bsc_m3d+i_iso_m3d+i_agg_m3d

!----------------------------------------------------------------
! js: sediment
!
      INTEGER, PARAMETER ::                                   &
     &          i_bsc_sed = 11,                               &
     &          jpowaic   =  1,                               &
     &          jpowaal   =  2,                               &
     &          jpowaph   =  3,                               &
     &          jpowaox   =  4,                               &
     &          jpown2    =  5,                               &
     &          jpowno3   =  6,                               &
     &          jpowasi   =  7,                               &
     &          jssso12   =  8,                               &
     &          jssssil   =  9,                               &
     &          jsssc12   = 10,                               &
     &          jssster   = 11          


      INTEGER, PARAMETER :: nbgct_sed = i_bsc_sed

!----------------------------------------------------------------

!      REAL, DIMENSION (:),       ALLOCATABLE :: stepspm      
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: bgct2d
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: bgcm2d
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: bgcm3d
     REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: bgct_sed 
!      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: bgct3d
      
      INTEGER :: meancnt_bgc_2D,meancnt_bgc_3D
      INTEGER :: mean_2D_day,mean_3D_month
      INTEGER :: meantime_2d, nmeantime_2d
      INTEGER :: meantime_3d, nmeantime_3d

      INTEGER :: n90depth, n1000depth, n2000depth
      
      INTEGER :: nc_2d_id, nc_bioz_id, nc_3d_id
 
      INTEGER :: nacc_bgc2d
      
      
      CONTAINS

      SUBROUTINE ALLOC_MEM_BGCMEAN(kpie,kpje,kpke)

      use mod_xc
      use mo_control_bgc
      use mo_param1_bgc 
      
      INTEGER :: kpie,kpje,kpke
      INTEGER :: errstat
      
!        IF (mnproc.eq.1) THEN
!        WRITE(io_stdo_bgc,*)'Memory allocation for variable stepspm ...'
!        WRITE(io_stdo_bgc,*)'First dimension    : ',nmeantime_2d
!        ENDIF

!        ALLOCATE (stepspm(nmeantime_2d))

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgct2d ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',nbgct2d
        ENDIF

        ALLOCATE (bgct2d(kpie,kpje,nbgct2d),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory bgct2d'

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgcm2d ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',nbgcm2d
        ENDIF

        ALLOCATE (bgcm2d(kpie,kpje,nbgcm2d),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory bgcm2d'

!        IF (mnproc.eq.1) THEN
!        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgcm3d ...'
!        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
!        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
!        WRITE(io_stdo_bgc,*)'Third dimension    : ',kwrbioz
!        WRITE(io_stdo_bgc,*)'Forth dimension    : ',nbgcm3d
!        ENDIF
!	
!        ALLOCATE (bgcm3d(kpie,kpje,kwrbioz,nbgcm3d))		

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgcm3d ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        WRITE(io_stdo_bgc,*)'Forth dimension    : ',nbgcm3d
        ENDIF

        ALLOCATE (bgcm3d(kpie,kpje,kpke,nbgcm3d),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory bgcm3d'

        IF (mnproc.eq.1) THEN
        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgctsed ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',ks
        WRITE(io_stdo_bgc,*)'Forth dimension    : ',nbgct_sed
        ENDIF

        ALLOCATE (bgct_sed(kpie,kpje,ks,nbgct_sed),stat=errstat)
        if(errstat.ne.0) stop 'not enough memory bgct_sed'

      END SUBROUTINE ALLOC_MEM_BGCMEAN

      END MODULE mo_bgcmean
