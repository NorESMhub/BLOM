      MODULE mo_param1_bgc


!***********************************************************************
!
!**** *MODULE mo_param1_bgc* - bgc tracer parameters.
!
!     Patrick Wetzel    *MPI-Met, HH*    01.09.03
!  
!     Purpose
!     -------
!     - declaration and memory allocation
!
!**********************************************************************
      implicit none
      
      INTEGER, PARAMETER :: ks=12,ksp=ks+1    ! ks: nb of sediment layers

      REAL,    PARAMETER :: dp_ez  = 100.0    ! depth of euphotic zone
      REAL,    PARAMETER :: dp_min = 1.0E-12  ! min layer thickness layers thinner than this are 
                                              ! ignored by HAMOCC
      REAL,    PARAMETER :: dp_min_sink = 1.0 ! min layer thickness for sinking (layers thinner than 
                                              ! this are ignored and set to the concentration of the 
                                              ! layer above). Note that the bottom layer index kbo(i,j)
                                              ! is defined as the lowermost layer thicker than dp_min_sink.

     INTEGER,  PARAMETER :: kmle   = 2        ! k-end index for layers that represent the mixed
                                              ! layer in MICOM

! advected tracers
      INTEGER, PARAMETER :: i_base_adv=17,                              &
     &                      isco212  =1,                                &
     &                      ialkali  =2,                                &
     &                      iphosph  =3,                                &
     &                      ioxygen  =4,                                &
     &                      igasnit  =5,                                &
     &                      iano3    =6,                                &
     &                      isilica  =7,                                &
     &                      idoc     =8,                                &
     &                      iphy     =9,                                &
     &                      izoo     =10,                               &
     &                      ian2o    =11,                               &
     &                      idms     =12,                               &
     &                      iiron    =13,                               &
     &                      ifdust   =14,                               &
     &                      iprefo2  =15,                               &
     &                      iprefpo4 =16,                               &
     &                      iprefalk =17  
#ifdef __c_isotopes
      INTEGER, PARAMETER ::                                             &
     &                      isco213  =i_base_adv+1,                     &
     &                      isco214  =i_base_adv+2,                     &
     &                      i_iso_adv=2                    
#else 
      INTEGER, PARAMETER ::                                             &
     &                      i_iso_adv=0
#endif

      INTEGER, PARAMETER ::                                             &
#ifdef CFC  
     &                      i_cfc_adv= 3,                               &
     &                      icfc11   = i_base_adv+i_iso_adv+1,          &
     &                      icfc12   = i_base_adv+i_iso_adv+2,          &
     &                      isf6     = i_base_adv+i_iso_adv+3           
#else 
     &                      i_cfc_adv= 0
#endif
      INTEGER, PARAMETER ::                                             &
#ifdef AGG
     &                      i_agg_adv= 2,                               &
     &                      inos     = i_base_adv+i_iso_adv+i_cfc_adv+1,&
     &                      iadust   = i_base_adv+i_iso_adv+i_cfc_adv+2
#else 
                            i_agg_adv= 0
#endif

! total number of advected tracers
      INTEGER, PARAMETER :: ntraad=i_base_adv+i_iso_adv+i_cfc_adv+i_agg_adv

! non-advected (fast sinking) tracers
      INTEGER, PARAMETER ::                                             &
     &                      idet     =ntraad+1,                         &
     &                      icalc    =ntraad+2,                         &
     &                      iopal    =ntraad+3,                         &
     &                      i_base   =3 
                            
      INTEGER, PARAMETER ::                                             &
#ifdef __c_isotopes
     &                      idet13   =ntraad+i_base+1,                  &
     &                      icalc13  =ntraad+i_base+2,                  &      
     &                      idet14   =ntraad+i_base+3,                  &
     &                      icalc14  =ntraad+i_base+4,                  &
     &                      i_iso    =4  
#else
     &                      i_iso    =0

#endif

     INTEGER, PARAMETER :: nocetra = ntraad+i_base+i_iso

! ATMOSPHERE
      INTEGER, PARAMETER :: iatmco2=1,                                  &
     &                      iatmo2 =2,                                  &
     &                      iatmn2 =3,                                  &
     &                      iatmn2o=4,                                  &
     &                      iatmdms=5,                                  &
     &                      i_base_atm=5
      INTEGER, PARAMETER ::                                             &
#ifdef __c_isotopes
     &                      iatmc13 = i_base_atm+1,                     &
     &                      iatmc14 = i_base_atm+2,                     &
     &                      i_iso_atm = 2
#else
     &                      i_iso_atm = 0
#endif

      INTEGER, PARAMETER ::  natm=i_base_atm+i_iso_atm

! sediment
#ifdef __c_isotopes
      INTEGER, PARAMETER :: issso12=1,                                  &
     &                      isssc12=2,                                  &
     &                      issssil=3,                                  &
     &                      issster=4,                                  &
     &                      issso13=5,                                  &
     &                      issso14=6,                                  &
     &                      isssc13=7,                                  &
     &                      isssc14=8,                                   &
     &                      nsedtra=8
     
! pore water tracers, index should be the same as for ocetra
      INTEGER, PARAMETER :: ipowaic=1,npowtra=9,                        &
     &                      ipowaal=2,                                  &
     &                      ipowaph=3,                                  &
     &                      ipowaox=4,                                  &
     &                      ipown2 =5,                                  &
     &                      ipowno3=6,                                  &
     &                      ipowasi=7,                                  &
     &                      ipowc13=8,                                  &
     &                      ipowc14=9
#else
      INTEGER, PARAMETER :: issso12=1,                                  &
     &                      isssc12=2,                                  &
     &                      issssil=3,                                  &
     &                      issster=4,                                  &
     &                      nsedtra=4

! pore water tracers, index should be the same as for ocetra
     INTEGER, PARAMETER :: ipowaic=1,npowtra=7,                        &
     &                      ipowaal=2,                                  &
     &                      ipowaph=3,                                  &
     &                      ipowaox=4,                                  &
     &                      ipown2 =5,                                  &
     &                      ipowno3=6,                                  &
     &                      ipowasi=7
#endif

      END MODULE mo_param1_bgc
