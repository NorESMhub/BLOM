      MODULE mo_timeser_bgc

!
!$Source: /scratch/local1/m212047/patrick/SRC_MPI/src_hamocc/RCS/mo_timeser_bgc.f90,v $\\
!$Revision: 1.1 $\\
!$Date: 2005/01/28 08:37:45 $\\
!$Name:  $\\
!
!***********************************************************************
!
!**** *MODULE mo_timeser_bgc* - Parameter and memory for time series.
!
!     S.Legutke,        *MPI-MaD, HH*    04.08.01
!
!     Modified
!     --------
!     
!     Purpose
!     -------
!     - declaration
!
!     *nfreqts1*     *INTEGER*  - sampling frequency of time series 1.
!     *nts1*         *INTEGER*  - no. of elements in time series 1.
!     *nvarts1*      *INTEGER*  - no. of variables sampled for each element.
!     *its1*         *INTEGER*  - 1st index of grid cell for samples in time series 1.
!     *jts1*         *INTEGER*  - 2nd index of grid cell for samples in time series 1.
!     *lts1*         *INTEGER*  - actual sample counter of time series 1.
!
!**********************************************************************
      implicit none

      INTEGER, PARAMETER :: nts=8

      REAL, DIMENSION (:,:,:), ALLOCATABLE :: ts1

      INTEGER :: lents1,nelets1

      INTEGER :: io_timeser_bgc = 27  ! logical unit number for time series output

      INTEGER :: its1(nts),jts1(nts),nts1,nfreqts1,lts1
      INTEGER :: k1ts1(nts),k2ts1(nts),k3ts1(nts)

      REAL :: rlonts1(nts),rlatts1(nts)
      REAL :: rdep1ts1(nts),rdep2ts1(nts),rdep3ts1(nts)


      INTEGER, PARAMETER ::                                    &  
     &          i_bgc_ts   =26,                                &
     &          itssco212  =1,                                 &
     &          itsphosy   =2,                                 &
     &          itsphosph  =3,                                 &
     &          itsoxygen  =4,                                 &
     &          itsgasnit  =5,                                 &
     &          itsano3    =6,                                 &
     &          itssilica  =7,                                 &
     &          itsdoc     =8,                                 &
     &          itsphy     =9,                                 &
     &          itszoo     =10,                                &
     &          itsdet     =11,                                &
     &          itscalc    =12,                                &
     &          itsopal    =13,                                &
     &          itsiron    =14,                                &
     &          its1fdet   =15,                                &
     &          its1fopa   =16,                                &
     &          its1fcal   =17,                                &
     &          its2fdet   =18,                                &
     &          its2fopa   =19,                                &
     &          its2fcal   =20,                                &
     &          its3fdet   =21,                                &
     &          its3fopa   =22,                                &
     &          its3fcal   =23,                                &
     &          itspco2    =24,                                &
     &          itsatm     =25,                                &
     &          itsco2f    =26

#ifdef AGG
!      INTEGER, PARAMETER ::                                    &  
!     &          i_agg_ts   =1,                                 &        
!     &          itsnos     =i_bgc_ts+i_ant_ts+1
      INTEGER, PARAMETER ::                                    &  
     &          i_agg_ts   =1
#else  
      INTEGER, PARAMETER :: i_agg_ts = 0 
#endif  

      INTEGER, PARAMETER :: nvarts1 =i_bgc_ts+i_agg_ts



      END MODULE mo_timeser_bgc
