! Copyright (C) 2003  P. Wetzel 
! Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger
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
                                              ! layer in BLOM

     REAL, PARAMETER :: safediv = 1.0e-25     ! added to the denominator of isotopic ratios (avoid div. by zero)
! advected tracers
      INTEGER, PARAMETER :: i_base=22,                                  &
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
     &                      idet     =11,                               &
     &                      icalc    =12,                               &
     &                      iopal    =13,                               &
     &                      ian2o    =14,                               &
     &                      idms     =15,                               &
     &                      iiron    =16,                               &
     &                      ifdust   =17,                               &
     &                      iprefo2  =18,                               &
     &                      iprefpo4 =19,                               &
     &                      iprefalk =20,                               &
     &                      iprefdic =21,                               &
     &                      idicsat  =22  
#ifdef cisonew
      INTEGER, PARAMETER :: i_iso=12,                                   &
     &                      isco213  = i_base+1,                        &
     &                      isco214  = i_base+2,                        &
     &                      idoc13   = i_base+3,                        &
     &                      idoc14   = i_base+4,                        &
     &                      iphy13   = i_base+5,                        &
     &                      iphy14   = i_base+6,                        &
     &                      izoo13   = i_base+7,                        &
     &                      izoo14   = i_base+8,                        &
     &                      idet13   = i_base+9,                        &
     &                      idet14   = i_base+10,                       &
     &                      icalc13  = i_base+11,                       &
     &                      icalc14  = i_base+12                                        
#else 
      INTEGER, PARAMETER :: i_iso=0
#endif
#ifdef CFC  
      INTEGER, PARAMETER :: i_cfc=3,                                    &
     &                      icfc11   = i_base+i_iso+1,                  &
     &                      icfc12   = i_base+i_iso+2,                  &
     &                      isf6     = i_base+i_iso+3           
#else 
      INTEGER, PARAMETER :: i_cfc=0
#endif
#ifdef AGG
      INTEGER, PARAMETER :: i_agg=2,                                    &
     &                      inos     = i_base+i_iso+i_cfc+1,            &
     &                      iadust   = i_base+i_iso+i_cfc+2
#else 
      INTEGER, PARAMETER :: i_agg=0
#endif
#ifdef natDIC
      INTEGER, PARAMETER :: i_nat_dic=3,                                &
     &                      inatsco212 = i_base+i_iso+i_cfc+i_agg+1,    &
     &                      inatalkali = i_base+i_iso+i_cfc+i_agg+2,    &
     &                      inatcalc   = i_base+i_iso+i_cfc+i_agg+3
#else 
      INTEGER, PARAMETER :: i_nat_dic=0
#endif

! total number of advected tracers
      INTEGER, PARAMETER :: nocetra=i_base+i_iso+i_cfc+i_agg+i_nat_dic


! ATMOSPHERE
      INTEGER, PARAMETER :: i_base_atm=5,                               &
     &                      iatmco2=1,                                  &
     &                      iatmo2 =2,                                  &
     &                      iatmn2 =3,                                  &
     &                      iatmn2o=4,                                  &
     &                      iatmdms=5

#ifdef cisonew
      INTEGER, PARAMETER :: i_iso_atm = 2,                              &
     &                      iatmc13 = i_base_atm+1,                     &
     &                      iatmc14 = i_base_atm+2
#else
      INTEGER, PARAMETER :: i_iso_atm = 0
#endif

#ifdef CFC
      INTEGER, PARAMETER :: i_cfc_atm = 3,                              &
     &                      iatmf11 = i_base_atm+i_iso_atm+1,           &
     &                      iatmf12 = i_base_atm+i_iso_atm+2,           &
     &                      iatmsf6 = i_base_atm+i_iso_atm+3                      
#else
      INTEGER, PARAMETER :: i_cfc_atm = 0
#endif

#ifdef natDIC
      INTEGER, PARAMETER :: i_ndic_atm = 1,                             &
     &                      iatmnco2 = i_base_atm+i_iso_atm+i_cfc_atm+1
#else
      INTEGER, PARAMETER :: i_ndic_atm = 0
#endif

! total number of atmosphere tracers
      INTEGER, PARAMETER :: natm=i_base_atm+i_iso_atm+i_cfc_atm+i_ndic_atm


! sediment
#ifdef cisonew
      INTEGER, PARAMETER :: issso12=1,                                  &
     &                      isssc12=2,                                  &
     &                      issssil=3,                                  &
     &                      issster=4,                                  &
     &                      issso13=5,                                  &
     &                      issso14=6,                                  &
     &                      isssc13=7,                                  &
     &                      isssc14=8,                                  &
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
