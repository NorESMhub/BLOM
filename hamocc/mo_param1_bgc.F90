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
!******************************************************************************
!
! MODULE mo_param1_bgc - bgc tracer parameters.
!
!  Patrick Wetzel    *MPI-Met, HH*    01.09.03
!
!  
!  Modified
!  --------
!  J.Schwinger,        *NORCE Climate, Bergen*    2020-05-26
!
!  - To facilitate easier use of 'only-lists' in use statements, make indices 
!    always defined also in case they are inside a #ifdef directive.
!    
!
!  Purpose
!  -------
!  - definition of indices in tracer arrays
!
!******************************************************************************
      implicit none
      
      INTEGER, PARAMETER :: ks=12,ksp=ks+1    ! ks: nb of sediment layers

      REAL,    PARAMETER :: safediv = 1.0e-25 ! added to the denominator of isotopic ratios (avoid div. by zero)

! Tracer indices
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
      INTEGER, PARAMETER :: i_iso=0,                                    &
     &                      isco213  = -1,                              &
     &                      isco214  = -1,                              &
     &                      idoc13   = -1,                              &
     &                      idoc14   = -1,                              &
     &                      iphy13   = -1,                              &
     &                      iphy14   = -1,                              &
     &                      izoo13   = -1,                              &
     &                      izoo14   = -1,                              &
     &                      idet13   = -1,                              &
     &                      idet14   = -1,                              &
     &                      icalc13  = -1,                              &
     &                      icalc14  = -1
#endif
#ifdef CFC  
      INTEGER, PARAMETER :: i_cfc=3,                                    &
     &                      icfc11   = i_base+i_iso+1,                  &
     &                      icfc12   = i_base+i_iso+2,                  &
     &                      isf6     = i_base+i_iso+3           
#else 
      INTEGER, PARAMETER :: i_cfc=0,                                    &
     &                      icfc11   = -1,                              &
     &                      icfc12   = -1,                              &
     &                      isf6     = -1
#endif
#ifdef AGG
      INTEGER, PARAMETER :: i_agg=2,                                    &
     &                      inos     = i_base+i_iso+i_cfc+1,            &
     &                      iadust   = i_base+i_iso+i_cfc+2
#else 
      INTEGER, PARAMETER :: i_agg=0,                                    &
     &                      inos     = -1,                              &
     &                      iadust   = -1
#endif
#ifdef natDIC
      INTEGER, PARAMETER :: i_nat_dic=3,                                &
     &                      inatsco212 = i_base+i_iso+i_cfc+i_agg+1,    &
     &                      inatalkali = i_base+i_iso+i_cfc+i_agg+2,    &
     &                      inatcalc   = i_base+i_iso+i_cfc+i_agg+3
#else 
      INTEGER, PARAMETER :: i_nat_dic=0,                                &
     &                      inatsco212 = -1,                            &
     &                      inatalkali = -1,                            &
     &                      inatcalc   = -1
#endif
#ifdef BROMO
      INTEGER, PARAMETER :: i_bromo=1,                                  &
     &                      ibromo=i_base+i_iso+i_cfc+i_agg+i_nat_dic+1
#else 
      INTEGER, PARAMETER :: i_bromo=0,                                  &
     &                      ibromo=-1
#endif
#ifdef extNcycle
      INTEGER, PARAMETER :: i_extn=2,                                   &
     &                      iano2 = i_base+i_iso+i_cfc+i_agg+i_nat_dic  &
     &                             +i_bromo+1,                          &
     &                      ianh4 = i_base+i_iso+i_cfc+i_agg+i_nat_dic  &
     &                             +i_bromo+2
#else
      INTEGER, PARAMETER :: i_extn=0,                                   &
     &                      iano2 = -1,                                 &
     &                      ianh4 = -1    
#endif
! total number of advected tracers
      INTEGER, PARAMETER :: nocetra=i_base+i_iso+i_cfc+i_agg+i_nat_dic  &
                                    +i_bromo+i_extn


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
      INTEGER, PARAMETER :: i_iso_atm = 0,                              &
     &                      iatmc13 = -1,                               &
     &                      iatmc14 = -1
#endif

#ifdef CFC
      INTEGER, PARAMETER :: i_cfc_atm = 3,                              &
     &                      iatmf11 = i_base_atm+i_iso_atm+1,           &
     &                      iatmf12 = i_base_atm+i_iso_atm+2,           &
     &                      iatmsf6 = i_base_atm+i_iso_atm+3                      
#else
      INTEGER, PARAMETER :: i_cfc_atm = 0,                              &
     &                      iatmf11 = -1,                               &
     &                      iatmf12 = -1,                               &
     &                      iatmsf6 = -1
#endif

#ifdef natDIC
      INTEGER, PARAMETER :: i_ndic_atm = 1,                             &
     &                      iatmnco2 = i_base_atm+i_iso_atm+i_cfc_atm+1
#else
      INTEGER, PARAMETER :: i_ndic_atm = 0,                             &
     &                      iatmnco2 = -1
#endif
#ifdef BROMO
      INTEGER, PARAMETER :: i_bromo_atm=1,                              &
     &                      iatmbromo=i_base_atm+i_iso_atm+i_cfc_atm+   &
     &                                i_ndic_atm+1
#else 
      INTEGER, PARAMETER :: i_bromo_atm=0,                              &
     &                      iatmbromo=-1
#endif
#ifdef extNcycle
      INTEGER, PARAMETER :: i_nh3_atm=1,                                &
     &                      iatmnh3=i_base_atm+i_iso_atm+i_cfc_atm+     &
     &                                i_ndic_atm+i_bromo_atm+1
#else 
      INTEGER, PARAMETER :: i_nh3_atm=0,                                &
     &                      iatmnh3=-1
#endif

! total number of atmosphere tracers
      INTEGER, PARAMETER :: natm=i_base_atm+i_iso_atm+i_cfc_atm+i_ndic_atm+i_bromo_atm+i_nh3_atm


! rivers
      integer, parameter :: nriv   =7                                     ! size of river input field
      integer, parameter :: irdin  =1,                                  & ! dissolved inorganic nitrogen
     &                      irdip  =2,                                  & ! dissolved inorganic phosphorous 
     &                      irsi   =3,                                  & ! dissolved silicate 
     &                      iralk  =4,                                  & ! alkalinity
     &                      iriron =5,                                  & ! dissolved bioavailable iron
     &                      irdoc  =6,                                  & ! dissolved organic carbon
     &                      irdet  =7                                     ! particulate carbon

      
! ---  sediment
     ! sediment solid components 
      INTEGER, PARAMETER :: i_sed_base = 4
      INTEGER, PARAMETER :: issso12=1,                                  &
     &                      isssc12=2,                                  &
     &                      issssil=3,                                  &
     &                      issster=4
#ifdef cisonew 
      INTEGER, PARAMETER :: i_sed_cisonew = 4
      INTEGER, PARAMETER :: issso13 = i_sed_base+1,                     &
     &                      issso14 = i_sed_base+2,                     &
     &                      isssc13 = i_sed_base+3,                     &
     &                      isssc14 = i_sed_base+4
#else
      INTEGER, PARAMETER :: i_sed_cisonew = 0
      INTEGER, PARAMETER :: issso13 = -1,                               &
     &                      issso14 = -1,                               &
     &                      isssc13 = -1,                               &
     &                      isssc14 = -1 
#endif
      INTEGER, PARAMETER :: nsedtra = i_sed_base + i_sed_cisonew


     ! sediment pore water components
      INTEGER, PARAMETER :: i_pow_base=7
      INTEGER, PARAMETER :: ipowaic=1,                                  &
     &                      ipowaal=2,                                  &
     &                      ipowaph=3,                                  &
     &                      ipowaox=4,                                  &
     &                      ipown2 =5,                                  &
     &                      ipowno3=6,                                  &
     &                      ipowasi=7
#ifdef cisonew
      INTEGER, PARAMETER :: i_pow_cisonew = 2
      INTEGER, PARAMETER :: ipowc13=i_pow_base + 1,                     &
     &                      ipowc14=i_pow_base + 2                       
#else
      INTEGER, PARAMETER :: i_pow_cisonew = 0
      INTEGER, PARAMETER :: ipowc13 = -1,                               &
     &                      ipowc14 = -1  
#endif
#ifdef extNcycle
      INTEGER, PARAMETER :: i_pow_extNcycle = 3 !indices not corresponding to ocetra
      INTEGER, PARAMETER :: ipownh4=i_pow_base + i_pow_cisonew+1,       &  
     &                      ipown2o=i_pow_base + i_pow_cisonew+2,       &
     &                      ipowno2=i_pow_base + i_pow_cisonew+3
#else
      INTEGER, PARAMETER :: i_pow_extNcycle = 0
      INTEGER, PARAMETER :: ipownh4 = -1,                               &
     &                      ipown2o = -1,                               &
     &                      ipowno2 = -1
#endif
      INTEGER, PARAMETER :: npowtra = i_pow_base + i_pow_cisonew + i_pow_extNcycle

     ! Mapping between pore water and ocean tracers needed for pore water diffusion
      INTEGER, SAVE      :: map_por2octra(npowtra)

      contains

      subroutine init_por2octra_mapping()
        
        map_por2octra(ipowaic) = isco212 
        map_por2octra(ipowaal) = ialkali 
        map_por2octra(ipowaph) = iphosph 
        map_por2octra(ipowaox) = ioxygen 
        map_por2octra(ipown2)  = igasnit 
        map_por2octra(ipowno3) = iano3 
        map_por2octra(ipowasi) = isilica 
       
        ! if statements for non-base tracers 
#ifdef cisonew
        map_por2octra(ipowc13) = isco213 
        map_por2octra(ipowc14) = isco214
#endif
#ifdef extNcycle
        map_por2octra(ipownh4) = ianh4
        map_por2octra(ipown2o) = ian2o
        map_por2octra(ipowno2) = iano2
#endif
      end subroutine init_por2octra_mapping

!******************************************************************************
      END MODULE mo_param1_bgc
