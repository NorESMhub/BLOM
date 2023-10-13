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
     !use mo_control_bgc, only: use_cisonew
      use mo_control_bgc

      implicit none
      public

      INTEGER, PARAMETER :: ks=12,ksp=ks+1    ! ks: nb of sediment layers
      REAL,    PARAMETER :: safediv = 1.0e-25 ! added to the denominator of isotopic ratios (avoid div. by zero)

      ! Tracer indices
      integer :: i_base
      integer :: isco212
      integer :: ialkali
      integer :: iphosph
      integer :: ioxygen
      integer :: igasnit
      integer :: iano3
      integer :: isilica
      integer :: idoc
      integer :: iphy
      integer :: izoo
      integer :: idet
      integer :: icalc
      integer :: iopal
      integer :: ian2o
      integer :: idms
      integer :: iiron
      integer :: ifdust
      integer :: iprefo2
      integer :: iprefpo4
      integer :: iprefalk
      integer :: iprefdic
      integer :: idicsat

      integer :: i_iso
      integer :: isco213
      integer :: isco214
      integer :: idoc13
      integer :: idoc14
      integer :: iphy13
      integer :: iphy14
      integer :: izoo13
      integer :: izoo14
      integer :: idet13
      integer :: idet14
      integer :: icalc13
      integer :: icalc14

      integer :: i_cfc
      integer :: icfc11
      integer :: icfc12
      integer :: isf6

      integer :: i_agg
      integer :: inos
      integer :: iadust

      integer :: i_nat_dic
      integer :: inatsco212
      integer :: inatalkali
      integer :: inatcalc

      integer :: i_bromo
      integer :: ibromo

      ! total number of advected tracers
      integer :: nocetra

      ! atmosphere
      integer :: i_base_atm
      integer :: iatmco2
      integer :: iatmo2
      integer :: iatmn2
      integer :: iatmn2o
      integer :: iatmdms

      integer :: i_iso_atm
      integer :: iatmc13
      integer :: iatmc14
      integer :: i_cfc_atm
      integer :: iatmf11
      integer :: iatmf12
      integer :: iatmsf6
      integer :: i_ndic_atm
      integer :: iatmnco2
      integer :: i_bromo_atm
      integer :: iatmbromo
      integer :: natm ! total number of atmosphere tracers

      ! rivers
      integer :: nriv   ! size of river input field
      integer :: irdin  ! dissolved inorganic nitrogen
      integer :: irdip  ! dissolved inorganic phosphorous
      integer :: irsi   ! dissolved silicate
      integer :: iralk  ! alkalinity
      integer :: iriron ! dissolved bioavailable iron
      integer :: irdoc  ! dissolved organic carbon
      integer :: irdet  ! particulate carbon

      ! ---  sediment
      ! sediment solid components
      integer :: i_sed_base
      integer :: issso12
      integer :: isssc12
      integer :: issssil
      integer :: issster

      integer :: i_sed_cisonew
      integer :: issso13
      integer :: issso14
      integer :: isssc13
      integer :: isssc14
      integer :: nsedtra

      ! sediment pore water components
      integer :: i_pow_base
      integer :: ipowaic
      integer :: ipowaal
      integer :: ipowaph
      integer :: ipowaox
      integer :: ipown2
      integer :: ipowno3
      integer :: ipowasi
      integer :: i_pow_cisonew
      integer :: ipowc13
      integer :: ipowc14
      integer :: npowtra

      ! Mapping between pore water and ocean tracers needed for pore water diffusion
      integer, allocatable :: map_por2octra(:)

    contains

      subroutine init_por2octra_mapping()
        map_por2octra(ipowaic) = isco212
        map_por2octra(ipowaal) = ialkali
        map_por2octra(ipowaph) = iphosph
        map_por2octra(ipowaox) = ioxygen
        map_por2octra(ipown2)  = igasnit
        map_por2octra(ipowno3) = iano3
        map_por2octra(ipowasi) = isilica
        if (use_cisonew) then
           map_por2octra(ipowc13) = isco213
           map_por2octra(ipowc14) = isco214
        end if
      end subroutine init_por2octra_mapping

      subroutine init_indices()

        ! Tracer indices
        i_base   = 22
        isco212  = 1
        ialkali  = 2
        iphosph  = 3
        ioxygen  = 4
        igasnit  = 5
        iano3    = 6
        isilica  = 7
        idoc     = 8
        iphy     = 9
        izoo     = 10
        idet     = 11
        icalc    = 12
        iopal    = 13
        ian2o    = 14
        idms     = 15
        iiron    = 16
        ifdust   = 17
        iprefo2  = 18
        iprefpo4 = 19
        iprefalk = 20
        iprefdic = 21
        idicsat  = 22
        if (use_cisonew) then
           i_iso    = 12
           isco213  = i_base+1
           isco214  = i_base+2
           idoc13   = i_base+3
           idoc14   = i_base+4
           iphy13   = i_base+5
           iphy14   = i_base+6
           izoo13   = i_base+7
           izoo14   = i_base+8
           idet13   = i_base+9
           idet14   = i_base+10
           icalc13  = i_base+11
           icalc14  = i_base+12
        else
           i_iso    =  0
           isco213  = -1
           isco214  = -1
           idoc13   = -1
           idoc14   = -1
           iphy13   = -1
           iphy14   = -1
           izoo13   = -1
           izoo14   = -1
           idet13   = -1
           idet14   = -1
           icalc13  = -1
           icalc14  = -1
        end if
        if (use_CFC) then
           i_cfc=3
           icfc11   = i_base+i_iso+1
           icfc12   = i_base+i_iso+2
           isf6     = i_base+i_iso+3
        else
           i_cfc=0
           icfc11   = -1
           icfc12   = -1
           isf6     = -1
        end if
        if (use_AGG) then
           i_agg=2
           inos     = i_base+i_iso+i_cfc+1
           iadust   = i_base+i_iso+i_cfc+2
        else
           i_agg=0
           inos     = -1
           iadust   = -1
        end if
        if (use_natDIC) then
           i_nat_dic=3
           inatsco212 = i_base+i_iso+i_cfc+i_agg+1
           inatalkali = i_base+i_iso+i_cfc+i_agg+2
           inatcalc   = i_base+i_iso+i_cfc+i_agg+3
        else
           i_nat_dic=0
           inatsco212 = -1
           inatalkali = -1
           inatcalc   = -1
        end if
        if (use_BROMO) then
           i_bromo=1
           ibromo=i_base+i_iso+i_cfc+i_agg+i_nat_dic+1
        else
           i_bromo=0
           ibromo=-1
        end if

        ! total number of advected tracers
        nocetra=i_base+i_iso+i_cfc+i_agg+i_nat_dic +i_bromo

        ! ATMOSPHERE
        i_base_atm=5
        iatmco2=1
        iatmo2 =2
        iatmn2 =3
        iatmn2o=4
        iatmdms=5
        if (use_cisonew) then
           i_iso_atm = 2
           iatmc13 = i_base_atm+1
           iatmc14 = i_base_atm+2
        else
           i_iso_atm = 0
           iatmc13 = -1
           iatmc14 = -1
        end if
        if (use_CFC) then
           i_cfc_atm = 3
           iatmf11 = i_base_atm+i_iso_atm+1
           iatmf12 = i_base_atm+i_iso_atm+2
           iatmsf6 = i_base_atm+i_iso_atm+3
        else
           i_cfc_atm = 0
           iatmf11 = -1
           iatmf12 = -1
           iatmsf6 = -1
        end if
        if (use_natDIC) then
           i_ndic_atm = 1
           iatmnco2 = i_base_atm+i_iso_atm+i_cfc_atm+1
        else
           i_ndic_atm = 0
           iatmnco2 = -1
        end if
        if (use_BROMO) then
           i_bromo_atm=1
           iatmbromo=i_base_atm+i_iso_atm+i_cfc_atm+ i_ndic_atm+1
        else
           i_bromo_atm=0
           iatmbromo=-1
        end if

        ! total number of atmosphere tracers
        natm=i_base_atm+i_iso_atm+i_cfc_atm+i_ndic_atm+i_bromo_atm

        ! rivers
        nriv   =7
        irdin  =1
        irdip  =2
        irsi   =3
        iralk  =4
        iriron =5
        irdoc  =6
        irdet  =7

        ! ---  sediment
        ! sediment solid components
        i_sed_base = 4
        issso12 = 1
        isssc12 = 2
        issssil = 3
        issster = 4
        if (use_cisonew) then
           i_sed_cisonew = 4
           issso13 = i_sed_base+1
           issso14 = i_sed_base+2
           isssc13 = i_sed_base+3
           isssc14 = i_sed_base+4
        else
           i_sed_cisonew = 0
           issso13 = -1
           issso14 = -1
           isssc13 = -1
           isssc14 = -1
        end if
        nsedtra = i_sed_base + i_sed_cisonew

        ! sediment pore water components
        i_pow_base =7
        ipowaic    =1
        ipowaal    =2
        ipowaph    =3
        ipowaox    =4
        ipown2     =5
        ipowno3    =6
        ipowasi    =7
        if (use_cisonew) then
           i_pow_cisonew = 2
           ipowc13=i_pow_base + 1
           ipowc14=i_pow_base + 2
        else
           i_pow_cisonew = 0
           ipowc13 = -1
           ipowc14 = -1
        end if
        npowtra = i_pow_base + i_pow_cisonew

        allocate(map_por2octra(-1:npowtra))

      end subroutine init_indices

!******************************************************************************
    END MODULE mo_param1_bgc
