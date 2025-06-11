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

module mo_param1_bgc

  !*************************************************************************************************
  ! Definition of indices in bgc tracer arrays
  !
  !  Patrick Wetzel    *MPI-Met, HH*    01.09.03
  !
  !  Modified
  !  J.Schwinger,        *NORCE Climate, Bergen*    2020-05-26
  !  T. Bourgeois,     *NORCE climate, Bergen*   2025-04-14
  !  - implement R2OMIP protocol
  !*************************************************************************************************

  use mo_control_bgc, only: use_BROMO, use_AGG, use_WLIN, use_natDIC, use_CFC,                     &
                            use_cisonew, use_PBGC_OCNP_TIMESTEP, use_PBGC_CK_TIMESTEP,             &
                            use_FB_BGC_OCE, use_BOXATM, use_sedbypass, use_extNcycle,              &
                            use_pref_tracers,use_sediment_quality
  implicit none
  public

  integer, parameter :: ks=12,ksp=ks+1    ! ks: nb of sediment layers
  real,    parameter :: safediv = 1.0e-25 ! added to the denominator of isotopic ratios (avoid div. by zero)

  ! ------------------
  ! Tracer indices
  ! ------------------

  integer :: i_base
  integer, protected :: isco212
  integer, protected :: ialkali
  integer, protected :: iphosph
  integer, protected :: ioxygen
  integer, protected :: igasnit
  integer, protected :: iano3
  integer, protected :: isilica
  integer, protected :: idoc
  integer, protected :: iphy
  integer, protected :: izoo
  integer, protected :: idet
  integer, protected :: icalc
  integer, protected :: iopal
  integer, protected :: ian2o
  integer, protected :: idms
  integer, protected :: iiron
  integer, protected :: ifdust
  integer, protected :: idicsat

  ! Indices for preformed tracers
  integer, protected :: i_pref
  integer, protected :: iprefo2
  integer, protected :: iprefpo4
  integer, protected :: iprefalk
  integer, protected :: iprefdic
  integer, protected :: iprefsilica

  ! Indices for C-isotope tracers
  integer, protected :: i_iso
  integer, protected :: isco213
  integer, protected :: isco214
  integer, protected :: idoc13
  integer, protected :: idoc14
  integer, protected :: iphy13
  integer, protected :: iphy14
  integer, protected :: izoo13
  integer, protected :: izoo14
  integer, protected :: idet13
  integer, protected :: idet14
  integer, protected :: icalc13
  integer, protected :: icalc14

  ! Indices for CFCs
  integer, protected :: i_cfc
  integer, protected :: icfc11
  integer, protected :: icfc12
  integer, protected :: isf6

  ! Indices for tracers related to aggregation scheme
  integer, protected :: i_agg
  integer, protected :: inos
  integer, protected :: iadust

  ! Indices for tracers related to natural DIC
  integer, protected :: i_nat_dic
  integer, protected :: inatsco212
  integer, protected :: inatalkali
  integer, protected :: inatcalc

  ! Indices for bromoform tracer
  integer, protected :: i_bromo
  integer, protected :: ibromo

  ! Indices for DOM tracer
  integer, protected :: i_dom
  integer, protected :: idocsl
  integer, protected :: idocr
  integer, protected :: idocsr
  integer, protected :: i_prefdom
  integer, protected :: iprefdoc
  integer, protected :: iprefdocsl
  integer, protected :: iprefdocsr
  integer, protected :: iprefdocr

  ! Indices for extended nitrogen cycle
  integer, protected :: i_extn
  integer, protected :: ianh4
  integer, protected :: iano2

  ! Indices for the age tracer for shelf water residence time
  integer, protected :: i_shelfage
  integer, protected :: ishelfage

  ! Indices for the R2O riverine terrestrial DOC
  integer, protected :: i_r2o
  integer, protected :: itdoc_lc
  integer, protected :: itdoc_hc
  integer, protected :: itdoc_lc13
  integer, protected :: itdoc_hc13
  integer, protected :: itdoc_lc14
  integer, protected :: itdoc_hc14

  ! total number of advected tracers (set by allocate_tracers in mod_tracers.F90)
  integer :: nocetra

  ! ------------------
  ! atmosphere
  ! ------------------

  integer, protected :: i_base_atm
  integer, protected :: iatmco2
  integer, protected :: iatmo2
  integer, protected :: iatmn2
  integer, protected :: iatmn2o
  integer, protected :: iatmdms

  ! Indices of C-isotope tracers in atm
  integer, protected :: i_iso_atm
  integer, protected :: iatmc13
  integer, protected :: iatmc14

  ! Indices of CFCs in atm
  integer, protected :: i_cfc_atm
  integer, protected :: iatmf11
  integer, protected :: iatmf12
  integer, protected :: iatmsf6

  ! Indices for tracers related to natDIC scheme in atm
  integer, protected :: i_ndic_atm
  integer, protected :: iatmnco2

  ! Indices for bromoform tracer in atm
  integer, protected :: i_bromo_atm
  integer, protected :: iatmbromo

  ! Indices for extended nitrogen tracer in atm
  integer, protected :: i_nh3_atm
  integer, protected :: iatmnh3

  integer, protected :: natm ! total number of atmosphere tracers

  ! --------------------
  ! Nitrogen deposition
  ! --------------------
  integer, protected :: nndep   ! size of N-deposition input field
  integer, protected :: idepnoy ! index for NOy deposition
  integer, protected :: idepnhx ! index for NHx deposition

  ! --------------------
  ! Dust/fe deposition
  ! --------------------
  integer, protected :: ndust   ! size of dust deposition input field
  integer, protected :: itdust  ! total atmospheric dust deposition
  integer, protected :: isfe    ! atmospheric deposition of soluble Fe

  ! ------------------
  ! rivers
  ! ------------------
  integer, protected :: nriv   ! size of river input field
  integer, protected :: irdin  ! dissolved inorganic nitrogen
  integer, protected :: irdip  ! dissolved inorganic phosphorous
  integer, protected :: irsi   ! dissolved silicate
  integer, protected :: iralk  ! alkalinity
  integer, protected :: iriron ! dissolved bioavailable iron
  integer, protected :: irdoc  ! dissolved organic carbon
  integer, protected :: irtdoc ! terrestrial semi-labile DOC
  integer, protected :: irdet  ! particulate carbon

  ! ------------------
  ! sediment
  ! ------------------
  ! Indices for solid sediment components
  integer, protected :: i_sed_base
  integer, protected :: issso12
  integer, protected :: isssc12
  integer, protected :: issssil
  integer, protected :: issster

  ! Indices for C-isotope tracers in sediment
  integer, protected :: i_sed_cisonew
  integer, protected :: issso13
  integer, protected :: issso14
  integer, protected :: isssc13
  integer, protected :: isssc14

  ! Indice for POC age
  integer, protected :: i_sed_age
  integer, protected :: issso12_age
  integer, protected :: nsedtra ! total number of solid sediment tracers
  integer, protected :: nsedtra_woage ! total number of solid sediment tracers without age tracer

  ! Indices for tracers in sediment pore water
  integer, protected :: i_pow_base
  integer, protected :: ipowaic
  integer, protected :: ipowaal
  integer, protected :: ipowaph
  integer, protected :: ipowaox
  integer, protected :: ipown2
  integer, protected :: ipowno3
  integer, protected :: ipowasi

  ! Indices for C-isotope tracers in sediment pore water
  integer, protected :: i_pow_cisonew
  integer, protected :: ipowc13
  integer, protected :: ipowc14
  integer, protected :: npowtra ! computed in init_indices

  ! Indices for extended nitrogen cycle
  integer, protected :: i_pow_extNcycle
  integer, protected :: ipownh4
  integer, protected :: ipown2o
  integer, protected :: ipowno2

  ! Mapping between pore water and ocean tracers needed for pore
  ! water diffusion
  integer, protected, allocatable :: map_por2octra(:)

contains

  ! ===========================================================================
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
    endif
    if (use_extNcycle) then
        map_por2octra(ipownh4) = ianh4
        map_por2octra(ipown2o) = ian2o
        map_por2octra(ipowno2) = iano2
    endif
  end subroutine init_por2octra_mapping

  ! ===========================================================================
  subroutine init_indices()

    use mod_xc        , only: lp, mnproc
    use mo_control_bgc, only: bgc_namelist,get_bgc_namelist, io_stdo_bgc
    use mo_control_bgc, only: use_BROMO,use_AGG,use_WLIN,use_natDIC,use_CFC,use_cisonew,           &
                              use_sedbypass,use_PBGC_OCNP_TIMESTEP,use_PBGC_CK_TIMESTEP,           &
                              use_FB_BGC_OCE, use_BOXATM,use_extNcycle,use_pref_tracers,           &
                              use_coupler_ndep,use_shelfsea_res_time,use_river2omip,use_DOMclasses

    integer :: iounit

    namelist / config_bgc / use_BROMO,use_AGG,use_WLIN,use_natDIC,use_CFC,use_cisonew,             &
                            use_sedbypass,use_PBGC_OCNP_TIMESTEP,use_PBGC_CK_TIMESTEP,             &
                            use_FB_BGC_OCE,use_BOXATM,use_extNcycle,use_pref_tracers,              &
                            use_coupler_ndep,use_shelfsea_res_time,use_sediment_quality,           &
                            use_river2omip,use_DOMclasses

    io_stdo_bgc = lp              !  standard out.

    if(.not. allocated(bgc_namelist)) call get_bgc_namelist()
    open (newunit=iounit, file=bgc_namelist, status='old', action='read')
    read (unit=iounit, nml=config_bgc)
    close (unit=iounit)

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)
      write(io_stdo_bgc,*) 'iHAMOCC: reading namelist CONFIG_BGC'
      write(io_stdo_bgc,nml=config_bgc)
    endif

    ! Tracer indices
    i_base   = 18
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
    idicsat  = 18
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
    endif
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
    endif
    if (use_AGG) then
      i_agg=2
      inos     = i_base+i_iso+i_cfc+1
      iadust   = i_base+i_iso+i_cfc+2
    else
      i_agg=0
      inos     = -1
      iadust   = -1
    endif
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
    endif
    if (use_BROMO) then
      i_bromo=1
      ibromo=i_base+i_iso+i_cfc+i_agg+i_nat_dic+1
    else
      i_bromo=0
      ibromo=-1
    endif
    if (use_extNcycle) then
      i_extn = 2
      iano2  = i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+1
      ianh4  = i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+2
    else
      i_extn = 0
      iano2  = -1
      ianh4  = -1
    endif
    if (use_pref_tracers) then
      i_pref      = 5
      iprefo2     = i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+i_extn+1
      iprefpo4    = i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+i_extn+2
      iprefalk    = i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+i_extn+3
      iprefdic    = i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+i_extn+4
      iprefsilica = i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+i_extn+5
    else
      i_pref      = 0
      iprefo2     = -1
      iprefpo4    = -1
      iprefalk    = -1
      iprefdic    = -1
      iprefsilica = -1
    endif
    if (use_shelfsea_res_time) then
      i_shelfage = 1
      ishelfage  = i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+i_extn+i_pref+1
    else
      i_shelfage = 0
      ishelfage  = -1
    endif
    if (use_river2omip) then
      itdoc_lc  = i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+i_extn+i_pref+i_shelfage+1
      itdoc_hc  = i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+i_extn+i_pref+i_shelfage+2
      if (use_cisonew) then
        i_r2o = 6
        itdoc_lc13  = i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+i_extn+i_pref+i_shelfage+3
        itdoc_hc13  = i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+i_extn+i_pref+i_shelfage+4
        itdoc_lc14  = i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+i_extn+i_pref+i_shelfage+5
        itdoc_hc14  = i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+i_extn+i_pref+i_shelfage+6
      else
        i_r2o = 2
        itdoc_lc13  = -1
        itdoc_hc13  = -1
        itdoc_lc14  = -1
        itdoc_hc14  = -1
      endif
    else
      i_r2o = 0
      itdoc_lc  = -1
      itdoc_hc  = -1
    endif
    if (use_DOMclasses) then
      i_dom=3
      idocsl=i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+i_extn+i_pref+i_shelfage+i_r2o+1
      idocsr=i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+i_extn+i_pref+i_shelfage+i_r2o+2
      idocr =i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+i_extn+i_pref+i_shelfage+i_r2o+3
    else
      i_dom=0
      idocsl=-1
      idocr =-1
      idocsr=-1
    endif
    if (use_DOMclasses .and. use_pref_tracers) then
      i_prefdom = 4
      iprefdoc=i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+i_extn+i_pref+i_shelfage+i_r2o+i_dom+1
      iprefdocsl=i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+i_extn+i_pref+i_shelfage+i_r2o+i_dom+2
      iprefdocsr=i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+i_extn+i_pref+i_shelfage+i_r2o+i_dom+3
      iprefdocr=i_base+i_iso+i_cfc+i_agg+i_nat_dic+i_bromo+i_extn+i_pref+i_shelfage+i_r2o+i_dom+4
    else
      i_prefdom =  0
      iprefdoc  = -1
      iprefdocsl= -1
      iprefdocsr= -1
      iprefdocr = -1
    endif

    ! total number of advected tracers
    nocetra=i_base+i_iso+i_cfc+i_agg+i_nat_dic +i_bromo+i_extn+i_pref+i_shelfage+i_r2o+i_dom+i_prefdom

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
    endif
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
    endif
    if (use_natDIC) then
      i_ndic_atm = 1
      iatmnco2 = i_base_atm+i_iso_atm+i_cfc_atm+1
    else
      i_ndic_atm = 0
      iatmnco2 = -1
    endif
    if (use_BROMO) then
      i_bromo_atm=1
      iatmbromo=i_base_atm+i_iso_atm+i_cfc_atm+ i_ndic_atm+1
    else
      i_bromo_atm=0
      iatmbromo=-1
    endif
    if (use_extNcycle) then
      i_nh3_atm = 1
      iatmnh3 = i_base_atm+i_iso_atm+i_cfc_atm+ i_ndic_atm+i_bromo_atm+1
    else
      i_nh3_atm = 0
      iatmnh3 = -1
    endif

    ! total number of atmosphere tracers
    natm=i_base_atm+i_iso_atm+i_cfc_atm+i_ndic_atm+i_bromo_atm+i_nh3_atm

    ! N-deposition
    if (use_extNcycle) then
      nndep   = 2
      idepnoy = 1
      idepnhx = 2
    else
      nndep   = 1
      idepnoy = 1
      idepnhx = -1
    endif

    ! Dust/fe
    ndust  = 2
    itdust = 1
    isfe   = 2

    ! rivers
    nriv   =8
    irdin  =1
    irdip  =2
    irsi   =3
    iralk  =4
    iriron =5
    irdoc  =6
    irtdoc =7
    irdet  =8

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
    endif

    !NOTE: The age tracer currently always needs to be the last sedlay tracer
    !      - otherwise issues in mo_sedshi.F90!
    if (use_sediment_quality) then
      i_sed_age   = 1
      issso12_age = i_sed_base + i_sed_cisonew +1
    else
      i_sed_age   = 0
      issso12_age = -1
    endif
    nsedtra_woage = i_sed_base + i_sed_cisonew ! needed in mo_sedshi.F90
    nsedtra       = nsedtra_woage + i_sed_age

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
    endif
    if (use_extNcycle) then
      i_pow_extNcycle = 3
      ipownh4 = i_pow_base + i_pow_cisonew+1
      ipown2o = i_pow_base + i_pow_cisonew+2
      ipowno2 = i_pow_base + i_pow_cisonew+3
    else
      i_pow_extNcycle = 0
      ipownh4 = -1
      ipown2o = -1
      ipowno2 = -1
    endif

    npowtra = i_pow_base + i_pow_cisonew+i_pow_extNcycle

    allocate(map_por2octra(-1:npowtra))

  end subroutine init_indices

end module mo_param1_bgc
