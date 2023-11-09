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

  !******************************************************************************
  ! bgc tracer parameters.
  ! - definition of indices in tracer arrays
  !
  !  Patrick Wetzel    *MPI-Met, HH*    01.09.03
  !  Modified
  !  J.Schwinger,        *NORCE Climate, Bergen*    2020-05-26
  !******************************************************************************

  use mo_control_bgc, only: use_BROMO, use_AGG, use_WLIN, use_natDIC, use_CFC,         &
                            use_cisonew, use_PBGC_OCNP_TIMESTEP, use_PBGC_CK_TIMESTEP, &
                            use_FB_BGC_OCE, use_BOXATM, use_sedbypass
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
  integer, protected :: iprefo2
  integer, protected :: iprefpo4
  integer, protected :: iprefalk
  integer, protected :: iprefdic
  integer, protected :: idicsat

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

  integer, protected :: natm ! total number of atmosphere tracers

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
  integer, protected :: nsedtra

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
  end subroutine init_por2octra_mapping

  ! ===========================================================================
  subroutine init_indices()

    use mod_xc        , only: lp, mnproc
    use mo_control_bgc, only: bgc_namelist,get_bgc_namelist, io_stdo_bgc
    use mo_control_bgc, only: use_BROMO, use_AGG, use_WLIN, use_natDIC, use_CFC,  &
         use_cisonew, use_sedbypass,                         &
         use_PBGC_OCNP_TIMESTEP, use_PBGC_CK_TIMESTEP,       &
         use_FB_BGC_OCE, use_BOXATM
    integer :: iounit

    namelist / config_bgc / use_BROMO, use_AGG, use_WLIN, &
         use_natDIC, use_CFC, use_cisonew, use_sedbypass, use_PBGC_OCNP_TIMESTEP, &
         use_PBGC_CK_TIMESTEP, use_FB_BGC_OCE, use_BOXATM

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
    endif
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
    endif
    npowtra = i_pow_base + i_pow_cisonew

    allocate(map_por2octra(-1:npowtra))

  end subroutine init_indices

end module mo_param1_bgc
