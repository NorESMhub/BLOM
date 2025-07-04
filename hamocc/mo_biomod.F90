! Copyright (C) 2001  S. Legutke
! Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger, I. Kriest
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


module mo_biomod
  !*************************************************************************************************
  !  Variables for marine biology (declaration and memory allocation).
  !
  ! S.Legutke,        *MPI-MaD, HH*    31.10.01
  !
  ! Modified
  ! I. Kriest, GEOMAR, 11.08.2016
  ! - included T-dependence of cyanobacteria growth
  ! - modified stoichiometry for denitrification
  ! J.Schwinger,      *Uni Research, Bergen*   2018-04-12
  ! - moved accumulation of all output fields to seperate subroutine,
  !   new global fields for output defined here
  !*************************************************************************************************

  use mo_kind,        only: rp

  implicit none
  private

  ! Routines

  public :: alloc_mem_biomod ! Allocate memory for biomod variables

  ! Module variables

  real(rp), dimension (:,:),   allocatable, public :: strahl
  real(rp), dimension (:,:),   allocatable, public :: expoor
  real(rp), dimension (:,:),   allocatable, public :: expoca
  real(rp), dimension (:,:),   allocatable, public :: exposi
  real(rp), dimension (:,:),   allocatable, public :: intphosy
  real(rp), dimension (:,:),   allocatable, public :: intdnit
  real(rp), dimension (:,:),   allocatable, public :: intnfix
  real(rp), dimension (:,:),   allocatable, public :: intdmsprod
  real(rp), dimension (:,:),   allocatable, public :: intdms_bac
  real(rp), dimension (:,:),   allocatable, public :: intdms_uv
  real(rp), dimension (:,:),   allocatable, public :: carflx0100
  real(rp), dimension (:,:),   allocatable, public :: carflx0500
  real(rp), dimension (:,:),   allocatable, public :: carflx1000
  real(rp), dimension (:,:),   allocatable, public :: carflx2000
  real(rp), dimension (:,:),   allocatable, public :: carflx4000
  real(rp), dimension (:,:),   allocatable, public :: carflx_bot
  real(rp), dimension (:,:),   allocatable, public :: bsiflx0100
  real(rp), dimension (:,:),   allocatable, public :: bsiflx0500
  real(rp), dimension (:,:),   allocatable, public :: bsiflx1000
  real(rp), dimension (:,:),   allocatable, public :: bsiflx2000
  real(rp), dimension (:,:),   allocatable, public :: bsiflx4000
  real(rp), dimension (:,:),   allocatable, public :: bsiflx_bot
  real(rp), dimension (:,:),   allocatable, public :: calflx0100
  real(rp), dimension (:,:),   allocatable, public :: calflx0500
  real(rp), dimension (:,:),   allocatable, public :: calflx1000
  real(rp), dimension (:,:),   allocatable, public :: calflx2000
  real(rp), dimension (:,:),   allocatable, public :: calflx4000
  real(rp), dimension (:,:),   allocatable, public :: calflx_bot
  real(rp), dimension (:,:),   allocatable, public :: dustflx0100
  real(rp), dimension (:,:),   allocatable, public :: dustflx0500
  real(rp), dimension (:,:),   allocatable, public :: dustflx1000
  real(rp), dimension (:,:),   allocatable, public :: dustflx2000
  real(rp), dimension (:,:),   allocatable, public :: dustflx4000
  real(rp), dimension (:,:),   allocatable, public :: dustflx_bot
  real(rp), dimension (:,:,:), allocatable, public :: phosy3d
  real(rp), dimension (:,:),   allocatable, public :: int_exudl
  real(rp), dimension (:,:),   allocatable, public :: int_exudsl
  real(rp), dimension (:,:),   allocatable, public :: int_excrl
  real(rp), dimension (:,:),   allocatable, public :: int_excrsl
  real(rp), dimension (:,:),   allocatable, public :: int_docl_rem
  real(rp), dimension (:,:),   allocatable, public :: int_docsl_rem
  real(rp), dimension (:,:),   allocatable, public :: int_docr_rem
  real(rp), dimension (:,:),   allocatable, public :: int_docsr_rem

  ! Variables for interactive phytoplanktion absorption (use_FB_BGC_OCE=.true.)
  real(rp), dimension (:,:,:), allocatable, public :: abs_oce

  ! Variables for aggregation scheme (use_AGG=.true.)
  real(rp), dimension (:,:,:), allocatable, public  :: wmass
  real(rp), dimension (:,:,:), allocatable, public  :: wnumb
  real(rp), dimension (:,:,:), allocatable, public  :: eps3d
  real(rp), dimension (:,:,:), allocatable, public  :: asize3d

  ! Variables for bromoform scheme (use_BROMO=.true.)
  real(rp), dimension (:,:),   allocatable, public  :: int_chbr3_prod
  real(rp), dimension (:,:),   allocatable, public  :: int_chbr3_uv

  real(rp), dimension (:,:,:), allocatable, public :: nitr_NH4
  real(rp), dimension (:,:,:), allocatable, public :: nitr_NO2
  real(rp), dimension (:,:,:), allocatable, public :: nitr_N2O_prod
  real(rp), dimension (:,:,:), allocatable, public :: nitr_NH4_OM
  real(rp), dimension (:,:,:), allocatable, public :: nitr_NO2_OM
  real(rp), dimension (:,:,:), allocatable, public :: denit_NO3
  real(rp), dimension (:,:,:), allocatable, public :: denit_NO2
  real(rp), dimension (:,:,:), allocatable, public :: denit_N2O
  real(rp), dimension (:,:,:), allocatable, public :: DNRA_NO2
  real(rp), dimension (:,:,:), allocatable, public :: anmx_N2_prod
  real(rp), dimension (:,:,:), allocatable, public :: anmx_OM_prod
  real(rp), dimension (:,:,:), allocatable, public :: phosy_NH4
  real(rp), dimension (:,:,:), allocatable, public :: phosy_NO3
  real(rp), dimension (:,:,:), allocatable, public :: remin_aerob
  real(rp), dimension (:,:,:), allocatable, public :: remin_sulf

CONTAINS

  subroutine alloc_mem_biomod(kpie,kpje,kpke)
    !******************************************************************************
    ! ALLOC_MEM_BIOMOD - Allocate variables in this module
    !******************************************************************************
    use mod_xc,         only: mnproc
    use mo_control_bgc, only: io_stdo_bgc
    use mo_control_bgc, only: use_FB_BGC_OCE,use_AGG,use_BROMO,use_extNcycle,use_DOMclasses

    ! Arguments
    integer, intent(in) :: kpie
    integer, intent(in) :: kpje
    integer, intent(in) :: kpke

    ! local variables
    integer :: errstat

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)' '
      write(io_stdo_bgc,*)'***************************************************'
      write(io_stdo_bgc,*)'Memory allocation for marine biology module :'
      write(io_stdo_bgc,*)' '
    endif

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable strahl ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif

    allocate (strahl(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory strahl'
    strahl(:,:) = 0.0_rp

    if (use_FB_BGC_OCE ) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable abs_oce'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif

      allocate (abs_oce(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory abs_oce'
      abs_oce(:,:,:) = 0.0_rp
    endif

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable expoor ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif

    allocate (expoor(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory expoor'
    expoor(:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable expoca ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif

    allocate (expoca(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory expoca'
    expoca(:,:) = 0.0_rp


    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable exposi ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif

    allocate (exposi(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory exposi'
    exposi(:,:) = 0.0_rp


    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable intphosy ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif

    allocate (intphosy(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory intphosy'
    intphosy(:,:) = 0.0_rp


    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable intdnit ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif

    allocate (intdnit(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory intdnit'
    intdnit(:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable intnfix ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif

    allocate (intnfix(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory intnfix'
    intnfix(:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable intdmsprod, intdms_bac, intdms_uv ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif

    allocate (intdmsprod(kpie,kpje),stat=errstat)
    allocate (intdms_bac(kpie,kpje),stat=errstat)
    allocate (intdms_uv(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory intdmsprod, intdms_bac, intdms_uv'
    intdmsprod(:,:) = 0.0_rp
    intdms_bac(:,:) = 0.0_rp
    intdms_uv(:,:)  = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable carflx* ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif

    allocate (carflx0100(kpie,kpje),stat=errstat)
    allocate (carflx0500(kpie,kpje),stat=errstat)
    allocate (carflx1000(kpie,kpje),stat=errstat)
    allocate (carflx2000(kpie,kpje),stat=errstat)
    allocate (carflx4000(kpie,kpje),stat=errstat)
    allocate (carflx_bot(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory carflx*'
    carflx0100(:,:) = 0.0_rp
    carflx0500(:,:) = 0.0_rp
    carflx1000(:,:) = 0.0_rp
    carflx2000(:,:) = 0.0_rp
    carflx4000(:,:) = 0.0_rp
    carflx_bot(:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable bsiflx* ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif

    allocate (bsiflx0100(kpie,kpje),stat=errstat)
    allocate (bsiflx0500(kpie,kpje),stat=errstat)
    allocate (bsiflx1000(kpie,kpje),stat=errstat)
    allocate (bsiflx2000(kpie,kpje),stat=errstat)
    allocate (bsiflx4000(kpie,kpje),stat=errstat)
    allocate (bsiflx_bot(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory bsiflx*'
    bsiflx0100(:,:) = 0.0_rp
    bsiflx0500(:,:) = 0.0_rp
    bsiflx1000(:,:) = 0.0_rp
    bsiflx2000(:,:) = 0.0_rp
    bsiflx4000(:,:) = 0.0_rp
    bsiflx_bot(:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable calflx* ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif

    allocate (calflx0100(kpie,kpje),stat=errstat)
    allocate (calflx0500(kpie,kpje),stat=errstat)
    allocate (calflx1000(kpie,kpje),stat=errstat)
    allocate (calflx2000(kpie,kpje),stat=errstat)
    allocate (calflx4000(kpie,kpje),stat=errstat)
    allocate (calflx_bot(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory calflx*'
    calflx0100(:,:) = 0.0_rp
    calflx0500(:,:) = 0.0_rp
    calflx1000(:,:) = 0.0_rp
    calflx2000(:,:) = 0.0_rp
    calflx4000(:,:) = 0.0_rp
    calflx_bot(:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable dustflx* ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif

    allocate (dustflx0100(kpie,kpje),stat=errstat)
    allocate (dustflx0500(kpie,kpje),stat=errstat)
    allocate (dustflx1000(kpie,kpje),stat=errstat)
    allocate (dustflx2000(kpie,kpje),stat=errstat)
    allocate (dustflx4000(kpie,kpje),stat=errstat)
    allocate (dustflx_bot(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory dustflx*'
    dustflx0100(:,:) = 0.0_rp
    dustflx0500(:,:) = 0.0_rp
    dustflx1000(:,:) = 0.0_rp
    dustflx2000(:,:) = 0.0_rp
    dustflx4000(:,:) = 0.0_rp
    dustflx_bot(:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable phosy3d ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
    endif

    allocate (phosy3d(kpie,kpje,kpke),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory phosy3d'
    phosy3d(:,:,:) = 0.0_rp

    if (use_AGG) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable wmass ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif

      allocate (wmass(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory eps3d'
      wmass(:,:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable wnumb ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif

      allocate (wnumb(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory eps3d'
      wnumb(:,:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable eps3d ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif

      allocate (eps3d(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory eps3d'
      eps3d(:,:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable asize3d ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif

      allocate (asize3d(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory asize3d'
      asize3d(:,:,:) = 0.0_rp
    endif

    if (use_BROMO) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable int_chbr3_prod, int_chbr3_uv ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      allocate (int_chbr3_prod(kpie,kpje),stat=errstat)
      allocate (int_chbr3_uv(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory int_chbr3_prod, int_chbr3_uv'
      int_chbr3_prod(:,:) = 0.0_rp
      int_chbr3_uv(:,:) = 0.0_rp
    endif

    if (use_extNcycle) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable of the extended nitrogen cycle ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif

      allocate (nitr_NH4(kpie,kpje,kpke),stat=errstat)
      allocate (nitr_NO2(kpie,kpje,kpke),stat=errstat)
      allocate (nitr_N2O_prod(kpie,kpje,kpke),stat=errstat)
      allocate (nitr_NH4_OM(kpie,kpje,kpke),stat=errstat)
      allocate (nitr_NO2_OM(kpie,kpje,kpke),stat=errstat)
      allocate (denit_NO3(kpie,kpje,kpke),stat=errstat)
      allocate (denit_NO2(kpie,kpje,kpke),stat=errstat)
      allocate (denit_N2O(kpie,kpje,kpke),stat=errstat)
      allocate (DNRA_NO2(kpie,kpje,kpke),stat=errstat)
      allocate (anmx_N2_prod(kpie,kpje,kpke),stat=errstat)
      allocate (anmx_OM_prod(kpie,kpje,kpke),stat=errstat)
      allocate (phosy_NH4(kpie,kpje,kpke),stat=errstat)
      allocate (phosy_NO3(kpie,kpje,kpke),stat=errstat)
      allocate (remin_aerob(kpie,kpje,kpke),stat=errstat)
      allocate (remin_sulf(kpie,kpje,kpke),stat=errstat)

      if(errstat.ne.0) stop 'not enough memory extended nitrogen cycle'
      nitr_NH4      = 0._rp
      nitr_NO2      = 0._rp
      nitr_N2O_prod = 0._rp
      nitr_NH4_OM   = 0._rp
      nitr_NO2_OM   = 0._rp
      denit_NO3     = 0._rp
      denit_NO2     = 0._rp
      denit_N2O     = 0._rp
      DNRA_NO2      = 0._rp
      anmx_N2_prod  = 0._rp
      anmx_OM_prod  = 0._rp
      phosy_NH4     = 0._rp
      phosy_NO3     = 0._rp
      remin_aerob   = 0._rp
      remin_sulf    = 0._rp
    endif

    if (use_DOMclasses) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable int_exudl ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      allocate (int_exudl(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory int_exudl'
      int_exudl(:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable int_exudsl ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      allocate (int_exudsl(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory int_exudsl'
      int_exudsl(:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable int_excrl ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      allocate (int_excrl(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory int_excrl'
      int_excrl(:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable int_excrsl ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      allocate (int_excrsl(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory intexcrsl'
      int_excrsl(:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable int_docl_rem ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      allocate (int_docl_rem(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory int_docl_rem'
      int_docl_rem(:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable int_docsl_rem ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      allocate (int_docsl_rem(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory int_docsl_rem'
      int_docsl_rem(:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable int_docsr_rem ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      allocate (int_docsr_rem(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory int_docsr_rem'
      int_docsr_rem(:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable int_docr_rem ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      allocate (int_docr_rem(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory int_docr_rem'
      int_docr_rem(:,:) = 0.0_rp
    endif

  end subroutine alloc_mem_biomod

end module mo_biomod
