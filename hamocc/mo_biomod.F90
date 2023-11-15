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

  implicit none
  private

  ! Routines

  public :: alloc_mem_biomod ! Allocate memory for biomod variables

  ! Module variables

  real, dimension (:,:),   allocatable, public :: strahl
  real, dimension (:,:),   allocatable, public :: expoor
  real, dimension (:,:),   allocatable, public :: expoca
  real, dimension (:,:),   allocatable, public :: exposi
  real, dimension (:,:),   allocatable, public :: intphosy
  real, dimension (:,:),   allocatable, public :: intdnit
  real, dimension (:,:),   allocatable, public :: intnfix
  real, dimension (:,:),   allocatable, public :: intdmsprod
  real, dimension (:,:),   allocatable, public :: intdms_bac
  real, dimension (:,:),   allocatable, public :: intdms_uv
  real, dimension (:,:),   allocatable, public :: carflx0100
  real, dimension (:,:),   allocatable, public :: carflx0500
  real, dimension (:,:),   allocatable, public :: carflx1000
  real, dimension (:,:),   allocatable, public :: carflx2000
  real, dimension (:,:),   allocatable, public :: carflx4000
  real, dimension (:,:),   allocatable, public :: carflx_bot
  real, dimension (:,:),   allocatable, public :: bsiflx0100
  real, dimension (:,:),   allocatable, public :: bsiflx0500
  real, dimension (:,:),   allocatable, public :: bsiflx1000
  real, dimension (:,:),   allocatable, public :: bsiflx2000
  real, dimension (:,:),   allocatable, public :: bsiflx4000
  real, dimension (:,:),   allocatable, public :: bsiflx_bot
  real, dimension (:,:),   allocatable, public :: calflx0100
  real, dimension (:,:),   allocatable, public :: calflx0500
  real, dimension (:,:),   allocatable, public :: calflx1000
  real, dimension (:,:),   allocatable, public :: calflx2000
  real, dimension (:,:),   allocatable, public :: calflx4000
  real, dimension (:,:),   allocatable, public :: calflx_bot
  real, dimension (:,:,:), allocatable, public :: phosy3d

  ! Variables for interactive phytoplanktion absorption (use_FB_BGC_OCE=.true.)
  real, dimension (:,:,:), allocatable, public :: abs_oce

  ! Variables for aggregation scheme (use_AGG=.true.)
  real, dimension (:,:,:), allocatable, public  :: wmass
  real, dimension (:,:,:), allocatable, public  :: wnumb
  real, dimension (:,:,:), allocatable, public  :: eps3d
  real, dimension (:,:,:), allocatable, public  :: asize3d

  ! Variables for bromoform scheme (use_BROMO=.true.)
  real, dimension (:,:),   allocatable, public  :: int_chbr3_prod
  real, dimension (:,:),   allocatable, public  :: int_chbr3_uv

  real, public :: growth_co2
  real, public :: bifr13_perm

CONTAINS

  subroutine alloc_mem_biomod(kpie,kpje,kpke)
    !******************************************************************************
    ! ALLOC_MEM_BIOMOD - Allocate variables in this module
    !******************************************************************************
    use mod_xc,         only: mnproc
    use mo_control_bgc, only: io_stdo_bgc
    use mo_control_bgc, only: use_FB_BGC_OCE,use_AGG,use_BROMO

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
    strahl(:,:) = 0.0

    if (use_FB_BGC_OCE ) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable abs_oce'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif

      allocate (abs_oce(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory abs_oce'
      abs_oce(:,:,:) = 0.0
    endif

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable expoor ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif

    allocate (expoor(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory expoor'
    expoor(:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable expoca ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif

    allocate (expoca(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory expoca'
    expoca(:,:) = 0.0


    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable exposi ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif

    allocate (exposi(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory exposi'
    exposi(:,:) = 0.0


    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable intphosy ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif

    allocate (intphosy(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory intphosy'
    intphosy(:,:) = 0.0


    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable intdnit ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif

    allocate (intdnit(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory intdnit'
    intdnit(:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable intnfix ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif

    allocate (intnfix(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory intnfix'
    intnfix(:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable intdmsprod, intdms_bac, intdms_uv ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif

    allocate (intdmsprod(kpie,kpje),stat=errstat)
    allocate (intdms_bac(kpie,kpje),stat=errstat)
    allocate (intdms_uv(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory intdmsprod, intdms_bac, intdms_uv'
    intdmsprod(:,:) = 0.0
    intdms_bac(:,:) = 0.0
    intdms_uv(:,:) = 0.0

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
    carflx0100(:,:) = 0.0
    carflx0500(:,:) = 0.0
    carflx1000(:,:) = 0.0
    carflx2000(:,:) = 0.0
    carflx4000(:,:) = 0.0
    carflx_bot(:,:) = 0.0

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
    bsiflx0100(:,:) = 0.0
    bsiflx0500(:,:) = 0.0
    bsiflx1000(:,:) = 0.0
    bsiflx2000(:,:) = 0.0
    bsiflx4000(:,:) = 0.0
    bsiflx_bot(:,:) = 0.0

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
    if(errstat.ne.0) stop 'not enough memory bsiflx*'
    calflx0100(:,:) = 0.0
    calflx0500(:,:) = 0.0
    calflx1000(:,:) = 0.0
    calflx2000(:,:) = 0.0
    calflx4000(:,:) = 0.0
    calflx_bot(:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable phosy3d ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
    endif

    allocate (phosy3d(kpie,kpje,kpke),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory phosy3d'
    phosy3d(:,:,:) = 0.0

    if (use_AGG) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable wmass ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif

      allocate (wmass(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory eps3d'
      wmass(:,:,:) = 0.0

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable wnumb ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif

      allocate (wnumb(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory eps3d'
      wnumb(:,:,:) = 0.0

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable eps3d ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif

      allocate (eps3d(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory eps3d'
      eps3d(:,:,:) = 0.0

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable asize3d ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif

      allocate (asize3d(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory asize3d'
      asize3d(:,:,:) = 0.0
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
      int_chbr3_prod(:,:) = 0.0
      int_chbr3_uv(:,:) = 0.0
    endif

  end subroutine alloc_mem_biomod

end module mo_biomod
