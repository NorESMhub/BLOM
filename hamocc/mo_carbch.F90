! Copyright (C) 2002  S. Legutke, P. Wetzel
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


module mo_carbch

  !*************************************************************************************************
  ! Variables for inorganic carbon cycle (declaration and memory allocation)
  !
  !  S.Legutke,        *MPI-MaD, HH*     31.10.01
  !
  !  Modified
  !  Patrick Wetzel    *MPI-Met, HH*     16.04.02
  !  - new: atm, atdifv, suppco2
  !  - changed: chemc(:,:,:) to chemcm(:,:,:,:)
  !  - new: bgcmean(:,:,:,:)
  !  J. Schwinger      *UiB-GfI, Bergen* 04.05.12
  !  - added initialisation of all vars after allocation
  !  J.Schwinger,      *Uni Research, Bergen*   2018-04-12
  !  - moved accumulation of all output fields to seperate subroutine,
  !    new global fields for output defined here
  !  - added OmegaA
  !*************************************************************************************************

  use mo_kind,        only: rp

  implicit none
  private

  ! Routines

  public :: alloc_mem_carbch ! Allocate memory for inorganic carbon variables

  ! Module variables

  real(rp), dimension (:,:,:,:), allocatable, public :: ocetra
  real(rp), dimension (:,:,:),   allocatable, public :: atm
  real(rp), dimension (:,:,:),   allocatable, public :: atmflx
  real(rp), dimension (:,:),     allocatable, public :: ndepnoyflx
  real(rp), dimension (:,:),     allocatable, public :: ndepnhxflx
  real(rp), dimension (:,:),     allocatable, public :: oalkflx
  real(rp), dimension (:,:,:),   allocatable, public :: dustflx
  real(rp), dimension (:,:,:),   allocatable, public :: rivinflx
  real(rp), dimension (:,:,:),   allocatable, public :: co3
  real(rp), dimension (:,:,:),   allocatable, public :: co2star
  real(rp), dimension (:,:,:),   allocatable, public :: hi
  real(rp), dimension (:,:,:),   allocatable, public :: omegaa
  real(rp), dimension (:,:,:),   allocatable, public :: omegac
  real(rp), dimension (:,:,:),   allocatable, public :: keqb

  real(rp), dimension (:,:,:),   allocatable, public :: satoxy
  real(rp), dimension (:,:),     allocatable, public :: satn2o
  real(rp), dimension (:,:),     allocatable, public :: pn2om
  real(rp), dimension (:,:),     allocatable, public :: pnh3
  real(rp), dimension (:,:),     allocatable, public :: atdifv
  real(rp), dimension (:,:),     allocatable, public :: suppco2
  real(rp), dimension (:,:,:),   allocatable, public :: sedfluxo
  real(rp), dimension (:,:,:),   allocatable, public :: sedfluxb
  real(rp), dimension (:,:,:,:), allocatable, public :: nutlim_diag
  real(rp), dimension (:,:,:),   allocatable, public :: zeu_nutlim_diag

  real(rp), dimension (:,:),     allocatable, public :: fco2
  real(rp), dimension (:,:),     allocatable, public :: pco2
  real(rp), dimension (:,:),     allocatable, public :: xco2
  real(rp), dimension (:,:),     allocatable, public :: pco2_gex
  real(rp), dimension (:,:),     allocatable, public :: kwco2sol
  real(rp), dimension (:,:),     allocatable, public :: kwco2a
  real(rp), dimension (:,:),     allocatable, public :: co2sol
  real(rp), dimension (:,:),     allocatable, public :: co2fxd
  real(rp), dimension (:,:),     allocatable, public :: co2fxu
  real(rp), dimension (:,:),     allocatable, public :: co213fxd
  real(rp), dimension (:,:),     allocatable, public :: co213fxu
  real(rp), dimension (:,:),     allocatable, public :: co214fxd
  real(rp), dimension (:,:),     allocatable, public :: co214fxu
  real(rp), dimension (:,:),     allocatable, public :: natpco2
  real(rp), dimension (:,:,:),   allocatable, public :: nathi
  real(rp), dimension (:,:,:),   allocatable, public :: natco3
  real(rp), dimension (:,:,:),   allocatable, public :: natomegaa
  real(rp), dimension (:,:,:),   allocatable, public :: natomegac

  real(rp), public :: atm_co2
  real(rp), public :: atm_cfc11_nh, atm_cfc11_sh
  real(rp), public :: atm_cfc12_nh, atm_cfc12_sh
  real(rp), public :: atm_sf6_nh, atm_sf6_sh

  ! Index for nutlim_diag
  integer, parameter, public :: inutlim_phosph = 1
  integer, parameter, public :: inutlim_n      = 2
  integer, parameter, public :: inutlim_fe     = 3

contains

  subroutine alloc_mem_carbch(kpie,kpje,kpke)

    !--------------------------------------------
    ! Allocate variables in this module
    !--------------------------------------------

    use mod_xc,         only: mnproc
    use mo_control_bgc, only: io_stdo_bgc
    use mo_param1_bgc,  only: nocetra,npowtra,nsedtra,natm,nriv
    use mo_control_bgc, only: use_natDIC,use_cisonew,use_extNcycle

    integer, intent(in) :: kpie
    integer, intent(in) :: kpje
    integer, intent(in) :: kpke

    ! Local variables
    integer :: errstat

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)' '
      write(io_stdo_bgc,*)'***************************************************'
      write(io_stdo_bgc,*)'Memory allocation for carbon chemistry module :'
      write(io_stdo_bgc,*)' '
    endif

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable ocetra ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
      write(io_stdo_bgc,*)'Forth dimension    : ',nocetra
    endif
    allocate (ocetra(kpie,kpje,kpke,nocetra),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory ocetra'
    ocetra(:,:,:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable hi ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
    endif
    allocate (hi(kpie,kpje,kpke),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory hi'
    hi(:,:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable co3 ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
    endif
    allocate (co3(kpie,kpje,kpke),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory co3'
    co3(:,:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable co2star ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
    endif
    allocate (co2star(kpie,kpje,kpke),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory co2star'
    co2star(:,:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable OmegaA, OmegaC ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
    endif
    allocate (OmegaA(kpie,kpje,kpke),stat=errstat)
    allocate (OmegaC(kpie,kpje,kpke),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory OmegaA, OmegaC'
    OmegaA(:,:,:) = 0.0_rp
    OmegaC(:,:,:) = 0.0_rp

    if (use_natDIC) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable natpco2 ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif
      allocate (natpco2(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory natpco2'
      natpco2(:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable nathi ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif
      allocate (nathi(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory nathi'
      nathi(:,:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable natco3 ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif
      allocate (natco3(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory natco3'
      natco3(:,:,:) = 0.0_rp

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable natOmegaA, natOmegaC ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif
      allocate (natOmegaA(kpie,kpje,kpke),stat=errstat)
      allocate (natOmegaC(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory natOmegaA, natOmegaC'
      natOmegaA(:,:,:) = 0.0_rp
      natOmegaC(:,:,:) = 0.0_rp
    endif

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable sedfluxo ..'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',npowtra
    endif
    allocate (sedfluxo(kpie,kpje,npowtra),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory sedfluxo'
    sedfluxo(:,:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable sedfluxb ..'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',nsedtra
    endif
    allocate (sedfluxb(kpie,kpje,nsedtra),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory sedfluxb'
    sedfluxb(:,:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable nutlim_diag ..'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
      write(io_stdo_bgc,*)'Fourth dimension   : ',3     ! number of potentially limiting nutrients
    endif
    allocate (nutlim_diag(kpie,kpje,kpke,3),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory nutlim_diag'
    nutlim_diag(:,:,:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable zeu_nutlim_diag ..'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',3     ! number of potentially limiting nutrients
    endif
    allocate (zeu_nutlim_diag(kpie,kpje,3),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory zeu_nutlim_diag'
    zeu_nutlim_diag(:,:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable satn2o ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (satn2o(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory satn2o'
    satn2o(:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable pn2om ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (pn2om(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory pn2om'
    pn2om(:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable keqb ...'
      write(io_stdo_bgc,*)'First dimension    : ',11
      write(io_stdo_bgc,*)'Second dimension   : ',kpie
      write(io_stdo_bgc,*)'Third dimension    : ',kpje
    endif
    allocate (keqb(11,kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory keqb'
    keqb(:,:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable satoxy ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
    endif
    allocate (satoxy(kpie,kpje,kpke),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory satoxy'
    satoxy(:,:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable atm ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',natm
    endif
    allocate (atm(kpie,kpje,natm),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory atm'
    atm(:,:,:) = 0.0_rp


    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable atmflx ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',natm
    endif
    allocate (atmflx(kpie,kpje,natm),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory atmflx'
    atmflx(:,:,:) = 0.0_rp

    ! Allocate field to hold N-deposition fluxes per timestep for
    ! inventory calculations and output
    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable ndepnoyflx ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (ndepnoyflx(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory ndepfnoylx'
    ndepnoyflx(:,:) = 0.0_rp

    ! Allocate field to hold OA alkalinity fluxes per timestep for
    ! inventory calculations and output
    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable oalkflx ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (oalkflx(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory oalkflx'
    oalkflx(:,:) = 0.0_rp

    ! Allocate field to hold riverine fluxes per timestep for inventory calculations
    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable rivinflx ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third  dimension   : ',nriv
    endif
    allocate(rivinflx(kpie,kpje,nriv),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory rivinflx'
    rivinflx(:,:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable fco2 ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (fco2(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory fco2'
    fco2(:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable pco2 ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (pco2(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory pco2'
    pco2(:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable xco2 ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (xco2(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory xco2'
    xco2(:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable pco2_gex ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (pco2_gex(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory pco2_gex'
    pco2_gex(:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable kwco2a ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (kwco2a(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory kwco2a'
    kwco2a(:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable kwco2sol ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (kwco2sol(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory co2fxd,co2fxu'
    kwco2sol(:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable co2sol ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (co2sol(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory co2sold'
    co2sol(:,:) = 0.0_rp

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable co2fxd, co2fxu ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (co2fxd(kpie,kpje),stat=errstat)
    allocate (co2fxu(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory co2fxd,co2fxu'
    co2fxd(:,:) = 0.0_rp
    co2fxu(:,:) = 0.0_rp

    if (use_cisonew) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable co213fxd,..., co214fxu ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif
      allocate (co213fxd(kpie,kpje),stat=errstat)
      allocate (co213fxu(kpie,kpje),stat=errstat)
      allocate (co214fxd(kpie,kpje),stat=errstat)
      allocate (co214fxu(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory co213fxd,..., co214fxu'
      co213fxd(:,:) = 0.0_rp
      co213fxu(:,:) = 0.0_rp
      co214fxd(:,:) = 0.0_rp
      co214fxu(:,:) = 0.0_rp
    endif

    if (use_extNcycle) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable pnh3 ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif
      allocate (pnh3(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory pnh3'
      pnh3(:,:) = 0.0_rp

      ! Allocate field to hold N-deposition NHx fluxes per timestep for inventory caluclations
      if (mnproc.eq.1) then
       write(io_stdo_bgc,*)'Memory allocation for variable ndepnhxflx ...'
       write(io_stdo_bgc,*)'First dimension    : ',kpie
       write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif
      allocate (ndepnhxflx(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory ndepnhxflx'
      ndepnhxflx(:,:) = 0.0_rp
    endif

  end subroutine alloc_mem_carbch
  !*************************************************************************************************

end module mo_carbch
