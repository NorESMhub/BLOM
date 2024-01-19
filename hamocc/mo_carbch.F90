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

  implicit none
  private

  ! Routines

  public :: alloc_mem_carbch ! Allocate memory for inorganic carbon variables

  ! Module variables

  real, dimension (:,:,:,:), allocatable, public :: ocetra
  real, dimension (:,:,:),   allocatable, public :: atm
  real, dimension (:,:,:),   allocatable, public :: atmflx
  real, dimension (:,:),     allocatable, public :: ndepnoyflx
  real, dimension (:,:),     allocatable, public :: ndepnhxflx
  real, dimension (:,:),     allocatable, public :: oalkflx
  real, dimension (:,:,:),   allocatable, public :: rivinflx
  real, dimension (:,:,:),   allocatable, public :: co3
  real, dimension (:,:,:),   allocatable, public :: co2star
  real, dimension (:,:,:),   allocatable, public :: hi
  real, dimension (:,:,:),   allocatable, public :: omegaa
  real, dimension (:,:,:),   allocatable, public :: omegac
  real, dimension (:,:,:),   allocatable, public :: keqb

  real, dimension (:,:,:),   allocatable, public :: satoxy
  real, dimension (:,:),     allocatable, public :: satn2o
  real, dimension (:,:),     allocatable, public :: pn2om
  real, dimension (:,:),     allocatable, public :: pnh3
  real, dimension (:,:),     allocatable, public :: atdifv
  real, dimension (:,:),     allocatable, public :: suppco2
  real, dimension (:,:,:),   allocatable, public :: sedfluxo
  real, dimension (:,:,:),   allocatable, public :: sedfluxb

  real, dimension (:,:),     allocatable, public :: pco2d
  real, dimension (:,:),     allocatable, public :: pco2m
  real, dimension (:,:),     allocatable, public :: kwco2sol
  real, dimension (:,:),     allocatable, public :: kwco2d
  real, dimension (:,:),     allocatable, public :: co2sold
  real, dimension (:,:),     allocatable, public :: co2solm
  real, dimension (:,:),     allocatable, public :: co2fxd
  real, dimension (:,:),     allocatable, public :: co2fxu
  real, dimension (:,:),     allocatable, public :: co213fxd
  real, dimension (:,:),     allocatable, public :: co213fxu
  real, dimension (:,:),     allocatable, public :: co214fxd
  real, dimension (:,:),     allocatable, public :: co214fxu
  real, dimension (:,:),     allocatable, public :: natpco2d
  real, dimension (:,:,:),   allocatable, public :: nathi
  real, dimension (:,:,:),   allocatable, public :: natco3
  real, dimension (:,:,:),   allocatable, public :: natomegaa
  real, dimension (:,:,:),   allocatable, public :: natomegac

  real, public :: atm_co2
  real, public :: atm_cfc11_nh, atm_cfc11_sh
  real, public :: atm_cfc12_nh, atm_cfc12_sh
  real, public :: atm_sf6_nh, atm_sf6_sh

contains

  subroutine alloc_mem_carbch(kpie,kpje,kpke)

    !--------------------------------------------
    ! Allocate variables in this module
    !--------------------------------------------

    use mod_xc,         only: mnproc
    use mo_control_bgc, only: io_stdo_bgc
    use mo_param1_bgc,  only: nocetra,npowtra,nsedtra,natm,nriv
    use mo_control_bgc, only: use_natDIC,use_cisonew

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
    ocetra(:,:,:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable hi ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
    endif
    allocate (hi(kpie,kpje,kpke),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory hi'
    hi(:,:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable co3 ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
    endif
    allocate (co3(kpie,kpje,kpke),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory co3'
    co3(:,:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable co2star ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
    endif
    allocate (co2star(kpie,kpje,kpke),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory co2star'
    co2star(:,:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable OmegaA, OmegaC ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
    endif
    allocate (OmegaA(kpie,kpje,kpke),stat=errstat)
    allocate (OmegaC(kpie,kpje,kpke),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory OmegaA, OmegaC'
    OmegaA(:,:,:) = 0.0
    OmegaC(:,:,:) = 0.0

    if (use_natDIC) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable natpco2d ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif
      allocate (natpco2d(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory natpco2d'
      natpco2d(:,:) = 0.0

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable nathi ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif
      allocate (nathi(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory nathi'
      nathi(:,:,:) = 0.0

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable natco3 ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif
      allocate (natco3(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory natco3'
      natco3(:,:,:) = 0.0

      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable natOmegaA, natOmegaC ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif
      allocate (natOmegaA(kpie,kpje,kpke),stat=errstat)
      allocate (natOmegaC(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory natOmegaA, natOmegaC'
      natOmegaA(:,:,:) = 0.0
      natOmegaC(:,:,:) = 0.0
    endif

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable sedfluxo ..'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',npowtra
    endif
    allocate (sedfluxo(kpie,kpje,npowtra),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory sedfluxo'
    sedfluxo(:,:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable sedfluxb ..'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',nsedtra
    endif
    allocate (sedfluxb(kpie,kpje,nsedtra),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory sedfluxb'
    sedfluxo(:,:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable satn2o ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (satn2o(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory satn2o'
    satn2o(:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable pn2om ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (pn2om(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory pn2om'
    pn2om(:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable keqb ...'
      write(io_stdo_bgc,*)'First dimension    : ',11
      write(io_stdo_bgc,*)'Second dimension   : ',kpie
      write(io_stdo_bgc,*)'Third dimension    : ',kpje
    endif
    allocate (keqb(11,kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory keqb'
    keqb(:,:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable satoxy ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
    endif
    allocate (satoxy(kpie,kpje,kpke),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory satoxy'
    satoxy(:,:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable atm ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',natm
    endif
    allocate (atm(kpie,kpje,natm),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory atm'
    atm(:,:,:) = 0.0


    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable atmflx ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',natm
    endif
    allocate (atmflx(kpie,kpje,natm),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory atmflx'
    atmflx(:,:,:) = 0.0

    ! Allocate field to hold N-deposition fluxes per timestep for
    ! inventory calculations and output
    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable ndepnoyflx ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (ndepnoyflx(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory ndepfnoylx'
    ndepnoyflx(:,:) = 0.0

    ! Allocate field to hold OA alkalinity fluxes per timestep for
    ! inventory calculations and output
    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable oalkflx ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (oalkflx(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory oalkflx'
    oalkflx(:,:) = 0.0

    ! Allocate field to hold riverine fluxes per timestep for inventory calculations
    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable rivinflx ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third  dimension   : ',nriv
    endif
    allocate(rivinflx(kpie,kpje,nriv),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory rivinflx'
    rivinflx(:,:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable pco2d ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (pco2d(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory pco2d'
    pco2d(:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable pco2m ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (pco2m(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory pco2m'
    pco2m(:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable kwco2d ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (kwco2d(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory kwco2d'
    kwco2d(:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable kwco2sol ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (kwco2sol(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory co2fxd,co2fxu'
    kwco2sol(:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable co2sold ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (co2sold(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory co2sold'
    co2sold(:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable co2solm ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (co2solm(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory co2solm'
    co2solm(:,:) = 0.0

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable co2fxd, co2fxu ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (co2fxd(kpie,kpje),stat=errstat)
    allocate (co2fxu(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory co2fxd,co2fxu'
    co2fxd(:,:) = 0.0
    co2fxu(:,:) = 0.0

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
      co213fxd(:,:) = 0.0
      co213fxu(:,:) = 0.0
      co214fxd(:,:) = 0.0
      co214fxu(:,:) = 0.0
    endif

#ifdef extNcycle
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable pnh3 ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif
      allocate (pnh3(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory pnh3'
      pnh3(:,:) = 0.0
      
      ! Allocate field to hold N-deposition NHx fluxes per timestep for inventory caluclations
      if (mnproc.eq.1) then
       write(io_stdo_bgc,*)'Memory allocation for variable ndepnhxflx ...'
       write(io_stdo_bgc,*)'First dimension    : ',kpie
       write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif
      allocate (ndepnhxflx(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory ndepnhxflx'
      ndepnhxflx(:,:) = 0.0
#endif

  end subroutine alloc_mem_carbch
  !*************************************************************************************************

end module mo_carbch
