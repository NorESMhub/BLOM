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


      MODULE mo_carbch
!***********************************************************************
!
! MODULE mo_carbch - Variables for inorganic carbon cycle.
!
!     S.Legutke,        *MPI-MaD, HH*     31.10.01
!
!     Modified
!     --------
!  
!     Patrick Wetzel    *MPI-Met, HH*     16.04.02
!     - new: atm, atdifv, suppco2
!     - changed: chemc(:,:,:) to chemcm(:,:,:,:)
!     - new: bgcmean(:,:,:,:)
!
!     J. Schwinger      *UiB-GfI, Bergen* 04.05.12
!     - added initialisation of all vars after allocation
!
!     J.Schwinger,      *Uni Research, Bergen*   2018-04-12
!     - moved accumulation of all output fields to seperate SUBROUTINE,
!       new global fields for output defined here
!     - added OmegaA
!
!     Purpose
!     -------
!     - declaration and memory allocation
!
!     Description:
!     ------------
!     Public routines and variable of this MODULE:
!
!     -SUBROUTINE alloc_mem_carbch
!        Allocate memory for inorganic carbon variables
!
!
!**********************************************************************
      implicit none

      real, dimension (:,:,:,:), allocatable :: ocetra
      real, dimension (:,:,:),   allocatable :: atm
      real, dimension (:,:,:),   allocatable :: atmflx
      real, dimension (:,:),     allocatable :: ndepflx
      real, dimension (:,:),     allocatable :: oalkflx
      real, dimension (:,:,:),   allocatable :: rivinflx
      real, dimension (:,:,:),   allocatable :: co3
      real, dimension (:,:,:),   allocatable :: co2star
      real, dimension (:,:,:),   allocatable :: hi
      real, dimension (:,:,:),   allocatable :: OmegaA
      real, dimension (:,:,:),   allocatable :: OmegaC
      real, dimension (:,:,:),   allocatable :: keqb

      real, dimension (:,:,:),   allocatable :: satoxy
      real, dimension (:,:),     allocatable :: satn2o
      real, dimension (:,:),     allocatable :: atdifv
      real, dimension (:,:),     allocatable :: suppco2
      real, dimension (:,:,:),   allocatable :: sedfluxo

      real, dimension (:,:),     allocatable :: pco2d
      real, dimension (:,:),     allocatable :: pco2m
      real, dimension (:,:),     allocatable :: kwco2sol
      real, dimension (:,:),     allocatable :: kwco2d
      real, dimension (:,:),     allocatable :: co2sold
      real, dimension (:,:),     allocatable :: co2solm
      real, dimension (:,:),     allocatable :: co2fxd
      real, dimension (:,:),     allocatable :: co2fxu
      real, dimension (:,:),     allocatable :: co213fxd
      real, dimension (:,:),     allocatable :: co213fxu
      real, dimension (:,:),     allocatable :: co214fxd
      real, dimension (:,:),     allocatable :: co214fxu
      real, dimension (:,:),     allocatable :: natpco2d
      real, dimension (:,:,:),   allocatable :: nathi
      real, dimension (:,:,:),   allocatable :: natco3
      real, dimension (:,:,:),   allocatable :: natOmegaA
      real, dimension (:,:,:),   allocatable :: natOmegaC
	  
      real :: atm_co2
      real :: atm_cfc11_nh,atm_cfc11_sh
      real :: atm_cfc12_nh,atm_cfc12_sh
      real :: atm_sf6_nh,atm_sf6_sh

      CONTAINS

      SUBROUTINE ALLOC_MEM_CARBCH(kpie,kpje,kpke)
!******************************************************************************
! ALLOC_MEM_CARBCH - Allocate variables in this MODULE
!******************************************************************************
      use mod_xc,         only: mnproc
      use mo_control_bgc, only: io_stdo_bgc
      use mo_param1_bgc,  only: nocetra,npowtra,natm,nriv
      use mo_control_bgc, only: use_natDIC,use_cisonew

      integer, intent(in) :: kpie,kpje,kpke
      integer             :: errstat

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

      ALLOCATE (ocetra(kpie,kpje,kpke,nocetra),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory ocetra'
      ocetra(:,:,:,:) = 0.0

      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable hi ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif

      ALLOCATE (hi(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory hi'
      hi(:,:,:) = 0.0

      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable co3 ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif

      ALLOCATE (co3(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory co3'
      co3(:,:,:) = 0.0

      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable co2star ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif

      ALLOCATE (co2star(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory co2star'
      co2star(:,:,:) = 0.0


      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable OmegaA, OmegaC ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif

      ALLOCATE (OmegaA(kpie,kpje,kpke),stat=errstat)
      ALLOCATE (OmegaC(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory OmegaA, OmegaC'
      OmegaA(:,:,:) = 0.0
      OmegaC(:,:,:) = 0.0

      if (use_natDIC) then
         if (mnproc.eq.1) then
         write(io_stdo_bgc,*)'Memory allocation for variable natpco2d ...'
         write(io_stdo_bgc,*)'First dimension    : ',kpie
         write(io_stdo_bgc,*)'Second dimension   : ',kpje
         endif

         ALLOCATE (natpco2d(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory natpco2d'
         natpco2d(:,:) = 0.0

         if (mnproc.eq.1) then
         write(io_stdo_bgc,*)'Memory allocation for variable nathi ...'
         write(io_stdo_bgc,*)'First dimension    : ',kpie
         write(io_stdo_bgc,*)'Second dimension   : ',kpje
         write(io_stdo_bgc,*)'Third dimension    : ',kpke
         endif

         ALLOCATE (nathi(kpie,kpje,kpke),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory nathi'
         nathi(:,:,:) = 0.0
         
         if (mnproc.eq.1) then
         write(io_stdo_bgc,*)'Memory allocation for variable natco3 ...'
         write(io_stdo_bgc,*)'First dimension    : ',kpie
         write(io_stdo_bgc,*)'Second dimension   : ',kpje
         write(io_stdo_bgc,*)'Third dimension    : ',kpke
         endif

         ALLOCATE (natco3(kpie,kpje,kpke),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory natco3'
         natco3(:,:,:) = 0.0
         
         if (mnproc.eq.1) then
         write(io_stdo_bgc,*)'Memory allocation for variable natOmegaA, natOmegaC ...'
         write(io_stdo_bgc,*)'First dimension    : ',kpie
         write(io_stdo_bgc,*)'Second dimension   : ',kpje
         write(io_stdo_bgc,*)'Third dimension    : ',kpke
         endif

         ALLOCATE (natOmegaA(kpie,kpje,kpke),stat=errstat)
         ALLOCATE (natOmegaC(kpie,kpje,kpke),stat=errstat)
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

      ALLOCATE (sedfluxo(kpie,kpje,npowtra),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory sedfluxo'
      sedfluxo(:,:,:) = 0.0

      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable satn2o ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (satn2o(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory satn2o'
      satn2o(:,:) = 0.0

      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable keqb ...'
      write(io_stdo_bgc,*)'First dimension    : ',11
      write(io_stdo_bgc,*)'Second dimension   : ',kpie
      write(io_stdo_bgc,*)'Third dimension    : ',kpje
      endif

      ALLOCATE (keqb(11,kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory keqb'
      keqb(:,:,:) = 0.0

      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable satoxy ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',kpke
      endif

      ALLOCATE (satoxy(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory satoxy'
      satoxy(:,:,:) = 0.0

      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable atm ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',natm
      endif

      ALLOCATE (atm(kpie,kpje,natm),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory atm'
      atm(:,:,:) = 0.0


      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable atmflx ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',natm
      endif

      ALLOCATE (atmflx(kpie,kpje,natm),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory atmflx'
      atmflx(:,:,:) = 0.0

      ! Allocate field to hold N-deposition fluxes per timestep for inventory calculations and output
      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable ndepflx ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (ndepflx(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory ndepflx'
      ndepflx(:,:) = 0.0

      ! Allocate field to hold OA alkalinity fluxes per timestep for inventory calculations and output
      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable oalkflx ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (oalkflx(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory oalkflx'
      oalkflx(:,:) = 0.0

      ! Allocate field to hold riverine fluxes per timestep for inventory calculations
      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable rivinflx ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third  dimension   : ',nriv
      endif
  
      ALLOCATE(rivinflx(kpie,kpje,nriv),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory rivinflx'
      rivinflx(:,:,:) = 0.0

      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable pco2d ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (pco2d(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory pco2d'
      pco2d(:,:) = 0.0

      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable pco2m ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (pco2m(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory pco2m'
      pco2m(:,:) = 0.0

      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable kwco2d ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (kwco2d(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory kwco2d'
      kwco2d(:,:) = 0.0

      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable kwco2sol ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif
      
      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable co2sold ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (co2sold(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory co2sold'
      co2sold(:,:) = 0.0
      
      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable co2solm ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (co2solm(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory co2solm'
      co2solm(:,:) = 0.0

      ALLOCATE (kwco2sol(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory co2fxd,co2fxu'
      kwco2sol(:,:) = 0.0


      if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable co2fxd, co2fxu ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif

      ALLOCATE (co2fxd(kpie,kpje),stat=errstat)
      ALLOCATE (co2fxu(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory co2fxd,co2fxu'
      co2fxd(:,:) = 0.0
      co2fxu(:,:) = 0.0

      if (use_cisonew) then
         if (mnproc.eq.1) then
         write(io_stdo_bgc,*)'Memory allocation for variable co213fxd,..., co214fxu ...'
         write(io_stdo_bgc,*)'First dimension    : ',kpie
         write(io_stdo_bgc,*)'Second dimension   : ',kpje
         endif

         ALLOCATE (co213fxd(kpie,kpje),stat=errstat)
         ALLOCATE (co213fxu(kpie,kpje),stat=errstat)
         ALLOCATE (co214fxd(kpie,kpje),stat=errstat)
         ALLOCATE (co214fxu(kpie,kpje),stat=errstat)
         if(errstat.ne.0) stop 'not enough memory co213fxd,..., co214fxu'
         co213fxd(:,:) = 0.0
         co213fxu(:,:) = 0.0
         co214fxd(:,:) = 0.0
         co214fxu(:,:) = 0.0
      endif

!******************************************************************************
      END SUBROUTINE ALLOC_MEM_CARBCH

      END MODULE mo_carbch
