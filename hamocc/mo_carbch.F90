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
!     - moved accumulation of all output fields to seperate subroutine,
!       new global fields for output defined here
!     - added OmegaA
!
!     Purpose
!     -------
!     - declaration and memory allocation
!
!
!**********************************************************************     
      implicit none
      
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: ocetra
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: atm      
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: atmflx
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: co3
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: co2star   
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: hi
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: OmegaA 
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: OmegaC 
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: keqb

      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: satoxy
      REAL, DIMENSION (:,:),     ALLOCATABLE :: satn2o
      REAL, DIMENSION (:,:),     ALLOCATABLE :: atdifv
      REAL, DIMENSION (:,:),     ALLOCATABLE :: suppco2
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: sedfluxo

      REAL, DIMENSION (:,:),     ALLOCATABLE :: pco2d
      REAL, DIMENSION (:,:),     ALLOCATABLE :: kwco2sol
      REAL, DIMENSION (:,:),     ALLOCATABLE :: co2fxd
      REAL, DIMENSION (:,:),     ALLOCATABLE :: co2fxu
#ifdef cisonew
      REAL, DIMENSION (:,:),     ALLOCATABLE :: co213fxd
      REAL, DIMENSION (:,:),     ALLOCATABLE :: co213fxu
      REAL, DIMENSION (:,:),     ALLOCATABLE :: co214fxd
      REAL, DIMENSION (:,:),     ALLOCATABLE :: co214fxu
#endif
      REAL :: dmspar(6)
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: pi_ph
#ifdef natDIC
      REAL                                   :: atm_co2_nat
      REAL, DIMENSION (:,:),     ALLOCATABLE :: natpco2d
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: nathi
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: natco3
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: natOmegaA
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: natOmegaC
#endif
      REAL :: atm_co2, atm_o2, atm_n2 
      REAL :: atm_c13, atm_c14  
#ifdef cisonew
      REAL :: c14_t_half, c14dec
#endif
#ifdef CFC
      REAL :: atm_cfc11_nh,atm_cfc11_sh
      REAL :: atm_cfc12_nh,atm_cfc12_sh
      REAL :: atm_sf6_nh,atm_sf6_sh
#endif

      CONTAINS

      SUBROUTINE ALLOC_MEM_CARBCH(kpie,kpje,kpke)
!******************************************************************************
! ALLOC_MEM_CARBCH - Allocate variables in this module
!******************************************************************************
      use mod_xc,         only: mnproc
      use mo_control_bgc, only: io_stdo_bgc
      use mo_param1_bgc,  only: nocetra,npowtra,natm

      INTEGER, intent(in) :: kpie,kpje,kpke
      INTEGER             :: errstat
      

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'***************************************************'
      WRITE(io_stdo_bgc,*)'Memory allocation for carbon chemistry module :'
      WRITE(io_stdo_bgc,*)' '
      ENDIF


      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable ocetra ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
      WRITE(io_stdo_bgc,*)'Forth dimension    : ',nocetra
      ENDIF

      ALLOCATE (ocetra(kpie,kpje,kpke,nocetra),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory ocetra'
      ocetra(:,:,:,:) = 0.0

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable hi ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
      ENDIF

      ALLOCATE (hi(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory hi'
      hi(:,:,:) = 0.0

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable co3 ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
      ENDIF

      ALLOCATE (co3(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory co3'
      co3(:,:,:) = 0.0

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable co2star ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
      ENDIF

      ALLOCATE (co2star(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory co2star'
      co2star(:,:,:) = 0.0


      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable OmegaA, OmegaC ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
      ENDIF

      ALLOCATE (OmegaA(kpie,kpje,kpke),stat=errstat)
      ALLOCATE (OmegaC(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory OmegaA, OmegaC'
      OmegaA(:,:,:) = 0.0
      OmegaC(:,:,:) = 0.0


#ifdef DMSPH
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable pi_ph ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
      ENDIF

      ALLOCATE (pi_ph(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory pi_ph'
      pi_ph(:,:,:) = 0.0
#endif
#ifdef natDIC
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable natpco2d ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      ENDIF

      ALLOCATE (natpco2d(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory natpco2d'
      natpco2d(:,:) = 0.0

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable nathi ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
      ENDIF

      ALLOCATE (nathi(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory nathi'
      nathi(:,:,:) = 0.0

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable natco3 ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
      ENDIF

      ALLOCATE (natco3(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory natco3'
      natco3(:,:,:) = 0.0

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable natOmegaA, natOmegaC ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
      ENDIF

      ALLOCATE (natOmegaA(kpie,kpje,kpke),stat=errstat)
      ALLOCATE (natOmegaC(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory natOmegaA, natOmegaC'
      natOmegaA(:,:,:) = 0.0
      natOmegaC(:,:,:) = 0.0
#endif

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable sedfluxo ..'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',npowtra
      ENDIF

      ALLOCATE (sedfluxo(kpie,kpje,npowtra),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory sedfluxo'
      sedfluxo(:,:,:) = 0.0

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable satn2o ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      ENDIF

      ALLOCATE (satn2o(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory satn2o'
      satn2o(:,:) = 0.0

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable keqb ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',11
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpie
      WRITE(io_stdo_bgc,*)'Third dimension    : ',kpje
      ENDIF

      ALLOCATE (keqb(11,kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory keqb'
      keqb(:,:,:) = 0.0

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable satoxy ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
      ENDIF

      ALLOCATE (satoxy(kpie,kpje,kpke),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory satoxy'
      satoxy(:,:,:) = 0.0

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable atm ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',natm
      ENDIF

      ALLOCATE (atm(kpie,kpje,natm),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory atm'
      atm(:,:,:) = 0.0


      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable atmflx ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',natm
      ENDIF

      ALLOCATE (atmflx(kpie,kpje,natm),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory atmflx'
      atmflx(:,:,:) = 0.0


      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable pco2d ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      ENDIF

      ALLOCATE (pco2d(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory pco2d'
      pco2d(:,:) = 0.0


      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable kwco2sol ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      ENDIF

      ALLOCATE (kwco2sol(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory co2fxd,co2fxu'
      kwco2sol(:,:) = 0.0


      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable co2fxd, co2fxu ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      ENDIF

      ALLOCATE (co2fxd(kpie,kpje),stat=errstat)
      ALLOCATE (co2fxu(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory co2fxd,co2fxu'
      co2fxd(:,:) = 0.0
      co2fxu(:,:) = 0.0

#ifdef cisonew
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable co213fxd,..., co214fxu ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      ENDIF

      ALLOCATE (co213fxd(kpie,kpje),stat=errstat)
      ALLOCATE (co213fxu(kpie,kpje),stat=errstat)
      ALLOCATE (co214fxd(kpie,kpje),stat=errstat)
      ALLOCATE (co214fxu(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory co213fxd,..., co214fxu'
      co213fxd(:,:) = 0.0
      co213fxu(:,:) = 0.0
      co214fxd(:,:) = 0.0
      co214fxu(:,:) = 0.0
#endif

!******************************************************************************
      END SUBROUTINE ALLOC_MEM_CARBCH

      END MODULE mo_carbch
