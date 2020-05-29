! Copyright (C) 2020  J. Schwinger
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


module mo_intfcblom
!******************************************************************************
!
! MODULE mo_intfcblom - Variables for BLOM-iHAMOCC interface
!
!  J.Schwinger,        *NORCE Climate, Bergen*    2020-05-19
!
!  Modified
!  --------
!
!  Purpose
!  -------
!   Declaration and memory allocation related to the BLOM-iHAMOCC interface.
!   This includes 2-time-level copies of sediment and amospheric fields.
!
!  To do
!  -----
!   Once BLOM has transitioned to free source format, the interface routines
!   blom2hamocc and hamocc2blom can be incorporated into this module
!
!
!   *nphys*      *INTEGER*  - number of bgc timesteps per ocean timestep.
!   *bgc_dx*     *REAL*     - size of grid cell (longitudinal) [m].
!   *bgc_dx*     *REAL*     - size of grid cell (latitudinal) [m].
!   *bgc_dp*     *REAL*     - size of grid cell (depth) [m].
!   *bgc_rho*    *REAL*     - sea water density [kg/m^3].
!   *omask*      *REAL*     - land ocean mask.
!
! The following arrays are used to keep a two time-level copy of sediment
! and prognostic atmosphere fields. These arrays are copied back and forth 
! in blom2hamocc.F and hamocc2blom.F in the same manner as the tracer field.
! Also, they written/read to and from restart files:
!
!   *sedlay2*    *REAL*     - two time-level copy of sedlay
!   *powtra2*    *REAL*     - two time-level copy of powtra
!   *burial2*    *REAL*     - two time-level copy of burial
!   *atm2*       *REAL*     - two time-level copy of atm
!
!******************************************************************************
  implicit none

  integer, parameter      :: nphys=2

  real, allocatable, save :: bgc_dx(:,:),bgc_dy(:,:)
  real, allocatable, save :: bgc_dp(:,:,:)
  real, allocatable, save :: bgc_rho(:,:,:)
  real, allocatable, save :: omask(:,:)

  ! Two time-level copy of sediment fields
  real, allocatable, save :: sedlay2(:,:,:,:)
  real, allocatable, save :: powtra2(:,:,:,:)
  real, allocatable, save :: burial2(:,:,:,:)

  ! Two time level copy of prognostic atmosphere field 
  ! used if BOXATM is activated
  real, allocatable, save :: atm2(:,:,:,:)      

contains
!******************************************************************************



subroutine alloc_mem_intfcblom(kpie,kpje,kpke)
!******************************************************************************
!
! ALLOC_MEM_VGRID - Allocate variables in this module
!
!  J.Schwinger            *NORCE Climate, Bergen*       2020-05-19
!
!******************************************************************************
  use mod_xc,         only: mnproc
  use mo_control_bgc, only: io_stdo_bgc
  use mo_param1_bgc,  only: ks,nsedtra,npowtra,natm

  INTEGER, intent(in) :: kpie,kpje,kpke
  INTEGER             :: errstat
      

  IF (mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)' '
    WRITE(io_stdo_bgc,*)'***************************************************'
    WRITE(io_stdo_bgc,*)'Memory allocation for module mo_intfcblom :'
    WRITE(io_stdo_bgc,*)' '
  ENDIF


  IF (mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)'Memory allocation for variable bgc_dx, bgc_dy ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
  ENDIF

  ALLOCATE (bgc_dx(kpie,kpje),stat=errstat)
  ALLOCATE (bgc_dy(kpie,kpje),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory bgc_dx, bgc_dy'
  bgc_dx(:,:) = 0.0
  bgc_dy(:,:) = 0.0


  IF (mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)'Memory allocation for variable bgc_dp ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
  ENDIF

  ALLOCATE (bgc_dp(kpie,kpje,kpke),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory bgc_dp'
  bgc_dp(:,:,:) = 0.0


  IF (mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)'Memory allocation for variable bgc_rho ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
  ENDIF

  ALLOCATE (bgc_rho(kpie,kpje,kpke),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory bgc_dp'
  bgc_rho(:,:,:) = 0.0


  IF(mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)'Memory allocation for variable omask ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
  ENDIF

  ALLOCATE(omask(kpie,kpje),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory omask'
  omask(:,:) = 0.0
  
#ifndef sedbypass
  IF(mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)'Memory allocation for variable sedlay2 ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    WRITE(io_stdo_bgc,*)'Third dimension    : ',2*ks
    WRITE(io_stdo_bgc,*)'Fourth dimension   : ',nsedtra
  ENDIF

  ALLOCATE (sedlay2(kpie,kpje,2*ks,nsedtra),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory sedlay2'
  sedlay2(:,:,:,:) = 0.0


  IF(mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)'Memory allocation for variable powtra2 ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    WRITE(io_stdo_bgc,*)'Third dimension    : ',2*ks
    WRITE(io_stdo_bgc,*)'Fourth dimension   : ',npowtra
  ENDIF

  ALLOCATE (powtra2(kpie,kpje,2*ks,npowtra),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory powtra2'
  powtra2(:,:,:,:) = 0.0


  IF(mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)'Memory allocation for variable burial2 ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    WRITE(io_stdo_bgc,*)'Third dimension    : ',2
    WRITE(io_stdo_bgc,*)'Fourth dimension   : ',nsedtra
  ENDIF

  ALLOCATE (burial2(kpie,kpje,2,nsedtra),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory burial2'
  burial2(:,:,:,:) = 0.0
#endif

#if defined(BOXATM)
  IF (mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)'Memory allocation for variable atm ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    WRITE(io_stdo_bgc,*)'Third dimension    : ',2
    WRITE(io_stdo_bgc,*)'Fourth dimension   : ',natm
  ENDIF

  ALLOCATE (atm2(kpie,kpje,2,natm),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory atm2'
  atm2(:,:,:,:) = 0.0
#endif

!******************************************************************************
end subroutine alloc_mem_intfcblom


!******************************************************************************
end module mo_intfcblom
