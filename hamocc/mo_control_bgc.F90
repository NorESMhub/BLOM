! Copyright (C) 2002  S. Legutke
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

MODULE mo_control_bgc
  !***********************************************************************
  !
  !**** *MODULE mo_control_bgc* - control variables for bgc modules.
  !
  !     S.Legutke,        *MPI-MaD, HH*    28.02.02
  !
  !     Modified
  !     --------
  !     J.Schwinger,      *Uni Research, Bergen*   2018-04-12
  !     - removed unused variables
  !
  !     Purpose
  !     -------
  !     - declaration
  !
  !
  !**********************************************************************
  implicit none

  ! Logical unit number for I/O.
  INTEGER :: io_stdo_bgc        !  standard out.

  ! File containing namelists
  CHARACTER(LEN=:), ALLOCATABLE, PROTECTED :: bgc_namelist

  ! Control variables
  REAL    :: dtbgc                    !  time step length [sec].
  REAL    :: dtb                      !  time step length [days].
  INTEGER :: ndtdaybgc                !  time steps per day.

  INTEGER :: ldtbgc                   !  time step number from bgc restart file
  INTEGER :: ldtrunbgc                !  actual time steps of run.

  INTEGER :: sedspin_yr_s = -1
  INTEGER :: sedspin_yr_e = -1
  INTEGER :: sedspin_ncyc = -1

  REAL    :: rmasks = 0.0             !  value at wet cells in sediment.
  REAL    :: rmasko = 99999.00        !  value at wet cells in ocean.

  ! Logical switches set via namelist
  LOGICAL :: l_3Dvarsedpor = .false.  ! apply lon-lat-depth variable sediment porosity via input file
  LOGICAL :: do_ndep     =.true.      ! apply n-deposition
  LOGICAL :: do_rivinpt  =.true.      ! apply riverine input
  LOGICAL :: do_sedspinup=.false.     ! apply sediment spin-up
  LOGICAL :: do_oalk     =.false.     ! apply ocean alkalinization
  logical :: with_dmsph  =.false.     ! apply DMS with pH dependence

  logical :: use_BROMO              = .false.
  logical :: use_AGG                = .false.
  logical :: use_WLIN               = .false.
  logical :: use_natDIC             = .false.
  logical :: use_CFC                = .false.
  logical :: use_cisonew            = .false.
  logical :: use_PBGC_OCNP_TIMESTEP = .false.
  logical :: use_PBGC_CK_TIMESTEP   = .false.
  logical :: use_FB_BGC_OCE         = .false.
  logical :: use_BOXATM             = .false.
  logical :: use_sedbypass          = .false.
  character(len=64) :: ocn_co2_type

contains

  subroutine get_bgc_namelist
    !-------------------------------------------------------------------------
    ! Get filename for namelist file
    !-------------------------------------------------------------------------
    use mod_config, only: inst_suffix
    use mod_xc,     only: xchalt

    implicit none

    logical :: exists

    if (.not. allocated(bgc_namelist)) then
      inquire (file='ocn_in'//trim(inst_suffix), exist=exists)
      if (exists) then
        allocate(character(len=len('ocn_in'//trim(inst_suffix))) :: bgc_namelist)
        bgc_namelist = 'ocn_in'//trim(inst_suffix)
      else
        inquire (file='limits', exist=exists)
        if (exists) then
          allocate(character(len=len('limits')) :: bgc_namelist)
          bgc_namelist = 'limits'
        else
          call xchalt('cannot find limits file')
          stop 'cannot find limits file'
        endif
      endif
    endif
  end subroutine get_bgc_namelist

END MODULE mo_control_bgc
