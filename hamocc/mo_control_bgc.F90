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
      INTEGER, save :: io_stdo_bgc        !  standard out.

! File containing namelists
      CHARACTER(LEN=:), ALLOCATABLE, PROTECTED :: bgc_namelist

! Control variables
      REAL,    save :: dtbgc              !  time step length [sec].
      REAL,    save :: dtb                !  time step length [days].
      INTEGER, save :: ndtdaybgc          !  time steps per day.

      INTEGER, save :: ldtbgc             !  time step number from bgc restart file
      INTEGER, save :: ldtrunbgc          !  actual time steps of run.

      INTEGER, save :: sedspin_yr_s = -1
      INTEGER, save :: sedspin_yr_e = -1
      INTEGER, save :: sedspin_ncyc = -1

      REAL,    save :: rmasks = 0.0       !  value at wet cells in sediment.
      REAL,    save :: rmasko = 99999.00  !  value at wet cells in ocean.
      
      LOGICAL, save :: l_sed_por  = .false.  ! apply lon-lat-variable sediment porosity via input file
! Logical switches set via namelist
      LOGICAL, save :: do_ndep     =.true.   ! apply n-deposition
      LOGICAL, save :: do_rivinpt  =.true.   ! apply riverine input
      LOGICAL, save :: do_sedspinup=.false.  ! apply sediment spin-up
      LOGICAL, save :: do_oalk     =.false.  ! apply ocean alkalinization
      logical, save :: with_dmsph  =.false.  ! apply DMS with pH dependence

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
              allocate(character(len=len('ocn_in'//trim(inst_suffix))) ::       &
                   bgc_namelist)
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
