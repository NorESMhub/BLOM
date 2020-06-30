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

      INTEGER :: io_stdo_bgc           !  standard out.
      INTEGER :: io_nml                !  namelist

! Control variables

      REAL    :: dtbgc            !  time step length [sec].
      REAL    :: dtb              !  time step length [days].
      INTEGER :: ndtdaybgc        !  time steps per day.

      INTEGER :: ldtbgc           !  time step number from bgc restart file
      INTEGER :: ldtrunbgc        !  actual time steps of run.

      INTEGER :: isac             !  acceleration factor for sediment, read from namelist


      REAL    :: rmasks = 0.0       !  value at wet cells in sediment.
      REAL    :: rmasko = 99999.00  !  value at wet cells in ocean.

! Logical switches
      LOGICAL, SAVE :: do_ndep=.true.    ! apply n-deposition   (set via namelist)
      LOGICAL, SAVE :: do_rivinpt=.true. ! apply riverine input (set via namelist)
      
      END MODULE mo_control_bgc
