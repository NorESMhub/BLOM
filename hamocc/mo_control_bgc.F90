      MODULE mo_control_bgc

!$Source: /scratch/local1/m212047/patrick/SRC_MPI/src_hamocc/RCS/mo_control_bgc.f90,v $\\
!$Revision: 1.1 $\\
!$Date: 2005/01/28 08:37:45 $\\

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
      INTEGER :: ndtrunbgc        !  total no. of time steps of run.

      INTEGER :: isac             !  acceleration factor for sediment, read from namelist


      REAL    :: rmasks = 0.0       !  value at wet cells in sediment.
      REAL    :: rmasko = 99999.00  !  value at wet cells in ocean.

! Logical switches
      LOGICAL, SAVE :: do_ndep=.true.    ! apply n-deposition   (set via namelist)
      LOGICAL, SAVE :: do_rivinpt=.true. ! apply riverine input (set via namelist)
      
      END MODULE mo_control_bgc
