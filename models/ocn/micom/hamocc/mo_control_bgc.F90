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
      INTEGER :: io_stdi_bgc           !  standard in.

      INTEGER :: io_rsti_bgc           !  restart in. 
      INTEGER :: io_rsto_bgc           !  restart out. 

      INTEGER :: io_nml                !  namelist

! Control variables

      REAL    :: dtbgc            !  time step length [sec].
      REAL    :: dtb              !  time step length [days].
      INTEGER :: ndtdaybgc        !  time steps per day.

      INTEGER :: ldtbgc           !  time step number from bgc restart file
      INTEGER :: ldtrunbgc        !  actual time steps of run.


      INTEGER :: icyclibgc        !  switch for cyclicity.
      INTEGER :: ndtrunbgc        !  total no. of time steps of run.

      INTEGER :: bgcstartyear            !  year of ocean restart file
      INTEGER :: bgcstartmonth           !  month of ocean restart file
      INTEGER :: bgcstartday             !  day of ocean restart file

! MPIOM is using variable lyear already !
!      INTEGER :: ldtoce           !  time step number from bgc ocean file

      INTEGER :: isac             !  acceleration factor for sediment, read from namelist


      REAL    :: rmasks = 0.0       !  value at wet cells in sediment.
      REAL    :: rmasko = 99999.00  !  value at wet cells in ocean.


      INTEGER :: kchck = 0          !  switch for extended print control (0/1). 
      
      END MODULE mo_control_bgc
