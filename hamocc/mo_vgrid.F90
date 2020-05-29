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


module mo_vgrid
!******************************************************************************
!
! MODULE mo_vgrid - Variables and routines related to vertical grid 
!                   structure
!
!  J.Schwinger,        *NORCE Climate, Bergen*    2020-05-19
!
!  Modified
!  --------
!
!  Purpose
!  -------
!   Declaration, memory allocation, and routines related to the
!   vertical grid structure. These have to be recalculated every
!   time step when iHAMOCC is coupled to BLOM.
!
!   *kbo*         *INTEGER*  - number of wet cells in column.
!   *kwrbioz*     *INTEGER*  - last k-index of euphotic zone.
!   *kxxxx*       *INTEGER*  - k-index of gridbox comprising xxxx m depth.
!   *bolay*       *REAL*     - height of bottom cell.
!   *ptiestu*     *REAL*     - depth of layer centres.
!   *ptiestw*     *REAL*     - depth of layer interfaces.
!
!******************************************************************************
  implicit none

  INTEGER, PARAMETER :: kmle   = 2        ! k-end index for layers that 
                                          ! represent the mixed layer in BLOM
  REAL,    PARAMETER :: dp_ez  = 100.0    ! depth of euphotic zone
  REAL,    PARAMETER :: dp_min = 1.0E-12  ! min layer thickness layers thinner 
                                          ! than this are ignored by HAMOCC
  REAL,    PARAMETER :: dp_min_sink = 1.0 ! min layer thickness for sinking (layers thinner than 
                                          ! this are ignored and set to the concentration of the 
                                          ! layer above). Note that the bottom layer index kbo(i,j)
                                          ! is defined as the lowermost layer thicker than dp_min_sink.

  INTEGER, DIMENSION(:,:),   ALLOCATABLE :: kbo
  INTEGER, DIMENSION(:,:),   ALLOCATABLE :: kwrbioz
  INTEGER, DIMENSION(:,:),   ALLOCATABLE :: k0100,k0500,k1000,k2000,k4000
  REAL,    DIMENSION(:,:),   ALLOCATABLE :: bolay
  REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ptiestu
  REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ptiestw

contains
!******************************************************************************



subroutine set_vgrid(kpie,kpje,kpke,pddpo)
!******************************************************************************
!
! SET_VGRID - Calculate variables related to the vertical grid structure. This
!             routine replaces calc_idepth and calc_bot.
!
!  J.Schwinger            *NORCE Climate, Bergen*       2020-05-19
!
!  Purpose
!  -------
!  -calculate depth of layer interfaces and centres based on layer thickness
!  -find lowest mass containing layer in the euphotic zone
!  -find k-index of 100,500,1000,2000, and 4000 m-surfaces
!
!  Parameter list:
!  ---------------
!     *INTEGER*   *kpie*    - 1st dimension of model grid.
!     *INTEGER*   *kpje*    - 2nd dimension of model grid.
!     *INTEGER*   *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*      *pddpo*   - size of grid cell (3rd dimension) [m].
!
!******************************************************************************
  INTEGER, intent(in) :: kpie,kpje,kpke
  REAL,    intent(in) :: pddpo(kpie,kpje,kpke)

  INTEGER             :: i,j,k


  ! --- set depth of surface interface to zero
  ptiestw(:,:,1)=0.
  ! --- depth of layer kpke+1 centre
  ptiestu(:,:,kpke+1)=9000.

!$OMP PARALLEL DO PRIVATE
  do k=1,kpke
    do j=1,kpje
    do i=1,kpie

        ! --- depth of layer interfaces
        ptiestw(i,j,k+1)=ptiestw(i,j,k)+pddpo(i,j,k)
        ! --- depth of layer centres
        ptiestu(i,j,k)=ptiestw(i,j,k)+0.5*pddpo(i,j,k)

    enddo
    enddo
  enddo
!$OMP END PARALLEL DO


  kbo(:,:)  =1
  bolay(:,:)=0.0

!$OMP PARALLEL DO
  DO j=1,kpje
  DO i=1,kpie

    DO k=kpke,1,-1
      IF(pddpo(i,j,k).GT.dp_min_sink) THEN
        bolay(i,j)=pddpo(i,j,k)
        kbo(i,j)=k
        exit
      ENDIF
    ENDDO

  ENDDO
  ENDDO
!$OMP END PARALLEL DO


!$OMP PARALLEL DO
  DO j=1,kpje
  DO i=1,kpie

    kwrbioz(i,j)=1
    DO k=2,kpke
      IF(pddpo(i,j,k) .gt. dp_min .and. ptiestw(i,j,k) .lt. dp_ez ) THEN
        kwrbioz(i,j)=k
       ENDIF
    END DO

  END DO
  END DO
!$OMP END PARALLEL DO


  k0100(:,:)=0
  k0500(:,:)=0
  k1000(:,:)=0
  k2000(:,:)=0
  k4000(:,:)=0

!$OMP PARALLEL DO
  DO j=1,kpje
  DO i=1,kpie

    DO k=2,kpke
      IF(pddpo(i,j,k) .gt. dp_min .and. ptiestw(i,j,k+1) .gt. 100.0 ) THEN
        k0100(i,j)=k
        exit
      ENDIF
    END DO

    DO k=2,kpke
      IF(pddpo(i,j,k) .gt. dp_min .and. ptiestw(i,j,k+1) .gt. 500.0 ) THEN
        k0500(i,j)=k
        exit
      ENDIF
    END DO

    DO k=2,kpke
      IF(pddpo(i,j,k) .gt. dp_min .and. ptiestw(i,j,k+1) .gt. 1000.0 ) THEN
        k1000(i,j)=k
        exit
      ENDIF
    END DO

    DO k=2,kpke
      IF(pddpo(i,j,k) .gt. dp_min .and. ptiestw(i,j,k+1) .gt. 2000.0 ) THEN
        k2000(i,j)=k
        exit
      ENDIF
    END DO

    DO k=2,kpke
      IF(pddpo(i,j,k) .gt. dp_min .and. ptiestw(i,j,k+1) .gt. 4000.0 ) THEN
        k4000(i,j)=k
        exit
      ENDIF
    END DO

  END DO
  END DO
!$OMP END PARALLEL DO

  RETURN

!******************************************************************************
end subroutine set_vgrid



subroutine alloc_mem_vgrid(kpie,kpje,kpke)
!******************************************************************************
!
! ALLOC_MEM_VGRID - Allocate variables in this module
!
!  J.Schwinger            *NORCE Climate, Bergen*       2020-05-19
!
!******************************************************************************
  use mod_xc,         only: mnproc
  use mo_control_bgc, only: io_stdo_bgc

  INTEGER, intent(in) :: kpie,kpje,kpke
  INTEGER             :: errstat

      
  IF (mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)' '
    WRITE(io_stdo_bgc,*)'***************************************************'
    WRITE(io_stdo_bgc,*)'Memory allocation for module mo_vgrid :'
    WRITE(io_stdo_bgc,*)' '
  ENDIF


  IF (mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)'Memory allocation for variable ptiestu ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke+1
  ENDIF

  ALLOCATE (ptiestu(kpie,kpje,kpke+1),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory ptiestu'
  ptiestu(:,:,:) = 0.0


  IF (mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)'Memory allocation for variable ptiestw ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke+1
  ENDIF

  ALLOCATE (ptiestw(kpie,kpje,kpke+1),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory ptiestw'
  ptiestw(:,:,:) = 0.0


  IF(mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)'Memory allocation for variable kbo ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
  ENDIF

  ALLOCATE(kbo(kpie,kpje),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory kbo'
  kbo(:,:) = 0


  IF(mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)'Memory allocation for variable kwrbioz...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
  ENDIF

  ALLOCATE(kwrbioz(kpie,kpje),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory kwrbioz'
  kwrbioz(:,:) = 0


  IF(mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)'Memory allocation for variables k0100, k0500, k1000, k2000 ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
  ENDIF

  ALLOCATE(k0100(kpie,kpje),stat=errstat)
  ALLOCATE(k0500(kpie,kpje),stat=errstat)
  ALLOCATE(k1000(kpie,kpje),stat=errstat)
  ALLOCATE(k2000(kpie,kpje),stat=errstat)
  ALLOCATE(k4000(kpie,kpje),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory k0100, k0500, k1000, k2000'
  k0100(:,:) = 0
  k0500(:,:) = 0
  k1000(:,:) = 0
  k2000(:,:) = 0
  k4000(:,:) = 0


  IF(mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)'Memory allocation for variable bolay ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
  ENDIF

  ALLOCATE (bolay(kpie,kpje),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory bolay'
  bolay(:,:) = 0.0


!******************************************************************************
end subroutine alloc_mem_vgrid


!******************************************************************************
end module mo_vgrid
