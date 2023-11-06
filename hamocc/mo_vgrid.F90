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


MODULE mo_vgrid
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
!  Description:
!  ------------
!  Public routines and variables of this MODULE:
!
!  -SUBROUTINE set_vgrid
!     Calculate variables related to the vertical grid structure.
!
!  -SUBROUTINE alloc_mem_vgrid
!     Allocate memory for vertical grid variables
!
!   *kbo*         *integer*  - number of wet cells in column.
!   *kwrbioz*     *integer*  - last k-index of euphotic zone.
!   *kxxxx*       *integer*  - k-index of gridbox comprising xxxx m depth.
!   *bolay*       *real*     - height of bottom cell.
!   *ptiestu*     *real*     - depth of layer centres.
!   *ptiestw*     *real*     - depth of layer interfaces.
!
!******************************************************************************
  IMPLICIT NONE

  integer, parameter :: kmle_static = 2   ! k-end index for layers that
                                          ! represent the mixed layer in BLOM.
                                          ! Default value used for isopycnic coordinates.
  real,    parameter :: dp_ez  = 100.0    ! depth of euphotic zone
  real,    parameter :: dp_min = 1.0E-12  ! min layer thickness layers thinner
                                          ! than this are ignored by HAMOCC
  real,    parameter :: dp_min_sink = 1.0 ! min layer thickness for sinking (layers thinner than 
                                          ! this are ignored and set to the concentration of the 
                                          ! layer above). Note that the bottom layer index kbo(i,j)
                                          ! is defined as the lowermost layer thicker than dp_min_sink.

  integer, dimension(:,:),   allocatable :: kmle
  integer, dimension(:,:),   allocatable :: kbo
  integer, dimension(:,:),   allocatable :: kwrbioz
  integer, dimension(:,:),   allocatable :: k0100,k0500,k1000,k2000,k4000
  real,    dimension(:,:),   allocatable :: bolay
  real,    dimension(:,:,:), allocatable :: ptiestu
  real,    dimension(:,:,:), allocatable :: ptiestw

CONTAINS
!******************************************************************************



SUBROUTINE set_vgrid(kpie,kpje,kpke,pddpo)
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
!     *integer*   *kpie*    - 1st dimension of model grid.
!     *integer*   *kpje*    - 2nd dimension of model grid.
!     *integer*   *kpke*    - 3rd (vertical) dimension of model grid.
!     *real*      *pddpo*   - size of grid cell (3rd dimension) [m].
!
!******************************************************************************
  integer, intent(in) :: kpie,kpje,kpke
  real,    intent(in) :: pddpo(kpie,kpje,kpke)

  integer             :: i,j,k


  ! --- set depth of surface interface to zero
  ptiestw(:,:,1)=0.
  ! --- depth of layer kpke+1 centre
  ptiestu(:,:,kpke+1)=9000.

!$OMP PARALLEL DO PRIVATE(j,i)
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

!$OMP PARALLEL DO PRIVATE(i,k)
  do j=1,kpje
  do i=1,kpie

    do k=kpke,1,-1
      if (pddpo(i,j,k).GT.dp_min_sink) then
        bolay(i,j)=pddpo(i,j,k)
        kbo(i,j)=k
        exit
      endif
    enddo

  enddo
  enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(i,k)
  do j=1,kpje
  do i=1,kpie

    kwrbioz(i,j)=1
    do k=2,kpke
      if (pddpo(i,j,k) .gt. dp_min .and. ptiestw(i,j,k) .lt. dp_ez ) then
        kwrbioz(i,j)=k
       endif
    END DO

  END DO
  END DO
!$OMP END PARALLEL DO


  k0100(:,:)=0
  k0500(:,:)=0
  k1000(:,:)=0
  k2000(:,:)=0
  k4000(:,:)=0

!$OMP PARALLEL DO PRIVATE(i,k)
  do j=1,kpje
  do i=1,kpie

    do k=2,kpke
      if (pddpo(i,j,k) .gt. dp_min .and. ptiestw(i,j,k+1) .gt. 100.0 ) then
        k0100(i,j)=k
        exit
      endif
    END DO

    do k=2,kpke
      if (pddpo(i,j,k) .gt. dp_min .and. ptiestw(i,j,k+1) .gt. 500.0 ) then
        k0500(i,j)=k
        exit
      endif
    END DO

    do k=2,kpke
      if (pddpo(i,j,k) .gt. dp_min .and. ptiestw(i,j,k+1) .gt. 1000.0 ) then
        k1000(i,j)=k
        exit
      endif
    END DO

    do k=2,kpke
      if (pddpo(i,j,k) .gt. dp_min .and. ptiestw(i,j,k+1) .gt. 2000.0 ) then
        k2000(i,j)=k
        exit
      endif
    END DO

    do k=2,kpke
      if (pddpo(i,j,k) .gt. dp_min .and. ptiestw(i,j,k+1) .gt. 4000.0 ) then
        k4000(i,j)=k
        exit
      endif
    END DO

  END DO
  END DO
!$OMP END PARALLEL DO

  return

!******************************************************************************
END SUBROUTINE set_vgrid



SUBROUTINE alloc_mem_vgrid(kpie,kpje,kpke)
!******************************************************************************
!
! ALLOC_MEM_VGRID - Allocate variables in this MODULE
!
!  J.Schwinger            *NORCE Climate, Bergen*       2020-05-19
!
!******************************************************************************
  use mod_xc,         only: mnproc
  use mo_control_bgc, only: io_stdo_bgc

  integer, intent(in) :: kpie,kpje,kpke
  integer             :: errstat

      
  if (mnproc.eq.1) then
    write(io_stdo_bgc,*)' '
    write(io_stdo_bgc,*)'***************************************************'
    write(io_stdo_bgc,*)'Memory allocation for MODULE mo_vgrid :'
    write(io_stdo_bgc,*)' '
  endif


  if (mnproc.eq.1) then
    write(io_stdo_bgc,*)'Memory allocation for variable ptiestu ...'
    write(io_stdo_bgc,*)'First dimension    : ',kpie
    write(io_stdo_bgc,*)'Second dimension   : ',kpje
    write(io_stdo_bgc,*)'Third dimension    : ',kpke+1
  endif

  ALLOCATE (ptiestu(kpie,kpje,kpke+1),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory ptiestu'
  ptiestu(:,:,:) = 0.0


  if (mnproc.eq.1) then
    write(io_stdo_bgc,*)'Memory allocation for variable ptiestw ...'
    write(io_stdo_bgc,*)'First dimension    : ',kpie
    write(io_stdo_bgc,*)'Second dimension   : ',kpje
    write(io_stdo_bgc,*)'Third dimension    : ',kpke+1
  endif

  ALLOCATE (ptiestw(kpie,kpje,kpke+1),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory ptiestw'
  ptiestw(:,:,:) = 0.0


  if (mnproc.eq.1) then
     write(io_stdo_bgc,*)'Memory allocation for variable kmle ...'
     write(io_stdo_bgc,*)'First dimension    : ',kpie
     write(io_stdo_bgc,*)'Second dimension   : ',kpje
  endif

  ALLOCATE(kmle(kpie,kpje),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory kmle'
  kmle(:,:) = kmle_static


  if (mnproc.eq.1) then
    write(io_stdo_bgc,*)'Memory allocation for variable kbo ...'
    write(io_stdo_bgc,*)'First dimension    : ',kpie
    write(io_stdo_bgc,*)'Second dimension   : ',kpje
  endif

  ALLOCATE(kbo(kpie,kpje),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory kbo'
  kbo(:,:) = 0


  if (mnproc.eq.1) then
    write(io_stdo_bgc,*)'Memory allocation for variable kwrbioz...'
    write(io_stdo_bgc,*)'First dimension    : ',kpie
    write(io_stdo_bgc,*)'Second dimension   : ',kpje
  endif

  ALLOCATE(kwrbioz(kpie,kpje),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory kwrbioz'
  kwrbioz(:,:) = 0


  if (mnproc.eq.1) then
    write(io_stdo_bgc,*)'Memory allocation for variables k0100, k0500, k1000, k2000 ...'
    write(io_stdo_bgc,*)'First dimension    : ',kpie
    write(io_stdo_bgc,*)'Second dimension   : ',kpje
  endif

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


  if (mnproc.eq.1) then
    write(io_stdo_bgc,*)'Memory allocation for variable bolay ...'
    write(io_stdo_bgc,*)'First dimension    : ',kpie
    write(io_stdo_bgc,*)'Second dimension   : ',kpje
  endif

  ALLOCATE (bolay(kpie,kpje),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory bolay'
  bolay(:,:) = 0.0


!******************************************************************************
END SUBROUTINE alloc_mem_vgrid


!******************************************************************************
END MODULE mo_vgrid
