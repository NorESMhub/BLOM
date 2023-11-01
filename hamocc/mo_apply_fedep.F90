! Copyright (C) 2022  S. J. Schwinger
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


module mo_apply_fedep
!********************************************************************************
!
!  J. Schwinger,     *NORCE climate, Bergen*   2022-05-19
!
!
! Purpose
! -------
!  - Routines for applying iron deposition data
!
!
! Description:
! ------------
!  Public routines and variable of this module:
!
!  -subroutine apply_fedep
!    apply iron deposition to the ocean tracer field
!
!  This module replaces code previously found inside the ocprod-routine and 
!  encapsulates it in a module.
!
!  
! Changes: 
! --------
!
!********************************************************************************
implicit none

private
public :: apply_fedep

!********************************************************************************
contains


subroutine apply_fedep(kpie,kpje,kpke,pddpo,omask,dust)
!--------------------------------------------------------------------------------
!
! Purpose:
! --------
!  Apply dust deposition input to oceanic tracer fields
!
! Description:
! ------------
! 
!
! Arguments:
! ----------
!  *INTEGER*   *kpie*    - 1st dimension of model grid.
!  *INTEGER*   *kpje*    - 2nd dimension of model grid.
!  *INTEGER*   *kpke*    - 3rd (vertical) dimension of model grid.
!  *REAL*      *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!  *REAL*      *omask*   - ocean mask
!  *REAL*      *dust*    - dust deposition flux [kg/m2/month].
!
!--------------------------------------------------------------------------------
  use mo_control_bgc, only: dtb
  use mo_param1_bgc,  only: ifdust,iiron
  use mo_param_bgc,   only: perc_diron
  use mo_carbch,      only: ocetra

  implicit none

  integer,intent(in) :: kpie,kpje,kpke
  real,   intent(in) :: pddpo(kpie,kpje,kpke)
  real,   intent(in) :: omask(kpie,kpje)
  real,   intent(in) :: dust(kpie,kpje)

  ! local variables
  integer            :: i,j
  real               :: dustinp

! dust flux from the atmosphere to the surface layer; dust fields are
! monthly mean values (kg/m2/month - assume 30 days per month here)
! dissolved iron is a fixed fraction (typically 3.5%), and immediately released

!$OMP PARALLEL DO PRIVATE(i,dustinp)
  do j = 1,kpje
  do i = 1,kpie
    if(omask(i,j) > 0.5) then
      dustinp = dust(i,j) / 30. * dtb / pddpo(i,j,1)
      ocetra(i,j,1,ifdust) = ocetra(i,j,1,ifdust) + dustinp
      ocetra(i,j,1,iiron) = ocetra(i,j,1,iiron) + dustinp * perc_diron
    endif
  enddo
  enddo
!$OMP END PARALLEL DO


!--------------------------------------------------------------------------------
end subroutine apply_fedep


!********************************************************************************
end module mo_apply_fedep

