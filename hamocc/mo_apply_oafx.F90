! Copyright (C) 2021  J. Schwinger
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


module mo_apply_oafx
!******************************************************************************
!
!   J.Schwinger             *NORCE Climate, Bergen*             2021-11-15
!
! Modified
! --------
!
! Purpose
! -------
!  -Routines for applying ocean alkalinization
!
!
! Description:
! ------------
!
!  -subroutine alkalinization
!     Apply alkalinization to the top-most model layer.
!
!
!******************************************************************************
  implicit none

  private
  public :: apply_oafx

!******************************************************************************
contains



subroutine apply_oafx(kpie,kpje,kpke,pddpo,omask,oafx)
!******************************************************************************
!
!     J. Schwinger            *NORCE Climate, Bergen*     2021-11-15
!
! Purpose
! -------
!  -apply alkalinization to the top-most model layer.
!
! Changes: 
! --------
!
!
! Parameter list:
! ---------------
!  *INTEGER*   *kpie*    - 1st dimension of model grid.
!  *INTEGER*   *kpje*    - 2nd dimension of model grid.
!  *REAL*      *pddpo*   - size of grid cell (depth) [m].
!  *REAL*      *omask*   - land/ocean mask (1=ocean)
!  *REAL*      *oafx*    - alkalinization field to apply [kmol m-2 yr-1]
!
!******************************************************************************
  use mo_control_bgc, only: dtb,do_oalk
  use mo_carbch,      only: ocetra
  use mo_param1_bgc,  only: ialkali,inatalkali

  implicit none

  integer, intent(in) :: kpie,kpje,kpke
  real,    intent(in) :: pddpo(kpie,kpje,kpke)
  real,    intent(in) :: omask(kpie,kpje)
  real,    intent(in) :: oafx(kpie,kpje)

  ! local variables 
  integer :: i,j

  if (.not. do_oalk) return 

  ! alkalinization in topmost layer 
  do j=1,kpje
  do i=1,kpie
    if (omask(i,j).gt.0.5) then
      ocetra(i,j,1,ialkali)=ocetra(i,j,1,ialkali)+oafx(i,j)*dtb/365./pddpo(i,j,1)
#ifdef natDIC
      ocetra(i,j,1,inatalkali)=ocetra(i,j,1,inatalkali)+oafx(i,j)*dtb/365./pddpo(i,j,1)
#endif
    endif
  enddo
  enddo

!******************************************************************************
end subroutine apply_oafx


!******************************************************************************
end module mo_apply_oafx
