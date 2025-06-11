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

  !*************************************************************************************************
  ! Routines for applying ocean alkalinization
  !
  ! J.Schwinger             *NORCE Climate, Bergen*             2021-11-15
  !*************************************************************************************************

  implicit none
  private

  public :: apply_oafx ! Apply alkalinization to the top-most model layer.

contains

  subroutine apply_oafx(kpie,kpje,kpke,pddpo,omask,oafx)

    !***********************************************************************************************
    ! Apply alkalinization to the top-most model layer.
    !
    ! J. Schwinger            *NORCE Climate, Bergen*     2021-11-15
    !***********************************************************************************************

    use mo_kind,        only: rp
    use mo_control_bgc, only: dtb,do_oalk
    use mo_param1_bgc,  only: ialkali
    use mo_carbch,      only: ocetra,oalkflx,OmegaA
    use mo_read_oafx,   only: thrh_omegaa

    ! Arguments
    integer, intent(in) :: kpie                     ! 1st dimension of model grid.
    integer, intent(in) :: kpje                     ! 2nd dimension of model grid.
    integer, intent(in) :: kpke                     ! 3rd (vertical) dimension of model grid.
    real(rp),intent(in) :: pddpo(kpie,kpje,kpke)    ! size of grid cell (depth) [m].
    real(rp),intent(in) :: omask(kpie,kpje)         ! land/ocean mask (1=ocean)
    real(rp),intent(in) :: oafx(kpie,kpje)          ! alkalinization field to apply [kmol m-2 yr-1]

    ! local variables
    integer :: i,j

    ! oalkflx stores the applied alaklinity flux for inventory calculations
    ! and output
    oalkflx(:,:)=0.0_rp

    if (.not. do_oalk) return

    ! alkalinization in topmost layer
    do j=1,kpje
      do i=1,kpie
        if (omask(i,j).gt.0.5_rp) then
          if (thrh_omegaa > 0.0_rp .and. OmegaA(i,j,1) > thrh_omegaa) cycle
          oalkflx(i,j) = oafx(i,j)*dtb/365._rp
          ocetra(i,j,1,ialkali)=ocetra(i,j,1,ialkali)+oalkflx(i,j)/pddpo(i,j,1)
        endif
      enddo
    enddo

  end subroutine apply_oafx

end module mo_apply_oafx
