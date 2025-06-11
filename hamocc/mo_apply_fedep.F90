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

  !************************************************************************************************
  !  Routines for applying total dust and iron deposition data
  !  This module replaces code previously found inside the ocprod-routine and
  !  encapsulates it in a module.
  !
  !  J. Schwinger,     *NORCE climate, Bergen*   2022-05-19
  !************************************************************************************************

  implicit none
  private

  public :: apply_fedep ! apply dust/iron deposition to the ocean tracer field

contains

  subroutine apply_fedep(kpie,kpje,kpke,pddpo,omask,dust)
    !--------------------------------------------------------------------------------
    !  Apply dust/iron deposition input to oceanic tracer fields
    !--------------------------------------------------------------------------------

    use mo_control_bgc, only: dtb
    use mo_kind,        only: rp
    use mo_param1_bgc,  only: ifdust,iiron,itdust,isfe,ndust
    use mo_param_bgc,   only: sec_per_day
    use mo_carbch,      only: ocetra,dustflx

    integer,intent(in) :: kpie                      ! 1st dimension of model grid.
    integer,intent(in) :: kpje                      ! 2nd dimension of model grid.
    integer,intent(in) :: kpke                      ! 3rd (vertical) dimension of model grid.
    real,   intent(in) :: pddpo(kpie,kpje,kpke)     ! size of scalar grid cell (3rd dimension) [m].
    real,   intent(in) :: omask(kpie,kpje)          ! ocean mask
    real,   intent(in) :: dust(kpie,kpje,ndust)     ! dust deposition flux [kg dust/m2/s] and [kmol sFe/m2/s].

    ! local variables
    integer :: i,j
    real    :: dustinp

    ! Total dust and soluble iron fluxes from the atmosphere to the surface layer; input
    ! fields are in kg/m2/s for total dust and in kmol fe/m2/s.
    !$OMP PARALLEL DO PRIVATE(i,dustinp)
    do j = 1,kpje
      do i = 1,kpie
        if(omask(i,j) > 0.5_rp) then
          dustflx(i,j,itdust)  = dust(i,j,itdust) * sec_per_day * dtb
          dustflx(i,j,isfe)    = dust(i,j,isfe)   * sec_per_day * dtb
          ocetra(i,j,1,ifdust) = ocetra(i,j,1,ifdust) + dustflx(i,j,itdust) / pddpo(i,j,1)
          ocetra(i,j,1,iiron)  = ocetra(i,j,1,iiron)  + dustflx(i,j,isfe)   / pddpo(i,j,1)
        endif
      enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine apply_fedep

end module mo_apply_fedep
