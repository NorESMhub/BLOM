! Copyright (C) 2020  S. Gao, I. Bethke, J. Tjiputra, J. Schwinger
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


module mo_apply_ndep

  !*************************************************************************************************
  ! Routine for applying the nitrogen deposition flux.
  !
  ! N-deposition is activated through a logical switch 'do_ndep' read from HAMOCC's bgcnml
  ! namelist. When coupled to NorESM, this is achieved by setting BLOM_N_DEPOSITION to
  ! TRUE in env_run.xml.
  !
  ! The routine n_deposition applies the nitrogen deposition flux to the top-most model layer.
  !
  ! S.Gao             *Gfi, Bergen*             2017-08-19
  !
  ! Modified:
  !  J. Tjiputra,      *Uni Research, Bergen*    2017-09-18
  !  -add 1 mol [H+], per mol [NO3] deposition, to alkalinity (minus 1 mol)
  !  J. Schwinger,     *NORCE climate, Bergen*   2022-05-18
  !  -seperate modules into one module that reads a specific data set, and this
  !   module that applies the n-deposition flux to the surface ocean
  !*************************************************************************************************

  implicit none
  private

  public :: apply_ndep ! Apply n-deposition to the top-most model layer.

contains

  subroutine apply_ndep(kpie,kpje,kpke,pddpo,omask,ndep)
    !***********************************************************************************************
    ! Apply n-deposition to the top-most model layer.
    !
    ! S. Gao               *Gfi, Bergen*    19.08.2017
    !
    ! Modified:
    !  Tjiputra (18.09.2017): add 1 mol [H+], per mol [NO3] deposition, to alkalinity (minus 1 mol)
    !***********************************************************************************************

    use mo_control_bgc, only: dtb,do_ndep,use_extNcycle
    use mo_carbch,      only: ocetra,ndepnoyflx,ndepnhxflx
    use mo_param1_bgc,  only: iano3,ialkali,inatalkali,nndep,idepnoy,ianh4,idepnhx
    use mo_control_bgc, only: use_natDIC

    ! Arguments
    integer, intent(in) :: kpie                  ! 1st dimension of model grid.
    integer, intent(in) :: kpje                  ! 2nd dimension of model grid.
    integer, intent(in) :: kpke                  ! 3rd (vertical) dimension of model grid.
    real,    intent(in) :: pddpo(kpie,kpje,kpke) ! size of grid cell (depth) [m].
    real,    intent(in) :: omask(kpie,kpje)      ! land/ocean mask (1=ocean)
    real,    intent(in) :: ndep(kpie,kpje,nndep) ! N-deposition field to apply

    ! local variables
    integer :: i,j

    ! ndepflx stores the applied n-deposition flux for inventory calculations and output
    ndepnoyflx(:,:)=0.0
    if (use_extNcycle) then
      ndepnhxflx(:,:)=0.0
    endif
    if (.not. do_ndep) return

    ! deposit N in topmost layer
    do j=1,kpje
      do i=1,kpie
        if (omask(i,j).gt.0.5) then
          ndepnoyflx(i,j) = ndep(i,j,idepnoy)*dtb/365.
          ocetra(i,j,1,iano3)=ocetra(i,j,1,iano3)+ndepnoyflx(i,j)/pddpo(i,j,1)
          ocetra(i,j,1,ialkali)=ocetra(i,j,1,ialkali)-ndepnoyflx(i,j)/pddpo(i,j,1)
          if (use_natDIC) then
            ocetra(i,j,1,inatalkali)=ocetra(i,j,1,inatalkali)-ndepnoyflx(i,j)/pddpo(i,j,1)
          endif
          if (use_extNcycle) then
            ndepnhxflx(i,j)       = ndep(i,j,idepnhx)*dtb/365.
            ocetra(i,j,1,ianh4)   = ocetra(i,j,1,ianh4)   + ndepnhxflx(i,j)/pddpo(i,j,1)
            ocetra(i,j,1,ialkali) = ocetra(i,j,1,ialkali) + ndepnhxflx(i,j)/pddpo(i,j,1)
          endif
        endif
      enddo
    enddo

  end subroutine apply_ndep
 
  !*************************************************************************************************
end module mo_apply_ndep
