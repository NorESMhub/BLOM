! Copyright (C) 2021-2024  j. maerz
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


module mo_shelfsea_restime

  !*************************************************************************************************
  ! Routine to calculate the residence time of shel-sea water according to:
  !
  !     Liu et al. 2019: Simulating Water Residence Time in the Coastal Ocean: A Global Perspective
  !
  ! j.maerz *UiB, Bergen* 2024-10-31
  !
  !*************************************************************************************************

  implicit none
  private

  public :: shelfsea_residence_time

contains

  subroutine shelfsea_residence_time(kpie,kpje,kpke,pddpo,shelfmask,omask)

    use mo_vgrid,       only : dp_min
    use mo_carbch,      only : ocetra
    use mo_param1_bgc,  only : ishelfage
    use mo_control_bgc, only : dtb,io_stdo_bgc

    ! Arguments
    integer, intent(in) :: kpie                   ! 1st dimension of model grid
    integer, intent(in) :: kpje                   ! 2nd dimension of model grid
    integer, intent(in) :: kpke                   ! 3rd dimension of model grid
    real,    intent(in) :: pddpo(kpie,kpje,kpke)  ! size of grid cell (3rd dimension) [m].
    logical, intent(in) :: shelfmask(kpie,kpje)   ! shelf-sea mask True: shelf, False: else
    real,    intent(in) :: omask(kpie,kpje)       ! land-ocean mask

    ! Local variables
    integer :: i,j,k

    !$OMP PARALLEL DO PRIVATE (i,j,k)
    do k=1,kpke
      do j=1,kpje
        do i=1,kpie
          if (pddpo(i,j,k)>dp_min .and. omask(i,j)>0.5) then
            ! Note that in Liu et al. 2019, min function is written,
            ! but a gradual decrease in residence time off the shelf should require max function
            ! to result in zero values in open ocean regions (and not negative values)
            ocetra(i,j,k,ishelfage) = merge(       ocetra(i,j,k,ishelfage) + dtb,                  &
                                            max(0.,ocetra(i,j,k,ishelfage) - dtb), shelfmask(i,j) )
          endif
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine shelfsea_residence_time

end module mo_shelfsea_restime
