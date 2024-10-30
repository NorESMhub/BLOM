module mo_shelfsea_restime
  implicit none
  private

  public :: shelfsea_residence_time

contains

  subroutine shelfsea_residence_time(kpie,kpje,kpke,pddpo,shelfmask,omask)
    use mo_vgrid,       only : dp_min
    use mo_carbch,      only : ocetra
    use mo_param1_bgc,  only : ishelfage
    use mo_control_bgc, only : dtb,io_stdo_bgc

    integer, intent(in) :: kpie,kpje,kpke
    real,    intent(in) :: pddpo(kpie,kpje,kpke)          ! size of grid cell (3rd dimension) [m].
    logical, intent(in) :: shelfmask(kpie,kpje)                      ! True: shelf, False: else
    real,    intent(in) :: omask(kpie,kpje)                          ! land-ocean mask

    integer :: i,j,k

    !$OMP PARALLEL DO PRIVATE (i,j,k)
    do k=1,kpke
      do j=1,kpje
        do i=1,kpie
          if (pddpo(i,j,k)>dp_min .and. omask(i,j)>0.5) then
            ! Note that in Liu et al. 2019, min function is written,
            ! but a gradual decrease in residence time off the shelf should requir max function
            ocetra(i,j,k,ishelfage) = merge(       ocetra(i,j,k,ishelfage) + dtb,                  &
                                            max(0.,ocetra(i,j,k,ishelfage) - dtb), shelfmask(i,j) )
          endif
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine shelfsea_residence_time

end module mo_shelfsea_restime
