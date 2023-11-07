! Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger, M. Bentsen
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

module mo_restart_hamoccwt

  implicit none
  private

  public :: restart_hamoccwt

contains

  subroutine restart_hamoccwt(rstfnm_ocn)
    !
    ! write restart for hamocc
    !
    use mod_time,      only: date,nstep
    use mod_xc,        only: idm,jdm,kdm
    use mod_tracers,   only: ntrbgc,ntr,itrbgc,trc
    use mo_intfcblom,  only: omask
    use mo_aufw_bgc,   only: aufw_bgc

    ! arguments
    character(len=*) :: rstfnm_ocn

    call aufw_bgc(idm,jdm,kdm,ntr,ntrbgc,itrbgc,trc,           &
                  date%year,date%month,date%day,nstep,omask,   &
                  rstfnm_ocn)

  end subroutine restart_hamoccwt

end module mo_restart_hamoccwt
