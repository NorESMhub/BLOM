! ------------------------------------------------------------------------------
! Copyright (C) 2004-2022 Mats Bentsen, Mehmet Ilicak

! This file is part of BLOM.

! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.

! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.

! You should have received a copy of the GNU Lesser General Public License
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

module mod_sfcstr_ben02

  use mod_xc,        only: xctilr, ii, jj, isu, ilu, ifu, isv, ifv, ilv, iu, iv, lp, &
                           halo_ps, mnproc
  use mod_forcing,   only: ztx, mty, taux, tauy
  use mod_seaice,    only: ficem, hicem, tauxice, tauyice
  use mod_checksum,  only: csdiag, chksummsk

  implicit none
  private

  public :: sfcstr_ben02

contains

  subroutine sfcstr_ben02()

    ! --- ------------------------------------------------------------------
    ! --- compute the surface stress
    ! --- ------------------------------------------------------------------

    integer :: i,j,l
    real :: facice

    call xctilr(ficem, 1,1, 1,1, halo_ps)
    call xctilr(hicem, 1,1, 1,1, halo_ps)

    !$omp parallel do private(l,i,facice)
    do j = 1,jj
      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
          facice = (ficem(i,j)+ficem(i-1,j)) * min(2.,hicem(i,j)+hicem(i-1,j))*.25
          taux(i,j) = (ztx(i,j)*(1.-facice)+tauxice(i,j)*facice)
        end do
      end do
      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          facice = (ficem(i,j)+ficem(i,j-1)) * min(2.,hicem(i,j)+hicem(i,j-1))*.25
          tauy(i,j) = (mty(i,j)*(1.-facice)+tauyice(i,j)*facice)
        end do
      end do
    end do
    !$omp end parallel do

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'sfcstr:'
      end if
      call chksummsk(taux,iu,1,'taux')
      call chksummsk(tauy,iv,1,'tauy')
    end if

  end subroutine sfcstr_ben02

end module mod_sfcstr_ben02
