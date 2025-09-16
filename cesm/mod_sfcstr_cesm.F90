! ------------------------------------------------------------------------------
! Copyright (C) 2015-2025 Mats Bentsen, Mehmet Ilicak, Mariana Vertenstein
!
! This file is part of BLOM.
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
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

module mod_sfcstr_cesm

  use mod_xc
  use mod_forcing,   only: ztx, mty, taux, tauy
  use mod_checksum,  only: csdiag, chksum

  implicit none
  private

  public :: sfcstr_cesm

contains

  subroutine sfcstr_cesm()

    ! ------------------------------------------------------------------
    ! Compute the surface stress. To be used when coupled to CESM
    ! ------------------------------------------------------------------

    ! Local variables
    integer :: i,j,l

    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
          taux(i,j) = ztx(i,j)
        end do
      end do
      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          tauy(i,j) = mty(i,j)
        end do
      end do
    end do
    !$omp end parallel do

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'sfcstr_cesm:'
      end if
      call chksum(taux, 1, halo_uv, 'taux')
      call chksum(tauy, 1, halo_vv, 'tauy')
    end if

  end subroutine sfcstr_cesm

end module mod_sfcstr_cesm
