! ------------------------------------------------------------------------------
! Copyright (C) 2010-2020 Mats Bentsen

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

module mod_idlage

  use mod_xc
  use mod_time,    only: nday_in_year, delt1
  use mod_tracers, only: itriag, trc

  implicit none
  private

  public :: idlage_init
  public :: idlage_step

contains

  subroutine idlage_init

    ! ------------------------------------------------------------------
    ! initialization of ideal age tracer
    ! ------------------------------------------------------------------

    integer :: i,j,k,l

    !$omp parallel do private(k,l,i)
    do j = 1-nbdy,jj+nbdy
      do k = 1,kk
        do l = 1,isp(j)
          do i = ifp(j,l),ilp(j,l)
            trc(i,j,k   ,itriag) = 0.
            trc(i,j,k+kk,itriag) = 0.
          end do
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine idlage_init

  subroutine idlage_step(m,n,mm,nn,k1m,k1n)

    ! ------------------------------------------------------------------
    ! update ideal age tracer
    ! ------------------------------------------------------------------

    ! Arguments
    integer, intent(in) :: m,n,mm,nn,k1m,k1n

    ! Local variables
    integer :: i,j,k,l,kn
    real :: q

    ! set ideal age to zero in the surface layer and increment the age
    ! in the subsurface layers
    !$omp parallel do private(l,i)
    do j = 1-nbdy,jj+nbdy
      do l = 1,isp(j)
        do i = ifp(j,l),ilp(j,l)
          trc(i,j,k1n,itriag) = 0.
        end do
      end do
    end do
    !$omp end parallel do

    q = delt1/(86400.*nday_in_year)

    !$omp parallel do private(k,kn,l,i)
    do j = 1-nbdy,jj+nbdy
      do k = 2,kk
        kn = k+nn
        do l = 1,isp(j)
          do i = ifp(j,l),ilp(j,l)
            trc(i,j,kn,itriag) = trc(i,j,kn,itriag)+q
          end do
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine idlage_step

end module mod_idlage
