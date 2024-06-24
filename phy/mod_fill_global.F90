! ------------------------------------------------------------------------------
! Copyright (C) 2004-2024 Mats Bentsen, Mariana Vertenstein
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

module mod_fill_global

  use mod_xc

  implicit none
  private

  public :: fill_global

contains

  subroutine fill_global(missing_value,fill_value,itype,field)

    ! ------------------------------------------------------------------
    ! Fill missing values by extrapolating values from neighbouring
    ! points. A value = fill_value will be modified.
    ! ------------------------------------------------------------------

    ! Arguments
    real,    intent(in) :: missing_value,fill_value
    integer, intent(in) :: itype
    real,    intent(inout), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: field

    ! Local variables
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tmp
    real :: vland_orig,sum
    integer :: i,j,mbdy,done,changed,num,il,jl

    vland_orig = vland
    vland = missing_value
    mbdy = 1
    done = 0
    changed = 1

    do while (changed == 1)

      if (mbdy == 1) then
        mbdy = nbdy
        call xctilr(field, 1,1, nbdy,nbdy, itype)
      else
        mbdy = mbdy-1
      end if

      done = 1
      changed = 0

      !$omp parallel do private(i,sum,num,jl,il) &
      !$omp reduction(min:done) reduction(max:changed)
      do j = 1-mbdy+1,jj+mbdy-1
        do i = 1-mbdy+1,ii+mbdy-1
          if (abs(field(i,j)) == abs(fill_value)) then
            done = min(0,done)
            sum = 0.
            num = 0
            do jl = j-1,j+1
              do il = i-1,i+1
                if (abs(field(il,jl)) /= abs(fill_value).and. &
                    abs(field(il,jl)) /= abs(missing_value)) then
                  sum = sum+field(il,jl)
                  num = num+1
                end if
              end do
            end do
            if (num > 0) then
              changed = max(1,changed)
              tmp(i,j) = sum/real(num)
            else
              tmp(i,j) = field(i,j)
            end if
          else
            tmp(i,j) = field(i,j)
          end if
        end do
      end do
      !$omp end parallel do

      !$omp parallel do private(i)
      do j = 1-mbdy+1,jj+mbdy-1
        do i = 1-mbdy+1,ii+mbdy-1
          field(i,j) = tmp(i,j)
        end do
      end do
      !$omp end parallel do

      call xcmin(done)
      call xcmax(changed)

    end do

    vland = vland_orig

    if (done == 0) then
      if (mnproc == 1) then
        write(lp,*) 'fill_global: filling failed!'
      end if
      return
    end if

  end subroutine fill_global

end module mod_fill_global
