! ------------------------------------------------------------------------------
! Copyright (C) 2004-2020 Mats Bentsen
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

subroutine sfcstr(m, n, mm, nn, k1m, k1n)
! ---------------------------------------------------------------------------
! Get surface stress.
! ---------------------------------------------------------------------------

    use mod_config, only: expcnf
    use mod_xc, only: lp, mnproc, xcstop

   implicit none

   integer, intent(in) :: m, n, mm, nn, k1m, k1n

   select case (trim(expcnf))
      case ('cesm')
         call sfcstr_cesm(m, n, mm, nn, k1m, k1n)
      case ('ben02clim', 'ben02syn', 'single_column')
         call sfcstr_ben02(m, n, mm, nn, k1m, k1n)
      case ('fuk95')
      case ('channel')
      case ('isomip1')
!        call sfcstr_isomip1(m, n, mm, nn, k1m, k1n)
      case ('isomip2')
!        call sfcstr_isomip2(m, n, mm, nn, k1m, k1n)
      case default
         if (mnproc == 1) then
            write (lp,'(3a)') ' sfcstr: expcnf = ', trim(expcnf), &
                              ' is unsupported!'
         endif
         call xcstop('(sfcstr)')
                stop '(sfcstr)'
   end select

end subroutine sfcstr
