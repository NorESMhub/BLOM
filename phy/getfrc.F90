! ------------------------------------------------------------------------------
! Copyright (C) 2015-2020 Mats Bentsen
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

subroutine getfrc
! ---------------------------------------------------------------------------
! Get forcing.
! ---------------------------------------------------------------------------

   use mod_config, only: expcnf
   use mod_xc, only: lp, mnproc, xcstop
   use mod_cesm, only: getfrc_cesm
   use mod_ben02, only: getfrc_ben02clim, getfrc_ben02syn

   implicit none

   select case (trim(expcnf))
      case ('cesm')
         call getfrc_cesm
      case ('ben02clim')
         call getfrc_ben02clim
      case ('ben02syn')
         call getfrc_ben02syn
      case ('channel')
      case ('fuk95')
      case ('isomip1')
!        call getfrc_isomip1
      case ('isomip2')
!        call getfrc_isomip2
      case ('single_column')
         call getfrc_ben02clim  
      case default
         if (mnproc == 1) then
            write (lp,'(3a)') ' getfrc: expcnf = ', trim(expcnf), &
                              ' is unsupported!'
         endif
         call xcstop('(getfrc)')
                stop '(getfrc)'
   end select

end subroutine getfrc
