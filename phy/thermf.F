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

      subroutine thermf(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Get surface forcing functions
c --- ------------------------------------------------------------------
c
      use mod_config, only: expcnf
      use mod_xc, only: lp, mnproc, xcstop
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
      select case (trim(expcnf))
        case ('cesm')
          call thermf_cesm(m,n,mm,nn,k1m,k1n)
        case ('ben02clim', 'ben02syn', 'single_column')
          call thermf_ben02(m,n,mm,nn,k1m,k1n)
        case ('channel')
          call thermf_channel(m,n,mm,nn,k1m,k1n)
        case ('fuk95')
        case ('isomip1')
c         call thermf_isomip1(m,n,mm,nn,k1m,k1n)
        case ('isomip2')
c         call thermf_isomip2(m,n,mm,nn,k1m,k1n)
        case default
          if (mnproc.eq.1) then
             write (lp,'(3a)') ' thermf: expcnf = ', trim(expcnf),
     .                         ' is unsupported!'
          endif
          call xcstop('(thermf)')
                 stop '(thermf)'
      end select
c
      return
      end
