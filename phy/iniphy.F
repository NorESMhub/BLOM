! ------------------------------------------------------------------------------
! Copyright (C) 2015-2021 Mats Bentsen, Mehmet Ilicak, Aleksi Nummelin
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

      subroutine iniphy
c
c --- ------------------------------------------------------------------
c --- Initialize physical parameterizations
c --- ------------------------------------------------------------------
c
      use mod_config, only: expcnf
      use mod_xc, only: lp, mnproc, xcstop
      use mod_vcoord, only: vcoord_type_tag, cntiso_hybrid
      use mod_tidaldissip, only: read_tidaldissip
      use mod_difest, only: init_difest
c
      implicit none
c
      if     (expcnf.eq.'cesm'.or.
     .        expcnf.eq.'ben02clim'.or.expcnf.eq.'ben02syn') then
        call read_tidaldissip
      elseif (expcnf.eq.'channel') then
      elseif (expcnf.eq.'fuk95') then
      elseif (expcnf.eq.'single_column') then
      elseif (expcnf.eq.'isomip1') then
      elseif (expcnf.eq.'isomip2') then
      else
        if (mnproc.eq.1) then
          write (lp,'(3a)') ' expcnf=',trim(expcnf),' is unsupported!'
        endif
        call xcstop('(iniphy)')
               stop '(iniphy)'
      endif
c
      if (vcoord_type_tag == cntiso_hybrid) then
        call init_difest
      endif
c
      return
      end
