! ------------------------------------------------------------------------------
! Copyright (C) 2007-2015 Mats Bentsen
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

      subroutine restart_trcwt(rstfnm_ocn)
c
c --- ------------------------------------------------------------------
c --- Write tracer state to restart files
c --- ------------------------------------------------------------------
c
      use mod_config, only: expcnf, runid, inst_suffix
      use mod_xc
c
      implicit none
c
      character rstfnm_ocn*(*)
c
      character(len=256) :: rstfnm_ocntrc
#ifdef HAMOCC
      character(len=256) :: rstfnm_hamocc
#endif
c
c --- ------------------------------------------------------------------
c --- Generate file name
c --- ------------------------------------------------------------------
c
      if (expcnf.eq.'cesm') then
         rstfnm_ocntrc = trim(runid)//
     .        '.blom.rtrc.'//rstfnm_ocn(len_trim(runid)+10:)
#ifdef HAMOCC
         rstfnm_hamocc = trim(runid)//
     .        '.blom.rbgc.'//rstfnm_ocn(len_trim(runid)+10:)
#endif
      else
         rstfnm_ocntrc = trim(runid)//
     .        '_resttrc_'//rstfnm_ocn(len_trim(runid)+10:)
#ifdef HAMOCC
         rstfnm_hamocc = trim(runid)//
     .        '_restbgc_'//rstfnm_ocn(len_trim(runid)+10:)
#endif
      endif
c
c --- ------------------------------------------------------------------
c --- Write restart files
c --- ------------------------------------------------------------------
c
      call restart_ocntrcwt(rstfnm_ocntrc)
c
#ifdef HAMOCC
      call restart_hamoccwt(rstfnm_hamocc)
#endif
c
      return
      end
