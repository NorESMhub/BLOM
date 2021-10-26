! ------------------------------------------------------------------------------
! Copyright (C) 2007-2021 Mats Bentsen, Tomas Torsvik
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
      use mod_config, only: expcnf
      use mod_xc
c
      implicit none
c
      character rstfnm_ocn*(*)
c
      logical :: error
      character(len=256) :: rstfnm_ocntrc
#ifdef HAMOCC
      character(len=256) :: rstfnm_hamocc
#endif
c
c --- ------------------------------------------------------------------
c --- Generate file name
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
         if (expcnf.eq.'cesm') then
            call restart_getfile(
     .           rstfnm_ocn, 'rtrc', rstfnm_ocntrc, error)
            if (error) then
               write(lp,*) 'restart_trcwt: '//
     .              'could not generate rstfnm_ocntrc file!'
               call xcstop('(restat_trcwt)')
               stop '(restart_trcwt)'
            endif
#ifdef HAMOCC
            call restart_getfile(
     .           rstfnm_ocn, 'rbgc', rstfnm_hamocc, error)
            if (error) then
               write(lp,*) 'restart_trcwt: '//
     .              'could not generate rstfnm_hamocc file!'
               call xcstop('(restat_trcwt)')
               stop '(restart_trcwt)'
            endif
#endif
         else
            call restart_getfile(
     .           rstfnm_ocn, 'resttrc', rstfnm_ocntrc, error)
            if (error) then
               write(lp,*) 'restart_trcwt: '//
     .              'could not generate rstfnm_ocntrc file!'
               call xcstop('(restat_trcwt)')
               stop '(restart_trcwt)'
            endif
#ifdef HAMOCC
            call restart_getfile(
     .           rstfnm_ocn, 'restbgc', rstfnm_hamocc, error)
            if (error) then
               write(lp,*) 'restart_trcwt: '//
     .              'could not generate rstfnm_hamocc file!'
               call xcstop('(restat_trcwt)')
               stop '(restart_trcwt)'
            endif
#endif
         endif
      endif

      call xcbcst(rstfn_ocntrc)
#ifdef HAMOCC
      call xcbcst(rstfn_hamocc)
#endif

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
