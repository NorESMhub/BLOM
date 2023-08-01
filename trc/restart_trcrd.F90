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

subroutine restart_trcrd(rstfnm_ocn)
!
! --- ------------------------------------------------------------------
! --- Read tracer state from restart files
! --- ------------------------------------------------------------------
!
  use mod_config, only: expcnf
  use mod_xc

  implicit none

  character(len=*), intent(in) :: rstfnm_ocn

  logical :: error
  character(len=256) :: rstfnm_ocntrc
  ! HAMOCC
  character(len=256) :: rstfnm_hamocc
#ifdef HAMOCC
  logical :: use_hamocc = .true.
#else
  logical :: use_hamocc = .false.
#endif
!
! --- ------------------------------------------------------------------
! --- Generate file name
! --- ------------------------------------------------------------------
!
  if (mnproc == 1) then
     if (expcnf == 'cesm') then
        call restart_getfile(rstfnm_ocn, 'rtrc', rstfnm_ocntrc, error)
        if (error) then
           write(lp,*) 'restart_trcrd: could not generate rstfnm_ocntrc file!'
           call xcstop('(restat_trcrd)')
           stop '(restart_trcrd)'
        endif
        if (use_HAMOCC) then
           call restart_getfile(rstfnm_ocn, 'rbgc', rstfnm_hamocc, error)
           if (error) then
              write(lp,*) 'restart_trcrd: could not generate rstfnm_hamocc file!'
              call xcstop('(restat_trcrd)')
              stop '(restart_trcrd)'
           endif
        end if
     else
        call restart_getfile(rstfnm_ocn, 'resttrc', rstfnm_ocntrc, error)
        if (error) then
           write(lp,*) 'restart_trcrd: could not generate rstfnm_ocntrc file!'
           call xcstop('(restat_trcrd)')
           stop '(restart_trcrd)'
        endif
        if (use_HAMOCC) then
           call restart_getfile(rstfnm_ocn, 'restbgc', rstfnm_hamocc, error)
           if (error) then
              write(lp,*) 'restart_trcrd: could not generate rstfnm_hamocc file!'
              call xcstop('(restat_trcrd)')
              stop '(restart_trcrd)'
           endif
        end if
     endif
  endif

  call xcbcst(rstfnm_ocntrc)
  if (use_HAMOCC) then
     call xcbcst(rstfnm_hamocc)
  end if

!
! --- ------------------------------------------------------------------
! --- Read restart files
! --- ------------------------------------------------------------------
!
  call restart_ocntrcrd(rstfnm_ocntrc)

  if (use_HAMOCC) then
     call hamocc_init(1,rstfnm_hamocc)
  end if

  return
end subroutine restart_trcrd
