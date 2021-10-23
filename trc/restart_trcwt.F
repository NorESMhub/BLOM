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
      use mod_config,   only: inst_suffix
      use mod_xc
c
      implicit none
c
      character rstfnm_ocn*(*)
c
      integer            :: i, leninrstfn
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
#ifdef CCSMCOUPLED
         leninrstfn = len('.blom'//trim(inst_suffix)//'.r.')-1
         i=1
         do while (rstfnm_ocn(i:i+leninrstfn) .ne.
     .        '.blom'//trim(inst_suffix)//'.r.')
            i=i+1
            if (i+leninrstfn .gt. len(rstfnm_ocn)) then
               write (lp,*)
     .              'Could not generate ocean tracer restart file name!'
               call xchalt('(restart_trcwt)')
               stop '(restart_trcwt)'
            endif
         enddo
         rstfnm_ocntrc = rstfnm_ocn(1:i-1)//'.blom'//trim(inst_suffix)//
     .        '.rtrc.'//rstfnm_ocn(i+leninrstfn+1:)
#ifdef HAMOCC
         rstfnm_hamocc = rstfnm_ocn(1:i-1)//'.blom'//trim(inst_suffix)//
     .        '.rbgc.'//rstfnm_ocn(i+leninrstfn+1:)
#endif
#else
         i=1
         do while (rstfnm_ocn(i:i+8) .ne. '_restphy_')
            i=i+1
            if (i+8 .gt. len(rstfnm_ocn)) then
               write (lp,*)
     .              'Could not generate ocean tracer restart file name!'
               call xchalt('(restart_trcwt)')
               stop '(restart_trcwt)'
            endif
         enddo
         rstfnm_ocntrc = rstfnm_ocn(1:i-1)//'_resttrc_'//
     .        rstfnm_ocn(i+9:)
#ifdef HAMOCC
         rstfnm_hamocc = rstfnm_ocn(1:i-1)//'_restbgc_'//
     .        rstfnm_ocn(i+9:)
#endif
#endif
      endif
c
c --- ------------------------------------------------------------------
c --- Write restart files
c --- ------------------------------------------------------------------
c
      call xcbcst(rstfnm_ocntrc)
      call restart_ocntrcwt(rstfnm_ocntrc)
c
#ifdef HAMOCC
      call xcbcst(rstfnm_hamocc)
      call restart_hamoccwt(rstfnm_hamocc)
#endif
c
      return
      end
