! ------------------------------------------------------------------------------
! Copyright (C) 2015-2021 Mats Bentsen, Alok Kumar Gupta, Tomas Torsvik
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

      subroutine restart_ocntrcrd(rstfnm)
c
c --- ------------------------------------------------------------------
c --- Read ocean tracer state from restart file
c --- ------------------------------------------------------------------
c
      use mod_calendar, only: date_type, operator(/=)
      use mod_time, only: date0, time0, time
      use mod_xc
      use mod_dia, only : iotype
      use mod_nctools
      use mod_tracers, only:  ntrocn, natr, ntr, trc
c
      implicit none
c
      character rstfnm*(*)
c
      type(date_type) :: date_rest
      integer nt,nat
      real time0r,timer
      logical fexist
      character(len=256) :: trcnm
c
c --- ------------------------------------------------------------------
c --- If no ocean tracers are defined, return
c --- ------------------------------------------------------------------
c
      if (ntrocn.eq.0) return
c
c --- ------------------------------------------------------------------
c --- Check for file existence
c --- ------------------------------------------------------------------
c
      inquire(file=rstfnm,exist=fexist)
c
      call xcbcst(fexist)
c
c --- ------------------------------------------------------------------
c --- If file exists, read tracer data from file
c --- ------------------------------------------------------------------
c
      if (.not.fexist) then
        if (mnproc.eq.1) then
          write (lp,*)
     .'Warning! No tracer restart file found. Calling ocntrc_init...'
        endif
        call ocntrc_init
      else
        if (mnproc.eq.1) then
          write (lp,'(2a)') ' reading ocean tracer restart file ',
     .                      trim(rstfnm)
        endif  
        call ncfopn(rstfnm,'r',' ',1,iotype)
        call ncgeti('nday0',date_rest%day)
        call ncgeti('nmonth0',date_rest%month)
        call ncgeti('nyear0',date_rest%year)
        call ncgetr('time0',time0r)
        call ncgetr('time',timer)
        if (mnproc.eq.1) then
          if (date_rest /= date0 .or. 
     .        time0r .ne. time0 .or. timer.ne.time) then
            write (lp,'(2a)')
     .        ' Warning! The time information of the model and',
     .        ' restart file is inconsistent'
            write (lp,'(a)') '                 model          file'
            write (lp,'(a,i4.4,2(i2.2),a,i4.4,2(i2.2))')
     .                       ' date0:      ',date0, '      ',date_rest
            write (lp,'(a,2f14.4)')
     .                       ' time0: ',time0,time0r
            write (lp,'(a,2f14.4)')
     .                       ' time:  ',time,timer
          endif
        endif
        do nt=1,ntrocn-natr
          write (trcnm,'(a,i3.3)') 'trc',nt
          call ncread(trim(trcnm),trc(1-nbdy,1-nbdy,1,nt),ip,1,0.)
        enddo
        do nt=1,natr
          nat=ntr-natr+nt
          write (trcnm,'(a,i3.3)') 'atrc',nt
          call ncread(trim(trcnm),trc(1-nbdy,1-nbdy,1,nat),ip,1,0.)
        enddo
        call ncfcls
      endif
c
      return
      end
