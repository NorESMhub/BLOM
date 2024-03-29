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

      subroutine restart_ocntrcwt(rstfnm)
c
c --- ------------------------------------------------------------------
c --- Write ocean tracer state to restart file
c --- ------------------------------------------------------------------
c
      use mod_time, only: date0, time0, time
      use mod_xc
      use mod_dia, only: iotype, rstfmt, rstcmp
      use mod_nctools
      use mod_tracers, only:  ntrocn, natr, ntr, trc
c
      implicit none
c
      character rstfnm*(*)
c
      integer nt,nat
      character(len=256) :: trcnm
c
c --- ------------------------------------------------------------------
c --- if no ocean tracers are defined, return
c --- ------------------------------------------------------------------
c
      if (ntrocn.eq.0) return
c
c --- ------------------------------------------------------------------
c --- Create file
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write (lp,'(2a)') ' saving ocean tracer restart file ',
     .                    trim(rstfnm)
      endif
      if (rstfmt.eq.1) then
        call ncfopn(rstfnm,'w','6',1,iotype)
      elseif (rstfmt.eq.2) then
        call ncfopn(rstfnm,'w','h',1,iotype)
      else
        call ncfopn(rstfnm,'w','c',1,iotype)
      endif
c
c --- ------------------------------------------------------------------
c --- Create attributes and dimensions
c --- ------------------------------------------------------------------
c
      call ncputi('nday0',date0%day)
      call ncputi('nmonth0',date0%month)
      call ncputi('nyear0',date0%year)
      call ncputr('time0',time0)
      call ncputr('time',time)
      if (rstcmp.eq.1) then 
        call ncdimc('pcomp',ip,0)
      else
        call ncdims('x',itdm)
        call ncdims('y',jtdm)
      endif 
      call ncdims('kk2',2*kk)
      call ncdims('time',1)
c
c --- ------------------------------------------------------------------
c --- Write tracer data to file
c --- ------------------------------------------------------------------
c
      do nt=1,ntrocn-natr
        write (trcnm,'(a,i3.3)') 'trc',nt
        if (rstcmp.eq.1) then
          call ncdefvar(trim(trcnm),'pcomp kk2 time',
     .                ndouble,8)
        else
          call ncdefvar(trim(trcnm),'x y kk2 time',
     .                ndouble,8)
        endif
      enddo
      do nt=1,natr
        nat=ntr-natr+nt
        write (trcnm,'(a,i3.3)') 'atrc',nt
        if (rstcmp.eq.1) then
          call ncdefvar(trim(trcnm),'pcomp kk2 time',
     .                ndouble,8)
        else
          call ncdefvar(trim(trcnm),'x y kk2 time',
     .                ndouble,8)
        endif
      enddo
c
      call ncedef
c
      do nt=1,ntrocn-natr
        write (trcnm,'(a,i3.3)') 'trc',nt
        if (rstcmp.eq.1) then
          call nccomp(trim(trcnm),'pcomp kk2 time',
     .                trc(1-nbdy,1-nbdy,1,nt),ip,1.,0.,8)
        else
          call ncwrtr(trim(trcnm),'x y kk2 time',
     .                trc(1-nbdy,1-nbdy,1,nt),ip,1,1.,0.,8)
        endif
      enddo
      do nt=1,natr
        nat=ntr-natr+nt
        write (trcnm,'(a,i3.3)') 'atrc',nt
        if (rstcmp.eq.1) then
          call nccomp(trim(trcnm),'pcomp kk2 time',
     .                trc(1-nbdy,1-nbdy,1,nat),ip,1.,0.,8)
        else
          call ncwrtr(trim(trcnm),'x y kk2 time',
     .                trc(1-nbdy,1-nbdy,1,nat),ip,1,1.,0.,8)
        endif
      enddo
c
      call ncfcls
c
      return
      end
