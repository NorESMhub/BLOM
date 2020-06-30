! ------------------------------------------------------------------------------
! Copyright (C) 2015-2020 Mats Bentsen, Alok Kumar Gupta
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

      subroutine restart_ocntrcwt(rstfnm_ocn)
c
c --- ------------------------------------------------------------------
c --- Write ocean tracer state to restart file
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mod_dia, only :iotype,rstfmt,rstcmp
      use mod_nctools
c
      implicit none
c
      character rstfnm_ocn*(*)
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "param_trc.h"
#include "common_trc.h"
c
      integer i,nt,nat
      character(len=256) :: rstfnm,trcnm
c
c --- ------------------------------------------------------------------
c --- if no ocean tracers are defined, return
c --- ------------------------------------------------------------------
c
      if (ntrocn.eq.0) return
c
c --- ------------------------------------------------------------------
c --- Generate file name and create file 
c --- ------------------------------------------------------------------
c
#ifdef CCSMCOUPLED
      i=1
      do while (rstfnm_ocn(i:i+7).ne.'.blom.r.')
        i=i+1
        if (i+7.gt.len(rstfnm_ocn)) then
          if (mnproc.eq.1)
     .      write (lp,*)
     .        'Could not generate ocean tracer restart file name!'
          call xchalt('(restart_ocntrcrd)')
          stop '(restart_ocntrcrd)'
        endif
      enddo
      rstfnm=rstfnm_ocn(1:i-1)//'.blom.rtrc.'//rstfnm_ocn(i+8:)
#else
      i=1
      do while (rstfnm_ocn(i:i+8).ne.'_restphy_')
        i=i+1
        if (i+8.gt.len(rstfnm_ocn)) then
          if (mnproc.eq.1)
     .      write (lp,*)
     .        'Could not generate ocean tracer restart file name!'
          call xchalt('(restart_ocntrcrd)')
          stop '(restart_ocntrcrd)'
        endif
      enddo
      rstfnm=rstfnm_ocn(1:i-1)//'_resttrc_'//rstfnm_ocn(i+9:)
#endif
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
      call ncputi('nday0',nday0)
      call ncputi('nmonth0',nmonth0)
      call ncputi('nyear0',nyear0)
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
#ifdef ATRC
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
      
#else
      do nt=1,ntrocn
        write (trcnm,'(a,i3.3)') 'trc',nt
        if (rstcmp.eq.1) then
          call ncdefvar(trim(trcnm),'pcomp kk2 time',
     .                ndouble,8)
        else
          call ncdefvar(trim(trcnm),'x y kk2 time',
     .                ndouble,8)
        endif
      enddo
#endif
          call ncedef
#ifdef ATRC
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
#else
      do nt=1,ntrocn
        call flush(lp)
        write (trcnm,'(a,i3.3)') 'trc',nt
        if (rstcmp.eq.1) then
          call nccomp(trim(trcnm),'pcomp kk2 time',
     .                trc(1-nbdy,1-nbdy,1,nt),ip,1.,0.,8)
        else
          call ncwrtr(trim(trcnm),'x y kk2 time',
     .                trc(1-nbdy,1-nbdy,1,nt),ip,1,1.,0.,8)
        endif
      enddo
#endif
c
      call ncfcls
c
      return
      end
