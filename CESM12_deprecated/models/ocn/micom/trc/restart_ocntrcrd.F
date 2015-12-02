      subroutine restart_ocntrcrd(rstfnm_ocn)
c
c --- ------------------------------------------------------------------
c --- Read ocean tracer state from restart file
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mod_dia, only : iotype
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
      integer i,ndayr,nmonthr,nyearr,nt,nat
      real time0r,timer
      logical fexist
      character(len=80) :: rstfnm,trcnm
c
c --- ------------------------------------------------------------------
c --- If no ocean tracers are defined, return
c --- ------------------------------------------------------------------
c
      if (ntrocn.eq.0) return
c
c --- ------------------------------------------------------------------
c --- Generate file name and check for file existence
c --- ------------------------------------------------------------------
c
#ifdef CCSMCOUPLED
        i=1
        do while (rstfnm_ocn(i:i+8).ne.'.micom.r.')
          i=i+1
          if (i+8.gt.len(rstfnm_ocn)) then
            if (mnproc.eq.1) 
     &       write (lp,*)
     &        'Could not generate ocean tracer restart file name!'
            call xchalt('(restart_ocntrcrd)')
            stop '(restart_ocntrcrd)'
          endif
        enddo
        rstfnm=rstfnm_ocn(1:i-1)//'.micom.rtrc.'//rstfnm_ocn(i+9:)
#else
        i=1
        do while (rstfnm_ocn(i:i+8).ne.'_restphy_')
          i=i+1
          if (i+8.gt.len(rstfnm_ocn)) then
            if (mnproc.eq.1) 
     &       write (lp,*)
     &        'Could not generate ocean tracer restart file name!'
            call xchalt('(restart_ocntrcrd)')
            stop '(restart_ocntrcrd)'
          endif
        enddo
        rstfnm=rstfnm_ocn(1:i-1)//'_resttrc_'//rstfnm_ocn(i+9:)
#endif
        inquire(file=path2(1:path2_len)//rstfnm,exist=fexist)
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
          call ncfopn(path2(1:path2_len)//rstfnm,'r',' ',1,iotype)
          call ncgeti('nday0',ndayr)
          call ncgeti('nmonth0',nmonthr)
          call ncgeti('nyear0',nyearr)
          call ncgetr('time0',time0r)
          call ncgetr('time',timer)
        if (mnproc.eq.1) then
          if (ndayr.ne.nday0.or.nmonthr.ne.nmonth0.or.
     .        nyearr.ne.nyear0.or.time0r.ne.time0.or.
     .        timer.ne.time) then
            write (lp,'(2a)')
     .        ' Warning! The time information of the model and',
     .        ' restart file is inconsistent'
            write (lp,'(a)') '                  model          file'
            write (lp,'(a,2i14)')   ' nday0:  ',nday0,ndayr
            write (lp,'(a,2i14)')   ' nmonth0:',nmonth0,nmonthr
            write (lp,'(a,2i14)')   ' nyear0: ',nyear0,nyearr
            write (lp,'(a,2f14.4)') ' time0:  ',time0,time0r
            write (lp,'(a,2f14.4)') ' time:   ',time,timer
          endif
        endif
#ifdef ATRC
        do nt=1,ntrocn-natr
          write (trcnm,'(a,i3.3)') 'trc',nt
          call ncread(trim(trcnm),trc(1-nbdy,1-nbdy,1,nt),ip,1,0.)
        enddo
        do nt=1,natr
          nat=ntr-natr+nt
          write (trcnm,'(a,i3.3)') 'atrc',nt
          call ncread(trim(trcnm),trc(1-nbdy,1-nbdy,1,nat),ip,1,0.)
        enddo
#else
        do nt=1,ntrocn
          write (trcnm,'(a,i3.3)') 'trc',nt
          call ncread(trim(trcnm),trc(1-nbdy,1-nbdy,1,nt),ip,1,0.)
        enddo
#endif
          call ncfcls
      endif
c
      return
      end
