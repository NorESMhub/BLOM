      subroutine restart_trcwt(rstfnm_ext)
c
c --- ------------------------------------------------------------------
c --- Write tracer state to restart file
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      character rstfnm_ext*(*)
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_dia.h"
#include "param_trc.h"
#include "common_trc.h"
c
      integer nt
      character*120 fname,trcname
c
      character startdate*20,timeunits*30
      real datenum
c
c --- set time information
      timeunits='hours since 0001-01-01 00:00'
      write(startdate,'(i4,a1,i2,a1,i2,a6)')
     .  nyear0,'-',nmonth0,'-',nday0,' 00:00'
      datenum=(time-.5*diag2dfq/nstep_in_day-time0)*24.
c
c --- set file name
      fname=path2(1:path2_len)//runid(1:runid_len)//'_resttrc'//
     .      rstfnm_ext
c 
c --- create file and save tracer data 
      if (mnproc.eq.1) 
     .  write (lp,'(a,a)') ' saving tracer restart file ',trim(fname)
      call ncfopn(trim(fname),'w')
c 
#ifdef COMPRES 
      call ncdimc('pcomp',ip,0)
#else 
      call ncdims('x',itdm)
      call ncdims('y',jtdm)
#endif 

      call nctime(datenum,calendar,timeunits,startdate)
c
      call ncputr('time',time,1)      
c
      do nt=1,ntr
        write(trcname,'(a3,i3.3)') 'trc',nt
#ifdef COMPRES
        call nccomp(trim(trcname),'pcomp time',trc(1-nbdy,1-nbdy,1,nt),
     .    ip)
#else
        call ncwrtr(trim(trcname),'x y time',trc(1-nbdy,1-nbdy,1,nt),
     .    ip,1)
#endif
      enddo
c 
      call ncfcls
c
      return
      end
