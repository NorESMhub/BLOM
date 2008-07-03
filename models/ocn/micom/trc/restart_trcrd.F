      subroutine restart_trcrd(nday1,rstfnm_ext)
c
c --- ------------------------------------------------------------------
c --- Read tracer initial conditions from restart file
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      character rstfnm_ext*(*)
      integer nday1
c
#include "common_blocks.h"
#include "param_trc.h"
#include "common_trc.h"
c
      integer nt
      character*120 fname,trcname
      logical fexist
      real rfexist
c
      rfexist=0.  
      fname=path2(1:path2_len)//runid(1:runid_len)//'_resttrc'//
     .      rstfnm_ext
      inquire (file=fname,exist=fexist)
      if (mnproc.eq.1.and.fexist) rfexist=1.
      call xcmaxr(rfexist)
      if (rfexist.gt.0.5) then
        call ncfopn(trim(fname),'r')
        if (mnproc.eq.1) then 
          write (lp,'(a,a)') ' reading tracer restart file ',trim(fname)
          call ncgetr('time',time)
          if (nint(time).ne.nday1) then
             write (lp,'(a,i6.6,a)')
     .         ' Warning! Integration day ',nint(time),
     .         ' in tracer restart file differs from'
             write (lp,'(a,i6.6,a)')
     .         ' start day ',nday1,' in limits file'
          end if
        end if   
        do nt=1,ntr
          write(trcname,'(a3,i3.3)') 'trc',nt 
          call ncread(trim(trcname),trc(1-nbdy,1-nbdy,1,nt),ip,1,0.)
        enddo
        call ncfcls 
      else
        if (mnproc.eq.1)
     .    write (lp,*) 
     .      'Warning! No tracer restart file found. Calling initrc...'
        call initrc
      endif
c
      return
      end
