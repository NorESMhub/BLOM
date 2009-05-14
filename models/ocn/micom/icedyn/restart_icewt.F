      subroutine restart_icewt(runid,runid_len,time,rstfnm_ext,
     .                         path2,path2_len)
c
c --- ------------------------------------------------------------------
c --- Write ice dynamics state to restart file
c --- ------------------------------------------------------------------
c
      use mod_xc, only: mnproc,nfu
c
      implicit none
c
      character runid*(*),rstfnm_ext*(*),path2*(*)
      real time
      integer runid_len,path2_len
c
#include "param_icedyn.h"
#include "common_icedyn.h"
c
      character fname*120
c
      fname=path2(1:path2_len)//runid(1:runid_len)//'_restice'//
     .      rstfnm_ext
c
      if (mnproc.eq.1) then
        write (*,'(a,a)') ' saving ice dynamics restart file ',
     .                    trim(fname)
        open (unit=nfu,file=fname,form='unformatted',status='unknown',
     .        action='write')
        write (nfu) time
        write (nfu) lold,lnew,u,v
        close (unit=nfu)
      endif
c
      return
      end
