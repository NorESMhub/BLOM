      subroutine restart_icerd(runid,runid_len,nday1,rstfnm_ext,
     .                         path2,path2_len)
c
c --- ------------------------------------------------------------------
c --- Read ice dynamics initial conditions from restart file
c --- ------------------------------------------------------------------
c
      use mod_xc, only: mnproc,nfu
c
      implicit none
c
      character runid*(*),rstfnm_ext*(*),path2*(*)
      integer runid_len,nday1,path2_len
c
#include "param_icedyn.h"
#include "common_icedyn.h"
c
      character fname*120
      logical fexist
      real time
c
      fname=path2(1:path2_len)//runid(1:runid_len)//'_restice'//
     .      rstfnm_ext
c
      inquire (file=fname,exist=fexist)
c
      if (fexist) then
        if (mnproc.eq.1) then
          open (unit=nfu,file=fname,form='unformatted')
          write (*,'(a,a)') ' reading ice dynamics restart file ',
     .                      trim(fname)
          read (nfu) time
          if (nint(time).ne.nday1) then
            write (*,'(a,i6.6,a)')
     .        ' Warning! Integration day ',nint(time),
     .        ' in ice dynamics restart file differs from'
            write (*,'(a,i6.6,a)')
     .        ' start day ',nday1,' in limits file'
          endif
          read (nfu) lold,lnew,u,v
          close (unit=nfu)
        endif
      else
        if (mnproc.eq.1) then
          write (*,*) 'Warning! No ice dynamics restart file found.'
        endif
      endif
c
      return
      end
