      subroutine restart_icewt(time,rstfnm_ocn)
c
c --- ------------------------------------------------------------------
c --- Write ice dynamics state to restart file
c --- ------------------------------------------------------------------
c
      use mod_xc, only: mnproc,nfu,xchalt
c
      implicit none
c
      real time
      character rstfnm_ocn*(*)
c
#include "param_icedyn.h"
#include "common_icedyn.h"
c
      integer i
      character rstfnm*256
c
      if (mnproc.eq.1) then

        i=1
        do while (rstfnm_ocn(i:i+8).ne.'_restphy_')
          i=i+1
          if (i+8.gt.len(rstfnm_ocn)) then
            write (lp,*) 'Could not generate restart file name!'
            call xchalt('(restart_icerd)')
            stop '(restart_icerd)'
          endif
        enddo
        rstfnm=rstfnm_ocn(1:i-1)//'_restice_'//rstfnm_ocn(i+9:)
c
        write (*,'(a,a)') ' saving ice dynamics restart file ',
     .                    trim(rstfnm)
        open (unit=nfu,file=rstfnm,
     .        form='unformatted',status='unknown',action='write')
        write (nfu) time
        write (nfu) lold,lnew,u,v,uwat2l,vwat2l,gradhx2l,gradhy2l
        close (unit=nfu)
      endif
c
      return
      end
