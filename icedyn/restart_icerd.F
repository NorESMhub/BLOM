      subroutine restart_icerd(time,rstfnm_ocn)
c
c --- ------------------------------------------------------------------
c --- Read ice dynamics initial conditions from restart file
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
      integer i,rolen
      character rstfnm*256
      logical fexist
      real timer
c
      if (mnproc.eq.1) then
c
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
        inquire (file=rstfnm,exist=fexist)
c
        if (fexist) then
          open (unit=nfu,file=rstfnm,form='unformatted')
          write (*,'(a,a)') ' reading ice dynamics restart file ',
     .                      trim(rstfnm)
          read (nfu) timer
          if (nint(timer).ne.nint(time)) then
            write (*,'(a,i6.6,a)')
     .        ' Warning! Integration day ',nint(timer),
     .        ' in ice dynamics restart file differs from'
            write (*,'(a,i6.6,a)')
     .        ' integration day ',nint(time),' in ocean model'
          endif
          read (nfu) lold,lnew,u,v,uwat2l,vwat2l,gradhx2l,gradhy2l
          close (unit=nfu)
        else
          write (*,*) 'Warning! No ice dynamics restart file found.'
        endif
c
      endif
c
      return
      end
