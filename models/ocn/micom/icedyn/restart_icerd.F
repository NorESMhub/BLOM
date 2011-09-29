      subroutine restart_icerd(nday1,rstfnm_ocn,path2,path2_len)
c
c --- ------------------------------------------------------------------
c --- Read ice dynamics initial conditions from restart file
c --- ------------------------------------------------------------------
c
      use mod_xc, only: mnproc,nfu,xchalt
c
      implicit none
c
      integer nday1,path2_len
      character rstfnm_ocn*(*),path2*(*)
c
#include "param_icedyn.h"
#include "common_icedyn.h"
c
      integer i,rolen
      character rstfnm*80
      logical fexist
      real time
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
        inquire (file=path2(1:path2_len)//rstfnm,exist=fexist)
c
        if (fexist) then
          open (unit=nfu,file=path2(1:path2_len)//rstfnm,
     .          form='unformatted')
          write (*,'(a,a)') ' reading ice dynamics restart file ',
     .                      trim(path2(1:path2_len)//rstfnm)
          read (nfu) time
          if (nint(time).ne.nday1) then
            write (*,'(a,i6.6,a)')
     .        ' Warning! Integration day ',nint(time),
     .        ' in ice dynamics restart file differs from'
            write (*,'(a,i6.6,a)')
     .        ' start day ',nday1,' in limits file'
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
