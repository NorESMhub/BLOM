      subroutine thermf(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- Get surface forcing functions
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "common_blocks.h"
c
      if     (expcnf.eq.'cesm') then
        call thermf_cesm(m,n,mm,nn,k1m,k1n)
      elseif (expcnf.eq.'ben02clim'.or.expcnf.eq.'ben02syn') then
        call thermf_ben02(m,n,mm,nn,k1m,k1n)
      elseif (expcnf.eq.'isomip1') then
c       call thermf_isomip1(m,n,mm,nn,k1m,k1n)
      elseif (expcnf.eq.'isomip2') then
c       call thermf_isomip2(m,n,mm,nn,k1m,k1n)
      else
        if (mnproc.eq.1) then
          write (lp,'(3a)') ' expcnf=',trim(expcnf),' is unsupported!'
        endif
        call xcstop('(thermf)')
               stop '(thermf)'
      endif
c
      return
      end
