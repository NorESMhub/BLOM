      subroutine iniphy
c
c --- ------------------------------------------------------------------
c --- Initialize physical parameterizations
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
c
      if     (expcnf.eq.'cesm'.or.
     .        expcnf.eq.'ben02clim'.or.expcnf.eq.'ben02syn') then
        call read_twedon
      elseif (expcnf.eq.'isomip1') then
      elseif (expcnf.eq.'isomip2') then
      else
        if (mnproc.eq.1) then
          write (lp,'(3a)') ' expcnf=',trim(expcnf),' is unsupported!'
        endif
        call xcstop('(iniphy)')
               stop '(iniphy)'
      endif
c
      return
      end
