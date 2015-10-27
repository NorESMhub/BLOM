      subroutine inifrc
c
c --- ------------------------------------------------------------------
c --- Initialize forcing
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mod_ben02, only: inifrc_ben02clim, inifrc_ben02syn
c
      implicit none
c
#include "common_blocks.h"
c
c --- initialize forcing arrays with huge values?
c
      if     (expcnf.eq.'cesm') then
        call inifrc_cesm
      elseif (expcnf.eq.'ben02clim') then
        call inifrc_ben02clim
      elseif (expcnf.eq.'ben02syn') then
        call inifrc_ben02syn
      elseif (expcnf.eq.'isomip1') then
c       call inifrc_isomip1
      elseif (expcnf.eq.'isomip2') then
c       call inifrc_isomip2
      else
        if (mnproc.eq.1) then
          write (lp,'(3a)') ' expcnf=',trim(expcnf),' is unsupported!'
        endif
        call xcstop('(inifrc)')
               stop '(inifrc)'
      endif
c
      return
      end
