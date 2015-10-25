      subroutine getfrc
c
c --- ------------------------------------------------------------------
c --- Get forcing
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mod_ben02, only: getfrc_ben02clim, getfrc_ben02syn
c
      implicit none
c
#include "common_blocks.h"
c
      if     (expcnf.eq.'cesm') then
        call getfrc_cesm
      elseif (expcnf.eq.'ben02clim') then
        call getfrc_ben02clim
      elseif (expcnf.eq.'ben02syn') then
        call getfrc_ben02syn
      elseif (expcnf.eq.'isomip1') then
c       call getfrc_isomip1
      elseif (expcnf.eq.'isomip2') then
c       call getfrc_isomip2
      else
        if (mnproc.eq.1) then
          write (lp,'(3a)') ' expcnf=',trim(expcnf),' is unsupported!'
        endif
        call xcstop('(getfrc)')
               stop '(getfrc)'
      endif
c
      return
      end
