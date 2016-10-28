      real function gettime()
c
c --- ------------------------------------------------------------------
c --- return time in seconds since last call to either initime of
c --- gettime
c --- ------------------------------------------------------------------
c
      use mod_xc, only: mnproc
      implicit none
c
#include "timing.h"
c
      if (mnproc.eq.1) then
        gettime=wtime()-wtimeold
        wtimeold=gettime+wtimeold
      else
        gettime=1.
      endif
c
      return
      end
