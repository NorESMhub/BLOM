      real function gettime()
c
c --- ------------------------------------------------------------------
c --- return time in seconds since last call to either initime of
c --- gettime
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "timing.h"
c
      start_count=end_count
      call system_clock(end_count)
      if (start_count.gt.end_count) then
        gettime=real(end_count-start_count+count_max)/real(count_rate)
      else
        gettime=real(end_count-start_count)/real(count_rate)
      endif
c
      return
      end
