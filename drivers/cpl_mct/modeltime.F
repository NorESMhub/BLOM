      subroutine modeltime(ymd, tod)

      ! Uses modules

      use mod_xc
      
      implicit none

      ! Input/output arguments

      integer, intent(out)    :: ymd, tod

      ! Include files

#include "common_blocks.h"
#include "common_clndr.h"

      ! Local variables


      ! External functions


      ymd = nyear*10000 + nmonth*100 + nday
      tod = nint(mod(nstep,nstep_in_day)*baclin)

      end subroutine modeltime
