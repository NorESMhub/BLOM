      subroutine getprecipfact_mct(lsend_precip_fact, precip_fact)

      ! Uses modules
      
      use types, only : r8
      use mod_xc

      implicit none

      ! Input/output arguments

      logical, intent(out)  :: lsend_precip_fact
      real(r8), intent(out) :: precip_fact

      ! Include files

#include "common_forc.h"

      lsend_precip_fact = sprfac
      precip_fact = prfac

      end subroutine getprecipfact_mct
