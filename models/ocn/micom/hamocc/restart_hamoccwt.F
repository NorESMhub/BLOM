      subroutine restart_hamoccwt(rstfnm_ocn)
c
cka write restart for HAMOCC      
c
      use mod_xc
c
      implicit none
c
#include "param_trc.h"
#include "common_bgc.h"
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_geo.h"
#include "common_trc.h"
c
      integer nt
      character rstfnm_ocn*(*)

      CALL AUFW_BGC(idm,jdm,kdm,bgc_dp,pgila,pgiph,bgc_pu            
     .             ,nyear,nmonth,nday,nstep/nphys,omask
     .             ,rstfnm_ocn,path2,path2_len)
c
      return
      end
