      subroutine restart_hamoccwt(rstfnm_ext)
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
      character rstfnm_ext*(*)

      CALL AUFW_BGC(idm,jdm,kdm,bgc_dp,pgila,pgiph,bgc_pu            
     .             ,nyear,nmonth,nday,nstep/nphys,omask
     .             ,runid,runid_len,rstfnm_ext,path2,path2_len)
c
      return
      end
