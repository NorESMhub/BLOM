      subroutine restart_hamoccwt(rstfnm_ocn)
c
c write restart for HAMOCC
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
      character(len=*) :: rstfnm_ocn

      CALL AUFW_BGC(idm,jdm,kdm,ntr,ntrbgc,itrbgc
     .             ,trc,bgc_dp,pglon,pglat            
     .             ,nyear,nmonth,nday,nstep/nphys,omask
     .             ,rstfnm_ocn,path2)
c
      return
      end
