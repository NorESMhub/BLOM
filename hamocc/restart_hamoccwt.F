      subroutine restart_hamoccwt(rstfnm_ocn)
c
c write restart for HAMOCC
c
      use mod_xc
      use mo_param1_bgc, only: ks,nsedtra,npowtra
c
      implicit none
c
#include "param_trc.h"
#include "common_trc.h"
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_bgc.h"
#ifndef sedbypass
#include "common_bgcs.h"
#endif
c
      character(len=*) :: rstfnm_ocn

      CALL AUFW_BGC(idm,jdm,kdm,ntr,ntrbgc,itrbgc,trc
#ifndef sedbypass
     .             ,sedlay2,powtra2,burial2
#endif
     .             ,nyear,nmonth,nday,nstep,omask
     .             ,rstfnm_ocn,path2)
c
      return
      end
