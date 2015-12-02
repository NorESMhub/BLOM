      subroutine restart_hamoccrd(nday1,rstfnm_ocn)
c
c initialisation of HAMOCC      
c
      use mod_xc
      use mo_param1_bgc, only: ks,nsedtra,npowtra
c
      implicit none
c
#include "param_trc.h"
#include "common_bgc.h"
#include "common_bgcs.h"
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_trc.h"
c
      integer          :: nday1,read_rest
      character(len=*) :: rstfnm_ocn
c
#ifdef RESTART_BGC
      read_rest=1
#else
      read_rest=0
#endif
c
c initialise arrays for MICOM-HAMOCC interface
      call inibgc
c
c initialise HAMOCC
c
      call micom2hamocc(2,1,kk,0)
c
      call INI_HAMOCC(read_rest,1,bgcdt,kpndtrun,idm,jdm,kdm,nbdy
     .               ,bgc_dp,bgc_t,bgc_s,bgc_rho
     .               ,bgc_dx,bgc_dy,bgc_pu,bgc_pw
     .               ,nyear,nmonth,nday,nstep1
     .               ,pglon,pglat,omask,ntr,ntrbgc,itrbgc
     .               ,trc,sedlay2,powtra2,burial2
     .               ,rstfnm_ocn,path,path2)
c

      return
      end
