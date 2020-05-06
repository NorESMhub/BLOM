c Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger, M. Bentsen
c
c This file is part of BLOM/iHAMOCC.
c
c BLOM is free software: you can redistribute it and/or modify it under the
c terms of the GNU Lesser General Public License as published by the Free 
c Software Foundation, either version 3 of the License, or (at your option) 
c any later version. 
c
c BLOM is distributed in the hope that it will be useful, but WITHOUT ANY 
c WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
c FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
c more details. 
c
c You should have received a copy of the GNU Lesser General Public License 
c along with BLOM. If not, see https://www.gnu.org/licenses/.


      subroutine restart_hamoccrd(rstfnm_ocn)
c
c initialisation of HAMOCC      
c
      use mod_xc
      use mo_param1_bgc, only: ks,nsedtra,npowtra
c
      implicit none
c
#include "param_trc.h"
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_trc.h"
#include "common_bgc.h"
#ifndef sedbypass
#include "common_bgcs.h"
#endif
c
      integer          :: read_rest
      character(len=*) :: rstfnm_ocn
c
#ifdef RESTART_BGC
      read_rest=1
#else
      read_rest=0
#endif
c
c initialise arrays for BLOM-HAMOCC interface
      call inibgc
c
c initialise HAMOCC
c
      call blom2hamocc(2,1,kk,0)
c
      call INI_HAMOCC(read_rest,bgcdt,kpndtrun,idm,jdm,kdm,nbdy
     .               ,bgc_dp,bgc_t,bgc_s,bgc_rho
     .               ,bgc_dx,bgc_dy,bgc_pu,bgc_pw
     .               ,nyear,nmonth,nday,nstep1
     .               ,pglon,pglat,omask,ntr,ntrbgc,itrbgc,trc
#ifndef sedbypass
     .               ,sedlay2,powtra2,burial2
#endif
     .               ,rstfnm_ocn,path)
c

      return
      end
