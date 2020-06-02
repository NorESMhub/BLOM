c Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger
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


      subroutine restart_hamoccwt(rstfnm_ocn)
c
c write restart for HAMOCC
c
      use mod_xc
      use mo_param1_bgc, only: ks,nsedtra,npowtra
      use mo_intfcblom,  only: omask
c
      implicit none
c
#include "param_trc.h"
#include "common_trc.h"
#include "common_blocks.h"
#include "common_clndr.h"
c
      character(len=*) :: rstfnm_ocn

      CALL AUFW_BGC(idm,jdm,kdm,ntr,ntrbgc,itrbgc,trc
     .             ,nyear,nmonth,nday,nstep,omask,rstfnm_ocn)
c
      return
      end
