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


      subroutine hamocc_init
c
c initialisation of HAMOCC      
c
      use mod_xc
      use mo_param1_bgc, only: ks,nsedtra,npowtra
#ifdef DIFFAT
      use mo_satm
#endif
c
      implicit none
c
#include "param_trc.h"
#include "common_trc.h"
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_geo.h"
#include "common_bgc.h"
#ifndef sedbypass
#include "common_bgcs.h"
#endif
c
      character(len=1) :: cdummy
      real             :: bgcdt         
c
c initialise arrays for BLOM-HAMOCC interface
      call inibgc
c
      bgcdt=nphys*baclin
c
c initialise HAMOCC
c
      call blom2hamocc(2,1,kk,0)
c
      call INI_HAMOCC(0,idm,jdm,kdm,nbdy,
     .                nyear,nmonth,nday,
     .                bgcdt,bgc_dx,bgc_dy,bgc_dp,bgc_rho,omask,
     .                plon,plat,ntr,ntrbgc,itrbgc,trc,
#ifndef sedbypass
     .                sedlay2,powtra2,burial2, 
#endif
     .                cdummy,path)
c
#ifdef DIFFAT
c     initialize slab atmosphere
      call satm_init(scp2,bgcdt)
#endif
c
      return
      end
