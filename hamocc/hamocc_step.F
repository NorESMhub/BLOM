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


      subroutine hamocc_step(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- perform one HAMOCC step     
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mo_bgcmean
      use mo_intfcblom, only: bgc_dx,bgc_dy,bgc_dp,bgc_rho,omask
      use mo_riverinpt, only: rivflx
      use mo_fedep,     only: get_fedep
      use mo_ndep,      only: get_ndep
c
      implicit none
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_forc.h"
#include "common_seaice.h"
#include "common_geo.h"
c
      integer, intent(in) :: m,n,mm,nn,k1m,k1n

      integer :: l,ldtday
      real    :: ndep(idm,jdm)
      real    :: dust(idm,jdm)
c
      call trc_limitc(nn)
c
      call blom2hamocc(m,n,mm,nn)
c
      ldtday = mod(nstep,nstep_in_day)
c
      do l=1,nbgc 
        bgcwrt(l)=0. 
        if (((diagann_bgc(l).and.nday_of_year.eq.1.or.diagmon_bgc(l)
     .    .and.nday.eq.1).and.mod(nstep,nstep_in_day).eq.0).or..not.
     .    (diagann_bgc(l).or.diagmon_bgc(l)).and.
     .    mod(nstep+.5,diagfq_bgc(l)).lt.1.) 
     .    bgcwrt(l)=1.
      enddo
c
      call get_fedep(idm,jdm,nmonth,dust)
      call get_ndep(idm,jdm,nyear,nmonth,omask,ndep)
c
      call hamocc4bcm(idm,jdm,kdm,nbdy,nyear,nmonth,nday,ldtday,
     .            bgc_dx,bgc_dy,bgc_dp,bgc_rho,plat,omask,
     .            dust,rivflx,ndep,
     .            swa,ficem,slp,abswnd,
     .            temp(1-nbdy,1-nbdy,1+nn),saln(1-nbdy,1-nbdy,1+nn),
     .            atmco2,flxco2,flxdms)
c
c --- accumulate fields and write output
c
      call accfields(idm,jdm,kdm,bgc_dx,bgc_dy,bgc_dp,omask)       

      call hamocc2blom(m,n,mm,nn)
c
      return
      end
