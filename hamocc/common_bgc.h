c Copyright (C) 2020  K. Assmann, M. Bentsen, I. Bethke
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


c--------------------------------------------------------------------
c arrays to feed BCM variables to HAMOCC
c
      c o m m o n  /bgc_hamocc/
     . bgc_dx  (idm,jdm)      ,bgc_dy  (idm,jdm)
     .,bgc_dp  (idm,jdm,kdm)
     .,bgc_pu  (idm,jdm,kdm+1),bgc_pw  (idm,jdm,kdm+1)
     .,bgc_rho (idm,jdm,kdm)
     .,bgc_t   (idm,jdm,kdm)  ,bgc_s   (idm,jdm,kdm)
     .,omask   (idm,jdm)
     .,bgc_swr (idm,jdm)      ,bgc_fice(idm,jdm)
     .,bgc_awnd(idm,jdm)      ,bgc_slp(idm,jdm) 
     .,bgc_atmco2(idm,jdm)    ,bgc_flxco2(idm,jdm)
     .,bgc_flxdms(idm,jdm)
c
      real bgc_dx,bgc_dy,bgc_dp
      real bgc_pu,bgc_pw,bgc_rho,bgc_t,bgc_s,omask
      real bgc_swr,bgc_fice,bgc_awnd,bgc_slp
      real bgc_atmco2,bgc_flxco2,bgc_flxdms
c
      c o m m o n /bgc_hamocc_b/
     . ldtday,ldtmonth,kpndtrun
c
      integer ldtday,ldtmonth,kpndtrun
c
c----------------------------------------------------------------------
      c o m m o n/bgcc/
     . bgcdt
c
      real bgcdt
c
c----------------------------------------------------------------------
c number of physical timesteps in every bgc timestep
c
      integer nphys
      parameter (nphys=2)
c
c----------------------------------------------------------------------

      c o m m o n /bgc_hamocc2/
     . pglon,pglat,bgc3dwrt,bgc2dwrt

      real pglon(idm,jdm),pglat(idm,jdm)
      real bgc3dwrt,bgc2dwrt

