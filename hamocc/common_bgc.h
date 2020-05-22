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
     .,bgc_rho (idm,jdm,kdm)
     .,omask   (idm,jdm)
c
      real bgc_dx,bgc_dy,bgc_dp
      real bgc_rho,omask
c
c----------------------------------------------------------------------
c number of physical timesteps in every bgc timestep
c
      integer nphys
      parameter (nphys=2)
c
c----------------------------------------------------------------------


