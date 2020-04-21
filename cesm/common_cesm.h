! ------------------------------------------------------------------------------
! Copyright (C) 2011-2015 Mats Bentsen
!
! This file is part of BLOM.
!
! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

c --- ------------------------------------------------------------------
c --- common blocks related to the application of CESM forcing fields
c --- ------------------------------------------------------------------
c
c --- daily average forcing fields
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::
     .  swa_da,        ! solar heat flux
     .  nsf_da,        ! non-solar heat flux
     .  hmlt_da,       ! heat flux due to melting
     .  lip_da,        ! liquid water flux
     .  sop_da,        ! solid precipitation
     .  eva_da,        ! evaporation
     .  rnf_da,        ! runoff, liquid
     .  rfi_da,        ! runoff, frozen
     .  fmltfz_da,     ! fresh water flux due to melting/freezing
     .  sfl_da,        ! salt flux
     .  ztx_da,        ! u component of wind stress
     .  mty_da,        ! v component of wind stress
     .  ustarw_da,     ! friction velocity for open water
     .  slp_da,        ! sea level pressure
     .  abswnd_da,     ! wind speed at measurement height -zu-
     .  atmco2_da,     ! atmospheric co2 concentration
     .  ficem_da       ! ice concentration
c
c --- time level indices for CESM fields
      integer l1ci,l2ci
c
      common /cesm/ swa_da,nsf_da,hmlt_da,lip_da,sop_da,eva_da,rnf_da,
     .              rfi_da,fmltfz_da,sfl_da,ztx_da,mty_da,ustarw_da,
     .              slp_da,abswnd_da,atmco2_da,ficem_da,l1ci,l2ci
