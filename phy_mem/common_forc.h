! ------------------------------------------------------------------------------
! Copyright (C) 2002-2020 Mats Bentsen, Jerry Tjiputra
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
c --- common blocks related to the application forcing fields that is
c --- shared among the various sources of forcing
c --- ------------------------------------------------------------------
c
c --- fields for diagnosed relaxation fluxes
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,48) ::
     .  tflxap,        ! heat flux to be applied
     .  sflxap,        ! fsalt flux to be applied
     .  tflxdi,        ! diagnosed heat flux
     .  sflxdi         ! diagnosed salt flux
c
c --- monthly climatological fields used in the computation of
c --- climatological fluxes and relaxation of sea surface temperature
c --- and salinity
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12) ::
     .  sstclm,        ! surface temperature
     .  ricclm,        ! ice concentration
     .  sssclm         ! sea surface salinity
c
c --- accumulation counter for diagnosed relaxation fluxes
      integer, dimension(48) :: nflxdi
c
c --- e-folding relaxation time scales
      real trxday,srxday
c
c --- maximum mixed layer depth for e-folding relaxation
      real trxdpt,srxdpt
c
c --- maximum absolute values for SST/SSS difference in relaxation
      real trxlim,srxlim
c
c --- interpolation parameters for monthly climatological fields
      real xmi
      integer l1mi,l2mi,l3mi,l4mi,l5mi
c
c --- flags concerning diagnosed heat and salt fluxes
      logical aptflx,apsflx,ditflx,disflx
c
c --- flags for balancing the SSS relaxation
      logical srxbal
c
c --- flag for smoothing of CESM forcing fields
      logical smtfrc
c
c --- flag for sending precipitation/runoff factor to CESM coupler
      logical sprfac
c
c --- name of file containing monthly SSS climatology
      character*256 scfile
c
      common /frc1/ tflxap,sflxap,tflxdi,sflxdi,sstclm,ricclm,sssclm,
     .              nflxdi,trxday,srxday,trxdpt,srxdpt,trxlim,srxlim,
     .              xmi,l1mi,l2mi,l3mi,l4mi,l5mi,
     .              aptflx,apsflx,ditflx,disflx,srxbal,smtfrc,sprfac,
     .              scfile
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
c
c --- various surface state and flux fields
     .  swa,           ! solar heat flux
     .  nsf,           ! non-solar heat flux
     .  hmltfz,        ! heat flux due to melting/freezing
     .  hmlt,          ! heat flux due to melting
     .  dfl,           ! derivative of non-solar heat flux by surface temp.
     .  lip,           ! liquid water flux
     .  sop,           ! solid precipitation
     .  eva,           ! evaporation
     .  rnf,           ! runoff, liquid
     .  rfi,           ! runoff, frozen
     .  fmltfz,        ! fresh water flux due to melting/freezing
     .  sfl,           ! salt flux
     .  ztx,           ! u component of wind stress
     .  mty,           ! v component of wind stress
     .  ustarw,        ! friction velocity for open water
     .  slp,           ! sea level pressure
     .  abswnd,        ! wind speed at measurement height -zu-
     .  albw,          ! daily mean open water albedo
     .  frzpot,        ! freezing potential
     .  mltpot,        ! melting potential
     .  atmco2,        ! atmospheric co2 concentration
     .  flxco2,        ! air-sea co2 flux
     .  flxdms,        ! sea-air dms flux
c
c --- albedo
     .  alb,         
c
c --- fields related to temporal smoothing of runoff
     .  rnfres,        ! runoff reservoar
     .  rnfflx,        ! liquid runoff freshwater flux taken out of the
                       ! reservoar
     .  rfiflx,        ! frozen runoff freshwater flux taken out of the
                       ! reservoar
c
c --- accumulation fields for balancing fresh water budget
     .  eiacc,         ! accumulation of evaporation and sea-ice melt./freezing
     .  pracc          ! accumulation of precipitation and runoff
c
c --- correction factor for precipitation and runoff
      real prfac
c
      common /frc2/ swa,nsf,hmltfz,hmlt,dfl,lip,sop,eva,rnf,rfi,fmltfz,
     .              sfl,ztx,mty,ustarw,slp,abswnd,albw,frzpot,mltpot,
     .              atmco2,flxco2,alb,rnfres,rnfflx,rfiflx,eiacc,pracc,
     .              prfac,flxdms
c
c --- constants set in 'frcdat'
      real albw_d,rhowat,sref
      integer nrfets
c
      common /frcpar/ albw_d,rhowat,sref,nrfets
