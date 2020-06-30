! ------------------------------------------------------------------------------
! Copyright (C) 2002-2020 Mats Bentsen
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

      module mod_ben02
c
c --- ------------------------------------------------------------------
c --- This module contains routines for generating surface fluxes
c --- from NCEP or ERA40 reanalysis following Bentsen and Drange (2002).
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mod_ben02func
      use mod_bulktf
c
      implicit none
c
      private
c
      real, allocatable, dimension(:,:) ::
     .  atm_lon,       ! longitudes
     .  atm_lat,       ! latitudes
     .  atm_topo       ! topography
c
      integer, allocatable, dimension(:,:) ::
     .  atm_mask       ! mask
c
      real*4, allocatable, dimension(:,:,:) ::
     .  atm_wgt,       ! interpolation weights
     .  rnf_wgt        ! interpolation weights
c
      integer*2, allocatable, dimension(:,:,:) ::
     .  atm_iwgt,      ! interpolation adresses
     .  atm_jwgt,      ! interpolation adresses
     .  rnf_ocdpi,     ! interpolation adresses
     .  rnf_ocdpj      ! interpolation adresses
c
      real, allocatable, dimension(:,:,:) ::
     .  taud,          ! wind stress
     .  tauxd,tauyd,   ! wind stress components
     .  dswrfl,        ! short-wave surface radiation
     .  nlwrfs,        ! net long-wave surface radiation
     .  shtflx,        ! sensible heat flux
     .  lhtflx,        ! latent heat flux
     .  precip,        ! precipitation
     .  clouds,        ! cloud cover
     .  slpres,        ! sea level pressure
     .  runoff,        ! runoff
     .  tmpsfc,        ! surface temperature
     .  ricec          ! ice concentration
c
      real, allocatable, dimension(:,:) ::
     .  cd_d,          ! transfer coefficient for momentum (data state)
     .  ch_d,          ! transfer coefficient for sensible heat (data
                       ! state)
     .  ce_d,          ! transfer coefficient for latent heat (data
                       ! state)
     .  wg2_d,         ! gustiness squared (data state)
     .  cd_m,          ! transfer coefficient for momentum (model state)
     .  ch_m,          ! transfer coefficient for sensible heat (model
                       ! state)
     .  ce_m,          ! transfer coefficient for latent heat (model
                       ! state)
     .  wg2_m,         ! gustiness squared (model state)
     .  rhoa           ! air density
c
      real, allocatable, dimension(:,:) ::
     .  tsi_tda,       ! accumulated snow/ice surface temperature
     .  tml_tda,       ! accumulated mixed layer temperature
     .  sml_tda,       ! accumulated mixed layer salinity
     .  alb_tda,       ! accumulated albedo
     .  fice_tda,      ! accumulated sea ice concentration
     .  tsi            ! averaged snow/ice surface temperature
c
      integer ::
     .  ntda           ! accumulation number
c
      integer, dimension(itdm,jtdm) ::
     .  itp            ! mask for the full domain
c
      real ::
     .  atm_ice_swgt,  ! smoothing weight for fields over ice
     .  atm_rnf_swgt,  ! smoothing weight for fields over ice
     .  atm_crnf,      ! runoff adjustment factor
     .  atm_cswa,      ! short-wave radiation adjustment factor
     .  xgi            ! interpolation parameter
c
      integer ::
     .  atm_idm,       ! zonal dimension of atmospheric grid 
     .  atm_jdm,       ! meridional dimension of atmospheric grid 
     .  atm_path_len,  ! length of path to daily atmospheric foring
     .  atm_ice_nsmt,  ! number of smoothing iterations over ice
     .  atm_rnf_nsmt,  ! number of smoothing iterations over ice
     .  l1gi,l2gi,     ! time level indexes for interpolation
     .  l3gi,l4gi,     ! time level indexes for interpolation
     .  l5gi           ! time level indexes for interpolation
c
      character*80 ::
     .  atm_path       ! path to daily atmospheric foring 
c
      real, parameter ::
     .  atm_mval=-9999.,   ! value of a point that should not have a
                           ! physical value
     .  atm_fval=-99999.   ! value of a point that do not have a
                           ! physical value but should have one
c
      integer, parameter ::
     .  atm_idm_ncep=192,  ! zonal dimension of NCEP atmospheric grid
     .  atm_jdm_ncep=94,   ! meridional dimension of NCEP atmospheric
                           ! grid 
     .  atm_idm_era=320,   ! zonal dimension of ERA atmospheric grid
     .  atm_jdm_era=160,   ! meridional dimension of ERA atmospheric
                           ! grid
     .  atm_nwgt=12,       ! number of neighbours used in the
                           ! interpolation
     .  atm_abdm=10        ! max. number of runoff discharge basins for
                           ! each atmospheric grid point
c
      real ::
     .  atm_ice_csmt_ncep, ! constant determining how much the atm.
                           ! fields are smoothed over ice covered
                           ! regions (NCEP)
     .  atm_rnf_csmt_ncep, ! constant determining how much the runoff is
                           ! smoothed at the coastal discharge points
                           ! (NCEP)
     .  atm_crnf_ncep,     ! runoff adjustment factor (NCEP)
     .  atm_cswa_ncep,     ! short-wave radiation adjustment factor
                           ! (NCEP)
     .  atm_ice_csmt_era,  ! constant determining how much the atm.
                           ! fields are smoothed over ice covered
                           ! regions (NCEP)
     .  atm_rnf_csmt_era,  ! constant determining how much the runoff is
                           ! smoothed at the coastal discharge points
                           ! (NCEP)
     .  atm_crnf_era,      ! runoff adjustment factor (NCEP)
     .  atm_cswa_era       ! short-wave radiation adjustment factor
                           ! (NCEP)
c
      data atm_ice_csmt_ncep,atm_rnf_csmt_ncep /2.e14,1.e13/,
     .     atm_crnf_ncep,atm_cswa_ncep /0.82073,0.88340/,
     .     atm_ice_csmt_era,atm_rnf_csmt_era /0.0,1.e13/,
     .     atm_crnf_era,atm_cswa_era /0.7234,0.9721/
c
      real ::
     . zu,                 ! measurement height of wind [m]
     . zt,                 ! measurement height of temperature [m]
     . zq,                 ! measurement height of specific humidity [m]
     . emiss,              ! emissivity of water []
     . cpair,              ! specific heat of dry air [J K-1 kg-1]
     . stefanb,            ! stefan-boltzman constant [W m-2 K-4]
     . cd_r,               ! reference transfer coefficient of momentum
     . ch_r,               ! reference transfer coefficient of sensible
                           ! heat
     . ce_r,               ! reference transfer coefficient of latent
                           ! heat
     . rhoa_r,             ! reference air density [g cm-3]
     . wg2_r               ! reference gustiness squared [cm2 s-2]
c
      data zu,zt,zq /10.,10.,10./,
     .     emiss    /.97/,
     .     cpair    /1004.7/,
     .     stefanb  /5.67e-8/,
     .     cd_r     /1.2e-3/,
     .     ch_r     /1.1e-3/,
     .     ce_r     /1.3e-3/,
     .     wg2_r    /1.e3/,
     .     rhoa_r   /1.225e-3/
c
      integer ::
     . tciter              ! iterations in the computation of transfer
                           ! coefficients
c
      public :: initai,fnlzai,rdcsic,rdctsf,inifrc_ben02clim,
     .          inifrc_ben02syn,
     .          getfrc_ben02clim,getfrc_ben02syn,atm_path,atm_path_len,
     .          cd_d,ch_d,ce_d,wg2_d,cd_m,ch_m,ce_m,wg2_m,rhoa,
     .          tsi_tda,tml_tda,sml_tda,alb_tda,fice_tda,tsi,ntda
c
      contains
c
c --- ------------------------------------------------------------------
c
      subroutine rdatm_dim(filename)
c
c --- ------------------------------------------------------------------
c --- Read dimensions of atmospheric grid
c --- ------------------------------------------------------------------
c
      use netcdf
c
      implicit none
c
      character filename*(*)
c
      integer status,ncid,dimid
c
c --- Open netcdf file.
      status=nf90_open(filename,nf90_nowrite,ncid)
      if (status.ne.nf90_noerr) then
        write(lp,'(4a)') 'nf90_open: ',trim(filename),': ',
     .                   nf90_strerror(status)
        call xchalt('(rdatm_dim)')
               stop '(rdatm_dim)'
      endif
c
c --- Get dimensions.
      status=nf90_inq_dimid(ncid,'lon',dimid)
      if (status.ne.nf90_noerr) then
        write(lp,'(2a)') 'nf90_inq_dimid: lon: ',nf90_strerror(status)
        call xchalt('(rdatm_dim)')
               stop '(rdatm_dim)'
      endif
      status=nf90_inquire_dimension(ncid,dimid,len=atm_idm)
      if (status.ne.nf90_noerr) then
        write(lp,'(2a)') 'nf90_inquire_dimension: lon: ',
     .                   nf90_strerror(status)
        call xchalt('(rdatm_dim)')
               stop '(rdatm_dim)'
      endif
      status=nf90_inq_dimid(ncid,'lat',dimid)
      if (status.ne.nf90_noerr) then
        write(lp,'(2a)') 'nf90_inq_dimid: lat: ',nf90_strerror(status)
        call xchalt('(rdatm_dim)')
               stop '(rdatm_dim)'
      endif
      status=nf90_inquire_dimension(ncid,dimid,len=atm_jdm)
      if (status.ne.nf90_noerr) then
        write(lp,'(2a)') 'nf90_inquire_dimension: lat: ',
     .                   nf90_strerror(status)
        call xchalt('(rdatm_dim)')
               stop '(rdatm_dim)'
      endif
c
c --- Close netcdf file.
      status=nf90_close(ncid)
      if (status.ne.nf90_noerr) then
        write(lp,'(4a)') 'nf90_close: ',trim(filename),': ',
     .                   nf90_strerror(status)
        call xchalt('(rdatm_dim)')
               stop '(rdatm_dim)'
      endif
c
      return
      end subroutine rdatm_dim
c
c --- ------------------------------------------------------------------
c
      subroutine rdatm_llm(filename)
c
c --- ------------------------------------------------------------------
c --- Read atmospheric lon/lat coordinates and land mask from netcdf
c --- file
c --- ------------------------------------------------------------------
c
      use netcdf
c
      implicit none
c
      character filename*(*)
c
      integer status,ncid,varid,i,j
      real*4 r4lon(atm_idm),r4lat(atm_jdm)
      integer*2 i2field(atm_idm,atm_jdm)
c
c --- Open netcdf file.
      status=nf90_open(filename,nf90_nowrite,ncid)
      if (status.ne.nf90_noerr) then
        write(lp,'(4a)') 'nf90_open: ',trim(filename),': ',
     .                   nf90_strerror(status)
        call xchalt('(rdatm_llm)')
               stop '(rdatm_llm)'
      endif
c
c --- Read longitudes.
      status=nf90_inq_varid(ncid,'lon',varid)
      if (status.ne.nf90_noerr) then
        write(lp,'(2a)') 'nf90_inq_varid: lon: ',nf90_strerror(status)
        call xchalt('(rdatm_llm)')
               stop '(rdatm_llm)'
      endif
      status=nf90_get_var(ncid,varid,r4lon)
      if (status.ne.nf90_noerr) then
        write(lp,'(2a)') 'nf90_get_var: lon: ',nf90_strerror(status)
        call xchalt('(rdatm_llm)')
               stop '(rdatm_llm)'
      endif
c
c --- Read latitudes.
      status=nf90_inq_varid(ncid,'lat',varid)
      if (status.ne.nf90_noerr) then
        write(lp,'(2a)') 'nf90_inq_varid: lat: ',nf90_strerror(status)
        call xchalt('(rdatm_llm)')
               stop '(rdatm_llm)'
      endif
      status=nf90_get_var(ncid,varid,r4lat)
      if (status.ne.nf90_noerr) then
        write(lp,'(2a)') 'nf90_get_var: lat: ',nf90_strerror(status)
        call xchalt('(rdatm_llm)')
               stop '(rdatm_llm)'
      endif
c
c --- Read land mask.
      status=nf90_inq_varid(ncid,'land',varid)
      if (status.ne.nf90_noerr) then
        write(lp,'(2a)') 'nf90_inq_varid: land: ',nf90_strerror(status)
        call xchalt('(rdatm_llm)')
               stop '(rdatm_llm)'
      endif
      status=nf90_get_var(ncid,varid,i2field)
      if (status.ne.nf90_noerr) then
        write(lp,'(2a)') 'nf90_get_var: land: ',nf90_strerror(status)
        call xchalt('(rdatm_llm)')
               stop '(rdatm_llm)'
      endif
c
c --- Close netcdf file.
      status=nf90_close(ncid)
      if (status.ne.nf90_noerr) then
        write(lp,'(4a)') 'nf90_close: ',trim(filename),': ',
     .                   nf90_strerror(status)
        call xchalt('(rdatm_dim)')
               stop '(rdatm_dim)'
      endif
c
c$OMP PARALLEL DO PRIVATE(i)
      do j=1,atm_jdm
        do i=1,atm_idm
          atm_lon(i,j)=r4lon(i)
          atm_lat(i,j)=r4lat(j)
          atm_mask(i,j)=1-i2field(i,j)
        enddo
      enddo
c$OMP END PARALLEL DO
c
      end subroutine rdatm_llm
c
c --- ------------------------------------------------------------------
c
      subroutine rdatm_topo(filename)
c
c --- ------------------------------------------------------------------
c --- Read atmospheric lon/lat coordinates and land mask from netcdf
c --- file
c --- ------------------------------------------------------------------
c
      use netcdf
c
      implicit none
c
      character filename*(*)
c
      integer status,ncid,varid,i,j
      real*4 r4field(atm_idm,atm_jdm)
c
c --- Open netcdf file.
      status=nf90_open(filename,nf90_nowrite,ncid)
      if (status.ne.nf90_noerr) then
        write(lp,'(4a)') 'nf90_open: ',trim(filename),': ',
     .                   nf90_strerror(status)
        call xchalt('(rdatm_topo)')
               stop '(rdatm_topo)'
      endif
c
c --- Read topography.
      status=nf90_inq_varid(ncid,'hgtsfc',varid)
      if (status.ne.nf90_noerr) then
        write(lp,'(2a)') 'nf90_inq_varid: hgtsfc: ',
     .                   nf90_strerror(status)
        call xchalt('(rdatm_topo)')
               stop '(rdatm_topo)'
      endif
      status=nf90_get_var(ncid,varid,r4field)
      if (status.ne.nf90_noerr) then
        write(lp,'(2a)') 'nf90_get_var: hgtsfc: ',nf90_strerror(status)
        call xchalt('(rdatm_topo)')
               stop '(rdatm_topo)'
      endif
c
c --- Close netcdf file.
      status=nf90_close(ncid)
      if (status.ne.nf90_noerr) then
        write(lp,'(4a)') 'nf90_close: ',trim(filename),': ',
     .                   nf90_strerror(status)
        call xchalt('(rdatm_topo)')
               stop '(rdatm_topo)'
      endif
c
c --- Convert field to real topography field
c$OMP PARALLEL DO PRIVATE(i)
      do j=1,atm_jdm
        do i=1,atm_idm
          atm_topo(i,j)=r4field(i,atm_jdm+1-j)
        enddo
      enddo
c$OMP END PARALLEL DO
c
      end subroutine rdatm_topo
c
c --- ------------------------------------------------------------------
c
      subroutine rdatm_ts(filename,fieldname,time_step,atm_field)
c
c --- ------------------------------------------------------------------
c --- Read field at selected time step from NetCDF file
c --- ------------------------------------------------------------------
c
      use netcdf
c
      implicit none
c
      character filename*(*),fieldname*(*)
      integer time_step
      real atm_field(atm_idm,atm_jdm)
c
      integer status,ncid,dimid,atm_idm_t,atm_jdm_t,timeid,varid,
     .        vartype,start(3),count(3),stride(3),i,j
      integer*2 i2field(atm_idm,atm_jdm),i2_mval
      real time
      real*4 offset,scale_factor
c
c --- open netcdf file
      status=nf90_open(filename,nf90_nowrite,ncid)
      if (status.ne.nf90_noerr) then
        write(lp,'(4a)') 'nf90_open: ',trim(filename),': ',
     .                   nf90_strerror(status)
        call xchalt('(rdatm_ts)')
               stop '(rdatm_ts)'
      endif
c
c --- check dimensions
      status=nf90_inq_dimid(ncid,'lon',dimid)
      if (status.ne.nf90_noerr) then
        write(lp,'(2a)') 'nf90_inq_dimid: lon: ',nf90_strerror(status)
        call xchalt('(rdatm_ts)')
               stop '(rdatm_ts)'
      endif
      status=nf90_inquire_dimension(ncid,dimid,len=atm_idm_t)
      if (status.ne.nf90_noerr) then
        write(lp,'(2a)') 'nf90_inquire_dimension: lon: ',
     .                   nf90_strerror(status)
        call xchalt('(rdatm_ts)')
               stop '(rdatm_ts)'
      endif
      status=nf90_inq_dimid(ncid,'lat',dimid)
      if (status.ne.nf90_noerr) then
        write(lp,'(2a)') 'nf90_inq_dimid: lat: ',nf90_strerror(status)
        call xchalt('(rdatm_ts)')
               stop '(rdatm_ts)'
      endif
      status=nf90_inquire_dimension(ncid,dimid,len=atm_jdm_t)
      if (status.ne.nf90_noerr) then
        write(lp,'(2a)') 'nf90_inquire_dimension: lat: ',
     .                   nf90_strerror(status)
        call xchalt('(rdatm_ts)')
               stop '(rdatm_ts)'
      endif
      if (atm_idm_t.ne.atm_idm.or.atm_jdm_t.ne.atm_jdm) then
        write(lp,'(a,i4,a,i4,a,i4,a,i4,3a)')
     .    'expected dimensions ',atm_idm,' X ',atm_jdm,', but found ',
     .    atm_idm_t,' X ',atm_jdm_t,' in file ',trim(filename),'!'
        call xchalt('(rdatm_ts)')
               stop '(rdatm_ts)'
      endif
c
c --- get id of time variable
      status=nf90_inq_varid(ncid,'time',timeid)
      if (status.ne.nf90_noerr) then
        write(lp,'(2a)') 'nf90_inq_varid: time: ',nf90_strerror(status)
        call xchalt('(rdatm_ts)')
               stop '(rdatm_ts)'
      endif
c
c --- read time variable at selected time level
      status=nf90_get_var(ncid,timeid,time,(/time_step/))
      if (status.ne.nf90_noerr) then
        write(lp,'(2a)') 'nf90_get_var: time: ',nf90_strerror(status)
        call xchalt('(rdatm_ts)')
               stop '(rdatm_ts)'
      endif
      if (time.eq.0.) then
        write(lp,'(2a)') 'rdatm_ts: ',trim(filename)
        write(lp,'(a)')
     .    '  Time variable is zero at selected time level!'
        write(lp,'(a)') '  Probable cause is a corrupted file.'
        call xchalt('(rdatm_ts)')
               stop '(rdatm_ts)'
      endif
c
c --- get id of field variable
      status=nf90_inq_varid(ncid,fieldname,varid)
      if (status.ne.nf90_noerr) then
        write(lp,'(4a)') 'nf90_inq_varid: ',trim(fieldname),': ',
     .                   nf90_strerror(status)
        call xchalt('(rdatm_ts)')
               stop '(rdatm_ts)'
      endif
c
c --- get data type of field variable
      status=nf90_inquire_variable(ncid,varid,xtype=vartype)
      if (status.ne.nf90_noerr) then
        write(lp,'(4a)') 'nf90_inquire_variable: ',trim(fieldname),
     .                   ': ',nf90_strerror(status)
        call xchalt('(rdatm_ts)')
               stop '(rdatm_ts)'
      endif
c
c --- read selected time step of field
c
      start(1)=1
      start(2)=1
      start(3)=time_step
      count(1)=atm_idm
      count(2)=atm_jdm
      count(3)=1
      stride(1)=1
      stride(2)=1
      stride(3)=1
c
      if (vartype.eq.NF90_SHORT) then
c
        status=nf90_get_var(ncid,varid,i2field,start,count,stride)
        if (status.ne.nf90_noerr) then
          write(lp,'(4a)') 'nf90_get_var: ',trim(fieldname),': ',
     .                     nf90_strerror(status)
          call xchalt('(rdatm_ts)')
                 stop '(rdatm_ts)'
        endif
c
c --- - read offset, scale factor, and value of no data
        status=nf90_get_att(ncid,varid,'add_offset',offset)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_att: add_offset: ',
     .                     nf90_strerror(status)
          call xchalt('(rdatm_ts)')
                 stop '(rdatm_ts)'
        endif
        status=nf90_get_att(ncid,varid,'scale_factor',scale_factor)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_att: scale_factor: ',
     .                     nf90_strerror(status)
          call xchalt('(rdatm_ts)')
                 stop '(rdatm_ts)'
        endif
        status=nf90_get_att(ncid,varid,'missing_value',i2_mval)
        if (status.ne.nf90_noerr) then
          write(lp,'(2a)') 'nf90_get_att: missing_value: ',
     .                     nf90_strerror(status)
          call xchalt('(rdatm_ts)')
                 stop '(rdatm_ts)'
        endif
c
c --- - scale and add offset to field
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,atm_jdm
          do i=1,atm_idm
            if (i2field(i,j).eq.i2_mval) then
              atm_field(i,j)=atm_mval
            else
              atm_field(i,j)=i2field(i,j)*scale_factor+offset
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
c
      else
c
        status=nf90_get_var(ncid,varid,atm_field,start,count,stride)
        if (status.ne.nf90_noerr) then
          write(lp,'(4a)') 'nf90_get_var: ',trim(fieldname),': ',
     .                     nf90_strerror(status)
          call xchalt('(rdatm_ts)')
                 stop '(rdatm_ts)'
        endif
c
      endif
c
c --- close netcdf file
      status=nf90_close(ncid)
      if (status.ne.nf90_noerr) then
        write(lp,'(4a)') 'nf90_close: ',trim(filename),': ',
     .                   nf90_strerror(status)
        call xchalt('(rdatm_ts)')
               stop '(rdatm_ts)'
      endif
c
      return
      end subroutine rdatm_ts
c
c --- ------------------------------------------------------------------
c
      subroutine initai
c
c --- ------------------------------------------------------------------
c --- Prepare interpolation of atmospheric surface fields
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_blocks.h"
#include "common_geo.h"
c
      real, dimension(itdm,jtdm) :: tmp2da,tmp2db
      real, dimension(0:atm_nwgt) :: r_wgt
      real piloc,min_d,d,r,l2
      integer atm_idm_t,atm_jdm_t,idm_t,jdm_t,atm_nwgt_t,nw_2,is,js,it,
     .        jt,iso,jso,m,n,i,j,iii,ism1,isp1,jsm1,jsp1,itm1,itp1,jtm1,
     .        jtp1
      logical dimok,exists
      character filename*120
c
c --- Get global coordinates of p-points.
      call xcaget(tmp2da,plat,1)
      call xcaget(tmp2db,plon,1)
c
      if (mnproc.eq.1) then
c
c --- - Get dimensions of atmospheric grid
        filename=atm_path(1:atm_path_len)//'land.sfc.gauss.nc'
        call rdatm_dim(filename)
        if ((atm_idm.ne.atm_idm_ncep.or.atm_jdm.ne.atm_jdm_ncep).and.
     .      (atm_idm.ne.atm_idm_era .or.atm_jdm.ne.atm_jdm_era )) then
          write (lp,*) 'initai: unknown atmospheric dimensions!'
          call xchalt('(iniai)')
                 stop '(iniai)'
        endif
        if     (atm_idm.eq.atm_idm_ncep.and.
     .          atm_jdm.eq.atm_jdm_ncep) then
          write (lp,*) 'initai: using NCEP reanalysis'
        elseif (atm_idm.eq.atm_idm_era .and.
     .          atm_jdm.eq.atm_jdm_era ) then
          write (lp,*) 'initai: using ERA40 reanalysis'
        endif
c
c --- - Allocate memory for the atmospheric interpolation.
        allocate(atm_lon(atm_idm,atm_jdm),
     .           atm_lat(atm_idm,atm_jdm),
     .           atm_mask(atm_idm,atm_jdm),
     .           atm_topo(atm_idm,atm_jdm),
     .           atm_wgt(atm_nwgt,itdm,jtdm),
     .           atm_iwgt(atm_nwgt,itdm,jtdm),
     .           atm_jwgt(atm_nwgt,itdm,jtdm))
c
c --- - Read atmospheric lon/lat coordinates and land mask.
        filename=atm_path(1:atm_path_len)//'land.sfc.gauss.nc'
        call rdatm_llm(filename)
c
c --- - Read atmospheric model topography.
        filename=atm_path(1:atm_path_len)//'hgt.sfc.nc'
        call rdatm_topo(filename)
c
c --- - Read interpolation weights if they exist, otherwise compute them
        dimok=.false.
        filename='atm_intwgt.uf'
        inquire(file=filename,exist=exists)
        if (exists) then
          open (unit=nfu,file=filename,form='unformatted')
          read (nfu) atm_idm_t,atm_jdm_t,idm_t,jdm_t,atm_nwgt_t
          if (atm_idm_t.eq.atm_idm.and.atm_jdm_t.eq.atm_jdm.and.
     .        idm_t.eq.itdm.and.jdm_t.eq.jtdm.and.
     .        atm_nwgt_t.eq.atm_nwgt) dimok=.true.
        endif
        if (exists.and.dimok) then
          read (nfu) atm_wgt,atm_iwgt,atm_jwgt
          close (unit=nfu)
        else
          if (exists) close (unit=nfu)

          write (lp,*)
     .      'initai: computing atmospheric interpolation weights...'

          piloc=4.*atan(1.)
          r_wgt(0)=-1.
          nw_2=atm_nwgt/2+1
c
          is=1
          js=1
c
          do jt=1,jtdm
            do it=1,itdm
c
              min_d=spherdist(1.,atm_lon(is,js),atm_lat(is,js),
     .                           tmp2db(it,jt),tmp2da(it,jt))
c
 100          iso=is
              jso=js
c
              i=mod(iso-2+atm_idm,atm_idm)+1
              d=spherdist(1.,atm_lon(i,jso),atm_lat(i,jso),
     .                       tmp2db(it,jt),tmp2da(it,jt))
              if (d.lt.min_d) then
                is=i
                js=jso
                min_d=d
              endif
              i=mod(iso,atm_idm)+1
              d=spherdist(1.,atm_lon(i,jso),atm_lat(i,jso),
     .                       tmp2db(it,jt),tmp2da(it,jt))
              if (d.lt.min_d) then
                is=i
                js=jso
                min_d=d
              endif
              j=max(jso-1,1)
              d=spherdist(1.,atm_lon(iso,j),atm_lat(iso,j),
     .                       tmp2db(it,jt),tmp2da(it,jt))
              if (d.lt.min_d) then
                is=iso
                js=j
                min_d=d
              endif
              j=min(jso+1,atm_jdm)
              d=spherdist(1.,atm_lon(iso,j),atm_lat(iso,j),
     .                       tmp2db(it,jt),tmp2da(it,jt))
              if (d.lt.min_d) then
                is=iso
                js=j
                min_d=d
              endif
c
              if (is.ne.iso.or.js.ne.jso) goto 100
c
              do m=1,atm_nwgt
                r_wgt(m)=999999.
              enddo
c
              do j=min(atm_jdm-nw_2*2,max(1,       js-nw_2)),
     .             min(atm_jdm       ,max(nw_2*2+1,js+nw_2))
                do iii=is-nw_2,is+nw_2
                  i=mod(iii-1+atm_idm,atm_idm)+1
                  r=spherdist(1.,atm_lon(i ,j ),atm_lat(i ,j ),
     .                           tmp2db(it,jt),tmp2da(it,jt))
                  m=atm_nwgt+1
 10               m=m-1
                  if (r.lt.r_wgt(m)) goto 10
                  m=m+1
                  if (m.le.atm_nwgt) then
                    do n=atm_nwgt-1,m,-1
                      r_wgt(n+1)=r_wgt(n)
                      atm_iwgt(n+1,it,jt)=atm_iwgt(n,it,jt)
                      atm_jwgt(n+1,it,jt)=atm_jwgt(n,it,jt)
                    enddo
                    r_wgt(m)=r
                    atm_iwgt(m,it,jt)=i
                    atm_jwgt(m,it,jt)=j
                  endif
                enddo
              enddo
c
              itm1=min(itdm-2,max(1,it-1))
              itp1=min(itdm  ,max(3,it+1))
              jtm1=mod(jt-2+jtdm,jtdm)+1
              jtp1=mod(jt,jtdm)+1
              ism1=mod(is-2+atm_idm,atm_idm)+1
              isp1=mod(is,atm_idm)+1
              jsm1=min(atm_jdm-2,max(1,js-1))
              jsp1=min(atm_jdm  ,max(3,js+1))
c
              l2=
     .  .25*max( spherdist(1.,tmp2db(itm1,jt),tmp2da(itm1,jt),
     .                        tmp2db(itp1,jt),tmp2da(itp1,jt))
     .          *spherdist(1.,tmp2db(it,jtm1),tmp2da(it,jtm1),
     .                        tmp2db(it,jtp1),tmp2da(it,jtp1)),
     .           spherdist(1.,atm_lon(ism1,js ),atm_lat(ism1,js ),
     .                        atm_lon(isp1,js ),atm_lat(isp1,js ))
     .          *spherdist(1.,atm_lon(is ,jsm1),atm_lat(is ,jsm1),
     .                        atm_lon(is ,jsp1),atm_lat(is ,jsp1)))
     .  /piloc
c
              do m=1,atm_nwgt
                atm_wgt(m,it,jt)=max(exp(-.5*r_wgt(m)*r_wgt(m)/l2),
     .                              1.e-9)
              enddo
c
            enddo
          enddo
c
          if (nreg.eq.2) then
            do it=1,itdm
              do m=1,atm_nwgt
                atm_wgt(m,it,jtdm)=atm_wgt(m,itdm-it+1,jtdm-1)
                atm_iwgt(m,it,jtdm)=atm_iwgt(m,itdm-it+1,jtdm-1)
                atm_jwgt(m,it,jtdm)=atm_jwgt(m,itdm-it+1,jtdm-1)
              enddo
            enddo
          endif
c
          open (unit=nfu,file=filename,form='unformatted')
          write (nfu) atm_idm,atm_jdm,itdm,jtdm,atm_nwgt
          write (nfu) atm_wgt,atm_iwgt,atm_jwgt
          close (unit=nfu)
c
        endif
c 
      endif
c
c --- Get mask for the full domain
c$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do i=1,ii
          util1(i,j)=ip(i,j)
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcaget(tmp2da,util1,1)
      if (mnproc.eq.1) then
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,jtdm
          do i=1,itdm
            itp(i,j)=nint(tmp2da(i,j))
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
      return
      end subroutine initai
c
c --- ------------------------------------------------------------------
c
      subroutine fnlzai
c
c --- ------------------------------------------------------------------
c --- Finalize interpolation of atmospheric surface fields (deallocate
c --- memory)
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        deallocate(atm_lon,
     .             atm_lat,
     .             atm_mask,
     .             atm_topo,
     .             atm_wgt,
     .             atm_iwgt,
     .             atm_jwgt)
      endif
c
      return
      end subroutine fnlzai
c
c --- ------------------------------------------------------------------
c
      subroutine inta2o(adata,odata)
c
c --- ------------------------------------------------------------------
c --- Gaussian interpolation of 2D fields from atmospheric grid to ocean
c --- grid using precomputed weights.
c --- ------------------------------------------------------------------
c
      implicit none
c
      real, dimension(atm_idm,atm_jdm) :: adata
      real, dimension(itdm,jtdm) :: odata
c
      real w_sum,d_sum,w
      integer it,jt,n,is,js
c
c --- interpolate
c$OMP PARALLEL DO PRIVATE(it,w_sum,d_sum,n,is,js,w)
      do jt=1,jtdm
        do it=1,itdm
          if (itp(it,jt).eq.0) then
            odata(it,jt)=atm_mval
          else
            w_sum=0.
            d_sum=0.
            do n=1,atm_nwgt
              is=atm_iwgt(n,it,jt)
              js=atm_jwgt(n,it,jt)
              if (atm_mask(is,js).eq.1.and.
     .            adata(is,js).ne.atm_mval) then
                w=atm_wgt(n,it,jt)
                w_sum=w_sum+w
                d_sum=d_sum+adata(is,js)*w
              endif
            enddo
            if (w_sum.eq.0.) then
              odata(it,jt)=atm_fval
            else
              odata(it,jt)=d_sum/w_sum
            endif
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end subroutine inta2o
c
c --- ------------------------------------------------------------------
c
      subroutine smtfld(fld,msk,nsmt,swgt)
c
c --- ------------------------------------------------------------------
c --- Smooth the field, fld, a scale independent and conservative
c --- manner. The smoothing will be iterated nsmt times with a smoothing
c --- weight swgt and the smoothing is applied where the mask, msk, is
c --- equal one.
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_blocks.h"
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: fld
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: msk
      integer nsmt
      real swgt
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: usmtfl,vsmtfl
      integer mbdy,n,i,j
c
      mbdy=1
c
      do n=1,nsmt
c
        if (mbdy.eq.1) then
          mbdy=nbdy
          call xctilr(fld,1,1,nbdy,nbdy,halo_ps)
        else
          mbdy=mbdy-1
        endif
c
c$OMP PARALLEL DO PRIVATE(i)
        do j=1-mbdy+1,jj+mbdy
          do i=1-mbdy+1,ii+mbdy
            if (msk(i-1,j).eq.1.and.msk(i,j).eq.1) then
              usmtfl(i,j)=scuy(i,j)*scuxi(i,j)*(fld(i-1,j)-fld(i,j))
            else
              usmtfl(i,j)=0.
            endif
            if (msk(i,j-1).eq.1.and.msk(i,j).eq.1) then
              vsmtfl(i,j)=scvx(i,j)*scvyi(i,j)*(fld(i,j-1)-fld(i,j))
            else
              vsmtfl(i,j)=0.
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(i)
        do j=1-mbdy+1,jj+mbdy-1
          do i=1-mbdy+1,ii+mbdy-1
            if (msk(i,j).eq.1)
     .        fld(i,j)=fld(i,j)
     .                -swgt*scp2i(i,j)*(usmtfl(i+1,j)-usmtfl(i,j)
     .                                 +vsmtfl(i,j+1)-vsmtfl(i,j))
          enddo
        enddo
c$OMP END PARALLEL DO
c
      enddo
c
      return
      end subroutine smtfld
c
c --- ------------------------------------------------------------------
c
      subroutine uvrotr2g(u,v)
c
c --- ------------------------------------------------------------------
c --- Rotate zonal/meridional vector components at C-grid scalar points
c --- to model grid components at C-grid velocity points
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_geo.h"
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: u,v
c
      real latlim,rad
      parameter (latlim=87.,rad=1.74532925199432958e-02)
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: utmp,vtmp
      real mlat,dlat,dlon,psi
      integer i,j
c
c --- ------------------------------------------------------------------
c --- rotate vector components
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(i,mlat,dlat,dlon,psi)
      do j=1,jj
        do i=1,ii
c
c --- --- u-component
          if (iu(i,j).eq.0) then
            utmp(i,j)=atm_mval
          else
            mlat=.5*(plat(i,j)+plat(i-1,j))
            if (abs(mlat).gt.latlim) then
              utmp(i,j)=atm_fval
            else
              dlat=plat(i,j)-plat(i-1,j)
              dlon=mod(plon(i  ,j)+360.,360.)
     .            -mod(plon(i-1,j)+360.,360.)
              if (abs(dlon+360.).lt.abs(dlon)) dlon=dlon+360.
              if (abs(dlon-360.).lt.abs(dlon)) dlon=dlon-360.
              psi=atan2(dlat,cos(mlat*rad)*dlon)
              utmp(i,j)=.5*((u(i,j)+u(i-1,j))*cos(psi)
     .                     +(v(i,j)+v(i-1,j))*sin(psi))
            endif
          endif
c
c --- --- v-component
          if (iv(i,j).eq.0) then
            vtmp(i,j)=atm_mval
          else
            mlat=.5*(plat(i,j)+plat(i,j-1))
            if (abs(mlat).gt.latlim) then
              vtmp(i,j)=atm_fval
            else
              dlat=plat(i,j)-plat(i,j-1)
              dlon=mod(plon(i,j  )+360.,360.)
     .            -mod(plon(i,j-1)+360.,360.)
              if (abs(dlon+360.).lt.abs(dlon)) dlon=dlon+360.
              if (abs(dlon-360.).lt.abs(dlon)) dlon=dlon-360.
              psi=atan2(dlat,cos(mlat*rad)*dlon)
              vtmp(i,j)=.5*((u(i,j)+u(i,j-1))*cos(psi)
     .                     +(v(i,j)+v(i,j-1))*sin(psi))
            endif
          endif
c
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do i=1,ii
          u(i,j)=utmp(i,j)
          v(i,j)=vtmp(i,j)
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- extrapolate values to velocity points near pole singularity
c --- ------------------------------------------------------------------
c
      call fill_global(atm_mval,atm_fval,halo_uv,u)
      call fill_global(atm_mval,atm_fval,halo_vv,v)
c
      return
      end subroutine uvrotr2g
c
c --- ------------------------------------------------------------------
c
      subroutine albw_eval(dangle,lat,cc,albw_d,albw)
c
c --- ------------------------------------------------------------------
c --- compute 24 hrs mean albedo at the marine surface layer
c --- ------------------------------------------------------------------
c
      implicit none
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: lat,cc,albw
      real dangle,albw_d
c
      real pi2,deg,eepsil,fraci,absh2o,s0,decli,sundv,sin2,cos2,stot,
     .     sads,hangle,cosz,srad,sdir,sdif,altdeg,cfac,ssurf,albdir
      integer ifrac,i,j,l,npart
c
c --- set various quantities
c
      pi2=8.*atan(1.)          !        2 times pi
      deg=360./pi2             !        convert from radians to degrees
      eepsil=1.e-9             !        small number
c
      ifrac=24                 !        split each 12 hrs day into ifrac parts
      fraci=1./ifrac           !        1 over ifrac
c
      absh2o=0.09              ! ---    absorption of water and ozone
      s0=1365.                 ! w/m^2  solar constant
c
c --- compute astronomic quantities
c
      decli=.006918+.070257*sin(dangle)   -.399912*cos(dangle)
     $             +.000907*sin(2.*dangle)-.006758*cos(2.*dangle)
     $             +.001480*sin(3.*dangle)-.002697*cos(3.*dangle)
c
      sundv=1.00011+.001280*sin(dangle)   +.034221*cos(dangle)
     $             +.000077*sin(2.*dangle)+.000719*cos(2.*dangle)
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+ l,i,sin2,cos2,stot,sads,npart,hangle,cosz,srad,sdir,sdif,altdeg,
c$OMP+ cfac,ssurf,albdir)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c --- --- compute astronomic quantities
c
          sin2=sin(lat(i,j)/deg)*sin(decli)
          cos2=cos(lat(i,j)/deg)*cos(decli)
c
c --- --- split each day into ifrac parts, and compute the solar
c --- --- radiance for each part. by assuming symmetry of the irradiance
c --- --- about noon, it is sufficient to compute the irradiance for the
c --- --- first 12 hrs of the (24 hrs) day (mean for the first 12 hrs
c --- --- equals then the mean for the last 12 hrs)
c
          stot=0.
          sads=0.
c
          do npart=1,ifrac
c
c --- ----- hour angle, in radians
            hangle=pi2*(npart-.5)*fraci*.5
c
c --- ----- cosine of the zenith angle
            cosz=min(1.,max(0.,sin2+cos2*cos(hangle)))
c
c --- ----- extraterrestrial radiation
            srad =s0*sundv*cosz
c
c --- ----- direct radiation component
            sdir=srad*0.7**min(100.,1./(cosz+eepsil))
c
c --- ----- diffusive radiation component
            sdif=((1.-absh2o)*srad-sdir)*.5
c
c --- ----- solar noon altitude in degrees
            altdeg=max(0.,asin(sin2+cos2))*deg
c
c --- ----- cloudiness correction
            cfac=(1.-0.62*cc(i,j)+0.0019*altdeg)
c
            ssurf=(sdir+sdif)*cfac+eepsil
            stot=stot+ssurf
c
c --- ----- albedo for direct light
            albdir=0.03*exp(0.742*acos(cosz)**2.866)
c
c --- ----- radiation weighted sum of direct albedo
            sads=sads+albdir*ssurf
c
          enddo
c
c --- --- daily mean albedo over water
          albw(i,j)=(1.-cc(i,j))*sads/stot+cc(i,j)*albw_d
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end subroutine albw_eval
c
c --- ------------------------------------------------------------------
c
      subroutine rdatm_syn
c
c --- ------------------------------------------------------------------
c --- Read and interpolate daily atmospheric forcing fields. The
c --- interpolated fields of heat fluxes, cloud cover, precipitation,
c --- and surface temperature are smoothed over ice covered regions.
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_geo.h"
c
      real, dimension(atm_idm,atm_jdm) :: atm_field,atm_skt,atm_tau
      real, dimension(itdm,jtdm) :: tmp2d
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: smtmsk
      integer i,j,l,nday_of_year1,nyear1,iii,jjj
      character cyear*4,filename*120,fieldname*5
c
c --- The actual fields to be read are 2 days ahead
      nday_of_year1=nday_of_year+2
      if (nday_of_year1.gt.nday_in_year) then
        nday_of_year1=nday_of_year1-nday_in_year
        nyear1=nyear+1
      else
        nyear1=nyear
      endif
c
c --- Do not go beyond the start of the atm. time series
      if     (atm_idm.eq.atm_idm_ncep.and.atm_jdm.eq.atm_jdm_ncep) then
        if (nyear1.lt.1948) then
          nday_of_year1=1
          nyear1=1948
        endif
      elseif (atm_idm.eq.atm_idm_era .and.atm_jdm.eq.atm_jdm_era ) then
        if (nyear1.lt.1958) then
          nday_of_year1=1
          nyear1=1958
        endif
      endif
c
      write(cyear,'(i4)') nyear1
c
c --- Rearrange the time level indexes
      i=l1gi
      l1gi=l2gi
      l2gi=l3gi
      l3gi=l4gi
      l4gi=l5gi
      l5gi=i
c
      if (mnproc.eq.1) then
        write (lp,'(a,i3,a,i4,a)')
     .    'Reading and interpolating atm. forcing fields for day ',
     .    nday_of_year1,' of year ',nyear1,'...'
        write (lp,'(a,5i2)') 'Time level indexes:',
     .                       l1gi,l2gi,l3gi,l4gi,l5gi
      endif
c
c --- ------------------------------------------------------------------
c --- Read/interpolate ice concentration [0-1]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'daily/ICECsfc/icec.sfc.gauss.'//cyear//'.nc'
        fieldname='icec'
        call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
        call inta2o(atm_field,tmp2d)
      endif
      call xcaput(tmp2d,ricec(1-nbdy,1-nbdy,l5gi),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 ricec(1-nbdy,1-nbdy,l5gi))
c
c --- create smoothing mask - smooth where ice conc. is above 0.5
c$OMP PARALLEL DO PRIVATE(i)
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          if (ricec(i,j,l5gi).gt..5.and.ip(i,j).eq.1) then
            smtmsk(i,j)=1
          else
            smtmsk(i,j)=0
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Read/interpolate short-wave radiation flux [W/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'daily/DSWRFsfc/dswrf.sfc.gauss.'//cyear//'.nc'
        fieldname='dswrf'
        call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
        call inta2o(atm_field,tmp2d)
      endif
      call xcaput(tmp2d,dswrfl(1-nbdy,1-nbdy,l5gi),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 dswrfl(1-nbdy,1-nbdy,l5gi))
      call smtfld(dswrfl(1-nbdy,1-nbdy,l5gi),smtmsk,
     .            atm_ice_nsmt,atm_ice_swgt)
c
c --- Adjust short-wave radiation field
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          dswrfl(i,j,l5gi)=dswrfl(i,j,l5gi)*atm_cswa
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Read/interpolate net long-wave radiation flux [W/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'daily/NLWRSsfc/nlwrs.sfc.gauss.'//cyear//'.nc'
        fieldname='nlwrs'
        call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
        call inta2o(atm_field,tmp2d)
      endif
      call xcaput(tmp2d,nlwrfs(1-nbdy,1-nbdy,l5gi),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 nlwrfs(1-nbdy,1-nbdy,l5gi))
      call smtfld(nlwrfs(1-nbdy,1-nbdy,l5gi),smtmsk,
     .            atm_ice_nsmt,atm_ice_swgt)
c
c --- ------------------------------------------------------------------
c --- Read/interpolate total cloud cover [0-100%]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'daily/TCDCeatm/tcdc.eatm.gauss.'//cyear//'.nc'
        fieldname='tcdc'
        call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
        call inta2o(atm_field,tmp2d)
      endif
      call xcaput(tmp2d,clouds(1-nbdy,1-nbdy,l5gi),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 clouds(1-nbdy,1-nbdy,l5gi))
      call smtfld(clouds(1-nbdy,1-nbdy,l5gi),smtmsk,
     .            atm_ice_nsmt,atm_ice_swgt)
c
c --- Convert range of cloudiness from 0-100 to 0-1
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          clouds(i,j,l5gi)=clouds(i,j,l5gi)*1.e-2
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Read/interpolate precipitation rate [kg/m^2/s]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'daily/PRATEsfc/prate.sfc.gauss.'//cyear//'.nc'
        fieldname='prate'
        call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
        call inta2o(atm_field,tmp2d)
      endif
      call xcaput(tmp2d,precip(1-nbdy,1-nbdy,l5gi),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 precip(1-nbdy,1-nbdy,l5gi))
      call smtfld(precip(1-nbdy,1-nbdy,l5gi),smtmsk,
     .            atm_ice_nsmt,atm_ice_swgt)
c
c --- ------------------------------------------------------------------
c --- Read/interpolate latent heat net flux [W/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'daily/LHTFLsfc/lhtfl.sfc.gauss.'//cyear//'.nc'
        fieldname='lhtfl'
        call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
        call inta2o(atm_field,tmp2d)
      endif
      call xcaput(tmp2d,lhtflx(1-nbdy,1-nbdy,l5gi),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 lhtflx(1-nbdy,1-nbdy,l5gi))
      call smtfld(lhtflx(1-nbdy,1-nbdy,l5gi),smtmsk,
     .            atm_ice_nsmt,atm_ice_swgt)
c
c --- ------------------------------------------------------------------
c --- Read/interpolate sensible heat net flux [W/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'daily/SHTFLsfc/shtfl.sfc.gauss.'//cyear//'.nc'
        fieldname='shtfl'
        call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
        call inta2o(atm_field,tmp2d)
      endif
      call xcaput(tmp2d,shtflx(1-nbdy,1-nbdy,l5gi),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 shtflx(1-nbdy,1-nbdy,l5gi))
      call smtfld(shtflx(1-nbdy,1-nbdy,l5gi),smtmsk,
     .            atm_ice_nsmt,atm_ice_swgt)
c
c --- ------------------------------------------------------------------
c --- Read/interpolate sea surface temperature [K]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'daily/SKTsfc/skt.sfc.gauss.'//cyear//'.nc'
        fieldname='skt'
        call rdatm_ts(filename,fieldname,nday_of_year1,atm_skt)
c
c --- - Compute sea level temperature
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,atm_jdm
          do i=1,atm_idm
            atm_skt(i,j)=atm_skt(i,j)+.0065*atm_topo(i,j)
          enddo
        enddo
c$OMP END PARALLEL DO
c
        call inta2o(atm_skt,tmp2d)
      endif
      call xcaput(tmp2d,tmpsfc(1-nbdy,1-nbdy,l5gi),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 tmpsfc(1-nbdy,1-nbdy,l5gi))
      call smtfld(tmpsfc(1-nbdy,1-nbdy,l5gi),smtmsk,
     .            atm_ice_nsmt,atm_ice_swgt)
c
c --- ------------------------------------------------------------------
c --- Read/interpolate surface pressure [Pa]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'daily/PRESsfc/pres.sfc.gauss.'//cyear//'.nc'
        fieldname='pres'
        call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
c
c --- - Compute sea level pressure
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,atm_jdm
          do i=1,atm_idm
            atm_field(i,j)=atm_field(i,j)
     .                      *exp(9.81*atm_topo(i,j)
     .                      /(287.*(atm_skt(i,j)-.00325*atm_topo(i,j))))
          enddo
        enddo
c$OMP END PARALLEL DO
c
        call inta2o(atm_field,tmp2d)
      endif
      call xcaput(tmp2d,slpres(1-nbdy,1-nbdy,l5gi),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 slpres(1-nbdy,1-nbdy,l5gi))
c
c --- ------------------------------------------------------------------
c --- Read/interpolate momentum flux [N/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'daily/UFLXsfc/uflx.sfc.gauss.'//cyear//'.nc'
        fieldname='uflx'
        call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
        call inta2o(atm_field,tmp2d)
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,atm_jdm
          do i=1,atm_idm
            atm_tau(i,j)=atm_field(i,j)*atm_field(i,j)
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
      call xcaput(tmp2d,tauxd(1-nbdy,1-nbdy,l5gi),1)
      call fill_global(atm_mval,atm_fval,halo_pv,
     .                 tauxd(1-nbdy,1-nbdy,l5gi))
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'daily/VFLXsfc/vflx.sfc.gauss.'//cyear//'.nc'
        fieldname='vflx'
        call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
        call inta2o(atm_field,tmp2d)
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,atm_jdm
          do i=1,atm_idm
            atm_tau(i,j)=sqrt(atm_tau(i,j)
     .                        +atm_field(i,j)*atm_field(i,j))
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
      call xcaput(tmp2d,tauyd(1-nbdy,1-nbdy,l5gi),1)
      call fill_global(atm_mval,atm_fval,halo_pv,
     .                 tauyd(1-nbdy,1-nbdy,l5gi))
      if (mnproc.eq.1) then
        call inta2o(atm_tau,tmp2d)
      endif
      call xcaput(tmp2d,taud(1-nbdy,1-nbdy,l5gi),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 taud(1-nbdy,1-nbdy,l5gi))
c
c --- Change sign of momentum flux
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          tauxd(i,j,l5gi)=-tauxd(i,j,l5gi)
          tauyd(i,j,l5gi)=-tauyd(i,j,l5gi)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Rotate the vector components
      call xctilr(tauxd(1-nbdy,1-nbdy,l5gi), 1,1, 1,1, halo_pv)
      call xctilr(tauyd(1-nbdy,1-nbdy,l5gi), 1,1, 1,1, halo_pv)
      call uvrotr2g(tauxd(1-nbdy,1-nbdy,l5gi),tauyd(1-nbdy,1-nbdy,l5gi))
c
c --- ------------------------------------------------------------------
c --- Read/distribute runoff rate [kg/m^2/s]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'daily/RUNOFsfc/runof.sfc.gauss.'//cyear//'.nc'
        fieldname='runof'
        call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
c
c --- - Place runoff at ocean discharge points
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,jtdm
          do i=1,itdm
            tmp2d(i,j)=0.
          enddo
        enddo
c$OMP END PARALLEL DO
        do j=1,atm_jdm
          do i=1,atm_idm
            do l=1,atm_abdm
              if (rnf_wgt(l,i,j).gt.0.) then
                iii=rnf_ocdpi(l,i,j)
                jjj=rnf_ocdpj(l,i,j)
                tmp2d(iii,jjj)=tmp2d(iii,jjj)
     .            +atm_field(i,j)*rnf_wgt(l,i,j)
              endif
            enddo
          enddo
        enddo
      endif
      call xcaput(tmp2d,runoff(1-nbdy,1-nbdy,l5gi),1)
c
c --- Multiply runoff by 2 for fields prior to May 1957 
      if ( nyear1.lt.1957.or.
     .    (nyear1.eq.1957.and.nday_of_year1.le.120)) then
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            runoff(i,j,l5gi)=runoff(i,j,l5gi)*2.
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
c --- Adjust runoff field
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          runoff(i,j,l5gi)=runoff(i,j,l5gi)*atm_crnf
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Smooth the runoff field
      call smtfld(runoff(1-nbdy,1-nbdy,l5gi),ip,
     .            atm_rnf_nsmt,atm_rnf_swgt)
c
c --- Convert unit of runoff from kg/m^2/day to kg/m^2/s
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          runoff(i,j,l5gi)=runoff(i,j,l5gi)/86400.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'rdatm_syn:'
        endif
        call chksummsk(ricec(1-nbdy,1-nbdy,l5gi),ip,1,'ricec')
        call chksummsk(dswrfl(1-nbdy,1-nbdy,l5gi),ip,1,'dswrfl')
        call chksummsk(nlwrfs(1-nbdy,1-nbdy,l5gi),ip,1,'nlwrfs')
        call chksummsk(clouds(1-nbdy,1-nbdy,l5gi),ip,1,'clouds')
        call chksummsk(precip(1-nbdy,1-nbdy,l5gi),ip,1,'precip')
        call chksummsk(lhtflx(1-nbdy,1-nbdy,l5gi),ip,1,'lhtflx')
        call chksummsk(shtflx(1-nbdy,1-nbdy,l5gi),ip,1,'shtflx')
        call chksummsk(tmpsfc(1-nbdy,1-nbdy,l5gi),ip,1,'tmpsfc')
        call chksummsk(slpres(1-nbdy,1-nbdy,l5gi),ip,1,'slpres')
        call chksummsk(tauxd(1-nbdy,1-nbdy,l5gi),ip,1,'tauxd')
        call chksummsk(tauyd(1-nbdy,1-nbdy,l5gi),ip,1,'tauyd')
        call chksummsk(taud(1-nbdy,1-nbdy,l5gi),ip,1,'taud')
        call chksummsk(runoff(1-nbdy,1-nbdy,l5gi),ip,1,'runoff')
      endif
c
      return
      end subroutine rdatm_syn
c
c --- ------------------------------------------------------------------
c
      subroutine asflux
c
c --- ------------------------------------------------------------------
c --- Compute air-sea fluxes. Same routine is used both for
c --- climatological and synoptic forcing
c --- ------------------------------------------------------------------
c
      use mod_swtfrz, only: swtfrz
c
      implicit none
c
#include "common_blocks.h"
#include "common_forc.h"
#include "common_geo.h"
#include "common_clndr.h"
c
c --- Parameters:
c ---   dtmax - maximum near surface temperature gradient [K]
c ---   dqmax - maximum near surface specific humidity gradient [kg/kg]
c
      real dtmax,dqmax
      parameter (dtmax=30.,dqmax=.05)
c
      integer i,j,l,niter
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: taufac,cc
      real rntda,tau_d,prcp,dswrf,nlwrs,shtfl,lhtfl,slpr,tsrf_d,rice,
     .     tml,sml,albedo,fice,tice_f,tml_d,tsi_d,qsrf_d,le,ua,sa,ta,qa,
     .     tsrf_m,qsrf_m,dqsrf_m,dangle
c
      real intp1d
      external intp1d
c
      rntda=1./real(ntda)
      ntda=0
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+ l,i,tau_d,prcp,dswrf,nlwrs,shtfl,lhtfl,slpr,tsrf_d,rice,tml,sml,
c$OMP+ albedo,fice,tice_f,tml_d,tsi_d,qsrf_d,le,ua,sa,ta,qa,niter,
c$OMP+ tsrf_m,qsrf_m,dqsrf_m)
c$OMP+SHARED(xgi,zu,zt,zq)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c --- ------------------------------------------------------------------
c --- --- Interpolate the climatological/synoptic fields
c --- ------------------------------------------------------------------
c
          tau_d =intp1d(taud  (i,j,l1gi),taud  (i,j,l2gi),
     .                  taud  (i,j,l3gi),taud  (i,j,l4gi),
     .                  taud  (i,j,l5gi),xgi)
          prcp  =intp1d(precip(i,j,l1gi),precip(i,j,l2gi),
     .                  precip(i,j,l3gi),precip(i,j,l4gi),
     .                  precip(i,j,l5gi),xgi)
          dswrf =intp1d(dswrfl(i,j,l1gi),dswrfl(i,j,l2gi),
     .                  dswrfl(i,j,l3gi),dswrfl(i,j,l4gi),
     .                  dswrfl(i,j,l5gi),xgi)
          nlwrs =intp1d(nlwrfs(i,j,l1gi),nlwrfs(i,j,l2gi),
     .                  nlwrfs(i,j,l3gi),nlwrfs(i,j,l4gi),
     .                  nlwrfs(i,j,l5gi),xgi)
          shtfl =intp1d(shtflx(i,j,l1gi),shtflx(i,j,l2gi),
     .                  shtflx(i,j,l3gi),shtflx(i,j,l4gi),
     .                  shtflx(i,j,l5gi),xgi)
          lhtfl =intp1d(lhtflx(i,j,l1gi),lhtflx(i,j,l2gi),
     .                  lhtflx(i,j,l3gi),lhtflx(i,j,l4gi),
     .                  lhtflx(i,j,l5gi),xgi)
          slpr  =intp1d(slpres(i,j,l1gi),slpres(i,j,l2gi),
     .                  slpres(i,j,l3gi),slpres(i,j,l4gi),
     .                  slpres(i,j,l5gi),xgi)
          if (expcnf.eq.'ben02clim') then
            tsrf_d=intp1d(sstclm(i,j,l1gi),sstclm(i,j,l2gi),
     .                    sstclm(i,j,l3gi),sstclm(i,j,l4gi),
     .                    sstclm(i,j,l5gi),xgi)
            rice  =intp1d(ricclm(i,j,l1gi),ricclm(i,j,l2gi),
     .                    ricclm(i,j,l3gi),ricclm(i,j,l4gi),
     .                    ricclm(i,j,l5gi),xgi)
          else
            tsrf_d=intp1d(tmpsfc(i,j,l1gi),tmpsfc(i,j,l2gi),
     .                    tmpsfc(i,j,l3gi),tmpsfc(i,j,l4gi),
     .                    tmpsfc(i,j,l5gi),xgi)
            rice  =intp1d(ricec (i,j,l1gi),ricec (i,j,l2gi),
     .                    ricec (i,j,l3gi),ricec (i,j,l4gi),
     .                    ricec (i,j,l5gi),xgi)
          endif
          rnf(i,j)=intp1d(runoff(i,j,l1gi),runoff(i,j,l2gi),
     .                    runoff(i,j,l3gi),runoff(i,j,l4gi),
     .                    runoff(i,j,l5gi),xgi)
c
          prcp =max(0.,prcp)
          dswrf=max(0.,dswrf)
          rice =max(0.,min(1.,rice))
c
c --- ------------------------------------------------------------------
c --- --- Get averaged quantities obtained using the previous surface
c --- --- fluxes
c --- ------------------------------------------------------------------
c
          tsi(i,j)=tsi_tda(i,j)*rntda
          tml=tml_tda(i,j)*rntda
          sml=sml_tda(i,j)*rntda
          albedo=alb_tda(i,j)*rntda
          fice=fice_tda(i,j)*rntda
c
          tsi_tda(i,j)=0.
          tml_tda(i,j)=0.
          sml_tda(i,j)=0.
          alb_tda(i,j)=0.
          fice_tda(i,j)=0.
c
c --- --- Freezing temperature of sea water
          tice_f=swtfrz(0.,sml)+t0deg
c
c --- ------------------------------------------------------------------
c --- --- Compute the atmospheric state by using the prescribed momentum
c --- --- and heat fluxes and the prescribed sea surface state
c --- ------------------------------------------------------------------
c
          tml_d=max(tsrf_d,tice_f)
          tsi_d=max(200.,(tsrf_d-(1.-rice)*tml_d)/max(rice,1.e-6))
          qsrf_d=rice*qsati(tsi_d,slpr)+(1.-rice)*qsatw(tml_d,slpr)
          le=(2.501-0.00237*(tsrf_d-273.15))*1.e6
c
c --- --- Make sure wind stress is not too small compared to latent and
c --- --- sensible heat fluxes
          sa=max(abs(shtfl)/(rhoa(i,j)*cpair*ch_d(i,j)*dtmax),
     .           abs(lhtfl)/(rhoa(i,j)*le*ce_d(i,j)*dqmax))
          tau_d=max(tau_d,rhoa(i,j)*cd_d(i,j)*sa*sa)
c
c --- --- First guess on the atmospheric state by using the transfer
c --- --- coefficients and density from the previous time step
          ua=sqrt(.5*(-wg2_d(i,j)
     .                +sqrt(wg2_d(i,j)*wg2_d(i,j)
     .                     +4.*(tau_d/(rhoa(i,j)*cd_d(i,j)))**2)))
          sa=sqrt(ua*ua+wg2_d(i,j))
          ta=tsrf_d-.0098*zt-shtfl/(rhoa(i,j)*cpair*ch_d(i,j)*sa)
          qa=qsrf_d-lhtfl/(rhoa(i,j)*le*ce_d(i,j)*sa)
          rhoa(i,j)=rhoair(ta,qa,slpr)
c
c --- --- Iteration loop for estimating transfer coefficients and
c --- --- atmospheric state
          do niter=1,tciter
c
c --- ----- update the transfer coefficients and gustiness
            call bulktf(ua,zu,ta,zt,qa,zq,tsrf_d,qsrf_d,rice,
     .                  cd_d(i,j),ch_d(i,j),ce_d(i,j),wg2_d(i,j))
c
c --- ----- update the atmospheric state
            ua=sqrt(.5*(-wg2_d(i,j)
     .                  +sqrt(wg2_d(i,j)*wg2_d(i,j)
     .                       +4.*(tau_d/(rhoa(i,j)*cd_d(i,j)))**2)))
            sa=sqrt(ua*ua+wg2_d(i,j))
            ta=tsrf_d-.0098*zt-shtfl/(rhoa(i,j)*cpair*ch_d(i,j)*sa)
            qa=qsrf_d-lhtfl/(rhoa(i,j)*le*ce_d(i,j)*sa)
            rhoa(i,j)=rhoair(ta,qa,slpr)
c
          enddo
c
c --- ------------------------------------------------------------------
c --- --- Update transfer coefficients and gustiness with the computed
c --- --- atmospheric state the models ocean state
c --- ------------------------------------------------------------------
c
          tsrf_m=fice*tsi(i,j)+(1.-fice)*tml
          qsrf_m=fice*qsati(tsi(i,j),slpr)+(1.-fice)*qsatw(tml,slpr)
c
          do niter=1,tciter
            call bulktf(ua,zu,ta,zt,qa,zq,tsrf_m,qsrf_m,fice,
     .                  cd_m(i,j),ch_m(i,j),ce_m(i,j),wg2_m(i,j))
          enddo
c
c --- ------------------------------------------------------------------
c --- --- Compute correction of the wind stress on the surface and wind
c --- --- generated tke [m/s]
c --- ------------------------------------------------------------------
c
          sa=sqrt(ua*ua+wg2_m(i,j))
          taufac(i,j)=rhoa(i,j)*cd_m(i,j)*sa*ua/tau_d
c
c --- --- Wind generated TKE
          ustarw(i,j)=sqrt(cd_m(i,j)*sa*ua*rhoa(i,j)/rhowat)
c
c --- ------------------------------------------------------------------
c --- --- Compute heat fluxes (positive downward) [W/m^2]
c --- ------------------------------------------------------------------
c
          swa(i,j)=dswrf*(1.-albedo)
          le=(2.501-.00237*(tsrf_m-273.15))*1.e6
          nsf(i,j)=rhoa(i,j)*cpair*ch_m(i,j)*sa*(ta+0.0098*zt-tsrf_m)
     .            +rhoa(i,j)*ce_m(i,j)*le*sa*(qa-qsrf_m)
     .            -nlwrs-4.*emiss*stefanb*ta**3*(tsrf_m-tsrf_d)
c
c --- ------------------------------------------------------------------
c --- --- Compute evaporation (positive downward) [kg/m^2/s]
c --- ------------------------------------------------------------------
c
          eva(i,j)=rhoa(i,j)*ce_m(i,j)*sa*(qa-qsrf_m)
c
c --- ------------------------------------------------------------------
c --- --- Compute derivative of non-solar flux by surface temperature
c --- --- [W/m^2/K]
c --- ------------------------------------------------------------------
c
          dqsrf_m=fice*dqsati(tsi(i,j),slpr)+(1.-fice)*dqsatw(tml,slpr)
          dfl(i,j)=-rhoa(i,j)*cpair*ch_m(i,j)*sa
     .             -rhoa(i,j)*ce_m(i,j)*le*sa*dqsrf_m
     .             -4.*emiss*stefanb*ta**3
c
c --- ------------------------------------------------------------------
c --- --- Split solid and liquid precipitation (positive downward)
c --- --- [kg/m^2/s]
c --- ------------------------------------------------------------------
c
          if (ta.lt.t0deg) then
            lip(i,j)=0.
            sop(i,j)=prcp
          else
            lip(i,j)=prcp
            sop(i,j)=0.
          endif
c
c --- ------------------------------------------------------------------
c --- --- Sea level pressure [Pa] and wind speed at measurement height
c --- --- [m/s]
c --- ------------------------------------------------------------------
c
          slp(i,j)=slpr
          abswnd(i,j)=sa
c
c --- ------------------------------------------------------------------
c --- --- Frozen runoff is not computed.
c --- ------------------------------------------------------------------
c
          rfi(i,j)=0.
c
c --- ------------------------------------------------------------------
c --- --- If requested, apply correction to precipitation and runoff to
c --- --- balance the fresh water budget
c --- ------------------------------------------------------------------
c
          if (sprfac) then
            lip(i,j)=lip(i,j)*prfac
            sop(i,j)=sop(i,j)*prfac
            rnf(i,j)=rnf(i,j)*prfac
            rfi(i,j)=rfi(i,j)*prfac
          endif
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Compute the surface wind stress [N/m^2]
c --- ------------------------------------------------------------------
c
      call xctilr(taufac,  1,1, 1,1, halo_ps)
c
c$OMP PARALLEL DO PRIVATE(l,i) SHARED(xgi)
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          ztx(i,j)=.5*(taufac(i,j)+taufac(i-1,j))
     .               *intp1d(tauxd(i,j,l1gi),tauxd(i,j,l2gi),
     .                       tauxd(i,j,l3gi),tauxd(i,j,l4gi),
     .                       tauxd(i,j,l5gi),xgi)
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          mty(i,j)=.5*(taufac(i,j)+taufac(i,j-1))
     .               *intp1d(tauyd(i,j,l1gi),tauyd(i,j,l2gi),
     .                       tauyd(i,j,l3gi),tauyd(i,j,l4gi),
     .                       tauyd(i,j,l5gi),xgi)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Compute open water albedo
c --- ------------------------------------------------------------------
c 
c --- Time interpolation of cloud cover
c$OMP PARALLEL DO PRIVATE(l,i) SHARED(xgi)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          cc(i,j)=intp1d(clouds(i,j,l1gi),clouds(i,j,l2gi),
     .                   clouds(i,j,l3gi),clouds(i,j,l4gi),
     .                   clouds(i,j,l5gi),xgi)
          cc(i,j)=max(0.,min(1.,cc(i,j)))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Obtain albedo for open water
      dangle=8.*atan(1.)*real(nday_of_year-1)/real(nday_in_year)
      call albw_eval(dangle,plat,cc,albw_d,albw)
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'asflux:'
        endif
        call chksummsk(swa,ip,1,'swa')
        call chksummsk(nsf,ip,1,'nsf')
        call chksummsk(dfl,ip,1,'dfl')
        call chksummsk(lip,ip,1,'lip')
        call chksummsk(sop,ip,1,'sop')
        call chksummsk(eva,ip,1,'eva')
        call chksummsk(ztx,iu,1,'ztx')
        call chksummsk(mty,iv,1,'mty')
        call chksummsk(rnf,ip,1,'rnf')
        call chksummsk(ustarw,ip,1,'ustarw')
        call chksummsk(tsi,ip,1,'tsi')
        call chksummsk(slp,ip,1,'slp')
        call chksummsk(abswnd,ip,1,'abswnd')
        call chksummsk(albw,ip,1,'albw')
      endif
c
      return
      end subroutine asflux
c
c --- ------------------------------------------------------------------
c
      subroutine rdcsic
c
c --- ------------------------------------------------------------------
c --- Read and interpolate climatological sea-ice concentration
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_blocks.h"
#include "common_forc.h"
c
      real, dimension(atm_idm,atm_jdm) :: atm_field
      real*4, dimension(atm_idm,atm_jdm) :: atm_field_r4
      real, dimension(itdm,jtdm) :: tmp2d
      integer i,j,k
c
      if (mnproc.eq.1) then
        write (lp,*) 'reading atm. climatological ice concentration...'
        open (unit=nfu,file=atm_path(1:atm_path_len)
     .                      //'clim/icec_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read(nfu) atm_field_r4
c$OMP PARALLEL DO PRIVATE(i)
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=atm_field_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call inta2o(atm_field,tmp2d)
        endif
        call xcaput(tmp2d,ricclm(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   ricclm(1-nbdy,1-nbdy,k))
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(ricclm, 1,12, nbdy,nbdy, halo_ps)
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'rdcsic:'
        endif
        call chksummsk(ricclm,ip,12,'ricclm')
      endif
c
      return
      end subroutine rdcsic
c
c --- ------------------------------------------------------------------
c
      subroutine rdctsf
c
c --- ------------------------------------------------------------------
c --- Read and interpolate climatological surface temperature
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_blocks.h"
#include "common_forc.h"
c
      real, dimension(atm_idm,atm_jdm) :: atm_field
      real*4, dimension(atm_idm,atm_jdm) :: atm_field_r4
      real, dimension(itdm,jtdm) :: tmp2d
      real dx2,dy2
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: smtmsk
      integer i,j,k,l
c
      if (mnproc.eq.1) then
        write (lp,*)
     .    'reading atm. climatological surface temperature...'
      endif
c
c --- Compute smooting weight atm_ice_swgt. For stability
c --- atm_ice_swgt < .5*dx^2*dy^2/(dx^2+dy^2).
      atm_ice_swgt=spval
c$OMP PARALLEL DO PRIVATE(l,i,dx2,dy2) REDUCTION(min:atm_ice_swgt)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          dx2=scpx(i,j)*scpx(i,j)
          dy2=scpy(i,j)*scpy(i,j)
          atm_ice_swgt=min(atm_ice_swgt,.5*dx2*dy2/(dx2+dy2))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcmin(atm_ice_swgt)
      atm_ice_swgt=.9*atm_ice_swgt
c
c --- Number of smoothing iterations is choosen to get a geographical
c --- extent of the smoothing independent of the grid resolution
      if (mnproc.eq.1) then
        if     (atm_idm.eq.atm_idm_ncep.and.
     .          atm_jdm.eq.atm_jdm_ncep) then
          atm_ice_nsmt=nint(atm_ice_csmt_ncep/atm_ice_swgt)
        elseif (atm_idm.eq.atm_idm_era .and.
     .          atm_jdm.eq.atm_jdm_era ) then
          atm_ice_nsmt=nint(atm_ice_csmt_era/atm_ice_swgt)
        endif
      endif
      call xcbcst(atm_ice_nsmt)
      if (mnproc.eq.1) then
        write (lp,*) 'atm_ice_swgt',atm_ice_swgt,
     .               'atm_ice_nsmt',atm_ice_nsmt
      endif
c
      if (mnproc.eq.1) then
        open (unit=nfu,file=atm_path(1:atm_path_len)
     .                      //'clim/skt_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) atm_field_r4
c$OMP PARALLEL DO PRIVATE(i)
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=atm_field_r4(i,j)+.0065*atm_topo(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call inta2o(atm_field,tmp2d)
        endif
        call xcaput(tmp2d,sstclm(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   sstclm(1-nbdy,1-nbdy,k))
c
c --- - create smoothing mask - smooth where ice conc. is above 0.5
c$OMP PARALLEL DO PRIVATE(i)
        do j=1-nbdy,jj+nbdy
          do i=1-nbdy,ii+nbdy
            if (ricclm(i,j,k).gt..5.and.ip(i,j).eq.1) then
              smtmsk(i,j)=1
            else
              smtmsk(i,j)=0
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
c
        call smtfld(sstclm(1-nbdy,1-nbdy,k),smtmsk,
     .              atm_ice_nsmt,atm_ice_swgt)
c
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(sstclm, 1,12, nbdy,nbdy, halo_ps)
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'rdctsf:'
        endif
        call chksummsk(sstclm,ip,12,'sstclm')
      endif
c
      return
      end subroutine rdctsf
c
c --- ------------------------------------------------------------------
c
      subroutine inifrc_ben02clim
c
c --- ------------------------------------------------------------------
c --- Initialize monthly climatological forcing fields
c --- ------------------------------------------------------------------
c
      use mod_thdysi, only: albs_f
c
      implicit none
c
#include "common_blocks.h"
#include "common_forc.h"
#include "common_geo.h"
c
      real, allocatable, dimension(:,:,:) :: atm_sktclm
      real, allocatable, dimension(:,:) :: atm_field
      real*4, allocatable, dimension(:,:) :: atm_field_r4
      real, dimension(itdm,jtdm) :: tmp2d
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12) :: smtmsk
      real dx2,dy2,prc_sum,eva_sum,rnf_sum,swa_sum,lwa_sum,lht_sum,
     .     sht_sum,fwf_fac,dangle,garea,le,albedo,fac,swa_ave,lwa_ave,
     .     lht_ave,sht_ave,crnf,cswa
      real*4 rw4
      integer i,j,k,l,il,jl
      integer*2 rn2,ri2,rj2
c
c --- Prepare interpolation of surface fields.
      call initai
c
c --- Read climatological sea-ice concentration.
      call rdcsic
c
c --- Read climatological surface temperature.
      call rdctsf
c
c --- Initialize diagnosing/application of relaxation fluxes
      call idarlx
c
c --- Allocate memory for additional monthly forcing fields.
      allocate(taud  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12),
     .         tauxd (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12),
     .         tauyd (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12),
     .         dswrfl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12),
     .         nlwrfs(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12),
     .         shtflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12),
     .         lhtflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12),
     .         precip(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12),
     .         clouds(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12),
     .         slpres(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12),
     .         runoff(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12))
c
c --- Allocate memory for temporary fields
      if (mnproc.eq.1) then
        allocate(atm_sktclm(atm_idm,atm_jdm,12),
     .           atm_field(atm_idm,atm_jdm),
     .           atm_field_r4(atm_idm,atm_jdm))
      endif
c
c --- Compute smoothing weights atm_ice_swgt and atm_rnf_swgt. For
c --- stability atm_ice_swgt,atm_rnf_swgt < .5*dx^2*dy^2/(dx^2+dy^2).
      atm_ice_swgt=spval
      atm_rnf_swgt=spval
c$OMP PARALLEL DO PRIVATE(l,i,dx2,dy2)
c$OMP+REDUCTION(min:atm_ice_swgt,atm_rnf_swgt)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          dx2=scpx(i,j)*scpx(i,j)
          dy2=scpy(i,j)*scpy(i,j)
          atm_ice_swgt=min(atm_ice_swgt,.5*dx2*dy2/(dx2+dy2))
          atm_rnf_swgt=min(atm_rnf_swgt,.5*dx2*dy2/(dx2+dy2))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcmin(atm_ice_swgt)
      atm_ice_swgt=.9*atm_ice_swgt
      call xcmin(atm_rnf_swgt)
      atm_rnf_swgt=.9*atm_rnf_swgt
c
c --- Number of smoothing iterations is choosen to get a geographical
c --- extent of the smoothing independent of the grid resolution
      if (mnproc.eq.1) then
        if     (atm_idm.eq.atm_idm_ncep.and.
     .          atm_jdm.eq.atm_jdm_ncep) then
          atm_ice_nsmt=nint(atm_ice_csmt_ncep/atm_ice_swgt)
          atm_rnf_nsmt=nint(atm_rnf_csmt_ncep/atm_rnf_swgt)
        elseif (atm_idm.eq.atm_idm_era .and.
     .          atm_jdm.eq.atm_jdm_era ) then
          atm_ice_nsmt=nint(atm_ice_csmt_era/atm_ice_swgt)
          atm_rnf_nsmt=nint(atm_rnf_csmt_era/atm_rnf_swgt)
        endif
      endif
      call xcbcst(atm_ice_nsmt)
      call xcbcst(atm_rnf_nsmt)
      if (mnproc.eq.1) then
        write (lp,*) 'atm_ice_swgt',atm_ice_swgt,
     .               'atm_ice_nsmt',atm_ice_nsmt
        write (lp,*) 'atm_rnf_swgt',atm_rnf_swgt,
     .               'atm_rnf_nsmt',atm_rnf_nsmt
      endif
c
c --- create smoothing mask - smooth where ice conc. is above 0.5
      do k=1,12
c$OMP PARALLEL DO PRIVATE(i)
        do j=1-nbdy,jj+nbdy
          do i=1-nbdy,ii+nbdy
            if (ricclm(i,j,k).gt..5.and.ip(i,j).eq.1) then
              smtmsk(i,j,k)=1
            else
              smtmsk(i,j,k)=0
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c
c --- ------------------------------------------------------------------
c --- read short-wave radiation flux [W/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write (lp,*)
     .    'reading atm. climatological short-wave radiation flux...'
        open (unit=nfu,file=atm_path(1:atm_path_len)
     .                      //'clim/dswrf_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) atm_field_r4
c$OMP PARALLEL DO PRIVATE(i)
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=atm_field_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call inta2o(atm_field,tmp2d)
        endif
        call xcaput(tmp2d,dswrfl(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   dswrfl(1-nbdy,1-nbdy,k))
        call smtfld(dswrfl(1-nbdy,1-nbdy,k),smtmsk(1-nbdy,1-nbdy,k),
     .              atm_ice_nsmt,atm_ice_swgt)
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(dswrfl, 1,12, nbdy,nbdy, halo_ps)
c
c --- ------------------------------------------------------------------
c --- read net long-wave radiation flux [W/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write (lp,*)
     .    'reading atm. climatological long-wave radiation flux...'
        open (unit=nfu,file=atm_path(1:atm_path_len)
     .                      //'clim/nlwrs_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) atm_field_r4
c$OMP PARALLEL DO PRIVATE(i)
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=atm_field_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call inta2o(atm_field,tmp2d)
        endif
        call xcaput(tmp2d,nlwrfs(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   nlwrfs(1-nbdy,1-nbdy,k))
        call smtfld(nlwrfs(1-nbdy,1-nbdy,k),smtmsk(1-nbdy,1-nbdy,k),
     .              atm_ice_nsmt,atm_ice_swgt)
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(nlwrfs, 1,12, nbdy,nbdy, halo_ps)
c
c --- ------------------------------------------------------------------
c --- read total cloud cover [0-100%]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write(lp,*) 'reading atm. climatological total cloud cover...'
        open (unit=nfu,file=atm_path(1:atm_path_len)
     .                      //'clim/tcdc_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) atm_field_r4
c$OMP PARALLEL DO PRIVATE(i)
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=atm_field_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call inta2o(atm_field,tmp2d)
        endif
        call xcaput(tmp2d,clouds(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   clouds(1-nbdy,1-nbdy,k))
        call smtfld(clouds(1-nbdy,1-nbdy,k),smtmsk(1-nbdy,1-nbdy,k),
     .              atm_ice_nsmt,atm_ice_swgt)
c
c --- - convert range of cloudiness from 0-100 to 0-1
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            clouds(i,j,k)=clouds(i,j,k)*1.e-2
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(clouds, 1,12, nbdy,nbdy, halo_ps)
c
c --- ------------------------------------------------------------------
c --- read precipitation rate [kg/m^2/s]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write(lp,*) 'reading atm. climatological precipitation rate...'
        open (unit=nfu,file=atm_path(1:atm_path_len)
     .                      //'clim/prate_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) atm_field_r4
c$OMP PARALLEL DO PRIVATE(i)
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=atm_field_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call inta2o(atm_field,tmp2d)
        endif
        call xcaput(tmp2d,precip(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   precip(1-nbdy,1-nbdy,k))
        call smtfld(precip(1-nbdy,1-nbdy,k),smtmsk(1-nbdy,1-nbdy,k),
     .              atm_ice_nsmt,atm_ice_swgt)
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(precip, 1,12, nbdy,nbdy, halo_ps)
c
c --- ------------------------------------------------------------------
c --- read latent heat net flux [W/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write(lp,*) 'reading atm. climatological latent heat flux...'
        open (unit=nfu,file=atm_path(1:atm_path_len)
     .                      //'clim/lhtfl_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) atm_field_r4
c$OMP PARALLEL DO PRIVATE(i)
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=atm_field_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call inta2o(atm_field,tmp2d)
        endif
        call xcaput(tmp2d,lhtflx(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   lhtflx(1-nbdy,1-nbdy,k))
        call smtfld(lhtflx(1-nbdy,1-nbdy,k),smtmsk(1-nbdy,1-nbdy,k),
     .              atm_ice_nsmt,atm_ice_swgt)
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(lhtflx, 1,12, nbdy,nbdy, halo_ps)
c
c --- ------------------------------------------------------------------
c --- read sensible heat net flux [W/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write(lp,*) 'reading atm. climatological sensible heat flux...'
        open (unit=nfu,file=atm_path(1:atm_path_len)
     .                      //'clim/shtfl_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) atm_field_r4
c$OMP PARALLEL DO PRIVATE(i)
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=atm_field_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call inta2o(atm_field,tmp2d)
        endif
        call xcaput(tmp2d,shtflx(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   shtflx(1-nbdy,1-nbdy,k))
        call smtfld(shtflx(1-nbdy,1-nbdy,k),smtmsk(1-nbdy,1-nbdy,k),
     .              atm_ice_nsmt,atm_ice_swgt)
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(shtflx, 1,12, nbdy,nbdy, halo_ps)
c
c --- ------------------------------------------------------------------
c --- Read surface temperature [K]. Needed for surface pressure
c --- correction.
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write (lp,*)
     .    'reading atm. climatological surface temperature...'
        open (unit=nfu,file=atm_path(1:atm_path_len)
     .                      //'clim/skt_1968-1996.uf',
     .        form='unformatted')
c
        do k=1,12
          read (nfu) atm_field_r4
c$OMP PARALLEL DO PRIVATE(i)
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_sktclm(i,j,k)=atm_field_r4(i,j)+.0065*atm_topo(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
        enddo
c
        close (unit=nfu)
      endif
c
c --- ------------------------------------------------------------------
c --- read surface pressure [Pa]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write(lp,*) 'reading atm. climatological surface pressure...'
        open (unit=nfu,file=atm_path(1:atm_path_len)
     .                      //'clim/pres_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) atm_field_r4
c$OMP PARALLEL DO PRIVATE(i)
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=
     .          atm_field_r4(i,j)*exp(9.81*atm_topo(i,j)
     .          /(287.*(atm_sktclm(i,j,k)-.00325*atm_topo(i,j))))
            enddo
          enddo
c$OMP END PARALLEL DO
          call inta2o(atm_field,tmp2d)
        endif
        call xcaput(tmp2d,slpres(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   slpres(1-nbdy,1-nbdy,k))
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(slpres, 1,12, nbdy,nbdy, halo_ps)
c
c --- ------------------------------------------------------------------
c --- read momentum flux [N/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write(lp,*) 'reading atm. climatological momentum flux...'
        open (unit=nfu,file=atm_path(1:atm_path_len)
     .                      //'clim/momfl_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) atm_field_r4
c$OMP PARALLEL DO PRIVATE(i)
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=atm_field_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call inta2o(atm_field,tmp2d)
        endif
        call xcaput(tmp2d,taud(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   taud(1-nbdy,1-nbdy,k))
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(taud, 1,12, nbdy,nbdy, halo_ps)
c
c --- ------------------------------------------------------------------
c --- read momentum flux components [N/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write(lp,*)
     .    'reading atm. climatological u-component of momentum flux...'
        open (unit=nfu,file=atm_path(1:atm_path_len)
     .                      //'clim/uflx_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) atm_field_r4
c$OMP PARALLEL DO PRIVATE(i)
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=atm_field_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call inta2o(atm_field,tmp2d)
        endif
        call xcaput(tmp2d,tauxd(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_pv,
     .                   tauxd(1-nbdy,1-nbdy,k))
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      if (mnproc.eq.1) then
        write(lp,*)
     .    'reading atm. climatological v-component of momentum flux...'
        open (unit=nfu,file=atm_path(1:atm_path_len)
     .                      //'clim/vflx_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) atm_field_r4
c$OMP PARALLEL DO PRIVATE(i)
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=atm_field_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call inta2o(atm_field,tmp2d)
        endif
        call xcaput(tmp2d,tauyd(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_pv,
     .                   tauyd(1-nbdy,1-nbdy,k))
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      do k=1,12
c
c --- - change sign of momentum flux
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tauxd(i,j,k)=-tauxd(i,j,k)
            tauyd(i,j,k)=-tauyd(i,j,k)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
c --- - rotate the vector components
        call xctilr(tauxd(1-nbdy,1-nbdy,k), 1,1, 1,1, halo_pv)
        call xctilr(tauyd(1-nbdy,1-nbdy,k), 1,1, 1,1, halo_pv)
        call uvrotr2g(tauxd(1-nbdy,1-nbdy,k),tauyd(1-nbdy,1-nbdy,k))
c
      enddo
c
      call xctilr(tauxd, 1,12, nbdy,nbdy, halo_uv)
      call xctilr(tauyd, 1,12, nbdy,nbdy, halo_vv)
c
c --- ------------------------------------------------------------------
c --- read runoff rate [kg/m^2/s]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write(lp,*) 'reading atm. climatological runoff rate...'
c
c --- - read file containing ocean discarge points/weights for land
c --- - areas
c
        allocate(rnf_wgt  (atm_abdm,atm_idm,atm_jdm),
     .           rnf_ocdpi(atm_abdm,atm_idm,atm_jdm),
     .           rnf_ocdpj(atm_abdm,atm_idm,atm_jdm))
c
        open (unit=nfu,file='runoffweights.uf',
     .        form='unformatted',status='old',action='read')
        do j=1,atm_jdm
          do i=1,atm_idm
            read (nfu) rn2
            if (rn2.gt.atm_abdm) then
              write (lp,*) '''atm_abdm'' too small!'
              call xchalt('(inifrc_ben02clim)')
                     stop '(inifrc_ben02clim)'
            endif
            do k=1,rn2
              read (nfu) rw4,ri2,rj2
              rnf_wgt(k,i,j)=rw4
              rnf_ocdpi(k,i,j)=ri2
              rnf_ocdpj(k,i,j)=rj2
            enddo
            do k=rn2+1,atm_abdm
              rnf_wgt(k,i,j)=0.
              rnf_ocdpi(k,i,j)=0
              rnf_ocdpj(k,i,j)=0
            enddo
          enddo
        enddo
        close (unit=nfu)
c
        open (unit=nfu,file=atm_path(1:atm_path_len)
     .                      //'clim/runof_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) atm_field_r4
c$OMP PARALLEL DO PRIVATE(i)
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=atm_field_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
c
c --- --- place runoff at ocean discharge points
c$OMP PARALLEL DO PRIVATE(i)
          do j=1,jtdm
            do i=1,itdm
              tmp2d(i,j)=0.
            enddo
          enddo
c$OMP END PARALLEL DO
          do j=1,atm_jdm
            do i=1,atm_idm
              do l=1,atm_abdm
                if (rnf_wgt(l,i,j).gt.0.) then
                  il=rnf_ocdpi(l,i,j)
                  jl=rnf_ocdpj(l,i,j)
                  tmp2d(il,jl)=tmp2d(il,jl)
     .              +atm_field(i,j)*rnf_wgt(l,i,j)
                endif
              enddo
            enddo
          enddo
        endif
c
        call xcaput(tmp2d,runoff(1-nbdy,1-nbdy,k),1)
c
c --- - smooth the runoff field
        call smtfld(runoff(1-nbdy,1-nbdy,k),ip,
     .              atm_rnf_nsmt,atm_rnf_swgt)
c
c --- - convert unit of runoff from kg/m^2/day to kg/m^2/s
c$OMP PARALLEL DO PRIVATE(l,i)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            runoff(i,j,k)=runoff(i,j,k)/86400.
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
        deallocate(rnf_wgt,
     .             rnf_ocdpi,
     .             rnf_ocdpj)
      endif
c
      call xctilr(runoff, 1,12, nbdy,nbdy, halo_ps)
c
c --- Deallocate memory used for interpolation of surface fields.
      if (mnproc.eq.1) then
        deallocate(atm_sktclm,
     .             atm_field,
     .             atm_field_r4)
      endif
      call fnlzai
c
c --- If SSS restoring is requested, read climatological sea surface
c --- salinity.
      if (srxday.gt.0.) call rdcsss
c
c --- Compute global surface budgets of climatological heat and
c --- freshwater fluxes
c
      fwf_fac=1.e-3 ! conversion factor kg/m^2/s -> m/s
c
c --- sum up freshwater fluxes
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=0.
          util2(i,j)=0.
          util3(i,j)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(k,l,i,garea,le)
      do j=1,jj
        do k=1,12
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            garea=scp2(i,j)*1.e-4 ! [m^2]
c
c --- ----- freshwater fluxes [m/s]
            util1(i,j)=util1(i,j)+precip(i,j,k)*fwf_fac*garea
            le=(2.501-.00237*(sstclm(i,j,k)-273.15))*1.e6
            util2(i,j)=util2(i,j)-(lhtflx(i,j,k)/le)*1.e-3*garea
            util3(i,j)=util3(i,j)+runoff(i,j,k)*fwf_fac*garea
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcsum(prc_sum,util1,ip)
      call xcsum(eva_sum,util2,ip)
      call xcsum(rnf_sum,util3,ip)
c
      fac=1.e-6/12.
      prc_sum=prc_sum*fac
      eva_sum=eva_sum*fac
      rnf_sum=rnf_sum*fac
c
c --- sum up heat fluxes
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=0.
          util2(i,j)=0.
          util3(i,j)=0.
          util4(i,j)=0.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO

      do k=1,12
c
c --- - compute albedo
        dangle=8.*atan(1.)*real(30*k-15)/360.
        call albw_eval(dangle,plat,clouds(1-nbdy,1-nbdy,k),albw_d,albw)
c$OMP PARALLEL DO PRIVATE(l,i,area,albedo)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            garea=scp2(i,j)*1.e-4 ! [m^2]
c
c --- ----- heat fluxes
            albedo=albs_f*ricclm(i,j,k)+albw(i,j)*(1.-ricclm(i,j,k))
            util1(i,j)=util1(i,j)+dswrfl(i,j,k)*(1.-albedo)*garea
            util2(i,j)=util2(i,j)-nlwrfs(i,j,k)*garea
            util3(i,j)=util3(i,j)-lhtflx(i,j,k)*garea
            util4(i,j)=util4(i,j)-shtflx(i,j,k)*garea
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c
      call xcsum(swa_sum,util1,ip)
      call xcsum(lwa_sum,util2,ip)
      call xcsum(lht_sum,util3,ip)
      call xcsum(sht_sum,util4,ip)
c
      fac=1.e4/(12.*area)
      swa_ave=swa_sum*fac
      lwa_ave=lwa_sum*fac
      lht_ave=lht_sum*fac
      sht_ave=sht_sum*fac
c
      if (mnproc.eq.1) then
      write (lp,*)
      write (lp,*) 'Global precipitation:         ',prc_sum,'Sv'
      write (lp,*) 'Global evaporation:           ',eva_sum,'Sv'
      write (lp,*) 'Global runoff:                ',rnf_sum,'Sv'
      write (lp,*) 'Global balance of freshwater: ',
     .            prc_sum+eva_sum+rnf_sum,'Sv'
      write (lp,*)
      write (lp,*) 'Global mean short-wave radiation: ',swa_ave,'W/m^2'
      write (lp,*) 'Global mean long-wave radiation:  ',lwa_ave,'W/m^2'
      write (lp,*) 'Global mean latent heat flux:     ',lht_ave,'W/m^2'
      write (lp,*) 'Global mean sensible heat flux:   ',sht_ave,'W/m^2'
      write (lp,*) 'Global balance of mean heat-flux: ',
     .           swa_ave+lwa_ave+lht_ave+sht_ave,'W/m^2'
      write (lp,*)
      endif
c
c --- Adjust runoff and short-wave radiation to balance freshwater and
c --- heat fluxes, respectively.
      crnf=-(prc_sum+eva_sum)/rnf_sum
      cswa=-(lwa_ave+lht_ave+sht_ave)/swa_ave
c$OMP PARALLEL DO PRIVATE(k,l,i)
      do j=1-nbdy,jj+nbdy
        do k=1,12
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            runoff(i,j,k)=runoff(i,j,k)*crnf
            dswrfl(i,j,k)=dswrfl(i,j,k)*cswa
          enddo
          enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      if (mnproc.eq.1) then
      write (lp,*) 'Runoff has been adjusted by a factor',crnf
      write (lp,*) 'Short-wave radiation has been adjusted by a factor',
     .              cswa
      write (lp,*)
      endif
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'chwbud:'
        endif
        call chksummsk(runoff,ip,12,'runoff')
        call chksummsk(dswrfl,ip,12,'dswrfl')
      endif
c
c --- Allocate memory for transfer coefficients, gustiness squared, and
c --- air density
      allocate(cd_d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         ch_d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         ce_d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         wg2_d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         cd_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         ch_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         ce_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         wg2_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         rhoa(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
c
c --- Initialize transfer coefficients, gustiness squared, and air density
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          cd_d(i,j)=cd_r
          ch_d(i,j)=ch_r
          ce_d(i,j)=ce_r
          wg2_d(i,j)=wg2_r*1.e-4
          cd_m(i,j)=cd_r
          ch_m(i,j)=ch_r
          ce_m(i,j)=ce_r
          wg2_m(i,j)=wg2_r*1.e-4
          rhoa(i,j)=rhoa_r*1.e3
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Number of iteration in the computation of transfer coefficients.
      tciter=5
c
c --- Allocate memory for accumulation variables
      allocate(tsi_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         tml_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         sml_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         alb_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         fice_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         tsi(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'inifrc_ben02clim:'
        endif
        call chksummsk(dswrfl,ip,12,'dswrfl')
        call chksummsk(nlwrfs,ip,12,'nlwrfs')
        call chksummsk(clouds,ip,12,'clouds')
        call chksummsk(precip,ip,12,'precip')
        call chksummsk(lhtflx,ip,12,'lhtflx')
        call chksummsk(shtflx,ip,12,'shtflx')
        call chksummsk(slpres,ip,12,'slpres')
        call chksummsk(taud,ip,12,'taud')
        call chksummsk(tauxd,iu,12,'tauxd')
        call chksummsk(tauyd,iv,12,'tauyd')
        call chksummsk(runoff,ip,12,'runoff')
      endif
c
      return
      end subroutine inifrc_ben02clim
c
c --- ------------------------------------------------------------------
c
      subroutine inifrc_ben02syn
c
c --- ------------------------------------------------------------------
c --- Prepare utilization of daily atmospheric forcing
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_blocks.h"
#include "common_forc.h"
#include "common_clndr.h"
c
      real dx2,dy2
      real*4 rw4
      integer i,j,k,l
      integer*2 rn2,ri2,rj2
c
c --- Prepare interpolation of surface fields.
      call initai
c
c --- Read climatological sea-ice concentration.
      call rdcsic
c
c --- Read climatological surface temperature.
      call rdctsf
c
c --- If SSS restoring is requested, read climatological sea surface
c --- salinity.
      if (srxday.gt.0.) call rdcsss
c
c --- Initialize diagnosing/application of relaxation fluxes
      call idarlx
c
c --- Allocate memory for daily forcing fields.
      allocate(taud  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5),
     .         tauxd (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5),
     .         tauyd (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5),
     .         dswrfl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5),
     .         nlwrfs(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5),
     .         shtflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5),
     .         lhtflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5),
     .         precip(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5),
     .         clouds(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5),
     .         slpres(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5),
     .         runoff(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5),
     .         tmpsfc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5),
     .         ricec (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5))
c
c --- Compute smoothing weights atm_ice_swgt and atm_rnf_swgt. For
c --- stability atm_ice_swgt,atm_rnf_swgt < .5*dx^2*dy^2/(dx^2+dy^2).
      atm_ice_swgt=spval
      atm_rnf_swgt=spval
c$OMP PARALLEL DO PRIVATE(l,i,dx2,dy2)
c$OMP+REDUCTION(min:atm_ice_swgt,atm_rnf_swgt)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          dx2=scpx(i,j)*scpx(i,j)
          dy2=scpy(i,j)*scpy(i,j)
          atm_ice_swgt=min(atm_ice_swgt,.5*dx2*dy2/(dx2+dy2))
          atm_rnf_swgt=min(atm_rnf_swgt,.5*dx2*dy2/(dx2+dy2))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcmin(atm_ice_swgt)
      atm_ice_swgt=.9*atm_ice_swgt
      call xcmin(atm_rnf_swgt)
      atm_rnf_swgt=.9*atm_rnf_swgt
c
c --- Number of smoothing iterations is choosen to get a geographical
c --- extent of the smoothing independent of the grid resolution
      if (mnproc.eq.1) then
        if     (atm_idm.eq.atm_idm_ncep.and.
     .          atm_jdm.eq.atm_jdm_ncep) then
          atm_ice_nsmt=nint(atm_ice_csmt_ncep/atm_ice_swgt)
          atm_rnf_nsmt=nint(atm_rnf_csmt_ncep/atm_rnf_swgt)
        elseif (atm_idm.eq.atm_idm_era .and.
     .          atm_jdm.eq.atm_jdm_era ) then
          atm_ice_nsmt=nint(atm_ice_csmt_era/atm_ice_swgt)
          atm_rnf_nsmt=nint(atm_rnf_csmt_era/atm_rnf_swgt)
        endif
      endif
      call xcbcst(atm_ice_nsmt)
      call xcbcst(atm_rnf_nsmt)
      if (mnproc.eq.1) then
        write (lp,*) 'atm_ice_swgt',atm_ice_swgt,
     .               'atm_ice_nsmt',atm_ice_nsmt
        write (lp,*) 'atm_rnf_swgt',atm_rnf_swgt,
     .               'atm_rnf_nsmt',atm_rnf_nsmt
      endif
c
c --- set runoff and short-wave radiation adjustment factors
      if (mnproc.eq.1) then
        if     (atm_idm.eq.atm_idm_ncep.and.
     .          atm_jdm.eq.atm_jdm_ncep) then
          atm_crnf=atm_crnf_ncep
          atm_cswa=atm_cswa_ncep
        elseif (atm_idm.eq.atm_idm_era .and.
     .          atm_jdm.eq.atm_jdm_era ) then
          atm_crnf=atm_crnf_era
          atm_cswa=atm_cswa_era
        endif
      endif
      call xcbcst(atm_crnf)
      call xcbcst(atm_cswa)
c
      if (mnproc.eq.1) then
c
c --- - read file containing ocean discarge points/weights for land
c --- - areas
c
        allocate(rnf_wgt  (atm_abdm,atm_idm,atm_jdm),
     .           rnf_ocdpi(atm_abdm,atm_idm,atm_jdm),
     .           rnf_ocdpj(atm_abdm,atm_idm,atm_jdm))
c
        open (unit=nfu,file='runoffweights.uf',
     .        form='unformatted',status='old',action='read')
        do j=1,atm_jdm
          do i=1,atm_idm
            read (nfu) rn2
            if (rn2.gt.atm_abdm) then
              write (lp,*) '''atm_abdm'' too small!'
              call xchalt('(inifrc_ben02syn)')
                     stop '(inifrc_ben02syn)'
            endif
            do k=1,rn2
              read (nfu) rw4,ri2,rj2
              rnf_wgt(k,i,j)=rw4
              rnf_ocdpi(k,i,j)=ri2
              rnf_ocdpj(k,i,j)=rj2
            enddo
            do k=rn2+1,atm_abdm
              rnf_wgt(k,i,j)=0.
              rnf_ocdpi(k,i,j)=0
              rnf_ocdpj(k,i,j)=0
            enddo
          enddo
        enddo
        close (unit=nfu)
c
      endif
c
c --- Allocate memory for transfer coefficients, gustiness squared, and
c --- air density
      allocate(cd_d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         ch_d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         ce_d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         wg2_d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         cd_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         ch_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         ce_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         wg2_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         rhoa(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
c
c --- Initialize transfer coefficients, gustiness squared, and air density
c$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          cd_d(i,j)=cd_r
          ch_d(i,j)=ch_r
          ce_d(i,j)=ce_r
          wg2_d(i,j)=wg2_r*1.e-4
          cd_m(i,j)=cd_r
          ch_m(i,j)=ch_r
          ce_m(i,j)=ce_r
          wg2_m(i,j)=wg2_r*1.e-4
          rhoa(i,j)=rhoa_r*1.e3
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Number of iteration in the computation of transfer coefficients.
      tciter=1
c
c --- Allocate memory for accumulation variables
      allocate(tsi_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         tml_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         sml_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         alb_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         fice_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     .         tsi(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
c
c --- Initial indexes of the 5 time levels
      l1gi=1
      l2gi=2
      l3gi=3
      l4gi=4
      l5gi=5
c
c --- Go back 4 days so that all necessary fields are read
      call pstdat(nyear,nmonth,nday,4)
      call inicln
      do i=1,4
        call rdatm_syn
        call updcln(0)
      enddo
c
      return
      end subroutine inifrc_ben02syn
c
c --- ------------------------------------------------------------------
c
      subroutine getfrc_ben02clim
c
c --- ------------------------------------------------------------------
c --- Get climatological forcing
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_forc.h"
c
c --- Set interpolation parameters
      xgi=xmi
      l1gi=l1mi
      l2gi=l2mi
      l3gi=l3mi
      l4gi=l4mi
      l5gi=l5mi
c
c --- ------------------------------------------------------------------
c --- Compute the air-sea fluxes
c --- ------------------------------------------------------------------
c
      call asflux
c
      return
      end subroutine getfrc_ben02clim
c
c --- ------------------------------------------------------------------
c
      subroutine getfrc_ben02syn
c
c --- ------------------------------------------------------------------
c --- Get synoptic forcing
c --- ------------------------------------------------------------------
c
      implicit none
c
#include "common_blocks.h"
#include "common_clndr.h"
c
c --- Set interpolation parameter
      xgi=real(mod(nstep-1,nstep_in_day)+1)/nstep_in_day
c
c --- ------------------------------------------------------------------
c --- The first time step of a new day, read new forcing fields
c --- ------------------------------------------------------------------
c
      if (mod(nstep,nstep_in_day).eq.1) call rdatm_syn
c
c --- ------------------------------------------------------------------
c --- Compute the air-sea fluxes
c --- ------------------------------------------------------------------
c
      call asflux
c
      return
      end subroutine getfrc_ben02syn
c
c --- ------------------------------------------------------------------
c
      end module mod_ben02
