! ------------------------------------------------------------------------------
! Copyright (C) 2002-2022 Mats Bentsen, Mehmet Ilicak

! This file is part of BLOM.

! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.

! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.

! You should have received a copy of the GNU Lesser General Public License
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

module mod_ben02

  ! --- ------------------------------------------------------------------
  ! --- This module contains routines for generating surface fluxes
  ! --- from NCEP or ERA40 reanalysis following Bentsen and Drange (2002).
  ! --- ------------------------------------------------------------------

  use mod_types,       only: i2, r4
  use mod_config,      only: expcnf
  use mod_constants,   only: t0deg, spval, L_mks2cgs
  use mod_calendar,    only: date_offset, calendar_noerr, &
                             calendar_errstr
  use mod_time,        only: date, calendar, nday_in_year, nday_of_year, &
                             nstep, nstep_in_day, set_day_of_year, &
                             xmi, l1mi, l2mi, l3mi, l4mi, l5mi
  use mod_xc
  use mod_grid,        only: scpx, scpy, scuy, scvx, scp2, scp2i, scuxi, &
                             scvyi, plon, plat, area
  use mod_forcing,     only: sprfac, srxday, sstclm, ricclm, prfac, swa, &
                             nsf, lip, sop, eva, rfi, ztx, mty, ustarw, &
                             slp, abswnd
  use mod_utility,     only: util1, util2, util3, util4
  use mod_thdysi,      only: albs_f
  use mod_swtfrz,      only: swtfrz
  use mod_ben02func
  use mod_bulktf
  use mod_checksum,    only: csdiag, chksummsk
  use netcdf
  use mod_rdcsss,      only: rdcsss
  use mod_fill_global, only: fill_global
  use mod_intp1d,      only: intp1d
  use mod_idarlx,      only: idarlx

  implicit none
  private

  real, allocatable, dimension(:,:) :: &
       atm_lon, &            ! longitudes
       atm_lat, &            ! latitudes
       atm_topo              ! topography

  integer, allocatable, dimension(:,:) :: &
       atm_mask              ! mask

  real(r4), allocatable, dimension(:,:,:) :: &
       atm_wgt, &            ! interpolation weights
       rnf_wgt               ! interpolation weights

  integer(i2), allocatable, dimension(:,:,:) :: &
       atm_iwgt, &           ! interpolation adresses
       atm_jwgt, &           ! interpolation adresses
       rnf_ocdpi, &          ! interpolation adresses
       rnf_ocdpj             ! interpolation adresses

  real, allocatable, dimension(:,:,:) :: &
       taud, &               ! wind stress
       tauxd,tauyd, &        ! wind stress components
       dswrfl, &             ! short-wave surface radiation
       nlwrfs, &             ! net long-wave surface radiation
       shtflx, &             ! sensible heat flux
       lhtflx, &             ! latent heat flux
       precip, &             ! precipitation
       clouds, &             ! cloud cover
       slpres, &             ! sea level pressure
       runoff, &             ! runoff
       tmpsfc, &             ! surface temperature
       ricec                 ! ice concentration

  real, allocatable, dimension(:,:) :: &
       cd_d, &               ! transfer coefficient for momentum (data state)
       ch_d, &               ! transfer coefficient for sensible heat (data state)
       ce_d, &               ! transfer coefficient for latent heat (data state)
       wg2_d, &              ! gustiness squared (data state)
       cd_m, &               ! transfer coefficient for momentum (model state)
       ch_m, &               ! transfer coefficient for sensible heat (model state)
       ce_m, &               ! transfer coefficient for latent heat (model state)
       wg2_m, &              ! gustiness squared (model state)
       rhoa                  ! air density

  real, allocatable, dimension(:,:) :: &
       tsi_tda, &            ! accumulated snow/ice surface temperature
       tml_tda, &            ! accumulated mixed layer temperature
       sml_tda, &            ! accumulated mixed layer salinity
       alb_tda, &            ! accumulated albedo
       fice_tda, &           ! accumulated sea ice concentration
       tsi                   ! averaged snow/ice surface temperature

  integer :: &
       ntda                  ! accumulation number

  real, allocatable, dimension(:,:) :: &
       dfl,    &             ! derivative of non-solar heat flux by surface
                             ! temperature
       albw,   &             ! daily mean open water albedo.
       alb,    &             ! total surface albedo.
       rnfins, &             ! instantaneous runoff flux
       rnfres                ! runoff reservoar

  integer :: &
       nrfets = 7            ! e-folding time scale for detrainment of runoff
                             ! reservoar [days]

  integer, dimension(itdm,jtdm) :: &
       itp                   ! mask for the full domain

  real :: &
       atm_ice_swgt, &       ! smoothing weight for fields over ice
       atm_rnf_swgt, &       ! smoothing weight for fields over ice
       atm_crnf, &           ! runoff adjustment factor
       atm_cswa, &           ! short-wave radiation adjustment factor
       xgi                   ! interpolation parameter

  integer :: &
       atm_idm, &            ! zonal dimension of atmospheric grid
       atm_jdm, &            ! meridional dimension of atmospheric grid
       atm_path_len, &       ! length of path to daily atmospheric foring
       atm_ice_nsmt, &       ! number of smoothing iterations over ice
       atm_rnf_nsmt, &       ! number of smoothing iterations over ice
       l1gi,l2gi, &          ! time level indexes for interpolation
       l3gi,l4gi, &          ! time level indexes for interpolation
       l5gi                  ! time level indexes for interpolation

  character*80 :: &
       atm_path              ! path to daily atmospheric foring

  real, parameter :: &
       atm_mval = -9999., &  ! value of a point that should not have a
                             ! physical value &
       atm_fval = -99999.    ! value of a point that do not have a
                             ! physical value but should have one

  integer, parameter :: &
       atm_idm_ncep = 192, & ! zonal dimension of NCEP atmospheric grid &
       atm_jdm_ncep = 94, &  ! meridional dimension of NCEP atmospheric
                             ! grid &
       atm_idm_era = 320, &  ! zonal dimension of ERA atmospheric grid &
       atm_jdm_era = 160, &  ! meridional dimension of ERA atmospheric
                             ! grid &
       atm_nwgt = 12, &      ! number of neighbours used in the
                             ! interpolation &
       atm_abdm = 10         ! max. number of runoff discharge basins for
                             ! each atmospheric grid point

  real :: &
       atm_ice_csmt_ncep, &  ! constant determining how much the atm.
                             ! fields are smoothed over ice covered
                             ! regions (NCEP) &
       atm_rnf_csmt_ncep, &  ! constant determining how much the runoff is
                             ! smoothed at the coastal discharge points
                             ! (NCEP) &
       atm_crnf_ncep, &      ! runoff adjustment factor (NCEP) &
       atm_cswa_ncep, &      ! short-wave radiation adjustment factor
                             ! (NCEP) &
       atm_ice_csmt_era, &   ! constant determining how much the atm.
                             ! fields are smoothed over ice covered
                             ! regions (NCEP) &
       atm_rnf_csmt_era, &   ! constant determining how much the runoff is
                             ! smoothed at the coastal discharge points
                             ! (NCEP) &
       atm_crnf_era, &       ! runoff adjustment factor (NCEP) &
       atm_cswa_era          ! short-wave radiation adjustment factor
                             ! (NCEP)

#ifdef MKS
  data atm_ice_csmt_ncep,atm_rnf_csmt_ncep /2.e10,1.e9/, &
       atm_crnf_ncep,atm_cswa_ncep /0.82073,0.88340/, &
       atm_ice_csmt_era,atm_rnf_csmt_era /0.0,1.e9/, &
       atm_crnf_era,atm_cswa_era /0.7234,0.9721/
#else
  data atm_ice_csmt_ncep,atm_rnf_csmt_ncep /2.e14,1.e13/, &
       atm_crnf_ncep,atm_cswa_ncep /0.82073,0.88340/, &
       atm_ice_csmt_era,atm_rnf_csmt_era /0.0,1.e13/, &
       atm_crnf_era,atm_cswa_era /0.7234,0.9721/
#endif

  real :: &
       zu, &                ! measurement height of wind [m]
       zt, &                ! measurement height of temperature [m]
       zq, &                ! measurement height of specific humidity [m]
       emiss, &             ! emissivity of water []
       cpair, &             ! specific heat of dry air [J K-1 kg-1]
       stefanb, &           ! stefan-boltzman constant [W m-2 K-4]
       cd_r, &              ! reference transfer coefficient of momentum
       ch_r, &              ! reference transfer coefficient of sensible
                            ! heat
       ce_r, &              ! reference transfer coefficient of latent
                            ! heat
       rhoa_r, &            ! reference air density [g cm-3]
       wg2_r, &             ! reference gustiness squared [cm2 s-2]
       albdif, &            ! diffusive light albedo over water []
       rhowat               ! approximare density of sea-surface water
                            ! [kg m-3].

  data zu,zt,zq /10.,10.,10./, &
       emiss    /.97/, &
       cpair    /1004.7/, &
       stefanb  /5.67e-8/, &
       cd_r     /1.2e-3/, &
       ch_r     /1.1e-3/, &
       ce_r     /1.3e-3/, &
       wg2_r    /1.e3/, &
       rhoa_r   /1.225e-3/, &
       albdif   /.065/, &
       rhowat   /1024./

  integer :: &
       tciter             ! iterations in the computation of transfer
                          ! coefficients

  public :: atm_path,atm_path_len, &
            cd_d,ch_d,ce_d,wg2_d,cd_m,ch_m,ce_m,wg2_m,rhoa, &
            tsi_tda,tml_tda,sml_tda,alb_tda,fice_tda,tsi,ntda, &
            dfl,albw,alb,rnfins,rnfres,nrfets, &
            cd_r,ch_r,ce_r,wg2_r,rhoa_r,albdif,rhowat

  ! Public routines
  public :: inifrc_ben02clim
  public :: inifrc_ben02syn
  public :: getfrc_ben02clim
  public :: getfrc_ben02syn
  public :: initai
  public :: fnlzai
  public :: rdcsic
  public :: rdctsf

  ! private routines
  private :: rdatm_dim
  private :: rdatm_llm
  private :: rdatm_topo
  private :: rdatm_ts
  private :: rdatm_syn
  private :: inta2o
  private :: smtfld
  private :: uvrotr2g
  private :: albw_eval
  private :: asflux

contains

  ! --- ------------------------------------------------------------------

  subroutine rdatm_dim(filename)

    ! --- ------------------------------------------------------------------
    ! --- Read dimensions of atmospheric grid
    ! --- ------------------------------------------------------------------

    character :: filename*(*)

    integer :: status,ncid,dimid

    ! --- Open netcdf file.
    status = nf90_open(filename,nf90_nowrite,ncid)
    if (status /= nf90_noerr) then
      write(lp,'(4a)') 'nf90_open: ',trim(filename),': ', &
           nf90_strerror(status)
      call xchalt('(rdatm_dim)')
      stop '(rdatm_dim)'
    end if

    ! --- Get dimensions.
    status = nf90_inq_dimid(ncid,'lon',dimid)
    if (status /= nf90_noerr) then
      write(lp,'(2a)') 'nf90_inq_dimid: lon: ',nf90_strerror(status)
      call xchalt('(rdatm_dim)')
      stop '(rdatm_dim)'
    end if
    status=nf90_inquire_dimension(ncid,dimid,len = atm_idm)
    if (status /= nf90_noerr) then
      write(lp,'(2a)') 'nf90_inquire_dimension: lon: ', &
           nf90_strerror(status)
      call xchalt('(rdatm_dim)')
      stop '(rdatm_dim)'
    end if
    status = nf90_inq_dimid(ncid,'lat',dimid)
    if (status /= nf90_noerr) then
      write(lp,'(2a)') 'nf90_inq_dimid: lat: ',nf90_strerror(status)
      call xchalt('(rdatm_dim)')
      stop '(rdatm_dim)'
    end if
    status=nf90_inquire_dimension(ncid,dimid,len = atm_jdm)
    if (status /= nf90_noerr) then
      write(lp,'(2a)') 'nf90_inquire_dimension: lat: ', &
           nf90_strerror(status)
      call xchalt('(rdatm_dim)')
      stop '(rdatm_dim)'
    end if

    ! --- Close netcdf file.
    status = nf90_close(ncid)
    if (status /= nf90_noerr) then
      write(lp,'(4a)') 'nf90_close: ',trim(filename),': ', &
           nf90_strerror(status)
      call xchalt('(rdatm_dim)')
      stop '(rdatm_dim)'
    end if

  end subroutine rdatm_dim

  ! --- ------------------------------------------------------------------

  subroutine rdatm_llm(filename)

    ! --- ------------------------------------------------------------------
    ! --- Read atmospheric lon/lat coordinates and land mask from netcdf
    ! --- file
    ! --- ------------------------------------------------------------------

    character :: filename*(*)

    integer :: status,ncid,varid,i,j
    real*4 r4lon(atm_idm),r4lat(atm_jdm)
    integer*2 i2field(atm_idm,atm_jdm)

    ! --- Open netcdf file.
    status = nf90_open(filename,nf90_nowrite,ncid)
    if (status /= nf90_noerr) then
      write(lp,'(4a)') 'nf90_open: ',trim(filename),': ', &
           nf90_strerror(status)
      call xchalt('(rdatm_llm)')
      stop '(rdatm_llm)'
    end if

    ! --- Read longitudes.
    status = nf90_inq_varid(ncid,'lon',varid)
    if (status /= nf90_noerr) then
      write(lp,'(2a)') 'nf90_inq_varid: lon: ',nf90_strerror(status)
      call xchalt('(rdatm_llm)')
      stop '(rdatm_llm)'
    end if
    status = nf90_get_var(ncid,varid,r4lon)
    if (status /= nf90_noerr) then
      write(lp,'(2a)') 'nf90_get_var: lon: ',nf90_strerror(status)
      call xchalt('(rdatm_llm)')
      stop '(rdatm_llm)'
    end if

    ! --- Read latitudes.
    status = nf90_inq_varid(ncid,'lat',varid)
    if (status /= nf90_noerr) then
      write(lp,'(2a)') 'nf90_inq_varid: lat: ',nf90_strerror(status)
      call xchalt('(rdatm_llm)')
      stop '(rdatm_llm)'
    end if
    status = nf90_get_var(ncid,varid,r4lat)
    if (status /= nf90_noerr) then
      write(lp,'(2a)') 'nf90_get_var: lat: ',nf90_strerror(status)
      call xchalt('(rdatm_llm)')
      stop '(rdatm_llm)'
    end if

    ! --- Read land mask.
    status = nf90_inq_varid(ncid,'land',varid)
    if (status /= nf90_noerr) then
      write(lp,'(2a)') 'nf90_inq_varid: land: ',nf90_strerror(status)
      call xchalt('(rdatm_llm)')
      stop '(rdatm_llm)'
    end if
    status = nf90_get_var(ncid,varid,i2field)
    if (status /= nf90_noerr) then
      write(lp,'(2a)') 'nf90_get_var: land: ',nf90_strerror(status)
      call xchalt('(rdatm_llm)')
      stop '(rdatm_llm)'
    end if

    ! --- Close netcdf file.
    status = nf90_close(ncid)
    if (status /= nf90_noerr) then
      write(lp,'(4a)') 'nf90_close: ',trim(filename),': ', &
           nf90_strerror(status)
      call xchalt('(rdatm_llm)')
      stop '(rdatm_llm)'
    end if

    !$OMP PARALLEL DO PRIVATE(i)
    do j = 1,atm_jdm
      do i = 1,atm_idm
        atm_lon(i,j) = r4lon(i)
        atm_lat(i,j) = r4lat(j)
        atm_mask(i,j) = 1-i2field(i,j)
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine rdatm_llm

  ! --- ------------------------------------------------------------------

  subroutine rdatm_topo(filename)

    ! --- ------------------------------------------------------------------
    ! --- Read atmospheric lon/lat coordinates and land mask from netcdf
    ! --- file
    ! --- ------------------------------------------------------------------

    character :: filename*(*)

    integer :: status,ncid,varid,i,j
    real*4 r4field(atm_idm,atm_jdm)

    ! --- Open netcdf file.
    status = nf90_open(filename,nf90_nowrite,ncid)
    if (status /= nf90_noerr) then
      write(lp,'(4a)') 'nf90_open: ',trim(filename),': ', &
           nf90_strerror(status)
      call xchalt('(rdatm_topo)')
      stop '(rdatm_topo)'
    end if

    ! --- Read topography.
    status = nf90_inq_varid(ncid,'hgtsfc',varid)
    if (status /= nf90_noerr) then
      write(lp,'(2a)') 'nf90_inq_varid: hgtsfc: ', &
           nf90_strerror(status)
      call xchalt('(rdatm_topo)')
      stop '(rdatm_topo)'
    end if
    status = nf90_get_var(ncid,varid,r4field)
    if (status /= nf90_noerr) then
      write(lp,'(2a)') 'nf90_get_var: hgtsfc: ',nf90_strerror(status)
      call xchalt('(rdatm_topo)')
      stop '(rdatm_topo)'
    end if

    ! --- Close netcdf file.
    status = nf90_close(ncid)
    if (status /= nf90_noerr) then
      write(lp,'(4a)') 'nf90_close: ',trim(filename),': ', &
           nf90_strerror(status)
      call xchalt('(rdatm_topo)')
      stop '(rdatm_topo)'
    end if

    ! --- Convert field to real topography field
    !$OMP PARALLEL DO PRIVATE(i)
    do j = 1,atm_jdm
      do i = 1,atm_idm
        atm_topo(i,j) = r4field(i,atm_jdm+1-j)
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine rdatm_topo

  ! --- ------------------------------------------------------------------

  subroutine rdatm_ts(filename,fieldname,time_step,atm_field)

    ! --- ------------------------------------------------------------------
    ! --- Read field at selected time step from NetCDF file
    ! --- ------------------------------------------------------------------

    character :: filename*(*),fieldname*(*)
    integer :: time_step
    real :: atm_field(atm_idm,atm_jdm)

    integer :: status,ncid,dimid,atm_idm_t,atm_jdm_t,timeid,varid, &
         vartype,start(3),count(3),stride(3),i,j
    integer*2 i2field(atm_idm,atm_jdm),i2_mval
    real :: time
    real*4 offset,scale_factor

    ! --- open netcdf file
    status = nf90_open(filename,nf90_nowrite,ncid)
    if (status /= nf90_noerr) then
      write(lp,'(4a)') 'nf90_open: ',trim(filename),': ', &
           nf90_strerror(status)
      call xchalt('(rdatm_ts)')
      stop '(rdatm_ts)'
    end if

    ! --- check dimensions
    status = nf90_inq_dimid(ncid,'lon',dimid)
    if (status /= nf90_noerr) then
      write(lp,'(2a)') 'nf90_inq_dimid: lon: ',nf90_strerror(status)
      call xchalt('(rdatm_ts)')
      stop '(rdatm_ts)'
    end if
    status=nf90_inquire_dimension(ncid,dimid,len = atm_idm_t)
    if (status /= nf90_noerr) then
      write(lp,'(2a)') 'nf90_inquire_dimension: lon: ', &
           nf90_strerror(status)
      call xchalt('(rdatm_ts)')
      stop '(rdatm_ts)'
    end if
    status = nf90_inq_dimid(ncid,'lat',dimid)
    if (status /= nf90_noerr) then
      write(lp,'(2a)') 'nf90_inq_dimid: lat: ',nf90_strerror(status)
      call xchalt('(rdatm_ts)')
      stop '(rdatm_ts)'
    end if
    status=nf90_inquire_dimension(ncid,dimid,len = atm_jdm_t)
    if (status /= nf90_noerr) then
      write(lp,'(2a)') 'nf90_inquire_dimension: lat: ', &
           nf90_strerror(status)
      call xchalt('(rdatm_ts)')
      stop '(rdatm_ts)'
    end if
    if (atm_idm_t /= atm_idm.or.atm_jdm_t /= atm_jdm) then
      write(lp,'(a,i4,a,i4,a,i4,a,i4,3a)') &
           'expected dimensions ',atm_idm,' X ',atm_jdm,', but found ', &
           atm_idm_t,' X ',atm_jdm_t,' in file ',trim(filename),'!'
      call xchalt('(rdatm_ts)')
      stop '(rdatm_ts)'
    end if

    ! --- get id of time variable
    status = nf90_inq_varid(ncid,'time',timeid)
    if (status /= nf90_noerr) then
      write(lp,'(2a)') 'nf90_inq_varid: time: ',nf90_strerror(status)
      call xchalt('(rdatm_ts)')
      stop '(rdatm_ts)'
    end if

    ! --- read time variable at selected time level
    status = nf90_get_var(ncid,timeid,time,(/time_step/))
    if (status /= nf90_noerr) then
      write(lp,'(2a)') 'nf90_get_var: time: ',nf90_strerror(status)
      call xchalt('(rdatm_ts)')
      stop '(rdatm_ts)'
    end if
    if (time == 0.) then
      write(lp,'(2a)') 'rdatm_ts: ',trim(filename)
      write(lp,'(a)') &
           '  Time variable is zero at selected time level!'
      write(lp,'(a)') '  Probable cause is a corrupted file.'
      call xchalt('(rdatm_ts)')
      stop '(rdatm_ts)'
    end if

    ! --- get id of field variable
    status = nf90_inq_varid(ncid,fieldname,varid)
    if (status /= nf90_noerr) then
      write(lp,'(4a)') 'nf90_inq_varid: ',trim(fieldname),': ', &
           nf90_strerror(status)
      call xchalt('(rdatm_ts)')
      stop '(rdatm_ts)'
    end if

    ! --- get data type of field variable
    status=nf90_inquire_variable(ncid,varid,xtype = vartype)
    if (status /= nf90_noerr) then
      write(lp,'(4a)') 'nf90_inquire_variable: ',trim(fieldname), &
           ': ',nf90_strerror(status)
      call xchalt('(rdatm_ts)')
      stop '(rdatm_ts)'
    end if

    ! --- read selected time step of field

    start(1) = 1
    start(2) = 1
    start(3) = time_step
    count(1) = atm_idm
    count(2) = atm_jdm
    count(3) = 1
    stride(1) = 1
    stride(2) = 1
    stride(3) = 1

    if (vartype == nf90_short) then

      status = nf90_get_var(ncid,varid,i2field,start,count,stride)
      if (status /= nf90_noerr) then
        write(lp,'(4a)') 'nf90_get_var: ',trim(fieldname),': ', &
             nf90_strerror(status)
        call xchalt('(rdatm_ts)')
        stop '(rdatm_ts)'
      end if

      ! --- - read offset, scale factor, and value of no data
      status = nf90_get_att(ncid,varid,'add_offset',offset)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') 'nf90_get_att: add_offset: ', &
             nf90_strerror(status)
        call xchalt('(rdatm_ts)')
        stop '(rdatm_ts)'
      end if
      status = nf90_get_att(ncid,varid,'scale_factor',scale_factor)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') 'nf90_get_att: scale_factor: ', &
             nf90_strerror(status)
        call xchalt('(rdatm_ts)')
        stop '(rdatm_ts)'
      end if
      status = nf90_get_att(ncid,varid,'missing_value',i2_mval)
      if (status /= nf90_noerr) then
        write(lp,'(2a)') 'nf90_get_att: missing_value: ', &
             nf90_strerror(status)
        call xchalt('(rdatm_ts)')
        stop '(rdatm_ts)'
      end if

      ! --- - scale and add offset to field
      !$OMP PARALLEL DO PRIVATE(i)
      do j = 1,atm_jdm
        do i = 1,atm_idm
          if (i2field(i,j) == i2_mval) then
            atm_field(i,j) = atm_mval
          else
            atm_field(i,j) = i2field(i,j)*scale_factor+offset
          end if
        end do
      end do
      !$OMP END PARALLEL DO

    else

      status = nf90_get_var(ncid,varid,atm_field,start,count,stride)
      if (status /= nf90_noerr) then
        write(lp,'(4a)') 'nf90_get_var: ',trim(fieldname),': ', &
             nf90_strerror(status)
        call xchalt('(rdatm_ts)')
        stop '(rdatm_ts)'
      end if

    end if

    ! --- close netcdf file
    status = nf90_close(ncid)
    if (status /= nf90_noerr) then
      write(lp,'(4a)') 'nf90_close: ',trim(filename),': ', &
           nf90_strerror(status)
      call xchalt('(rdatm_ts)')
      stop '(rdatm_ts)'
    end if

  end subroutine rdatm_ts

  ! --- ------------------------------------------------------------------

  subroutine initai

    ! --- ------------------------------------------------------------------
    ! --- Prepare interpolation of atmospheric surface fields
    ! --- ------------------------------------------------------------------

    real, dimension(itdm,jtdm) :: tmp2da,tmp2db
    real, dimension(0:atm_nwgt) :: r_wgt
    real :: piloc,min_d,d,r,l2
    integer :: atm_idm_t,atm_jdm_t,idm_t,jdm_t,atm_nwgt_t,nw_2,is,js,it, &
               jt,iso,jso,m,n,i,j,iii,ism1,isp1,jsm1,jsp1,itm1,itp1,jtm1, &
               jtp1,nfu
    logical :: dimok,exists
    character :: filename*120

    ! --- Get global coordinates of p-points.
    call xcaget(tmp2da,plat,1)
    call xcaget(tmp2db,plon,1)

    if (mnproc == 1) then

      ! --- - Get dimensions of atmospheric grid
      filename = atm_path(1:atm_path_len)//'land.sfc.gauss.nc'
      call rdatm_dim(filename)
      if ((atm_idm /= atm_idm_ncep.or.atm_jdm /= atm_jdm_ncep).and. &
           (atm_idm /= atm_idm_era .or.atm_jdm /= atm_jdm_era )) then
        write (lp,*) 'initai: unknown atmospheric dimensions!'
        call xchalt('(iniai)')
        stop '(iniai)'
      end if
      if     (atm_idm == atm_idm_ncep.and. &
           atm_jdm == atm_jdm_ncep) then
        write (lp,*) 'initai: using NCEP reanalysis'
      else if (atm_idm == atm_idm_era .and. &
           atm_jdm == atm_jdm_era ) then
        write (lp,*) 'initai: using ERA40 reanalysis'
      end if

      ! --- - Allocate memory for the atmospheric interpolation.
      allocate(atm_lon(atm_idm,atm_jdm), &
               atm_lat(atm_idm,atm_jdm), &
               atm_mask(atm_idm,atm_jdm), &
               atm_topo(atm_idm,atm_jdm), &
               atm_wgt(atm_nwgt,itdm,jtdm), &
               atm_iwgt(atm_nwgt,itdm,jtdm), &
               atm_jwgt(atm_nwgt,itdm,jtdm))

      ! --- - Read atmospheric lon/lat coordinates and land mask.
      filename = atm_path(1:atm_path_len)//'land.sfc.gauss.nc'
      call rdatm_llm(filename)

      ! --- - Read atmospheric model topography.
      filename = atm_path(1:atm_path_len)//'hgt.sfc.nc'
      call rdatm_topo(filename)

      ! --- - Read interpolation weights if they exist, otherwise compute them
      dimok = .false.
      filename = 'atm_intwgt.uf'
      inquire(file=filename,exist = exists)
      if (exists) then
        open (newunit=nfu,file=filename,form = 'unformatted')
        read (nfu) atm_idm_t,atm_jdm_t,idm_t,jdm_t,atm_nwgt_t
        if (atm_idm_t == atm_idm.and.atm_jdm_t == atm_jdm.and. &
             idm_t == itdm.and.jdm_t == jtdm.and. &
             atm_nwgt_t == atm_nwgt) dimok = .true.
      end if
      if (exists.and.dimok) then
        read (nfu) atm_wgt,atm_iwgt,atm_jwgt
        close (unit = nfu)
      else
        if (exists) close (unit = nfu)

        write (lp,*) &
             'initai: computing atmospheric interpolation weights...'

        piloc = 4.*atan(1.)
        r_wgt(0) = -1.
        nw_2 = atm_nwgt/2+1

        is = 1
        js = 1

        do jt = 1,jtdm
          do it = 1,itdm

            min_d = spherdist(1.,atm_lon(is,js),atm_lat(is,js), &
                 tmp2db(it,jt),tmp2da(it,jt))

100         iso = is
            jso = js

            i = mod(iso-2+atm_idm,atm_idm)+1
            d = spherdist(1.,atm_lon(i,jso),atm_lat(i,jso), &
                 tmp2db(it,jt),tmp2da(it,jt))
            if (d < min_d) then
              is = i
              js = jso
              min_d = d
            end if
            i = mod(iso,atm_idm)+1
            d = spherdist(1.,atm_lon(i,jso),atm_lat(i,jso), &
                 tmp2db(it,jt),tmp2da(it,jt))
            if (d < min_d) then
              is = i
              js = jso
              min_d = d
            end if
            j = max(jso-1,1)
            d = spherdist(1.,atm_lon(iso,j),atm_lat(iso,j), &
                 tmp2db(it,jt),tmp2da(it,jt))
            if (d < min_d) then
              is = iso
              js = j
              min_d = d
            end if
            j = min(jso+1,atm_jdm)
            d = spherdist(1.,atm_lon(iso,j),atm_lat(iso,j), &
                 tmp2db(it,jt),tmp2da(it,jt))
            if (d < min_d) then
              is = iso
              js = j
              min_d = d
            end if

            if (is /= iso.or.js /= jso) goto 100

            do m = 1,atm_nwgt
              r_wgt(m) = 999999.
            end do

            do j = min(atm_jdm-nw_2*2,max(1,       js-nw_2)), &
                 min(atm_jdm       ,max(nw_2*2+1,js+nw_2))
              do iii = is-nw_2,is+nw_2
                i = mod(iii-1+atm_idm,atm_idm)+1
                r = spherdist(1.,atm_lon(i ,j ),atm_lat(i ,j ), &
                     tmp2db(it,jt),tmp2da(it,jt))
                m = atm_nwgt+1
10              m = m-1
                if (r < r_wgt(m)) goto 10
                m = m+1
                if (m <= atm_nwgt) then
                  do n = atm_nwgt-1,m,-1
                    r_wgt(n+1) = r_wgt(n)
                    atm_iwgt(n+1,it,jt) = atm_iwgt(n,it,jt)
                    atm_jwgt(n+1,it,jt) = atm_jwgt(n,it,jt)
                  end do
                  r_wgt(m) = r
                  atm_iwgt(m,it,jt) = int(i,i2)
                  atm_jwgt(m,it,jt) = int(j,i2)
                end if
              end do
            end do

            itm1 = max(min(itdm-2,max(1,it-1)),1)
            itp1 = min(itdm  ,max(3,it+1))
            jtm1 = mod(jt-2+jtdm,jtdm)+1
            jtp1 = mod(jt,jtdm)+1
            ism1 = mod(is-2+atm_idm,atm_idm)+1
            isp1 = mod(is,atm_idm)+1
            jsm1 = min(atm_jdm-2,max(1,js-1))
            jsp1 = min(atm_jdm  ,max(3,js+1))

            l2= &
                 .25*max( spherdist(1.,tmp2db(itm1,jt),tmp2da(itm1,jt), &
                 tmp2db(itp1,jt),tmp2da(itp1,jt)) &
                 *spherdist(1.,tmp2db(it,jtm1),tmp2da(it,jtm1), &
                 tmp2db(it,jtp1),tmp2da(it,jtp1)), &
                 spherdist(1.,atm_lon(ism1,js ),atm_lat(ism1,js ), &
                 atm_lon(isp1,js ),atm_lat(isp1,js )) &
                 *spherdist(1.,atm_lon(is ,jsm1),atm_lat(is ,jsm1), &
                 atm_lon(is ,jsp1),atm_lat(is ,jsp1))) &
                 /piloc

            do m = 1,atm_nwgt
              atm_wgt(m,it,jt)= &
                   real(max(exp(-.5*r_wgt(m)*r_wgt(m)/l2),1.e-9),r4)
            end do

          end do
        end do

        if (nreg == 2) then
          do it = 1,itdm
            do m = 1,atm_nwgt
              atm_wgt(m,it,jtdm) = atm_wgt(m,itdm-it+1,jtdm-1)
              atm_iwgt(m,it,jtdm) = atm_iwgt(m,itdm-it+1,jtdm-1)
              atm_jwgt(m,it,jtdm) = atm_jwgt(m,itdm-it+1,jtdm-1)
            end do
          end do
        end if

        open (newunit=nfu,file=filename,form = 'unformatted')
        write (nfu) atm_idm,atm_jdm,itdm,jtdm,atm_nwgt
        write (nfu) atm_wgt,atm_iwgt,atm_jwgt
        close (unit = nfu)

      end if

    end if

    ! --- Get mask for the full domain
    !$OMP PARALLEL DO PRIVATE(i)
    do j = 1,jj
      do i = 1,ii
        util1(i,j) = ip(i,j)
      end do
    end do
    !$OMP END PARALLEL DO
    call xcaget(tmp2da,util1,1)
    if (mnproc == 1) then
      !$OMP PARALLEL DO PRIVATE(i)
      do j = 1,jtdm
        do i = 1,itdm
          itp(i,j) = nint(tmp2da(i,j))
        end do
      end do
      !$OMP END PARALLEL DO
    end if

  end subroutine initai

  ! --- ------------------------------------------------------------------

  subroutine fnlzai

    ! --- ------------------------------------------------------------------
    ! --- Finalize interpolation of atmospheric surface fields (deallocate
    ! --- memory)
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      deallocate(atm_lon, &
                 atm_lat, &
                 atm_mask, &
                 atm_topo, &
                 atm_wgt, &
                 atm_iwgt, &
                 atm_jwgt)
    end if

  end subroutine fnlzai

  ! --- ------------------------------------------------------------------

  subroutine inta2o(adata,odata)

    ! --- ------------------------------------------------------------------
    ! --- Gaussian interpolation of 2D fields from atmospheric grid to ocean
    ! --- grid using precomputed weights.
    ! --- ------------------------------------------------------------------

    real, dimension(atm_idm,atm_jdm) :: adata
    real, dimension(itdm,jtdm) :: odata

    real :: w_sum,d_sum,w
    integer :: it,jt,n,is,js

    ! --- interpolate
    !$OMP PARALLEL DO PRIVATE(it,w_sum,d_sum,n,is,js,w)
    do jt = 1,jtdm
      do it = 1,itdm
        if (itp(it,jt) == 0) then
          odata(it,jt) = atm_mval
        else
          w_sum = 0.
          d_sum = 0.
          do n = 1,atm_nwgt
            is = atm_iwgt(n,it,jt)
            js = atm_jwgt(n,it,jt)
            if (atm_mask(is,js) == 1.and. &
                 adata(is,js) /= atm_mval) then
              w = atm_wgt(n,it,jt)
              w_sum = w_sum+w
              d_sum = d_sum+adata(is,js)*w
            end if
          end do
          if (w_sum == 0.) then
            odata(it,jt) = atm_fval
          else
            odata(it,jt) = d_sum/w_sum
          end if
        end if
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine inta2o

  ! --- ------------------------------------------------------------------

  subroutine smtfld(fld,msk,nsmt,swgt)

    ! --- ------------------------------------------------------------------
    ! --- Smooth the field, fld, a scale independent and conservative
    ! --- manner. The smoothing will be iterated nsmt times with a smoothing
    ! --- weight swgt and the smoothing is applied where the mask, msk, is
    ! --- equal one.
    ! --- ------------------------------------------------------------------

    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: fld
    integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: msk
    integer :: nsmt
    real :: swgt

    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: usmtfl,vsmtfl
    integer :: mbdy,n,i,j

    mbdy = 1

    do n = 1,nsmt

      if (mbdy == 1) then
        mbdy = nbdy
        call xctilr(fld,1,1,nbdy,nbdy,halo_ps)
      else
        mbdy = mbdy-1
      end if

      !$OMP PARALLEL DO PRIVATE(i)
      do j = 1-mbdy+1,jj+mbdy
        do i = 1-mbdy+1,ii+mbdy
          if (msk(i-1,j) == 1.and.msk(i,j) == 1) then
            usmtfl(i,j) = scuy(i,j)*scuxi(i,j)*(fld(i-1,j)-fld(i,j))
          else
            usmtfl(i,j) = 0.
          end if
          if (msk(i,j-1) == 1.and.msk(i,j) == 1) then
            vsmtfl(i,j) = scvx(i,j)*scvyi(i,j)*(fld(i,j-1)-fld(i,j))
          else
            vsmtfl(i,j) = 0.
          end if
        end do
      end do
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(i)
      do j = 1-mbdy+1,jj+mbdy-1
        do i = 1-mbdy+1,ii+mbdy-1
          if (msk(i,j) == 1) &
               fld(i,j) = fld(i,j) &
               -swgt*scp2i(i,j)*(usmtfl(i+1,j)-usmtfl(i,j) &
               +vsmtfl(i,j+1)-vsmtfl(i,j))
        end do
      end do
      !$OMP END PARALLEL DO

    end do

  end subroutine smtfld

  ! --- ------------------------------------------------------------------

  subroutine uvrotr2g(u,v)

    ! --- ------------------------------------------------------------------
    ! --- Rotate zonal/meridional vector components at C-grid scalar points
    ! --- to model grid components at C-grid velocity points
    ! --- ------------------------------------------------------------------

    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: u,v

    real :: latlim,rad
    parameter (latlim=87.,rad = 1.74532925199432958e-02)

    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: utmp,vtmp
    real :: mlat,dlat,dlon,psi
    integer :: i,j

    ! --- ------------------------------------------------------------------
    ! --- rotate vector components
    ! --- ------------------------------------------------------------------

    !$OMP PARALLEL DO PRIVATE(i,mlat,dlat,dlon,psi)
    do j = 1,jj
      do i = 1,ii

        ! --- --- u-component
        if (iu(i,j) == 0) then
          utmp(i,j) = atm_mval
        else
          mlat = .5*(plat(i,j)+plat(i-1,j))
          if (abs(mlat) > latlim) then
            utmp(i,j) = atm_fval
          else
            dlat = plat(i,j)-plat(i-1,j)
            dlon = mod(plon(i  ,j)+360.,360.) &
                 -mod(plon(i-1,j)+360.,360.)
            if (abs(dlon+360.) < abs(dlon)) dlon = dlon+360.
            if (abs(dlon-360.) < abs(dlon)) dlon = dlon-360.
            psi = atan2(dlat,cos(mlat*rad)*dlon)
            utmp(i,j) = .5*((u(i,j)+u(i-1,j))*cos(psi) &
                 +(v(i,j)+v(i-1,j))*sin(psi))
          end if
        end if

        ! --- --- v-component
        if (iv(i,j) == 0) then
          vtmp(i,j) = atm_mval
        else
          mlat = .5*(plat(i,j)+plat(i,j-1))
          if (abs(mlat) > latlim) then
            vtmp(i,j) = atm_fval
          else
            dlat = plat(i,j)-plat(i,j-1)
            dlon = mod(plon(i,j  )+360.,360.) &
                 -mod(plon(i,j-1)+360.,360.)
            if (abs(dlon+360.) < abs(dlon)) dlon = dlon+360.
            if (abs(dlon-360.) < abs(dlon)) dlon = dlon-360.
            psi = atan2(dlat,cos(mlat*rad)*dlon)
            vtmp(i,j) = .5*((u(i,j)+u(i,j-1))*cos(psi) &
                 +(v(i,j)+v(i,j-1))*sin(psi))
          end if
        end if

      end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i)
    do j = 1,jj
      do i = 1,ii
        u(i,j) = utmp(i,j)
        v(i,j) = vtmp(i,j)
      end do
    end do
    !$OMP END PARALLEL DO

    ! --- ------------------------------------------------------------------
    ! --- extrapolate values to velocity points near pole singularity
    ! --- ------------------------------------------------------------------

    call fill_global(atm_mval,atm_fval,halo_uv,u)
    call fill_global(atm_mval,atm_fval,halo_vv,v)

  end subroutine uvrotr2g

  ! --- ------------------------------------------------------------------

  subroutine albw_eval(dangle,lat,cc,albw_d,albw)

    ! --- ------------------------------------------------------------------
    ! --- compute 24 hrs mean albedo at the marine surface layer
    ! --- ------------------------------------------------------------------

    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: lat,cc,albw
    real :: dangle,albw_d

    real :: pi2,deg,eepsil,fraci,absh2o,s0,decli,sundv,sin2,cos2,stot, &
         sads,hangle,cosz,srad,sdir,sdif,altdeg,cfac,ssurf,albdir
    integer :: ifrac,i,j,l,npart

    ! --- set various quantities

    pi2 = 8.*atan(1.)          !        2 times pi
    deg = 360./pi2             !        convert from radians to degrees
    eepsil = 1.e-9             !        small number

    ifrac = 24                 !        split each 12 hrs day into ifrac parts
    fraci = 1./ifrac           !        1 over ifrac

    absh2o = 0.09              ! ---    absorption of water and ozone
    s0 = 1365.                 ! w/m^2  solar constant

    ! --- compute astronomic quantities

    decli = .006918+.070257*sin(dangle)   -.399912*cos(dangle) &
         +.000907*sin(2.*dangle)-.006758*cos(2.*dangle) &
         +.001480*sin(3.*dangle)-.002697*cos(3.*dangle)

    sundv = 1.00011+.001280*sin(dangle)   +.034221*cos(dangle) &
         +.000077*sin(2.*dangle)+.000719*cos(2.*dangle)

    !$omp parallel do private( &
    !$omp l,i,sin2,cos2,stot,sads,npart,hangle,cosz,srad,sdir,sdif,altdeg, &
    !$omp cfac,ssurf,albdir)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))

          ! --- --- compute astronomic quantities

          sin2 = sin(lat(i,j)/deg)*sin(decli)
          cos2 = cos(lat(i,j)/deg)*cos(decli)

          ! --- --- split each day into ifrac parts, and compute the solar
          ! --- --- radiance for each part. by assuming symmetry of the irradiance
          ! --- --- about noon, it is sufficient to compute the irradiance for the
          ! --- --- first 12 hrs of the (24 hrs) day (mean for the first 12 hrs
          ! --- --- equals then the mean for the last 12 hrs)

          stot = 0.
          sads = 0.

          do npart = 1,ifrac

            ! --- ----- hour angle, in radians
            hangle = pi2*(npart-.5)*fraci*.5

            ! --- ----- cosine of the zenith angle
            cosz = min(1.,max(0.,sin2+cos2*cos(hangle)))

            ! --- ----- extraterrestrial radiation
            srad =s0*sundv*cosz

            ! --- ----- direct radiation component
            sdir = srad*0.7**min(100.,1./(cosz+eepsil))

            ! --- ----- diffusive radiation component
            sdif = ((1.-absh2o)*srad-sdir)*.5

            ! --- ----- solar noon altitude in degrees
            altdeg = max(0.,asin(sin2+cos2))*deg

            ! --- ----- cloudiness correction
            cfac = (1.-0.62*cc(i,j)+0.0019*altdeg)

            ssurf = (sdir+sdif)*cfac+eepsil
            stot = stot+ssurf

            ! --- ----- albedo for direct light
            albdir = 0.03*exp(0.742*acos(cosz)**2.866)

            ! --- ----- radiation weighted sum of direct albedo
            sads = sads+albdir*ssurf

          end do

          ! --- --- daily mean albedo over water
          albw(i,j) = (1.-cc(i,j))*sads/stot+cc(i,j)*albw_d

        end do
      end do
    end do
    !$omp end parallel do

  end subroutine albw_eval

  ! --- ------------------------------------------------------------------

  subroutine rdatm_syn

    ! --- ------------------------------------------------------------------
    ! --- Read and interpolate daily atmospheric forcing fields. The
    ! --- interpolated fields of heat fluxes, cloud cover, precipitation,
    ! --- and surface temperature are smoothed over ice covered regions.
    ! --- ------------------------------------------------------------------

    real, dimension(atm_idm,atm_jdm) :: atm_field,atm_skt,atm_tau
    real, dimension(itdm,jtdm) :: tmp2d
    integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: smtmsk
    integer :: i,j,l,nday_of_year1,nyear1,iii,jjj
    character :: cyear*4,filename*120,fieldname*5

    ! --- The actual fields to be read are 2 days ahead
    nday_of_year1 = nday_of_year+2
    if (nday_of_year1 > nday_in_year) then
      nday_of_year1 = nday_of_year1-nday_in_year
      nyear1 = date%year+1
    else
      nyear1 = date%year
    end if

    ! --- Do not go beyond the start of the atm. time series
    if     (atm_idm == atm_idm_ncep.and.atm_jdm == atm_jdm_ncep) then
      if (nyear1 < 1948) then
        nday_of_year1 = 1
        nyear1 = 1948
      end if
    else if (atm_idm == atm_idm_era .and.atm_jdm == atm_jdm_era ) then
      if (nyear1 < 1958) then
        nday_of_year1 = 1
        nyear1 = 1958
      end if
    end if

    write(cyear,'(i4)') nyear1

    ! --- Rearrange the time level indexes
    i = l1gi
    l1gi = l2gi
    l2gi = l3gi
    l3gi = l4gi
    l4gi = l5gi
    l5gi = i

    if (mnproc == 1) then
      write (lp,'(a,i3,a,i4,a)') &
           'Reading and interpolating atm. forcing fields for day ', &
           nday_of_year1,' of year ',nyear1,'...'
      write (lp,'(a,5i2)') 'Time level indexes:', &
           l1gi,l2gi,l3gi,l4gi,l5gi
    end if

    ! --- ------------------------------------------------------------------
    ! --- Read/interpolate ice concentration [0-1]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      filename = atm_path(1:atm_path_len) &
           //'daily/ICECsfc/icec.sfc.gauss.'//cyear//'.nc'
      fieldname = 'icec'
      call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
      call inta2o(atm_field,tmp2d)
    end if
    call xcaput(tmp2d,ricec(1-nbdy,1-nbdy,l5gi),1)
    call fill_global(atm_mval,atm_fval,halo_ps, &
         ricec(1-nbdy,1-nbdy,l5gi))

    ! --- create smoothing mask - smooth where ice conc. is above 0.5
    !$omp parallel do private(i)
    do j = 1-nbdy,jj+nbdy
      do i = 1-nbdy,ii+nbdy
        if (ricec(i,j,l5gi) > .5.and.ip(i,j) == 1) then
          smtmsk(i,j) = 1
        else
          smtmsk(i,j) = 0
        end if
      end do
    end do
    !$omp end parallel do

    ! --- ------------------------------------------------------------------
    ! --- Read/interpolate short-wave radiation flux [W/m^2]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      filename = atm_path(1:atm_path_len) &
           //'daily/DSWRFsfc/dswrf.sfc.gauss.'//cyear//'.nc'
      fieldname = 'dswrf'
      call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
      call inta2o(atm_field,tmp2d)
    end if
    call xcaput(tmp2d,dswrfl(1-nbdy,1-nbdy,l5gi),1)
    call fill_global(atm_mval,atm_fval,halo_ps, &
         dswrfl(1-nbdy,1-nbdy,l5gi))
    call smtfld(dswrfl(1-nbdy,1-nbdy,l5gi),smtmsk, &
         atm_ice_nsmt,atm_ice_swgt)

    ! --- Adjust short-wave radiation field
    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          dswrfl(i,j,l5gi) = dswrfl(i,j,l5gi)*atm_cswa
        end do
      end do
    end do
    !$omp end parallel do

    ! --- ------------------------------------------------------------------
    ! --- Read/interpolate net long-wave radiation flux [W/m^2]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      filename = atm_path(1:atm_path_len) &
           //'daily/NLWRSsfc/nlwrs.sfc.gauss.'//cyear//'.nc'
      fieldname = 'nlwrs'
      call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
      call inta2o(atm_field,tmp2d)
    end if
    call xcaput(tmp2d,nlwrfs(1-nbdy,1-nbdy,l5gi),1)
    call fill_global(atm_mval,atm_fval,halo_ps, &
         nlwrfs(1-nbdy,1-nbdy,l5gi))
    call smtfld(nlwrfs(1-nbdy,1-nbdy,l5gi),smtmsk, &
         atm_ice_nsmt,atm_ice_swgt)

    ! --- ------------------------------------------------------------------
    ! --- Read/interpolate total cloud cover [0-100%]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      filename = atm_path(1:atm_path_len) &
           //'daily/TCDCeatm/tcdc.eatm.gauss.'//cyear//'.nc'
      fieldname = 'tcdc'
      call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
      call inta2o(atm_field,tmp2d)
    end if
    call xcaput(tmp2d,clouds(1-nbdy,1-nbdy,l5gi),1)
    call fill_global(atm_mval,atm_fval,halo_ps, &
         clouds(1-nbdy,1-nbdy,l5gi))
    call smtfld(clouds(1-nbdy,1-nbdy,l5gi),smtmsk, &
         atm_ice_nsmt,atm_ice_swgt)

    ! --- Convert range of cloudiness from 0-100 to 0-1
    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          clouds(i,j,l5gi) = clouds(i,j,l5gi)*1.e-2
        end do
      end do
    end do
    !$omp end parallel do

    ! --- ------------------------------------------------------------------
    ! --- Read/interpolate precipitation rate [kg/m^2/s]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      filename = atm_path(1:atm_path_len) &
           //'daily/PRATEsfc/prate.sfc.gauss.'//cyear//'.nc'
      fieldname = 'prate'
      call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
      call inta2o(atm_field,tmp2d)
    end if
    call xcaput(tmp2d,precip(1-nbdy,1-nbdy,l5gi),1)
    call fill_global(atm_mval,atm_fval,halo_ps, &
         precip(1-nbdy,1-nbdy,l5gi))
    call smtfld(precip(1-nbdy,1-nbdy,l5gi),smtmsk, &
         atm_ice_nsmt,atm_ice_swgt)

    ! --- ------------------------------------------------------------------
    ! --- Read/interpolate latent heat net flux [W/m^2]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      filename = atm_path(1:atm_path_len) &
           //'daily/LHTFLsfc/lhtfl.sfc.gauss.'//cyear//'.nc'
      fieldname = 'lhtfl'
      call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
      call inta2o(atm_field,tmp2d)
    end if
    call xcaput(tmp2d,lhtflx(1-nbdy,1-nbdy,l5gi),1)
    call fill_global(atm_mval,atm_fval,halo_ps, &
         lhtflx(1-nbdy,1-nbdy,l5gi))
    call smtfld(lhtflx(1-nbdy,1-nbdy,l5gi),smtmsk, &
         atm_ice_nsmt,atm_ice_swgt)

    ! --- ------------------------------------------------------------------
    ! --- Read/interpolate sensible heat net flux [W/m^2]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      filename = atm_path(1:atm_path_len) &
           //'daily/SHTFLsfc/shtfl.sfc.gauss.'//cyear//'.nc'
      fieldname = 'shtfl'
      call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
      call inta2o(atm_field,tmp2d)
    end if
    call xcaput(tmp2d,shtflx(1-nbdy,1-nbdy,l5gi),1)
    call fill_global(atm_mval,atm_fval,halo_ps, &
         shtflx(1-nbdy,1-nbdy,l5gi))
    call smtfld(shtflx(1-nbdy,1-nbdy,l5gi),smtmsk, &
         atm_ice_nsmt,atm_ice_swgt)

    ! --- ------------------------------------------------------------------
    ! --- Read/interpolate sea surface temperature [K]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      filename = atm_path(1:atm_path_len) &
           //'daily/SKTsfc/skt.sfc.gauss.'//cyear//'.nc'
      fieldname = 'skt'
      call rdatm_ts(filename,fieldname,nday_of_year1,atm_skt)

      ! --- - Compute sea level temperature
      !$omp parallel do private(i)
      do j = 1,atm_jdm
        do i = 1,atm_idm
          atm_skt(i,j) = atm_skt(i,j)+.0065*atm_topo(i,j)
        end do
      end do
      !$omp end parallel do

      call inta2o(atm_skt,tmp2d)
    end if
    call xcaput(tmp2d,tmpsfc(1-nbdy,1-nbdy,l5gi),1)
    call fill_global(atm_mval,atm_fval,halo_ps, &
         tmpsfc(1-nbdy,1-nbdy,l5gi))
    call smtfld(tmpsfc(1-nbdy,1-nbdy,l5gi),smtmsk, &
         atm_ice_nsmt,atm_ice_swgt)

    ! --- ------------------------------------------------------------------
    ! --- Read/interpolate surface pressure [Pa]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      filename = atm_path(1:atm_path_len) &
           //'daily/PRESsfc/pres.sfc.gauss.'//cyear//'.nc'
      fieldname = 'pres'
      call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)

      ! --- - compute sea level pressure
      !$omp parallel do private(i)
      do j = 1,atm_jdm
        do i = 1,atm_idm
          atm_field(i,j) = atm_field(i,j) &
               *exp(9.81*atm_topo(i,j) &
               /(287.*(atm_skt(i,j)-.00325*atm_topo(i,j))))
        end do
      end do
      !$omp end parallel do

      call inta2o(atm_field,tmp2d)
    end if
    call xcaput(tmp2d,slpres(1-nbdy,1-nbdy,l5gi),1)
    call fill_global(atm_mval,atm_fval,halo_ps, &
         slpres(1-nbdy,1-nbdy,l5gi))

    ! --- ------------------------------------------------------------------
    ! --- Read/interpolate momentum flux [N/m^2]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      filename = atm_path(1:atm_path_len) &
           //'daily/UFLXsfc/uflx.sfc.gauss.'//cyear//'.nc'
      fieldname = 'uflx'
      call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
      call inta2o(atm_field,tmp2d)
      !$omp parallel do private(i)
      do j = 1,atm_jdm
        do i = 1,atm_idm
          atm_tau(i,j) = atm_field(i,j)*atm_field(i,j)
        end do
      end do
      !$omp end parallel do
    end if
    call xcaput(tmp2d,tauxd(1-nbdy,1-nbdy,l5gi),1)
    call fill_global(atm_mval,atm_fval,halo_pv, &
         tauxd(1-nbdy,1-nbdy,l5gi))
    if (mnproc == 1) then
      filename = atm_path(1:atm_path_len) &
           //'daily/VFLXsfc/vflx.sfc.gauss.'//cyear//'.nc'
      fieldname = 'vflx'
      call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)
      call inta2o(atm_field,tmp2d)
      !$omp parallel do private(i)
      do j = 1,atm_jdm
        do i = 1,atm_idm
          atm_tau(i,j) = sqrt(atm_tau(i,j) &
               +atm_field(i,j)*atm_field(i,j))
        end do
      end do
      !$omp end parallel do
    end if
    call xcaput(tmp2d,tauyd(1-nbdy,1-nbdy,l5gi),1)
    call fill_global(atm_mval,atm_fval,halo_pv, &
         tauyd(1-nbdy,1-nbdy,l5gi))
    if (mnproc == 1) then
      call inta2o(atm_tau,tmp2d)
    end if
    call xcaput(tmp2d,taud(1-nbdy,1-nbdy,l5gi),1)
    call fill_global(atm_mval,atm_fval,halo_ps, &
         taud(1-nbdy,1-nbdy,l5gi))

    ! --- Change sign of momentum flux
    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          tauxd(i,j,l5gi) = -tauxd(i,j,l5gi)
          tauyd(i,j,l5gi) = -tauyd(i,j,l5gi)
        end do
      end do
    end do
    !$omp end parallel do

    ! --- rotate the vector components
    call xctilr(tauxd(1-nbdy,1-nbdy,l5gi), 1,1, 1,1, halo_pv)
    call xctilr(tauyd(1-nbdy,1-nbdy,l5gi), 1,1, 1,1, halo_pv)
    call uvrotr2g(tauxd(1-nbdy,1-nbdy,l5gi),tauyd(1-nbdy,1-nbdy,l5gi))

    ! --- ------------------------------------------------------------------
    ! --- Read/distribute runoff rate [kg/m^2/s]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      filename = atm_path(1:atm_path_len) &
           //'daily/RUNOFsfc/runof.sfc.gauss.'//cyear//'.nc'
      fieldname = 'runof'
      call rdatm_ts(filename,fieldname,nday_of_year1,atm_field)

      ! --- - Place runoff at ocean discharge points
      !$omp parallel do private(i)
      do j = 1,jtdm
        do i = 1,itdm
          tmp2d(i,j) = 0.
        end do
      end do
      !$omp end parallel do
      do j = 1,atm_jdm
        do i = 1,atm_idm
          do l = 1,atm_abdm
            if (rnf_wgt(l,i,j) > 0.) then
              iii = rnf_ocdpi(l,i,j)
              jjj = rnf_ocdpj(l,i,j)
              tmp2d(iii,jjj) = tmp2d(iii,jjj) &
                   +atm_field(i,j)*rnf_wgt(l,i,j)
            end if
          end do
        end do
      end do
    end if
    call xcaput(tmp2d,runoff(1-nbdy,1-nbdy,l5gi),1)

    ! --- Multiply runoff by 2 for fields prior to May 1957
    if ( nyear1 < 1957.or. &
         (nyear1 == 1957.and.nday_of_year1 <= 120)) then
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            runoff(i,j,l5gi) = runoff(i,j,l5gi)*2.
          end do
        end do
      end do
      !$omp end parallel do
    end if

    ! --- Adjust runoff field
    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          runoff(i,j,l5gi) = runoff(i,j,l5gi)*atm_crnf
        end do
      end do
    end do
    !$omp end parallel do

    ! --- Smooth the runoff field
    call smtfld(runoff(1-nbdy,1-nbdy,l5gi),ip, &
         atm_rnf_nsmt,atm_rnf_swgt)

    ! --- Convert unit of runoff from kg/m^2/day to kg/m^2/s
    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          runoff(i,j,l5gi) = runoff(i,j,l5gi)/86400.
        end do
      end do
    end do
    !$omp end parallel do

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'rdatm_syn:'
      end if
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
    end if

  end subroutine rdatm_syn

  ! --- ------------------------------------------------------------------

  subroutine asflux

    ! --- ------------------------------------------------------------------
    ! --- Compute air-sea fluxes. Same routine is used both for
    ! --- climatological and synoptic forcing
    ! --- ------------------------------------------------------------------

    ! --- Parameters:
    ! ---   dtmax - maximum near surface temperature gradient [K]
    ! ---   dqmax - maximum near surface specific humidity gradient [kg/kg]

    real :: dtmax,dqmax
    parameter (dtmax=30.,dqmax = .05)

    integer :: i,j,l,niter
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: taufac,cc
    real :: rntda,tau_d,prcp,dswrf,nlwrs,shtfl,lhtfl,slpr,tsrf_d,rice, &
         tml,sml,albedo,fice,tice_f,tml_d,tsi_d,qsrf_d,le,ua,sa,ta,qa, &
         tsrf_m,qsrf_m,dqsrf_m,dangle


    rntda = 1./real(ntda)
    ntda = 0

    !$omp parallel do private( &
    !$omp l,i,tau_d,prcp,dswrf,nlwrs,shtfl,lhtfl,slpr,tsrf_d,rice,tml,sml, &
    !$omp albedo,fice,tice_f,tml_d,tsi_d,qsrf_d,le,ua,sa,ta,qa,niter, &
    !$omp tsrf_m,qsrf_m,dqsrf_m) &
    !$omp shared(xgi,zu,zt,zq)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))

          ! --- ------------------------------------------------------------------
          ! --- --- Interpolate the climatological/synoptic fields
          ! --- ------------------------------------------------------------------

          tau_d =intp1d(taud  (i,j,l1gi),taud  (i,j,l2gi), &
               taud  (i,j,l3gi),taud  (i,j,l4gi), &
               taud  (i,j,l5gi),xgi)
          prcp  =intp1d(precip(i,j,l1gi),precip(i,j,l2gi), &
               precip(i,j,l3gi),precip(i,j,l4gi), &
               precip(i,j,l5gi),xgi)
          dswrf =intp1d(dswrfl(i,j,l1gi),dswrfl(i,j,l2gi), &
               dswrfl(i,j,l3gi),dswrfl(i,j,l4gi), &
               dswrfl(i,j,l5gi),xgi)
          nlwrs =intp1d(nlwrfs(i,j,l1gi),nlwrfs(i,j,l2gi), &
               nlwrfs(i,j,l3gi),nlwrfs(i,j,l4gi), &
               nlwrfs(i,j,l5gi),xgi)
          shtfl =intp1d(shtflx(i,j,l1gi),shtflx(i,j,l2gi), &
               shtflx(i,j,l3gi),shtflx(i,j,l4gi), &
               shtflx(i,j,l5gi),xgi)
          lhtfl =intp1d(lhtflx(i,j,l1gi),lhtflx(i,j,l2gi), &
               lhtflx(i,j,l3gi),lhtflx(i,j,l4gi), &
               lhtflx(i,j,l5gi),xgi)
          slpr  =intp1d(slpres(i,j,l1gi),slpres(i,j,l2gi), &
               slpres(i,j,l3gi),slpres(i,j,l4gi), &
               slpres(i,j,l5gi),xgi)
          if (expcnf == 'ben02clim') then
            tsrf_d = intp1d(sstclm(i,j,l1gi),sstclm(i,j,l2gi), &
                 sstclm(i,j,l3gi),sstclm(i,j,l4gi), &
                 sstclm(i,j,l5gi),xgi)
            rice  =intp1d(ricclm(i,j,l1gi),ricclm(i,j,l2gi), &
                 ricclm(i,j,l3gi),ricclm(i,j,l4gi), &
                 ricclm(i,j,l5gi),xgi)
          else if (expcnf == 'single_column') then
            tsrf_d = intp1d(sstclm(i,j,l1gi),sstclm(i,j,l2gi), &
                 sstclm(i,j,l3gi),sstclm(i,j,l4gi), &
                 sstclm(i,j,l5gi),xgi)
            rice  =intp1d(ricclm(i,j,l1gi),ricclm(i,j,l2gi), &
                 ricclm(i,j,l3gi),ricclm(i,j,l4gi), &
                 ricclm(i,j,l5gi),xgi)
          else
            tsrf_d = intp1d(tmpsfc(i,j,l1gi),tmpsfc(i,j,l2gi), &
                 tmpsfc(i,j,l3gi),tmpsfc(i,j,l4gi), &
                 tmpsfc(i,j,l5gi),xgi)
            rice  =intp1d(ricec (i,j,l1gi),ricec (i,j,l2gi), &
                 ricec (i,j,l3gi),ricec (i,j,l4gi), &
                 ricec (i,j,l5gi),xgi)
          end if
          rnfins(i,j) = intp1d(runoff(i,j,l1gi),runoff(i,j,l2gi), &
               runoff(i,j,l3gi),runoff(i,j,l4gi), &
               runoff(i,j,l5gi),xgi)

          prcp =max(0.,prcp)
          dswrf = max(0.,dswrf)
          rice =max(0.,min(1.,rice))

          ! --- ------------------------------------------------------------------
          ! --- --- Get averaged quantities obtained using the previous surface
          ! --- --- fluxes
          ! --- ------------------------------------------------------------------

          tsi(i,j) = tsi_tda(i,j)*rntda
          tml = tml_tda(i,j)*rntda
          sml = sml_tda(i,j)*rntda
          albedo = alb_tda(i,j)*rntda
          fice = fice_tda(i,j)*rntda

          tsi_tda(i,j) = 0.
          tml_tda(i,j) = 0.
          sml_tda(i,j) = 0.
          alb_tda(i,j) = 0.
          fice_tda(i,j) = 0.

          ! --- --- Freezing temperature of sea water
          tice_f = swtfrz(0.,sml)+t0deg

          ! --- ------------------------------------------------------------------
          ! --- --- Compute the atmospheric state by using the prescribed momentum
          ! --- --- and heat fluxes and the prescribed sea surface state
          ! --- ------------------------------------------------------------------

          tml_d = max(tsrf_d,tice_f)
          tsi_d = max(200.,(tsrf_d-(1.-rice)*tml_d)/max(rice,1.e-6))
          qsrf_d = rice*qsati(tsi_d,slpr)+(1.-rice)*qsatw(tml_d,slpr)
          le = (2.501-0.00237*(tsrf_d-273.15))*1.e6

          ! --- --- Make sure wind stress is not too small compared to latent and
          ! --- --- sensible heat fluxes
          sa = max(abs(shtfl)/(rhoa(i,j)*cpair*ch_d(i,j)*dtmax), &
               abs(lhtfl)/(rhoa(i,j)*le*ce_d(i,j)*dqmax))
          tau_d = max(tau_d,rhoa(i,j)*cd_d(i,j)*sa*sa)

          ! --- --- First guess on the atmospheric state by using the transfer
          ! --- --- coefficients and density from the previous time step
          ua = sqrt(.5*(-wg2_d(i,j) &
               +sqrt(wg2_d(i,j)*wg2_d(i,j) &
               +4.*(tau_d/(rhoa(i,j)*cd_d(i,j)))**2)))
          sa = sqrt(ua*ua+wg2_d(i,j))
          ta = tsrf_d-.0098*zt-shtfl/(rhoa(i,j)*cpair*ch_d(i,j)*sa)
          qa = qsrf_d-lhtfl/(rhoa(i,j)*le*ce_d(i,j)*sa)
          rhoa(i,j) = rhoair(ta,qa,slpr)

          ! --- --- Iteration loop for estimating transfer coefficients and
          ! --- --- atmospheric state
          do niter = 1,tciter

            ! --- ----- update the transfer coefficients and gustiness
            call bulktf(ua,zu,ta,zt,qa,zq,tsrf_d,qsrf_d,rice, &
                 cd_d(i,j),ch_d(i,j),ce_d(i,j),wg2_d(i,j))

            ! --- ----- update the atmospheric state
            ua = sqrt(.5*(-wg2_d(i,j) &
                 +sqrt(wg2_d(i,j)*wg2_d(i,j) &
                 +4.*(tau_d/(rhoa(i,j)*cd_d(i,j)))**2)))
            sa = sqrt(ua*ua+wg2_d(i,j))
            ta = tsrf_d-.0098*zt-shtfl/(rhoa(i,j)*cpair*ch_d(i,j)*sa)
            qa = qsrf_d-lhtfl/(rhoa(i,j)*le*ce_d(i,j)*sa)
            rhoa(i,j) = rhoair(ta,qa,slpr)

          end do

          ! --- ------------------------------------------------------------------
          ! --- --- Update transfer coefficients and gustiness with the computed
          ! --- --- atmospheric state the models ocean state
          ! --- ------------------------------------------------------------------

          tsrf_m = fice*tsi(i,j)+(1.-fice)*tml
          qsrf_m = fice*qsati(tsi(i,j),slpr)+(1.-fice)*qsatw(tml,slpr)

          do niter = 1,tciter
            call bulktf(ua,zu,ta,zt,qa,zq,tsrf_m,qsrf_m,fice, &
                 cd_m(i,j),ch_m(i,j),ce_m(i,j),wg2_m(i,j))
          end do

          ! --- ------------------------------------------------------------------
          ! --- --- Compute correction of the wind stress on the surface and wind
          ! --- --- generated tke [m/s]
          ! --- ------------------------------------------------------------------

          sa = sqrt(ua*ua+wg2_m(i,j))
          taufac(i,j) = rhoa(i,j)*cd_m(i,j)*sa*ua/tau_d

          ! --- --- Wind generated TKE
          ustarw(i,j) = sqrt(cd_m(i,j)*sa*ua*rhoa(i,j)/rhowat)

          ! --- ------------------------------------------------------------------
          ! --- --- Compute heat fluxes (positive downward) [W/m^2]
          ! --- ------------------------------------------------------------------

          swa(i,j) = dswrf*(1.-albedo)
          le = (2.501-.00237*(tsrf_m-273.15))*1.e6
          nsf(i,j) = rhoa(i,j)*cpair*ch_m(i,j)*sa*(ta+0.0098*zt-tsrf_m) &
               +rhoa(i,j)*ce_m(i,j)*le*sa*(qa-qsrf_m) &
               -nlwrs-4.*emiss*stefanb*ta**3*(tsrf_m-tsrf_d)

          ! --- ------------------------------------------------------------------
          ! --- --- Compute evaporation (positive downward) [kg/m^2/s]
          ! --- ------------------------------------------------------------------

          eva(i,j) = rhoa(i,j)*ce_m(i,j)*sa*(qa-qsrf_m)

          ! --- ------------------------------------------------------------------
          ! --- --- Compute derivative of non-solar flux by surface temperature
          ! --- --- [W/m^2/K]
          ! --- ------------------------------------------------------------------

          dqsrf_m = fice*dqsati(tsi(i,j),slpr)+(1.-fice)*dqsatw(tml,slpr)
          dfl(i,j) = -rhoa(i,j)*cpair*ch_m(i,j)*sa &
               -rhoa(i,j)*ce_m(i,j)*le*sa*dqsrf_m &
               -4.*emiss*stefanb*ta**3

          ! --- ------------------------------------------------------------------
          ! --- --- Split solid and liquid precipitation (positive downward)
          ! --- --- [kg/m^2/s]
          ! --- ------------------------------------------------------------------

          if (ta < t0deg) then
            lip(i,j) = 0.
            sop(i,j) = prcp
          else
            lip(i,j) = prcp
            sop(i,j) = 0.
          end if

          ! --- ------------------------------------------------------------------
          ! --- --- Sea level pressure [Pa] and wind speed at measurement height
          ! --- --- [m/s]
          ! --- ------------------------------------------------------------------

          slp(i,j) = slpr
          abswnd(i,j) = sa

          ! --- ------------------------------------------------------------------
          ! --- --- Frozen runoff is not computed.
          ! --- ------------------------------------------------------------------

          rfi(i,j) = 0.

          ! --- ------------------------------------------------------------------
          ! --- --- If requested, apply correction to precipitation and runoff to
          ! --- --- balance the fresh water budget
          ! --- ------------------------------------------------------------------

          if (sprfac) then
            lip(i,j) = lip(i,j)*prfac
            sop(i,j) = sop(i,j)*prfac
            rnfins(i,j) = rnfins(i,j)*prfac
            rfi(i,j) = rfi(i,j)*prfac
          end if

        end do
      end do
    end do
    !$omp end parallel do

    ! --- ------------------------------------------------------------------
    ! --- compute the surface wind stress [n/m^2]
    ! --- ------------------------------------------------------------------

    call xctilr(taufac,  1,1, 1,1, halo_ps)

    !$omp parallel do private(l,i) shared(xgi)
    do j = 1,jj
      do l = 1,isu(j)
        do i = max(1,ifu(j,l)),min(ii,ilu(j,l))
          ztx(i,j) = .5*(taufac(i,j)+taufac(i-1,j)) &
               *intp1d(tauxd(i,j,l1gi),tauxd(i,j,l2gi), &
               tauxd(i,j,l3gi),tauxd(i,j,l4gi), &
               tauxd(i,j,l5gi),xgi)
        end do
      end do
      do l = 1,isv(j)
        do i = max(1,ifv(j,l)),min(ii,ilv(j,l))
          mty(i,j) = .5*(taufac(i,j)+taufac(i,j-1)) &
               *intp1d(tauyd(i,j,l1gi),tauyd(i,j,l2gi), &
               tauyd(i,j,l3gi),tauyd(i,j,l4gi), &
               tauyd(i,j,l5gi),xgi)
        end do
      end do
    end do
    !$omp end parallel do

    ! --- ------------------------------------------------------------------
    ! --- compute open water albedo
    ! --- ------------------------------------------------------------------

    ! --- Time interpolation of cloud cover
    !$omp parallel do private(l,i) shared(xgi)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          cc(i,j) = intp1d(clouds(i,j,l1gi),clouds(i,j,l2gi), &
               clouds(i,j,l3gi),clouds(i,j,l4gi), &
               clouds(i,j,l5gi),xgi)
          cc(i,j) = max(0.,min(1.,cc(i,j)))
        end do
      end do
    end do
    !$omp end parallel do

    ! --- Obtain albedo for open water
    dangle = 8.*atan(1.)*real(nday_of_year-1)/real(nday_in_year)
    call albw_eval(dangle,plat,cc,albdif,albw)

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'asflux:'
      end if
      call chksummsk(swa,ip,1,'swa')
      call chksummsk(nsf,ip,1,'nsf')
      call chksummsk(dfl,ip,1,'dfl')
      call chksummsk(lip,ip,1,'lip')
      call chksummsk(sop,ip,1,'sop')
      call chksummsk(eva,ip,1,'eva')
      call chksummsk(ztx,iu,1,'ztx')
      call chksummsk(mty,iv,1,'mty')
      call chksummsk(rnfins,ip,1,'rnfins')
      call chksummsk(ustarw,ip,1,'ustarw')
      call chksummsk(tsi,ip,1,'tsi')
      call chksummsk(slp,ip,1,'slp')
      call chksummsk(abswnd,ip,1,'abswnd')
      call chksummsk(albw,ip,1,'albw')
    end if

  end subroutine asflux

  ! --- ------------------------------------------------------------------

  subroutine rdcsic

    ! --- ------------------------------------------------------------------
    ! --- Read and interpolate climatological sea-ice concentration
    ! --- ------------------------------------------------------------------

    real, dimension(atm_idm,atm_jdm) :: atm_field
    real*4, dimension(atm_idm,atm_jdm) :: atm_field_r4
    real, dimension(itdm,jtdm) :: tmp2d
    integer :: i,j,k,nfu

    if (mnproc == 1) then
      write (lp,*) 'reading atm. climatological ice concentration...'
      open (newunit=nfu,file = atm_path(1:atm_path_len) &
           //'clim/icec_1968-1996.uf', &
           form = 'unformatted')
    end if

    do k = 1,12
      if (mnproc == 1) then
        read(nfu) atm_field_r4
        !$omp parallel do private(i)
        do j = 1,atm_jdm
          do i = 1,atm_idm
            atm_field(i,j) = atm_field_r4(i,j)
          end do
        end do
        !$omp end parallel do
        call inta2o(atm_field,tmp2d)
      end if
      call xcaput(tmp2d,ricclm(1-nbdy,1-nbdy,k),1)
      call fill_global(atm_mval,atm_fval,halo_ps, &
           ricclm(1-nbdy,1-nbdy,k))
    end do

    if (mnproc == 1) then
      close (unit = nfu)
    end if

    call xctilr(ricclm, 1,12, nbdy,nbdy, halo_ps)

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'rdcsic:'
      end if
      call chksummsk(ricclm,ip,12,'ricclm')
    end if

  end subroutine rdcsic

  ! --- ------------------------------------------------------------------

  subroutine rdctsf

    ! --- ------------------------------------------------------------------
    ! --- Read and interpolate climatological surface temperature
    ! --- ------------------------------------------------------------------

    real, dimension(atm_idm,atm_jdm) :: atm_field
    real*4, dimension(atm_idm,atm_jdm) :: atm_field_r4
    real, dimension(itdm,jtdm) :: tmp2d
    real :: dx2,dy2
    integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: smtmsk
    integer :: i,j,k,l,nfu

    if (mnproc == 1) then
      write (lp,*) &
           'reading atm. climatological surface temperature...'
    end if

    ! --- compute smooting weight atm_ice_swgt. for stability
    ! --- atm_ice_swgt < .5*dx^2*dy^2/(dx^2+dy^2).
    atm_ice_swgt = spval
    !$omp parallel do private(l,i,dx2,dy2) reduction(min:atm_ice_swgt)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          dx2 = scpx(i,j)*scpx(i,j)
          dy2 = scpy(i,j)*scpy(i,j)
          atm_ice_swgt = min(atm_ice_swgt,.5*dx2*dy2/(dx2+dy2))
        end do
      end do
    end do
    !$omp end parallel do
    call xcmin(atm_ice_swgt)
    atm_ice_swgt = .9*atm_ice_swgt

    ! --- Number of smoothing iterations is choosen to get a geographical
    ! --- extent of the smoothing independent of the grid resolution
    if (mnproc == 1) then
      if     (atm_idm == atm_idm_ncep.and. &
           atm_jdm == atm_jdm_ncep) then
        atm_ice_nsmt = nint(atm_ice_csmt_ncep/atm_ice_swgt)
      else if (atm_idm == atm_idm_era .and. &
           atm_jdm == atm_jdm_era ) then
        atm_ice_nsmt = nint(atm_ice_csmt_era/atm_ice_swgt)
      end if
    end if
    call xcbcst(atm_ice_nsmt)
    if (mnproc == 1) then
      write (lp,*) 'atm_ice_swgt',atm_ice_swgt, &
           'atm_ice_nsmt',atm_ice_nsmt
    end if

    if (mnproc == 1) then
      open (newunit=nfu,file = atm_path(1:atm_path_len) &
           //'clim/skt_1968-1996.uf', &
           form = 'unformatted')
    end if

    do k = 1,12
      if (mnproc == 1) then
        read (nfu) atm_field_r4
        !$omp parallel do private(i)
        do j = 1,atm_jdm
          do i = 1,atm_idm
            atm_field(i,j) = atm_field_r4(i,j)+.0065*atm_topo(i,j)
          end do
        end do
        !$omp end parallel do
        call inta2o(atm_field,tmp2d)
      end if
      call xcaput(tmp2d,sstclm(1-nbdy,1-nbdy,k),1)
      call fill_global(atm_mval,atm_fval,halo_ps, &
           sstclm(1-nbdy,1-nbdy,k))

      ! --- - create smoothing mask - smooth where ice conc. is above 0.5
      !$omp parallel do private(i)
      do j = 1-nbdy,jj+nbdy
        do i = 1-nbdy,ii+nbdy
          if (ricclm(i,j,k) > .5.and.ip(i,j) == 1) then
            smtmsk(i,j) = 1
          else
            smtmsk(i,j) = 0
          end if
        end do
      end do
      !$omp end parallel do

      call smtfld(sstclm(1-nbdy,1-nbdy,k),smtmsk, &
           atm_ice_nsmt,atm_ice_swgt)

    end do

    if (mnproc == 1) then
      close (unit = nfu)
    end if

    call xctilr(sstclm, 1,12, nbdy,nbdy, halo_ps)

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'rdctsf:'
      end if
      call chksummsk(sstclm,ip,12,'sstclm')
    end if

  end subroutine rdctsf

  ! --- ------------------------------------------------------------------

  subroutine inifrc_ben02clim

    ! --- ------------------------------------------------------------------
    ! --- Initialize monthly climatological forcing fields
    ! --- ------------------------------------------------------------------

    real, allocatable, dimension(:,:,:) :: atm_sktclm
    real, allocatable, dimension(:,:) :: atm_field
    real*4, allocatable, dimension(:,:) :: atm_field_r4
    real, dimension(itdm,jtdm) :: tmp2d
    integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12) :: smtmsk
    real :: dx2,dy2,prc_sum,eva_sum,rnf_sum,swa_sum,lwa_sum,lht_sum, &
         sht_sum,fwf_fac,dangle,garea,le,albedo,fac,swa_ave,lwa_ave, &
         lht_ave,sht_ave,crnf,cswa,A_cgs2mks
    real*4 :: rw4
    integer :: i,j,k,l,il,jl,nfu
    integer*2 :: rn2,ri2,rj2

    A_cgs2mks = 1./(L_mks2cgs**2)

    ! --- Allocate memory for additional monthly forcing fields.
    allocate(taud  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12), &
             tauxd (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12), &
             tauyd (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12), &
             dswrfl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12), &
             nlwrfs(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12), &
             shtflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12), &
             lhtflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12), &
             precip(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12), &
             clouds(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12), &
             slpres(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12), &
             runoff(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12))

    ! --- Allocate memory for transfer coefficients, gustiness squared, and
    ! --- air density
    allocate(cd_d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             ch_d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             ce_d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             wg2_d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             cd_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             ch_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             ce_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             wg2_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             rhoa(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))

    ! --- Allocate memory for accumulation variables
    allocate(tsi_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             tml_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             sml_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             alb_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             fice_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             tsi(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))

    ! --- Allocate memory for derivative of non-solar heat flux by surface
    ! --- temperature, albedos and instantaneous runoff flux and runoff
    ! --- reservoar
    allocate(dfl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             albw(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             alb(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             rnfins(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             rnfres(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))

    ! --- number of iteration in the computation of transfer coefficients.
    tciter = 5

    ! --- Prepare interpolation of surface fields.
    call initai

    ! --- Read climatological sea-ice concentration.
    call rdcsic

    ! --- Read climatological surface temperature.
    call rdctsf

    ! --- Initialize diagnosing/application of relaxation fluxes
    call idarlx

    ! --- Allocate memory for temporary fields
    if (mnproc == 1) then
      allocate(atm_sktclm(atm_idm,atm_jdm,12), &
               atm_field(atm_idm,atm_jdm), &
               atm_field_r4(atm_idm,atm_jdm))
    end if

    ! --- Compute smoothing weights atm_ice_swgt and atm_rnf_swgt. For
    ! --- stability atm_ice_swgt,atm_rnf_swgt < .5*dx^2*dy^2/(dx^2+dy^2).
    atm_ice_swgt = spval
    atm_rnf_swgt = spval

    !$omp parallel do private(l,i,dx2,dy2) &
    !$omp reduction(min:atm_ice_swgt,atm_rnf_swgt)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          dx2 = scpx(i,j)*scpx(i,j)
          dy2 = scpy(i,j)*scpy(i,j)
          atm_ice_swgt = min(atm_ice_swgt,.5*dx2*dy2/(dx2+dy2))
          atm_rnf_swgt = min(atm_rnf_swgt,.5*dx2*dy2/(dx2+dy2))
        end do
      end do
    end do
    !$omp end parallel do

    call xcmin(atm_ice_swgt)
    atm_ice_swgt = .9*atm_ice_swgt
    call xcmin(atm_rnf_swgt)
    atm_rnf_swgt = .9*atm_rnf_swgt

    ! --- Number of smoothing iterations is choosen to get a geographical
    ! --- extent of the smoothing independent of the grid resolution
    if (mnproc == 1) then
      if     (atm_idm == atm_idm_ncep.and. &
           atm_jdm == atm_jdm_ncep) then
        atm_ice_nsmt = nint(atm_ice_csmt_ncep/atm_ice_swgt)
        atm_rnf_nsmt = nint(atm_rnf_csmt_ncep/atm_rnf_swgt)
      else if (atm_idm == atm_idm_era .and. &
           atm_jdm == atm_jdm_era ) then
        atm_ice_nsmt = nint(atm_ice_csmt_era/atm_ice_swgt)
        atm_rnf_nsmt = nint(atm_rnf_csmt_era/atm_rnf_swgt)
      end if
    end if
    call xcbcst(atm_ice_nsmt)
    call xcbcst(atm_rnf_nsmt)
    if (mnproc == 1) then
      write (lp,*) 'atm_ice_swgt',atm_ice_swgt, &
           'atm_ice_nsmt',atm_ice_nsmt
      write (lp,*) 'atm_rnf_swgt',atm_rnf_swgt, &
           'atm_rnf_nsmt',atm_rnf_nsmt
    end if

    ! --- create smoothing mask - smooth where ice conc. is above 0.5
    do k = 1,12
      !$omp parallel do private(i)
      do j = 1-nbdy,jj+nbdy
        do i = 1-nbdy,ii+nbdy
          if (ricclm(i,j,k) > .5.and.ip(i,j) == 1) then
            smtmsk(i,j,k) = 1
          else
            smtmsk(i,j,k) = 0
          end if
        end do
      end do
      !$omp end parallel do
    end do

    ! --- ------------------------------------------------------------------
    ! --- read short-wave radiation flux [W/m^2]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      write (lp,*) &
           'reading atm. climatological short-wave radiation flux...'
      open (newunit=nfu,file = atm_path(1:atm_path_len) &
           //'clim/dswrf_1968-1996.uf', &
           form = 'unformatted')
    end if

    do k = 1,12
      if (mnproc == 1) then
        read (nfu) atm_field_r4
        !$omp parallel do private(i)
        do j = 1,atm_jdm
          do i = 1,atm_idm
            atm_field(i,j) = atm_field_r4(i,j)
          end do
        end do
        !$omp end parallel do
        call inta2o(atm_field,tmp2d)
      end if
      call xcaput(tmp2d,dswrfl(1-nbdy,1-nbdy,k),1)
      call fill_global(atm_mval,atm_fval,halo_ps, &
           dswrfl(1-nbdy,1-nbdy,k))
      call smtfld(dswrfl(1-nbdy,1-nbdy,k),smtmsk(1-nbdy,1-nbdy,k), &
           atm_ice_nsmt,atm_ice_swgt)
    end do

    if (mnproc == 1) then
      close (unit = nfu)
    end if

    call xctilr(dswrfl, 1,12, nbdy,nbdy, halo_ps)

    ! --- ------------------------------------------------------------------
    ! --- read net long-wave radiation flux [W/m^2]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      write (lp,*) &
           'reading atm. climatological long-wave radiation flux...'
      open (newunit=nfu,file = atm_path(1:atm_path_len) &
           //'clim/nlwrs_1968-1996.uf', &
           form = 'unformatted')
    end if

    do k = 1,12
      if (mnproc == 1) then
        read (nfu) atm_field_r4
        !$omp parallel do private(i)
        do j = 1,atm_jdm
          do i = 1,atm_idm
            atm_field(i,j) = atm_field_r4(i,j)
          end do
        end do
        !$omp end parallel do
        call inta2o(atm_field,tmp2d)
      end if
      call xcaput(tmp2d,nlwrfs(1-nbdy,1-nbdy,k),1)
      call fill_global(atm_mval,atm_fval,halo_ps, &
           nlwrfs(1-nbdy,1-nbdy,k))
      call smtfld(nlwrfs(1-nbdy,1-nbdy,k),smtmsk(1-nbdy,1-nbdy,k), &
           atm_ice_nsmt,atm_ice_swgt)
    end do

    if (mnproc == 1) then
      close (unit = nfu)
    end if

    call xctilr(nlwrfs, 1,12, nbdy,nbdy, halo_ps)

    ! --- ------------------------------------------------------------------
    ! --- read total cloud cover [0-100%]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      write(lp,*) 'reading atm. climatological total cloud cover...'
      open (newunit=nfu,file = atm_path(1:atm_path_len) &
           //'clim/tcdc_1968-1996.uf', &
           form = 'unformatted')
    end if

    do k = 1,12
      if (mnproc == 1) then
        read (nfu) atm_field_r4
        !$omp parallel do private(i)
        do j = 1,atm_jdm
          do i = 1,atm_idm
            atm_field(i,j) = atm_field_r4(i,j)
          end do
        end do
        !$omp end parallel do
        call inta2o(atm_field,tmp2d)
      end if
      call xcaput(tmp2d,clouds(1-nbdy,1-nbdy,k),1)
      call fill_global(atm_mval,atm_fval,halo_ps, &
           clouds(1-nbdy,1-nbdy,k))
      call smtfld(clouds(1-nbdy,1-nbdy,k),smtmsk(1-nbdy,1-nbdy,k), &
           atm_ice_nsmt,atm_ice_swgt)

      ! --- - convert range of cloudiness from 0-100 to 0-1
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            clouds(i,j,k) = clouds(i,j,k)*1.e-2
          end do
        end do
      end do
      !$omp end parallel do
    end do

    if (mnproc == 1) then
      close (unit = nfu)
    end if

    call xctilr(clouds, 1,12, nbdy,nbdy, halo_ps)

    ! --- ------------------------------------------------------------------
    ! --- read precipitation rate [kg/m^2/s]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      write(lp,*) 'reading atm. climatological precipitation rate...'
      open (newunit=nfu,file = atm_path(1:atm_path_len) &
           //'clim/prate_1968-1996.uf', &
           form = 'unformatted')
    end if

    do k = 1,12
      if (mnproc == 1) then
        read (nfu) atm_field_r4
        !$omp parallel do private(i)
        do j = 1,atm_jdm
          do i = 1,atm_idm
            atm_field(i,j) = atm_field_r4(i,j)
          end do
        end do
        !$omp end parallel do
        call inta2o(atm_field,tmp2d)
      end if
      call xcaput(tmp2d,precip(1-nbdy,1-nbdy,k),1)
      call fill_global(atm_mval,atm_fval,halo_ps, &
           precip(1-nbdy,1-nbdy,k))
      call smtfld(precip(1-nbdy,1-nbdy,k),smtmsk(1-nbdy,1-nbdy,k), &
           atm_ice_nsmt,atm_ice_swgt)
    end do

    if (mnproc == 1) then
      close (unit = nfu)
    end if

    call xctilr(precip, 1,12, nbdy,nbdy, halo_ps)

    ! --- ------------------------------------------------------------------
    ! --- read latent heat net flux [W/m^2]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      write(lp,*) 'reading atm. climatological latent heat flux...'
      open (newunit=nfu,file = atm_path(1:atm_path_len) &
           //'clim/lhtfl_1968-1996.uf', &
           form = 'unformatted')
    end if

    do k = 1,12
      if (mnproc == 1) then
        read (nfu) atm_field_r4
        !$omp parallel do private(i)
        do j = 1,atm_jdm
          do i = 1,atm_idm
            atm_field(i,j) = atm_field_r4(i,j)
          end do
        end do
        !$omp end parallel do
        call inta2o(atm_field,tmp2d)
      end if
      call xcaput(tmp2d,lhtflx(1-nbdy,1-nbdy,k),1)
      call fill_global(atm_mval,atm_fval,halo_ps, &
           lhtflx(1-nbdy,1-nbdy,k))
      call smtfld(lhtflx(1-nbdy,1-nbdy,k),smtmsk(1-nbdy,1-nbdy,k), &
           atm_ice_nsmt,atm_ice_swgt)
    end do

    if (mnproc == 1) then
      close (unit = nfu)
    end if

    call xctilr(lhtflx, 1,12, nbdy,nbdy, halo_ps)

    ! --- ------------------------------------------------------------------
    ! --- read sensible heat net flux [W/m^2]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      write(lp,*) 'reading atm. climatological sensible heat flux...'
      open (newunit=nfu,file = atm_path(1:atm_path_len) &
           //'clim/shtfl_1968-1996.uf', &
           form = 'unformatted')
    end if

    do k = 1,12
      if (mnproc == 1) then
        read (nfu) atm_field_r4
        !$omp parallel do private(i)
        do j = 1,atm_jdm
          do i = 1,atm_idm
            atm_field(i,j) = atm_field_r4(i,j)
          end do
        end do
        !$omp end parallel do
        call inta2o(atm_field,tmp2d)
      end if
      call xcaput(tmp2d,shtflx(1-nbdy,1-nbdy,k),1)
      call fill_global(atm_mval,atm_fval,halo_ps, &
           shtflx(1-nbdy,1-nbdy,k))
      call smtfld(shtflx(1-nbdy,1-nbdy,k),smtmsk(1-nbdy,1-nbdy,k), &
           atm_ice_nsmt,atm_ice_swgt)
    end do

    if (mnproc == 1) then
      close (unit = nfu)
    end if

    call xctilr(shtflx, 1,12, nbdy,nbdy, halo_ps)

    ! --- ------------------------------------------------------------------
    ! --- Read surface temperature [K]. Needed for surface pressure
    ! --- correction.
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      write (lp,*) &
           'reading atm. climatological surface temperature...'
      open (newunit=nfu,file = atm_path(1:atm_path_len) &
           //'clim/skt_1968-1996.uf', &
           form = 'unformatted')

      do k = 1,12
        read (nfu) atm_field_r4
        !$omp parallel do private(i)
        do j = 1,atm_jdm
          do i = 1,atm_idm
            atm_sktclm(i,j,k) = atm_field_r4(i,j)+.0065*atm_topo(i,j)
          end do
        end do
        !$omp end parallel do
      end do

      close (unit = nfu)
    end if

    ! --- ------------------------------------------------------------------
    ! --- read surface pressure [Pa]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      write(lp,*) 'reading atm. climatological surface pressure...'
      open (newunit=nfu,file = atm_path(1:atm_path_len) &
           //'clim/pres_1968-1996.uf', &
           form = 'unformatted')
    end if

    do k = 1,12
      if (mnproc == 1) then
        read (nfu) atm_field_r4
        !$omp parallel do private(i)
        do j = 1,atm_jdm
          do i = 1,atm_idm
            atm_field(i,j)= &
                 atm_field_r4(i,j)*exp(9.81*atm_topo(i,j) &
                 /(287.*(atm_sktclm(i,j,k)-.00325*atm_topo(i,j))))
          end do
        end do
        !$omp end parallel do
        call inta2o(atm_field,tmp2d)
      end if
      call xcaput(tmp2d,slpres(1-nbdy,1-nbdy,k),1)
      call fill_global(atm_mval,atm_fval,halo_ps, &
           slpres(1-nbdy,1-nbdy,k))
    end do

    if (mnproc == 1) then
      close (unit = nfu)
    end if

    call xctilr(slpres, 1,12, nbdy,nbdy, halo_ps)

    ! --- ------------------------------------------------------------------
    ! --- read momentum flux [N/m^2]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      write(lp,*) 'reading atm. climatological momentum flux...'
      open (newunit=nfu,file = atm_path(1:atm_path_len) &
           //'clim/momfl_1968-1996.uf', &
           form = 'unformatted')
    end if

    do k = 1,12
      if (mnproc == 1) then
        read (nfu) atm_field_r4
        !$omp parallel do private(i)
        do j = 1,atm_jdm
          do i = 1,atm_idm
            atm_field(i,j) = atm_field_r4(i,j)
          end do
        end do
        !$omp end parallel do
        call inta2o(atm_field,tmp2d)
      end if
      call xcaput(tmp2d,taud(1-nbdy,1-nbdy,k),1)
      call fill_global(atm_mval,atm_fval,halo_ps, &
           taud(1-nbdy,1-nbdy,k))
    end do

    if (mnproc == 1) then
      close (unit = nfu)
    end if

    call xctilr(taud, 1,12, nbdy,nbdy, halo_ps)

    ! --- ------------------------------------------------------------------
    ! --- read momentum flux components [n/m^2]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      write(lp,*) &
           'reading atm. climatological u-component of momentum flux...'
      open (newunit=nfu,file = atm_path(1:atm_path_len) &
           //'clim/uflx_1968-1996.uf', &
           form = 'unformatted')
    end if

    do k = 1,12
      if (mnproc == 1) then
        read (nfu) atm_field_r4
        !$omp parallel do private(i)
        do j = 1,atm_jdm
          do i = 1,atm_idm
            atm_field(i,j) = atm_field_r4(i,j)
          end do
        end do
        !$omp end parallel do
        call inta2o(atm_field,tmp2d)
      end if
      call xcaput(tmp2d,tauxd(1-nbdy,1-nbdy,k),1)
      call fill_global(atm_mval,atm_fval,halo_pv, &
           tauxd(1-nbdy,1-nbdy,k))
    end do

    if (mnproc == 1) then
      close (unit = nfu)
    end if

    if (mnproc == 1) then
      write(lp,*) &
           'reading atm. climatological v-component of momentum flux...'
      open (newunit=nfu,file = atm_path(1:atm_path_len) &
           //'clim/vflx_1968-1996.uf', &
           form = 'unformatted')
    end if

    do k = 1,12
      if (mnproc == 1) then
        read (nfu) atm_field_r4
        !$omp parallel do private(i)
        do j = 1,atm_jdm
          do i = 1,atm_idm
            atm_field(i,j) = atm_field_r4(i,j)
          end do
        end do
        !$omp end parallel do
        call inta2o(atm_field,tmp2d)
      end if
      call xcaput(tmp2d,tauyd(1-nbdy,1-nbdy,k),1)
      call fill_global(atm_mval,atm_fval,halo_pv, &
           tauyd(1-nbdy,1-nbdy,k))
    end do

    if (mnproc == 1) then
      close (unit = nfu)
    end if

    do k = 1,12

      ! --- - change sign of momentum flux
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            tauxd(i,j,k) = -tauxd(i,j,k)
            tauyd(i,j,k) = -tauyd(i,j,k)
          end do
        end do
      end do
      !$omp end parallel do

      ! --- - rotate the vector components
      call xctilr(tauxd(1-nbdy,1-nbdy,k), 1,1, 1,1, halo_pv)
      call xctilr(tauyd(1-nbdy,1-nbdy,k), 1,1, 1,1, halo_pv)
      call uvrotr2g(tauxd(1-nbdy,1-nbdy,k),tauyd(1-nbdy,1-nbdy,k))

    end do

    call xctilr(tauxd, 1,12, nbdy,nbdy, halo_uv)
    call xctilr(tauyd, 1,12, nbdy,nbdy, halo_vv)

    ! --- ------------------------------------------------------------------
    ! --- read runoff rate [kg/m^2/s]
    ! --- ------------------------------------------------------------------

    if (mnproc == 1) then
      write(lp,*) 'reading atm. climatological runoff rate...'

      ! --- - read file containing ocean discarge points/weights for land
      ! --- - areas

      allocate(rnf_wgt  (atm_abdm,atm_idm,atm_jdm), &
               rnf_ocdpi(atm_abdm,atm_idm,atm_jdm), &
               rnf_ocdpj(atm_abdm,atm_idm,atm_jdm))

      if (expcnf == 'single_column') then
        rnf_wgt = 0.
        rnf_ocdpi = 0
        rnf_ocdpj = 0
      else
        open (newunit=nfu,file = 'runoffweights.uf', &
             form='unformatted',status='old',action = 'read')
        do j = 1,atm_jdm
          do i = 1,atm_idm
            read (nfu) rn2
            if (rn2 > atm_abdm) then
              write (lp,*) '''atm_abdm'' too small!'
              call xchalt('(inifrc_ben02clim)')
              stop '(inifrc_ben02clim)'
            end if
            do k = 1,rn2
              read (nfu) rw4,ri2,rj2
              rnf_wgt(k,i,j) = rw4
              rnf_ocdpi(k,i,j) = ri2
              rnf_ocdpj(k,i,j) = rj2
            end do
            do k = rn2+1,atm_abdm
              rnf_wgt(k,i,j) = 0.
              rnf_ocdpi(k,i,j) = 0
              rnf_ocdpj(k,i,j) = 0
            end do
          end do
        end do
        close (unit = nfu)
      end if

      open (newunit=nfu,file = atm_path(1:atm_path_len) &
           //'clim/runof_1968-1996.uf', &
           form = 'unformatted')
    end if

    do k = 1,12
      if (mnproc == 1) then
        read (nfu) atm_field_r4
        !$omp parallel do private(i)
        do j = 1,atm_jdm
          do i = 1,atm_idm
            atm_field(i,j) = atm_field_r4(i,j)
          end do
        end do
        !$omp end parallel do

        ! --- --- place runoff at ocean discharge points
        !$omp parallel do private(i)
        do j = 1,jtdm
          do i = 1,itdm
            tmp2d(i,j) = 0.
          end do
        end do
        !$omp end parallel do
        do j = 1,atm_jdm
          do i = 1,atm_idm
            do l = 1,atm_abdm
              if (rnf_wgt(l,i,j) > 0.) then
                il = rnf_ocdpi(l,i,j)
                jl = rnf_ocdpj(l,i,j)
                tmp2d(il,jl) = tmp2d(il,jl) &
                     +atm_field(i,j)*rnf_wgt(l,i,j)
              end if
            end do
          end do
        end do
      end if

      call xcaput(tmp2d,runoff(1-nbdy,1-nbdy,k),1)

      ! --- - smooth the runoff field
      call smtfld(runoff(1-nbdy,1-nbdy,k),ip, &
           atm_rnf_nsmt,atm_rnf_swgt)

      ! --- - convert unit of runoff from kg/m^2/day to kg/m^2/s
      !$omp parallel do private(l,i)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            runoff(i,j,k) = runoff(i,j,k)/86400.
          end do
        end do
      end do
      !$omp end parallel do

    end do

    if (mnproc == 1) then
      close (unit = nfu)
      deallocate(rnf_wgt, &
                 rnf_ocdpi, &
                 rnf_ocdpj)
    end if

    call xctilr(runoff, 1,12, nbdy,nbdy, halo_ps)

    ! --- Deallocate memory used for interpolation of surface fields.
    if (mnproc == 1) then
      deallocate(atm_sktclm, &
                 atm_field, &
                 atm_field_r4)
    end if
    call fnlzai

    ! --- If SSS restoring is requested, read climatological sea surface
    ! --- salinity.
    if (srxday > 0.) call rdcsss

    ! --- Compute global surface budgets of climatological heat and
    ! --- freshwater fluxes

    fwf_fac = 1.e-3 ! conversion factor kg/m^2/s -> m/s

    ! --- sum up freshwater fluxes
    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j) = 0.
          util2(i,j) = 0.
          util3(i,j) = 0.
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(k,l,i,garea,le)
    do j = 1,jj
      do k = 1,12
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            garea = scp2(i,j)*A_cgs2mks ! [m^2]

            ! --- ----- freshwater fluxes [m/s]
            util1(i,j) = util1(i,j)+precip(i,j,k)*fwf_fac*garea
            le = (2.501-.00237*(sstclm(i,j,k)-273.15))*1.e6
            util2(i,j) = util2(i,j)-(lhtflx(i,j,k)/le)*1.e-3*garea
            util3(i,j) = util3(i,j)+runoff(i,j,k)*fwf_fac*garea
          end do
        end do
      end do
    end do
    !$omp end parallel do
    call xcsum(prc_sum,util1,ip)
    call xcsum(eva_sum,util2,ip)
    call xcsum(rnf_sum,util3,ip)

    fac = 1.e-6/12.
    prc_sum = prc_sum*fac
    eva_sum = eva_sum*fac
    rnf_sum = rnf_sum*fac

    ! --- sum up heat fluxes
    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j) = 0.
          util2(i,j) = 0.
          util3(i,j) = 0.
          util4(i,j) = 0.
        end do
      end do
    end do
    !$omp end parallel do

    do k = 1,12

      ! --- - compute albedo
      dangle = 8.*atan(1.)*real(30*k-15)/360.
      call albw_eval(dangle,plat,clouds(1-nbdy,1-nbdy,k),albdif,albw)
      !$omp parallel do private(l,i,garea,albedo)
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            garea = scp2(i,j)*A_cgs2mks ! [m^2]

            ! --- ----- heat fluxes
            albedo = albs_f*ricclm(i,j,k)+albw(i,j)*(1.-ricclm(i,j,k))
            util1(i,j) = util1(i,j)+dswrfl(i,j,k)*(1.-albedo)*garea
            util2(i,j) = util2(i,j)-nlwrfs(i,j,k)*garea
            util3(i,j) = util3(i,j)-lhtflx(i,j,k)*garea
            util4(i,j) = util4(i,j)-shtflx(i,j,k)*garea
          end do
        end do
      end do
      !$omp end parallel do
    end do

    call xcsum(swa_sum,util1,ip)
    call xcsum(lwa_sum,util2,ip)
    call xcsum(lht_sum,util3,ip)
    call xcsum(sht_sum,util4,ip)

    fac = (L_mks2cgs*L_mks2cgs)/(12.*area)
    swa_ave = swa_sum*fac
    lwa_ave = lwa_sum*fac
    lht_ave = lht_sum*fac
    sht_ave = sht_sum*fac

    if (mnproc == 1) then
      write (lp,*)
      write (lp,*) 'Global precipitation:         ',prc_sum,'Sv'
      write (lp,*) 'Global evaporation:           ',eva_sum,'Sv'
      write (lp,*) 'Global runoff:                ',rnf_sum,'Sv'
      write (lp,*) 'Global balance of freshwater: ', &
           prc_sum+eva_sum+rnf_sum,'Sv'
      write (lp,*)
      write (lp,*) 'Global mean short-wave radiation: ',swa_ave,'W/m^2'
      write (lp,*) 'Global mean long-wave radiation:  ',lwa_ave,'W/m^2'
      write (lp,*) 'Global mean latent heat flux:     ',lht_ave,'W/m^2'
      write (lp,*) 'Global mean sensible heat flux:   ',sht_ave,'W/m^2'
      write (lp,*) 'Global balance of mean heat-flux: ', &
           swa_ave+lwa_ave+lht_ave+sht_ave,'W/m^2'
      write (lp,*)
    end if

    ! --- Adjust runoff and short-wave radiation to balance freshwater and
    ! --- heat fluxes, respectively.
    if (expcnf == 'single_column') then
      crnf = 0.0
    else
      crnf = -(prc_sum+eva_sum)/rnf_sum
    end if
    cswa = -(lwa_ave+lht_ave+sht_ave)/swa_ave
    !$omp parallel do private(k,l,i)
    do j = 1-nbdy,jj+nbdy
      do k = 1,12
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            runoff(i,j,k) = runoff(i,j,k)*crnf
            dswrfl(i,j,k) = dswrfl(i,j,k)*cswa
          end do
        end do
      end do
    end do
    !$omp end parallel do
    if (mnproc == 1) then
      write (lp,*) 'Runoff has been adjusted by a factor',crnf
      write (lp,*) 'Short-wave radiation has been adjusted by a factor', &
           cswa
      write (lp,*)
    end if

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'inifrc_ben02clim:'
      end if
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
    end if

  end subroutine inifrc_ben02clim

  ! --- ------------------------------------------------------------------

  subroutine inifrc_ben02syn

    ! --- ------------------------------------------------------------------
    ! --- Prepare utilization of daily atmospheric forcing
    ! --- ------------------------------------------------------------------

    real :: dx2,dy2
    real*4 :: rw4
    integer :: errstat,i,j,k,l,nfu
    integer*2 :: rn2,ri2,rj2

    ! --- Allocate memory for daily forcing fields.
    allocate(taud  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5), &
             tauxd (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5), &
             tauyd (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5), &
             dswrfl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5), &
             nlwrfs(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5), &
             shtflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5), &
             lhtflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5), &
             precip(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5), &
             clouds(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5), &
             slpres(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5), &
             runoff(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5), &
             tmpsfc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5), &
             ricec (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,5))

    ! --- Allocate memory for transfer coefficients, gustiness squared, and
    ! --- air density
    allocate(cd_d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             ch_d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             ce_d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             wg2_d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             cd_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             ch_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             ce_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             wg2_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             rhoa(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))

    ! --- Allocate memory for accumulation variables
    allocate(tsi_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             tml_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             sml_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             alb_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             fice_tda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             tsi(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))

    ! --- Allocate memory for derivative of non-solar heat flux by surface
    ! --- temperature, albedos and instantaneous runoff flux and runoff
    ! --- reservoar
    allocate(dfl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             albw(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             alb(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             rnfins(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             rnfres(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))

    ! --- Number of iteration in the computation of transfer coefficients.
    tciter = 1

    ! --- Prepare interpolation of surface fields.
    call initai

    ! --- Read climatological sea-ice concentration.
    call rdcsic

    ! --- Read climatological surface temperature.
    call rdctsf

    ! --- If SSS restoring is requested, read climatological sea surface
    ! --- salinity.
    if (srxday > 0.) call rdcsss

    ! --- Initialize diagnosing/application of relaxation fluxes
    call idarlx

    ! --- Compute smoothing weights atm_ice_swgt and atm_rnf_swgt. For
    ! --- stability atm_ice_swgt,atm_rnf_swgt < .5*dx^2*dy^2/(dx^2+dy^2).
    atm_ice_swgt = spval
    atm_rnf_swgt = spval
    !$omp parallel do private(l,i,dx2,dy2) &
    !$omp reduction(min:atm_ice_swgt,atm_rnf_swgt)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          dx2 = scpx(i,j)*scpx(i,j)
          dy2 = scpy(i,j)*scpy(i,j)
          atm_ice_swgt = min(atm_ice_swgt,.5*dx2*dy2/(dx2+dy2))
          atm_rnf_swgt = min(atm_rnf_swgt,.5*dx2*dy2/(dx2+dy2))
        end do
      end do
    end do
    !$omp end parallel do
    call xcmin(atm_ice_swgt)
    atm_ice_swgt = .9*atm_ice_swgt
    call xcmin(atm_rnf_swgt)
    atm_rnf_swgt = .9*atm_rnf_swgt

    ! --- Number of smoothing iterations is choosen to get a geographical
    ! --- extent of the smoothing independent of the grid resolution
    if (mnproc == 1) then
      if     (atm_idm == atm_idm_ncep.and. &
           atm_jdm == atm_jdm_ncep) then
        atm_ice_nsmt = nint(atm_ice_csmt_ncep/atm_ice_swgt)
        atm_rnf_nsmt = nint(atm_rnf_csmt_ncep/atm_rnf_swgt)
      else if (atm_idm == atm_idm_era .and. &
           atm_jdm == atm_jdm_era ) then
        atm_ice_nsmt = nint(atm_ice_csmt_era/atm_ice_swgt)
        atm_rnf_nsmt = nint(atm_rnf_csmt_era/atm_rnf_swgt)
      end if
    end if
    call xcbcst(atm_ice_nsmt)
    call xcbcst(atm_rnf_nsmt)
    if (mnproc == 1) then
      write (lp,*) 'atm_ice_swgt',atm_ice_swgt, &
           'atm_ice_nsmt',atm_ice_nsmt
      write (lp,*) 'atm_rnf_swgt',atm_rnf_swgt, &
           'atm_rnf_nsmt',atm_rnf_nsmt
    end if

    ! --- set runoff and short-wave radiation adjustment factors
    if (mnproc == 1) then
      if     (atm_idm == atm_idm_ncep.and. &
           atm_jdm == atm_jdm_ncep) then
        atm_crnf = atm_crnf_ncep
        atm_cswa = atm_cswa_ncep
      else if (atm_idm == atm_idm_era .and. &
           atm_jdm == atm_jdm_era ) then
        atm_crnf = atm_crnf_era
        atm_cswa = atm_cswa_era
      end if
    end if
    call xcbcst(atm_crnf)
    call xcbcst(atm_cswa)

    if (mnproc == 1) then

      ! --- - read file containing ocean discarge points/weights for land
      ! --- - areas

      allocate(rnf_wgt  (atm_abdm,atm_idm,atm_jdm), &
               rnf_ocdpi(atm_abdm,atm_idm,atm_jdm), &
               rnf_ocdpj(atm_abdm,atm_idm,atm_jdm))

      open (newunit=nfu,file = 'runoffweights.uf', &
           form='unformatted',status='old',action = 'read')
      do j = 1,atm_jdm
        do i = 1,atm_idm
          read (nfu) rn2
          if (rn2 > atm_abdm) then
            write (lp,*) '''atm_abdm'' too small!'
            call xchalt('(inifrc_ben02syn)')
            stop '(inifrc_ben02syn)'
          end if
          do k = 1,rn2
            read (nfu) rw4,ri2,rj2
            rnf_wgt(k,i,j) = rw4
            rnf_ocdpi(k,i,j) = ri2
            rnf_ocdpj(k,i,j) = rj2
          end do
          do k = rn2+1,atm_abdm
            rnf_wgt(k,i,j) = 0.
            rnf_ocdpi(k,i,j) = 0
            rnf_ocdpj(k,i,j) = 0
          end do
        end do
      end do
      close (unit = nfu)

    end if

    ! --- Initial indexes of the 5 time levels
    l1gi = 1
    l2gi = 2
    l3gi = 3
    l4gi = 4
    l5gi = 5

    ! --- Go back 4 days so that all necessary fields are read
    errstat = date_offset(calendar,date,-4)
    if (errstat /= calendar_noerr) then
      if (mnproc == 1) then
        write (lp, '(2a)') ' inifrc_ben02syn: date_offset error: ', &
             trim(calendar_errstr(errstat))
      end if
      call xcstop('(inifrc_ben02syn)')
      stop '(inifrc_ben02syn)'
    end if
    call set_day_of_year
    do i = 1,4
      call rdatm_syn
      errstat = date_offset(calendar,date,1)
      if (errstat /= calendar_noerr) then
        if (mnproc == 1) then
          write (lp, '(2a)') ' inifrc_ben02syn: date_offset error: ', &
               trim(calendar_errstr(errstat))
        end if
        call xcstop('(inifrc_ben02syn)')
        stop '(inifrc_ben02syn)'
      end if
      call set_day_of_year
    end do

  end subroutine inifrc_ben02syn

  ! --- ------------------------------------------------------------------

  subroutine getfrc_ben02clim

    ! --- ------------------------------------------------------------------
    ! --- Get climatological forcing
    ! --- ------------------------------------------------------------------

    ! --- Set interpolation parameters
    xgi = xmi
    l1gi = l1mi
    l2gi = l2mi
    l3gi = l3mi
    l4gi = l4mi
    l5gi = l5mi

    ! --- ------------------------------------------------------------------
    ! --- Compute the air-sea fluxes
    ! --- ------------------------------------------------------------------

    call asflux

  end subroutine getfrc_ben02clim

  ! --- ------------------------------------------------------------------

  subroutine getfrc_ben02syn

    ! --- ------------------------------------------------------------------
    ! --- Get synoptic forcing
    ! --- ------------------------------------------------------------------

    ! --- Set interpolation parameter
    xgi = real(mod(nstep-1,nstep_in_day)+1)/nstep_in_day

    ! --- ------------------------------------------------------------------
    ! --- The first time step of a new day, read new forcing fields
    ! --- ------------------------------------------------------------------

    if (mod(nstep,nstep_in_day) == 1) call rdatm_syn

    ! --- ------------------------------------------------------------------
    ! --- Compute the air-sea fluxes
    ! --- ------------------------------------------------------------------

    call asflux

  end subroutine getfrc_ben02syn

end module mod_ben02
