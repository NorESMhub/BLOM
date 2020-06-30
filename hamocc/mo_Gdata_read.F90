! Copyright (C) 2020  J. Schwinger
!
! This file is part of BLOM/iHAMOCC.
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
! along with BLOM. If not, see https://www.gnu.org/licenses/.


module mo_Gdata_read

!********************************************************************************
!     J.Schwinger,        *Gfi, Bergen*           2011-05-19
!
! Modified
! --------
!    J.Schwinger,         *Uni Climate, BCCR*     2017-07-07
!    - adapted this module to read the initial conditions for OMIP-BGC.
!
!     J.Schwinger,      *Uni Research, Bergen*   2018-04-12
!     - adaptions for reading c-isotope initial values as d13C and d14C
!     
! Purpose
! -------
!  - Routines for reading initial condition files for OMIP-BGC, which are based
!    on WOA 2013 and GLODAPv2 gridded data netCDF files
!
! Description:
! ------------
!  Public routines and variable of this module:
!
!  -subroutine set_Gdata
!     Initialise global varibles and read in one data set. Must be
!     called before the processing of one data set starts.
!
!  -subroutine clean_Gdata
!     Deallocate global fields of this module and reset all global variables.
!     Should be called each time, the processing of one data set is finished.
!
!  -subroutine get_profile
!     Returns one profile from the currently open data set (opened by a
!     previous call to set_Gdata). See header of get profile for details.
!
!  -function get_region
!     Returns the index of the region a given point belongs to. If no region
!     is found get_region returns 0, which is the index of the 'global region'.
!     Note that the regions are defined below in the module header.
!
!  -nz_woa
!     Number of z-levels in the WOA data files.
!
!  -nz_glo
!     Number of z-levels in the GLODAP data files.
!
!  -nzmax
!     Max nuber of z-levels (=nzwoa)
!
!  -zlev
!     Depth of each z-level [m] in the current data file.
!
!
!********************************************************************************

use netcdf
use mod_xc,         only: mnproc,xchalt
use mo_control_bgc, only: io_stdo_bgc,io_nml

implicit none

private

public :: set_Gdata,clean_Gdata,get_profile,get_region,nzmax,nz,zlev,zlev_bnds,fillval
public :: inidic,inialk,inipo4,inioxy,inino3,inisil,inid13c,inid14c

! Number of latitudes, longitudes, and z-levels in the WOA and GLODAP data
integer, parameter          :: nlon   = 360
integer, parameter          :: nlat   = 180
integer, parameter          :: nz_woa = 102
integer, parameter          :: nz_glo = 33
integer, parameter          :: nzmax  = nz_woa
! Resolution of data in degree
real, parameter             :: dres = 1.0


! Max number of gridpoints to select around the center for averaging in 
! longitude direction
integer, parameter          :: dnmax = 100.0


! Fill value used in this module, original fill values of data files are
! replaced by this fill value during read
real,             parameter :: fillval = -1.e+32

! Input file names (incl. full path) set through namelist
character(len=256), save                    :: inidic  = ''
character(len=256), save                    :: inialk  = ''
character(len=256), save                    :: inipo4  = ''
character(len=256), save                    :: inioxy  = ''
character(len=256), save                    :: inino3  = ''
character(len=256), save                    :: inisil  = ''
character(len=256), save                    :: inid13c = ''
character(len=256), save                    :: inid14c = ''
character(len=256), save                    :: inic13 = ''   ! currently not used
character(len=256), save                    :: inic14 = ''   ! currently not used

! Variables set by call to Gdata_set
integer, save                               :: nz
real, save                                  :: cfac, ddeg
real, save, dimension(:),       allocatable :: lon,lat,zlev
real, save, dimension(:,:),     allocatable :: zlev_bnds
real, save, dimension(:, :, :), allocatable :: rvar,gdata
character(len=16),  save                    :: var,ncname
character(len=3) ,  save                    :: dsrc
character(len=256), save                    :: infile

logical, save :: lset  = .false.


!-----------------------------------------
! Definitions for regional mean profiles:
!-----------------------------------------
type region
   character(len=64)  :: name        ! Region name
   integer            :: idx         ! Region index
   integer            :: npts(nzmax) ! nb of valid data points at each level
   real               :: clon, clat  ! center longitude and latitude
   real               :: dlon, dlat  ! latitude and longitude extent
   real               :: mprf(nzmax) ! mean profile for region
   logical            :: global      ! global extent T/F
end type region

integer, parameter    :: nreg=10
type(region)          :: rg(0:nreg)


! Set regions for fall-back profiles

! Global profile;
data  rg(0)%idx,  rg(0)%name / 0, 'global'  /
data  rg(0)%clon, rg(0)%clat /   0.0,   0.0 /
data  rg(0)%dlon, rg(0)%dlat / 360.0, 180.0 /
data  rg(0)%global           / .true.       /

! Indian Ocean
data  rg(1)%idx,  rg(1)%name / 1, 'Indian Ocean'  /
data  rg(1)%clon, rg(1)%clat /  65.0,-10.0 /
data  rg(1)%dlon, rg(1)%dlat /  90.0, 80.0 /
data  rg(1)%global           / .false.      /

! North Atlantic
data  rg(2)%idx,  rg(2)%name / 2, 'North Atlantic'  /
data  rg(2)%clon, rg(2)%clat /   0.0,  70.0 /
data  rg(2)%dlon, rg(2)%dlat / 180.0,  40.0 /
data  rg(2)%global           / .false.      /

! northern subtropical Atlantic
data  rg(3)%idx,  rg(3)%name / 3, 'Northern subtropical Atlantic'  /
data  rg(3)%clon, rg(3)%clat / 330.0,  35.0 /
data  rg(3)%dlon, rg(3)%dlat / 140.0,  30.0 /
data  rg(3)%global           / .false.      /

! Tropical Atlantic
data  rg(4)%idx,  rg(4)%name / 4, 'Tropical Atlantic'  /
data  rg(4)%clon, rg(4)%clat / 335.0,   0.0 /
data  rg(4)%dlon, rg(4)%dlat /  90.0,  40.0 /
data  rg(4)%global           / .false.      /

! Southern subtropical Atlantic
data  rg(5)%idx,  rg(5)%name / 5, 'Southern subtropical Atlantic'  /
data  rg(5)%clon, rg(5)%clat / 335.0, -35.0 /
data  rg(5)%dlon, rg(5)%dlat /  90.0,  30.0 /
data  rg(5)%global           / .false.      /

! North Pacific
data  rg(6)%idx,  rg(6)%name / 6, 'North Pacific'  /
data  rg(6)%clon, rg(6)%clat / 180.0,  70.0 /
data  rg(6)%dlon, rg(6)%dlat / 180.0,  40.0 /
data  rg(6)%global           / .false.      /

! northern subtropical Pacific
data  rg(7)%idx,  rg(7)%name / 7, 'Northern subtropical Pacific'  /
data  rg(7)%clon, rg(7)%clat / 185.0,  35.0 /
data  rg(7)%dlon, rg(7)%dlat / 150.0,  30.0 /
data  rg(7)%global           / .false.      /

! Tropical Pacific
data  rg(8)%idx,  rg(8)%name / 8, 'Tropical Pacific'  /
data  rg(8)%clon, rg(8)%clat / 200.0,   0.0 /
data  rg(8)%dlon, rg(8)%dlat / 180.0,  40.0 /
data  rg(8)%global           / .false.      /

! Southern subtropical Pacific
data  rg(9)%idx,  rg(9)%name / 9, 'Southern subtropical Pacific'  /
data  rg(9)%clon, rg(9)%clat / 200.0, -35.0 /
data  rg(9)%dlon, rg(9)%dlat / 180.0,  30.0 /
data  rg(9)%global           / .false.      /


! Southern Ocean
data  rg(10)%idx,  rg(10)%name / 10, 'Southern Ocean'  /
data  rg(10)%clon, rg(10)%clat / 180.0, -70.0 /
data  rg(10)%dlon, rg(10)%dlat / 360.0,  40.0 /
data  rg(10)%global            / .false.      /


!********************************************************************************
contains


subroutine set_Gdata(vname,inddeg)
!--------------------------------------------------------------------------------
!
! Purpose:
! --------
!  Initialise global varibles and read data set specified by vname. Must be
!  called before the first call to any routine of this module.
!
! Description:
! ------------
! 
!
! Arguments:
! ----------
!  vname:        data set name to read in; valid names are
!                'pho' - WOA phosphate
!                'nit' - WOA nitrate
!                'sil' - WOA silicate
!                'oxy' - WOA dissolved oxygen
!                'alk' - GLODAP alkalinity
!                'dic' - GLODAP dissolved inorganic carbon
!                'C13' - Dissolved inorganic 13C carbon isotope
!                'd13' - delta13C of dissolved inorganic carbon
!                'C14' - Dissolved inorganic 14C carbon isotope
!                'd14' - delta14C of dissolved inorganic carbon
!  inddeg:       extent (in degrees) of region used for averaging
!
!--------------------------------------------------------------------------------
character(len=*), intent(in) :: vname
real,             intent(in) :: inddeg

! Local variables
character(len=*), parameter  :: routinestr = 'set_Gdata'


if( allocated(lon)       ) deallocate( lon  ) 
if( allocated(lat)       ) deallocate( lat  ) 
if( allocated(zlev)      ) deallocate( zlev ) 
if( allocated(zlev_bnds) ) deallocate( zlev_bnds ) 
if( allocated(rvar)      ) deallocate( rvar )
if( allocated(gdata)     ) deallocate( gdata )

! Select settings specific to each variable
select case (vname)

case ('pho') ! phosphate
   infile = inipo4
   ncname = 'po4'
   dsrc   = 'WOA'
   cfac   = 1.0e-6  ! data in mumol/L -> kmol/m3

case ('nit') ! nitrate
   infile = inino3
   ncname = 'no3'
   dsrc   = 'WOA'
   cfac   = 1.0e-6  ! data in mumol/L -> kmol/m3

case ('sil') ! silicate
   infile = inisil
   ncname = 'si'
   dsrc   = 'WOA'
   cfac   = 1.0e-6  ! data in mumol/L -> kmol/m3

case ('oxy') ! oxygen
   infile = inioxy
   ncname = 'o2'
   dsrc   = 'WOA'
   cfac   = 44.661*1.0e-6  ! conversion ml/L -> mumol/L -> kmol/m3

case ('alk') ! alkalinity
   infile = inialk
   ncname = 'At'
   dsrc   = 'GLO'
   cfac   = 1.0e-6  ! data in mumol/kg -> mol/kg

case ('dic') ! DIC
   infile   = inidic
   ncname = 'Ct_preind'
   dsrc   = 'GLO'
   cfac   = 1.0e-6  ! data in mumol/kg -> mol/kg

case ('C13') ! natural 13C [micromoles/kg]
   infile = inic13
   ncname = 'C13'
   dsrc   = 'ISO'
   cfac   = 1.0e-6  ! data in mumol/kg -> mol/kg

case ('d13') ! natural delta13C [permil]
   infile = inid13c
   ncname = 'd13C'
   dsrc   = 'ISO'
   cfac   = 1.0

case ('C14') ! natural 14C [micromoles/kg]
   infile = inic14
   ncname = 'C14'
   dsrc   = 'ISO'
   cfac   = 1.0e-6  ! data in mumol/kg -> mol/kg

case ('d14') ! natural delta14C [permil]
   infile = inid14c
   ncname = 'd14C'
   dsrc   = 'ISO'
   cfac   = 1.0

case default
   call moderr(routinestr,'Invalid vname')  

end select

var  = vname
ddeg = inddeg

if(mnproc == 1) write(io_stdo_bgc,*) 'iHAMOCC: initialising ', trim(vname)

call read_Gdata()

! extend data array by +/-dnmax data points in longitude
allocate( gdata(-dnmax:nlon+dnmax,nlat,nz) )
gdata(:,:,:)                 = 0.0
gdata( 1:nlon,          :,:) = rvar(:,:,:)
gdata(-dnmax:0,         :,:) = rvar(nlon-dnmax:nlon,:,:)
gdata(nlon+1:nlon+dnmax,:,:) = rvar(1:dnmax,:,:)

lset  = .true.

call set_regional_profiles()


!--------------------------------------------------------------------------------
end subroutine set_Gdata



subroutine get_profile(clon,clat,prf)
!--------------------------------------------------------------------------------
!
! Purpose:
! --------
!  Return a profile suitable for initialisation of HAMCC at point clon/clat.
!
! Description:
! ------------
!  A mean profile is calculated by calling calc_mean_profile with the settings
!  defined by a previous call to set_Gdata. If no valid data is found for the
!  point clon/clat, it is tried to obtain a mean regional profile (e.g. for the
!  north atlantic area). These mean profiles are initialised as part of 
!  set_Gdata.
! 
!
! Arguments:
! ----------
!  clon, clat:   center lon/lat of mean profile
!  prf:          mean profile for initialisation
!
!--------------------------------------------------------------------------------
real,              intent(in)  :: clon, clat  
real,              intent(out) :: prf(nzmax)

! Local variables
integer                        :: idx, npts(nzmax)
real                           :: clon_tmp,clat_tmp
character(len=*), parameter    :: routinestr = 'mo_Gdata_read, get_profile'


if( .not. lset ) call moderr(routinestr, ' Module not initialised yet')


if( clon < 0 ) then
   clon_tmp=clon+360.0
   clat_tmp=clat
else
   clon_tmp=clon
   clat_tmp=clat
end if

! Try to obtain a mean profile for a region centered at clon/clat
call calc_mean_profile(clon_tmp,clat_tmp,ddeg,ddeg,prf,npts)


! Fall back to regional profile if number of valid data points is smaller
! than 3 for the surface layer. A global mean profile is used if
! get_region returns 0.
if( npts(1) < 3 ) then

   idx = get_region(clon_tmp,clat_tmp)
   prf = rg(idx)%mprf
   !write(*,*) 'Region is ', rg(idx)%name, clon, clat

end if


!--------------------------------------------------------------------------------
end subroutine get_profile



function get_region(clon,clat)
!--------------------------------------------------------------------------------
!
! Purpose:
! --------
!  Return index of region the point clon/clat belongs to
!
! Description:
! ------------
!  The rectangular regions as defined in the module header (and stored in the
!  data type 'rg') are searched. If point clon/clat belongs to region i, the 
!  index i is the result of this function. If no region is found, get_region
!  returns 0, which is the index of the 'global' region defined in the header.
!
! Arguments:
! ----------
!  clon, clat:   lon/lat of point
!
!--------------------------------------------------------------------------------
real, intent(in)               :: clon,clat
integer                        :: get_region

! Local variables
integer                        :: i
real                           :: ll_lon, ur_lon
real                           :: ll_lat, ur_lat
logical                        :: boundwithin, found
character(len=*), parameter    :: routinestr = 'mo_Gdata_read, get_region'


if( clon < 0     ) call moderr(routinestr, ' clon must be in the range [0,360]')
if( clon > 360.0 ) call moderr(routinestr, ' clon must be in the range [0,360]')

found       = .false.

do i=1,nreg

   boundwithin = .false.

   ll_lon = rg(i)%clon-rg(i)%dlon/2.0
   ur_lon = rg(i)%clon+rg(i)%dlon/2.0
   ll_lat = rg(i)%clat-rg(i)%dlat/2.0
   ur_lat = rg(i)%clat+rg(i)%dlat/2.0

   if( ll_lon < 0.0    ) ll_lon = ll_lon+360.0
   if( ur_lon > 360.0  ) ur_lon = ur_lon-360.0

   if( ll_lon > ur_lon ) boundwithin = .true.

   if( clat < ll_lat .or. clat > ur_lat ) cycle

   if( boundwithin ) then
 
      if( clon < ll_lon .and. clon > ur_lon ) cycle

   else

      if( clon < ll_lon .or. clon > ur_lon ) cycle

   end if

   found = .true.
   exit

end do

if( found ) then
   get_region = rg(i)%idx
else
   get_region = 0
end if


!--------------------------------------------------------------------------------
end function get_region



subroutine set_regional_profiles()
!--------------------------------------------------------------------------------
!
! Calculate the mean profiles in regions as defined in the module header
!
!--------------------------------------------------------------------------------

! Local variables
integer                     :: i
character(len=*), parameter :: routinestr = 'mo_Gdata_read, set_regional_profiles'

if( .not. lset ) call moderr(routinestr, ' Module not initialised yet')

do i=0,nreg

   call calc_mean_profile(rg(i)%clon,rg(i)%clat,rg(i)%dlon,rg(i)%dlat, &
                          rg(i)%mprf,rg(i)%npts,rg(i)%global)


   !write(*,*) 'Calculated mean profile for ', rg(i)%name
   !write(*,*) '==============='
   !write(*,*) rg(i)%mprf
   !write(*,*) '==============='

end do

!--------------------------------------------------------------------------------
end subroutine set_regional_profiles



subroutine read_Gdata()
!--------------------------------------------------------------------------------
!
! Read the WOA or GLODAP data into variables lon/lat/zlev and rvar
!
!--------------------------------------------------------------------------------

! Local variables
integer                     :: ncId, vId, dId 
integer                     :: numlon, numlat, numlev
integer                     :: i, ndim, natts
integer                     :: dimid(7)
integer                     :: status
real                        :: fval
character(len=16)           :: lonstr,latstr,depthstr,depthbndsstr,fvalstr
character(len=*), parameter :: routinestr = 'mo_Gdata_read, read_Gdata'


lonstr   = 'lon'
latstr   = 'lat'
fvalstr  = '_FillValue'

select case (dsrc)

case ('WOA')
   nz = nz_woa
   depthstr='depth'
   depthbndsstr='depth_bnds'
case ('GLO')
   nz = nz_glo
   depthstr='depthz'
   depthbndsstr='depthz_bnds'
case ('ISO')
   nz = nz_glo
   depthstr='depthz'
   depthbndsstr='depthz_bnds'
case default
   call moderr(routinestr,'Invalid dsrc')  

end select


! Open file
if(mnproc == 1) write(io_stdo_bgc,*) 'Reading ', trim(infile)
status = nf90_open(infile,nf90_nowrite,ncid); call ncerr(status)


! Get dimensions
status = nf90_inq_dimid(ncid, trim(lonstr), dId)
call ncerr(status)
status = nf90_inquire_dimension(ncid, dID, len=numlon)
call ncerr(status)

status = nf90_inq_dimid(ncid, trim(latstr), dId)
call ncerr(status)
status = nf90_inquire_dimension(ncid, dID, len=numlat)
call ncerr(status)

status = nf90_inq_dimid(ncid, trim(depthstr), dId)
call ncerr(status)
status = nf90_inquire_dimension(ncid, dId, len=numlev)
call ncerr(status)

if( numlon /= nlon .or. numlat /= nlat .or. numlev /= nz ) &
     call moderr(routinestr,'Unexpected nb of elements in data file')  

allocate( lon(nlon), lat(nlat), zlev(nz), zlev_bnds(2,nz) )
allocate( rvar(nlon,nlat,nz) )

! Get lon, lat, and lev
status = nf90_inq_varid(ncid, trim(lonstr), vId)
call ncerr(status)
status = nf90_get_var(ncid, vId, lon)
call ncerr(status)

status = nf90_inq_varid(ncid, trim(latstr), vId)
call ncerr(status)
status = nf90_get_var(ncid, vId, lat)
call ncerr(status)

status = nf90_inq_varid(ncid, trim(depthstr), vId)
call ncerr(status)
status = nf90_get_var(ncid, vId, zlev)
call ncerr(status)

status = nf90_inq_varid(ncid, trim(depthbndsstr), vId)
call ncerr(status)
status = nf90_get_var(ncid, vId, zlev_bnds)
call ncerr(status)

! Get varid and fill value
status = nf90_inq_varid(ncid, ncname, vId)
call ncerr(status)
status = nf90_inquire_variable(ncid, vId, ndims=ndim, dimids=dimid, nAtts=natts)
call ncerr(status)
status =  nf90_get_att(ncid, vid, trim(fvalstr), fval)
call ncerr(status)

! GetRead the data 
status = nf90_get_var(ncid, vId, rvar)
call ncerr(status)

! arrange data to correspond to [0,360] in longitude
select case (dsrc)

case ('WOA')
   lon  = cshift(lon, -180)
   rvar = cshift(rvar,-180,1)
case ('GLO')
   lon  = cshift(lon, -20)
   rvar = cshift(rvar,-20,1)
case ('ISO')
   lon  = cshift(lon, -180)
   rvar = cshift(rvar,-180,1)
end select

do i=1,nlon
   if(lon(i)<  0.0) lon(i)=lon(i)+360.0
   if(lon(i)>360.0) lon(i)=lon(i)-360.0
end do

! Fillvalues are assumed to be < 0 currently, otherwise the below code would fail
if(fval > 0.0) call moderr(routinestr,'FillValue > 0 found in data')  

where( rvar < fval*0.1 )
   ! Replace fill values:
   rvar = fillval
elsewhere
   ! unit conversion
   rvar = rvar*cfac
end where


! Close data file
status = nf90_close(ncid)
call ncerr(status)


!--------------------------------------------------------------------------------
end subroutine read_Gdata



subroutine calc_mean_profile(clon,clat,dlon,dlat,prf,npts,global)
!--------------------------------------------------------------------------------
!
! Purpose:
! --------
!  Return mean profile around the center point clon/clat.
!
! Description:
! ------------
!  The mean profile is calculated from valid data points in the square defined
!  by clon+/-dlon/2 clat+/-dlat/2. The number of valid data points per depth 
!  level is returned in npts. By setting the optional argument global to true,
!  all valid data points are used to calculate a global mean profile. clon, clat,
!  dlon, and dlat are ignored in this case.
!
! Arguments:
! ----------
!  clon, clat:   center lon/lat of mean profile
!  dlon, dlat:   lon/lat extent of region to select for averaging
!  prf:          mean profile calculated from all data in selected region
!  npts:         nb of valid data points found for each depth level
!  global:       if set to true, calculate mean over the whole data set
!
!--------------------------------------------------------------------------------
real,              intent(in)  :: clon, clat  
real,              intent(in)  :: dlon, dlat 
real,              intent(out) :: prf(nzmax)
integer,           intent(out) :: npts(nzmax)
logical, optional, intent(in)  :: global

! Local variables
integer                        :: ilonc, ilons, ilone, dnlon
integer                        :: ilatc, ilats, ilate, dnlat
integer                        :: l, nelmlon,nelmlat
logical                        :: gl = .false.
character(len=*), parameter    :: routinestr = 'mo_Gdata_read, calc_mean_profile'


if( .not. lset   ) call moderr(routinestr, ' Module not initialised yet')
if( clon < 0     ) call moderr(routinestr, ' clon must be in the range [0,360]')
if( clon > 360.0 ) call moderr(routinestr, ' clon must be in the range [0,360]')

prf(:)  = fillval
npts(:) = 0.0

if( present(global) ) gl=global

if( gl ) then

   ilons=1
   ilone=nlon
   ilats=1
   ilate=nlat

else

   ! Find index of nearest gridpoint (not exact but okay for this purpose)
   do ilonc=1,nlon
      if( clon < lon(ilonc) ) exit
   end do
   if( ilonc > nlon ) ilonc = nlon
   if( lon(ilonc)-clon > dres/2.0 ) ilonc=ilonc-1
   if( ilonc < 1 ) ilonc = 1
   
   do ilatc=1,nlat
      if( clat < lat(ilatc) ) exit
   end do
   if( ilatc > nlat ) ilatc = nlat
   if( lat(ilatc)-clat > dres/2.0 ) ilatc=ilatc-1
   if( ilatc < 1 ) ilatc = 1
   
   dnlon  = int(dlon/2.0*dres) ! Nb of gridpoints to select around the center lon 
   dnlat  = int(dlat/2.0*dres) ! Nb of gridpoints to select around the center lat
   
   nelmlon  = 2*dnlon+1
   nelmlat  = 2*dnlat+1
   
   ! Start idices of rectangle:
   ilons = ilonc-dnlon
   ilats = ilatc-dnlat
   
   ! There is no "wrap-around" if southern/northen boundary of rectangle
   ! goes beyond the pole. Instead rectangle is adjusted such that boundaries
   ! are aligned with northernmost/sothernmost data gridpoint
   if(ilats <= 0            ) ilats=1
   if(ilats > nlat-nelmlat+1) ilats= nlat-nelmlat+1
   
   ! End indices of rectangle:
   ilone = ilons+nelmlon-1
   ilate = ilats+nelmlat-1

   if( ilons < -dnmax      ) call moderr(routinestr,'error: data array too small')  
   if( ilone >  dnmax+nlon ) call moderr(routinestr,'error: data array too small')  

end if


! Calculate mean profile:
do l=1,nz

   npts(l) = count(gdata(ilons:ilone,ilats:ilate,l) > fillval*0.1)
   prf(l)  = sum(gdata(ilons:ilone,ilats:ilate,l), mask=gdata(ilons:ilone,ilats:ilate,l) > fillval*0.1)
   if( npts(l) > 0) then
      prf(l) = prf(l)/npts(l)
   else
      prf(l) = fillval
   end if

end do


!write(*,*) '================'
!if( gl ) then
!   write(*,*) 'global'
!else
!   write(*,*) dnlon,dnlat
!   write(*,*) ilonc,ilons,ilone,lon(ilonc)
!   write(*,*) ilatc,ilats,ilate,lat(ilatc)
!end if
!write(*,*) '================'

!--------------------------------------------------------------------------------
end subroutine calc_mean_profile



subroutine clean_Gdata()
!--------------------------------------------------------------------------------
! Deallocate fields and reset global variables
!--------------------------------------------------------------------------------

if( allocated(lon)        ) deallocate( lon   ) 
if( allocated(lat)        ) deallocate( lat   ) 
if( allocated(zlev)       ) deallocate( zlev  ) 
if( allocated(zlev_bnds)  ) deallocate( zlev_bnds  ) 
if( allocated(rvar)       ) deallocate( rvar  ) 
if( allocated(gdata)      ) deallocate( gdata ) 

infile = ''
ncname = ''
var    = ''
dsrc   = ''
cfac   = 1.0
ddeg   = 0.0
nz     = 0
lset   = .false.

!--------------------------------------------------------------------------------
end subroutine clean_Gdata



subroutine ncerr(status)
!--------------------------------------------------------------------------------
! Handle netCDF-errors
!--------------------------------------------------------------------------------
integer, intent(in) :: status

if(status == nf90_NoErr) return

write(io_stdo_bgc,*) 'NetCDF error: ',nf90_strerror(status)
write(io_stdo_bgc,*) 'Abort... '
call flush(io_stdo_bgc)
call xchalt('(Module mo_Gdata_read, ncerr)')
stop '(Module mo_Gdata_read, ncerr)'


!--------------------------------------------------------------------------------
end subroutine ncerr


subroutine moderr(routinestr,errstr)
!--------------------------------------------------------------------------------
! Handle errors, which occur in this module
!--------------------------------------------------------------------------------
character(len=*), intent(in) :: routinestr,errstr


write(io_stdo_bgc,'(/3a)') routinestr, ': ', errstr
write(io_stdo_bgc,*) 'Abort... '
call flush(io_stdo_bgc)
call xchalt('(Module mo_Gdata_read)')
stop '(Module mo_Gdata_read)'


!--------------------------------------------------------------------------------
end subroutine moderr


!********************************************************************************
end module mo_Gdata_read
