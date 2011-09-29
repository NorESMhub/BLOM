module mo_Gdata_read

!********************************************************************************
!     J.Schwinger,        *Gfi, Bergen*    19.05.2011
!
! Modified
! --------
!     
! Purpose
! -------
!  - Routines for reading WOA and GLODAP gridded data netCDF files for
!    initialisation of HAMOCC fields
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
!  -nz
!     Number of z-levels in the data files.
!
!  -zlev
!     Depth of each z-level [m] in the current data file.
!
!
!********************************************************************************

use netcdf
use mod_xc,         only: mnproc,lp,xchalt
use mo_control_bgc, only: io_stdo_bgc

implicit none

private

public:: set_Gdata,clean_Gdata,get_profile,get_region,nz,zlev

! Number of latitudes, longitudes, and z-levels in the WOA and GLODAP data
integer, parameter          :: nlon = 360
integer, parameter          :: nlat = 180
integer, parameter          :: nz   = 33

! Resolution of data in degree
real, parameter             :: dres = 1.0


! Max number of gridpoints to select around the center for averaging in 
! longitude direction
integer, parameter          :: dnmax = 100.0

! Subdirectories for WOA and GLODAP data within global input-data directory:
character(len=*), parameter :: WOA_dir = 'WOA_2009'
character(len=*), parameter :: GLP_dir = 'GLODAP'

! Fill value used in this module, original fill values of data files are
! replaced by this fill value during read
real,             parameter :: fillval = -999.9

! Variables set by call to Gdata_set
character(len=16)                     :: var,ncname
character(len=3)                      :: dsrc
character(len=256)                    :: file
real                                  :: cfac, ddeg
real, dimension(:),       allocatable :: lon,lat,zlev
real, dimension(:, :, :), allocatable :: rvar,data

logical, save :: lset  = .false.


!-----------------------------------------
! Definitions for regional mean profiles:
!-----------------------------------------
type region
   character(len=64)  :: name       ! Region name
   integer            :: idx        ! Region index
   integer            :: npts(nz)   ! nb of valid data points at each level
   real               :: clon, clat ! center longitude and latitude
   real               :: dlon, dlat ! latitude and longitude extent
   real               :: mprf(nz)   ! mean profile for region
   logical            :: global     ! global extent T/F
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


subroutine set_Gdata(path,vname,inddeg)
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
!  path:         path to the data directory; it is assumed that WOA and GLODAP
!                data reside in subdirectories WOA_dir and GLP_dir (as
!                defined in the header of this module) 
!  vname:        data set name to read in; valid names are
!                'pho' - WOA phosphate
!                'nit' - WOA nitrate
!                'sil' - WOA silicate
!                'oxy' - WOA dissolved oxygen
!                'alk' - GLODAP alkalinity
!                'dic' - GLODAP dissolved inorganic carbon
!  inddeg:       extent (in degrees) of region used for averaging
!
!--------------------------------------------------------------------------------
character(len=*), intent(in) :: path,vname
real,             intent(in) :: inddeg

! Local variables
character(len=*), parameter  :: routinestr = 'set_Gdata'


if( allocated(lon)  ) deallocate( lon  ) 
if( allocated(lat)  ) deallocate( lat  ) 
if( allocated(zlev) ) deallocate( zlev ) 
if( allocated(rvar) ) deallocate( rvar )
if( allocated(data) ) deallocate( data )

! Select settings specific to each variable
select case (vname)

case ('pho') ! phosphate
   file   = trim(path)//'/'//trim(WOA_dir)//'/'//'phosphate_annual_1deg.nc'
   ncname = 'p_an'
   dsrc   = 'WOA'
   cfac   = 1.0  ! no unit conversion within this module

case ('nit') ! nitrate
   file   = trim(path)//'/'//trim(WOA_dir)//'/'//'nitrate_annual_1deg.nc'
   ncname = 'n_an'
   dsrc   = 'WOA'
   cfac   = 1.0  ! no unit conversion within this module

case ('sil') ! silicate
   file   = trim(path)//'/'//trim(WOA_dir)//'/'//'silicate_annual_1deg.nc'
   ncname = 'i_an'
   dsrc   = 'WOA'
   cfac   = 1.0  ! no unit conversion within this module

case ('oxy') ! oxygen
   file   = trim(path)//'/'//trim(WOA_dir)//'/'//'dissolved_oxygen_annual_1deg.nc'
   ncname = 'o_an'
   dsrc   = 'WOA'
   cfac   = 44.661  ! conversion ml/l -> micromoles/l

case ('alk') ! alkalinity
   file   = trim(path)//'/'//trim(GLP_dir)//'/'//'Alk.nc'
   ncname = 'Alk'
   dsrc   = 'GLP'
   cfac   = 1.0  ! no unit conversion within this module

case ('dic') ! DIC
   file   = trim(path)//'/'//trim(GLP_dir)//'/'//'TCO2.nc'
   ncname = 'TCO2'
   dsrc   = 'GLP'
   cfac   = 1.0  ! no unit conversion within this module

case default
   call moderr(routinestr,'Invalid vname')  

end select

var  = vname
ddeg = inddeg

if(mnproc == 1) write(io_stdo_bgc,*) 'HAMOCC: initialising ', trim(vname)

call read_Gdata()

! extend data array by +/-dnmax data points in longitude
allocate( data(-dnmax:nlon+dnmax,nlat,nz) )
data(:,:,:)                 = 0.0
data( 1:nlon,          :,:) = rvar(:,:,:)
data(-dnmax:0,         :,:) = rvar(nlon-dnmax:nlon,:,:)
data(nlon+1:nlon+dnmax,:,:) = rvar(1:dnmax,:,:)

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
real,              intent(out) :: prf(nz)

! Local variables
integer                        :: idx, npts(nz)
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
! than 10 for the surface layer. A global mean profile is used if
! get_region returns 0.
if( npts(1) < 10 ) then

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
integer                     :: ncId, varId 
integer                     :: LonDimId,   lonId,   numlon
integer                     :: LatDimId,   latId,   numlat
integer                     :: DepthDimId, depthId, numlev
integer                     :: ndim, natts
integer                     :: i,l
integer                     :: dimid(7)
integer                     :: status
real, allocatable           :: tmp(:,:,:)
real                        :: fval
character(len=16)           :: lonstr,latstr,fvalstr
character(len=*), parameter :: routinestr = 'mo_Gdata_read, read_Gdata'
logical                     :: lonfirst


select case (dsrc)

case ('WOA')
  lonstr   = 'lon'
  latstr   = 'lat'
  fvalstr  = '_FillValue'
  lonfirst = .true.

case ('GLP')
  lonstr   = 'longitude'
  latstr   = 'latitude'
  fvalstr  = 'missing_value'
  lonfirst = .false.

case default
   call moderr(routinestr,'Invalid dsrc')  

end select


! Open file
if(mnproc == 1) write(io_stdo_bgc,*) 'Reading ', trim(file)
status = nf90_open(file,nf90_nowrite,ncid); call ncerr(status)


! Get dimensions
status = nf90_inq_dimid(ncid, trim(lonstr), LonDimId)
call ncerr(status)
status = nf90_inquire_dimension(ncid, LonDimID, len=numlon)
call ncerr(status)

status = nf90_inq_dimid(ncid, trim(latstr), LatDimId)
call ncerr(status)
status = nf90_inquire_dimension(ncid, LatDimID, len=numlat)
call ncerr(status)

status = nf90_inq_dimid(ncid, "depth", DepthDimId)
call ncerr(status)
status = nf90_inquire_dimension(ncid, DepthDimID, len=numlev)
call ncerr(status)


if( numlon /= nlon .or. numlat /= nlat .or. numlev /= nz ) &
     call moderr(routinestr,'Unexpected nb of elements in data file')  

allocate( lon(nlon), lat(nlat), zlev(nz) )
allocate( rvar(nlon,nlat,nz) )

! Get lon, lat, and lev
status = nf90_inq_varid(ncid, trim(lonstr), lonId)
call ncerr(status)
status = nf90_get_var(ncid, lonId, lon)
call ncerr(status)

status = nf90_inq_varid(ncid, trim(latstr), latId)
call ncerr(status)
status = nf90_get_var(ncid, latId, lat)
call ncerr(status)

status = nf90_inq_varid(ncid, "depth", depthId)
call ncerr(status)
status = nf90_get_var(ncid, depthId, zlev)
call ncerr(status)

! Get varid and fill value
status = nf90_inq_varid(ncid, ncname, varId)
call ncerr(status)
status = nf90_inquire_variable(ncid, varId, ndims=ndim, dimids=dimid, nAtts=natts)
call ncerr(status)
status =  nf90_get_att(ncid, varid, trim(fvalstr), fval)
call ncerr(status)

! Get data values, switch order of lat and long for GLODAP data
select case (dsrc)

case ('WOA') 
   ! just read the data
   status = nf90_get_var(ncid, varId, rvar)
   call ncerr(status)

   ! Replace fill values:
   where( rvar > fval*0.1 ) rvar = fillval

   ! unit conversion
   where( rvar > 0.0 ) rvar = rvar*cfac

case ('GLP')
   ! read data into temporary field
   allocate( tmp(nlat,nlon,nz) )
   status = nf90_get_var(ncid, varId, tmp)
   call ncerr(status)
   ! switch order of lat/lon
   do l=1,nz
      rvar(:,:,l) = transpose( tmp(:,:,l) )
   end do
   deallocate( tmp )
   ! arrange data to correspond to [0,360] in longitude
   lon  = cshift(lon, -180)
   rvar = cshift(rvar,-180,1)
   do i=1,nlon
      if(lon(i)<0.0) lon(i)=lon(i)+360.0
   end do

   ! Replace fill values:
   where( rvar < fval*0.1 ) rvar = fillval

   ! unit conversion
   where( rvar > 0.0 ) rvar = rvar*cfac

end select

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
real,              intent(out) :: prf(nz)
integer,           intent(out) :: npts(nz)
logical, optional, intent(in)  :: global

! Local variables
integer                        :: ilonc, ilons, ilone, dnlon
integer                        :: ilatc, ilats, ilate, dnlat
integer                        :: k,l, nelmlon,nelmlat
logical                        :: gl = .false.
character(len=*), parameter    :: routinestr = 'mo_Gdata_read, calc_mean_profile'


if( .not. lset   ) call moderr(routinestr, ' Module not initialised yet')
if( clon < 0     ) call moderr(routinestr, ' clon must be in the range [0,360]')
if( clon > 360.0 ) call moderr(routinestr, ' clon must be in the range [0,360]')


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
   ! goes beyond the pole. Insted rectangle is adjusted such that boundaries
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

   npts(l) = count(data(ilons:ilone,ilats:ilate,l)>0.0)
   prf(l)  = sum(data(ilons:ilone,ilats:ilate,l), mask=data(ilons:ilone,ilats:ilate,l)>0.0)
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

if( allocated(lon)  ) deallocate( lon  ) 
if( allocated(lat)  ) deallocate( lat  ) 
if( allocated(zlev) ) deallocate( zlev ) 
if( allocated(rvar) ) deallocate( rvar ) 
if( allocated(data) ) deallocate( data ) 

file   = ''
ncname = ''
var    = ''
dsrc   = ''
cfac   = 1.0
ddeg   = 0.0
lset   = .false.

!--------------------------------------------------------------------------------
end subroutine clean_Gdata



subroutine ncerr(status)
!--------------------------------------------------------------------------------
! Handle netCDF-errors
!--------------------------------------------------------------------------------
integer, intent(in) :: status

if(status == nf90_NoErr) return

write(lp,*) 'NetCDF error: ',nf90_strerror(status)
call flush(lp)
call xchalt('(Module mo_Gdata_read, ncerr)')
stop '(Module mo_Gdata_read, ncerr)'


!--------------------------------------------------------------------------------
end subroutine ncerr


subroutine moderr(routinestr,errstr)
!--------------------------------------------------------------------------------
! Handle errors, which occur in this module
!--------------------------------------------------------------------------------
character(len=*), intent(in) :: routinestr,errstr


write(lp,'(/3a)') routinestr, ': ', errstr
call flush(lp)
call xchalt('(Module mo_Gdata_read)')
stop '(Module mo_Gdata_read)'


!--------------------------------------------------------------------------------
end subroutine moderr


!********************************************************************************
end module mo_Gdata_read
