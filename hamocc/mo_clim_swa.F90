! Copyright (C) 2021  J. Tjiputra 
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


MODULE mo_clim_swa
!******************************************************************************
!
! MODULE mo_clim_swa - Variables and routines for climatology short-wave fields
!
!  J.Tjiputra,        *NORCE Climate, Bergen*    2021-04-15
!
!  Modified
!  --------
!
!  Purpose
!  -------
!   Declaration, memory allocation, and routines related to swa_clim fields
!
!******************************************************************************
  implicit none

  private
  public :: ini_swa_clim, swaclimfile, swa_clim

  ! File name (incl. full path) for input data, set through namelist 
  ! in hamocc_init.F
  character(len=512), save :: swaclimfile=''
  ! Array to store swa flux after reading from file
  real, allocatable,  save :: swa_clim(:,:,:)


CONTAINS
!******************************************************************************



SUBROUTINE ini_swa_clim(kpie,kpje,omask)
!******************************************************************************
!
! INI_SWA_CLIM - initialise the climatology SWA field MODULE.
!
!  J.Tjiputra             *NORCE Climate, Bergen*       2021-04-15
!
!  Purpose
!  -------
!   Initialise the climatology swa MODULE, read in the swa (short-wave radiation) data set.
!
!  Parameter list:
!  ---------------
!   *integer*   *kpie*    - 1st dimension of model grid.
!   *integer*   *kpje*    - 2nd dimension of model grid.
!   *real*      *omask*   - land/ocean mask (1=ocean)
!
!******************************************************************************
  use netcdf,         only: nf90_noerr,nf90_nowrite,nf90_close,nf90_open 
  use mod_xc,         only: mnproc,xchalt
  use mo_control_bgc, only: io_stdo_bgc

  implicit none

  integer,          intent(in) :: kpie,kpje
  real,             intent(in) :: omask(kpie,kpje)

  integer :: i,j
  integer :: ncid,ncstat,ncvarid,errstat


  ! allocate field to hold swa fields
  if (mnproc.eq.1) then
    write(io_stdo_bgc,*)' '
    write(io_stdo_bgc,*)'***************************************************'
    write(io_stdo_bgc,*)'iHAMOCC: Initialization of module mo_clim_swa:'
    write(io_stdo_bgc,*)' '
  endif

  if (mnproc.eq.1) then
    write(io_stdo_bgc,*)'Memory allocation for variable swa_clim ...'
    write(io_stdo_bgc,*)'First dimension    : ',kpie
    write(io_stdo_bgc,*)'Second dimension   : ',kpje
  endif
   
  ALLOCATE (swa_clim(kpie,kpje,1),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory swa_clim'
  swa_clim(:,:,1) = 0.0

  ! Open netCDF data file     
  if (mnproc==1) then
    ncstat = NF90_OPEN(trim(swaclimfile),NF90_NOWRITE, ncid)
    if (ncstat.NE.NF90_NOERR ) then
      call xchalt('(ini_swa_clim: Problem with netCDF1)')
             stop '(ini_swa_clim: Problem with netCDF1)'
    endif
  endif
  
  ! Read  data
  call read_netcdf_var(ncid,'swa',swa_clim(1,1,1),1,1,0)

  ! Close file
  if (mnproc==1) then
    ncstat = NF90_CLOSE(ncid)
    if ( ncstat .NE. NF90_NOERR ) then
      call xchalt('(ini_swa_clim: Problem with netCDF200)')
             stop '(ini_swa_clim: Problem with netCDF200)'
    endif
  endif

  if (mnproc.eq.1) then 
    write(io_stdo_bgc,*) ''
    write(io_stdo_bgc,*) 'ini_swa_clim: Using climatology swa file '//trim(swaclimfile) 
  endif

  ! set flux to zero over land
  do j=1,kpje
    do i=1,kpie

      if(omask(i,j).lt.0.5) swa_clim(i,j,1) = 0.0
    
    enddo
  enddo


  return

!******************************************************************************
END SUBROUTINE ini_swa_clim


!******************************************************************************
END MODULE mo_clim_swa

