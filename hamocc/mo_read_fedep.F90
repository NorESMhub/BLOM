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


MODULE mo_read_fedep
!******************************************************************************
!
! MODULE mo_read_fedep - routines for reading iron deposition data
!
!
!  J.Schwinger,      *NORCE Climate, Bergen*   2020-05-27
!
!  Modified
!  --------
!  J. Schwinger,     *NORCE climate, Bergen*   2022-06-02
!  -revise structure of this MODULE, split into a MODULE for reading the 
!   data (mo_read_fedep) and a MODULE that applies the fluxes in core 
!   hamocc (mo_apply_fedep)
!
!  Purpose
!  -------
!   Declaration, memory allocation, and routines related to reading iron 
!   deposition input data
!
!  Description:
!  ------------
!  Public routines and variable of this MODULE:
!
!  -SUBROUTINE ini_read_fedep
!     Initialise the MODULE for reading iron deposition data
!
!  -SUBROUTINE get_fedep
!     Get the iron (dust) deposition for a given month
!
!
!******************************************************************************
  IMPLICIT NONE

  private
  public :: ini_read_fedep,get_fedep,fedepfile

  ! File name (incl. full path) for input data, set through namelist 
  ! in hamocc_init.F
  character(len=512), save :: fedepfile=''
  ! Array to store dust deposition flux after reading from file
  real, allocatable,  save :: dustflx(:,:,:)


CONTAINS
!******************************************************************************



SUBROUTINE ini_read_fedep(kpie,kpje,omask)
!******************************************************************************
!
! INI_FEDEP - initialise the iron deposition MODULE.
!
!
!  J.Schwinger            *NORCE Climate, Bergen*       2020-05-19
!
!  Purpose
!  -------
!   Initialise the iron deposition MODULE, read in the iron (dust) data set.
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

  IMPLICIT NONE

  integer, intent(in) :: kpie,kpje
  real,    intent(in) :: omask(kpie,kpje)

  integer             :: i,j,l
  integer             :: ncid,ncstat,ncvarid,errstat


  ! allocate field to hold iron deposition fluxes
  if (mnproc.eq.1) then
    write(io_stdo_bgc,*)' '
    write(io_stdo_bgc,*)'***************************************************'
    write(io_stdo_bgc,*)'iHAMOCC: Initialization of module mo_fedep:'
    write(io_stdo_bgc,*)' '
  endif

  if (mnproc.eq.1) then
    write(io_stdo_bgc,*)'Memory allocation for variable dustflx ...'
    write(io_stdo_bgc,*)'First dimension    : ',kpie
    write(io_stdo_bgc,*)'Second dimension   : ',kpje
    write(io_stdo_bgc,*)'Third dimension    :  12'
  endif
   
  ALLOCATE (dustflx(kpie,kpje,12),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory dustflx'
  dustflx(:,:,:) = 0.0

  ! Open netCDF data file     
  if (mnproc==1) then
    ncstat = NF90_OPEN(trim(fedepfile),NF90_NOWRITE, ncid)
    if (ncstat.NE.NF90_NOERR ) then
      call xchalt('(get_dust: Problem with netCDF1)')
             stop '(get_dust: Problem with netCDF1)'
    endif
  endif

  ! Read  data
  call read_netcdf_var(ncid,'DUST',dustflx(1,1,1),12,0,0)

  ! Close file
  if (mnproc==1) then
    ncstat = NF90_CLOSE(ncid)
    if ( ncstat .NE. NF90_NOERR ) then
      call xchalt('(get_dust: Problem with netCDF200)')
             stop '(get_dust: Problem with netCDF200)'
    endif
  endif

  if (mnproc.eq.1) then 
    write(io_stdo_bgc,*) ''
    write(io_stdo_bgc,*) 'ini_fedep: Using dust deposition file '//trim(fedepfile) 
  endif

  ! set flux to zero over land
  do l=1,12
    do j=1,kpje
    do i=1,kpie

      if(omask(i,j).lt.0.5) dustflx(i,j,l) = 0.0
    
    enddo
    enddo
  enddo


  return

!******************************************************************************
END SUBROUTINE ini_read_fedep


SUBROUTINE get_fedep(kpie,kpje,kplmon,dust)
!******************************************************************************
!
! GET_FEDEP - get iron (dust) deposition for current month
!
!
!  J.Schwinger            *NORCE Climate, Bergen*       2020-05-19
!
!  Purpose
!  -------
!   Initialise the iron deposition MODULE, read in the iron (dust) data set.
!
!  Parameter list:
!  ---------------
!   *integer*   *kpie*    - 1st dimension of model grid.
!   *integer*   *kpje*    - 2nd dimension of model grid.
!   *integer*   *kplmon*  - current month.
!   *real*      *dust*    - dust flux for current month

!
!******************************************************************************
  integer, intent(in)  :: kpie,kpje,kplmon
  real,    intent(out) :: dust(kpie,kpje)

  dust = dustflx(:,:,kplmon)


!******************************************************************************
END SUBROUTINE get_fedep


!******************************************************************************
END MODULE mo_read_fedep
