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


module mo_fedep
!******************************************************************************
!
! MODULE mo_fedep - Variables and routines for iron deposition
!
!
!  J.Schwinger,        *NORCE Climate, Bergen*    2020-05-27
!
!  Modified
!  --------
!
!  Purpose
!  -------
!   Declaration, memory allocation, and routines related to iron deposition
!
!
!******************************************************************************
  implicit none

  private
  public :: ini_fedep, get_fedep,fedepfile

  ! File name (incl. full path) for input data, set through namelist 
  ! in hamocc_init.F
  character(len=512), save :: fedepfile=''
  ! Array to store dust deposition flux after reading from file
  real, allocatable,  save :: dustflx(:,:,:)


contains
!******************************************************************************



subroutine ini_fedep(kpie,kpje,omask)
!******************************************************************************
!
! INI_FEDEP - initialise the iron deposition module.
!
!
!  J.Schwinger            *NORCE Climate, Bergen*       2020-05-19
!
!  Purpose
!  -------
!   Initialise the iron deposition module, read in the iron (dust) data set.
!
!  Parameter list:
!  ---------------
!   *INTEGER*   *kpie*    - 1st dimension of model grid.
!   *INTEGER*   *kpje*    - 2nd dimension of model grid.
!   *REAL*      *omask*   - land/ocean mask (1=ocean)
!
!******************************************************************************
  use netcdf
  use mod_xc,         only: mnproc,xchalt
  use mo_control_bgc, only: io_stdo_bgc

  implicit none

  integer,          intent(in) :: kpie,kpje
  real,             intent(in) :: omask(kpie,kpje)

  integer :: i,j,l
  integer :: ncid,ncstat,ncvarid,errstat


  ! allocate field to hold iron deposition fluxes
  IF (mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)' '
    WRITE(io_stdo_bgc,*)'***************************************************'
    WRITE(io_stdo_bgc,*)'iHAMOCC: Initialization of module mo_fedep:'
    WRITE(io_stdo_bgc,*)' '
  ENDIF

  IF (mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)'Memory allocation for variable dustflx ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    WRITE(io_stdo_bgc,*)'Third dimension    :  12'
  ENDIF
   
  ALLOCATE (dustflx(kpie,kpje,12),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory dustflx'
  dustflx(:,:,:) = 0.0

  ! Open netCDF data file     
  IF(mnproc==1) THEN
    ncstat = NF90_OPEN(trim(fedepfile),NF90_NOWRITE, ncid)
    IF (ncstat.NE.NF90_NOERR ) THEN
      CALL xchalt('(get_dust: Problem with netCDF1)')
             stop '(get_dust: Problem with netCDF1)'
    END IF
  END IF

  ! Read  data
  call read_netcdf_var(ncid,'DUST',dustflx(1,1,1),12,0,0)

  ! Close file
  IF(mnproc==1) THEN
    ncstat = NF90_CLOSE(ncid)
    IF ( ncstat .NE. NF90_NOERR ) THEN
      CALL xchalt('(get_dust: Problem with netCDF200)')
             stop '(get_dust: Problem with netCDF200)'
    END IF
  END IF

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


  RETURN

!******************************************************************************
end subroutine ini_fedep


subroutine get_fedep(kpie,kpje,kplmon,dust)
!******************************************************************************
!
! GET_FEDEP - get iron (dust) deposition for current month
!
!
!  J.Schwinger            *NORCE Climate, Bergen*       2020-05-19
!
!  Purpose
!  -------
!   Initialise the iron deposition module, read in the iron (dust) data set.
!
!  Parameter list:
!  ---------------
!   *INTEGER*   *kpie*    - 1st dimension of model grid.
!   *INTEGER*   *kpje*    - 2nd dimension of model grid.
!   *INTEGER*   *kplmon*  - current month.
!   *REAL*      *dust*    - dust flux for current month

!
!******************************************************************************
  integer,          intent(in)  :: kpie,kpje,kplmon
  real,             intent(out) :: dust(kpie,kpje)

  dust = dustflx(:,:,kplmon)


!******************************************************************************
end subroutine get_fedep


!******************************************************************************
end module mo_fedep
