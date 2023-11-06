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


module mo_clim_swa
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


contains
  !******************************************************************************



  subroutine ini_swa_clim(kpie,kpje,omask)
    !******************************************************************************
    !
    ! INI_SWA_CLIM - initialise the climatology SWA field module.
    !
    !  J.Tjiputra             *NORCE Climate, Bergen*       2021-04-15
    !
    !  Purpose
    !  -------
    !   Initialise the climatology swa module, read in the swa (short-wave radiation) data set.
    !
    !  Parameter list:
    !  ---------------
    !   *INTEGER*   *kpie*    - 1st dimension of model grid.
    !   *INTEGER*   *kpje*    - 2nd dimension of model grid.
    !   *REAL*      *omask*   - land/ocean mask (1=ocean)
    !
    !******************************************************************************
    use netcdf,         only: nf90_noerr,nf90_nowrite,nf90_close,nf90_open
    use mod_xc,         only: mnproc,xchalt
    use mo_control_bgc, only: io_stdo_bgc
    use mo_read_netcdf_var, only: read_netcdf_var

    implicit none

    ! Arguments
    integer, intent(in) :: kpie,kpje
    real,    intent(in) :: omask(kpie,kpje)

    ! Local variables
    integer :: i,j
    integer :: ncid,ncstat,ncvarid,errstat


    ! allocate field to hold swa fields
    IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'***************************************************'
      WRITE(io_stdo_bgc,*)'iHAMOCC: Initialization of module mo_clim_swa:'
      WRITE(io_stdo_bgc,*)' '
    ENDIF

    IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable swa_clim ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    ENDIF

    ALLOCATE (swa_clim(kpie,kpje,1),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory swa_clim'
    swa_clim(:,:,1) = 0.0

    ! Open netCDF data file
    IF(mnproc==1) THEN
      ncstat = NF90_OPEN(trim(swaclimfile),NF90_NOWRITE, ncid)
      IF (ncstat.NE.NF90_NOERR ) THEN
        CALL xchalt('(ini_swa_clim: Problem with netCDF1)')
        stop '(ini_swa_clim: Problem with netCDF1)'
      END IF
    END IF

    ! Read  data
    call read_netcdf_var(ncid,'swa',swa_clim(1,1,1),1,1,0)

    ! Close file
    IF(mnproc==1) THEN
      ncstat = NF90_CLOSE(ncid)
      IF ( ncstat .NE. NF90_NOERR ) THEN
        CALL xchalt('(ini_swa_clim: Problem with netCDF200)')
        stop '(ini_swa_clim: Problem with netCDF200)'
      END IF
    END IF

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


    RETURN

    !******************************************************************************
  end subroutine ini_swa_clim


  !******************************************************************************
end module mo_clim_swa
