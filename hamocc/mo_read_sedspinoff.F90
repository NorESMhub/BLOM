! Copyright (C) 2020  S. Gao, I. Bethke, J. Tjiputra, J. Schwinger
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

module mo_read_sedspinoff
  !*************************************************************************************************
  ! Routine for reading climatological sedimentation fluxes for the offline sediment spinup
  !
  !*************************************************************************************************

  implicit none
  private

  ! Routintes
  public :: read_sedspinoff ! read sediment porosity file

  ! Module variables
  character(len=512), public :: sedspinoff_file = ''

  real,dimension(:,:), allocatable, protected, public :: clim_prcaca
  real,dimension(:,:), allocatable, protected, public :: clim_prorca
  real,dimension(:,:), allocatable, protected, public :: clim_produs
  real,dimension(:,:), allocatable, protected, public :: clim_silpro

contains

  subroutine read_sedspinoff(kpie,kpje,omask)

    use mod_xc,             only: mnproc,xchalt
    use mo_control_bgc,     only: io_stdo_bgc,offline_sediment_spinup
    use netcdf,             only: nf90_noerr,nf90_nowrite,nf90_close,nf90_open
    use mo_netcdf_bgcrw,    only: read_netcdf_var

    implicit none

    ! Arguments
    integer, intent(in)    :: kpie
    integer, intent(in)    :: kpje
    real,    intent(in)    :: omask(kpie,kpje)

    ! Local variables
    integer :: i,j
    logical :: file_exists = .false.
    integer :: ncid,ncstat,errstat
    real    :: tmp_prorca(kpie,kpje)
    real    :: tmp_prcaca(kpie,kpje)
    real    :: tmp_produs(kpie,kpje)
    real    :: tmp_silpro(kpie,kpje)

    ! Return if offline_sediment_spinup is turned off
    if (.not. offline_sediment_spinup) then
      if (mnproc == 1) then
        write(io_stdo_bgc,*) ''
        write(io_stdo_bgc,*) 'read_sedspinoff: offline sediment spinup is not activated'
      endif
      return
    endif

    ! Check if offline-sediment-spinup file  file exists. If not, abort.
    inquire(file=sedspinoff_file,exist=file_exists)
    if (.not. file_exists .and. mnproc == 1) then
      write(io_stdo_bgc,*) ''
      write(io_stdo_bgc,*) 'read_sedspinoff: Cannot find offline-sediment-spinup file ... '
      call xchalt('(read_sedspinoff)')
      stop        '(read_sedspinoff)'
    endif

    ! read offline-sediment-spinup forcing info from file
    if (mnproc == 1) then
      write(io_stdo_bgc,*) ''
      write(io_stdo_bgc,*) 'read_sedspinoff: read offline-sediment-spinup information from ',      &
                                             trim(sedspinoff_file)
    endif

    ! Open netCDF data file
    if (mnproc == 1) then
      ncstat = NF90_OPEN(trim(sedspinoff_file),NF90_NOWRITE, ncid)
      if (ncstat /= NF90_NOERR ) then
        call xchalt('(read_sedspinoff: Problem with netCDF1)')
        stop        '(read_sedspinoff: Problem with netCDF1)'
      end if
    end if

    if (mnproc == 1) then
      write(io_stdo_bgc,*)'Memory allocation for variable clim_prorca ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (clim_prorca(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory clim_prorca'
    clim_prorca(:,:) = 0.0

    if (mnproc == 1) then
      write(io_stdo_bgc,*)'Memory allocation for variable clim_prcaca ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (clim_prcaca(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory clim_prcaca'
    clim_prcaca(:,:) = 0.0

    if (mnproc == 1) then
      write(io_stdo_bgc,*)'Memory allocation for variable clim_produs ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (clim_produs(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory clim_produs'
    clim_produs(:,:) = 0.0

    if (mnproc == 1) then
      write(io_stdo_bgc,*)'Memory allocation for variable clim_silpro ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate (clim_silpro(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory clim_silpro'
    clim_silpro(:,:) = 0.0

    ! Read  data - expecting data to be in format:
    !   - carflx_bot:  mol C m-2 s-1
    !   - calflx_bot:  mol Ca m-2 s-1
    !   - dustflx_bot: g m-2 s-1
    !   - bsiflx_bot:  mol Si m-2 s-1
    ! as usually provided as output by iHAMOCC
    call read_netcdf_var(ncid,'carflx_bot',tmp_prorca,1,1,0)
    call read_netcdf_var(ncid,'calflx_bot',tmp_prcaca,1,1,0)
    call read_netcdf_var(ncid,'dustflx_bot',tmp_produs,1,1,0)
    call read_netcdf_var(ncid,'bsiflx_bot',tmp_silpro,1,1,0)

    ! Close file
    if (mnproc == 1) then
      ncstat = NF90_CLOSE(ncid)
      if ( ncstat /=  NF90_NOERR ) then
        call xchalt('(read_sedspinoff: Problem with netCDF2)')
        stop        '(read_sedspinoff: Problem with netCDF2)'
      end if
    end if

    do j=1,kpje
      do i=1,kpie
        if (omask(i,j) > 0.5) then
          clim_prorca(i,j) = tmp_prorca(i,j)
          clim_prcaca(i,j) = tmp_prcaca(i,j)
          clim_produs(i,j) = tmp_produs(i,j)
          clim_silpro(i,j) = tmp_silpro(i,j)
        endif
      enddo
    enddo


  end subroutine read_sedspinoff

end module mo_read_sedspinoff
