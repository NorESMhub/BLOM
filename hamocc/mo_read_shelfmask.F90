! Copyright (C) 2021-2024  j. maerz
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


module mo_read_shelfmask

  !*************************************************************************************************
  ! Routine to read and/or initialize an ocean shelf-sea mask either from netcdf file or from
  ! internal bathymetry values, if use_shelfsea_res_time = .true.
  !
  ! If a file is read, it should hold real values with 1. for shelf-sea regions and 0 elsewhere
  ! and the variable name is 'shelfmask'
  !
  ! The shelf-sea mask is currently used to calculate the shelf-sea water residence time
  !
  ! j.maerz *UiB, Bergen* 2024-10-31
  !
  !*************************************************************************************************

  implicit none
  private

  public :: ini_read_shelfmask

  character(len=512),               public :: shelfsea_maskfile =''
  logical, allocatable,  protected, public :: shelfmask(:,:)

contains

  subroutine ini_read_shelfmask(kpie,kpje,kbnd,pbath,omask)

    use mod_xc,         only: mnproc,xchalt
    use mod_dia,        only: iotype
    use mo_control_bgc, only: use_shelfsea_res_time,io_stdo_bgc
    use mo_param_bgc,   only: shelfbreak_depth
    use netcdf,         only: nf90_open,nf90_close,nf90_nowrite
    use mo_netcdf_bgcrw,only: read_netcdf_var

    ! Arguments
    integer, intent(in) :: kpie                                     ! 1st dimension of model grid
    integer, intent(in) :: kpje                                     ! 2nd dimension of model grid
    integer, intent(in) :: kbnd                                     ! number of halo grid points
    real,    intent(in) :: pbath(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd) ! bathymetry fields - water depth [m]
    real,    intent(in) :: omask(kpie,kpje)                         ! land/ocean mask.

    ! Local variables
    logical :: file_exists=.false.
    integer :: i,j,errstat,ncid,ncstat
    real,allocatable  :: mask(:,:)

    ! Always allocate field to hold shelfsea mask
    if(mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable shelfmask ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate(shelfmask(kpie,kpje),stat=errstat)
    allocate(mask(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory shelfmask'
    shelfmask(:,:) = .false.
    mask = 0.

    ! Check, if we are going to run with shelf-sea water residence time tracers
    if (.not.use_shelfsea_res_time) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*) ''
        write(io_stdo_bgc,*) '******************************************************'
        write(io_stdo_bgc,*) 'read_shelfmask: shelf sea age tracer is not activated.'
      endif
      return
    endif

    ! Check if shelfsea mask file exists. If not, run in default mode.
    inquire(file=shelfsea_maskfile,exist=file_exists)
    if (.not. file_exists .and. mnproc.eq.1) then
      write(io_stdo_bgc,*) ''
      write(io_stdo_bgc,*) '************************************************************'
      write(io_stdo_bgc,*) 'read_shelfmask: Cannot find (shelf sea) region mask file.   '
      write(io_stdo_bgc,*) 'Fallback to default ...                                     '
      write(io_stdo_bgc,*) '... using internal bathymetry data to reconstruct the mask  '
    endif

    if (file_exists) then
      ! read shelf sea mask from file
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*) ''
        write(io_stdo_bgc,'(a)') 'read_shelfmask: read mask from ',trim(shelfsea_maskfile)
      endif
      ncstat=nf90_open(trim(shelfsea_maskfile),nf90_nowrite,ncid)
      call read_netcdf_var(ncid,'shelfmask',mask,1,1,0)
      ncstat=nf90_close(ncid)
    else
      ! reconstruct shelf sea mask from internal bathymetry
      !$OMP DO PARALLEL PRIVATE (i,j)
      do j=1,kpje
        do i=1,kpie
          if((omask(i,j) > 0.5) .and. (pbath(i,j) <= shelfbreak_depth)) then
            mask(i,j) = 1.
          endif
        enddo
      enddo
      !$OMP END PARALLEL DO
    endif

    !$OMP DO PARALLEL PRIVATE (i,j)
    ! Eventually fill the logical shelfsea mask field
    do j = 1,kpje
      do i = 1,kpie
        if (1 == nint(mask(i,j))) then
          shelfmask(i,j) = .true.
        else
          shelfmask(i,j) = .false.
        endif
      enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine ini_read_shelfmask

end module mo_read_shelfmask
