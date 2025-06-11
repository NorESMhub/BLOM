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

module mo_read_sedqual
  !*************************************************************************************************
  ! Routine for reading sediment POC age from netcdf file use_sediment_quality must be set to true
  !
  ! The model attempts to read the 3D sediment age information and the climatological prorca
  ! from the input file 'SEDQUALFILE' (incl. full path)
  !
  ! sed_POCage_init holds then the age that can be applied later in mo_sedmnt (ini_sed_qual)
  !*************************************************************************************************

  use mo_kind, only: bgc_fnmlen

  implicit none
  private

  ! Routintes
  public :: read_sedqual ! read sediment porosity file

  ! Module variables
  character(len=bgc_fnmlen), public :: sedqualfile = ''

contains

  subroutine read_sedqual(kpie,kpje,ks,omask,sed_POCage_init,prorca_mavg_init)

    use mod_xc,             only: mnproc,xchalt
    use mo_kind,            only: rp
    use mo_control_bgc,     only: io_stdo_bgc,use_sediment_quality
    use netcdf,             only: nf90_noerr,nf90_nowrite,nf90_close,nf90_open
    use mo_netcdf_bgcrw,    only: read_netcdf_var

    implicit none

    ! Arguments
    integer, intent(in)    :: kpie
    integer, intent(in)    :: kpje
    integer, intent(in)    :: ks
    real(rp),intent(in)    :: omask(kpie,kpje)
    real(rp),intent(inout) :: sed_POCage_init(kpie,kpje,ks)
    real(rp),intent(inout) :: prorca_mavg_init(kpie,kpje)

    ! Local variables
    integer :: i,j,k
    real(rp):: sed_age(kpie,kpje,ks)
    real(rp):: mavg_prorca(kpie,kpje)
    logical :: file_exists = .false.
    integer :: ncid,ncstat

    ! Return if use_sediment_quality is turned off
    if (.not. use_sediment_quality) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*) ''
        write(io_stdo_bgc,*) 'read_sedqual: sediment quality-based remineralization is not activated'
      endif
      return
    endif

    ! Check if sediment quality file exists. If not, abort.
    inquire(file=sedqualfile,exist=file_exists)
    if (.not. file_exists .and. mnproc.eq.1) then
      write(io_stdo_bgc,*) ''
      write(io_stdo_bgc,*) 'read_sedqual: Cannot find sediment quality file... '
      call xchalt('(read_sedqual)')
      stop        '(read_sedqual)'
    endif

    ! read sediment quality initialization info from file
    if (mnproc.eq.1) then
      write(io_stdo_bgc,*) ''
      write(io_stdo_bgc,*) 'read_sedqual: read sediment quality information from ',trim(sedqualfile)
    endif

    ! Open netCDF data file
    if (mnproc==1) then
      ncstat = NF90_OPEN(trim(sedqualfile),NF90_NOWRITE, ncid)
      if (ncstat /= NF90_NOERR ) then
        call xchalt('(read_sedqual: Problem with netCDF1)')
        stop        '(read_sedqual: Problem with netCDF1)'
      end if
    end if

    ! Read  data
    call read_netcdf_var(ncid,'sedPOCage',sed_age(1,1,1),ks,0,0)
    call read_netcdf_var(ncid,'prorca_mavg',mavg_prorca,1,1,0)

    ! Close file
    if (mnproc==1) then
      ncstat = NF90_CLOSE(ncid)
      if ( ncstat /=  NF90_NOERR ) then
        call xchalt('(read_sedqual: Problem with netCDF2)')
        stop        '(read_sedqual: Problem with netCDF2)'
      end if
    end if

    do k=1,ks
      do j=1,kpje
        do i=1,kpie
          if(omask(i,j) > 0.5_rp)then
            sed_POCage_init(i,j,k) = sed_age(i,j,k)
          else
            sed_POCage_init(i,j,k) = 0._rp
          endif
        enddo
      enddo
    enddo

    do j=1,kpje
      do i=1,kpie
        if (omask(i,j) > 0.5_rp) then
          prorca_mavg_init(i,j) = mavg_prorca(i,j)
        else
          prorca_mavg_init(i,j) = 0._rp
        endif
      enddo
    enddo



  end subroutine read_sedqual

end module mo_read_sedqual
