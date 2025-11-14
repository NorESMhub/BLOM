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

module mo_read_sedpor

  !*************************************************************************************************
  ! Routine for reading sediment porosity from netcdf file L_SED_POR must be set to true in nml
  ! to activate lon-lat variable sediment porosity.
  !
  ! The model attempts to read lon-lat-sediment depth variable porosity
  ! from the input file 'SEDPORFILE' (incl. full path)
  !
  ! sed_por holds then the porosity that can be applied later in mo_sedmnt (ini_sedmnt_por)
  !*************************************************************************************************

  use mo_kind, only: bgc_fnmlen,rp

  implicit none
  private

  ! Routines
  public :: read_sedpor ! read sediment porosity file

  ! Module variables
  character(len=bgc_fnmlen), public :: sedporfile = ''

  real(rp), dimension(:,:,:), allocatable, public :: sed_por

contains

  subroutine read_sedpor(kpie,kpje,ks,omask)

    use mod_xc,             only: mnproc,xchalt
    use mo_kind,            only: rp
    use mo_control_bgc,     only: io_stdo_bgc,l_3Dvarsedpor
    use netcdf,             only: nf90_noerr,nf90_nowrite,nf90_close,nf90_open
    use mo_netcdf_bgcrw,    only: read_netcdf_var

    implicit none

    ! Arguments
    integer, intent(in)    :: kpie
    integer, intent(in)    :: kpje
    integer, intent(in)    :: ks
    real(rp),intent(in)    :: omask(kpie,kpje)

    !local variables
    integer :: i,j,k
    real(rp):: sed_por_in(kpie,kpje,ks)
    logical :: file_exists = .false.
    integer :: ncid,ncstat,errstat

    ! Return if l_3Dvarsedpor is turned off
    if (.not. l_3Dvarsedpor) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*) ''
        write(io_stdo_bgc,*) 'read_sedpor: spatially variable sediment porosity is not activated.'
      endif
      allocate (sed_por(1,1,1),stat=errstat) ! Dummy allocation to pass allocation check in hamocc_init
      if (mnproc==1 .and. errstat.ne.0) then
        write(io_stdo_bgc,*) 'Dummy memory allocation failed for sed_por in mo_read_sedpor'
      endif
      return
    endif

    allocate (sed_por(kpie,kpje,ks),stat=errstat)
    if (mnproc==1 .and. errstat.ne.0) then
      write(io_stdo_bgc,*) 'Memory allocation failed for sed_por in mo_read_sedpor'
    endif
    sed_por(:,:,:)         = 0._rp

    ! Check if sediment porosity file exists. If not, abort.
    inquire(file=sedporfile,exist=file_exists)
    if (.not. file_exists .and. mnproc.eq.1) then
      write(io_stdo_bgc,*) ''
      write(io_stdo_bgc,*) 'read_sedpor: Cannot find sediment porosity file... '
      call xchalt('(read_sedpor)')
      stop        '(read_sedpor)'
    endif

    ! read sediment porosity from file
    if (mnproc.eq.1) then
      write(io_stdo_bgc,*) ''
      write(io_stdo_bgc,*) 'read_sedpor: read sediment porosity from ',trim(sedporfile)
    endif

    ! Open netCDF data file
    if (mnproc==1) then
      ncstat = NF90_OPEN(trim(sedporfile),NF90_NOWRITE, ncid)
      if (ncstat /= NF90_NOERR ) then
        call xchalt('(read_sedpor: Problem with netCDF1)')
        stop        '(read_sedpor: Problem with netCDF1)'
      end if
    end if

    ! Read  data
    call read_netcdf_var(ncid,'sedpor',sed_por_in(1,1,1),ks,0,0)

    ! Close file
    if (mnproc==1) then
      ncstat = NF90_CLOSE(ncid)
      if ( ncstat /=  NF90_NOERR ) then
        call xchalt('(read_sedpor: Problem with netCDF2)')
        stop        '(read_sedpor: Problem with netCDF2)'
      end if
    end if

    do k=1,ks
      do j=1,kpje
        do i=1,kpie
          if(omask(i,j).gt. 0.5_rp)then
            sed_por(i,j,k)=sed_por_in(i,j,k)
          else
            sed_por(i,j,k)=0._rp
          endif
        enddo
      enddo
    enddo

  end subroutine read_sedpor

end module mo_read_sedpor
