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
!*****************************************************************************
! Purpose
! -------
!   - Routine for reading sediment porosity from netcdf file
!
! Description
! -----------
! Public routines and variable of this module:
!
!   - subroutine ini_read_sedpor
!        read sediment porosity file
!
!   L_SED_POR must be set to true in nml to activate 
!   lon-lat variable sediment porosity.
!
!   The model attempts to read lon-lat-sediment depth variable porosity
!   from the input file 'SEDPORFILE' (incl. full path)
!  
!   sed_por holds then the porosity that can be applied later 
!   via mo_apply_sedpor
!
!*****************************************************************************

implicit none

private

public :: read_sedpor,sedporfile

character(len=512),save :: sedporfile = ''

contains

subroutine read_sedpor(kpie,kpje,ks,omask,sed_por)
  use mod_xc,         only: mnproc,xchalt
  use mod_dia,        only: iotype
  use mo_control_bgc, only: io_stdo_bgc,l_3Dvarsedpor
  use mod_nctools,    only: ncfopn,ncread,ncfcls

  implicit none

  integer, intent(in) :: kpie,kpje,ks
  real,    intent(in) :: omask(kpie,kpje)
  real,    intent(inout) :: sed_por(kpie,kpje,ks)

  !local variables
  integer :: i,j,k,errstat,dummymask(2)
  real    :: sed_por_in(kpie,kpje,ks)
  logical :: file_exists = .false.

  ! Return if l_3Dvarsedpor is turned off
  if (.not. l_3Dvarsedpor) then
    if (mnproc.eq.1) then
      write(io_stdo_bgc,*) ''
      write(io_stdo_bgc,*) 'read_sedpor: spatially variable sediment porosity is not activated.'
    endif
    return
  endif

  ! Check if sediment porosity file exists. If not, abort. 
  inquire(file=sedporfile,exist=file_exists)
  if (.not. file_exists .and. mnproc.eq.1) then 
    write(io_stdo_bgc,*) ''
    write(io_stdo_bgc,*) 'read_sedpor: Cannot find sediment porosity file... '
    call xchalt('(read_sedpor)')
    stop '(read_sedpor)' 
  endif 

  ! read sediment porosity from file
  if (mnproc.eq.1) then
    write(io_stdo_bgc,*) ''
    write(io_stdo_bgc,*) 'read_sedpor: read sediment porosity from ',       & 
                          trim(sedporfile)
  endif
  call ncfopn(trim(sedporfile),'r',' ',1,iotype)
  call ncread('sedpor',sed_por_in,dummymask,0,0.)
  call ncfcls

  do k=1,ks
  do j=1,kpje
  do i=1,kpie
     if(omask(i,j).gt. 0.5)then
       sed_por(i,j,k)=sed_por_in(i,j,k)
     endif
  enddo
  enddo
  enddo

end subroutine read_sedpor
end module mo_read_sedpor
