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


module mo_ndep
!********************************************************************************
!
!  S.Gao             *Gfi, Bergen*             2017-08-19
!
! Modified
! --------
!  J. Tjiputra,      *Uni Research, Bergen*    2017-09-18
!  -add 1 mol [H+], per mol [NO3] deposition, to alkalinity (minus 1 mol)
!
!  J. Schwinger,     *Uni Research, Bergen*    2018-04-12
!  -re-organised this module into an initialisation routine and a routine that 
!   does the deposition; introduced logical switch to activate N deposition.
!   
!
! Purpose
! -------
!  -Routines for reading and applying nitrogen deposition fluxes 
!
!
! Description:
! ------------
!  
!  The routine n_deposition reads nitrogen deposition from file and applies it 
!  to the top-most model layer. 
!
!  N deposition is activated through a logical switch 'do_ndep' read from 
!  HAMOCC's bgcnml namelist. If N deposition is acitvated, a valid filename
!  needs to be provided via HAMOCC's bgcnml namelist (variable ndepfname). If 
!  the input file is not found, an error will be issued.  
! 
!  The input data must be already pre-interpolated to the ocean grid and stored 
!  in the same folder with BLOM's grid information.   
!
!
!********************************************************************************
implicit none

character(len=512), save :: ndepfname='ndep_CMIP6_1850.nc'
character(len=512), save :: ndepfile
integer,            save :: startyear,endyear  
real, allocatable,  save :: ndep(:,:)
logical,            save :: lini = .false.

!********************************************************************************
contains



subroutine ini_ndep(kpie,kpje,path)
!********************************************************************************
!
!     S. Gao               *Gfi, Bergen*    19.08.2017 
!
! Purpose
! -------
!  -Initialise the n-deposition module.
!
! Changes: 
! --------
!
!********************************************************************************
use mod_xc,         only: mnproc,xchalt
use mo_control_bgc, only: io_stdo_bgc,do_ndep
use mod_dia,        only: iotype
use mod_nctools,    only: ncfopn,ncgeti,ncfcls

implicit none 

integer,          intent(in) :: kpie,kpje
character(len=*), intent(in) :: path

logical                      :: file_exists=.false.

! Return if N deposition is turned off
if (.not. do_ndep) then
  if (mnproc.eq.1) then
    write(io_stdo_bgc,*) ''
    write(io_stdo_bgc,*) 'ini_ndep: N deposition is not activated.'
  endif
  return
end if

! Initialise the module
if (.not. lini) then 

  ! Check if nitrogen deposition specified. If not, abort. 
  ndepfile=trim(path)//trim(ndepfname)
  inquire(file=ndepfile,exist=file_exists)
  if (.not. file_exists .and. mnproc.eq.1) then 
    write(io_stdo_bgc,*) ''
    write(io_stdo_bgc,*) 'ini_ndep: Cannot find N deposition file... '
    call xchalt('(ini_ndep)')
    stop '(ini_ndep)' 
  endif 

  ! read start and end year of n-deposition file
  call ncfopn(trim(ndepfile),'r',' ',1,iotype)
  call ncgeti('startyear',startyear)
  call ncgeti('endyear',endyear)
  call ncfcls

  ! allocate the field to hold one month of n-deposition data
  allocate(ndep(kpie,kpje))
  ndep(:,:) = 0.0
  if (mnproc.eq.1) then 
    write(io_stdo_bgc,*) ''
    write(io_stdo_bgc,*) 'ini_ndep: Using N deposition file '//trim(ndepfile) 
  endif

endif

lini=.true.


end



subroutine n_deposition(kpie,kpje,kpke,kplyear,kplmon,pddpo,omask)
!********************************************************************************
!
!     S. Gao               *Gfi, Bergen*    19.08.2017 
!
! Purpose
! -------
!  -Read n-deposition data and apply to the top-most model layer.
!
! Changes: 
! --------
!  Tjiputra (18.09.2017): add 1 mol [H+], per mol [NO3] deposition, to 
!    alkalinity (minus 1 mol)
!
!********************************************************************************
use mod_xc,         only: mnproc
use netcdf,         only: nf90_open,nf90_close,nf90_nowrite
use mo_control_bgc, only: io_stdo_bgc,dtb,do_ndep
use mo_carbch,      only: ocetra
#ifdef natDIC
use mo_param1_bgc,  only: iano3,ialkali,inatalkali
#else
use mo_param1_bgc,  only: iano3,ialkali
#endif

implicit none

integer :: kpie,kpje,kpke,kplyear,kplmon,kplday 
real, dimension(kpie,kpje,kpke) :: pddpo
real, dimension(kpie,kpje) :: omask

! local variables 
integer :: i,j,month_in_file,ncstat,ncid
integer, save :: oldmonth=0 

if (.not. do_ndep) return 

  
! read ndep data from file 
if (kplmon.ne.oldmonth) then 
  month_in_file=(max(startyear,min(endyear,kplyear))-startyear)*12+kplmon 
  if (mnproc.eq.1) then
    write(io_stdo_bgc,*) 'Read N deposition month ',month_in_file, & 
      & ' from file ',trim(ndepfile) 
  endif 
  ncstat=nf90_open(trim(ndepfile),nf90_nowrite,ncid)
  call read_netcdf_var(ncid,'ndep',ndep,1,month_in_file,0) 
  ncstat=nf90_close(ncid)   
  oldmonth=kplmon 
endif 

! deposite N in topmost layer 
do j=1,kpje
  do i=1,kpie
    if (omask(i,j).gt.0.5) then
      ocetra(i,j,1,iano3)=ocetra(i,j,1,iano3)+ndep(i,j)*dtb/365./pddpo(i,j,1)
      ocetra(i,j,1,ialkali)=ocetra(i,j,1,ialkali)-ndep(i,j)*dtb/365./pddpo(i,j,1)
#ifdef natDIC
      ocetra(i,j,1,inatalkali)=ocetra(i,j,1,inatalkali)-ndep(i,j)*dtb/365./pddpo(i,j,1)
#endif
   endif
  enddo
enddo

end subroutine n_deposition 

!********************************************************************************
end module mo_ndep 
