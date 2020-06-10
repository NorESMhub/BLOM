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
!******************************************************************************
!
!   S.Gao             *Gfi, Bergen*             2017-08-19
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
!  J. Schwinger,     *NORCE climate, Bergen*   2020-05-27
!  -put reading of a time-slice of n-deposition data into own subroutine
!  -removed default file name
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
!  needs to be provided via HAMOCC's bgcnml namelist (variable ndepfile). If 
!  the input file is not found, an error will be issued.  
! 
!  The input data must be already pre-interpolated to the ocean grid and stored 
!  in the same folder with BLOM's grid information.   
!
!
!******************************************************************************
  implicit none

  private
  public :: ini_ndep,get_ndep,n_deposition,ndepfile

  character(len=512), save :: ndepfile=''
  real,  allocatable, save :: ndepread(:,:)
  integer,            save :: startyear,endyear  
  logical,            save :: lini = .false.

!******************************************************************************
contains



subroutine ini_ndep(kpie,kpje)
!******************************************************************************
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
! Parameter list:
! ---------------
!  *INTEGER* *kpie*       - 1st dimension of model grid.
!  *INTEGER* *kpje*       - 2nd dimension of model grid.
!
!******************************************************************************
  use mod_xc,         only: mnproc,xchalt
  use mo_control_bgc, only: io_stdo_bgc,do_ndep
  use mod_dia,        only: iotype
  use mod_nctools,    only: ncfopn,ncgeti,ncfcls

  implicit none 

  integer, intent(in) :: kpie,kpje

  integer :: errstat
  logical :: file_exists=.false.


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

    ! Check if nitrogen deposition file exists. If not, abort. 
    inquire(file=ndepfile,exist=file_exists)
    if (.not. file_exists .and. mnproc.eq.1) then 
      write(io_stdo_bgc,*) ''
      write(io_stdo_bgc,*) 'ini_ndep: Cannot find N deposition file... '
      call xchalt('(ini_ndep)')
      stop '(ini_ndep)' 
    endif 

    IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'***************************************************'
      WRITE(io_stdo_bgc,*)'iHAMOCC: Initialization of module mo_ndep:'
      WRITE(io_stdo_bgc,*)' '
    ENDIF

    ! Allocate field to hold N-deposition fluxes
    IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable ndepread ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    ENDIF
   
    ALLOCATE (ndepread(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory ndep'
    ndepread(:,:) = 0.0

    ! read start and end year of n-deposition file
    call ncfopn(trim(ndepfile),'r',' ',1,iotype)
    call ncgeti('startyear',startyear)
    call ncgeti('endyear',endyear)
    call ncfcls

    if (mnproc.eq.1) then 
      write(io_stdo_bgc,*) ''
      write(io_stdo_bgc,*) 'ini_ndep: Using N deposition file '//trim(ndepfile) 
    endif

    lini=.true.

  endif


!******************************************************************************
end subroutine ini_ndep


subroutine get_ndep(kpie,kpje,kplyear,kplmon,omask,ndep)
!******************************************************************************
!
!     S. Gao               *Gfi, Bergen*    19.08.2017 
!
! Purpose
! -------
!  -Read and return n-deposition data for a given month.
!
! Parameter list:
! ---------------
!  *INTEGER*   *kpie*    - 1st dimension of model grid.
!  *INTEGER*   *kpje*    - 2nd dimension of model grid.
!  *INTEGER*   *kplyear* - current year.
!  *INTEGER*   *kplmon*  - current month.
!  *REAL*      *omask*   - land/ocean mask (1=ocean)
!  *REAL*      *ndep*    - N-deposition field for current year and month
!
!******************************************************************************
  use mod_xc,         only: mnproc
  use netcdf,         only: nf90_open,nf90_close,nf90_nowrite
  use mo_control_bgc, only: io_stdo_bgc,do_ndep

  implicit none

  integer, intent(in)  :: kpie,kpje,kplyear,kplmon
  real,    intent(in)  :: omask(kpie,kpje)
  real,    intent(out) :: ndep(kpie,kpje)

  ! local variables 
  integer              :: month_in_file,ncstat,ncid
  integer, save        :: oldmonth=0 


  ! if N-deposition is switched off set ndep to zero and return
  if (.not. do_ndep) then
    ndep(:,:) = 0.0
    return 
  endif


  ! read ndep data from file 
  if (kplmon.ne.oldmonth) then 
    month_in_file=(max(startyear,min(endyear,kplyear))-startyear)*12+kplmon 
    if (mnproc.eq.1) then
      write(io_stdo_bgc,*) 'Read N deposition month ',month_in_file, & 
                           ' from file ',trim(ndepfile) 
    endif 
    ncstat=nf90_open(trim(ndepfile),nf90_nowrite,ncid)
    call read_netcdf_var(ncid,'ndep',ndepread,1,month_in_file,0) 
    ncstat=nf90_close(ncid)   
    oldmonth=kplmon 
  endif

  ndep(:,:) = ndepread

!******************************************************************************
end subroutine get_ndep



subroutine n_deposition(kpie,kpje,kpke,pddpo,omask,ndep)
!******************************************************************************
!
!     S. Gao               *Gfi, Bergen*    19.08.2017 
!
! Purpose
! -------
!  -apply n-deposition to the top-most model layer.
!
! Changes: 
! --------
!  Tjiputra (18.09.2017): add 1 mol [H+], per mol [NO3] deposition, to 
!    alkalinity (minus 1 mol)
!
! Parameter list:
! ---------------
!  *INTEGER*   *kpie*    - 1st dimension of model grid.
!  *INTEGER*   *kpje*    - 2nd dimension of model grid.
!  *INTEGER*   *kpke*    - 3rd (vertical) dimension of model grid.
!  *REAL*      *pddpo*   - size of grid cell (depth) [m].
!  *REAL*      *omask*   - land/ocean mask (1=ocean)
!  *REAL*      *ndep*    - N-deposition field to apply
!
!******************************************************************************
  use mod_xc,         only: mnproc
  use mo_control_bgc, only: io_stdo_bgc,dtb,do_ndep
  use mo_carbch,      only: ocetra
  use mo_param1_bgc,  only: iano3,ialkali,inatalkali

  implicit none

  integer, intent(in) :: kpie,kpje,kpke
  real,    intent(in) :: pddpo(kpie,kpje,kpke)
  real,    intent(in) :: omask(kpie,kpje)
  real,    intent(in) :: ndep(kpie,kpje)

  ! local variables 
  integer :: i,j

  if (.not. do_ndep) return 

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

!******************************************************************************
end subroutine n_deposition 


!******************************************************************************
end module mo_ndep 
