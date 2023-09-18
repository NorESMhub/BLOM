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


module mo_read_ndep
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
!  J. Schwinger,     *NORCE climate, Bergen*   2022-06-02
!  -revise structure of this module, split into a module for reading the 
!   data (mo_read_ndep) and a module that applies the fluxes in core 
!   hamocc (mo_apply_ndep)
!
!
! Purpose
! -------
!  -Routines for reading nitrogen deposition fluxes from netcdf files
!
!
! Description:
! ------------
!
!  The routine get_ndep reads nitrogen deposition from file. The n-deposition
!  field is then passed to hamocc4bcm where it is applied to the top-most model 
!  layer by a call to apply_ndep (mo_apply_ndep).
!
!  N deposition is activated through a logical switch 'do_ndep' read from
!  HAMOCC's bgcnml namelist. If N deposition is acitvated, a valid filename
!  (including the full path) needs to be provided via HAMOCC's bgcnml namelist 
!  (variable ndepfile). If the input file is not found, an error will be issued.
!  The input data must be already pre-interpolated to the ocean grid.
!
!  -subroutine ini_read_ndep
!     Initialise the module
!
!  -subroutine get_ndep
!     Read and return n-deposition data for a given month.
!
!******************************************************************************
  implicit none

  private
  public :: ini_read_ndep,get_ndep,ndepfile

  character(len=512), save :: ndepfile=''
  real,  allocatable, save :: ndepread(:,:)
  real,  allocatable, save :: noydepread(:,:)
  real,  allocatable, save :: nhxdepread(:,:)
  integer,            save :: startyear,endyear  
  logical,            save :: lini = .false.

!******************************************************************************
contains



subroutine ini_read_ndep(kpie,kpje)
!******************************************************************************
!
!     S. Gao               *Gfi, Bergen*    19.08.2017 
!
! Purpose
! -------
!  -Initialise the module, check existence of input file, allocate array
!   for reading the data
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
  use mo_control_bgc, only: io_stdo_bgc,do_ndep,do_ndep_coupled
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
      write(io_stdo_bgc,*) 'ini_read_ndep: N deposition is not activated.'
    endif
    return
  end if
  if (do_ndep_coupled) then
    if (mnproc.eq.1) then
      write(io_stdo_bgc,*) ''
      write(io_stdo_bgc,*) 'ini_read_ndep: N deposition in interactive mode.'
    endif
    return
  end if

  ! Initialise the module
  if (.not. lini) then 

    IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'***************************************************'
      WRITE(io_stdo_bgc,*)'iHAMOCC: Initialization of module mo_read_ndep:'
      WRITE(io_stdo_bgc,*)' '
    ENDIF

    ! Check if nitrogen deposition file exists. If not, abort. 
    inquire(file=ndepfile,exist=file_exists)
    if (.not. file_exists .and. mnproc.eq.1) then 
      write(io_stdo_bgc,*) ''
      write(io_stdo_bgc,*) 'ini_read_ndep: Cannot find N deposition file... '
      call xchalt('(ini_read_ndep)')
      stop '(ini_read_ndep)' 
    endif 

#ifdef extNcycle
    ! Allocate field to hold N-deposition fluxes
    IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable nhxdepread ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    ENDIF
   
    ALLOCATE (nhxdepread(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory nhxdepread'
    nhxdepread(:,:) = 0.0

    IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable noydepread ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    ENDIF
   
    ALLOCATE (noydepread(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory noydepread'
    noydepread(:,:) = 0.0
#else
    ! Allocate field to hold N-deposition fluxes
    IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)'Memory allocation for variable ndepread ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    ENDIF
   
    ALLOCATE (ndepread(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory ndep'
    ndepread(:,:) = 0.0
#endif

    ! read start and end year of n-deposition file
    call ncfopn(trim(ndepfile),'r',' ',1,iotype)
    call ncgeti('startyear',startyear)
    call ncgeti('endyear',endyear)
    call ncfcls

    if (mnproc.eq.1) then 
      write(io_stdo_bgc,*) ''
      write(io_stdo_bgc,*) 'ini_read_ndep: Using N deposition file '//trim(ndepfile) 
    endif

    lini=.true.

  endif


!******************************************************************************
end subroutine ini_read_ndep


subroutine get_ndep(kpie,kpje,kplyear,kplmon,omask,ndep)
!******************************************************************************
!
!     S. Gao               *Gfi, Bergen*    19.08.2017 
!
! Purpose
! -------
!  -Read and return CMIP6 n-deposition data for a given month.
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
  use mo_param1_bgc,  only: nndep,idepnoy
#ifdef extNcycle
  use mo_param1_bgc, only: idepnhx
#endif  

  implicit none

  integer, intent(in)  :: kpie,kpje,kplyear,kplmon
  real,    intent(in)  :: omask(kpie,kpje)
  real,    intent(out) :: ndep(kpie,kpje,nndep)

  ! local variables 
  integer              :: month_in_file,ncstat,ncid,i,j
  integer, save        :: oldmonth=0 


  ! if N-deposition is switched off set ndep to zero and return
  if (.not. do_ndep) then
    ndep(:,:,:) = 0.0
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
#ifdef extNcycle
    call read_netcdf_var(ncid,'nhxdep',nhxdepread,1,month_in_file,0) 
    call read_netcdf_var(ncid,'noydep',noydepread,1,month_in_file,0) 
#else
    call read_netcdf_var(ncid,'ndep',ndepread,1,month_in_file,0) 
#endif
    ncstat=nf90_close(ncid)   
    oldmonth=kplmon 
  endif

!$OMP PARALLEL DO PRIVATE(i)
  ! 1 = NO3; 2 = NH4
  DO  j=1,kpje
  DO  i=1,kpie
#ifdef extNcycle
        ndep(i,j,idepnoy) = noydepread(i,j)  
        ndep(i,j,idepnhx) = nhxdepread(i,j)
#else
        ndep(i,j,idepnoy) = ndepread(i,j)  
#endif
  ENDDO
  ENDDO
!$OMP END PARALLEL DO 

!******************************************************************************
end subroutine get_ndep



!******************************************************************************
end module mo_read_ndep
