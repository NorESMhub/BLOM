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

  !*************************************************************************************************
  !  Routines for reading nitrogen deposition fluxes from netcdf files
  !
  !  N-deposition is activated through a logical switch 'do_ndep' read from HAMOCC's bgcnml
  !  namelist. When coupled to NorESM, this is achieved by setting BLOM_N_DEPOSITION to
  !  TRUE in env_run.xml.
  !
  !  The routine get_ndep reads nitrogen deposition from file. The n-deposition
  !  field is then passed to hamocc4bcm where it is applied to the top-most model
  !  layer by a call to apply_ndep (mo_apply_ndep). If N deposition is acitvated, a
  !  valid filename (including the full path) needs to be provided via HAMOCC's bgcnml
  !  namelist (variable ndepfile). If the input file is not found, an error will be issued.
  !  The input data must be already pre-interpolated to the ocean grid.
  !
  ! S.Gao             *Gfi, Bergen*             2017-08-19
  !
  ! Modified
  !  J. Tjiputra,      *Uni Research, Bergen*    2017-09-18
  !   - add 1 mol [H+], per mol [NO3] deposition, to alkalinity (minus 1 mol)
  !  J. Schwinger,     *Uni Research, Bergen*    2018-04-12
  !   - re-organised this module into an initialisation routine and a routine that
  !     does the deposition; introduced logical switch to activate N deposition.
  !  J. Schwinger,     *NORCE climate, Bergen*   2020-05-27
  !   - put reading of a time-slice of n-deposition data into own subroutine
  !   - removed default file name
  !  J. Schwinger,     *NORCE climate, Bergen*   2022-06-02
  !   - revise structure of this module, split into a module for reading the data (mo_read_ndep)
  !     and a module that applies the fluxes in core hamocc (mo_apply_ndep)
  !
  !*************************************************************************************************

  implicit none
  private

  public :: ini_read_ndep ! Initialise the module
  public :: get_ndep      ! Read and return n-deposition data for a given month.
  public :: ndepfile

  character(len=512)  :: ndepfile=''
  real,  allocatable  :: noydepread(:,:)
  real,  allocatable  :: nhxdepread(:,:)
  integer             :: startyear,endyear
  logical             :: lini = .false.
  integer             :: oldmonth=0

contains

  subroutine ini_read_ndep(kpie,kpje)

    !***********************************************************************************************
    ! Initialise the module, check existence of input file, allocate array for reading the data
    !
    ! S. Gao   *Gfi, Bergen*    19.08.2017
    !***********************************************************************************************

    use mod_xc,             only: mnproc,xchalt
    use mo_control_bgc,     only: io_stdo_bgc,do_ndep,use_extNcycle
    use mod_dia,            only: iotype
    use mod_nctools,        only: ncfopn,ncgeti,ncfcls
    use mo_netcdf_bgcrw,    only: read_netcdf_var

    ! Arguments
    integer, intent(in) :: kpie ! 1st dimension of model grid.
    integer, intent(in) :: kpje ! 2nd dimension of model grid.

    ! Local variables
    integer :: errstat
    logical :: file_exists=.false.

    ! Return if N deposition is turned off
    if (.not. do_ndep) then
      if (mnproc == 1) then
        write(io_stdo_bgc,*) ''
        write(io_stdo_bgc,*) 'ini_read_ndep: N deposition is not activated.'
      endif
      return
    endif

    ! Initialise the module
    if (.not. lini) then

      if (mnproc == 1) then
        write(io_stdo_bgc,*)' '
        write(io_stdo_bgc,*)'***************************************************'
        write(io_stdo_bgc,*)'iHAMOCC: Initialization of module mo_read_ndep:'
        write(io_stdo_bgc,*)' '
      endif

      ! Check if nitrogen deposition file exists. If not, abort.
      inquire(file=ndepfile,exist=file_exists)
      if (.not. file_exists .and. mnproc == 1) then
        write(io_stdo_bgc,*) ''
        write(io_stdo_bgc,*) 'ini_read_ndep: Cannot find N deposition file... '
        call xchalt('(ini_read_ndep)')
        stop '(ini_read_ndep)'
      endif

      ! Allocate field to hold N-deposition fluxes
      if (mnproc == 1) then
        write(io_stdo_bgc,*)'Memory allocation for variable nhxdepread ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif
      allocate (nhxdepread(kpie,kpje),stat=errstat)
      if(errstat /= 0) stop 'not enough memory nhxdepread'
      nhxdepread(:,:) = 0.0

      if (mnproc == 1) then
        write(io_stdo_bgc,*)'Memory allocation for variable noydepread ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif
      allocate (noydepread(kpie,kpje),stat=errstat)
      if(errstat /= 0) stop 'not enough memory noydepread'
      noydepread(:,:) = 0.0

      ! read start and end year of n-deposition file
      call ncfopn(trim(ndepfile),'r',' ',1,iotype)
      call ncgeti('startyear',startyear)
      call ncgeti('endyear',endyear)
      call ncfcls

      if (mnproc == 1) then
        write(io_stdo_bgc,*) ''
        write(io_stdo_bgc,*) 'ini_read_ndep: Using N deposition file '//trim(ndepfile)
      endif

      lini=.true.

    endif

  end subroutine ini_read_ndep


  subroutine get_ndep(kpie,kpje,kbnd,kplyear,kplmon,omask,ndep,patmnhxdep,patmnoydep)

    !***********************************************************************************************
    ! Read and return CMIP6 n-deposition data for a given month or use atmosphere input
    !
    ! S. Gao               *Gfi, Bergen*    19.08.2017
    !***********************************************************************************************

    use mod_xc,             only: mnproc
    use netcdf,             only: nf90_open,nf90_close,nf90_nowrite
    use mo_control_bgc,     only: io_stdo_bgc, do_ndep, use_extNcycle, use_coupler_ndep
    use mo_netcdf_bgcrw,    only: read_netcdf_var
    use mo_param1_bgc,      only: nndep,idepnoy,idepnhx
    use mo_chemcon,         only: mw_nitrogen

    ! Arguments
    integer, intent(in)  :: kpie              ! 1st dimension of model grid.
    integer, intent(in)  :: kpje              ! 2nd dimension of model grid.
    integer, intent(in)  :: kbnd              !
    integer, intent(in)  :: kplyear           ! current year.
    integer, intent(in)  :: kplmon            ! current month.
    real,    intent(in)  :: omask(kpie,kpje)  ! land/ocean mask (1=ocean)
    real,    intent(out) :: ndep(kpie,kpje,nndep) ! N-deposition field for current year and month
    real,    intent(in)  :: patmnhxdep(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)   ! Atmospheric NHx deposition [kgN m-2 s-1]
    real,    intent(in)  :: patmnoydep(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)   ! Atmospheric NOy deposition [kgN m-2 s-1]


    ! local variables
    integer  :: month_in_file, ncstat, ncid, i, j
    real     :: fatmndep
    logical  :: first_call = .true.

    ndep(:,:,:) = 0.0

    if (.not. do_ndep) then
      ! if N-deposition is switched off return
      return
    endif

    if (use_coupler_ndep) then
      ! If  use_coupler_ndep, nitrogen deposition is ALWAYS obtained from the
      ! mct coupler or nuopc mediator (CAM or CDEPS)
      if (mnproc == 1 .and. first_call) then
        write (io_stdo_bgc,*) 'iHAMOCC: getting NOy and NHx deposition from atm'
      endif

      ! convert from kgN/m2/s to climatological input file units: kmolN/m2/yr
      fatmndep = 365.*86400./mw_nitrogen

      if (use_extNcycle) then
        !$omp parallel do private(i)
        do j=1,kpje
          do i=1,kpie
            if (patmnoydep(i,j) > 0.) then
              ndep(i,j,idepnoy) = patmnoydep(i,j)*fatmndep
            endif
            if (patmnhxdep(i,j) > 0.) then
              ndep(i,j,idepnhx) = patmnhxdep(i,j)*fatmndep
            endif
          enddo
        enddo
        !$omp end parallel do
      else
        !$omp parallel do private(i)
        do j=1,kpje
          do i=1,kpie
            if (patmnoydep(i,j) > 0. .and.  patmnhxdep(i,j) > 0.) then
              ! reduced and oxidized forms all enter the NO3 pool
              ndep(i,j,idepnoy) = (patmnoydep(i,j)+patmnhxdep(i,j))*fatmndep
            endif
          enddo
        enddo
        !$omp end parallel do
      end if

    else
      ! NOTE: No online coupling of the extended nitrogen cycling through MCT
      !       - only input files are read!
      ! read ndep data from file
      if (kplmon /= oldmonth) then
        month_in_file=(max(startyear,min(endyear,kplyear))-startyear)*12+kplmon
        if (mnproc == 1) then
          write(io_stdo_bgc,*) 'Read N deposition month ',month_in_file,' from file ',trim(ndepfile)
        endif
        ncstat=nf90_open(trim(ndepfile),nf90_nowrite,ncid)
        call read_netcdf_var(ncid,'nhxdep',nhxdepread,1,month_in_file,0)
        call read_netcdf_var(ncid,'noydep',noydepread,1,month_in_file,0)
        ncstat=nf90_close(ncid)
        oldmonth=kplmon
      endif

      !$omp parallel do private(i)
      ! 1 = NO3; 2 = NH4
      do j=1,kpje
        do i=1,kpie
          if (use_extNcycle) then
            ndep(i,j,idepnoy) = noydepread(i,j)
            ndep(i,j,idepnhx) = nhxdepread(i,j)
          else
            ndep(i,j,idepnoy) = noydepread(i,j) + nhxdepread(i,j)
          endif
        enddo
      enddo
      !$omp end parallel do

    endif
    first_call = .false.
  end subroutine get_ndep

end module mo_read_ndep
