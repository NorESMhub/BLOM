! Copyright (C) 2021-2022  J. Schwinger
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


module mo_read_oafx
!******************************************************************************
!
!   J.Schwinger             *NORCE Climate, Bergen*             2022-08-24
!
! Modified
! --------
!
! Purpose
! -------
!  -Routines for reading ocean alkalinization fluxes from netcdf files
!
!
! Description:
! ------------
!  The routine get_oafx reads a fluxs of alkalinity from file (or, for simple
!  cases, constructs an alkalinity flux field from scratch). The alkalinity
!  flux is then passed to hamocc4bcm where it is applied to the top-most model 
!  layer by a call to apply_oafx (mo_apply_oafx).
!
!  Ocean alkalinization is activated through a logical switch 'do_oalk' read from
!  HAMOCC's bgcnml namelist. If ocean alkalinization is acitvated, a valid 
!  name of an alkalinisation scenario (defined in this module, see below) and 
!  the file name (including the full path) of the corresponding OA-scenario 
!  input file needs to be provided via HAMOCC's bgcnml namelist (variables 
!  oascenario and oafxfile). If the input file is not found, an error will be 
!  issued. The input data must be already pre-interpolated to the ocean grid.
!
!  Currently available ocean alkalinisation scenarios:
!  (no input file needed, flux and latitude range can be defined in the
!   namelist, default values are defined):
!    -'const':        constant alkalinity flux applied to the surface ocean 
!                     between two latitudes.
!    -'ramp':         ramping-up alkalinity flux from 0 Pmol yr-1 to a maximum
!                     value between two specified years and kept constant
!                     onward, applied to the surface ocean between two
!                     latitudes.
!
!  -subroutine ini_read_oafx
!     Initialise the module
!
!  -subroutine get_oafx
!     Gets the alkalinity flux to apply at a given time.
!
!
!******************************************************************************
  implicit none

  private
  public :: ini_read_oafx,get_oafx,oalkscen,oalkfile

  real,allocatable, save :: oalkflx(:,:)
   
  character(len=128), save :: oalkscen   =''
  character(len=512), save :: oalkfile   =''
  real, parameter          :: Pmol2kmol  = 1.0e12
  
  ! Parameter used in the definition of alkalinization scenarios. The following 
  ! scenarios are defined in this module:
  !
  !  const         Constant homogeneous addition of alkalinity between latitude
  !                cdrmip_latmin and latitude cdrmip_latmax
  !  ramp          Linear increase of homogeneous addition from 0 to addalk
  !                Pmol ALK/yr-1 from year ramp_start to year ramp_end between
  !                latitude cdrmip_latmin and latitude cdrmip_latmax
  !
  real, protected :: addalk    = 0.56  ! Pmol alkalinity/yr added in the
                                           ! scenarios. Read from namelist file
                                           ! to overwrite default value.
  real, protected :: cdrmip_latmax =  70.0 ! Min and max latitude where
  real, protected :: cdrmip_latmin = -60.0 ! alkalinity is added according
                                           ! to the CDRMIP protocol. Read from
                                           ! namelist file to overwrite default
                                           ! value.
  integer, protected :: ramp_start = 2025  ! In 'ramp' scenario, start at
  integer, protected :: ramp_end   = 2035  ! 0 Pmol/yr in ramp_start, and max
                                           ! addalk Pmol/yr in ramp_end.
                                           ! Read from namelist file to
                                           ! overwrite default value.

  logical,   save :: lini = .false.

!******************************************************************************
contains



subroutine ini_read_oafx(kpie,kpje,pdlxp,pdlyp,pglat,omask)
!******************************************************************************
!
!     J.Schwinger               *NORCE Climate, Bergen*         2021-11-15
!
! Purpose
! -------
!  -Initialise the alkalinization module.
!
! Changes: 
! --------
!
! Parameter list:
! ---------------
!  *INTEGER* *kpie*       - 1st dimension of model grid.
!  *INTEGER* *kpje*       - 2nd dimension of model grid.
!  *REAL*    *pdlxp*      - size of grid cell (longitudinal) [m].
!  *REAL*    *pdlyp*      - size of grid cell (latitudinal) [m].
!  *REAL*    *pglat*      - latitude grid cell centres [degree N].
!  *REAL*    *omask*      - land/ocean mask.
!
!******************************************************************************
  use mod_xc,         only: xcsum,xchalt,mnproc,nbdy,ips
  use mo_control_bgc, only: io_stdo_bgc,do_oalk,bgc_namelist,get_bgc_namelist

  implicit none

  integer, intent(in) :: kpie,kpje
  real,    intent(in) :: pdlxp(kpie,kpje), pdlyp(kpie,kpje)
  real,    intent(in) :: pglat(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy)
  real,    intent(in) :: omask(kpie,kpje)

  integer :: i,j,errstat
  integer :: iounit
  real    :: avflx,ztotarea
  real    :: ztmp1(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy)

  namelist /bgcoafx/ do_oalk,oalkscen,oalkfile,addalk,cdrmip_latmax,          &
       &             cdrmip_latmin,ramp_start,ramp_end

  ! Read parameters for alkalinization fluxes from namelist file
  if(.not. allocated(bgc_namelist)) call get_bgc_namelist
  open (newunit=iounit, file=bgc_namelist, status='old'                   &
       &   ,action='read')
  read (unit=iounit, nml=BGCOAFX)
  close (unit=iounit)

  ! Return if alkalinization is turned off
  if (.not. do_oalk) then
    if (mnproc.eq.1) then
      write(io_stdo_bgc,*) ''
      write(io_stdo_bgc,*) 'ini_read_oafx: ocean alkalinization is not activated.'
    endif
    return
  end if

  ! Initialise the module
  if(.not. lini) then 

    if(mnproc.eq.1) then
      write(io_stdo_bgc,*)' '
      write(io_stdo_bgc,*)'***************************************************'
      write(io_stdo_bgc,*)'iHAMOCC: Initialization of module mo_read_oafx:'
      write(io_stdo_bgc,*)' '
    endif

    if( trim(oalkscen)=='const' .or. trim(oalkscen)=='ramp' ) then

      if(mnproc.eq.1) then
        write(io_stdo_bgc,*)'Using alkalinization scenario ', trim(oalkscen)
        write(io_stdo_bgc,*)' '
      endif

      ! Allocate field to hold alkalinization fluxes
      if(mnproc.eq.1) then
        write(io_stdo_bgc,*)'Memory allocation for variable oalkflx ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
      endif
   
      allocate(oalkflx(kpie,kpje),stat=errstat)
      if(errstat.ne.0) stop 'not enough memory oalkflx'
      oalkflx(:,:) = 0.0

      ! Calculate total ocean area 
      ztmp1(:,:)=0.0
      do j=1,kpje
      do i=1,kpie
        if( omask(i,j).gt.0.5 .and. pglat(i,j)<cdrmip_latmax                  &
                              .and. pglat(i,j)>cdrmip_latmin ) then
          ztmp1(i,j)=ztmp1(i,j)+pdlxp(i,j)*pdlyp(i,j)
        endif
      enddo
      enddo

      call xcsum(ztotarea,ztmp1,ips)
      
      ! Calculate alkalinity flux (kmol m^2 yr-1) to be applied
      avflx = addalk/ztotarea*Pmol2kmol
      if(mnproc.eq.1) then
        write(io_stdo_bgc,*)' '
        write(io_stdo_bgc,*)' applying alkalinity flux of ', avflx, ' kmol m-2 yr-1'
        write(io_stdo_bgc,*)'             over an area of ', ztotarea , ' m2'
        if( trim(oalkscen)=='ramp' ) then
          write(io_stdo_bgc,*)'             ramping-up from ', ramp_start, ' to ', ramp_end
        endif
      endif

      do j=1,kpje
      do i=1,kpie
        if( omask(i,j).gt.0.5 .and. pglat(i,j)<cdrmip_latmax                  &
                              .and. pglat(i,j)>cdrmip_latmin ) then
          oalkflx(i,j) = avflx
        endif
      enddo
      enddo

      lini=.true.

    !--------------------------------
    ! No valid scenario specified
    !--------------------------------
    else
    
      write(io_stdo_bgc,*) ''
      write(io_stdo_bgc,*) 'ini_read_oafx: invalid alkalinization scenario... '
      call xchalt('(ini_read_oafx)')
      stop '(ini_read_oafx)' 
      
    endif

  endif ! not lini


!******************************************************************************
end subroutine ini_read_oafx


subroutine get_oafx(kpie,kpje,kplyear,kplmon,omask,oafx)
!******************************************************************************
!
!     J. Schwinger            *NORCE Climate, Bergen*     2021-11-15
!
! Purpose
! -------
!  -return ocean alkalinization flux.
!
! Changes: 
! --------
!
!
! Parameter list:
! ---------------
!  *INTEGER*   *kpie*    - 1st dimension of model grid.
!  *INTEGER*   *kpje*    - 2nd dimension of model grid.
!  *INTEGER*   *kplyear* - current year.
!  *INTEGER*   *kplmon*  - current month.
!  *REAL*      *omask*   - land/ocean mask (1=ocean)
!  *REAL*      *oaflx*   - alkalinization flux [kmol m-2 yr-1]
!
!******************************************************************************
  use mod_xc,         only: xchalt,mnproc
  use mo_control_bgc, only: io_stdo_bgc,do_oalk
  use mod_time,       only: nday_of_year

  implicit none

  integer, intent(in)  :: kpie,kpje,kplyear,kplmon
  real,    intent(in)  :: omask(kpie,kpje)
  real,    intent(out) :: oafx(kpie,kpje)
  integer              :: current_day

  ! local variables 
  integer :: i,j

  if (.not. do_oalk) then
    oafx(:,:) = 0.0
    return 
  endif
  
  !--------------------------------
  ! Scenarios of constant fluxes
  !--------------------------------
  if( trim(oalkscen)=='const' ) then
      
    oafx(:,:) = oalkflx(:,:)

  !--------------------------------
  ! Scenario of ramping-up fluxes
  !--------------------------------
  elseif(trim(oalkscen)=='ramp' ) then

    if(kplyear.lt.ramp_start ) then
      oafx(:,:) = 0.0
    elseif(kplyear.ge.ramp_end ) then
      oafx(:,:) = oalkflx(:,:)
    else
      current_day = (kplyear-ramp_start)*365.+nday_of_year
      oafx(:,:) = oalkflx(:,:) * current_day / ((ramp_end-ramp_start)*365.)
    endif

    if(mnproc.eq.138 ) then
      write(io_stdo_bgc,*) 'get_oafx: oafx (kmol m-2 yr-1) ', oafx
    endif

  else
    
    write(io_stdo_bgc,*) ''
    write(io_stdo_bgc,*) 'get_oafx: invalid alkalinization scenario... '
    call xchalt('(get_oafx)')
    stop '(get_oafx)' 
      
  endif

!******************************************************************************
end subroutine get_oafx



!******************************************************************************
end module mo_read_oafx
