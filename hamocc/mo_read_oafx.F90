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
!    -'const_0p14':   constant alkalinity flux of 0.14 Pmol yr-1 applied to the 
!                     surface ocean between 60S and 70N (no input file needed)
!    -'const_0p56':   constant alkalinity flux of 0.56 Pmol yr-1 applied to the 
!                     surface ocean between 60S and 70N (no input file needed)
!    -'ramp':         ramping-up alkalinity flux from 0 Pmol yr-1 in 2025 to
!                     0.135 Pmol yr-1 in 2035 and onward, applied to the surface
!                     ocean between 60S and 70N (no input file needed)
!                     From G.Tran: 4279324154000 umol/s *3600 *24 *365 *1e-15
!                     *1e-6 = 0.135 Pmol yr-1
!
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
   
  character(len=128), save :: oalkscen=''
  character(len=512), save :: oalkfile=''
  real, parameter          :: Pmol2kmol  = 1.0e12
  
  ! Parameter used in the definition of alkalinization scenarios. The following 
  ! scenarios are defined in this module:
  !
  !  const_0p14    Homogeneous addition of 0.14 Pmol ALK/yr-1 over the ice-free
  !                surface ocean (assumed to be between 60S and 70N)
  !  const_0p56    Homogeneous addition of 0.56 Pmol ALK/yr-1 over the ice-free
  !                surface ocean (assumed to be between 60S and 70N)
  !  ramp          Linear increase of homogeneous addition of 0 to 0.135 Pmol
  !                ALK/yr-1 from 2025 to 2035 over the ice-free surface ocean  
  !                (assumed to be between 60S and 70N)
  !
  real, parameter :: addalk_0p14   = 0.14  ! Pmol alkalinity/yr added in the
  real, parameter :: addalk_0p56   = 0.56  ! 'const_0p14' and 'const_0p56' 
                                           ! scenarios
  real, parameter :: cdrmip_latmax =  70.0 ! Min and max latitude where
  real, parameter :: cdrmip_latmin = -60.0 ! alkalinity is added according
                                           ! to the CDRMIP protocol
  real, parameter :: addalk_ramp   = 0.135 ! Max Pmol alkalinity/yr added
  integer, parameter :: ramp_start = 2025  ! in 2035 in the 'ramp' scenario,
  integer, parameter :: ramp_end   = 2035  ! starting at 0 Pmol/yr in 2025.

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
  use mo_control_bgc, only: io_stdo_bgc,do_oalk

  implicit none 

  integer, intent(in) :: kpie,kpje
  real,    intent(in) :: pdlxp(kpie,kpje), pdlyp(kpie,kpje)
  real,    intent(in) :: pglat(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy)  
  real,    intent(in) :: omask(kpie,kpje)

  integer :: i,j,errstat
  real    :: avflx,ztotarea,addalk_tot
  real    :: ztmp1(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy)

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

    if( trim(oalkscen)=='const_0p14' .or. trim(oalkscen)=='const_0p56' .or.   &
        trim(oalkscen)=='ramp' ) then

      if(mnproc.eq.1) then
        write(io_stdo_bgc,*)'Using alkalinization scenario ', trim(oalkscen)
        write(io_stdo_bgc,*)' '
      endif

      ! Allocate field to hold constant alkalinization fluxes
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
      
      if( trim(oalkscen)=='const_0p14') then
        addalk_tot = addalk_0p14
      else if( trim(oalkscen)=='const_0p56') then
        addalk_tot = addalk_0p56
      else
        addalk_tot = addalk_ramp
      endif
    
      ! Calculate alkalinity flux (kmol m^2 yr-1) to be applied
      avflx = addalk_tot/ztotarea*Pmol2kmol
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
  if( trim(oalkscen)=='const_0p14' .or. trim(oalkscen)=='const_0p56' ) then
      
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
      current_day = (kplyear-ramp_start)*365+nday_of_year
      oafx(:,:) = oalkflx(:,:) * current_day / ((ramp_end-ramp_start)*365)
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
