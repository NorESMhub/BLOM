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


module mo_read_rivin

  !*************************************************************************************************
  ! Routines for reading riverine nutrient and carbon input data
  !
  ! Riverine carbon and nutrient input is activated through a logical switch 'do_rivinpt' read
  ! from HAMOCC's bgcnml namelist. When coupled to NorESM, this is achieved by setting
  ! BLOM_RIVER_NUTRIENTS to TRUE in env_run.xml.
  !
  ! The model attempts to read nutrient fluxes from a NetCDF file
  ! derived from the GNEWS 2000 data base, which is specified through the
  ! namelist. The nutrient fluxes in the file are pre-interpolated to the
  ! ocean grid.
  !
  ! The nutrient discharge is distributed on the ocean grid in manner that is
  ! consistent with how model distributes its freshwater runoff.
  ! This has been achieved by using the mapping file used to interpolate the
  ! runoff also to interpolate the GNEWS nutrient fluxes to the ocean grid.
  !
  ! Since only alkalinity is available from measurements, DIC is updated using
  ! the assumtions that a_t=a_c+a_n and DIC=a_c (a_t: total alkalinity,
  ! a_c: carbonate alkalinity, a_n: contribution of nutrients to a_t).
  !
  ! S. Gao,              *Gfi, Bergen*    19.08.2017
  !
  ! Changes:
  !  J. Schwinger,     *NORCE climate, Bergen*   2020-05-27
  !   - re-structured this module such that riverine input can be passed as an
  !     argument to iHAMOCC's main routine
  !  J. Schwinger,     *NORCE climate, Bergen*   2022-05-18
  !   - re-structured and renamed this module such that reading and application of
  !     data are seperated into two distinct modules
  !*************************************************************************************************

  use dimensions,  only: idm,jdm
  use mod_xc ,     only: nbdy
  use mod_forcing, only: rivflx_stream ! here will holds input data as read from file

  implicit none
  private

  ! Routines
  public :: ini_read_rivin ! read gnews riverine nutrient and carbon data

  ! File name (incl. full path) for input data, set through namelist in mo_hamocc_init
  character(len=256), public :: rivinfile = ''

  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: riv_DIN2d
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: riv_DIP2d
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: riv_DSI2d
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: riv_DIC2d
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: riv_DFe2d
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: riv_idet2d
  real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: riv_idoc2d

contains

  subroutine ini_read_rivin(omask)
    !***********************************************************************************************
    !  Initialise reading of riverine input data (GNEWS 2000)
    !***********************************************************************************************

    use mod_xc,         only: mnproc,idm,jdm
    use mod_dia,        only: iotype
    use mod_nctools,    only: ncfopn,ncread,ncfcls
    use mo_control_bgc, only: io_stdo_bgc,do_rivinpt
    use mo_param1_bgc,  only: nriv,irdin,irdip,irsi,iralk,iriron,irdoc,irdet

    ! Arguments
    real, intent(in) :: omask(idm,jdm) ! ocean mask

    ! local variables
    integer :: i,j,errstat,dummymask(2)

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)' '
      write(io_stdo_bgc,*)'***************************************************'
      write(io_stdo_bgc,*)'iHAMOCC: Initialization of module mo_read_rivin:'
      write(io_stdo_bgc,*)' '
    endif

    ! Allocate field to hold river fluxes
    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable rivflx ...'
      write(io_stdo_bgc,*)'First dimension    : ',idm
      write(io_stdo_bgc,*)'Second dimension   : ',jdm
      write(io_stdo_bgc,*)'Third  dimension   : ',nriv
    endif

    allocate (rivflx_stream(idm,jdm,nriv),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory rivflx'
    rivflx_stream(:,:,:) = 0.0

    ! Return if riverine input is turned off
    if (.not. do_rivinpt) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*) ''
        write(io_stdo_bgc,*) 'ini_read_rivin: riverine input is not activated.'
      endif
      return
    endif

    ! read riverine nutrient fluxes from file
    if (mnproc.eq.1) then
      write(io_stdo_bgc,*) ''
      write(io_stdo_bgc,'(a)') 'ini_read_rivin: read riverine nutrients from ',trim(rivinfile)
    endif

    call ncfopn(trim(rivinfile),'r',' ',1,iotype)
    call ncread('DIN',riv_DIN2d,dummymask,0,0.)
    call ncread('DIP',riv_DIP2d,dummymask,0,0.)
    call ncread('DSi',riv_DSI2d,dummymask,0,0.)
    call ncread('DIC',riv_DIC2d,dummymask,0,0.) ! It is actually alkalinity that is observed
    call ncread('Fe' ,riv_DFe2d,dummymask,0,0.)
    call ncread('DOC',riv_idoc2d,dummymask,0,0.)
    call ncread('DET',riv_idet2d,dummymask,0,0.)
    call ncfcls

    do j=1,jdm
      do i=1,idm
        if (omask(i,j) > 0.5) then
          rivflx_stream(i,j,irdin)  = riv_DIN2d(i,j)
          rivflx_stream(i,j,irdip)  = riv_DIP2d(i,j)
          rivflx_stream(i,j,irsi)   = riv_DSI2d(i,j)
          rivflx_stream(i,j,iralk)  = riv_DIC2d(i,j)
          rivflx_stream(i,j,iriron) = riv_DFe2d(i,j)
          rivflx_stream(i,j,irdoc)  = riv_idoc2d(i,j)
          rivflx_stream(i,j,irdet)  = riv_idet2d(i,j)
        endif
      enddo
    enddo

  end subroutine ini_read_rivin

end module mo_read_rivin
