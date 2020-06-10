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


module mo_riverinpt
!********************************************************************************
!
!     S. Gao,              *Gfi, Bergen*    19.08.2017
!
! Purpose
! -------
!  - Routines for reading and applying riverine nutrient and carbon input data
!
!
! Description:
! ------------
!  Public routines and variable of this module:
!
!  -subroutine ini_riverinpt
!    read gnews riverine nutrient and cabon data 
!
!  -subroutine riverinpt
!    apply riverine input to the ocean tracer fields
!
!  BLOM_RIVER_NUTRIENTS must be set to TRUE in env_run.xml to activate 
!  riverine nutrients.
!
!  The model attempts to read nutrient fluxes from the NetCDF file 
!  river_nutrients_GNEWS2000.nc, which has to be located in the same directory  
!  with BLOM's grid information. The nutrient fluxes in the file are 
!  pre-interpolated to the ocean grid. 
! 
!  The nutrient discharge is distributed on the ocean grid in manner that is 
!  consistent with how model distributes its freshwater runoff. 
!  This has been achieved by using the mapping file used to interpolate the 
!  runoff also to interpolate the GNEWS nutrient fluxes to the ocean grid.     
!  
!  Since only alkalinity is available from measurements, DIC is updated using
!  the assumtions that a_t=a_c+a_n and DIC=a_c (a_t: total alkalinity, 
!  a_c: carbonate alkalinity, a_n: contribution of nutrients to a_t).
!  
! Changes: 
! --------
!  J. Schwinger,     *NORCE climate, Bergen*   2020-05-27
!  - re-structured this moodule such that riverine input can be passed as an 
!    argument to iHAMOCC's main routine
! 
!********************************************************************************
use dimensions,     only: idm,jdm
use mod_xc ,        only: mnproc,nbdy

implicit none

private
public :: ini_riverinpt,riverinpt,nriv,rivflx,rivinfile
public :: irdin,irdip,irsi,iralk,iriron,irdoc,irdet

integer,         parameter :: nriv     = 7    ! size of river input field
integer,         parameter :: irdin    = 1, & ! dissolved inorganic nitrogen
                              irdip    = 2, & ! dissolved inorganic phosphorous 
                              irsi     = 3, & ! dissolved silicate 
                              iralk    = 4, & ! alkalinity
                              iriron   = 5, & ! dissolved bioavailable iron
                              irdoc    = 6, & ! dissolved organic carbon
                              irdet    = 7    ! particulate carbon
real,save,allocatable      :: rivflx(:,:,:)

! File name (incl. full path) for input data, set through namelist 
! in hamocc_init.F 
character(len=256),save    :: rivinfile = ''

! arrays for reading riverine inputs on the model grid
real,save,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: riv_DIN2d, riv_DIP2d,   &
                                                        riv_DSI2d, riv_DIC2d,   &
                                                        riv_idet2d,riv_idoc2d,  &
                                                        riv_DFe2d

!********************************************************************************
contains



subroutine ini_riverinpt(kpie,kpje,omask)
!--------------------------------------------------------------------------------
!
! Purpose:
! --------
!  Initialise riverine input to oceanic tracer fields
!
!
! Arguments:
! ----------
!  *INTEGER*     *kpie*    - 1st dimension of model grid.
!  *INTEGER*     *kpje*    - 2nd dimension of model grid.
!  *REAL*        *omask*   - ocean mask
!
!--------------------------------------------------------------------------------
  use mod_dia,        only: iotype
  use mod_nctools,    only: ncfopn,ncread,ncfcls 
  use mo_control_bgc, only: io_stdo_bgc,do_rivinpt

  implicit none

  integer,          intent(in) :: kpie,kpje
  real,             intent(in) :: omask(kpie,kpje)

  ! local variables
  integer :: i,j,errstat,dummymask(2)


  IF (mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)' '
    WRITE(io_stdo_bgc,*)'***************************************************'
    WRITE(io_stdo_bgc,*)'iHAMOCC: Initialization of module mo_riverinpt:'
    WRITE(io_stdo_bgc,*)' '
  ENDIF

  ! Allocate field to hold river fluxes
  IF (mnproc.eq.1) THEN
    WRITE(io_stdo_bgc,*)'Memory allocation for variable rivflx ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    WRITE(io_stdo_bgc,*)'Third  dimension   : ',nriv
  ENDIF

  ALLOCATE (rivflx(kpie,kpje,nriv),stat=errstat)
  if(errstat.ne.0) stop 'not enough memory rivflx'
  rivflx(:,:,:) = 0.0

  ! Return if riverine input is turned off
  if (.not. do_rivinpt) then
    if (mnproc.eq.1) then
      write(io_stdo_bgc,*) ''
      write(io_stdo_bgc,*) 'ini_riverinpt: riverine input is not activated.'
    endif
    return
  endif

  ! read riverine nutrient fluxes from file
  if (mnproc.eq.1) then
    write(io_stdo_bgc,*) ''
    write(io_stdo_bgc,*) 'ini_riverinpt: read riverine nutrients from ',       & 
                          trim(rivinfile)
  endif
  call ncfopn(trim(rivinfile),'r',' ',1,iotype)
  call ncread('DIN',riv_DIN2d,dummymask,0,0.)
  call ncread('DIP',riv_DIP2d,dummymask,0,0.)
  call ncread('DSi',riv_DSI2d,dummymask,0,0.)
  call ncread('DIC',riv_DIC2d,dummymask,0,0.) ! It is actually alkalinity that is observed
  call ncread('Fe',riv_DFe2d,dummymask,0,0.)
  call ncread('DOC',riv_idoc2d,dummymask,0,0.)
  call ncread('DET',riv_idet2d,dummymask,0,0.)
  call ncfcls



  DO j=1,kpje
  DO i=1,kpie
    IF(omask(i,j).GT.0.5) THEN

    rivflx(i,j,irdin)    = riv_DIN2d(i,j)
    rivflx(i,j,irdip)    = riv_DIP2d(i,j)
    rivflx(i,j,irsi)     = riv_DSI2d(i,j)
    rivflx(i,j,iralk)    = riv_DIC2d(i,j)
    rivflx(i,j,iriron)   = riv_DFe2d(i,j)
    rivflx(i,j,irdoc)    = riv_idoc2d(i,j)
    rivflx(i,j,irdet)    = riv_idet2d(i,j)

    ENDIF
  ENDDO
  ENDDO

!--------------------------------------------------------------------------------
end subroutine ini_riverinpt 



subroutine riverinpt(kpie,kpje,kpke,pddpo,omask,rivin)
!--------------------------------------------------------------------------------
!
! Purpose:
! --------
!  Apply riverine input to oceanic tracer fields
!
! Description:
! ------------
! 
!
! Arguments:
! ----------
!  *INTEGER*   *kpie*    - 1st dimension of model grid.
!  *INTEGER*   *kpje*    - 2nd dimension of model grid.
!  *INTEGER*   *kpke*    - 3rd (vertical) dimension of model grid.
!  *REAL*      *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!  *REAL*      *omask*   - ocean mask
!  *REAL*      *rivin*   - riverine input field [kmol m-2 yr-1]
!
!--------------------------------------------------------------------------------
  use mo_param1_bgc,  only: iano3,iphosph,isilica,isco212,iiron,idoc,idet,      &
                            ialkali,inatsco212,inatalkali
  use mo_control_bgc, only: dtb,do_rivinpt
  use mo_vgrid,       only: kmle
  use mo_carbch,      only: ocetra

  implicit none

  integer,intent(in) :: kpie,kpje,kpke
  real,   intent(in) :: pddpo(kpie,kpje,kpke)
  real,   intent(in) :: omask(kpie,kpje)
  real,   intent(in) :: rivin(kpie,kpje,nriv)

  ! local variables
  integer            :: i,j,k
  real               :: fdt,volij


  if (.not. do_rivinpt) return

!$OMP PARALLEL DO PRIVATE(fdt,volij)
  DO j=1,kpje
  DO i=1,kpie
    IF(omask(i,j).GT.0.5) THEN

      fdt   = dtb/365.

      ! Distribute riverine inputs over the model mixed layer
      volij = 0.
      DO k=1,kmle
        volij=volij+pddpo(i,j,k)
      ENDDO

      ! Notes:
      ! 1) DIC is updated using the assumtions that a_t=a_c+a_n and DIC=a_c (a_t: total 
      ! alkalinity, a_c: carbonate alkalinity, a_n: contribution of nutrients to a_t). 
      ! 2) Approx. 80-99% of dFe input is lost to the particulate phase in estuaries at 
      ! low salinities [Boyle et al., 1977; Chester, 1990; Dai and Martin, 1995; Lohan 
      ! and Bruland, 2006; Sholkovitz, 1978]. Below we assume 99% loss.
      ocetra(i,j,1:kmle,iano3)      = ocetra(i,j,1:kmle,iano3)      + rivin(i,j,irdin)*fdt/volij
      ocetra(i,j,1:kmle,iphosph)    = ocetra(i,j,1:kmle,iphosph)    + rivin(i,j,irdip)*fdt/volij
      ocetra(i,j,1:kmle,isilica)    = ocetra(i,j,1:kmle,isilica)    + rivin(i,j,irsi) *fdt/volij
      ocetra(i,j,1:kmle,isco212)    = ocetra(i,j,1:kmle,isco212)    + rivin(i,j,iralk)*fdt/volij  &
                                                                    + rivin(i,j,irdin)*fdt/volij  &
                                                                    + rivin(i,j,irdip)*fdt/volij
      ocetra(i,j,1:kmle,ialkali)    = ocetra(i,j,1:kmle,ialkali)    + rivin(i,j,iralk)*fdt/volij
#ifdef natDIC
      ocetra(i,j,1:kmle,inatsco212) = ocetra(i,j,1:kmle,inatsco212) + rivin(i,j,iralk)*fdt/volij  &
                                                                    + rivin(i,j,irdin)*fdt/volij  &
                                                                    + rivin(i,j,irdip)*fdt/volij
      ocetra(i,j,1:kmle,inatalkali) = ocetra(i,j,1:kmle,inatalkali) + rivin(i,j,iralk)*fdt/volij
#endif
      ocetra(i,j,1:kmle,iiron)      = ocetra(i,j,1:kmle,iiron)      + rivin(i,j,iriron)*fdt/volij*0.01
      ocetra(i,j,1:kmle,idoc)       = ocetra(i,j,1:kmle,idoc)       + rivin(i,j,irdoc)*fdt/volij
      ocetra(i,j,1:kmle,idet)       = ocetra(i,j,1:kmle,idet)       + rivin(i,j,irdet)*fdt/volij

  ENDIF
ENDDO
ENDDO
!$OMP END PARALLEL DO

!--------------------------------------------------------------------------------
end subroutine riverinpt



!********************************************************************************
end module mo_riverinpt
