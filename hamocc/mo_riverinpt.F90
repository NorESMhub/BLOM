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
!     C. Bernard
!     J. Schwinger,        *Gfi, Bergen*    28.09.2012
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
! 
! Changes: 
! --------
!  19.08.2017
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
!********************************************************************************

use dimensions,     only: idm,jdm
use mod_xc ,        only: mnproc,nbdy

implicit none

private
public :: ini_riverinpt, riverinpt

character(len=*),parameter :: infile = 'river_nutrients_GNEWS2000.nc'

! arrays for riverine inputs on the model grid
real,save,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: riv_DIN2d, riv_DIP2d,   &
                                                        riv_DSI2d, riv_DIC2d,   &
                                                        riv_idet2d,riv_idoc2d,  &
                                                        riv_DFe2d

!********************************************************************************
contains



subroutine ini_riverinpt(path)
!--------------------------------------------------------------------------------
!
! Purpose:
! --------
!  Initialise riverine input to oceanic tracer fields
!
!--------------------------------------------------------------------------------
use mod_dia,        only: iotype
use mod_nctools,    only: ncfopn,ncread,ncfcls 
use mo_control_bgc, only: io_stdo_bgc,do_rivinpt

implicit none

character(len=*)   :: path 

! local variables
integer            :: dummymask(2)


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
  write(io_stdo_bgc,*) 'ini_riverinpt: read riverine nutrients from ',trim(path)//trim(infile)
endif
call ncfopn(trim(path)//trim(infile),'r',' ',1,iotype)
call ncread('DIN',riv_DIN2d,dummymask,0,0.)
call ncread('DIP',riv_DIP2d,dummymask,0,0.)
call ncread('DSi',riv_DSI2d,dummymask,0,0.)
call ncread('DIC',riv_DIC2d,dummymask,0,0.)
call ncread('Fe',riv_DFe2d,dummymask,0,0.)
call ncread('DOC',riv_idoc2d,dummymask,0,0.)
call ncread('DET',riv_idet2d,dummymask,0,0.)
call ncfcls

!--------------------------------------------------------------------------------
end subroutine ini_riverinpt 



subroutine riverinpt(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,omask)
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
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st dimension) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (2nd dimension) [m].
!     *REAL*    *omask*   - ocean mask
!
!--------------------------------------------------------------------------------
#ifdef natDIC
use mo_param1_bgc,  only: iano3,iphosph,isilica,isco212,iiron,idoc,idet,ialkali,&
                          inatsco212,inatalkali
#else
use mo_param1_bgc,  only: iano3,iphosph,isilica,isco212,iiron,idoc,idet,ialkali
#endif
use mo_control_bgc, only: dtb,do_rivinpt
use mo_vgrid,       only: kmle
use mo_carbch,      only: ocetra

implicit none

integer,intent(in) :: kpie,kpje,kpke
real,   intent(in) :: pddpo(kpie,kpje,kpke)
real,   intent(in) :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
real,   intent(in) :: omask(kpie,kpje)

! local variables
integer            :: i,j,k,dummymask(2)
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

    ! Inorganic elements
    ocetra(i,j,1:kmle,iano3)      = ocetra(i,j,1:kmle,iano3)      + riv_DIN2d(i,j)*fdt/volij
    ocetra(i,j,1:kmle,iphosph)    = ocetra(i,j,1:kmle,iphosph)    + riv_DIP2d(i,j)*fdt/volij
    ocetra(i,j,1:kmle,isilica)    = ocetra(i,j,1:kmle,isilica)    + riv_DSI2d(i,j)*fdt/volij
    ocetra(i,j,1:kmle,isco212)    = ocetra(i,j,1:kmle,isco212)    + riv_DIC2d(i,j)*fdt/volij  &
                                                                  + riv_DIN2d(i,j)*fdt/volij  &
                                                                  + riv_DIP2d(i,j)*fdt/volij
    ocetra(i,j,1:kmle,ialkali)    = ocetra(i,j,1:kmle,ialkali)    + riv_DIC2d(i,j)*fdt/volij
#ifdef natDIC
    ocetra(i,j,1:kmle,inatsco212) = ocetra(i,j,1:kmle,inatsco212) + riv_DIC2d(i,j)*fdt/volij  &
                                                                  + riv_DIN2d(i,j)*fdt/volij  &
                                                                  + riv_DIP2d(i,j)*fdt/volij
    ocetra(i,j,1:kmle,inatalkali) = ocetra(i,j,1:kmle,inatalkali) + riv_DIC2d(i,j)*fdt/volij
#endif
    ocetra(i,j,1:kmle,iiron)      = ocetra(i,j,1:kmle,iiron)      + riv_DFe2d(i,j)*fdt/volij*0.01
!SG: Approx. 80-99% of dFe input is lost to the particulate phase in estuaries at low salinities 
!    [Boyle et al., 1977; Chester, 1990; Dai and Martin, 1995; Lohan and Bruland, 2006; Sholkovitz, 1978] 

    ! Dissolved organic matter
    ocetra(i,j,1:kmle,idoc)       = ocetra(i,j,1:kmle,idoc)       + riv_idoc2d(i,j)*fdt/volij

    ! Particulate organic matter
    ocetra(i,j,1:kmle,idet)       = ocetra(i,j,1:kmle,idet)       + riv_idet2d(i,j)*fdt/volij

  ENDIF
ENDDO
ENDDO
!$OMP END PARALLEL DO

!--------------------------------------------------------------------------------
end subroutine riverinpt



!********************************************************************************
end module mo_riverinpt
