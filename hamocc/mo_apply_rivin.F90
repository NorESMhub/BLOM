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


module mo_apply_rivin
!********************************************************************************
!
!     S. Gao,              *Gfi, Bergen*    19.08.2017
!
! Purpose
! -------
!  - Routines for applying riverine nutrient and carbon input data
!
!
! Description:
! ------------
!  Public routines and variable of this module:
!
!  -subroutine apply_rivin
!    apply riverine input to the ocean tracer field
!
!  BLOM_RIVER_NUTRIENTS must be set to TRUE in env_run.xml to activate 
!  riverine nutrients.
!
!  
! Changes: 
! --------
!  J. Schwinger,     *NORCE climate, Bergen*   2020-05-27
!  - re-structured this module such that riverine input can be passed as an 
!    argument to iHAMOCC's main routine
! 
!  J. Schwinger,     *NORCE climate, Bergen*   2022-05-18
!  - re-structured and renamed this module such that reading and application of 
!    data are seperated into two distinct modules
!
!********************************************************************************
implicit none

private
public :: apply_rivin

! Approx. 80-99% of dFe riverine input is lost to the particulate phase in 
! estuaries at low salinities [Boyle et al., 1977; Chester, 1990; Dai and 
! Martin, 1995; Lohan and Bruland, 2006; Sholkovitz, 1978]. dFe_frac is the
! fraction of dissolved iron that enters the costal ocean.
real, parameter :: dFe_frac = 0.01  ! assume 99% loss of dissolved iron


!********************************************************************************
contains



subroutine apply_rivin(kpie,kpje,kpke,pddpo,omask,rivin)
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
  use mo_control_bgc, only: dtb,do_rivinpt
  use mo_param1_bgc,  only: nriv,irdin,irdip,irsi,iralk,iriron,irdoc,irdet,     &
                            iano3,iphosph,isilica,isco212,iiron,idoc,idet,      &
                            ialkali,inatsco212,inatalkali
#ifdef cisonew
  use mo_param1_bgc,  only: idet13,idet14,idoc13,idoc14,isco213,isco214,safediv
#endif

  use mo_vgrid,       only: kmle
  use mo_carbch,      only: ocetra,rivinflx

  implicit none

  integer,intent(in) :: kpie,kpje,kpke
  real,   intent(in) :: pddpo(kpie,kpje,kpke)
  real,   intent(in) :: omask(kpie,kpje)
  real,   intent(in) :: rivin(kpie,kpje,nriv)

  ! local variables
  integer            :: i,j,k
  real               :: fdt,volij


  ! rivinflx stores the applied n-deposition flux for inventory calculations 
  ! and output
  rivinflx(:,:,:) = 0.0

  if (.not. do_rivinpt) return
  
  fdt = dtb/365.

!$OMP PARALLEL DO PRIVATE(i,k,volij)
  DO j=1,kpje
  DO i=1,kpie
    IF(omask(i,j).GT.0.5) THEN

      ! Distribute riverine inputs over the model mixed layer
      volij = 0.
      DO k=1,kmle(i,j)
        volij=volij+pddpo(i,j,k)
      ENDDO

#ifdef cisonew
      ocetra(i,j,1:kmle(i,j),isco213)    = ocetra(i,j,1:kmle(i,j),isco213)    +                   &
            ocetra(i,j,1:kmle(i,j),isco213)/(ocetra(i,j,1:kmle(i,j),isco212)+safediv)*            &
            (rivin(i,j,iralk)*fdt/volij + rivin(i,j,irdin)*fdt/volij + rivin(i,j,irdip)*fdt/volij)
      ocetra(i,j,1:kmle(i,j),isco214)    = ocetra(i,j,1:kmle(i,j),isco214)    +                   &
            ocetra(i,j,1:kmle(i,j),isco214)/(ocetra(i,j,1:kmle(i,j),isco212)+safediv)*            &
            (rivin(i,j,iralk)*fdt/volij + rivin(i,j,irdin)*fdt/volij + rivin(i,j,irdip)*fdt/volij)

      ocetra(i,j,1:kmle(i,j),idoc13)     = ocetra(i,j,1:kmle(i,j),idoc13)     +                   &
            ocetra(i,j,1:kmle(i,j),idoc13)/(ocetra(i,j,1:kmle(i,j),idoc)+safediv)*                &
            rivin(i,j,irdoc)*fdt/volij
      ocetra(i,j,1:kmle(i,j),idoc14)     = ocetra(i,j,1:kmle(i,j),idoc14)     +                   &
            ocetra(i,j,1:kmle(i,j),idoc14)/(ocetra(i,j,1:kmle(i,j),idoc)+safediv)*                &
            rivin(i,j,irdoc)*fdt/volij

      ocetra(i,j,1:kmle(i,j),idet13)     = ocetra(i,j,1:kmle(i,j),idet13)     +                   &
            ocetra(i,j,1:kmle(i,j),idet13)/(ocetra(i,j,1:kmle(i,j),idet)+safediv)*                &
            rivin(i,j,irdet)*fdt/volij
      ocetra(i,j,1:kmle(i,j),idet14)     = ocetra(i,j,1:kmle(i,j),idet14)     +                   &
            ocetra(i,j,1:kmle(i,j),idet14)/(ocetra(i,j,1:kmle(i,j),idet)+safediv)*                &
            rivin(i,j,irdet)*fdt/volij
#endif
      ! DIC is updated using the assumtions that a_t=a_c+a_n and DIC=a_c (a_t: total 
      ! alkalinity, a_c: carbonate alkalinity, a_n: contribution of nutrients to a_t). 
      ocetra(i,j,1:kmle(i,j),iano3)      = ocetra(i,j,1:kmle(i,j),iano3)      + rivin(i,j,irdin)*fdt/volij
      ocetra(i,j,1:kmle(i,j),iphosph)    = ocetra(i,j,1:kmle(i,j),iphosph)    + rivin(i,j,irdip)*fdt/volij
      ocetra(i,j,1:kmle(i,j),isilica)    = ocetra(i,j,1:kmle(i,j),isilica)    + rivin(i,j,irsi) *fdt/volij
      ocetra(i,j,1:kmle(i,j),isco212)    = ocetra(i,j,1:kmle(i,j),isco212)    + rivin(i,j,iralk)*fdt/volij  &
                                                                    + rivin(i,j,irdin)*fdt/volij  &
                                                                    + rivin(i,j,irdip)*fdt/volij
      ocetra(i,j,1:kmle(i,j),ialkali)    = ocetra(i,j,1:kmle(i,j),ialkali)    + rivin(i,j,iralk)*fdt/volij
#ifdef natDIC
      ocetra(i,j,1:kmle(i,j),inatsco212) = ocetra(i,j,1:kmle(i,j),inatsco212) + rivin(i,j,iralk)*fdt/volij  &
                                                                    + rivin(i,j,irdin)*fdt/volij  &
                                                                    + rivin(i,j,irdip)*fdt/volij
      ocetra(i,j,1:kmle(i,j),inatalkali) = ocetra(i,j,1:kmle(i,j),inatalkali) + rivin(i,j,iralk)*fdt/volij
#endif
      ocetra(i,j,1:kmle(i,j),iiron)      = ocetra(i,j,1:kmle(i,j),iiron)      + rivin(i,j,iriron)*fdt/volij*dFe_frac
      ocetra(i,j,1:kmle(i,j),idoc)       = ocetra(i,j,1:kmle(i,j),idoc)       + rivin(i,j,irdoc)*fdt/volij
      ocetra(i,j,1:kmle(i,j),idet)       = ocetra(i,j,1:kmle(i,j),idet)       + rivin(i,j,irdet)*fdt/volij

      rivinflx(i,j,irdin)    = rivin(i,j,irdin)*fdt
      rivinflx(i,j,irdip)    = rivin(i,j,irdip)*fdt
      rivinflx(i,j,irsi)     = rivin(i,j,irsi)*fdt
      rivinflx(i,j,iralk)    = rivin(i,j,iralk)*fdt
      rivinflx(i,j,iriron)   = rivin(i,j,iriron)*fdt*dFe_frac
      rivinflx(i,j,irdoc)    = rivin(i,j,irdoc)*fdt
      rivinflx(i,j,irdet)    = rivin(i,j,irdet)*fdt 
    ENDIF
  ENDDO
  ENDDO
!$OMP END PARALLEL DO

!--------------------------------------------------------------------------------
end subroutine apply_rivin


!********************************************************************************
end module mo_apply_rivin

