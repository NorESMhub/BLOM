! Copyright (C) 2001  Ernst Maier-Reimer, S. Legutke
! Copyright (C) 2020  J. Schwinger, I. Kriest
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


SUBROUTINE CYANO(kpie,kpje,kpke,kbnd,pddpo,omask,ptho)
!**********************************************************************
!
!**** *CYANO* -  .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!     - included : surface reduction of gaseous nitrogen
!
!     I.Kriest,         *GEOMAR, Kiel*           2016-08-11
!     - included T-dependence of cyanobacteria growth
!     - modified oxygen stoichiometry for N2-Fixation
!
!     J.Schwinger,      *Uni Research, Bergen*   2018-04-12
!     - moved accumulation of all output fields to seperate subroutine,
!       related code-restructuring
!     - added reduction of alkalinity through N-fixation
!
!     Purpose
!     -------
!     Nitrogen-fixation by cyano bacteria, followed by remineralisation 
!     and nitrification
!
!     Method:
!     ------
!
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *INTEGER* *kbnd*    - nb of halo grid points
!     *REAL*    *ptho*    - potential temperature.
!
!     Externals
!     ---------
!     .
!**********************************************************************

  use mo_carbch,     only: ocetra
  use mo_biomod,     only: bluefix,intnfix,rnit,tf0,tf1,tf2,tff
  use mo_param1_bgc, only: ialkali,iano3,igasnit,iphosph,ioxygen
  use mo_vgrid,      only: kmle,kwrbioz
  use mo_control_bgc,only: leuphotic_cya
#ifdef natDIC
  use mo_param1_bgc, only: inatalkali
#endif
#ifdef extNcycle
  use mo_param1_bgc, only: ianh4
#endif

  implicit none

  INTEGER, intent(in) :: kpie,kpje,kpke,kbnd
  REAL,    intent(in) :: pddpo(kpie,kpje,kpke)
  REAL,    intent(in) :: omask(kpie,kpje)
  REAL,    intent(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)

  ! Local variables
  INTEGER :: i,j,k
  REAL :: oldocetra,anavail,dansp,dox,dalk
  REAL :: ttemp,nfixtfac

  intnfix(:,:)=0.0

!
! N-fixation by cyano bacteria (followed by remineralisation and nitrification,
! or, for the extended nitrogen cycle only by remin to NH4), 
! it is assumed here that this process is limited to the mixed layer
!
  DO j=1,kpje
  DO i=1,kpie
  IF(omask(i,j).gt.0.5) THEN
     DO k=1,merge(kwrbioz(i,j),kmle(i,j),leuphotic_cya) ! if leuphotic_cya=.true., do bluefix only in euphotic zone
#ifdef extNcycle
        ! assuming nitrate and ammonium required for cyanobacteria growth (as bulk PP)
        anavail = ocetra(i,j,k,iano3)+ocetra(i,j,k,ianh4)
#else
        anavail = ocetra(i,j,k,iano3)
#endif 
        IF(anavail.LT.(rnit*ocetra(i,j,k,iphosph))) THEN

            ttemp = min(40.,max(-3.,ptho(i,j,k)))

            ! Temperature dependence of nitrogen fixation, Kriest and Oschlies 2015.
            nfixtfac = MAX(0.0,tf2*ttemp*ttemp + tf1*ttemp + tf0)/tff

#ifndef extNcycle
            oldocetra = ocetra(i,j,k,iano3)
            ocetra(i,j,k,iano3)=ocetra(i,j,k,iano3)*(1.-bluefix*nfixtfac) &
     &                      +bluefix*nfixtfac*rnit*ocetra(i,j,k,iphosph)
            dansp=ocetra(i,j,k,iano3)-oldocetra
            ! Note: to fix one mole N2 requires: N2+H2O+y*O2 = 2* HNO3 <-> y=2.5 mole O2.
            ! I.e., to release one mole HNO3 = H+ + NO3- requires 1.25 mole O2
            dox  = -dansp*1.25
            ! Nitrogen fixation followed by remineralisation and nitrification decreases
            ! alkalinity by 1 mole per mole nitrogen fixed (Wolf-Gladrow et al. 2007)
            dalk = -dansp 
#else
            oldocetra = ocetra(i,j,k,ianh4)
            ocetra(i,j,k,ianh4)=ocetra(i,j,k,ianh4)*(1.-bluefix*nfixtfac) &
     &                      +bluefix*nfixtfac*rnit*ocetra(i,j,k,iphosph)
            dansp=ocetra(i,j,k,ianh4)-oldocetra
            dox  = dansp*0.75
            dalk = dansp 

#endif
            ocetra(i,j,k,igasnit)=ocetra(i,j,k,igasnit)-dansp*(1./2.)

            ocetra(i,j,k,ioxygen)=ocetra(i,j,k,ioxygen)+dox

            ocetra(i,j,k,ialkali)=ocetra(i,j,k,ialkali)+dalk
#ifdef natDIC
            ocetra(i,j,k,inatalkali)=ocetra(i,j,k,inatalkali)+dalk
#endif

            intnfix(i,j) = intnfix(i,j) + dansp*pddpo(i,j,k)

        ENDIF
     ENDDO
  ENDIF
  ENDDO
  ENDDO

END SUBROUTINE CYANO
