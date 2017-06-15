       SUBROUTINE CYANO(kpie,kpje,kpke,ptho,pddpo,omask)
!**********************************************************************
!
!**** *CYANO* -  .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,             *MPI-MaD, HH*    10.04.01
!     - included : surface reduction of gaseous nitrogen
!
!     I. Kriest, GEOMAR, 11.08.2016
!     - included T-dependence of cyanobacteria growth
!     - modified oxygen stoichiometry for N2-Fixation
!
!     Purpose
!     -------
!     Nitrate reduction by cyano bacteria (2NO3 + O2 => N2O + O2).
!
!     Method:
!     ------
!
!     *CALL*       *CYANO(kpie,kpje,kpke,pddpo)*
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *ptho*    - potential temperature.
!
!     Externals
!     ---------
!     .
!**********************************************************************

      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      USE mo_control_bgc
      use mo_param1_bgc 
      USE mo_bgcmean, only: jintnfix,accsrf


      implicit none

      INTEGER :: kpie,kpje,kpke
      REAL :: ptho(kpie,kpje,kpke)
      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: omask(kpie,kpje)

      INTEGER :: i,j,k
      REAL :: oldocetra
      REAL :: ttemp,nfixtfac
      REAL :: aux2d_nfix(kpie,kpje)

      aux2d_nfix(:,:)=0.0
      
!
!  N-fixation by cyano bacteria (this is not a surface flux!)
!
      DO k=1,kmle
!$OMP PARALLEL DO PRIVATE(oldocetra) 
      DO j=1,kpje
      DO i=1,kpie
        IF(omask(i,j).gt.0.5) THEN
          IF(ocetra(i,j,k,iano3).LT.(rnit*ocetra(i,j,k,iphosph))) THEN

            oldocetra = ocetra(i,j,k,iano3)
            ttemp = ptho(i,j,k)

! Temperature dependence of nitrogen fixation, Kriest and Oschlies 2015.
            nfixtfac = MAX(0.0,tf2*ttemp*ttemp + tf1*ttemp + tf0)/tff

            ocetra(i,j,k,iano3)=ocetra(i,j,k,iano3)*(1-bluefix*nfixtfac) &
     &                      +bluefix*nfixtfac*rnit*ocetra(i,j,k,iphosph)

            ocetra(i,j,k,igasnit)=ocetra(i,j,k,igasnit)-                 &
     &         (ocetra(i,j,k,iano3)-oldocetra)*(1./2.)

! Note: to fix one mole N2 requires: N2+H2O+y*O2 = 2* HNO3 <-> y=2.5 mole O2.
! I.e., to release one mole HNO3 = H+ + NO3- requires 1.25 mole O2
            ocetra(i,j,k,ioxygen)=ocetra(i,j,k,ioxygen)-                 &
     &         (ocetra(i,j,k,iano3)-oldocetra)*1.25


            aux2d_nfix(i,j) = aux2d_nfix(i,j) +                          &
     &         (ocetra(i,j,k,iano3)-oldocetra)*pddpo(i,j,k)

            ENDIF  
         ENDIF  
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      ENDDO


! Accumulate 2d diagnostics 
     call accsrf(jintnfix,aux2d_nfix,omask,0)

      RETURN
      END
