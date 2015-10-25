      SUBROUTINE CALC_BOT(kpie,kpje,kpke,pddpo)
!
!**********************************************************************
!
!**** *CALC_SED* - .
!
!     Karen Assmann          *BCCR*           14.11.05
! 
!     set lowest mass containing layer and its thickness
!     since both kbo and bolay need to be recalculated every
!     time step in MICOM
!
!     *CALL*       *BODENSED
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!
! ******************************************************************

      USE mo_biomod
      USE mo_control_bgc
      USE mo_param1_bgc, only: dp_min_sink
      USE mod_xc

      implicit none

      INTEGER :: kpie,kpje,kpke
      REAL    :: pddpo(kpie,kpje,kpke)

      INTEGER :: i,j,k

      kbo(:,:)  =1
      bolay(:,:)=0.0

!$OMP PARALLEL DO
      DO j=1,kpje
      DO i=1,kpie
      DO k=kpke,1,-1
         IF(pddpo(i,j,k).GT.dp_min_sink) THEN
            bolay(i,j)=pddpo(i,j,k)
            kbo(i,j)=k
            exit
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END
