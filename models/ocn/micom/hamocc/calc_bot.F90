      SUBROUTINE CALC_BOT(kpie,kpje,kpke,pddpo)
!
!**********************************************************************
!
!**** *CALC_SED* - .
!
!     Karen Assmann          *BCCR*           14.11.05
! 
!     set lowest mass containing layer and its thickness
!     extracted from subroutine bodensed
!     since both kbot and bolay need to be recalculated every
!     time step in MICOM
!
!     *CALL*       *BODENSED
!
!     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!     *COMMON*     *PARAM1_BGC.h* - declaration of ocean/sediment tracer.
!     *COMMON*     *COMMO1_BGC.h* - ocean/sediment tracer arrays.
!     *COMMON*     *UNITS_BGC.h*  - std I/O logical units.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!

      USE mo_biomod
      USE mo_control_bgc
      USE mod_xc

      implicit none

      REAL :: pddpo(kpie,kpje,kpke)
      INTEGER :: kpie,kpje,kpke,i,j,k,fbl

! ******************************************************************

      k=kpke
!$OMP PARALLEL DO
      DO 1321 j=1,kpje
      DO 1321 i=1,kpie
         kbo(i,j)=1
         bolay(i,j)=0.
         IF(pddpo(i,j,k).GT.0.5) THEN
            bolay(i,j)=pddpo(i,j,k)
            kbo(i,j)=k
         ENDIF
1321  CONTINUE
!$OMP END PARALLEL DO

! evaluate min depth of last layer 
      bolaymin=8000.
      DO 1322 j=1,kpje
      DO 1322 i=1,kpie
      fbl=0.
      DO 1322 k=kpke-1,1,-1
         IF(fbl.lt.0.5) THEN
         IF(pddpo(i,j,k).GT.0.5.AND.pddpo(i,j,k+1).LE.0.5) THEN
            bolay(i,j)=pddpo(i,j,k)
            kbo(i,j)=k
            bolaymin = min(bolaymin,bolay(i,j))
            fbl=1.
         ENDIF
         ENDIF
1322  CONTINUE
      CALL xcminr(bolaymin)
!      WRITE(io_stdo_bgc,*)  'bolaymin=', bolaymin

      RETURN
      END
