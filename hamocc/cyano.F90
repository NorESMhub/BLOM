       SUBROUTINE CYANO(kpie,kpje,kpke,pddpo,omask)
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
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
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


      implicit none

      REAL :: pddpo(kpie,kpje,kpke),omask(kpie,kpje)
      INTEGER :: kpie,kpje,kpke, i,j,k
      REAL :: oldocetra,contppm,osum,nsum
      
      contppm=1/0.35e-3
!
!  N-fixation by cyano bacteria (this is not a surface flux!)
!
!$OMP PARALLEL DO PRIVATE(oldocetra) 
      DO j=1,kpje
      DO i=1,kpie
        IF(pddpo(i,j,1).GT.1.e-6.and.omask(i,j).gt.0.5) THEN
          IF(ocetra(i,j,1,iano3).LT.(rnit*ocetra(i,j,1,iphosph))) THEN

!          nsum=ocetra(i,j,1,iano3)+2.*ocetra(i,j,1,igasnit)
!          osum=ocetra(i,j,1,ioxygen)+1.5*ocetra(i,j,1,iano3)

            oldocetra = ocetra(i,j,1,iano3)
            ocetra(i,j,1,iano3)= ocetra(i,j,1,iano3)*(1-bluefix)       &
     &                          +bluefix*rnit*ocetra(i,j,1,iphosph)

            ocetra(i,j,1,igasnit)=ocetra(i,j,1,igasnit)-               &
     &      (ocetra(i,j,1,iano3)-oldocetra)*(1./2.)

            ocetra(i,j,1,ioxygen)=ocetra(i,j,1,ioxygen)-               &
     &      (ocetra(i,j,1,iano3)-oldocetra)*1.5

!       IF(abs(ocetra(i,j,1,iano3)+2.*ocetra(i,j,1,igasnit)-nsum).gt.1.e-205)THEN
!          write(io_stdo_bgc,*) 'CYANO',i,j
!          write(io_stdo_bgc,*) 'N',ocetra(i,j,1,iano3)+2.*ocetra(i,j,1,igasnit), &
!     &    nsum
!       ENDIF
!       IF(abs(ocetra(i,j,1,ioxygen)+1.5*ocetra(i,j,1,iano3)-osum).gt.1.e-12)THEN
!          write(io_stdo_bgc,*) 'CYANO',i,j
!          write(io_stdo_bgc,*) 'O',ocetra(i,j,1,ioxygen)+1.5*ocetra(i,j,1,iano3),&
!     &    osum 
!       ENDIF

!#ifdef DIFFAT
!            atm(i,j,iatmn2)=atm(i,j,iatmn2) -                          &
!     &      (ocetra(i,j,1,iano3)-oldocetra)/2.*contppm*pddpo(i,j,1)
!#endif
            ENDIF  
         ENDIF  
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END
