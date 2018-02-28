      SUBROUTINE PREFTRC(kpie,kpje,omask)

!****************************************************************
!
!**** *PREFTRC* - update preformed tracers in the mixed layer.
!
!     J. Tjiputra, J.Schwinger,    *BCCR, Bergen*   2015-01-23
!
!     Modified
!     --------
!
!
!     Method
!     -------
!     Preformed tracers are set to the value of their full counterparts
!     in the mixed layer.
!
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!
!**************************************************************************

      USE mo_carbch
      use mo_param1_bgc

      implicit none

      INTEGER :: kpie,kpje
      REAL    :: omask(kpie,kpje)

      INTEGER :: i,j,k

      do k=1,kmle
!$OMP PARALLEL DO
      do j=1,kpje
      do i=1,kpie
        if (omask(i,j) .gt. 0.5 ) then
          ocetra(i,j,k,iprefo2) =ocetra(i,j,k,ioxygen)
          ocetra(i,j,k,iprefpo4)=ocetra(i,j,k,iphosph)
          ocetra(i,j,k,iprefalk)=ocetra(i,j,k,ialkali)
          ocetra(i,j,k,iprefdic)=ocetra(i,j,k,isco212)
        endif
      enddo
      enddo
!$OMP END PARALLEL DO
      enddo


      END SUBROUTINE PREFTRC
