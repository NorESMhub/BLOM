! Copyright (C) 2020  J. Schwinger
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


      SUBROUTINE calc_idepth(kpie,kpje,kpke,pddpo,ptiestw)
!**********************************************************************
!
!**** *CALC_IDEPTH * - .
!
!     J.Schwinger            *GFI, UiB*       2013-05-06
!     - replaces calc_ez with enhanced functionality
!
!     Purpose
!     -------
!     find lowest mass containing layer in the euphotic zone
!     find k-index of 100,500,1000,2000, and 4000 m-surfaces
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!     *INTEGER*   *kpie*    - 1st dimension of model grid.
!     *INTEGER*   *kpje*    - 2nd dimension of model grid.
!     *INTEGER*   *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*      *pddpo*   - size of grid cell (3rd dimension) [m].
!     *REAL*      *ptiestw* - depth of grid cell upper interface [m].
!
!**********************************************************************

      USE mo_biomod
      USE mo_param1_bgc, only: dp_min, dp_ez

      implicit none

      REAL    :: pddpo(kpie,kpje,kpke),ptiestw(kpie,kpje,kpke+1)
      INTEGER :: kpie,kpje,kpke,i,j,k


!$OMP PARALLEL DO
      DO j=1,kpje
      DO i=1,kpie
        
        kwrbioz(i,j)=1
        DO k=2,kpke
          IF(pddpo(i,j,k) .gt. dp_min .and. ptiestw(i,j,k) .lt. dp_ez ) THEN
            kwrbioz(i,j)=k
          ENDIF
        END DO

      END DO
      END DO
!$OMP END PARALLEL DO

      k0100(:,:)=0
      k0500(:,:)=0
      k1000(:,:)=0
      k2000(:,:)=0
      k4000(:,:)=0

!$OMP PARALLEL DO
      DO j=1,kpje
      DO i=1,kpie
        
        DO k=2,kpke
          IF(pddpo(i,j,k) .gt. dp_min .and. ptiestw(i,j,k+1) .gt. 100.0 ) THEN
            k0100(i,j)=k
            exit
          ENDIF
        END DO

        DO k=2,kpke
          IF(pddpo(i,j,k) .gt. dp_min .and. ptiestw(i,j,k+1) .gt. 500.0 ) THEN
            k0500(i,j)=k
            exit
          ENDIF
        END DO

        DO k=2,kpke
          IF(pddpo(i,j,k) .gt. dp_min .and. ptiestw(i,j,k+1) .gt. 1000.0 ) THEN
            k1000(i,j)=k
            exit
          ENDIF
        END DO

        DO k=2,kpke
          IF(pddpo(i,j,k) .gt. dp_min .and. ptiestw(i,j,k+1) .gt. 2000.0 ) THEN
            k2000(i,j)=k
            exit
          ENDIF
        END DO

        DO k=2,kpke
          IF(pddpo(i,j,k) .gt. dp_min .and. ptiestw(i,j,k+1) .gt. 4000.0 ) THEN
            k4000(i,j)=k
            exit
          ENDIF
        END DO

      END DO
      END DO
!$OMP END PARALLEL DO

     RETURN
     END
