! Copyright (C) 2020  J. Tjiputra, J. Schwinger
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


     SUBROUTINE PREFTRC(kpie,kpje,omask)
!****************************************************************
!
!**** *PREFTRC* - update preformed tracers in the mixed layer.
!
!     J. Tjiputra, J.Schwinger,    *BCCR, Bergen*   2015-01-23
!
!     Modified
!     --------
!     J.Tjiputra,       *Uni Research, Bergen*   2018-04-12
!     - added preformed DIC tracer
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
      use mo_vgrid, only: kmle

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
