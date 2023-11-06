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

MODULE MO_PREFTRC

  implicit none
  private

  public :: PREFTRC

CONTAINS

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

    use mo_carbch,     only: ocetra
    use mo_param1_bgc, only: ialkali,ioxygen,iphosph,iprefalk,iprefdic,iprefo2,iprefpo4,isco212
    use mo_vgrid,      only: kmle

    ! Arguments
    integer :: kpie,kpje
    real    :: omask(kpie,kpje)

    ! Local variables
    integer :: i,j

    do j=1,kpje
      do i=1,kpie
        if (omask(i,j) .gt. 0.5 ) then
          ocetra(i,j,1:kmle(i,j),iprefo2)  = ocetra(i,j,1:kmle(i,j),ioxygen)
          ocetra(i,j,1:kmle(i,j),iprefpo4) = ocetra(i,j,1:kmle(i,j),iphosph)
          ocetra(i,j,1:kmle(i,j),iprefalk) = ocetra(i,j,1:kmle(i,j),ialkali)
          ocetra(i,j,1:kmle(i,j),iprefdic) = ocetra(i,j,1:kmle(i,j),isco212)
        endif
      enddo
    enddo


  END SUBROUTINE PREFTRC

END MODULE MO_PREFTRC
