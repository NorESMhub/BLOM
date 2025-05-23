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

module mo_preftrc

  implicit none
  private

  public :: preftrc

contains

  subroutine preftrc(kpie,kpje,omask)

    !***********************************************************************************************
    ! Update preformed tracers in the mixed layer.
    !
    ! Preformed tracers are set to the value of their full counterparts in the mixed layer.
    !
    ! J. Tjiputra, J.Schwinger,    *BCCR, Bergen*   2015-01-23
    !
    ! Modified
    ! J.Tjiputra,       *Uni Research, Bergen*   2018-04-12
    !  - added preformed DIC tracer
    !***********************************************************************************************

    use mo_carbch,     only: ocetra
    use mo_param1_bgc, only: ialkali,ioxygen,iphosph,isilica,iprefalk,iprefdic,iprefo2,iprefpo4,isco212,iprefsilica,&
                             idoc,idocsl,idocsr,idocr,iprefdoc,iprefdocsl,iprefdocsr,iprefdocr
    use mo_vgrid,      only: kmle
    use mo_control_bgc,only: use_DOMclasses

    ! Arguments
    integer, intent(in) :: kpie ! 1st dimension of model grid.
    integer, intent(in) :: kpje ! 2nd dimension of model grid.
    real,    intent(in) :: omask(kpie,kpje) ! land-ocean mask

    ! Local variables
    integer :: i,j

    do j=1,kpje
      do i=1,kpie
        if (omask(i,j) > 0.5 ) then
          ocetra(i,j,1:kmle(i,j),iprefo2)  = ocetra(i,j,1:kmle(i,j),ioxygen)
          ocetra(i,j,1:kmle(i,j),iprefpo4) = ocetra(i,j,1:kmle(i,j),iphosph)
          ocetra(i,j,1:kmle(i,j),iprefsilica)= ocetra(i,j,1:kmle(i,j),isilica)
          ocetra(i,j,1:kmle(i,j),iprefalk) = ocetra(i,j,1:kmle(i,j),ialkali)
          ocetra(i,j,1:kmle(i,j),iprefdic) = ocetra(i,j,1:kmle(i,j),isco212)
          if (use_DOMclasses) then
            ocetra(i,j,1:kmle(i,j),iprefdoc)   = ocetra(i,j,1:kmle(i,j),idoc)
            ocetra(i,j,1:kmle(i,j),iprefdocsl) = ocetra(i,j,1:kmle(i,j),idocsl)
            ocetra(i,j,1:kmle(i,j),iprefdocsr) = ocetra(i,j,1:kmle(i,j),idocsr)
            ocetra(i,j,1:kmle(i,j),iprefdocr)  = ocetra(i,j,1:kmle(i,j),idocr)
          endif
        endif
      enddo
    enddo

  end subroutine preftrc

end module mo_preftrc
