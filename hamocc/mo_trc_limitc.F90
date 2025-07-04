! Copyright (C) 2020  J. Schwinger, M. Bentsen
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

module mo_trc_limitc

  implicit none
  private

  public :: trc_limitc

contains

  subroutine trc_limitc(nn)

    !***********************************************************************************************
    ! Remove negative tracer values.
    !
    ! Remove negative tracer values in the first layer in a mass
    ! conservative fashion (i.e. the mass deficit removed is
    ! transfered to non-negative points by a multiplicative
    ! correction). This is done since the virtual tracer fluxes
    ! (applied in mxlayr.F directly before HAMOCC is called) can
    ! cause negative tracer values in regions with low concentration
    ! and strong precipitation.
    !
    ! J. Schwinger      *GFI, UiB        initial version, 2014-06-17
    !
    ! Modified
    ! J.Schwinger,      *Uni Research, Bergen*   2018-04-12
    !  - fixed a bug related to the 2 time-level scheme
    !***********************************************************************************************

    use mod_xc,      only: ii,jj,ips,ifp,isp,ilp,xcsum
    use mod_grid,    only: scp2
    use mod_state,   only: dp
    use mod_tracers, only: ntrbgc, itrbgc, trc
    use mod_utility, only: util1
    use mo_kind,     only: rp

    ! Arguments
    integer :: nn

    ! Local variables
    integer :: i,j,l,nt,kn
    real(rp):: trbudo(ntrbgc),trbudn,q

    ! --- ------------------------------------------------------------------
    ! --- - compute tracer budgets before removing negative values
    ! --- ------------------------------------------------------------------

    kn=1+nn

    do nt=1,ntrbgc

      util1(:,:)=0._rp

      !$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j) = util1(i,j)+trc(i,j,kn,itrbgc+nt-1)*dp(i,j,kn)*scp2(i,j)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

      call xcsum(trbudo(nt),util1,ips)

    enddo

    ! --- ------------------------------------------------------------------
    ! --- - remove negative tracer values in the surface layer
    ! --- ------------------------------------------------------------------

    !$OMP PARALLEL DO PRIVATE(j,l,i)
    do nt=itrbgc,itrbgc+ntrbgc-1
      do j=1,jj
        do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            trc(i,j,kn,nt) = max(trc(i,j,kn,nt),0.0_rp)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    ! --- ------------------------------------------------------------------
    ! --- - recalculate and correct tracer budgets
    ! --- ------------------------------------------------------------------

    do nt=1,ntrbgc

      util1(:,:)=0._rp

      !$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            util1(i,j) = util1(i,j)+trc(i,j,kn,itrbgc+nt-1)*dp(i,j,kn)*scp2(i,j)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

      call xcsum(trbudn,util1,ips)
      q = trbudo(nt)/max(1.e-14_rp,trbudn)

      !$OMP PARALLEL DO PRIVATE(l,i)
      do j=1,jj
        do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            trc(i,j,kn,itrbgc+nt-1) = trc(i,j,kn,itrbgc+nt-1)*q
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

    enddo

  end subroutine trc_limitc

end module mo_trc_limitc
