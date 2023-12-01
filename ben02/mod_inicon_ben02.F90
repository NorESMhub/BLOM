! ------------------------------------------------------------------------------
! Copyright (C) 2015-2020 Mats Bentsen

! This file is part of BLOM.

! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.

! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.

! You should have received a copy of the GNU Lesser General Public License
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

module mod_inicon_ben02

  use mod_constants, only: t0deg
  use mod_xc
  use mod_state,     only: temp, saln
  use mod_forcing,   only: sstclm, ricclm
  use mod_ben02,     only: cd_d, ch_d, ce_d, wg2_d, cd_m, ch_m, ce_m, &
                           wg2_m, rhoa, tsi_tda, tml_tda, sml_tda, &
                           alb_tda, fice_tda, ntda, rnfres, &
                           cd_r,ch_r,ce_r,wg2_r,rhoa_r,albdif
  use mod_thdysi,    only: tsrfm, ticem, albs_f, fice_max
  use mod_seaice,    only: ficem, hicem, hsnwm, iagem

  implicit none
  private

  public inicon_ben02

contains

  subroutine inicon_ben02

    ! --- ------------------------------------------------------------------
    ! --- Define initial conditions specifically when using the Bentsen and
    ! --- Drange (2002) forcing and thermodynamic sea ice
    ! --- ------------------------------------------------------------------

    integer :: i,j,l

    ! --- ------------------------------------------------------------------
    ! --- Set initial conditions for:
    ! --- * thermodynamic sea ice model
    ! --- * accumulation variables
    ! --- * transfer coefficients, gustiness squared, and air density
    ! --- * runoff reservoar
    ! --- ------------------------------------------------------------------

    !$OMP PARALLEL DO PRIVATE(l,i)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))

          if (ricclm(i,j,12) < .05) then
            hicem(i,j) = 0.
            ficem(i,j) = 0.
            hsnwm(i,j) = 0.
            tsrfm(i,j) = temp(i,j,1)+t0deg
            tsi_tda(i,j) = temp(i,j,1)+t0deg
          else
            hicem(i,j) = 2.
            ficem(i,j) = min(fice_max,ricclm(i,j,12)) ! assuming start in January
            hsnwm(i,j) = 0.1
            tsrfm(i,j) = sstclm(i,j,1)
            tsi_tda(i,j) = (sstclm(i,j,1) &
                 -(1.-ficem(i,j))*(temp(i,j,1)+t0deg)) &
                 /ficem(i,j)
            iagem(i,j) = 0.
          end if
          ticem(i,j) = tsrfm(i,j)
          tml_tda(i,j) = temp(i,j,1)+t0deg
          sml_tda(i,j) = saln(i,j,1)
          alb_tda(i,j) = albs_f*ficem(i,j)+albdif*(1.-ficem(i,j))
          fice_tda(i,j) = ficem(i,j)

          cd_d(i,j) = cd_r
          ch_d(i,j) = ch_r
          ce_d(i,j) = ce_r
          wg2_d(i,j) = wg2_r*1.e-4
          cd_m(i,j) = cd_r
          ch_m(i,j) = ch_r
          ce_m(i,j) = ce_r
          wg2_m(i,j) = wg2_r*1.e-4
          rhoa(i,j) = rhoa_r*1.e3

          rnfres(i,j) = 0.

        end do
      end do
    end do
    !$OMP END PARALLEL DO

    ntda = 1

  end subroutine inicon_ben02

end module mod_inicon_ben02
