! ------------------------------------------------------------------------------
! Copyright (C) 2021-2025 Mats Bentsen, Mariana Vertenstein, Mehmet Ilicak
!
! This file is part of BLOM.
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
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

module mod_ale_forcing
! ------------------------------------------------------------------------------
! This module contains procedures for computing penetration factors and buoyancy
! flux at model layer interfaces when the ALE method is used.
! ------------------------------------------------------------------------------

  use mod_types,     only: r8
  use mod_constants, only: grav, spcifh, alpha0, onem, onecm, onemu
  use mod_xc
  use mod_eos,       only: dsigdt0, dsigds0
  use mod_state,     only: dp, temp, saln, p
  use mod_swabs,     only: swamxd, swfc1, swfc2, swal1, swal2
  use mod_forcing,   only: surflx, sswflx, salflx, brnflx, buoyfl, &
                           t_sw_nonloc, s_br_nonloc
  use mod_cmnfld,    only: mlts
  use mod_checksum,  only: csdiag, chksum

  implicit none
  private

  public :: ale_forcing

contains

  subroutine ale_forcing(m, n, mm, nn, k1m, k1n)
  ! ----------------------------------------------------------------------------
  ! Compute penetration factors for shortwave and brine flux and compute
  ! interface buoyancy flux.
  ! ----------------------------------------------------------------------------

    ! Arguments
    integer, intent(in) :: m, n, mm, nn, k1m, k1n

    ! Local variables
    ! Numeric constants for brine absorption profile.
    real(r8), parameter :: cbra1 = 2._r8**(1._r8/3._r8)
    real(r8), parameter :: cbra2 = cbra1*cbra1/12._r8
    ! real(r8), parameter :: cbra1 = 2._r8**(1._r8/3._r8), &
    ! real(r8), parameter :: cbra2 = cbra1*cbra1/288._r8
    real(r8) :: cpi, gaa, pmax, lei1, lei2, lei, q, q3, pmaxi, nlbot, &
                dsgdt, dsgds
    real(r8) :: hf, hfsw, hfns, sf, sfbr, sfnb
    integer  :: i, j, k, l, kmax, kn

    ! Set some constants:
    cpi = 1._r8/spcifh      ! Multiplicative inverse of specific heat capacity.
    gaa = grav*alpha0*alpha0

    ! --------------------------------------------------------------------------
    ! Compute shortwave flux penetration factors.
    ! --------------------------------------------------------------------------

    ! Maximum pressure of shortwave absorption.
    pmax = swamxd*onem

    !$omp parallel do private(l, i, lei1, lei2, kmax, k, kn, pmaxi, nlbot)
    do j = 1, jj
      do l = 1, isp(j)
        do i = max(1, ifp(j,l)), min(ii, ilp(j,l))

          ! Penetration factors at layer interfaces.
          lei1 = 1._r8/(swal1(i,j)*onem)
          lei2 = 1._r8/(swal2(i,j)*onem)
          kmax = 1
          t_sw_nonloc(i,j,1) = 1._r8
          do k = 1, kk
            kn = k + nn
            if (dp(i,j,kn) > onemu) then
              t_sw_nonloc(i,j,k+1) = &
                   swfc1(i,j)*exp( - lei1*min(pmax, p(i,j,k+1))) &
                 + swfc2(i,j)*exp( - lei2*min(pmax, p(i,j,k+1)))
              kmax = k
            else
              t_sw_nonloc(i,j,k+1) = t_sw_nonloc(i,j,k)
            endif
            if (p(i,j,k+1) > pmax) exit
          enddo

          ! Modify penetration factors so that fluxes destined to penetrate
          ! below the lowest model layer are evenly absorbed in the water
          ! column.
          pmaxi = 1._r8/min(pmax, p(i,j,kmax+1))
          nlbot = t_sw_nonloc(i,j,kmax+1)
          do k = kmax+1, kk+1
            t_sw_nonloc(i,j,k) = 0._r8
          enddo
          do k = kmax, 2, -1
            kn = k + nn
            if (dp(i,j,kn) > onemu) then
              t_sw_nonloc(i,j,k) = t_sw_nonloc(i,j,k) - nlbot*p(i,j,k)*pmaxi
            else
              t_sw_nonloc(i,j,k) = t_sw_nonloc(i,j,k+1)
            endif
          enddo

        enddo
      enddo
    enddo
    !$omp end parallel do

    ! --------------------------------------------------------------------------
    ! Compute brine flux penetration factors.
    ! --------------------------------------------------------------------------

    !$omp parallel do private(l, i, lei, pmax, kmax, k, kn, q, q3, pmaxi, nlbot)
    do j = 1, jj
      do l = 1, isp(j)
        do i = max(1, ifp(j,l)), min(ii, ilp(j,l))

          ! Penetration factors at layer interfaces.
          lei = 1._r8/(mlts(i,j)*onem)
          pmax = cbra1*mlts(i,j)*onem
          kmax = 1
          s_br_nonloc(i,j,1) = 1._r8
          do k = 1, kk
            kn = k + nn
            if (dp(i,j,kn) > onemu) then
              q = min(cbra1, lei*p(i,j,k+1))
              q3 = q*q*q
              ! s_br_nonloc(i,j,k+1) = &
              !    1._r8 - cbra2*q*q3*q3*(q3*(35._r8*q3 - 182._r8) + 260._r8)
              s_br_nonloc(i,j,k+1) = 1._r8 - cbra2*q*q3*(7._r8-2._r8*q3)
              kmax = k
            else
              s_br_nonloc(i,j,k+1) = s_br_nonloc(i,j,k)
            endif
            if (p(i,j,k+1) > pmax) exit
          enddo

          ! Modify penetration factors so that fluxes destined to penetrate
          ! below the lowest model layer are evenly absorbed in the water
          ! column.
          pmaxi = 1._r8/min(pmax, p(i,j,kmax+1))
          nlbot = s_br_nonloc(i,j,kmax+1)
          do k = kmax+1, kk+1
            s_br_nonloc(i,j,k) = 0._r8
          enddo
          do k = kmax, 2, -1
            kn = k + nn
            if (dp(i,j,kn) > onemu) then
              s_br_nonloc(i,j,k) = s_br_nonloc(i,j,k) - nlbot*p(i,j,k)*pmaxi
            else
              s_br_nonloc(i,j,k) = s_br_nonloc(i,j,k+1)
            endif
          enddo

        enddo
      enddo
    enddo
    !$omp end parallel do

    ! --------------------------------------------------------------------------
    ! Compute buoyancy flux.
    ! --------------------------------------------------------------------------

    !$omp parallel do private(l, i, dsgdt, dsgds, hf, hfsw, hfns, sf, sfbr, sfnb, k)
    do j = 1, jj
      do l = 1, isp(j)
        do i = max(1, ifp(j,l)), min(ii, ilp(j,l))

          ! Derivatives of potential density referenced at the surface.
          dsgdt = dsigdt0(temp(i,j,k1n), saln(i,j,k1n))
          dsgds = dsigds0(temp(i,j,k1n), saln(i,j,k1n))

          ! Surface heat fluxes.
          hf   = surflx(i,j) ! Total.
          hfsw = sswflx(i,j) ! Shortwave.
          hfns = hf - hfsw   ! Non-shortwave.

          ! Surface salt fluxes.
          sf   = salflx(i,j) ! Total.
          sfbr = brnflx(i,j) ! Brine.
          sfnb = sf - sfbr   ! Non-brine.

          ! Surface buoyancy flux [m2 s-3].
          buoyfl(i,j,1) = - (dsgdt*hf*cpi + dsgds*sf)*gaa

          ! Buoyancy flux at subsurface layer interfaces [m2 s-3].
          do k = 2, kk+1
            buoyfl(i,j,k) = - ( dsgdt*t_sw_nonloc(i,j,k)*hfsw*cpi &
                              + dsgds*s_br_nonloc(i,j,k)*sfbr)*gaa
          enddo

        enddo
      enddo
    enddo
    !$omp end parallel do

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'ale_forcing:'
      endif
      call chksum(t_sw_nonloc, kk+1, halo_ps, 't_sw_nonloc')
      call chksum(s_br_nonloc, kk+1, halo_ps, 's_br_nonloc')
      call chksum(buoyfl     , kk+1, halo_ps, 'buoyfl'     )
    endif

  end subroutine ale_forcing

end module mod_ale_forcing
