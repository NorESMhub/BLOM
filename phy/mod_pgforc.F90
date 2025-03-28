! ------------------------------------------------------------------------------
! Copyright (C) 2005-2025 Mats Bentsen, Mehmet Ilicak, Mariana Vertenstein
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

module mod_pgforc
! ------------------------------------------------------------------------------
! This module contains variables and procedures related to time integration of
! the baroclinic momentum equation.
! ------------------------------------------------------------------------------

  use dimensions,    only: idm, jdm, kdm
  use mod_xc,        only: xctilr, nbdy, ii, jj, kk, &
                           isp, ifp, ilp, isu, ifu, ilu, isv, ifv, ilv, &
                           halo_ps, mnproc, lp, ip, iu, iv, xcstop
  use mod_types,     only: r8
  use mod_constants, only: grav, epsilp, onemm, spval
  use mod_state,     only: dp, dpu, dpv, temp, saln, p, pu, pv, phi, &
                           pb_p, pbu_p, pbv_p, sealv
  use mod_eos,       only: pref, alp, p_alpha, delphi, dalpdt, dalpds, &
                           dynh_derivatives
  use mod_checksum,  only: csdiag, chksummsk

  implicit none
  private

  ! Variables to be set in namelist:
  character(len = 80) :: pgfmth

  ! Parameters.
  real(r8), parameter :: &
       wpgf = .25_r8, & ! Weight for time averaging of pressure gradient force
                        ! [].
       p0_dynh = 0._r8  ! Reference pressure for dynamic enthalpy [kg m-1 s-2].

  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2*kdm) :: &
       pgfx, &    ! x-component of baroclinic pressure gradient force [m2 s-2].
       pgfy       ! y-component of baroclinic pressure gradient force [m2 s-2].

  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: &
       pgfx_o, &  ! 'pgfx' at old time level [m2 s-2]. &
       pgfy_o     ! 'pgfy' at old time level [m2 s-2].

  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: &
       pgfxm, &   ! x-component of barotropic pressure gradient force, not
                  ! dependent on bottom pressure [m2 s-2]
       pgfym, &   ! y-component of barotropic pressure gradient force, not
                  ! dependent on bottom pressure [m2 s-2]
       xixp, &    ! Dependeny of x-component of barotropic pressure gradient
                  ! force on bottom pressure at (i,j) [m3 g-1].
       xixm, &    ! Dependeny of x-component of barotropic pressure gradient
                  ! force on bottom pressure at (i-1,j) [m3 g-1].
       xiyp, &    ! Dependeny of y-component of barotropic pressure gradient
                  ! force on bottom pressure at (i,j) [m3 g-1].
       xiym       ! Dependeny of y-component of barotropic pressure gradient
                  ! force on bottom pressure at (i,j-1) [m3 g-1].

  real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
       pgfxm_o, & ! 'pgfxm' at old time level [m2 s-2]
       pgfym_o, & ! 'pgfym' at old time level [m2 s-2]
       xixp_o, &  ! 'xixp' at old time level [m3 g-1]
       xixm_o, &  ! 'xixm' at old time level [m3 g-1]
       xiyp_o, &  ! 'xiyp' at old time level [m3 g-1]
       xiym_o     ! 'xiym' at old time level [m3 g-1].

  ! Public variables.
  public :: pgfmth, wpgf, pgfx, pgfy, pgfx_o, pgfy_o, pgfxm, pgfym, &
            xixp, xixm, xiyp, xiym, pgfxm_o, pgfym_o, &
            xixp_o, xixm_o, xiyp_o, xiym_o

  ! Public routines.
  public :: inivar_pgforc, pgforc

contains

  ! ----------------------------------------------------------------------------
  ! Private procedures.
  ! ----------------------------------------------------------------------------

  subroutine pgforc_geopotential(m, n, mm, nn, k1m, k1n)
  ! ----------------------------------------------------------------------------
  ! Compute the pressure gradient force (PGF) as the gradient of the
  ! geopotential on pressure surfaces.
  ! ----------------------------------------------------------------------------

    ! Arguments.
    integer, intent(in) :: m, n, mm, nn, k1m, k1n

    ! Local variables.
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) :: phip
    real :: dphi, alpl, alpu, prs, dphip, dphim, alplp, alpup, alplm, alpum, &
            cp, cm, phi_p, phi_m, q
    integer, dimension(idm) :: kup, kum, kvp, kvm
    integer :: i, j, k, l, kn

    !$omp parallel do private(l,i,k,kn,dphi,alpu,alpl)
    do j = 0, jj
      do l = 1, isp(j)
        do i = max(0, ifp(j,l)), min(ii, ilp(j,l))
          phip(i,j,kk+1) = 0._r8
        end do
      end do
      do k = kk, 1, -1
        kn = k + nn
        do l = 1, isp(j)
          do i = max(0, ifp(j,l)), min(ii, ilp(j,l))
            if (dp(i,j,kn) < epsilp) then
              phi (i,j,k) = phi (i,j,k+1)
              phip(i,j,k) = phip(i,j,k+1)
            else
              call delphi(p(i,j,k), p(i,j,k+1), temp(i,j,kn), saln(i,j,kn), &
                          dphi, alpu, alpl)
              phi (i,j,k) = phi (i,j,k+1) - dphi
              phip(i,j,k) = phip(i,j,k+1) + p(i,j,k+1)*alpl - p(i,j,k)*alpu
            end if
          end do
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private( &
    !$omp l,i,kup,kum,kvp,kvm,k,kn,prs,dphip,alpup,alplp,dphim,alpum,alplm, &
    !$omp cp,cm,q,phi_p,phi_m)
    do j = 1, jj

      do l = 1, isu(j)
        do i = max(1, ifu(j,l)), min(ii, ilu(j,l))
          kup(i) = kk
          kum(i) = kk
          xixp(i,j,n) = 0._r8
          xixm(i,j,n) = 0._r8
          pgfxm(i,j,n) = 0._r8
        end do
      end do

      do l = 1, isv(j)
        do i = max(1, ifv(j,l)), min(ii, ilv(j,l))
          kvp(i) = kk
          kvm(i) = kk
          xiyp(i,j,n) = 0._r8
          xiym(i,j,n) = 0._r8
          pgfym(i,j,n) = 0._r8
        end do
      end do

      do k = kk, 1, -1
        kn = k + nn

        do l = 1, isu(j)
          do i = max(1, ifu(j,l)), min(ii, ilu(j,l))

            prs = pu(i,j,k+1) - .5_r8*dpu(i,j,kn)

            do while (p(i  ,j,kup(i)) > prs)
              kup(i) = kup(i) - 1
            end do

            do while (p(i-1,j,kum(i)) > prs)
              kum(i) = kum(i) - 1
            end do

            call delphi(prs, p(i,j,kup(i)+1), &
                        temp(i,j,kup(i)+nn), saln(i,j,kup(i)+nn), &
                        dphip, alpup, alplp)

            call delphi(prs, p(i-1,j,kum(i)+1), &
                        temp(i-1,j,kum(i)+nn), saln(i-1,j,kum(i)+nn), &
                        dphim, alpum, alplm)

            cp = .25_r8*(p(i  ,j,k+1) + p(i  ,j,k))
            cm = .25_r8*(p(i-1,j,k+1) + p(i-1,j,k))
            q = prs/(cp + cm)
            cp = q*cp
            cm = q*cm

            phi_p = phi(i  ,j,kup(i)+1) - dphip
            xixp(i,j,n) = xixp(i,j,n) &
                        + ( phip(i  ,j,kup(i)+1) &
                          + p(i  ,j,kup(i)+1)*alplp - cp*(alpup-alpum)) &
                          *dpu(i,j,kn)

            phi_m = phi(i-1,j,kum(i)+1) - dphim
            xixm(i,j,n) = xixm(i,j,n) &
                        + ( phip(i-1,j,kum(i)+1) &
                          + p(i-1,j,kum(i)+1)*alplm - cm*(alpum-alpup)) &
                          *dpu(i,j,kn)

            pgfx(i,j,kn) = - (phi_p - phi_m)
            pgfxm(i,j,n) = pgfxm(i,j,n) + pgfx(i,j,kn)*dpu(i,j,kn)

          end do
        end do

        do l = 1, isv(j)
          do i = max(1, ifv(j,l)), min(ii, ilv(j,l))

            prs = pv(i,j,k+1) - .5_r8*dpv(i,j,kn)

            do while (p(i,j  ,kvp(i)) > prs)
              kvp(i) = kvp(i)-1
            end do

            do while (p(i,j-1,kvm(i)) > prs)
              kvm(i) = kvm(i)-1
            end do

            call delphi(prs, p(i,j  ,kvp(i)+1), &
                        temp(i,j  ,kvp(i)+nn), saln(i,j  ,kvp(i)+nn), &
                        dphip, alpup, alplp)

            call delphi(prs, p(i,j-1,kvm(i)+1), &
                        temp(i,j-1,kvm(i)+nn), saln(i,j-1,kvm(i)+nn), &
                        dphim, alpum, alplm)

            cp = .25_r8*(p(i,j  ,k+1) + p(i,j  ,k))
            cm = .25_r8*(p(i,j-1,k+1) + p(i,j-1,k))
            q = prs/(cp + cm)
            cp = q*cp
            cm = q*cm

            phi_p = phi(i,j  ,kvp(i)+1) - dphip
            xiyp(i,j,n) = xiyp(i,j,n) &
                        + ( phip(i,j  ,kvp(i)+1) &
                          + p(i,j  ,kvp(i)+1)*alplp - cp*(alpup-alpum)) &
                          *dpv(i,j,kn)

            phi_m = phi(i,j-1,kvm(i)+1) - dphim
            xiym(i,j,n) = xiym(i,j,n) &
                        + ( phip(i,j-1,kvm(i)+1) &
                          + p(i,j-1,kvm(i)+1)*alplm - cm*(alpum-alpup)) &
                          *dpv(i,j,kn)

            pgfy(i,j,kn) = - (phi_p - phi_m)
            pgfym(i,j,n) = pgfym(i,j,n) + pgfy(i,j,kn)*dpv(i,j,kn)

          end do
        end do

      end do

    end do
    !$omp end parallel do

  end subroutine pgforc_geopotential

  subroutine pgforc_dynamic_enthalpy(m, n, mm, nn, k1m, k1n)
  ! ----------------------------------------------------------------------------
  ! Compute the pressure gradient force (PGF) using layer-wise gradients.
  ! ----------------------------------------------------------------------------

    ! Arguments.
    integer, intent(in) :: m, n, mm, nn, k1m, k1n

    ! Local variables.
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: &
      pot_dynh, pot_dynh_pb, dynh_a, dynh_t, alpha_r
    real :: dynh_ts_t, dynh_ts_s
    integer :: i, j, k, l, kn

    !$omp parallel do private(kn, l, i, k, dynh_ts_t, dynh_ts_s)
    do j = 0, jj

      ! For each layer, obtain the potential, pot_dynh, which is the sum of
      ! dynamic enthalpy and the geopotential. Also get the linearized response
      ! of this potential to barotropic pressure excursions, pot_dynh_pb, and
      ! the geopotential, phi, at layer interfaces.

      kn = kk + nn
      do l = 1, isp(j)
        do i = max(0, ifp(j,l)), min(ii, ilp(j,l))
          pot_dynh(i,j,kk) = phi(i,j,kk+1) &
           + p_alpha(p0_dynh, p(i,j,kk+1), temp(i,j,kn), saln(i,j,kn))
          pot_dynh_pb(i,j,kk) = &
             alp(p(i,j,kk+1), temp(i,j,kn),saln(i,j,kn))*p(i,j,kk+1)
          phi(i,j,kk) = phi(i,j,kk+1) &
           + p_alpha(p(i,j,kk), p(i,j,kk+1), temp(i,j,kn), saln(i,j,kn))
        end do
      end do
      do k = kk-1, 1, -1
        kn = k + nn
        do l = 1, isp(j)
          do i = max(0, ifp(j,l)), min(ii, ilp(j,l))
            pot_dynh(i,j,k) = pot_dynh(i,j,k+1) &
             + p_alpha(p0_dynh, p(i,j,k+1), temp(i,j,kn  ),saln(i,j,kn  )) &
             - p_alpha(p0_dynh, p(i,j,k+1), temp(i,j,kn+1),saln(i,j,kn+1))
            pot_dynh_pb(i,j,k) = pot_dynh_pb(i,j,k+1) &
             + ( alp(p(i,j,k+1), temp(i,j,kn  ),saln(i,j,kn  )) &
               - alp(p(i,j,k+1), temp(i,j,kn+1),saln(i,j,kn+1)))*p(i  ,j,k+1)
            phi(i,j,k) = phi(i,j,k+1) &
             + p_alpha(p(i,j,k), p(i,j,k+1), temp(i,j,kn), saln(i,j,kn))
          end do
        end do
      end do

      ! For each layer, find the derivatives of dynamic enthalpy with respect to
      ! potential temperature and specific volume with a reference pressure
      ! (reciprocal of potential density).

      do k = 1, kk
        kn = k + nn
        do l = 1, isp(j)
          do i = max(0, ifp(j,l)), min(ii, ilp(j,l))
            if (dp(i,j,kn) < onemm) then
              dynh_a(i,j,k) = 0._r8
              dynh_t(i,j,k) = 0._r8
            else
              call dynh_derivatives(p0_dynh, p(i,j,k), p(i,j,k+1), &
                                    temp(i,j,kn), saln(i,j,kn), &
                                    dynh_ts_t, dynh_ts_s)
              dynh_a(i,j,k) = &
                 dynh_ts_s/dalpds(pref, temp(i,j,kn), saln(i,j,kn))
              dynh_t(i,j,k) = &
                 dynh_ts_t &
               - dynh_a(i,j,k)*dalpdt(pref, temp(i,j,kn), saln(i,j,kn))
            end if
            alpha_r(i,j,k) = alp(pref, temp(i,j,kn), saln(i,j,kn))
          end do
        end do
      end do
    end do

    ! Compute the pressure gradient force as the gradient on vertical coordinate
    ! surfaces of pot_dynh with additional contributions by potential
    ! temperature and specific volume gradients. Also accumulate layer pressure
    ! thickness weighted PGF and linearized response of PGF to barotropic
    ! pressure excursions.

    !$omp parallel do private(l, i, k, kn)
    do j = 1, jj

      do l = 1, isu(j)
        do i = max(1, ifu(j,l)), min(ii, ilu(j,l))
          xixp(i,j,n) = 0._r8
          xixm(i,j,n) = 0._r8
          pgfxm(i,j,n) = 0._r8
        end do
      end do

      do l = 1, isv(j)
        do i = max(1, ifv(j,l)), min(ii, ilv(j,l))
          xiyp(i,j,n) = 0._r8
          xiym(i,j,n) = 0._r8
          pgfym(i,j,n) = 0._r8
        end do
      end do

      do k = kk, 1, -1
        kn = k + nn

        do l = 1, isu(j)
          do i = max(1, ifu(j,l)), min(ii, ilu(j,l))

            pgfx(i,j,kn) = - (pot_dynh(i,j,k) - pot_dynh(i-1,j,k))
            if (dp(i-1,j,kn) >= onemm .and. dp(i,j,kn) >= onemm) then
              pgfx(i,j,kn) = pgfx(i,j,kn) &
                           + .5_r8*( (dynh_t(i-1,j,k) + dynh_t(i,j,k)) &
                                     *(temp   (i,j,k) - temp   (i-1,j,k)) &
                                   + (dynh_a(i-1,j,k) + dynh_a(i,j,k)) &
                                     *(alpha_r(i,j,k) - alpha_r(i-1,j,k)))
            endif

            pgfxm(i,j,n) = pgfxm(i,j,n) + pgfx(i,j,kn)*dpu(i,j,kn)
            xixm(i,j,n) = xixm(i,j,n) + pot_dynh_pb(i-1,j,k)*dpu(i,j,kn)
            xixp(i,j,n) = xixp(i,j,n) + pot_dynh_pb(i  ,j,k)*dpu(i,j,kn)

          end do
        end do

        do l = 1, isv(j)
          do i = max(1, ifv(j,l)), min(ii, ilv(j,l))

            pgfy(i,j,kn) = - (pot_dynh(i,j,k) - pot_dynh(i,j-1,k))
            if (dp(i,j-1,kn) >= onemm .and. dp(i,j,kn) >= onemm) then
              pgfy(i,j,kn) = pgfy(i,j,kn) &
                           + .5_r8*( (dynh_t(i,j-1,k) + dynh_t(i,j,k)) &
                                     *(temp   (i,j,k) - temp   (i,j-1,k)) &
                                   + (dynh_a(i,j-1,k) + dynh_a(i,j,k)) &
                                     *(alpha_r(i,j,k) - alpha_r(i,j-1,k)))
            endif

            pgfym(i,j,n) = pgfym(i,j,n) + pgfy(i,j,kn)*dpv(i,j,kn)
            xiym(i,j,n) = xiym(i,j,n) + pot_dynh_pb(i,j-1,k)*dpv(i,j,kn)
            xiyp(i,j,n) = xiyp(i,j,n) + pot_dynh_pb(i,j  ,k)*dpv(i,j,kn)

          end do
        end do

      end do

    end do

  end subroutine pgforc_dynamic_enthalpy

  ! ----------------------------------------------------------------------------
  ! Public procedures.
  ! ----------------------------------------------------------------------------

  subroutine inivar_pgforc()
  ! ----------------------------------------------------------------------------
  ! Initialize arrays.
  ! ----------------------------------------------------------------------------

    ! Local variables
    integer :: i,j,k

    !$omp parallel do private(i,k)
    do j = 1-nbdy, jj+nbdy
      do k = 1, 2*kk
        do i = 1-nbdy, ii+nbdy
          pgfx(i,j,k) = spval
          pgfy(i,j,k) = spval
        end do
      end do
      do k = 1, kk
        do i = 1-nbdy, ii+nbdy
          pgfx_o(i,j,k) = spval
          pgfy_o(i,j,k) = spval
        end do
      end do
      do k = 1, 2
        do i = 1-nbdy, ii+nbdy
          pgfxm(i,j,k) = spval
          pgfym(i,j,k) = spval
          xixp(i,j,k) = spval
          xixm(i,j,k) = spval
          xiyp(i,j,k) = spval
          xiym(i,j,k) = spval
        end do
      end do
      do i = 1-nbdy, ii+nbdy
        pgfxm_o(i,j) = spval
        pgfym_o(i,j) = spval
        xixp_o(i,j) = spval
        xixm_o(i,j) = spval
        xiyp_o(i,j) = spval
        xiym_o(i,j) = spval
      end do
    end do
    !$omp end parallel do

  end subroutine inivar_pgforc

  subroutine pgforc(m, n, mm, nn, k1m, k1n)
  ! ----------------------------------------------------------------------------
  ! Compute the pressure gradient force (PGF).
  ! ----------------------------------------------------------------------------

    ! Arguments.
    integer, intent(in) :: m, n, mm, nn, k1m, k1n

    ! Local variables.
    real :: q
    integer :: i, j, k, l, kn

    ! Compute new -dpu,dpv- fields.
    !$omp parallel do private(k,kn,l,i)
    do j = -2, jj+2
      do k = 1, kk
        kn = k + nn
        do l = 1, isp(j)
          do i = max(-2, ifp(j,l)), min(ii+2, ilp(j,l))
            p(i,j,k+1) = p(i,j,k) + dp(i,j,kn)
          end do
        end do
      end do
    end do
    !$omp end parallel do
    !$omp parallel do private(k,kn,l,i,q)
    do j = -1, jj+2
      do k = 1, kk
        kn = k + nn
        do l = 1, isu(j)
          do i = max(-1, ifu(j,l)), min(ii+2, ilu(j,l))
            q = min(p(i,j,kk+1), p(i-1,j,kk+1))
            dpu(i,j,kn) = .5_r8*( (min(q, p(i-1,j,k+1)) - min(q, p(i-1,j,k))) &
                                + (min(q, p(i  ,j,k+1)) - min(q, p(i  ,j,k))))
            pu(i,j,k+1) = pu(i,j,k) + dpu(i,j,kn)
          end do
        end do
        do l = 1, isv(j)
          do i = max(-1, ifv(j,l)), min(ii+2, ilv(j,l))
            q = min(p(i,j,kk+1), p(i,j-1,kk+1))
            dpv(i,j,kn) = .5_r8*( (min(q, p(i,j-1,k+1)) - min(q, p(i,j-1,k))) &
                                + (min(q, p(i,j  ,k+1)) - min(q, p(i,j  ,k))))
            pv(i,j,k+1) = pv(i,j,k) + dpv(i,j,kn)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    ! Copy old PGF fields.
    !$omp parallel do private(l,i)
    do j = -1, jj+2
      do l = 1, isu(j)
        do i = max(0, ifu(j,l)), min(ii+1, ilu(j,l))
          xixp_o(i,j) = xixp(i,j,n)
          xixm_o(i,j) = xixm(i,j,n)
          pgfxm_o(i,j) = pgfxm(i,j,n)
        end do
      end do
      do l = 1, isv(j)
        do i = max(0, ifv(j,l)), min(ii+1, ilv(j,l))
          xiyp_o(i,j) = xiyp(i,j,n)
          xiym_o(i,j) = xiym(i,j,n)
          pgfym_o(i,j) = pgfym(i,j,n)
        end do
      end do
    end do
    !$omp end parallel do
    !$omp parallel do private(k, kn, l, i)
    do j = 1, jj
      do k = kk, 1, -1
        kn = k + nn
        do l = 1, isu(j)
          do i = max(1, ifu(j,l)), min(ii, ilu(j,l))
            pgfx_o(i,j,k) = pgfx(i,j,kn)
          end do
        end do
        do l = 1, isv(j)
          do i = max(1, ifv(j,l)), min(ii, ilv(j,l))
            pgfy_o(i,j,k) = pgfy(i,j,kn)
          end do
        end do
      end do
    end do
    !$omp end parallel do

    ! Compute new PGF fields.
    if (pgfmth == 'geopotential') then
      call pgforc_geopotential(m, n, mm, nn, k1m, k1n)
    elseif (pgfmth == 'dynamic enthalpy') then
      call pgforc_dynamic_enthalpy(m, n, mm, nn, k1m, k1n)
    else
      if (mnproc == 1) then
        write (lp,'(3a)') ' pgfmth = ',trim(pgfmth),' is unsupported!'
      end if
      call xcstop('(pgforc)')
      stop '(pgforc)'
    endif

    ! Finalise vertical average PGF fields and the dependency fields of
    ! barotropic PGF on bottom pressure.

    call xctilr(pb_p, 1, 1, 1, 1, halo_ps)

    !$omp parallel do private(l,i,q,k,kn)
    do j = 1,jj

      do l = 1, isu(j)
        do i = max(1, ifu(j,l)), min(ii, ilu(j,l))
          q = 1._r8/pbu_p(i,j)
          pgfxm(i,j,n) = pgfxm(i,j,n)*q
          xixp(i,j,n) = xixp(i,j,n)*q
          xixm(i,j,n) = xixm(i,j,n)*q
        end do
      end do
      do l = 1, isv(j)
        do i = max(1, ifv(j,l)), min(ii, ilv(j,l))
          q = 1._r8/pbv_p(i,j)
          pgfym(i,j,n) = pgfym(i,j,n)*q
          xiyp(i,j,n) = xiyp(i,j,n)*q
          xiym(i,j,n) = xiym(i,j,n)*q
        end do
      end do

      do k = 1, kk
        kn = k + nn
        do l = 1, isu(j)
          do i = max(1, ifu(j,l)), min(ii, ilu(j,l))
            pgfx(i,j,kn) = pgfx(i,j,kn) - pgfxm(i,j,n)
          end do
        end do
        do l = 1, isv(j)
          do i = max(1, ifv(j,l)), min(ii, ilv(j,l))
            pgfy(i,j,kn) = pgfy(i,j,kn) - pgfym(i,j,n)
          end do
        end do
      end do

      do l = 1, isu(j)
        do i = max(1, ifu(j,l)), min(ii, ilu(j,l))
          pgfxm(i,j,n) = pgfxm(i,j,n) + xixp(i,j,n) - xixm(i,j,n)
          xixp(i,j,n) = xixp(i,j,n)/pb_p(i  ,j)
          xixm(i,j,n) = xixm(i,j,n)/pb_p(i-1,j)
        end do
      end do
      do l = 1, isv(j)
        do i = max(1, ifv(j,l)), min(ii, ilv(j,l))
          pgfym(i,j,n) = pgfym(i,j,n) + xiyp(i,j,n) - xiym(i,j,n)
          xiyp(i,j,n) = xiyp(i,j,n)/pb_p(i,j  )
          xiym(i,j,n) = xiym(i,j,n)/pb_p(i,j-1)
        end do
      end do

      do l = 1, isp(j)
        do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
          sealv(i,j) = phi(i,j,1)/grav
        end do
      end do

    end do
    !$omp end parallel do

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'pgforc:'
      end if
      call chksummsk(phi   ,ip, kk+1, 'phi')
      call chksummsk(pgfx  ,iu, 2*kk, 'pgfx')
      call chksummsk(pgfy  ,iv, 2*kk, 'pgfy')
      call chksummsk(pgfxm ,iu, 2,    'pgfxm')
      call chksummsk(pgfym ,iv, 2,    'pgfym')
      call chksummsk(xixp  ,iu, 2,    'xixp')
      call chksummsk(xixm  ,iu, 2,    'xixm')
      call chksummsk(xiyp  ,iv, 2,    'xiyp')
      call chksummsk(xiym  ,iv, 2,    'xiym')
    end if

  end subroutine pgforc

end module mod_pgforc
