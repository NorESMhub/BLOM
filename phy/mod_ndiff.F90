! ------------------------------------------------------------------------------
! Copyright (C) 2024 Mats Bentsen, Mehmet Ilicak, Mariana Vertenstein
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

module mod_ndiff
! ------------------------------------------------------------------------------
! This module contains procedures for solving vertical diffusion equations.
! ------------------------------------------------------------------------------

  use mod_types,     only: r8
  use mod_constants, only: grav, alpha0, epsilp, onemm, P_mks2cgs, R_mks2cgs
  use mod_time,      only: delt1
  use mod_xc
  use mod_grid,      only: scuy, scvx, scp2, scuxi, scvyi
  use mod_eos,       only: drhodt, drhods, rho
  use mod_state,     only: dp, temp, saln, utflx, vtflx, usflx, vsflx, pu, pv
  use mod_diffusion, only: difiso, utflld, vtflld, usflld, vsflld
  use mod_cmnfld,    only: nslpx, nslpy
  use mod_tracers,   only: trc

  implicit none
  private

  real(r8), parameter :: &
       ndiff_dstsnp_fac = .01_r8, &
       rho_eps = 1.e-5_r8*R_mks2cgs, &
       dp_eps = 1.e-5_r8*P_mks2cgs
  integer, parameter :: &
       it = 1, &
       is = 2

  integer :: ntr_loc

  ! Public routines
  public :: ndiff_prep_jslice
  public :: ndiff_uflx_jslice
  public :: ndiff_vflx_jslice
  public :: ndiff_update_trc_jslice

contains

  ! ----------------------------------------------------------------------------
  ! Private procedures.
  ! ----------------------------------------------------------------------------

  pure function peval(pc, x) result(f)
  ! ----------------------------------------------------------------------------
  ! Evaluate the polynomial with coefficients pc at x.
  ! ----------------------------------------------------------------------------

    real(r8), dimension(:), intent(in) :: pc
    real(r8), intent(in) :: x

    real(r8) :: f

    f = (((pc(5)*x + pc(4))*x + pc(3))*x + pc(2))*x + pc(1)

  end function peval

  pure function pmeval(pc, x0, x1) result(f)
  ! ----------------------------------------------------------------------------
  ! Evaluate the mean of a polynomial with coefficients pc over the interval
  ! (x0, x1).
  ! ----------------------------------------------------------------------------

    real(r8), dimension(:), intent(in) :: pc
    real(r8), intent(in) :: x0, x1

    real(r8) :: f

    real(r8), parameter :: &
         c1_2 = 1._r8/2._r8, &
         c1_3 = 1._r8/3._r8, &
         c1_4 = 1._r8/4._r8, &
         c1_5 = 1._r8/5._r8

    real(r8) :: b1, b2, b3, b4, b5

    b5 =         c1_5*pc(5)
    b4 = b5*x1 + c1_4*pc(4)
    b3 = b4*x1 + c1_3*pc(3)
    b2 = b3*x1 + c1_2*pc(2)
    b1 = b2*x1 +      pc(1)
    f = (((b5*x0 + b4)*x0 + b3)*x0 + b2)*x0 + b1

  end function pmeval

  pure function drhoroot(tpc, spc, tf, sf, &
                         drhodt_l, drhodt_u, drhods_l, drhods_u) result(x)

    real(r8), dimension(:), intent(in) :: tpc, spc
    real(r8), intent(in) :: tf, sf, drhodt_l, drhodt_u, drhods_l, drhods_u

    real(r8) :: x

    real(r8), parameter :: &
         c1_2 = 1._r8/2._r8, &
         c0 = 0._r8, &
         c1 = 1._r8, &
         c2 = 2._r8, &
         c3 = 3._r8, &
         c4 = 4._r8, &
         eps = 1.e-14_r8, &
         x_tol = 1.e-4_r8

    real(r8) :: ddrdtdx, ddrdsdx, dt, ds, drdt, drds, dtdx, dsdx, dr, ddrdx, &
                x_old
    integer :: n

    x = c1_2
    ddrdtdx = drhodt_l - drhodt_u
    ddrdsdx = drhods_l - drhods_u

    do n = 1, 10

      dt = tf - (tpc(1) + (tpc(2) + (tpc(3) + (tpc(4) + tpc(5)*x)*x)*x)*x)
      ds = sf - (spc(1) + (spc(2) + (spc(3) + (spc(4) + spc(5)*x)*x)*x)*x)
      drdt = drhodt_l*x + drhodt_u*(c1 - x)
      drds = drhods_l*x + drhods_u*(c1 - x)
      dtdx = - (tpc(2) + (c2*tpc(3) + (c3*tpc(4) + c4*tpc(5)*x)*x)*x)
      dsdx = - (spc(2) + (c2*spc(3) + (c3*spc(4) + c4*spc(5)*x)*x)*x)

      dr = drdt*dt + drds*ds
      ddrdx = ddrdtdx*dt + drdt*dtdx + ddrdsdx*ds + drds*dsdx

      x_old = x
      x = max(c0, min(c1, x_old - dr/sign(max(eps, abs(ddrdx)), ddrdx)))
      if (abs(x - x_old) < x_tol) return

    enddo

  end function drhoroot

  pure function drho(t1, s1, t2, s2, drhodt, drhods) result(dr)

    real(r8), intent(in) :: t1, s1, t2, s2, drhodt, drhods

    real(r8) :: dr

    dr = drhodt*(t2 - t1) + drhods*(s2 - s1)

  end function drho

  pure subroutine ndiff_flx(p_srcdi_m, t_srcdi_m, tpc_src_m, &
                            drhodt_srcdi_m, drhods_srcdi_m, &
                            p_dst_m, ksmx_m, kdmx_m, &
                            p_srcdi_p, t_srcdi_p, tpc_src_p, &
                            drhodt_srcdi_p, drhods_srcdi_p, &
                            p_dst_p, ksmx_p, kdmx_p, &
                            cdiff, cnslp, puv, flxconv_js, &
                            uvtflld, uvsflld, uvtflx, uvsflx, nslpxy, &
                            ntr_loc, i_m, j_m, i_p, j_p, js_m, js_p, mm, nn)

    real(r8), dimension(:,:), intent(in) :: &
         p_srcdi_m, drhodt_srcdi_m, drhods_srcdi_m, &
         p_srcdi_p, drhodt_srcdi_p, drhods_srcdi_p
    real(r8), dimension(:,:,:), intent(in) :: &
         t_srcdi_m, tpc_src_m, t_srcdi_p, tpc_src_p
    real(r8), dimension(:), intent(in) :: &
         p_dst_m, p_dst_p
    real(r8), dimension(1-nbdy:,1-nbdy:,:), intent(in) :: &
         puv
    real(r8), dimension(:,:,1-nbdy:,:), intent(inout) :: &
         flxconv_js
    real(r8), dimension(1-nbdy:,1-nbdy:,:), intent(inout) :: &
         uvtflld, uvsflld, uvtflx, uvsflx
    real(r8), dimension(1-nbdy:,1-nbdy:,:), intent(out) :: &
         nslpxy
    real(r8), intent(in) :: cdiff, cnslp
    integer, intent(in) :: &
         ksmx_m, ksmx_p, kdmx_m, kdmx_p, ntr_loc, i_m, j_m, i_p, j_p, &
         js_m, js_p, mm, nn

    real(r8), dimension(4*(kk+1)) :: nslp_src, p_nslp_src
    real(r8), dimension(2,kk) :: p_ni_srcdi_m, p_ni_srcdi_p
    real(r8), dimension(kk+1) :: p_dstsnp_m, p_dstsnp_p
    real(r8), dimension(ntr_loc,2) :: t_ni_m, t_ni_p
    real(r8), dimension(ntr_loc) :: t_nl_m, t_nl_p
    real(r8), dimension(2) :: x_ni_m, x_ni_p, p_ni_m, p_ni_p
    real(r8) :: drho_curr, p_ni_m_prev, p_ni_p_prev, &
         drhodt_x0, drhodt_x1, drhods_x0, drhods_x1, &
         x_ni, p_ni, drho_prev, dp_dst_u, dp_dst_l, pu_m, pl_m, pu_p, pl_p, &
         pp1, pp2, dp_ni_m, dp_ni_p, dp_ni, q, dt, ds, &
         tflx, sflx, p_ni_up, p_ni_lo, dp_ni_i, mlfrac, p_nslp_dst
    integer :: nns, is_m, is_p, ks_m, ks_p, k, ks_m_prev, ks_p_prev, &
         kd_m, kd_p, isn_m, isn_p, ksn_m, ksn_p, &
         nip, nic, kuv, case_m, case_p, nt, kuvm, kd, ks
    logical, dimension(kk) :: stab_src_m, stab_src_p
    logical :: drho_neg, drho_pos, drho_zero, &
         advance_src_m, advance_src_p, advance_dst_m, advance_dst_p, &
         found_ni

    real(r8), parameter :: mval = 1.e30_r8

    ! ------------------------------------------------------------------------
    ! Search the source columns, starting from the surface, to identify
    ! neutral interfaces anchored at layer interfaces. Store information
    ! about whether layers are stably stratified with current reconstruction
    ! and neutral slopes.
    ! ------------------------------------------------------------------------

    p_ni_srcdi_m(:,:) = mval
    p_ni_srcdi_p(:,:) = mval
    stab_src_m(:) = .false.
    stab_src_p(:) = .false.

    nns = 0

    is_m = 1
    ks_m = 1
    ks_p = 1
    is_p = 1
    drho_curr = drho(t_srcdi_m(is_m,ks_m,it), &
                     t_srcdi_m(is_m,ks_m,is), &
                     t_srcdi_p(is_p,ks_p,it), &
                     t_srcdi_p(is_p,ks_p,is), &
                    .5_r8*( drhodt_srcdi_m(is_m,ks_m) &
                          + drhodt_srcdi_p(is_p,ks_p)), &
                    .5_r8*( drhods_srcdi_m(is_m,ks_m) &
                          + drhods_srcdi_p(is_p,ks_p)))
    p_ni_m_prev = p_srcdi_m(1,1)
    p_ni_p_prev = p_srcdi_p(1,1)

    search_loop1: do

      drho_neg = drho_curr <= - rho_eps
      drho_pos = drho_curr >=   rho_eps
      drho_zero = .not. (drho_neg .or. drho_pos)

      if (is_m + ks_m > 2 .and. is_p + ks_p > 2) then
        if (drho_neg) then
          if (is_m == 2) then
            drhodt_x0 = .5_r8*( drhodt_srcdi_m(1   ,ks_m) &
                              + drhodt_srcdi_p(is_p,ks_p))
            drhodt_x1 = .5_r8*( drhodt_srcdi_m(2   ,ks_m) &
                              + drhodt_srcdi_p(is_p,ks_p))
            drhods_x0 = .5_r8*( drhods_srcdi_m(1   ,ks_m) &
                              + drhods_srcdi_p(is_p,ks_p))
            drhods_x1 = .5_r8*( drhods_srcdi_m(2   ,ks_m) &
                             + drhods_srcdi_p(is_p,ks_p))
            x_ni = drhoroot(tpc_src_m(:,ks_m,it), tpc_src_m(:,ks_m,is), &
                            t_srcdi_p(is_p,ks_p,it), t_srcdi_p(is_p,ks_p,is), &
                            drhodt_x1, drhodt_x0, drhods_x1, drhods_x0)
            p_ni = p_srcdi_m(2,ks_m)*x_ni + p_srcdi_m(1,ks_m)*(1._r8 - x_ni)
            if (p_ni > p_ni_m_prev) then
              p_ni_m_prev = p_ni
              p_ni_srcdi_p(is_p,ks_p) = p_ni
              nns = nns + 1
              nslp_src(nns) = - cnslp*(p_srcdi_p(is_p,ks_p) - p_ni)
              p_nslp_src(nns) = .5_r8*(p_srcdi_p(is_p,ks_p) + p_ni)
            endif
          endif
        elseif (drho_pos) then
          if (is_p == 2) then
            drhodt_x0 = .5_r8*( drhodt_srcdi_m(is_m,ks_m) &
                              + drhodt_srcdi_p(1   ,ks_p))
            drhodt_x1 = .5_r8*( drhodt_srcdi_m(is_m,ks_m) &
                              + drhodt_srcdi_p(2   ,ks_p))
            drhods_x0 = .5_r8*( drhods_srcdi_m(is_m,ks_m) &
                              + drhods_srcdi_p(1   ,ks_p))
            drhods_x1 = .5_r8*( drhods_srcdi_m(is_m,ks_m) &
                              + drhods_srcdi_p(2   ,ks_p))
            x_ni = drhoroot(tpc_src_p(:,ks_p,it), tpc_src_p(:,ks_p,is), &
                            t_srcdi_m(is_m,ks_m,it), t_srcdi_m(is_m,ks_m,is), &
                            drhodt_x1, drhodt_x0, drhods_x1, drhods_x0)
            p_ni = p_srcdi_p(2,ks_p)*x_ni + p_srcdi_p(1,ks_p)*(1._r8 - x_ni)
            if (p_ni > p_ni_p_prev) then
              p_ni_p_prev = p_ni
              p_ni_srcdi_m(is_m,ks_m) = p_ni
              nns = nns + 1
              nslp_src(nns) = - cnslp*(p_ni - p_srcdi_m(is_m,ks_m))
              p_nslp_src(nns) = .5_r8*(p_ni + p_srcdi_m(is_m,ks_m))
            endif
          endif
        else
          p_ni_srcdi_p(is_p,ks_p) = p_srcdi_m(is_m,ks_m)
          p_ni_srcdi_m(is_m,ks_m) = p_srcdi_p(is_p,ks_p)
          nns = nns + 1
          nslp_src(nns) = - cnslp*( p_srcdi_p(is_p,ks_p) &
                                  - p_srcdi_m(is_m,ks_m))
          p_nslp_src(nns) = .5_r8*( p_srcdi_p(is_p,ks_p) &
                                  + p_srcdi_m(is_m,ks_m))
        endif
      endif

      if (drho_zero .or. drho_pos) then
        do
          drho_prev = drho_curr
          if (is_m == 1) then
            is_m = 2
          else
            ks_m = ks_m + 1
            if (ks_m > ksmx_m) exit search_loop1
            is_m = 1
          endif
          drho_curr = drho(t_srcdi_m(is_m,ks_m,it), &
                           t_srcdi_m(is_m,ks_m,is), &
                           t_srcdi_p(is_p,ks_p,it), &
                           t_srcdi_p(is_p,ks_p,is), &
                           .5_r8*( drhodt_srcdi_m(is_m,ks_m) &
                                 + drhodt_srcdi_p(is_p,ks_p)), &
                           .5_r8*( drhods_srcdi_m(is_m,ks_m) &
                                 + drhods_srcdi_p(is_p,ks_p)))
          if (drho_prev - drho_curr > rho_eps) then
            if (is_m == 2 .and. p_srcdi_m(2,ks_m) - p_srcdi_m(1,ks_m) > onemm) &
              stab_src_m(ks_m) = .true.
            exit
          endif
          if (is_m == 1) then
            p_ni_srcdi_m(is_m,ks_m) = p_ni_srcdi_m(2,ks_m-1)
          endif
        enddo
      endif

      if (drho_zero .or. drho_neg) then
        do
          drho_prev = drho_curr
          if (is_p == 1) then
            is_p = 2
          else
            ks_p = ks_p + 1
            if (ks_p > ksmx_p) exit search_loop1
            is_p = 1
          endif
          drho_curr = drho(t_srcdi_m(is_m,ks_m,it), &
                           t_srcdi_m(is_m,ks_m,is), &
                           t_srcdi_p(is_p,ks_p,it), &
                           t_srcdi_p(is_p,ks_p,is), &
                           .5_r8*( drhodt_srcdi_m(is_m,ks_m) &
                                 + drhodt_srcdi_p(is_p,ks_p)), &
                           .5_r8*( drhods_srcdi_m(is_m,ks_m) &
                                 + drhods_srcdi_p(is_p,ks_p)))
          if (drho_curr - drho_prev > rho_eps) then
            if (is_p == 2 .and. p_srcdi_p(2,ks_p) - p_srcdi_p(1,ks_p) > onemm) &
              stab_src_p(ks_p) = .true.
            exit
          endif
          if (is_p == 1) then
            p_ni_srcdi_p(is_p,ks_p) = p_ni_srcdi_p(2,ks_p-1)
          endif
        enddo
      endif

    enddo search_loop1

    ! ------------------------------------------------------------------------
    ! Do another search from the surface, this time including destination
    ! interfaces, to identify neutral layers and compute fluxes that are added
    ! to a flux convergence for the target layers.
    ! ------------------------------------------------------------------------

    ! For the search, reset destination interfaces to their corresponding source
    ! interfaces (snap) if the interface difference is less than a specified
    ! fraction of adjacent layer thicknesses.

    p_dstsnp_m(1) = p_dst_m(1)
    dp_dst_u = p_dst_m(2) - p_dst_m(1)
    do k = 2, min(ksmx_m, kdmx_m)
      dp_dst_l = p_dst_m(k+1) - p_dst_m(k)
      if (abs(p_dst_m(k) - p_srcdi_m(1,k)) &
          < min(dp_dst_u, dp_dst_l)*ndiff_dstsnp_fac) then
        p_dstsnp_m(k) = p_srcdi_m(1,k)
      else
        p_dstsnp_m(k) = p_dst_m(k)
      endif
      dp_dst_u = dp_dst_l
    enddo
    do k = min(ksmx_m, kdmx_m) + 1, kdmx_m + 1
      p_dstsnp_m(k) = p_dst_m(k)
    enddo

    p_dstsnp_p(1) = p_dst_p(1)
    dp_dst_u = p_dst_m(2) - p_dst_m(1)
    dp_dst_u = p_dst_p(2) - p_dst_p(1)
    do k = 2, min(ksmx_p, kdmx_p)
      dp_dst_l = p_dst_p(k+1) - p_dst_p(k)
      if (abs(p_dst_p(k) - p_srcdi_p(1,k)) &
          < min(dp_dst_u, dp_dst_l)*ndiff_dstsnp_fac) then
        p_dstsnp_p(k) = p_srcdi_p(1,k)
      else
        p_dstsnp_p(k) = p_dst_p(k)
      endif
      dp_dst_u = dp_dst_l
    enddo
    do k = min(ksmx_p, kdmx_p) + 1, kdmx_p + 1
      p_dstsnp_p(k) = p_dst_p(k)
    enddo

    is_m = 2
    ks_m = 0
    is_p = 2
    ks_p = 0
    kd_m = 0
    kd_p = 0
    advance_src_m = .true.
    advance_src_p = .true.
    advance_dst_m = .true.
    advance_dst_p = .true.
    ks_m_prev = 0
    ks_p_prev = 0
    nip = 1
    nic = 2
    p_ni_m(nip) = - mval
    p_ni_p(nip) = - mval
    kuv = 1

    search_loop2: do

      ! Advance source interface indices as requested. When source interfaces
      ! indices are advanced, keep seperate indices for next interface and next
      ! interface that is anchoring a neutral interface.

      if (advance_src_m) then
        do
          if (is_m == 1) then
            is_m = 2
            if (stab_src_m(ks_m)) exit
          else
            ks_m = ks_m + 1
            if (ks_m > ksmx_m) exit search_loop2
            is_m = 1
            if (stab_src_m(ks_m) .and. p_ni_srcdi_m(is_m,ks_m) /= mval) &
                 exit
          endif
        enddo
        isn_m = is_m
        ksn_m = ks_m
        do while (p_ni_srcdi_m(isn_m,ksn_m) == mval)
          if (isn_m == 1) then
            isn_m = 2
          else
            if (ksn_m == ksmx_m) exit
            ksn_m = ksn_m + 1
            isn_m = 1
          endif
        enddo
      endif

      if (advance_src_p) then
        do
          if (is_p == 1) then
            is_p = 2
            if (stab_src_p(ks_p)) exit
          else
            ks_p = ks_p + 1
            if (ks_p > ksmx_p) exit search_loop2
            is_p = 1
            if (stab_src_p(ks_p) .and. p_ni_srcdi_p(is_p,ks_p) /= mval) then
              exit
            end if
          endif
        enddo
        isn_p = is_p
        ksn_p = ks_p
        do while (p_ni_srcdi_p(isn_p,ksn_p) == mval)
          if (isn_p == 1) then
            isn_p = 2
          else
            if (ksn_p == ksmx_p) exit
            ksn_p = ksn_p + 1
            isn_p = 1
          endif
        enddo
      endif

      ! If previous neutral interface is unset, set it to the first neutral
      ! interface.
      if (p_ni_m(nip) == - mval) then
        if ( (p_ni_srcdi_m(isn_m,ksn_m) - p_srcdi_p(isn_p,ksn_p)) &
           < (p_ni_srcdi_p(isn_p,ksn_p) - p_srcdi_m(isn_m,ksn_m))) then
          p_ni_m(nip) = p_srcdi_m(isn_m,ksn_m)
          p_ni_p(nip) = p_ni_srcdi_m(isn_m,ksn_m)
        else
          p_ni_m(nip) = p_ni_srcdi_p(isn_p,ksn_p)
          p_ni_p(nip) = p_srcdi_p(isn_p,ksn_p)
        endif
      endif

      ! Advance destination interface indices as requested.

      if (advance_dst_m) then
        kd_m = kd_m + 1
        if (kd_m > kdmx_m) exit search_loop2
      endif

      if (advance_dst_p) then
        kd_p = kd_p + 1
        if (kd_p > kdmx_p) exit search_loop2
      endif

      do while (p_dstsnp_m(kd_m+1) &
           <= max(p_srcdi_m(1,ks_m), p_ni_m(nip)))
        kd_m = kd_m + 1
        if (kd_m > kdmx_m) exit search_loop2
      enddo

      do while (p_dstsnp_p(kd_p+1) &
           <= max(p_srcdi_p(1,ks_p), p_ni_p(nip)))
        kd_p = kd_p + 1
        if (kd_p > kdmx_p) exit search_loop2
      enddo

      advance_src_m = .false.
      advance_src_p = .false.
      advance_dst_m = .false.
      advance_dst_p = .false.

      ! By considering current destination interface, source interface and
      ! neutral interface anchored by the neighbour column, find which of
      ! those are the shallowes for each column. Whichever interfaces are
      ! minimums defines cases that are considered to find the next neutral
      ! interface between the columns.

      case_m = 3
      if (p_srcdi_m(is_m,ks_m) <= p_ni_srcdi_p(isn_p,ksn_p)) then
        if (p_srcdi_m(is_m,ks_m) <= p_dstsnp_m(kd_m+1)) case_m = 1
      elseif (p_ni_srcdi_p(isn_p,ksn_p) <= p_dstsnp_m(kd_m+1)) then
        case_m = 2
      endif
      case_p = 3
      if (p_srcdi_p(is_p,ks_p) <= p_ni_srcdi_m(isn_m,ksn_m)) then
        if (p_srcdi_p(is_p,ks_p) <= p_dstsnp_p(kd_p+1)) case_p = 1
      elseif (p_ni_srcdi_m(isn_m,ksn_m) <= p_dstsnp_p(kd_p+1)) then
        case_p = 2
      endif

      found_ni = .false.

      if     (case_m == 3 .and. case_p == 3) then

        if (is_p == 2 .and. is_m == 2) then
          p_ni_m(nic) = p_dstsnp_m(kd_m+1)
          p_ni_p(nic) = p_dstsnp_p(kd_p+1)
          pu_m = p_ni_m(nip)
          pu_p = p_ni_p(nip)
          if ( (p_ni_srcdi_m(isn_m,ksn_m) - p_srcdi_p(isn_p,ksn_p)) &
             < (p_ni_srcdi_p(isn_p,ksn_p) - p_srcdi_m(isn_m,ksn_m))) then
            pl_m = p_srcdi_m(isn_m,ksn_m)
            pl_p = p_ni_srcdi_m(isn_m,ksn_m)
          else
            pl_m = p_ni_srcdi_p(isn_p,ksn_p)
            pl_p = p_srcdi_p(isn_p,ksn_p)
          endif
          pp1 = (p_ni_m(nic) - pu_m)*(pl_p - pu_p)
          pp2 = (p_ni_p(nic) - pu_p)*(pl_m - pu_m)
          if ( abs(pp1 - pp2) &
             < dp_eps*max(dp_eps, pl_m - pu_m + pl_p - pu_p)) then
            advance_dst_m = .true.
            advance_dst_p = .true.
          elseif (pp1 < pp2) then
            p_ni_p(nic) = pu_p + pp1/(pl_m - pu_m)
            advance_dst_m = .true.
          else
            p_ni_m(nic) = pu_m + pp2/(pl_p - pu_p)
            advance_dst_p = .true.
          endif
          if (p_ni_m(nic) >= p_srcdi_m(1,ks_m) .and. &
              p_ni_m(nic) <= p_srcdi_m(2,ks_m) .and. &
              p_ni_p(nic) >= p_srcdi_p(1,ks_p) .and. &
              p_ni_p(nic) <= p_srcdi_p(2,ks_p)) then
            x_ni_m(nic) = (p_ni_m(nic)        - p_srcdi_m(1,ks_m)) &
                          /(p_srcdi_m(2,ks_m) - p_srcdi_m(1,ks_m))
            x_ni_p(nic) = (p_ni_p(nic)        - p_srcdi_p(1,ks_p)) &
                          /(p_srcdi_p(2,ks_p) - p_srcdi_p(1,ks_p))
            do nt = 1, ntr_loc
              t_ni_m(nt,nic) = peval(tpc_src_m(:,ks_m,nt), x_ni_m(nic))
              t_ni_p(nt,nic) = peval(tpc_src_p(:,ks_p,nt), x_ni_p(nic))
            enddo
            found_ni = .true.
          endif
        else
          if (is_p /= 2) advance_dst_m = .true.
          if (is_m /= 2) advance_dst_p = .true.
        endif

      elseif (case_m == 3) then

        if (is_p == 2) then
          p_ni_m(nic) = p_dstsnp_m(kd_m+1)
          if (case_p == 1) then 
            p_ni_p(nic) = p_ni_p(nip) &
                        + (p_ni_m(nic) - p_ni_m(nip)) &
                          *(p_srcdi_p(isn_p,ksn_p) - p_ni_p(nip)) &
                          /(p_ni_srcdi_p(isn_p,ksn_p) - p_ni_m(nip))
          else
            p_ni_p(nic) = p_ni_p(nip) &
                        + (p_ni_m(nic) - p_ni_m(nip)) &
                          *(p_ni_srcdi_m(isn_m,ksn_m) - p_ni_p(nip)) &
                          /(p_srcdi_m(isn_m,ksn_m) - p_ni_m(nip))
          endif
          if (p_ni_p(nic) >= p_srcdi_p(1,ks_p) .and. &
              p_ni_p(nic) <= p_srcdi_p(2,ks_p)) then
            x_ni_m(nic) = (p_dstsnp_m(kd_m+1) - p_srcdi_m(1,ks_m)) &
                          /(p_srcdi_m(2,ks_m) - p_srcdi_m(1,ks_m))
            x_ni_p(nic) = (p_ni_p(nic)        - p_srcdi_p(1,ks_p)) &
                          /(p_srcdi_p(2,ks_p) - p_srcdi_p(1,ks_p))
            do nt = 1, ntr_loc
              t_ni_m(nt,nic) = peval(tpc_src_m(:,ks_m,nt), x_ni_m(nic))
              t_ni_p(nt,nic) = peval(tpc_src_p(:,ks_p,nt), x_ni_p(nic))
            enddo
            found_ni = .true.
            advance_dst_m = .true.
          else
            if (case_p == 1 .and. p_ni_srcdi_p(is_p,ks_p) == mval) then
              advance_src_p = .true.
            else
              advance_dst_m = .true.
            endif
          endif
        else
          advance_dst_m = .true.
        endif

      elseif (case_p == 3) then

        if (is_m == 2) then
          p_ni_p(nic) = p_dstsnp_p(kd_p+1)
          if (case_m == 1) then
            p_ni_m(nic) = p_ni_m(nip) &
                        + (p_ni_p(nic) - p_ni_p(nip)) &
                          *(p_srcdi_m(isn_m,ksn_m) - p_ni_m(nip)) &
                          /(p_ni_srcdi_m(isn_m,ksn_m) - p_ni_p(nip))
          else
            p_ni_m(nic) = p_ni_m(nip) &
                        + (p_ni_p(nic) - p_ni_p(nip)) &
                          *(p_ni_srcdi_p(isn_p,ksn_p) - p_ni_m(nip)) &
                          /(p_srcdi_p(isn_p,ksn_p) - p_ni_p(nip))
          endif
          if (p_ni_m(nic) >= p_srcdi_m(1,ks_m) .and. &
              p_ni_m(nic) <= p_srcdi_m(2,ks_m)) then
            x_ni_p(nic) = (p_dstsnp_p(kd_p+1) - p_srcdi_p(1,ks_p)) &
                          /(p_srcdi_p(2,ks_p) - p_srcdi_p(1,ks_p))
            x_ni_m(nic) = (p_ni_m(nic)        - p_srcdi_m(1,ks_m)) &
                          /(p_srcdi_m(2,ks_m) - p_srcdi_m(1,ks_m))
            do nt = 1, ntr_loc
              t_ni_m(nt,nic) = peval(tpc_src_m(:,ks_m,nt), x_ni_m(nic))
              t_ni_p(nt,nic) = peval(tpc_src_p(:,ks_p,nt), x_ni_p(nic))
            enddo
            found_ni = .true.
            advance_dst_p = .true.
          else
            if (case_m == 1 .and. p_ni_srcdi_m(is_m,ks_m) == mval) then
              advance_src_m = .true.
            else
              advance_dst_p = .true.
            endif
          endif
        else
          advance_dst_p = .true.
        endif

      elseif (case_m == 1 .and. case_p == 1) then

        if (p_ni_srcdi_m(is_m,ks_m) /= mval .and. &
            p_ni_srcdi_p(is_p,ks_p) /= mval) then
          x_ni_m(nic) = real(is_m - 1, r8)
          p_ni_m(nic) = p_srcdi_m(is_m,ks_m)
          x_ni_p(nic) = real(is_p - 1, r8)
          p_ni_p(nic) = p_srcdi_p(is_p,ks_p)
          do nt = 1, ntr_loc
            t_ni_m(nt,nic) = t_srcdi_m(is_m,ks_m,nt)
            t_ni_p(nt,nic) = t_srcdi_p(is_p,ks_p,nt)
          enddo
          found_ni = .true.
          advance_src_m = .true.
          advance_src_p = .true.
        else
          if (p_ni_srcdi_m(is_m,ks_m) == mval) advance_src_m = .true.
          if (p_ni_srcdi_p(is_p,ks_p) == mval) advance_src_p = .true.
        endif

      elseif (case_m == 1) then

        if (p_ni_srcdi_m(is_m,ks_m) /= mval .and. &
            p_ni_srcdi_m(is_m,ks_m) >= p_srcdi_p(1,ks_p)) then
          x_ni_m(nic) = real(is_m - 1, r8)
          p_ni_m(nic) = p_srcdi_m(is_m,ks_m)
          p_ni_p(nic) = p_ni_srcdi_m(is_m,ks_m)
          x_ni_p(nic) = (p_ni_p(nic)       - p_srcdi_p(1,ks_p)) &
                       /(p_srcdi_p(2,ks_p) - p_srcdi_p(1,ks_p))
          do nt = 1, ntr_loc
            t_ni_m(nt,nic) = t_srcdi_m(is_m,ks_m,nt)
            t_ni_p(nt,nic) = peval(tpc_src_p(:,ks_p,nt), x_ni_p(nic))
          enddo
          found_ni = .true.
        endif
        advance_src_m = .true.

      elseif (case_p == 1) then

        if (p_ni_srcdi_p(is_p,ks_p) /= mval .and. &
            p_ni_srcdi_p(is_p,ks_p) >= p_srcdi_m(1,ks_m)) then
          x_ni_p(nic) = real(is_p - 1, r8)
          p_ni_p(nic) = p_srcdi_p(is_p,ks_p)
          p_ni_m(nic) = p_ni_srcdi_p(is_p,ks_p)
          x_ni_m(nic) = (p_ni_m(nic)       - p_srcdi_m(1,ks_m)) &
                       /(p_srcdi_m(2,ks_m) - p_srcdi_m(1,ks_m))
          do nt = 1, ntr_loc
            t_ni_p(nt,nic) = t_srcdi_p(is_p,ks_p,nt)
            t_ni_m(nt,nic) = peval(tpc_src_m(:,ks_m,nt), x_ni_m(nic))
          enddo
          found_ni = .true.
        endif
        advance_src_p = .true.

      else
        advance_src_m = .true.
        advance_src_p = .true.
      endif

      if (found_ni) then

        ! If a neutral interface is found, check 1) whether the current and
        ! previous neutral interfaces are between same source and destination
        ! layers and 2) whether the neutral layer thickness is above a small
        ! threshold. If so, a neutral layer, suitable for diffusive flux
        ! computations, has been found.

        dp_ni_m = min(p_ni_m(nic) - p_ni_m(nip), &
                      p_dst_m(kd_m+1) - p_dst_m(kd_m))
        dp_ni_p = min(p_ni_p(nic) - p_ni_p(nip), &
                      p_dst_p(kd_p+1) - p_dst_p(kd_p))
        dp_ni = 2._r8*dp_ni_m*dp_ni_p/max(dp_ni_m + dp_ni_p, 2._r8*dp_eps)

        if (ks_m == ks_m_prev .and. ks_p == ks_p_prev .and. &
            p_ni_m(nip) >= p_dstsnp_m(kd_m) .and. &
            p_ni_m(nic) <= p_dstsnp_m(kd_m+1) .and. &
            p_ni_p(nip) >= p_dstsnp_p(kd_p) .and. &
            p_ni_p(nic) <= p_dstsnp_p(kd_p+1) .and. &
            dp_ni > 2._r8*dp_eps) then

          do nt = 1, ntr_loc
            t_nl_m(nt) = pmeval(tpc_src_m(:,ks_m,nt), x_ni_m(nip), x_ni_m(nic))
            t_nl_p(nt) = pmeval(tpc_src_p(:,ks_p,nt), x_ni_p(nip), x_ni_p(nic))
          enddo

          q = .5_r8*cdiff*(difiso(i_m,j_m,ks_m) + difiso(i_p,j_p,ks_p))*dp_ni

          dt = t_nl_m(it) - t_nl_p(it)
          ds = t_nl_m(is) - t_nl_p(is)

          if (dt*( temp(i_m,j_m,ks_m+nn) &
                 - temp(i_p,j_p,ks_p+nn)) >= 0._r8 .and. &
              dt*( t_ni_m(it,nip) - t_ni_p(it,nip)) >= 0._r8 .and. &
              dt*( t_ni_m(it,nic) - t_ni_p(it,nic)) >= 0._r8 .and. &
              ds*( saln(i_m,j_m,ks_m+nn) &
                 - saln(i_p,j_p,ks_p+nn)) >= 0._r8 .and. &
              ds*( t_ni_m(is,nip) - t_ni_p(is,nip)) >= 0._r8 .and. &
              ds*( t_ni_m(is,nic) - t_ni_p(is,nic)) >= 0._r8) then
            tflx = q*dt
            flxconv_js(kd_m,it,i_m,js_m) = flxconv_js(kd_m,it,i_m,js_m) + tflx
            flxconv_js(kd_p,it,i_p,js_p) = flxconv_js(kd_p,it,i_p,js_p) - tflx
            sflx = q*ds
            flxconv_js(kd_m,is,i_m,js_m) = flxconv_js(kd_m,is,i_m,js_m) + sflx
            flxconv_js(kd_p,is,i_p,js_p) = flxconv_js(kd_p,is,i_p,js_p) - sflx
            p_ni_up = .5_r8*(p_ni_m(nip) + p_ni_p(nip))
            p_ni_lo = .5_r8*(p_ni_m(nic) + p_ni_p(nic))
            dp_ni_i = 1._r8/max(epsilp, p_ni_lo - p_ni_up)
            do while (kuv <= kk)
              kuvm = kuv + mm
              if (puv(i_p,j_p,kuv+1) < p_ni_lo) then
                mlfrac = max(0._r8, puv(i_p,j_p,kuv+1) &
                       - max(p_ni_up, puv(i_p,j_p,kuv))) *dp_ni_i
                uvtflld(i_p,j_p,kuvm) = uvtflld(i_p,j_p,kuvm) + tflx*mlfrac
                uvsflld(i_p,j_p,kuvm) = uvsflld(i_p,j_p,kuvm) + sflx*mlfrac
                uvtflx(i_p,j_p,kuvm) = uvtflx(i_p,j_p,kuvm) + tflx*mlfrac
                uvsflx(i_p,j_p,kuvm) = uvsflx(i_p,j_p,kuvm) + sflx*mlfrac
                kuv = kuv + 1
              else
                mlfrac = (p_ni_lo - max(p_ni_up, puv(i_p,j_p,kuv))) * dp_ni_i
                uvtflld(i_p,j_p,kuvm) = uvtflld(i_p,j_p,kuvm) + tflx*mlfrac
                uvsflld(i_p,j_p,kuvm) = uvsflld(i_p,j_p,kuvm) + sflx*mlfrac
                uvtflx(i_p,j_p,kuvm) = uvtflx(i_p,j_p,kuvm) + tflx*mlfrac
                uvsflx(i_p,j_p,kuvm) = uvsflx(i_p,j_p,kuvm) + sflx*mlfrac
                exit
              endif
            enddo
          endif

          do nt = 3, ntr_loc
            dt = t_nl_m(nt) - t_nl_p(nt)
            if (dt*( trc(i_m,j_m,ks_m+nn,nt-2) &
                 - trc(i_p,j_p,ks_p+nn,nt-2)) >= 0._r8 .and. &
                 dt*( t_ni_m(nt,nip) - t_ni_p(nt,nip)) >= 0._r8 .and. &
                 dt*( t_ni_m(nt,nic) - t_ni_p(nt,nic)) >= 0._r8) then
              tflx = q*dt
              flxconv_js(kd_m,nt,i_m,js_m) = flxconv_js(kd_m,nt,i_m,js_m) + tflx
              flxconv_js(kd_p,nt,i_p,js_p) = flxconv_js(kd_p,nt,i_p,js_p) - tflx
            endif
          enddo

        endif

        ks_m_prev = ks_m
        ks_p_prev = ks_p
        nip = 3 - nip
        nic = 3 - nic

      endif

    enddo search_loop2

    ! Linearly interpolate the neutral slope estimates from the source data to
    ! destination interfaces.
    if (nns == 0) then
      nslpxy(i_p,j_p,:) = 0._r8
    else
      do kd = 1, kk
        p_nslp_dst = .5_r8*(p_dst_m(kd) + p_dst_p(kd))
        if (p_nslp_dst > p_nslp_src(1)) exit
        nslpxy(i_p,j_p,kd) = nslp_src(1)
      enddo
      ks = 1
      interp_loop: do
        do while (p_nslp_dst > p_nslp_src(ks))
          if (ks == nns) exit interp_loop
          ks = ks + 1
        enddo
        q = (p_nslp_src(ks) - p_nslp_dst) &
             /max(p_nslp_src(ks) - p_nslp_src(ks-1), epsilp)
        nslpxy(i_p,j_p,kd) = q*nslp_src(ks-1) + (1._r8 - q)*nslp_src(ks)
        kd = kd + 1
        if (kd > kk) exit
        p_nslp_dst = .5_r8*(p_dst_m(kd) + p_dst_p(kd))
      enddo interp_loop
      do kd = kd, kk
        nslpxy(i_p,j_p,kd) = nslp_src(nns)
      enddo
    endif

  end subroutine ndiff_flx

  ! ----------------------------------------------------------------------------
  ! Public procedures.
  ! ----------------------------------------------------------------------------

  subroutine ndiff_prep_jslice(p_src_js, ksmx_js, &
                               tpc_src_js, t_srcdi_js, &
                               p_dst_js, kdmx_js, p_srcdi_js, &
                               drhodt_srcdi_js, drhods_srcdi_js, &
                               flxconv_js, &
                               ilb, iub, j, js, mm)

    real(r8), dimension(:,1-nbdy:,:), intent(in) :: p_src_js, p_dst_js
    integer, dimension(1-nbdy:,:), intent(in) :: ksmx_js
    integer, dimension(1-nbdy:,:), intent(out) :: kdmx_js
    real(r8), dimension(:,:,:,1-nbdy:,:), intent(in) :: tpc_src_js, t_srcdi_js
    real(r8), dimension(:,:,1-nbdy:,:), intent(out) :: &
      p_srcdi_js, drhodt_srcdi_js, drhods_srcdi_js, flxconv_js
    integer, intent(in) :: ilb, iub, j, js, mm

    integer :: l, i, nt, k, km, errstat

    do l = 1, isp(j)
      do i = max(ilb, ifp(j, l)), min(iub, ilp(j, l))

        ! Find index of deepest destination layer with non-zero thickness.
        kdmx_js(i,js) = kk
        do k = kk, 1, -1
          if (p_dst_js(k,i,js) == p_dst_js(kk+1,i,js)) &
            kdmx_js(i,js) = k - 1
        enddo

        ! Store variables in dual interface arrays with with values
        ! corresponding to upper and lower interface of each layer.
        do k = 1, ksmx_js(i,js)
          p_srcdi_js(1,k,i,js) = p_src_js(k  ,i,js)
          p_srcdi_js(2,k,i,js) = p_src_js(k+1,i,js)
          drhodt_srcdi_js(1,k,i,js) = drhodt(p_srcdi_js(1,k   ,i,js), &
                                             t_srcdi_js(1,k,it,i,js), &
                                             t_srcdi_js(1,k,is,i,js))
          drhodt_srcdi_js(2,k,i,js) = drhodt(p_srcdi_js(2,k   ,i,js), &
                                             t_srcdi_js(2,k,it,i,js), &
                                             t_srcdi_js(2,k,is,i,js))
          drhods_srcdi_js(1,k,i,js) = drhods(p_srcdi_js(1,k   ,i,js), &
                                             t_srcdi_js(1,k,it,i,js), &
                                             t_srcdi_js(1,k,is,i,js))
          drhods_srcdi_js(2,k,i,js) = drhods(p_srcdi_js(2,k   ,i,js), &
                                             t_srcdi_js(2,k,it,i,js), &
                                             t_srcdi_js(2,k,is,i,js))
        enddo

        flxconv_js(:,:,i,js) = 0._r8

      enddo
    enddo

    do k = 1, kk
      km = k + mm
      do l = 1, isu(j)
        do i = max(ilb, ifu(j, l)), min(iub, ilu(j, l))
          utflld(i,j,km) = 0._r8
          usflld(i,j,km) = 0._r8
        enddo
      enddo
      do l = 1, isv(j)
        do i = max(ilb, ifv(j, l)), min(iub, ilv(j, l))
          vtflld(i,j,km) = 0._r8
          vsflld(i,j,km) = 0._r8
        enddo
      enddo
    enddo

  end subroutine ndiff_prep_jslice

  subroutine ndiff_uflx_jslice(ksmx_js, tpc_src_js, t_srcdi_js, &
                               p_dst_js, kdmx_js, p_srcdi_js, &
                               drhodt_srcdi_js, drhods_srcdi_js, &
                               flxconv_js, &
                               ntr_loc, ilb, iub, j, js, mm, nn)

    integer, dimension(1-nbdy:,:), intent(in) :: ksmx_js, kdmx_js
    real(r8), dimension(:,:,:,1-nbdy:,:), target, intent(in) :: &
         tpc_src_js, t_srcdi_js
    real(r8), dimension(:,1-nbdy:,:), target, intent(in) :: p_dst_js
    real(r8), dimension(:,:,1-nbdy:,:), target, intent(in) :: &
         p_srcdi_js, drhodt_srcdi_js, drhods_srcdi_js
    real(r8), dimension(:,:,1-nbdy:,:), intent(inout) :: flxconv_js
    integer, intent(in) :: ntr_loc, ilb, iub, j, js, mm, nn

    real(r8), dimension(:,:,:), pointer :: &
         t_srcdi_m, tpc_src_m, t_srcdi_p, tpc_src_p
    real(r8), dimension(:,:), pointer :: &
         p_srcdi_m, drhodt_srcdi_m, drhods_srcdi_m, &
         p_srcdi_p, drhodt_srcdi_p, drhods_srcdi_p
    real(r8), dimension(:), pointer :: &
         p_dst_m, p_dst_p
    real(r8) :: cdiff, cnslp
    integer :: l, i, ksmx_m, ksmx_p, kdmx_m, kdmx_p

    do l = 1, isu(j)
      do i = max(ilb, ifu(j, l)), min(iub, ilu(j, l))

        p_srcdi_m => p_srcdi_js(:,:,i-1,js)
        p_srcdi_p => p_srcdi_js(:,:,i  ,js)
        t_srcdi_m => t_srcdi_js(:,:,:,i-1,js)
        t_srcdi_p => t_srcdi_js(:,:,:,i  ,js)
        tpc_src_m => tpc_src_js(:,:,:,i-1,js)
        tpc_src_p => tpc_src_js(:,:,:,i  ,js)
        drhodt_srcdi_m => drhodt_srcdi_js(:,:,i-1,js)
        drhodt_srcdi_p => drhodt_srcdi_js(:,:,i  ,js)
        drhods_srcdi_m => drhods_srcdi_js(:,:,i-1,js)
        drhods_srcdi_p => drhods_srcdi_js(:,:,i  ,js)
        p_dst_m => p_dst_js(:,i-1,js)
        p_dst_p => p_dst_js(:,i  ,js)
        ksmx_m = ksmx_js(i-1,js)
        ksmx_p = ksmx_js(i  ,js)
        kdmx_m = kdmx_js(i-1,js)
        kdmx_p = kdmx_js(i  ,js)
        cdiff = delt1*scuy(i,j)*scuxi(i,j)
        cnslp = alpha0*scuxi(i,j)/grav

        call ndiff_flx(p_srcdi_m, t_srcdi_m, tpc_src_m, &
                       drhodt_srcdi_m, drhods_srcdi_m, &
                       p_dst_m, ksmx_m, kdmx_m, &
                       p_srcdi_p, t_srcdi_p, tpc_src_p, &
                       drhodt_srcdi_p, drhods_srcdi_p, &
                       p_dst_p, ksmx_p, kdmx_p, &
                       cdiff, cnslp, pu, flxconv_js, &
                       utflld, usflld, utflx, usflx, nslpx, &
                       ntr_loc, i-1, j, i, j, js, js, mm, nn)

      enddo
    enddo

  end subroutine ndiff_uflx_jslice

  subroutine ndiff_vflx_jslice(ksmx_js, tpc_src_js, t_srcdi_js, &
                               p_dst_js, kdmx_js, p_srcdi_js, &
                               drhodt_srcdi_js, drhods_srcdi_js, &
                               flxconv_js, &
                               ntr_loc, ilb, iub, j, js_m, js_p, mm, nn)

    integer, dimension(1-nbdy:,:), intent(in) :: ksmx_js, kdmx_js
    real(r8), dimension(:,:,:,1-nbdy:,:), target, intent(in) :: &
         tpc_src_js, t_srcdi_js
    real(r8), dimension(:,1-nbdy:,:), target, intent(in) :: p_dst_js
    real(r8), dimension(:,:,1-nbdy:,:), target, intent(in) :: &
         p_srcdi_js, drhodt_srcdi_js, drhods_srcdi_js
    real(r8), dimension(:,:,1-nbdy:,:), intent(inout) :: flxconv_js
    integer, intent(in) :: ntr_loc, ilb, iub, j, js_m, js_p, mm, nn

    real(r8), dimension(:,:,:), pointer :: &
         t_srcdi_m, tpc_src_m, t_srcdi_p, tpc_src_p
    real(r8), dimension(:,:), pointer :: &
         p_srcdi_m, drhodt_srcdi_m, drhods_srcdi_m, &
         p_srcdi_p, drhodt_srcdi_p, drhods_srcdi_p
    real(r8), dimension(:), pointer :: &
         p_dst_m, p_dst_p
    real(r8) :: cdiff, cnslp
    integer :: l, i, ksmx_m, ksmx_p, kdmx_m, kdmx_p

    do l = 1, isv(j)
      do i = max(ilb, ifv(j, l)), min(iub, ilv(j, l))

        p_srcdi_m => p_srcdi_js(:,:,i,js_m)
        p_srcdi_p => p_srcdi_js(:,:,i,js_p)
        t_srcdi_m => t_srcdi_js(:,:,:,i,js_m)
        t_srcdi_p => t_srcdi_js(:,:,:,i,js_p)
        tpc_src_m => tpc_src_js(:,:,:,i,js_m)
        tpc_src_p => tpc_src_js(:,:,:,i,js_p)
        drhodt_srcdi_m => drhodt_srcdi_js(:,:,i,js_m)
        drhodt_srcdi_p => drhodt_srcdi_js(:,:,i,js_p)
        drhods_srcdi_m => drhods_srcdi_js(:,:,i,js_m)
        drhods_srcdi_p => drhods_srcdi_js(:,:,i,js_p)
        p_dst_m => p_dst_js(:,i,js_m)
        p_dst_p => p_dst_js(:,i,js_p)
        ksmx_m = ksmx_js(i,js_m)
        ksmx_p = ksmx_js(i,js_p)
        kdmx_m = kdmx_js(i,js_m)
        kdmx_p = kdmx_js(i,js_p)
        cdiff = delt1*scvx(i,j)*scvyi(i,j)
        cnslp = alpha0*scvyi(i,j)/grav

        call ndiff_flx(p_srcdi_m, t_srcdi_m, tpc_src_m, &
                       drhodt_srcdi_m, drhods_srcdi_m, &
                       p_dst_m, ksmx_m, kdmx_m, &
                       p_srcdi_p, t_srcdi_p, tpc_src_p, &
                       drhodt_srcdi_p, drhods_srcdi_p, &
                       p_dst_p, ksmx_p, kdmx_p, &
                       cdiff, cnslp, pv, flxconv_js, &
                       vtflld, vsflld, vtflx, vsflx, nslpy, &
                       ntr_loc, i, j-1, i, j, js_m, js_p, mm, nn)

      enddo
    enddo

  end subroutine ndiff_vflx_jslice

  pure subroutine ndiff_update_trc_jslice(p_dst_js, flxconv_js, trc_rm, &
                                          ntr_loc, ilb, iub, j, js)

    real(r8), dimension(:,1-nbdy:,:), intent(in) :: p_dst_js
    real(r8), dimension(:,:,1-nbdy:,:), intent(in) :: flxconv_js
    real(r8), dimension(:,:,1-nbdy:), intent(inout) :: trc_rm
    integer, intent(in) :: ntr_loc, ilb, iub, j, js

    real(r8) :: q
    integer :: k, l, i, nt

    do l = 1, isp(j)
      do i = max(ilb, ifp(j, l)), min(iub, ilp(j, l))
        do k = 1, kk
          q = 1._r8/(scp2(i,j)*max( p_dst_js(k+1,i,js) &
                                  - p_dst_js(k  ,i,js), dp_eps))
          do nt = 1, ntr_loc
            trc_rm(k,nt,i) = trc_rm(k,nt,i) - q*flxconv_js(k,nt,i,js)
          enddo
        enddo
      enddo
    enddo

  end subroutine ndiff_update_trc_jslice

end module mod_ndiff
