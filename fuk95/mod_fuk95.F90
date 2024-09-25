! ------------------------------------------------------------------------------
! Copyright (C) 2021-2024 Mats Bentsen
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

module mod_fuk95
  ! ------------------------------------------------------------------------------
  ! This module contains variables and routines related to a periodic channel
  ! experiment with configuration and initial conditions from Fukamachi et al.
  ! (1995) and further elaborated in Eldevik and Dysthe (2002).
  ! ------------------------------------------------------------------------------

  use mod_types,     only: r8
  use mod_constants, only: g, rearth, rho0, pi, radian, epsilz, &
                           L_mks2cgs, R_mks2cgs
  use mod_xc
  use mod_vcoord,    only: vcoord_tag, vcoord_isopyc_bulkml, sigmar
  use mod_grid,      only: qclon, qclat, pclon, pclat, uclon, uclat, vclon, vclat, &
                           scqx, scqy, scpx, scpy, scux, scuy, scvx, scvy, &
                           scq2, scp2, scu2, scv2, &
                           qlon, qlat, plon, plat, ulon, ulat, vlon, vlat, &
                           depths, corioq, coriop, betafp, angle, cosang, sinang, &
                           nwp
  use mod_eos,       only: tofsig
  use mod_forcing,   only: surflx, surrlx, sswflx, salflx, brnflx, salrlx, &
                           taux, tauy
  use mod_mxlayr,    only: mltmin
  use mod_state,     only: v, temp, saln, sigma, phi
  use mod_checksum,  only: csdiag, chksummsk

  implicit none

  private

  !  real(r8), parameter :: &
  !     u0     = 30._r8, &     ! Maximum jet velocity [cm s-1].
  !     h1     = 1.e4_r8, &    ! Depth of active layer [cm].
  !     h0     = 2.e4_r8, &    ! Depth of water column [cm].
  !     l0     = 2.e6_r8, &    ! Half-width of the jet [cm].
  !     drho   = 0.19e-3_r8, & ! Active layer density difference [g cm-3].
  !     rhoc   = 1.0259_r8, &  ! Density at the center of active layer [g cm-3].
  !     rhob   = 1.0270_r8, &  ! Density of water beneath active layer [g cm-3].
  !     f      = 1.e-4_r8, &   ! Coriolis parameter [1 s-1].
  !     lat0   = 45._r8, &     ! Center latitude of grid domain [deg].
  !     lambda = 20.8e5, &     ! Channel length [cm].
  !     mindz  = 1.e2_r8, &    ! Minimum interior layer thickness [cm].
  !     saln0  = 35._r8        ! Constant salinity value [g kg-1].
  real(r8), parameter :: &
       u0     = .3_r8*L_mks2cgs, &     ! Maximum jet velocity [m s-1].
       h1     = 1.e2_r8*L_mks2cgs, &   ! Depth of active layer [m].
       h0     = 2.e2_r8*L_mks2cgs, &   ! Depth of water column [m].
       l0     = 2.e4_r8*L_mks2cgs, &   ! Half-width of the jet [m].
       drho   = 0.19_r8*R_mks2cgs, &   ! Active layer density difference [kg m-3].
       rhoc   = 1025.9_r8*R_mks2cgs, & ! Density at the center of active layer [kg m-3].
       rhob   = 1027.0_r8*R_mks2cgs, & ! Density of water beneath active layer [kg m-3].
       f      = 1.e-4_r8, &            ! Coriolis parameter [1 s-1].
       lat0   = 45._r8, &              ! Center latitude of grid domain [deg].
       lambda = 20.8e3*L_mks2cgs, &    ! Channel length [m].
       mindz  = 1._r8*L_mks2cgs, &     ! Minimum interior layer thickness [m].
       saln0  = 35._r8                 ! Constant salinity value [g kg-1].

  public :: geoenv_fuk95, inifrc_fuk95, ictsz_fuk95

contains

  ! ---------------------------------------------------------------------------
  ! Private procedures.
  ! ---------------------------------------------------------------------------

  pure real(r8) function x_nudge(ri, rj)
    ! ---------------------------------------------------------------------------
    ! Perturbed x-position as a function of grid location.
    ! ---------------------------------------------------------------------------

    real(r8), intent(in) :: ri, rj

    x_nudge = ( ri + i0 - itdm/2 - .5_r8 &
         + .1_r8*sin(2._r8*(rj + j0 - 1)*pi/jtdm))*lambda/jtdm

  end function x_nudge

  pure real(r8) function psi(x)
    ! ---------------------------------------------------------------------------
    ! Shape function of geostrophic jet.
    ! ---------------------------------------------------------------------------

    real(r8), intent(in) :: x

    if (abs(x) >= l0) then
      psi = 0._r8
    else
      psi = .5_r8*(1._r8 + cos(pi*x/l0))
    endif

  end function psi

  pure real(r8) function x_psi(x)
    ! ---------------------------------------------------------------------------
    ! Integral of shape function with respect to x.
    ! ---------------------------------------------------------------------------

    real(r8), intent(in) :: x

    if     (x <= -l0) then
      x_psi = -.5_r8*l0
    elseif (x >= l0) then
      x_psi = .5_r8*l0
    else
      x_psi = .5_r8*(x + l0/pi*sin(pi*x/l0))
    endif

  end function x_psi

  ! ---------------------------------------------------------------------------
  ! Public procedures.
  ! ---------------------------------------------------------------------------

  subroutine geoenv_fuk95
    ! ---------------------------------------------------------------------------
    ! Define bathymetry, grid specification and Coriolis parameter for the
    ! Fukamachi et al. (1995) experiment.
    ! ---------------------------------------------------------------------------

    real(r8), dimension(itdm,jtdm) :: tmpg
    real(r8) :: gs, dlat, dlon
    integer i, j

    ! Set water depth.
    if (mnproc == 1) then
      !$omp parallel do private(i)
      do j = 1, jtdm
        tmpg(1   , j) = 0._r8
        tmpg(itdm, j) = 0._r8
        do i = 2, itdm - 1
          tmpg(i, j) = h0*L_mks2cgs**(-1)
        enddo
      enddo
      !$omp end parallel do
    endif
    call xcaput(tmpg, depths, 1)

    ! Number of wet points.
    nwp = (itdm - 2)*jtdm

    ! Grid spacing.
    gs = lambda/jtdm

    ! Set grid coordinates, Coriolis parameter and local grid angle.

    dlat = gs*radian/rearth
    dlon = dlat*sin(lat0/radian)

    !$omp parallel do private(i)
    do j = 1, jj
      do i = 1, ii

        qlon(i, j) = (j + j0        )*dlon
        qlat(i, j) = (i + i0 - itdm/2 - .5_r8)*dlat + lat0

        plon(i, j) = (j + j0 + .5_r8)*dlon
        plat(i, j) = (i + i0 - itdm/2        )*dlat + lat0

        ulon(i, j) = (j + j0        )*dlon
        ulat(i, j) = (i + i0 - itdm/2        )*dlat + lat0

        vlon(i, j) = (j + j0 + .5_r8)*dlon
        vlat(i, j) = (i + i0 - itdm/2 - .5_r8)*dlat + lat0

        qclon(i, j, 1) = (j + j0 - .5_r8)*dlon
        qclon(i, j, 2) = (j + j0 + .5_r8)*dlon
        qclon(i, j, 3) = (j + j0 + .5_r8)*dlon
        qclon(i, j, 4) = (j + j0 - .5_r8)*dlon
        qclat(i, j, 1) = (i + i0 - itdm/2 - 1._r8)*dlat + lat0
        qclat(i, j, 2) = (i + i0 - itdm/2 - 1._r8)*dlat + lat0
        qclat(i, j, 3) = (i + i0 - itdm/2        )*dlat + lat0
        qclat(i, j, 4) = (i + i0 - itdm/2        )*dlat + lat0

        pclon(i, j, 1) = (j + j0        )*dlon
        pclon(i, j, 2) = (j + j0 + 1._r8)*dlon
        pclon(i, j, 3) = (j + j0 + 1._r8)*dlon
        pclon(i, j, 4) = (j + j0        )*dlon
        pclat(i, j, 1) = (i + i0 - itdm/2 - .5_r8)*dlat + lat0
        pclat(i, j, 2) = (i + i0 - itdm/2 - .5_r8)*dlat + lat0
        pclat(i, j, 3) = (i + i0 - itdm/2 + .5_r8)*dlat + lat0
        pclat(i, j, 4) = (i + i0 - itdm/2 + .5_r8)*dlat + lat0

        uclon(i, j, 1) = (j + j0 - .5_r8)*dlon
        uclon(i, j, 2) = (j + j0 + .5_r8)*dlon
        uclon(i, j, 3) = (j + j0 + .5_r8)*dlon
        uclon(i, j, 4) = (j + j0 - .5_r8)*dlon
        uclat(i, j, 1) = (i + i0 - itdm/2 - .5_r8)*dlat + lat0
        uclat(i, j, 2) = (i + i0 - itdm/2 - .5_r8)*dlat + lat0
        uclat(i, j, 3) = (i + i0 - itdm/2 + .5_r8)*dlat + lat0
        uclat(i, j, 4) = (i + i0 - itdm/2 + .5_r8)*dlat + lat0

        vclon(i, j, 1) = (j + j0        )*dlon
        vclon(i, j, 2) = (j + j0 + 1._r8)*dlon
        vclon(i, j, 3) = (j + j0 + 1._r8)*dlon
        vclon(i, j, 4) = (j + j0        )*dlon
        vclat(i, j, 1) = (i + i0 - itdm/2 - 1._r8)*dlat + lat0
        vclat(i, j, 2) = (i + i0 - itdm/2 - 1._r8)*dlat + lat0
        vclat(i, j, 3) = (i + i0 - itdm/2        )*dlat + lat0
        vclat(i, j, 4) = (i + i0 - itdm/2        )*dlat + lat0

        scqx(i, j) = gs
        scqy(i, j) = gs
        scpx(i, j) = gs
        scpy(i, j) = gs
        scux(i, j) = gs
        scuy(i, j) = gs
        scvx(i, j) = gs
        scvy(i, j) = gs
        scq2(i, j) = gs*gs
        scp2(i, j) = gs*gs
        scu2(i, j) = gs*gs
        scv2(i, j) = gs*gs

        corioq(i, j) = f
        coriop(i, j) = f
        betafp(i, j) = f/(tan(lat0/radian)*rearth)

        angle(i, j) = 0._r8
        cosang(i, j) = 1._r8
        sinang(i, j) = 0._r8

      enddo
    enddo
    !$omp end parallel do

  end subroutine geoenv_fuk95

  subroutine inifrc_fuk95

    integer :: i, j, l

    !$omp parallel do private(l, i)
    do j = 1, jj
      do l = 1, isp(j)
        do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
          surflx(i, j) = 0._r8
          sswflx(i, j) = 0._r8
          salflx(i, j) = 0._r8
          brnflx(i, j) = 0._r8
          surrlx(i, j) = 0._r8
          salrlx(i, j) = 0._r8
        enddo
      enddo
      do l = 1, isu(j)
        do i = max(1, ifu(j, l)), min(ii, ilu(j, l))
          taux(i, j) = 0._r8
        enddo
      enddo
      do l = 1, isv(j)
        do i = max(1, ifv(j, l)), min(ii, ilv(j, l))
          tauy(i, j) = 0._r8
        enddo
      enddo
    enddo
    !$omp end parallel do

  end subroutine inifrc_fuk95

  subroutine ictsz_fuk95

    real(r8), dimension(1 - nbdy:idm + nbdy, &
         1 - nbdy:jdm + nbdy, kdm + 1) :: z
    real(r8), dimension(kdm) :: sigref
    real(r8) :: drhojet, dsig, x, sigm, sigi, s0, s1, v0, v1, zu, zl
    integer :: i, j, k, l


    ! Set reference potential density, interface depths and layer salinity and
    ! temperature.

    select case (vcoord_tag)

      case (vcoord_isopyc_bulkml)

        ! For vertical coordinate featuring bulk surface mixed with
        ! isopycnic layers below, set layer reference potential densities
        ! and corresponding isopycnic layer structure. The bulk mixed layer
        ! is set to the minimum mixed layer thickness.

        drhojet = rhoc*f*u0*l0/(g*h1)
        dsig = (drho + drhojet)/(kk - 4)
        sigref(kk) = rhob - rho0
        sigref(kk - 1) = rhoc + .5_r8*(drho + drhojet) - rho0
        do k = kk - 2, 1, -1
          sigref(k) = sigref(k + 1) - dsig
        enddo

        !$omp parallel do private(k, l, i)
        do j = 1, jj
          do k = 1, kk
            do l = 1, isp(j)
              do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
                sigmar(i, j, k) = sigref(k)
                sigma(i, j, k) = sigref(k)
                saln(i, j, k) = saln0
                temp(i, j, k) = tofsig(sigma(i, j, k), saln(i, j, k))
              enddo
            enddo
          enddo
        enddo
        !$omp end parallel do

        !$omp parallel do private(k, l, i, x, sigm, sigi)
        do j = 1, jj
          do l = 1, isp(j)
            do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
              x = x_nudge(real(i, r8), real(j, r8))
              z(i, j, 1) = 0._r8
              z(i, j, 2) = .5_r8*mltmin*L_mks2cgs
              z(i, j, 3) = mltmin*L_mks2cgs
              z(i, j, kk    ) = h1
              z(i, j, kk + 1) = h0
              sigm = rhoc*(1._r8 + f*u0*x_psi(x)/(g*h1)) - rho0
              sigma(i, j, 1) = sigm &
                   + .5_r8*drho*(z(i, j, 2) + z(i, j, 1) - h1)/h1
              sigma(i, j, 2) = sigm &
                   + .5_r8*drho*(z(i, j, 3) + z(i, j, 2) - h1)/h1
              temp(i, j, 1) = tofsig(sigma(i, j, 1), saln(i, j, 1))
              temp(i, j, 2) = tofsig(sigma(i, j, 2), saln(i, j, 2))
            enddo
          enddo
          do k = 4, kk - 1
            do l = 1, isp(j)
              do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
                x = x_nudge(real(i, r8), real(j, r8))
                sigm = rhoc*(1._r8 + f*u0*x_psi(x)/(g*h1)) - rho0
                sigi = .5_r8*(sigref(k - 1) + sigref(k))
                z(i, j, k) = ((sigi - sigm)/drho + .5_r8)*h1
                z(i, j, k) = min(z(i, j, kk) - mindz*(kk - k), &
                     max(z(i, j, 3), z(i, j, k)))
              enddo
            enddo
          enddo
        enddo
        !$omp end parallel do

      case default

        ! For hybrid vertical coordinate featuring pressure coordinates
        ! towards the surface and continous isopycnal below, set layer
        ! interface reference potential densities. Initially the lowest
        ! model layer occupy everything below the active layer, while the
        ! active layer is distributed equally among the remaining model
        ! layers using constant z-level interfaces.

        !           drhojet = rhoc*f*u0*l0/(g*h1)
        !           dsig = (drho + drhojet)/(kk - 4)
        !           sigref(kk) = .5_r8*(rhob + rhoc) + .25_r8*(drho + drhojet) - rho0
        !           sigref(kk - 1) = rhoc + .5_r8*(drho + drhojet - dsig) - rho0
        !           do k = kk - 2, 1, -1
        !              sigref(k) = sigref(k + 1) - dsig
        !           enddo
        drhojet = rhoc*f*u0*l0/(g*h1)
        dsig = (drho + drhojet)/(kk - 5)
        sigref(kk - 2) = rhoc + .5_r8*(drho + drhojet - dsig) - rho0
        do k = kk - 3, 1, -1
          sigref(k) = sigref(k + 1) - dsig
        enddo
        sigref(kk    ) = rhob - rho0
        sigref(kk - 1) = (2._r8*sigref(kk - 2) + sigref(kk))/3._r8
        sigref(kk    ) = (sigref(kk - 2) + 2._r8*sigref(kk))/3._r8

        !$omp parallel do private(k, l, i)
        do j = 1, jj
          do k = 1, kk
            do l = 1, isp(j)
              do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
                sigmar(i, j, k) = sigref(k)
                saln(i, j, k) = saln0
                z(i, j, k) = real(k - 1, r8)*h0/real(kk, r8)
              enddo
            enddo
          enddo
          do l = 1, isp(j)
            do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
              z(i, j, kk + 1) = h0
            enddo
          enddo
        enddo
        !$omp end parallel do

        s0 = rhob - rho0
        !$omp parallel do private(k, l, i, x, s1)
        do j = 1, jj
          do k = 1, kk
            do l = 1, isp(j)
              do i = max(1, ifp(j, l)), min(ii, ilp(j, l))
                x = x_nudge(real(i, r8), real(j, r8))
                s1 = rhoc*(1._r8 + f*u0*x_psi(x)/(g*h1)) - rho0 &
                     + .5_r8*drho*(z(i, j, k + 1) + z(i, j, k) - h1)/h1
                sigma(i, j, k) = &
                     ( s1*max(0._r8, min(z(i, j, k + 1), h1) - z(i, j, k)) &
                     + s0*max(0._r8, z(i, j, k + 1) - max(z(i, j, k), h1))) &
                     /(z(i, j, k + 1) - z(i, j, k))
                temp(i, j, k) = tofsig(sigma(i, j, k), saln(i, j, k))
              enddo
            enddo
          enddo
        enddo
        !$omp end parallel do

    end select

    ! Set layer velocity.
    call xctilr(z, 1, kk + 1, 1, 1, halo_ps)
    v0 = 0._r8
    !$omp parallel do private(k, l, i, x, zu, zl, v1)
    do j = 1, jj
      do k = 1, kk - 1
        do l = 1, isv(j)
          do i = max(1, ifv(j, l)), min(ii, ilv(j, l))
            x = x_nudge(real(i, r8), real(j, r8) - .5_r8)
            zu = .5_r8*(z(i, j - 1, k    ) + z(i, j, k    ))
            zl = .5_r8*(z(i, j - 1, k + 1) + z(i, j, k + 1))
            v1 = u0*psi(x)*(h1 - .5*(zu + zl))/h1
            v1 = 0._r8
            if (abs(zl - zu) < epsilz) then
              v(i, j, k) = v1
            else
              v(i, j, k) = ( v1*max(0._r8, min(zl, h1) - zu) &
                   + v0*max(0._r8, zl - max(zu, h1)))/(zl - zu)
            endif
          enddo
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! Compute layer interface geopotential.
    !$omp parallel do private(k, l, i)
    do j = 1, jj
      do k = 1, kk + 1
        do l = 1, isp(j)
          do i = max(1, ifp(j, l)),min(ii, ilp(j, l))
            phi(i, j, k) = - g*z(i, j, k)
          enddo
        enddo
      enddo
    enddo
    !$omp end parallel do

  end subroutine ictsz_fuk95

end module mod_fuk95
