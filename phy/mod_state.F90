! ------------------------------------------------------------------------------
! Copyright (C) 2020-2025 Mats Bentsen
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

module mod_state
! ------------------------------------------------------------------------------
! This module contains declaration and initialization of state variables.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_constants, only: spval
   use mod_xc
   use mod_checksum, only: csdiag, chksum

   implicit none

   private

   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, 2*kdm) :: &
      u,         & ! u-component of baroclinic velocity [m s-1].
      v,         & ! v-component of baroclinic velocity [m s-1].
      dp,        & ! Layer pressure thickness [kg m-1 s-2].
      dpu,       & ! Layer pressure thickness at u-point [kg m-1 s-2].
      dpv,       & ! Layer pressure thickness at v-point [kg m-1 s-2].
      temp,      & ! Potential temperature [deg C].
      saln,      & ! Salinity [g kg-1].
      sigma,     & ! Potential density [kg m-3].
      uflx,      & ! u-component of mass flux [kg m s-2].
      vflx,      & ! v-component of mass flux [kg m s-2].
      utflx,     & ! u-component of heat flux [K kg m s-2].
      vtflx,     & ! v-component of heat flux [K kg m s-2].
      usflx,     & ! u-component of heat flux [g m s-2].
      vsflx        ! v-component of heat flux [g m s-2].

   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, kdm + 1) :: &
      p,         & ! Layer interface pressure [kg m-1 s-2].
      pu,        & ! Layer interface pressure at u-points [kg m-1 s-2].
      pv,        & ! Layer interface pressure at v-points [kg m-1 s-2].
      phi          ! Layer interface geopotential [m2 s-2].

   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, kdm) :: &
      cau,       & ! u-component of flux area [m2].
      cav          ! v-component of flux area [m2].

   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, 3) :: &
      ubflxs,    & ! u-component of barotropic mass flux sum [kg m s-3].
      vbflxs       ! v-component of barotropic mass flux sum [kg m s-3].

   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, 2) :: &
      ub,        & ! u-component of barotropic velocity [m s-1].
      vb,        & ! v-component of barotropic velocity [m s-1].
      pb,        & ! Bottom pressure [kg m-1 s-2].
      pbu,       & ! Bottom pressure at u-point [kg m-1 s-2].
      pbv,       & ! Bottom pressure at v-point [kg m-1 s-2].
      ubflxs_p,  & ! u-component of predicted barotropic mass flux sum
                   ! [kg m s-3].
      vbflxs_p     ! v-component of predicted barotropic mass flux sum
                   ! [kg m s-3].

   real(r8), dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy) :: &
      pb_p,      & ! Predicted bottom pressure [kg m-1 s-2].
      pbu_p,     & ! Predicted bottom pressure at u-point [kg m-1 s-2].
      pbv_p,     & ! Predicted bottom pressure at v-point [kg m-1 s-2].
      ubcors_p,  & ! u-component of predicted sum of barotropic coriolis term
                   ! [m s-2].
      vbcors_p,  & ! v-component of predicted sum of barotropic coriolis term
                   ! [m s-2].
      sealv        ! Sea surface height [m].

   integer, dimension(1 - nbdy:idm + nbdy, 1 - nbdy:jdm + nbdy, 2) :: &
      kfpla        ! Index of first mass containing layer below the mixed layer.

   public :: u, v, dp, dpu, dpv, temp, saln, sigma, &
             uflx, vflx, utflx, vtflx, usflx, vsflx, &
             p, pu, pv, phi, cau, cav, ubflxs, vbflxs, &
             ub, vb, pb, pbu, pbv, ubflxs_p, vbflxs_p, &
             pb_p, pbu_p, pbv_p, ubcors_p, vbcors_p, sealv, kfpla, &
             inivar_state, init_fluxes

contains

   subroutine inivar_state
   ! ---------------------------------------------------------------------------
   ! Initialize state variables.
   ! ---------------------------------------------------------------------------

      integer :: i, j, k, l

      pb_p(:,:) = spval
      pbu_p(:,:) = spval
      pbv_p(:,:) = spval
      ubcors_p(:,:) = spval
      vbcors_p(:,:) = spval
      sealv(:,:) = spval
      pb(:,:,:) = spval
      ub(:,:,:) = spval
      vb(:,:,:) = spval
      ubflxs_p(:,:,:) = spval
      vbflxs_p(:,:,:) = spval
      pbu(:,:,:) = spval
      pbv(:,:,:) = spval
      ubflxs(:,:,:) = spval
      vbflxs(:,:,:) = spval
      p(:,:,:) = spval
      pu(:,:,:) = spval
      pv(:,:,:) = spval
      phi(:,:,:) = spval
      u(:,:,:) = spval
      v(:,:,:) = spval
      uflx(:,:,:) = spval
      utflx(:,:,:) = spval
      usflx(:,:,:) = spval
      vflx(:,:,:) = spval
      vtflx(:,:,:) = spval
      vsflx(:,:,:) = spval
      dp(:,:,:) = spval
      dpu(:,:,:) = spval
      dpv(:,:,:) = spval
      temp(:,:,:) = spval
      saln(:,:,:) = spval
      sigma(:,:,:) = spval

      cau(:,:,:) = 0._r8
      cav(:,:,:) = 0._r8

   !$omp parallel do private(l, i, k)
      do j = 1, jj + 1
         do l = 1, isq(j)
         do i = max(1, ifq(j, l)), min(ii + 1, ilq(j, l))
            do k = 1, 2
               pb(i    , j    , k) = 0._r8
               pb(i - 1, j    , k) = 0._r8
               pb(i    , j - 1, k) = 0._r8
               pb(i - 1, j - 1, k) = 0._r8
            enddo
            pb_p(i    , j    ) = 0._r8
            pb_p(i - 1, j    ) = 0._r8
            pb_p(i    , j - 1) = 0._r8
            pb_p(i - 1, j - 1) = 0._r8
            p(i    , j    , 1) = 0._r8
            p(i - 1, j    , 1) = 0._r8
            p(i    , j - 1, 1) = 0._r8
            p(i - 1, j - 1, 1) = 0._r8
            do k = 1, 2*kk
               dp(i    , j    , k) = 0._r8
               dp(i - 1, j    , k) = 0._r8
               dp(i    , j - 1, k) = 0._r8
               dp(i - 1, j - 1, k) = 0._r8
            enddo
         enddo
         enddo
      enddo
   !$omp end parallel do

      call xctilr(pb,   1,    2, nbdy, nbdy, halo_ps)
      call xctilr(pb_p, 1,    1, nbdy, nbdy, halo_ps)
      call xctilr(p,    1,    1, nbdy, nbdy, halo_ps)
      call xctilr(dp,   1, 2*kk, nbdy, nbdy, halo_ps)

   !$omp parallel do private(l, i)
      do j = 1 - nbdy, jj + nbdy
         do l = 1, isp(j)
         do i = ifp(j, l), ilp(j, l)
            p(i, j, 1) = 0._r8
         enddo
         enddo
      enddo
   !$omp end parallel do

   ! Initialize 'u', 'ub', 'uflx' at points located upstream and downstream (in
   ! i-direction) of p-points. Initialize 'pbu', 'dpu' upstream and downstream
   ! of p-points as well as at lateral neighbors of interior u-points.
   !$omp parallel do private(l, i, k)
      do j = 0, jj + 1
         do l = 1, isu(j)
         do i = max(1, ifu(j, l)), min(ii, ilu(j, l))
            pu(i, j, 1) = 0._r8
            do k = 1, 2
               pbu(i, j - 1, k) = 0._r8
               pbu(i, j + 1, k) = 0._r8
            enddo
            pbu_p(i, j - 1) = 0._r8
            pbu_p(i, j + 1) = 0._r8
            do k = 1, 2*kk
               dpu(i, j - 1, k) = 0._r8
               dpu(i, j + 1, k) = 0._r8
            enddo
         enddo
         enddo
      enddo
   !$omp end parallel do
   !$omp parallel do private(l, i, k)
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j, l)), min(ii, ilp(j, l) + 1)
            do k = 1, 3
               ubflxs(i, j, k) = 0._r8
            enddo
            do k = 1, 2
               ub(i, j, k) = 0._r8
               ubflxs_p(i, j, k) = 0._r8
               pbu(i, j, k) = 0._r8
            enddo
            pbu_p(i, j) = 0._r8
            ubcors_p(i, j) = 0._r8
            do k = 1, 2*kk
               dpu  (i, j, k) = 0._r8
               uflx (i, j, k) = 0._r8
               utflx(i, j, k) = 0._r8
               usflx(i, j, k) = 0._r8
               u    (i, j, k) = 0._r8
            enddo
         enddo
         enddo
      enddo
   !$omp end parallel do

      call xctilr(pu,       1,    1, nbdy, nbdy, halo_us)
      call xctilr(pbu,      1,    2, nbdy, nbdy, halo_us)
      call xctilr(pbu_p,    1,    1, nbdy, nbdy, halo_us)
      call xctilr(dpu,      1, 2*kk, nbdy, nbdy, halo_us)
      call xctilr(ub,       1,    2, nbdy, nbdy, halo_us)
      call xctilr(ubflxs,   1,    3, nbdy, nbdy, halo_us)
      call xctilr(ubflxs_p, 1,    2, nbdy, nbdy, halo_us)
      call xctilr(ubcors_p, 1,    1, nbdy, nbdy, halo_us)
      call xctilr(uflx,     1, 2*kk, nbdy, nbdy, halo_us)
      call xctilr(utflx,    1, 2*kk, nbdy, nbdy, halo_us)
      call xctilr(usflx,    1, 2*kk, nbdy, nbdy, halo_us)
      call xctilr(u,        1, 2*kk, nbdy, nbdy, halo_us)

   ! Initialize 'v', 'vb', 'vflx' at points located upstream and downstream (in
   ! j-direction) of p-points. Initialize 'pbv', 'dpv' upstream and downstream
   ! of p-points as well as at lateral neighbors of interior v-points.
   !$omp parallel do private(l, j, k)
      do i = 0, ii + 1
         do l = 1, jsv(i)
         do j = max(1, jfv(i, l)), min(jj, jlv(i, l))
            pv(i, j, 1) = 0._r8
            do k = 1, 2
               pbv(i - 1, j, k) = 0._r8
               pbv(i + 1, j, k) = 0._r8
            enddo
            pbv_p(i - 1, j) = 0._r8
            pbv_p(i + 1, j) = 0._r8
            do k = 1, 2*kk
               dpv(i - 1, j, k) = 0._r8
               dpv(i + 1, j, k) = 0._r8
            enddo
         enddo
         enddo
      enddo
   !$omp end parallel do
   !$omp parallel do private(l, j, k)
      do i = 1, ii
         do l = 1, jsp(i)
         do j = max(1, jfp(i, l)), min(jj, jlp(i, l) + 1)
            do k = 1, 3
               vbflxs(i, j, k) = 0._r8
            enddo
            do k = 1, 2
               vb(i, j, k) = 0._r8
               vbflxs_p(i, j, k) = 0._r8
               pbv(i, j, k) = 0._r8
            enddo
            pbv_p(i, j) = 0._r8
            vbcors_p(i, j) = 0._r8
            do k = 1, 2*kk
               dpv  (i, j, k) = 0._r8
               vflx (i, j, k) = 0._r8
               vtflx(i, j, k) = 0._r8
               vsflx(i, j, k) = 0._r8
               v    (i, j, k) = 0._r8
            enddo
         enddo
         enddo
      enddo
   !$omp end parallel do

      call xctilr(pv,       1,    1, nbdy, nbdy, halo_vs)
      call xctilr(pbv,      1,    2, nbdy, nbdy, halo_vs)
      call xctilr(pbv_p,    1,    1, nbdy, nbdy, halo_vs)
      call xctilr(dpv,      1, 2*kk, nbdy, nbdy, halo_vs)
      call xctilr(vb,       1,    2, nbdy, nbdy, halo_vs)
      call xctilr(vbflxs,   1,    3, nbdy, nbdy, halo_vs)
      call xctilr(vbflxs_p, 1,    2, nbdy, nbdy, halo_vs)
      call xctilr(vbcors_p, 1,    1, nbdy, nbdy, halo_vs)
      call xctilr(vflx,     1, 2*kk, nbdy, nbdy, halo_vs)
      call xctilr(vtflx,    1, 2*kk, nbdy, nbdy, halo_vs)
      call xctilr(vsflx,    1, 2*kk, nbdy, nbdy, halo_vs)
      call xctilr(v,        1, 2*kk, nbdy, nbdy, halo_vs)

      if (csdiag) then
         if (mnproc == 1) then
            write (lp, *) 'inivar_state:'
         endif
         call chksum(p       , kk + 1, halo_ps, 'p'       )
         call chksum(pu      , kk + 1, halo_us, 'pu'      )
         call chksum(pv      , kk + 1, halo_vs, 'pv'      )
         call chksum(pb      , 2     , halo_ps, 'pb'      )
         call chksum(pb_p    , 1     , halo_ps, 'pb_p'    )
         call chksum(ub      , 2     , halo_uv, 'ub'      )
         call chksum(vb      , 2     , halo_vv, 'vb'      )
         call chksum(pbu     , 2     , halo_us, 'pbu'     )
         call chksum(pbu_p   , 1     , halo_us, 'pbu_p'   )
         call chksum(pbv     , 2     , halo_vs, 'pbv'     )
         call chksum(pbv_p   , 1     , halo_vs, 'pbv_p'   )
         call chksum(ubflxs  , 3     , halo_uv, 'ubflxs'  )
         call chksum(ubflxs_p, 2     , halo_uv, 'ubflxs_p')
         call chksum(vbflxs  , 3     , halo_vv, 'vbflxs'  )
         call chksum(vbflxs_p, 2     , halo_vv, 'vbflxs_p')
         call chksum(ubcors_p, 1     , halo_uv, 'ubcors_p')
         call chksum(vbcors_p, 1     , halo_vv, 'vbcors_p')
         call chksum(u       , 2*kk  , halo_uv, 'u'       )
         call chksum(v       , 2*kk  , halo_vv, 'v'       )
         call chksum(cau     , kk    , halo_uv, 'cau'     )
         call chksum(cav     , kk    , halo_vv, 'cav'     )
         call chksum(uflx    , 2*kk  , halo_uv, 'uflx'    )
         call chksum(vflx    , 2*kk  , halo_vv, 'vflx'    )
         call chksum(dp      , 2*kk  , halo_ps, 'dp'      )
         call chksum(dpu     , 2*kk  , halo_us, 'dpu'     )
         call chksum(dpv     , 2*kk  , halo_vs, 'dpv'     )
      endif

   end subroutine inivar_state

   subroutine init_fluxes(m, n, mm, nn, k1m, k1n, update_flux_halos)
   ! ---------------------------------------------------------------------------
   ! Reset fluxes to be accumulated over a model time step and update flux
   ! halos if requested.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n
      logical, intent(in) :: update_flux_halos

      integer :: i, j, k, km, l

    !$omp parallel do private(k, km, l, i)
        do j = 0, jj+2
           do k = 1, kk
              km = k + mm
              do l=1, isu(j)
              do i=max(0, ifu(j,l)), min(ii+2, ilu(j,l))
                 uflx(i,j,km) = 0._r8
                 utflx(i,j,km) = 0._r8
                 usflx(i,j,km) = 0._r8
              enddo
              enddo
              do l=1, isv(j)
              do i=max(0, ifv(j,l)), min(ii+2, ilv(j,l))
                 vflx(i,j,km) = 0._r8
                 vtflx(i,j,km) = 0._r8
                 vsflx(i,j,km) = 0._r8
              enddo
              enddo
           enddo
        enddo
    !$omp end parallel do

      if (update_flux_halos) then
         call xctilr(uflx (1-nbdy,1-nbdy,k1n), 1, kk, 1, 1, halo_uv)
         call xctilr(utflx(1-nbdy,1-nbdy,k1n), 1, kk, 1, 1, halo_uv)
         call xctilr(usflx(1-nbdy,1-nbdy,k1n), 1, kk, 1, 1, halo_uv)
         call xctilr(vflx (1-nbdy,1-nbdy,k1n), 1, kk, 1, 1, halo_vv)
         call xctilr(vtflx(1-nbdy,1-nbdy,k1n), 1, kk, 1, 1, halo_vv)
         call xctilr(vsflx(1-nbdy,1-nbdy,k1n), 1, kk, 1, 1, halo_vv)
      endif

   end subroutine init_fluxes

end module mod_state
