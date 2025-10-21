! ------------------------------------------------------------------------------
! Copyright (C) 2021-2025 Mats Bentsen, Mehmet Ilicak
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

module mod_ale_regrid_remap
! ------------------------------------------------------------------------------
! This module contains parameter, variables and procedures related to the
! regridding and remapping needed by the ALE method.
! ------------------------------------------------------------------------------

   use mod_types,     only: r8
   use mod_config,    only: inst_suffix
   use mod_constants, only: grav, epsilp, onem
   use mod_time,      only: delt1
   use mod_xc
   use mod_grid,      only: scuy, scvx, scp2, scuxi, scvyi, scp2i
   use mod_vcoord,    only: vcoord_tag, vcoord_isopyc_bulkml, &
                            vcoord_cntiso_hybrid, vcoord_plevel, &
                            sigmar, plevel
   use mod_eos,       only: sig, dsigdt, dsigds
   use mod_state,     only: u, v, dp, dpu, dpv, temp, saln, sigma, p, pu, pv, &
                            ub, vb, utflx, vtflx, usflx, vsflx
   use mod_hor3map,   only: recon_grd_struct, recon_src_struct, remap_struct, &
                            hor3map_plm, hor3map_ppm, hor3map_pqm, &
                            hor3map_monotonic, hor3map_non_oscillatory, &
                            hor3map_non_oscillatory_posdef, &
                            initialize_rcgs, initialize_rcss, initialize_rms, &
                            prepare_reconstruction, reconstruct, &
                            extract_polycoeff, regrid, &
                            prepare_remapping, remap, &
                            hor3map_noerr, hor3map_errstr
   use mod_diffusion, only: ltedtp_opt, ltedtp_neutral, difmxp, &
                            utflld, vtflld, usflld, vsflld
   use mod_ndiff,     only: ndiff_prep_jslice, ndiff_uflx_jslice, &
                            ndiff_vflx_jslice, ndiff_update_trc_jslice
   use mod_dia,       only: ddm, nphy, alarm_phy, &
                            depthslev_bnds, pbath, ubath, vbath, phylvl, &
                            acc_templvl, acc_salnlvl, &
                            acc_uvellvl, acc_vvellvl, &
                            acc_idlagelvl
   use mod_checksum,  only: csdiag, chksum
   use mod_tracers,   only: ntr, itriag, trc
   use mod_ifdefs,    only: use_TRC, use_IDLAGE
   use mod_utility,   only: util1

   implicit none
   private

   ! Options with default values, modifiable by namelist.
   character(len = 80) :: &
        reconstruction_method = 'ppm', &
        density_limiting = 'monotonic', &
        tracer_limiting = 'non_oscillatory', &
        velocity_limiting = 'non_oscillatory', &
        regrid_method = 'nudge'
   logical :: &
        density_pc_upper_bndr = .false., &
        density_pc_lower_bndr = .false., &
        tracer_pc_upper_bndr = .true., &
        tracer_pc_lower_bndr = .false., &
        velocity_pc_upper_bndr = .true., &
        velocity_pc_lower_bndr = .false.
   real(r8) :: &
        dpmin_interior         = .1_r8, &
        regrid_nudge_ts        = 86400._r8, &
        stab_fac_limit         = .75_r8, &
        dpvar_fac              = .75_r8, &
        smooth_diff_max        = 50000._r8
   integer :: &
        upper_bndr_ord = 6, &
        lower_bndr_ord = 4, &
        k_range_plevel = 1, &
        dktzu = 4, &
        dktzl = 2

   ! Options derived from string options.
   integer :: &
        reconstruction_method_tag, &
        density_limiting_tag, &
        tracer_limiting_tag, &
        velocity_limiting_tag, &
        regrid_method_tag

   real(r8), parameter :: &
        bfsq_min      = 1.e-7_r8, &   ! Minimum buoyancy frequency squared in
                                      ! monotonized potential density to be used
                                      ! in regridding [s-2].
        regrid_mval   = - 1.e33_r8, & ! Missing value for regridding.
        x_eps         = 1.e-14_r8     ! Small non-dimensional value used in the
                                      ! construction of Bezier curves.
   integer, parameter :: &
        regrid_method_direct = 1, &   ! Regrid method (vcoord_tag ==
                                      ! vcoord_cntiso_hybrid): On the basis of
                                      ! reconstructed potential density, regrid
                                      ! interface pressures so interface
                                      ! potential densities match target values.
        regrid_method_nudge  = 2      ! Regrid method (vcoord_tag ==
                                      ! vcoord_cntiso_hybrid): Nudge interface
                                      ! pressures to reduce the deviation from
                                      ! the interface reference potential
                                      ! density.

   integer :: ntr_loc

   type(recon_grd_struct) :: rcgs
   type(recon_src_struct) :: d_rcss, v_rcss
   type(recon_src_struct), allocatable, dimension(:) :: trc_rcss
   type(remap_struct) :: rms, rms_diazlv

   public :: regrid_method_tag, regrid_method_direct, &
             readnml_ale_regrid_remap, init_ale_regrid_remap, &
             ale_regrid_remap, ale_remap_diazlv

contains

   ! ---------------------------------------------------------------------------
   ! Private procedures.
   ! ---------------------------------------------------------------------------

   pure function peval0(pc) result(f)

      real(r8), dimension(:), intent(in) :: pc

      real(r8) :: f

      f = pc(1)

   end function peval0

   pure function peval1(pc) result(f)

      real(r8), dimension(:), intent(in) :: pc

      real(r8) :: f

      f = pc(1) + pc(2) + pc(3) + pc(4) + pc(5)

   end function peval1

   pure function dpeval0(pc) result(f)

      real(r8), dimension(:), intent(in) :: pc

      real(r8) :: f

      f = pc(2)

   end function dpeval0

   pure function dpeval1(pc) result(f)

      real(r8), dimension(:), intent(in) :: pc

      real(r8) :: f

      real(r8), parameter :: &
           c2 = 2._r8, &
           c3 = 3._r8, &
           c4 = 4._r8

      f = pc(2) + c2*pc(3) + c3*pc(4) + c4*pc(5)

   end function dpeval1

   subroutine reconstruct_trc_jslice(p_src, ksmx, tpc_src, t_srcdi, &
                                     ilb, iub, j, js, nn)
   ! ---------------------------------------------------------------------------
   ! Vertically reconstruct temperature, salinity and additional tracers along a
   ! j-slice of the model data.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:,1-nbdy:), intent(out) :: p_src
      integer, dimension(1-nbdy:), intent(out) :: ksmx
      real(r8), intent(out), dimension(:,:,:,1-nbdy:) :: tpc_src, t_srcdi
      integer, intent(in) :: ilb, iub, j, js, nn

      real(r8), dimension(kdm,ntr_loc) :: trc_1d
      integer :: l, i, k, kn, nt, errstat


      do l = 1, isp(j)
         do i = max(ilb, ifp(j,l)), min(iub, ilp(j,l))

            ! Compute source layer interface pressure and copy variables into 1D
            ! arrays.
            p_src(1,i) = p(i,j,1)
            do k = 1, kk
               kn = k + nn
               p_src(k+1,i) = p_src(k,i) + dp(i,j,k+nn)
               trc_1d(k,1) = temp(i,j,kn)
               trc_1d(k,2) = saln(i,j,kn)
               do nt = 1, ntr
                  trc_1d(k,nt+2) = trc(i,j,kn,nt)
               enddo
            enddo

            ! Find index of deepest source layer with non-zero thickness.
            ksmx(i) = kk
            do k = kk, 1, -1
               if (p_src(k,i) == p_src(kk+1,i)) ksmx(i) = k - 1
            enddo

            errstat = prepare_reconstruction(rcgs, p_src(:,i), i, js)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(prep_recon_jslice)')
               stop '(prep_recon_jslice)'
            endif

            ! Reconstruct tracers.
            do nt = 1, ntr_loc
               errstat = reconstruct(rcgs, trc_rcss(nt), trc_1d(:,nt), i, js)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(recon_trc_jslice)')
                  stop '(recon_trc_jslice)'
               endif
            enddo

            ! Extract polynomial coefficients of the reconstructions and store
            ! tracer variables in dual interface arrays with with values
            ! corresponding to upper and lower interface of each layer.
            do nt = 1, ntr_loc
               errstat = extract_polycoeff(trc_rcss(nt), tpc_src(:,:,nt,i), &
                                           i, js)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(recon_trc_jslice)')
                  stop '(recon_trc_jslice)'
               endif
               do k = 1, ksmx(i)
                  t_srcdi(1,k,nt,i) = peval0(tpc_src(:,k,nt,i))
                  t_srcdi(2,k,nt,i) = peval1(tpc_src(:,k,nt,i))
               enddo
            enddo

         enddo
      enddo

   end subroutine reconstruct_trc_jslice

   subroutine regrid_plevel_jslice(p_src, p_dst, ilb, iub, j)
   ! ---------------------------------------------------------------------------
   ! For vcoord == 'plevel', regrid interface pressures to specified pressure
   ! levels.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:,1-nbdy:), intent(in) :: p_src
      real(r8), dimension(:,1-nbdy:), intent(out) :: p_dst
      integer, intent(in) :: ilb, iub, j

      integer :: l, i, k

      do l = 1, isp(j)
         do i = max(ilb, ifp(j,l)), min(iub, ilp(j,l))
            do k = 1, kk
               p_dst(k,i) = min(plevel(k) + p_src(1,i), p_src(kk+1,i))
            enddo
            p_dst(kk+1,i) = p_src(kk+1,i)
         enddo
      enddo

   end subroutine regrid_plevel_jslice

   subroutine regrid_cntiso_hybrid_direct_jslice(p_src, p_dst, &
                                                 ilb, iub, j, js, nn)
   ! ---------------------------------------------------------------------------
   ! For vcoord == 'cntiso_hybrid' and regrid_method = 'direct', regrid
   ! interface pressures so interface potential densities match target values
   ! except where minimum layer thickness towards the surface must be
   ! maintained.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:,1-nbdy:), intent(in) :: p_src
      real(r8), dimension(:,1-nbdy:), intent(out) :: p_dst
      integer, intent(in) :: ilb, iub, j, js, nn

      real(r8), dimension(kdm+1) :: sig_trg
      real(r8), dimension(kdm) :: sig_src
      real(r8) :: beta, sdpsum, smean, dpmin, pku, pku_test, pmin, dpt, &
                  pt, ptu1, ptl1, ptu2, ptl2, w1, x
      integer :: l, i, k, kn, ks, ke, kl, ku, errstat
      logical :: thin_layers, layer_added

      ! Minimum potential density difference with respect to pressure for
      ! potential density to be used in regridding.
      beta = bfsq_min/(grav*grav)

      do l = 1, isp(j)
         do i = max(ilb, ifp(j,l)), min(iub, ilp(j,l))

            ! Copy source and target potential densities into 1D arrays.
            do k = 1, kk
               kn = k + nn
               sig_src(k) = sigma(i,j,kn)
               sig_trg(k) = sigmar(i,j,k)
            enddo
            sig_trg(kk+1) = sig_trg(kk)

            ! Make sure potential density to be used in regridding is
            ! monotonically increasing with depth.
            kl = kk
            ku = kl - 1
            do while (ku > 0)
               thin_layers = p_src(kl+1,i) - p_src(ku,i) < epsilp
               if (thin_layers .or. &
                   sig_src(kl) - sig_src(ku) &
                   < .5_r8*beta*(p_src(kl+1,i) - p_src(ku,i))) then
                  sdpsum = sig_src(ku)*(p_src(ku+1,i) - p_src(ku,i)) &
                         + sig_src(kl)*(p_src(kl+1,i) - p_src(kl,i))
                  if (.not. thin_layers) &
                     smean = sdpsum/(p_src(kl+1,i) - p_src(ku,i))
                  do
                     layer_added = .false.
                     if (ku > 1) then
                        if (thin_layers) then
                           ku = ku - 1
                           sdpsum = sdpsum &
                                  + sig_src(ku)*(p_src(ku+1,i) - p_src(ku,i))
                           thin_layers = p_src(kl+1,i) - p_src(ku,i) < epsilp
                           if (.not. thin_layers) &
                              smean = sdpsum/(p_src(kl+1,i) - p_src(ku,i))
                           layer_added = .true.
                        else
                           if (smean - sig_src(ku-1) &
                              < .5_r8*beta*(p_src(kl+1,i) - p_src(ku-1,i))) then
                              ku = ku - 1
                              sdpsum = sdpsum &
                                     + sig_src(ku)*(p_src(ku+1,i) - p_src(ku,i))
                              smean = sdpsum/(p_src(kl+1,i) - p_src(ku,i))
                              layer_added = .true.
                           endif
                        endif
                     endif
                     if (kl < kk) then
                        if (thin_layers) then
                           kl = kl + 1
                           sdpsum = sdpsum &
                                  + sig_src(kl)*(p_src(kl+1,i) - p_src(kl,i))
                           thin_layers = p_src(kl+1,i) - p_src(ku,i) < epsilp
                           if (.not. thin_layers) &
                              smean = sdpsum/(p_src(kl+1,i) - p_src(ku,i))
                           layer_added = .true.
                        else
                           if (sig_src(kl+1) - smean &
                              < .5_r8*beta*(p_src(kl+2,i) - p_src(ku,i))) then
                              kl = kl + 1
                              sdpsum = sdpsum &
                                     + sig_src(kl)*(p_src(kl+1,i) - p_src(kl,i))
                              smean = sdpsum/(p_src(kl+1,i) - p_src(ku,i))
                              layer_added = .true.
                           endif
                        endif
                     endif
                     if (.not. layer_added) exit
                  enddo
                  do k = ku, kl
                     sig_src(k) = smean &
                                + .5_r8*beta*( p_src(k ,i) + p_src(k +1,i) &
                                             - p_src(ku,i) - p_src(kl+1,i))
                  enddo
               endif
               kl = ku
               ku = kl - 1
            enddo

            ! Monotonically reconstruct potential density.
            errstat = reconstruct(rcgs, d_rcss, sig_src, i, js)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(cntiso_regrid_direct_jslice)')
               stop '(cntiso_regrid_direct_jslice)'
            endif

            ! On the basis of the reconstructed potential density, regrid
            ! interface pressures so interface potential densities match target
            ! values.
            errstat = regrid(d_rcss, sig_trg, p_dst(:,i), regrid_mval, i, js)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(cntiso_regrid_direct_jslice)')
               stop '(cntiso_regrid_direct_jslice)'
            endif

            ! Modify regridded interface pressures to ensure the water column is
            ! properly bounded.
            k = 1
            do
               ks = k
               if (p_dst(k,i) /= regrid_mval) exit
               p_dst(k,i) = p_src(1,i)
               if (k > kk) exit
               k = k + 1
            enddo
            k = kk + 1
            do
               ke = k
               if (p_dst(k,i) /= regrid_mval) exit
               p_dst(k,i) = p_src(kk+1,i)
               if (k == 1) exit
               k = k - 1
            enddo
            p_dst(1,i) = p_src(1,i)
            p_dst(kk+1,i) = p_src(kk+1,i)

            ! If no regrid interface is found in the water column, try to place
            ! all water in the layer with potential density bounds that include
            ! the column mean potential density.
            if (ks == ke) then
               sdpsum = 0._r8
               do k = 1, kk
                  sdpsum = sdpsum + sig_src(k)*(p_src(k+1,i) - p_src(k,i))
               enddo
               smean = sdpsum/(p_src(kk+1,i) - p_src(1,i))
               ks = 2
               do while (ks <= kk)
                  if (smean < sig_trg(ks)) exit
                  ks = ks + 1
               enddo
               do k = ks, kk
                  p_dst(k,i) = p_src(kk+1,i)
               enddo
               ke = ks - 1
            endif

            ! Modify interface pressures so that layer thicknesses are
            ! above a specified threshold.
            dpmin = min(plevel(2) - plevel(1), dpmin_interior)
            ks = max(2, ks)
            ke = min(kk, ke)
            k = ks
            do while (k <= ke)
               if (p_dst(k+1,i) - p_dst(k,i) < dpmin) then
                  if (k == ke) then
                     p_dst(k,i) = p_dst(ke+1,i)
                  else
                     ku = k
                     kl = k + 1
                     pku = .5_r8*(p_dst(kl,i) + p_dst(ku,i) - dpmin)
                     do
                        layer_added = .false.
                        kl = kl + 1
                        pku_test = ((pku - dpmin)*(kl - ku) + p_dst(kl,i)) &
                                   /(kl - ku + 1)
                        if (pku_test + (kl - ku)*dpmin > p_dst(kl,i)) then
                           if (kl == ke + 1) exit
                           pku = pku_test
                           layer_added = .true.
                        else
                           kl = kl - 1
                        endif
                        ku = ku - 1
                        pku_test = ((pku - dpmin)*(kl - ku) + p_dst(ku,i)) &
                                   /(kl - ku + 1)
                        if (pku_test < p_dst(ku,i)) then
                           if (ku == 1) exit
                           pku = pku_test
                           layer_added = .true.
                        else
                           ku = ku + 1
                        endif
                        if (.not. layer_added) exit
                     enddo
                     if     (ku == 1) then
                        do k = 2, kl
                           p_dst(k,i) = min(p_dst(ke+1,i), &
                                            p_dst(k -1,i) + dpmin)
                        enddo
                        do k = kl+1, ke
                           p_dst(k,i) = &
                              min(p_dst(ke+1,i), &
                              max(p_dst(k,i), p_dst(1,i) + dpmin*(k - 1)))
                        enddo
                     elseif (kl == ke + 1) then
                        do k = ku, kl
                           p_dst(k,i) = p_dst(ke+1,i)
                        enddo
                     else
                        p_dst(ku,i) = pku
                        do k = ku+1, kl
                           p_dst(k,i) = p_dst(k-1,i) + dpmin
                        enddo
                     endif
                     k = kl
                  endif
               endif
               k = k + 1
            enddo

            ! Modify regridded interface pressures to ensure that a minimum
            ! layer thickness towards the surface is maintained. A smooth
            ! transition between modified and unmodified interfaces is sought.
            do k = 2, k_range_plevel
               p_dst(k,i) = min(p_dst(kk+1,i), plevel(k) + p_src(1,i))
            enddo
            dpt = plevel(k_range_plevel+1) - plevel(k_range_plevel)
            do k = k_range_plevel + 1, ke
               pmin = plevel(k) + p_src(1,i)
               dpt = max(p_dst(k+1,i) - p_dst(k,i), dpt, &
                         plevel(min(k,kk-1)+1) - plevel(min(k,kk-1)))
               pt = max(p_dst(k,i), pmin)
               ptu1 = pmin - dpt
               ptl1 = pmin + dpt
               ptu2 = pmin
               ptl2 = pmin + 2._r8*dpt
               w1 = min(1._r8,(p_dst(k,i) - p_src(1,i))/(pmin - p_src(1,i)))
               if (p_dst(k,i) > ptu1 .and. p_dst(k,i) < ptl1) then
                  x = .5_r8*(p_dst(k,i) - ptu1)/dpt
                  pt = pmin + dpt*x*x
               endif
               if (p_dst(k+1,i) > ptu2 .and. p_dst(k+1,i) < ptl2) then
                  x = .5_r8*(p_dst(k+1,i) - ptu2)/dpt
                  pt = w1*pt + (1._r8 - w1)*(pmin + dpt*x*x)
               endif
               p_dst(k,i) = min(p_dst(ke+1,i), max(p_dst(k-1,i) + dpmin, pt))
            enddo

         enddo
      enddo

   end subroutine regrid_cntiso_hybrid_direct_jslice

   subroutine regrid_cntiso_hybrid_nudge_jslice( &
                 p_src, ksmx, tpc_src, t_srcdi, p_dst, smooth_fac, ilb, iub, j)
   ! ---------------------------------------------------------------------------
   ! For vcoord == 'cntiso_hybrid' and regrid_method = 'nudge', nudge the
   ! interface pressures to reduce the deviation from interface target potential
   ! densities, while maintaining minimum layer thicknesses towards the surface.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:,1-nbdy:), intent(in) :: p_src
      integer, dimension(1-nbdy:), intent(in) :: ksmx
      real(r8), dimension(:,:,:,1-nbdy:), intent(in) :: tpc_src, t_srcdi
      real(r8), dimension(:,1-nbdy:), intent(out) :: p_dst, smooth_fac
      integer, intent(in) :: ilb, iub, j

      integer, parameter :: &
           it = 1, &
           is = 2

      real(r8), dimension(2,kdm) :: sig_srcdi
      integer, dimension(1-nbdy:idm+nbdy) :: kdmx

      real(r8), dimension(kdm+1) :: sig_trg, sig_pmin
      real(r8), dimension(kdm) :: dsig_trg, pmin, dpmin
      real(r8) :: sig_max, ckt, sig_up, sig_lo, dk, dki, &
                  dsigdx_up, dsigdx_lo, x, xi, si, t, nudge_fac, &
                  dsig, dsigdx, stab_fac, dp_up, dp_lo, sig_intrp, &
                  dpmin_sum, dp_sum, dpmin_sum_test, dp_sum_test
      integer :: l, i, k, kt, kl, ktzmin, ktzmax, ks, ke, ku
      logical :: layer_added

      do l = 1, isp(j)
         do i = max(ilb, ifp(j,l)), min(iub, ilp(j,l))

            ! Store density in a dual interface array with with values
            ! corresponding to upper and lower interface of each layer.  Also
            ! find the maximum lower interface potential density of the
            ! reconstructed column.
            sig_max = 0._r8
            do k = 1, ksmx(i)
               sig_srcdi(1,k) = sig(t_srcdi(1,k,it,i), t_srcdi(1,k,is,i))
               sig_srcdi(2,k) = sig(t_srcdi(2,k,it,i), t_srcdi(2,k,is,i))
               sig_max = max(sig_max, sig_srcdi(2,k))
            enddo

            ! Copy variables into 1D arrays.
            do k = 1, kk
               sig_trg(k) = sigmar(i,j,k)
            enddo
            sig_trg(kk+1) = sig_trg(kk)
            do k = 1, kk-1
               dsig_trg(k) = sig_trg(k+1) - sigmar(i,j,k)
            enddo
            dsig_trg(kk) = dsig_trg(kk-1)

            ! Find the index of the first layer which lower interface reference
            ! potential density is denser than the maximum lower interface
            ! potential density of the reconstructed column.
            do k = kk, 1, -1
               if (sig_trg(k) < sig_max) exit
            enddo
            kdmx(i) = max(1, k)

            ! Set minimum interface pressure.
            do k = 1, kk
               pmin(k) = min(plevel(k) + p_src(1,i), p_src(kk+1,i))
            enddo

            ! Set non-dimensional nudging factor.
            nudge_fac = delt1/regrid_nudge_ts

            ! Enforce or nudge towards minimum interface pressure for layer
            ! interface indices 1 to k_range_plevel.
            kl = 1
            sig_pmin(1) = sig_srcdi(1,1)
            p_dst(1,i) = pmin(1)
            smooth_fac(1,i) = 1._r8
            do k = 2, k_range_plevel
               do while (p_src(kl+1,i) < pmin(k))
                  kl = kl + 1
               enddo
               sig_pmin(k) = ( (p_src(kl+1,i) - pmin(k))*sig_srcdi(1,kl) &
                             + (pmin(k) - p_src(kl,i))*sig_srcdi(2,kl)) &
                             /(p_src(kl+1,i) - p_src(kl,i))
               p_dst(k,i) = p_src(k,i) + nudge_fac*(pmin(k) - p_src(k,i))
               p_dst(k,i) = min(max(p_dst(k,i), pmin(k), &
                                    p_dst(k-1,i) + dpmin_interior), &
                                p_src(kk+1,i))
               smooth_fac(k,i) = 1._r8
            enddo

            ! Find the index of the first interface with potential density at
            ! minimum interface pressure smaller than the reference potential
            ! density of this transition interface. A layer range above and
            ! below the transition interface may be specified, making a
            ! transition zone where interface reference potential densities are
            ! adjusted to achieve a more gradual change from pressure level to
            ! isopycnic interfaces.
            kt = k_range_plevel + 1
            do while (kt <= kdmx(i))
               do while (p_src(kl+1,i) < pmin(kt))
                  kl = kl + 1
               enddo
               sig_pmin(kt) = ( (p_src(kl+1,i) - pmin(kt))*sig_srcdi(1,kl) &
                              + (pmin(kt) - p_src(kl,i))*sig_srcdi(2,kl)) &
                              /(p_src(kl+1,i) - p_src(kl,i))
               if (sig_trg(kt) > sig_pmin(kt)) then
                  ktzmin = max(k_range_plevel + 2, kt - dktzu)
                  ktzmax = min(kk - 1, kt + dktzl)
                  if (ktzmin < kt .and. ktzmax - ktzmin > 1) then
                     ! For a smooth transition in layer reference potential
                     ! densities, try to construct a quadratic Bezier curve
                     ! specified by the density and density gradients at the
                     ! boundary of the transition zone. If construction of a
                     ! Bezier curve fails, use a linear change of reference
                     ! potential densities in the transition zone.
                     ckt = (sig_trg(kt) - sig_pmin(kt)) &
                           /( sig_trg(kt) - sig_trg(kt-1) &
                            - sig_pmin (kt) + sig_pmin (kt-1))
                     sig_up = sig_pmin (ktzmin-1)*ckt &
                            + sig_pmin (ktzmin  )*(1._r8 - ckt)
                     sig_lo = sig_trg(ktzmax-1)*ckt &
                            + sig_trg(ktzmax  )*(1._r8 - ckt)
                     dk = real(ktzmax - ktzmin, r8)
                     dki = 1._r8/dk
                     dsigdx_up = .5*( ( sig_pmin (ktzmin  ) &
                                      - sig_pmin (ktzmin-2))*ckt &
                                    + ( sig_pmin (ktzmin+1) &
                                      - sig_pmin (ktzmin-1))*(1. - ckt))*dk
                     dsigdx_lo = .5*( ( sig_trg(ktzmax  ) &
                                      - sig_trg(ktzmax-2))*ckt &
                                    + ( sig_trg(ktzmax+1) &
                                      - sig_trg(ktzmax-1))*(1. - ckt))*dk
                     dsigdx_up = max(0._r8, dsigdx_up)
                     if (dsigdx_lo <= dsigdx_up .or. &
                         sig_up - sig_lo <= - dsigdx_lo .or. &
                         sig_up - sig_lo >= - dsigdx_up) then
                        do k = ktzmin, ktzmax - 1
                           x = (k - ktzmin + ckt)*dki
                           sig_trg(k) = sig_up*(1._r8 - x) + sig_lo*x
                        enddo
                     else
                        xi = (sig_up - sig_lo + dsigdx_lo) &
                             /(dsigdx_lo - dsigdx_up)
                        si = ( dsigdx_lo*(sig_up + dsigdx_up) &
                             - dsigdx_up*sig_lo)/(dsigdx_lo - dsigdx_up)
                        if (abs(xi - .5_r8) < x_eps) then
                           do k = ktzmin, ktzmax-1
                              t = (k - ktzmin + ckt)*dki
                              sig_trg(k) = &
                                 (1._r8 - t)*((1._r8 - t)*sig_up + 2._r8*t*si) &
                               + t*t*sig_lo
                           enddo
                        else
                           do k = ktzmin, ktzmax-1
                              x = (k - ktzmin + ckt)*dki
                              t = (sqrt(xi*(xi - 2_r8*x) + x) - xi) &
                                  /(1._r8 - 2_r8*xi)
                              sig_trg(k) = &
                                 (1._r8 - t)*((1._r8 - t)*sig_up + 2._r8*t*si) &
                               + t*t*sig_lo
                           enddo
                        endif
                     endif
                     kt = ktzmin
                  endif
                  exit
               endif
               p_dst(kt,i) = p_src(kt,i) + nudge_fac*(pmin(kt) - p_src(kt,i))
               p_dst(kt,i) = min(max(p_dst(kt,i), pmin(kt), &
                                     p_dst(kt-1,i) + dpmin_interior), &
                                 p_src(kk+1,i))
               smooth_fac(kt,i) = 1._r8
               kt = kt + 1
            enddo

            ! Starting at the transition interface, nudge the interface
            ! pressures to reduce the deviation from the interface reference
            ! potential density.

            do k = kt, kk+1
               p_dst(k,i) = p_src(kk+1,i)
               smooth_fac(k,i) = 0._r8
            enddo

            do k = kt, min(ksmx(i), kdmx(i))
               if      (sig_trg(k) < sig_srcdi(2,k-1) .and. &
                        sig_trg(k) < sig_srcdi(1,k  )) then
                  dsig = sig_trg(k) - sig_srcdi(2,k-1)
                  dsigdx = dsigdt(t_srcdi(2,k-1,it,i), t_srcdi(2,k-1,is,i)) &
                           *dpeval1(tpc_src(:,k-1,it,i)) &
                         + dsigds(t_srcdi(2,k-1,it,i), t_srcdi(2,k-1,is,i)) &
                           *dpeval1(tpc_src(:,k-1,is,i))
                  stab_fac = dsigdx/dsig_trg(k-1)
                  dsigdx = dsig_trg(k-1)*max(stab_fac, stab_fac_limit)
                  p_dst(k,i) = p_src(k,i) &
                             + max(- .5_r8, dsig*nudge_fac/dsigdx) &
                               *(p_src(k,i) - p_src(k-1,i))
               elseif  (sig_trg(k) > sig_srcdi(2,k-1) .and. &
                        sig_trg(k) > sig_srcdi(1,k  )) then
                  dsig = sig_trg(k) - sig_srcdi(1,k)
                  dsigdx = dsigdt(t_srcdi(1,k,it,i), t_srcdi(1,k,is,i)) &
                           *dpeval0(tpc_src(:,k,it,i)) &
                         + dsigds(t_srcdi(1,k,it,i), t_srcdi(1,k,is,i)) &
                           *dpeval0(tpc_src(:,k,is,i))
                  stab_fac = dsigdx/dsig_trg(k)
                  dsigdx = dsig_trg(k)*max(stab_fac, stab_fac_limit)
                  p_dst(k,i) = p_src(k,i) &
                             + min(.5_r8, dsig*nudge_fac/dsigdx) &
                               *(p_src(k+1,i) - p_src(k,i))
               else
                  dsigdx_up = dsigdt(t_srcdi(2,k-1,it,i), t_srcdi(2,k-1,is,i)) &
                              *dpeval1(tpc_src(:,k-1,it,i)) &
                            + dsigds(t_srcdi(2,k-1,it,i), t_srcdi(2,k-1,is,i)) &
                              *dpeval1(tpc_src(:,k-1,is,i))
                  dsigdx_lo = dsigdt(t_srcdi(1,k,it,i), t_srcdi(1,k,is,i)) &
                              *dpeval0(tpc_src(:,k,it,i)) &
                            + dsigds(t_srcdi(1,k,it,i), t_srcdi(1,k,is,i)) &
                              *dpeval0(tpc_src(:,k,is,i))
                  dp_up = max(p_src(k  ,i) - p_src(k-1,i), epsilp)
                  dp_lo = max(p_src(k+1,i) - p_src(k  ,i), epsilp)
                  sig_intrp = ( (sig_srcdi(1,k  ) + .5_r8*dsigdx_lo)*dp_up &
                              + (sig_srcdi(2,k-1) - .5_r8*dsigdx_up)*dp_lo) &
                              /(dp_up + dp_lo)
                  sig_intrp = max(min(sig_srcdi(2,k-1),sig_srcdi(1,k)), &
                                  min(max(sig_srcdi(2,k-1),sig_srcdi(1,k)), &
                                          sig_intrp))
                  dsig = sig_trg(k) - sig_intrp
                  if (dsig < 0._r8) then
                     dsigdx = dsigdx_up + 2._r8*(sig_intrp - sig_srcdi(2,k-1))
                     stab_fac = dsigdx/dsig_trg(k-1)
                     dsigdx = dsig_trg(k-1)*max(stab_fac, stab_fac_limit)
                     p_dst(k,i) = p_src(k,i) &
                                + max(- .5_r8, dsig*nudge_fac/dsigdx) &
                                  *(p_src(k,i) - p_src(k-1,i))
                  else
                     dsigdx = dsigdx_lo + 2._r8*(sig_srcdi(1,k  ) - sig_intrp)
                     stab_fac = dsigdx/dsig_trg(k)
                     dsigdx = dsig_trg(k)*max(stab_fac, stab_fac_limit)
                     p_dst(k,i) = p_src(k,i) &
                                + min(.5_r8, dsig*nudge_fac/dsigdx) &
                                  *(p_src(k+1,i) - p_src(k,i))
                  endif
               endif
               p_dst(k,i) = min(max(p_dst(k,i), pmin(k), &
                                    p_dst(k-1,i) + dpmin_interior), &
                                p_src(kk+1,i))
               smooth_fac(k,i) = &
                  max(0._r8, min(1._r8, (stab_fac_limit - stab_fac) &
                                        /stab_fac_limit))
            enddo

            do k = max(kt, min(ksmx(i), kdmx(i))) + 1, kdmx(i)
               if (sig_trg(k) < sig_srcdi(2,ksmx(i))) then
                  dsig = sig_trg(k) - sig_srcdi(2,ksmx(i))
                  dsigdx = dsigdt(t_srcdi(2,ksmx(i),it,i), &
                                  t_srcdi(2,ksmx(i),is,i)) &
                           *dpeval1(tpc_src(:,ksmx(i),it,i)) &
                         + dsigds(t_srcdi(2,ksmx(i),it,i), &
                                  t_srcdi(2,ksmx(i),is,i)) &
                           *dpeval1(tpc_src(:,ksmx(i),is,i))
                  stab_fac = dsigdx/dsig_trg(ksmx(i)-1)
                  dsigdx = dsig_trg(ksmx(i)-1)*max(stab_fac, stab_fac_limit)
                  p_dst(k,i) = p_src(kk+1,i) &
                             + max(- .5_r8, dsig*nudge_fac/dsigdx) &
                               *(p_src(kk+1,i) - p_src(ksmx(i),i))
                  p_dst(k,i) = min(max(p_dst(k,i), pmin(k), &
                                       p_dst(k-1,i) + dpmin_interior), &
                                   p_src(kk+1,i))
                  smooth_fac(k,i) = &
                     max(0._r8, min(1._r8, (stab_fac_limit - stab_fac) &
                                           /stab_fac_limit))
               endif
            enddo

            ! Limit the local vertical layer thickness variation. The overall
            ! goal is that layer thickness of layer k has a lower bound of:
            !
            !    dpvar_fac*(dp(k-1) + dp(k) + dp(k+1))/3

            ks = kt
            ke = kk
            do k = kk, 1, -1
               if (p_dst(k,i) == p_dst(kk+1,i)) ke = k - 1
            enddo

            do k = ks, ke - 1
               dpmin(k) = &
                  min(2._r8*p_dst(ke+1,i) - p_dst(k+1,i) - p_dst(k,i), &
                      max(dpmin_interior, &
                          dpvar_fac*(p_dst(k+2,i) - p_dst(k-1,i))/3._r8))
            enddo

            k = ks
            do while (k < ke)
               if (p_dst(k+1,i) - p_dst(k,i) < dpmin(k)) then
                  ku = k
                  kl = k + 1
                  dpmin_sum = dpmin(ku)
                  dp_sum = p_dst(k+1,i) - p_dst(k,i)
                  layer_added = .true.
                  do while (layer_added)
                     layer_added = .false.
                     if (kl + 1 < ke) then
                        dpmin_sum_test = dpmin_sum + dpmin(kl)
                        dp_sum_test = dp_sum + p_dst(kl+1,i) - p_dst(kl,i)
                        if (dpmin_sum_test > dp_sum_test) then
                           dpmin_sum = dpmin_sum_test
                           dp_sum = dp_sum_test
                           kl = kl + 1
                           layer_added = .true.
                        endif
                     endif
                     if (ku > ks) then
                        dpmin_sum_test = dpmin_sum + dpmin(ku-1)
                        dp_sum_test = dp_sum + p_dst(ku,i) - p_dst(ku-1,i)
                        if (dpmin_sum_test > dp_sum_test) then
                           dpmin_sum = dpmin_sum_test
                           dp_sum = dp_sum_test
                           ku = ku - 1
                           layer_added = .true.
                        endif
                     endif
                  enddo
                  if (ku == ks) then
                     do k = ks, kl - 1
                        p_dst(k+1,i) = min(p_dst(ke+1,i), p_dst(k,i) + dpmin(k))
                     enddo
                     do k = kl, ke - 1
                        p_dst(k+1,i) = max(p_dst(k+1,i), p_dst(k,i))
                     enddo
                     ks = kl - 1
                     k = ks
                  else
                     dp_up = p_dst(ku,i) - p_dst(ku-1,i)
                     dp_lo = p_dst(kl+1,i) - p_dst(kl,i)
                     p_dst(ku,i) = &
                        max(p_dst(ku-1,i), &
                            p_dst(ku,i) - (dpmin_sum - dp_sum)*dp_up &
                                          /max(epsilp, dp_up + dp_lo))
                     do k = ku, kl - 1
                        p_dst(k+1,i) = min(p_dst(ke+1,i), p_dst(k,i) + dpmin(k))
                     enddo
                     k = kl
                  endif
               endif
               k = k + 1
            enddo

         enddo
      enddo

   end subroutine regrid_cntiso_hybrid_nudge_jslice

   subroutine regrid_jslice(p_src, ksmx, tpc_src, t_srcdi, p_dst, smooth_fac, &
                            ilb, iub, j, js, nn)
   ! ---------------------------------------------------------------------------
   ! Carry out regridding layer interfaces.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:,1-nbdy:), intent(in) :: p_src
      integer, dimension(1-nbdy:), intent(in) :: ksmx
      real(r8), dimension(:,:,:,1-nbdy:), intent(in) :: tpc_src, t_srcdi
      real(r8), dimension(:,1-nbdy:), intent(out) :: p_dst, smooth_fac
      integer, intent(in) :: ilb, iub, j, js, nn

      if (vcoord_tag == vcoord_plevel) then
         call regrid_plevel_jslice(p_src, p_dst, ilb, iub, j)
      else
         if (regrid_method_tag == regrid_method_direct) then
            call regrid_cntiso_hybrid_direct_jslice(p_src, p_dst, &
                                                    ilb, iub, j, js, nn)
         else
            call regrid_cntiso_hybrid_nudge_jslice(p_src, ksmx, tpc_src, &
                                                   t_srcdi, p_dst, smooth_fac, &
                                                   ilb, iub, j)
         endif
      endif

   end subroutine regrid_jslice

   subroutine regrid_smooth_jslice(p_dst_js, smooth_fac_js, smtflxconv_js, &
                                   ilb, iub, j, js2, js3)
   ! ---------------------------------------------------------------------------
   ! For vcoord == 'cntiso_hybrid' and regrid_method == 'nudge', apply lateral
   ! smoothing of the regridded interfaces when a vertical stability factor is
   ! below a specified threshold.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:,1-nbdy:,:), intent(inout) :: &
         p_dst_js, smooth_fac_js, smtflxconv_js
      integer, intent(in) :: ilb, iub, j, js2, js3

      real(r8) :: cdiff, difmx, flxhi, flxlo, flx, sdiff
      integer :: l, i, k

      smtflxconv_js(:,:,js3) = 0._r8

      do l = 1, isu(j+1)
         do i = max(ilb, ifu(j+1,l)), min(iub+1, ilu(j+1,l))
            cdiff = delt1*scuy(i,j+1)*scuxi(i,j+1)
            difmx = .5_r8*(difmxp(i-1,j+1) + difmxp(i,j+1))
            do k = 2, kk
               flxhi =   .125_r8*min(( p_dst_js(k,  i-1,js3) &
                                     - p_dst_js(k-1,i-1,js3))*scp2(i-1,j+1), &
                                     ( p_dst_js(k+1,i  ,js3) &
                                     - p_dst_js(k  ,i  ,js3))*scp2(i  ,j+1))
               flxlo = - .125_r8*min(( p_dst_js(k  ,i  ,js3) &
                                     - p_dst_js(k-1,i  ,js3))*scp2(i  ,j+1), &
                                     ( p_dst_js(k+1,i-1,js3) &
                                     - p_dst_js(k  ,i-1,js3))*scp2(i-1,j+1))
               sdiff = min(.5_r8*( smooth_fac_js(k,i-1,js3) &
                                 + smooth_fac_js(k,i  ,js3))*smooth_diff_max, &
                           difmx)
               flx = min(flxhi, max(flxlo, cdiff*sdiff*( p_dst_js(k,i-1,js3) &
                                                       - p_dst_js(k,i  ,js3))))
               smtflxconv_js(k,i-1,js3) = smtflxconv_js(k,i-1,js3) + flx
               smtflxconv_js(k,i  ,js3) = smtflxconv_js(k,i  ,js3) - flx
            enddo
         enddo
      enddo

      do l = 1, isv(j+1)
         do i = max(ilb, ifv(j+1,l)), min(iub, ilv(j+1,l))
            cdiff = delt1*scvx(i,j+1)*scvyi(i,j+1)
            difmx = .5_r8*(difmxp(i,j) + difmxp(i,j+1))
            do k = 2, kk
               flxhi =   .125_r8*min(( p_dst_js(k,  i,js2) &
                                     - p_dst_js(k-1,i,js2))*scp2(i,j  ), &
                                     ( p_dst_js(k+1,i,js3) &
                                     - p_dst_js(k  ,i,js3))*scp2(i,j+1))
               flxlo = - .125_r8*min(( p_dst_js(k  ,i,js3) &
                                     - p_dst_js(k-1,i,js3))*scp2(i,j+1), &
                                     ( p_dst_js(k+1,i,js2) &
                                     - p_dst_js(k  ,i,js2))*scp2(i,j  ))
               sdiff = min(.5_r8*( smooth_fac_js(k,i,js2) &
                                 + smooth_fac_js(k,i,js3))*smooth_diff_max, &
                           difmx)
               flx = min(flxhi, max(flxlo, cdiff*sdiff*( p_dst_js(k,i,js2) &
                                                       - p_dst_js(k,i,js3))))
               smtflxconv_js(k,i,js2) = smtflxconv_js(k,i,js2) + flx
               smtflxconv_js(k,i,js3) = smtflxconv_js(k,i,js3) - flx
            enddo
         enddo
      enddo

      do l = 1, isp(j)
         do i = max(ilb, ifp(j,l)), min(iub, ilp(j,l))
            do k = 2, kk
               p_dst_js(k,i,js2) = p_dst_js(k,i,js2) &
                                 - smtflxconv_js(k,i,js2)*scp2i(i,j)
            enddo
         enddo
      enddo

   end subroutine regrid_smooth_jslice

   subroutine remap_trc_jslice(p_dst, trc_rm, ilb, iub, j, js)
   ! ---------------------------------------------------------------------------
   ! Remap tracers from source to destination layer structure.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:,1-nbdy:), intent(in) :: p_dst
      real(r8), dimension(:,:,1-nbdy:), intent(out) :: trc_rm
      integer, intent(in) :: ilb, iub, j, js

      integer :: l, i, nt, errstat

      do l = 1, isp(j)
         do i = max(ilb, ifp(j,l)), min(iub, ilp(j,l))

            ! Prepare remapping to destination layers.
            errstat = prepare_remapping(rcgs, rms, p_dst(:,i), i, js)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(remap_trc_jslice)')
               stop '(remap_trc_jslice)'
            endif

            ! Remap tracers.
            do nt = 1, ntr_loc
               errstat = remap(trc_rcss(nt), rms, trc_rm(:,nt,i), i, js)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(remap_trc_jslice)')
                  stop '(remap_trc_jslice)'
               endif
            enddo

         enddo
      enddo

   end subroutine remap_trc_jslice

   subroutine remap_trc_diazlv_jslice(p_src, ilb, iub, j, js, &
                                      do_acc_templvl, do_acc_salnlvl, &
                                      do_acc_idlagelvl)
   ! ---------------------------------------------------------------------------
   ! Remap tracers from source to diagnostic z-levels.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:,1-nbdy:), intent(in) :: p_src
      integer, intent(in) :: ilb, iub, j, js
      logical, intent(in) :: do_acc_templvl, do_acc_salnlvl, do_acc_idlagelvl

      real(r8), dimension(ddm+1) :: p_dst
      real(r8), dimension(ddm) :: trc_rm
      real(r8) :: q
      integer :: l, i, kd, errstat, iogrp

      if (.not. (do_acc_templvl .or. do_acc_salnlvl .or. &
                 do_acc_idlagelvl)) return

      do l = 1, isp(j)
         do i = max(ilb, ifp(j,l)), min(iub, ilp(j,l))

            ! Prepare remapping to destination z-levels.
            p_dst(1) = p_src(1,i)
            q = 1._r8/pbath(i,j)
            do kd = 2, ddm
               p_dst(kd) = min(depthslev_bnds(1,kd)*q, 1._r8) &
                           *(p_src(kk+1,i) - p_src(1,i)) + p_src(1,i)
            enddo
            p_dst(ddm+1) = p_src(kk+1,i)
            errstat = prepare_remapping(rcgs, rms_diazlv, p_dst, i, js)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(remap_trc_diazlv_jslice)')
               stop '(remap_trc_diazlv_jslice)'
            endif

            ! Remap tracers.

            if (do_acc_templvl) then
               errstat = remap(trc_rcss(1), rms_diazlv, trc_rm, i, js)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(remap_trc_diazlv_jslice)')
                  stop '(remap_trc_diazlv_jslice)'
               endif
               do iogrp = 1, nphy
                  if (acc_templvl(iogrp) /= 0) then
                     do kd = 1, ddm
                        phylvl(i,j,kd,acc_templvl(iogrp)) = &
                           phylvl(i,j,kd,acc_templvl(iogrp)) + trc_rm(kd)
                     enddo
                  endif
               enddo
            endif

            if (do_acc_salnlvl) then
               errstat = remap(trc_rcss(2), rms_diazlv, trc_rm, i, js)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(remap_trc_diazlv_jslice)')
                  stop '(remap_trc_diazlv_jslice)'
               endif
               do iogrp = 1, nphy
                  if (acc_salnlvl(iogrp) /= 0) then
                     do kd = 1, ddm
                        phylvl(i,j,kd,acc_salnlvl(iogrp)) = &
                           phylvl(i,j,kd,acc_salnlvl(iogrp)) + trc_rm(kd)
                     enddo
                  endif
               enddo
            endif

            if (do_acc_idlagelvl) then
               errstat = remap(trc_rcss(itriag+2), rms_diazlv, trc_rm, i, js)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(remap_trc_diazlv_jslice)')
                  stop '(remap_trc_diazlv_jslice)'
               endif
               do iogrp = 1, nphy
                  if (acc_idlagelvl(iogrp) /= 0) then
                     do kd = 1, ddm
                        phylvl(i,j,kd,acc_idlagelvl(iogrp)) = trc_rm(kd)
                     enddo
                  endif
               enddo
            endif

         enddo
      enddo

   end subroutine remap_trc_diazlv_jslice

   subroutine copy_jslice_to_3d(p_dst, trc_rm, ilb, iub, j, nn)

      real(r8), dimension(:,1-nbdy:), intent(in) :: p_dst
      real(r8), dimension(:,:,1-nbdy:), intent(in) :: trc_rm

      integer, intent(in) :: ilb, iub, j, nn

      integer :: l, i, k, kn, nt

      do l = 1, isp(j)
         do i = max(ilb, ifp(j,l)), min(iub, ilp(j,l))

            do k = 1, kk
               kn = k + nn
               temp(i,j,kn) = trc_rm(k,1,i)
               saln(i,j,kn) = trc_rm(k,2,i)
               dp(i,j,kn) = p_dst(k+1,i) - p_dst(k,i)
               sigma(i,j,kn) = sig(trc_rm(k,1,i), trc_rm(k,2,i))
               do nt = 1, ntr
                  trc(i,j,kn,nt) = trc_rm(k,nt+2,i)
               enddo
            enddo

         enddo
      enddo

   end subroutine copy_jslice_to_3d

   ! ---------------------------------------------------------------------------
   ! Public procedures.
   ! ---------------------------------------------------------------------------

   subroutine readnml_ale_regrid_remap
   ! ---------------------------------------------------------------------------
   ! Read variables in the namelist group 'ale_regrid_remap' and resolve
   ! options.
   ! ---------------------------------------------------------------------------

      character(len = 80) :: nml_fname
      integer :: nfu, ios
      logical :: fexist

      namelist /ale_regrid_remap/ &
         reconstruction_method, upper_bndr_ord, lower_bndr_ord, &
         density_limiting, tracer_limiting, velocity_limiting, &
         density_pc_upper_bndr, density_pc_lower_bndr, &
         tracer_pc_upper_bndr, tracer_pc_lower_bndr, &
         velocity_pc_upper_bndr, velocity_pc_lower_bndr, dpmin_interior, &
         regrid_method, k_range_plevel, regrid_nudge_ts, stab_fac_limit, &
         dpvar_fac, smooth_diff_max, dktzu, dktzl

      ! Return if ALE method is not required.
      if (vcoord_tag == vcoord_isopyc_bulkml) return

      ! Read variables in the namelist group 'ale_regrid_remap'.
      if (mnproc == 1) then
         nml_fname = 'ocn_in'//trim(inst_suffix)
         inquire(file = nml_fname, exist = fexist)
         if (fexist) then
            open (newunit = nfu, file = nml_fname, status = 'old', &
                  action = 'read')
         else
            nml_fname = 'limits'//trim(inst_suffix)
            inquire(file = nml_fname, exist = fexist)
            if (fexist) then
               open (newunit = nfu, file = nml_fname, status = 'old', &
                     action = 'read')
            else
               write (lp,*) &
                  'readnml_ale_regrid_remap: could not find namelist file!'
               call xchalt('(readnml_ale_regrid_remap)')
               stop '(readnml_ale_regrid_remap)'
            endif
         endif
         read (unit = nfu, nml = ale_regrid_remap, iostat = ios)
         close (unit = nfu)
      endif
      call xcbcst(ios)
      if (ios /= 0) then
         if (mnproc == 1) &
              write (lp,*) &
                 'readnml_ale_regrid_remap: No vertical coordinate variable group found in namelist. Using defaults.'
      else
         call xcbcst(reconstruction_method)
         call xcbcst(upper_bndr_ord)
         call xcbcst(lower_bndr_ord)
         call xcbcst(density_limiting)
         call xcbcst(tracer_limiting)
         call xcbcst(velocity_limiting)
         call xcbcst(density_pc_upper_bndr)
         call xcbcst(density_pc_lower_bndr)
         call xcbcst(tracer_pc_upper_bndr)
         call xcbcst(tracer_pc_lower_bndr)
         call xcbcst(velocity_pc_upper_bndr)
         call xcbcst(velocity_pc_lower_bndr)
         call xcbcst(dpmin_interior)
         call xcbcst(regrid_method)
         call xcbcst(k_range_plevel)
         call xcbcst(regrid_nudge_ts)
         call xcbcst(stab_fac_limit)
         call xcbcst(dpvar_fac)
         call xcbcst(smooth_diff_max)
         call xcbcst(dktzu)
         call xcbcst(dktzl)
      endif
      if (mnproc == 1) then
         write (lp,*) 'readnml_ale_regrid_remap: ALE regrid-remap variables:'
         write (lp,*) '  reconstruction_method =  ', trim(reconstruction_method)
         write (lp,*) '  upper_bndr_ord =         ', upper_bndr_ord
         write (lp,*) '  lower_bndr_ord =         ', lower_bndr_ord
         write (lp,*) '  density_limiting =       ', trim(density_limiting)
         write (lp,*) '  tracer_limiting =        ', trim(tracer_limiting)
         write (lp,*) '  velocity_limiting =      ', trim(velocity_limiting)
         write (lp,*) '  density_pc_upper_bndr =  ', density_pc_upper_bndr
         write (lp,*) '  density_pc_lower_bndr =  ', density_pc_lower_bndr
         write (lp,*) '  tracer_pc_upper_bndr =   ', tracer_pc_upper_bndr
         write (lp,*) '  tracer_pc_lower_bndr =   ', tracer_pc_lower_bndr
         write (lp,*) '  velocity_pc_upper_bndr = ', velocity_pc_upper_bndr
         write (lp,*) '  velocity_pc_lower_bndr = ', velocity_pc_lower_bndr
         write (lp,*) '  dpmin_interior =         ', dpmin_interior
         write (lp,*) '  regrid_method =          ', trim(regrid_method)
         write (lp,*) '  k_range_plevel =         ', k_range_plevel
         write (lp,*) '  regrid_nudge_ts =        ', regrid_nudge_ts
         write (lp,*) '  stab_fac_limit =         ', stab_fac_limit
         write (lp,*) '  dpvar_fac =              ', dpvar_fac
         write (lp,*) '  smooth_diff_max =        ', smooth_diff_max
         write (lp,*) '  dktzu =                  ', dktzu
         write (lp,*) '  dktzl =                  ', dktzl
      endif

      ! Resolve options.
      select case (trim(reconstruction_method))
         case ('plm')
            reconstruction_method_tag = hor3map_plm
         case ('ppm')
            reconstruction_method_tag = hor3map_ppm
         case ('pqm')
            reconstruction_method_tag = hor3map_pqm
         case default
            if (mnproc == 1) &
                 write (lp,'(3a)') &
                 ' readnml_ale_regrid_remap: reconstruction_method = ', &
                 trim(reconstruction_method), ' is unsupported!'
            call xcstop('(readnml_ale_regrid_remap)')
            stop '(readnml_ale_regrid_remap)'
      end select
      select case (trim(density_limiting))
         case ('monotonic')
            density_limiting_tag = hor3map_monotonic
         case default
            if (mnproc == 1) &
                 write (lp,'(3a)') &
                 ' readnml_ale_regrid_remap: density_limiting = ', &
                 trim(density_limiting), ' is unsupported!'
            call xcstop('(readnml_ale_regrid_remap)')
            stop '(readnml_ale_regrid_remap)'
      end select
      select case (trim(tracer_limiting))
         case ('monotonic')
            tracer_limiting_tag = hor3map_monotonic
         case ('non_oscillatory')
            tracer_limiting_tag = hor3map_non_oscillatory
         case default
            if (mnproc == 1) &
                 write (lp,'(3a)') &
                 ' readnml_ale_regrid_remap: tracer_limiting = ', &
                 trim(tracer_limiting), ' is unsupported!'
            call xcstop('(readnml_ale_regrid_remap)')
            stop '(readnml_ale_regrid_remap)'
      end select
      select case (trim(velocity_limiting))
         case ('monotonic')
            velocity_limiting_tag = hor3map_monotonic
         case ('non_oscillatory')
            velocity_limiting_tag = hor3map_non_oscillatory
         case default
            if (mnproc == 1) &
                 write (lp,'(3a)') &
                 ' readnml_ale_regrid_remap: velocity_limiting = ', &
                 trim(velocity_limiting), ' is unsupported!'
            call xcstop('(readnml_ale_regrid_remap)')
            stop '(readnml_ale_regrid_remap)'
      end select
      if (vcoord_tag == vcoord_cntiso_hybrid) then
         select case (trim(regrid_method))
            case ('direct')
               regrid_method_tag = regrid_method_direct
            case ('nudge')
               regrid_method_tag = regrid_method_nudge
            case default
               if (mnproc == 1) &
                    write (lp,'(3a)') &
                    ' readnml_ale_regrid_remap: regrid_method = ', &
                    trim(regrid_method), ' is unsupported!'
               call xcstop('(readnml_ale_regrid_remap)')
               stop '(readnml_ale_regrid_remap)'
         end select
      endif

      ! Change units from [m] to [kg m-1 s-2] of depth interval variables.
      dpmin_interior = dpmin_interior*onem

   end subroutine readnml_ale_regrid_remap

   subroutine init_ale_regrid_remap
   ! ---------------------------------------------------------------------------
   ! Initialize arrays and data structures.
   ! ---------------------------------------------------------------------------

      integer :: i, j, k, nt, errstat

      ! Only do the initialization if the vertical coordinate require the ALE
      ! approach.
      if (vcoord_tag == vcoord_isopyc_bulkml) return

      ! Local number of tracers where temperature and salinity is added to the
      ! ntr parameter.
      ntr_loc = ntr + 2

      ! Allocate reconstruction data structures for tracer source data.
      allocate(trc_rcss(ntr_loc), stat = errstat)
      if (errstat /= 0) then
         write(lp,*) 'Failed to allocate trc_rcss!'
         call xchalt('(init_ale_regrid_remap)')
         stop '(init_ale_regrid_remap)'
      endif

      ! Configuration of the reconstruction data structure that only depends on
      ! the source grid.
      rcgs%n_src = kk
      rcgs%i_ubound = ii
      if (ltedtp_opt == ltedtp_neutral) then
         ! Neutral diffusion is requested so increase the index range of the
         ! reconstruction data structure.
         rcgs%i_lbound = rcgs%i_lbound - 1
         rcgs%i_ubound = rcgs%i_ubound + 1
         rcgs%j_ubound = rcgs%j_ubound + 1
      endif
      if (vcoord_tag == vcoord_cntiso_hybrid .and. &
          regrid_method_tag == regrid_method_nudge .and. &
          smooth_diff_max > 0._r8) then
         ! Lateral smoothing of interfaces after regridding is requested so
         ! increase the index range of the reconstruction data structure.
         rcgs%i_lbound = rcgs%i_lbound - 1
         rcgs%i_ubound = rcgs%i_ubound + 1
         rcgs%j_ubound = rcgs%j_ubound + 1
      endif
      rcgs%method = reconstruction_method_tag
      rcgs%left_bndr_ord = upper_bndr_ord
      rcgs%right_bndr_ord = lower_bndr_ord

      ! Configuration of reconstruction data structures that is specific to
      ! various source data.

      d_rcss%limiting = density_limiting_tag
      d_rcss%pc_left_bndr = density_pc_upper_bndr
      d_rcss%pc_right_bndr = density_pc_lower_bndr

      trc_rcss(1)%limiting = tracer_limiting_tag
      trc_rcss(1)%pc_left_bndr = tracer_pc_upper_bndr
      trc_rcss(1)%pc_right_bndr = tracer_pc_lower_bndr
      if (tracer_limiting_tag == hor3map_non_oscillatory) then
         do nt = 2, ntr_loc
            trc_rcss(nt)%limiting = hor3map_non_oscillatory_posdef
            trc_rcss(nt)%pc_left_bndr = tracer_pc_upper_bndr
            trc_rcss(nt)%pc_right_bndr = tracer_pc_lower_bndr
         enddo
      else
         do nt = 2, ntr_loc
            trc_rcss(nt)%limiting = tracer_limiting_tag
            trc_rcss(nt)%pc_left_bndr = tracer_pc_upper_bndr
            trc_rcss(nt)%pc_right_bndr = tracer_pc_lower_bndr
         enddo
      endif

      v_rcss%limiting = velocity_limiting_tag
      v_rcss%pc_left_bndr = velocity_pc_upper_bndr
      v_rcss%pc_right_bndr = velocity_pc_lower_bndr

      ! Configuration of remapping data structure.
      rms%n_dst = kk

      ! Configuration of remapping data structure for diagnostic z-levels.
      rms_diazlv%n_dst = ddm

      ! Initialize reconstruction and remapping data structures.

      errstat = initialize_rcgs(rcgs)
      if (errstat /= hor3map_noerr) then
         write(lp,*) trim(hor3map_errstr(errstat))
         call xchalt('(init_ale_regrid_remap)')
         stop '(init_ale_regrid_remap)'
      endif

      errstat = initialize_rcss(rcgs, d_rcss)
      if (errstat /= hor3map_noerr) then
         write(lp,*) trim(hor3map_errstr(errstat))
         call xchalt('(init_ale_regrid_remap)')
         stop '(init_ale_regrid_remap)'
      endif

      do nt = 1, ntr_loc
         errstat = initialize_rcss(rcgs, trc_rcss(nt))
         if (errstat /= hor3map_noerr) then
            write(lp,*) trim(hor3map_errstr(errstat))
            call xchalt('(init_ale_regrid_remap)')
            stop '(init_ale_regrid_remap)'
         endif
      enddo

      errstat = initialize_rcss(rcgs, v_rcss)
      if (errstat /= hor3map_noerr) then
         write(lp,*) trim(hor3map_errstr(errstat))
         call xchalt('(init_ale_regrid_remap)')
         stop '(init_ale_regrid_remap)'
      endif

      errstat = initialize_rms(rcgs, rms)
      if (errstat /= hor3map_noerr) then
         write(lp,*) trim(hor3map_errstr(errstat))
         call xchalt('(init_ale_regrid_remap)')
         stop '(init_ale_regrid_remap)'
      endif

      errstat = initialize_rms(rcgs, rms_diazlv)
      if (errstat /= hor3map_noerr) then
         write(lp,*) trim(hor3map_errstr(errstat))
         call xchalt('(init_ale_regrid_remap)')
         stop '(init_ale_regrid_remap)'
      endif

   end subroutine init_ale_regrid_remap

   subroutine ale_regrid_remap(m, n, mm, nn, k1m, k1n)
   ! ---------------------------------------------------------------------------
   ! Regrid, remap and carry out additional operations that make use of the
   ! reconstructed vertical profiles, such as neutral diffusion and remapping to
   ! diagnostic z-levels.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      integer, parameter :: p_ord = 4

      real(r8), dimension(kdm,1-nbdy:idm+nbdy,3) :: smtflxconv_js
      real(r8), dimension(kdm+1,1-nbdy:idm+nbdy,3) :: &
         p_src_js, p_dst_js, smooth_fac_js
      real(r8), dimension(p_ord+1,kdm,ntr_loc,1-nbdy:idm+nbdy,3) :: tpc_src_js
      real(r8), dimension(2,kdm,ntr_loc,1-nbdy:idm+nbdy,3) :: t_srcdi_js
      real(r8), dimension(2,kdm,1-nbdy:idm+nbdy,3) :: &
         p_srcdi_js, drhodt_srcdi_js, drhods_srcdi_js
      real(r8), dimension(kdm,ntr_loc,1-nbdy:idm+nbdy,3) :: flxconv_js
      real(r8), dimension(kdm,ntr_loc,1-nbdy:idm+nbdy) :: trc_rm
      real(r8), dimension(kdm+1) :: p_src, p_dst
      real(r8), dimension(kdm) :: v_src, v_rm
      real(r8), dimension(ddm+1) :: p_dst_diazlv
      real(r8), dimension(ddm) :: v_rm_diazlv
      real(r8) :: q
      integer, dimension(1-nbdy:idm+nbdy,3) :: ksmx_js, kdmx_js
      integer :: ilb1, ilb2, ilb3, iub1, iub2, iub3, jofs2, jofs3, &
                 jlb_regrid_smooth, jlb_ndiff_prep, jlb_ndiff_uflx, &
                 jlb_ndiff_vflx, jlb_ndiff_update_trc, &
                 js1, js2, js3, j, nt, i, k, l, kn, kd, iogrp, errstat
      logical :: do_regrid_smooth, do_ndiff, do_acc_templvl, do_acc_salnlvl, &
                 do_acc_uvellvl, do_acc_vvellvl, do_acc_idlagelvl
      character(len = 2) cnt

      ! ------------------------------------------------------------------------
      ! Check if accumulation of diagnostic variables should be done.
      ! ------------------------------------------------------------------------

      do_acc_templvl = sum(acc_templvl(1:nphy)) /= 0 .and. &
                       sum(acc_templvl(1:nphy)*alarm_phy(1:nphy)) == 0
      do_acc_salnlvl = sum(acc_salnlvl(1:nphy)) /= 0 .and. &
                       sum(acc_salnlvl(1:nphy)*alarm_phy(1:nphy)) == 0
      do_acc_uvellvl = sum(acc_uvellvl(1:nphy)) /= 0 .and. &
                       sum(acc_uvellvl(1:nphy)*alarm_phy(1:nphy)) == 0
      do_acc_vvellvl = sum(acc_vvellvl(1:nphy)) /= 0 .and. &
                       sum(acc_vvellvl(1:nphy)*alarm_phy(1:nphy)) == 0
      do_acc_idlagelvl = .false.

      ! ------------------------------------------------------------------------
      ! Regrid and remap tracers. Also carry out neutral diffusion if requested.
      ! The operations are carried out on temporary arrays with minimal j-index
      ! range, j-slices, to minimize the memory needed to hold the
      ! reconstruction data structures.
      ! ------------------------------------------------------------------------

      ! Check if lateral smoothing of interfaces after regridding should be
      ! applied.
      if (vcoord_tag == vcoord_cntiso_hybrid .and. &
          regrid_method_tag == regrid_method_nudge .and. &
          smooth_diff_max > 0._r8) then
         do_regrid_smooth = .true.
         smtflxconv_js(:,:,:) = 0._r8
      else
         do_regrid_smooth = .false.
      endif

      ! Check if neutral diffusion is requested.
      if (ltedtp_opt == ltedtp_neutral) then
         do_ndiff = .true.
      else
         do_ndiff = .false.
      endif

      ! Set j-slice offsets.
      if (do_ndiff) then
         jofs2 = 1
      else
         jofs2 = 0
      endif
      if (do_regrid_smooth) then
         jofs3 = jofs2 + 1
      else
         jofs3 = jofs2
      endif

      ! Set lower j-index bound of the various operations. If an operation is
      ! not requested, the bound is set so the operation is never called.
      jlb_regrid_smooth = jj + 1
      jlb_ndiff_prep = jj + 1
      jlb_ndiff_uflx = jj + 1
      jlb_ndiff_vflx = jj + 1
      jlb_ndiff_update_trc = jj + 1
      if (do_regrid_smooth) then
         jlb_regrid_smooth = 1 - 2*jofs3 + 1
      endif
      if (do_ndiff) then
         jlb_ndiff_prep = - 1
         jlb_ndiff_uflx = 1
         jlb_ndiff_vflx = 0
         jlb_ndiff_update_trc = 1
      endif

      ! Set lower and upper bounds of i-indices corresponding to j-slice
      ! offsets.
      ilb1 =  1
      ilb2 =  1 - jofs2
      ilb3 =  1 - jofs3
      iub1 = ii
      iub2 = ii + jofs2
      iub3 = ii + jofs3

      ! Update halos as needed.
      if (jofs3 > 0) then
         call xctilr(temp (1-nbdy,1-nbdy,k1n), 1, kk, jofs3, jofs3, halo_ps)
         call xctilr(saln (1-nbdy,1-nbdy,k1n), 1, kk, jofs3, jofs3, halo_ps)
         call xctilr(sigma(1-nbdy,1-nbdy,k1n), 1, kk, jofs3, jofs3, halo_ps)
      endif
      if (do_ndiff) then
         do nt = 1, ntr
            call xctilr(trc(1-nbdy,1-nbdy,k1n,nt), 1, kk, 1, 1, halo_ps)
         enddo
      end if

      ! Initial j-slice indices.
      js1 = 1
      js2 = js1 + jofs2
      js3 = js1 + jofs3

      do j = 1 - 2*jofs3, jj

         ! Update j-slice indices.
         js1 = mod(js1            , jofs3 + 1) + 1
         js2 = mod(js1 + jofs2 - 1, jofs3 + 1) + 1
         js3 = mod(js1 + jofs3 - 1, jofs3 + 1) + 1

         ! Vertically reconstruct tracers.
         call reconstruct_trc_jslice(p_src_js(:,:,js3), ksmx_js(:,js3), &
                                     tpc_src_js(:,:,:,:,js3), &
                                     t_srcdi_js(:,:,:,:,js3), &
                                     ilb3, iub3, j+jofs3, js3, nn)

         ! Regrid.
         call regrid_jslice(p_src_js(:,:,js3), ksmx_js(:,js3), &
                            tpc_src_js(:,:,:,:,js3), &
                            t_srcdi_js(:,:,:,:,js3), &
                            p_dst_js(:,:,js3), &
                            smooth_fac_js(:,:,js3), &
                            ilb3, iub3, j+jofs3, js3, nn)

         ! If requested, apply lateral smoothing of the interfaces after
         ! regridding.
         if (j >= jlb_regrid_smooth) &
            call regrid_smooth_jslice(p_dst_js, smooth_fac_js, smtflxconv_js, &
                                      ilb2, iub2, j+jofs2, js2, js3)

         ! If requested, prepare neutral diffusion.
         if (j >= jlb_ndiff_prep) &
            call ndiff_prep_jslice(p_src_js, ksmx_js, &
                                   tpc_src_js, t_srcdi_js, &
                                   p_dst_js, kdmx_js, p_srcdi_js, &
                                   drhodt_srcdi_js, drhods_srcdi_js, &
                                   flxconv_js, &
                                   ilb2, iub2, j+jofs2, js2, mm)

         ! If requested, compute the contribution of u-component fluxes to the
         ! flux convergence of neutral diffusion.
         if (j >= jlb_ndiff_uflx) &
            call ndiff_uflx_jslice(ksmx_js, tpc_src_js, t_srcdi_js, &
                                   p_dst_js, kdmx_js, p_srcdi_js, &
                                   drhodt_srcdi_js, drhods_srcdi_js, &
                                   flxconv_js, &
                                   ntr_loc, ilb1, iub2, j, js1, mm, nn)

         ! If requested, compute the contribution of v-component fluxes to the
         ! flux convergence of neutral diffusion.
         if (j >= jlb_ndiff_vflx) &
            call ndiff_vflx_jslice(ksmx_js, tpc_src_js, t_srcdi_js, &
                                   p_dst_js, kdmx_js, p_srcdi_js, &
                                   drhodt_srcdi_js, drhods_srcdi_js, &
                                   flxconv_js, &
                                   ntr_loc, ilb1, iub1, j+jofs2, &
                                   js1, js2, mm, nn)

         ! Remap tracers to regridded layers and diagnostic z-levels.
         if (j >= 1) then
            call remap_trc_jslice(p_dst_js(:,:,js1), trc_rm, &
                                  ilb1, iub1, j, js1)
            call remap_trc_diazlv_jslice(p_src_js(:,:,js1), &
                                         ilb1, iub1, j, js1, &
                                         do_acc_templvl, do_acc_salnlvl, &
                                         do_acc_idlagelvl)
         endif

         ! If requested, update the tracers by applying the neutral diffusion
         ! flux convergence.
         if (j >= jlb_ndiff_update_trc) &
            call ndiff_update_trc_jslice(p_dst_js, flxconv_js, trc_rm, &
                                         ntr_loc, ilb1, iub1, j, js1)

         ! Copy from the j-slice array to the full tracer array.
         if (j >= 1) &
            call copy_jslice_to_3d(p_dst_js(:,:,js1), trc_rm, &
                                   ilb1, iub1, j, nn)
      enddo

      ! ------------------------------------------------------------------------
      ! Remap velocity.
      ! ------------------------------------------------------------------------

      !$omp parallel do private(k, kn, l, i)
      do j = 1, jj
         do k = 1, kk
            kn = k + nn
            do l = 1, isu(j)
               do i = max(1, ifu(j,l)), min(ii, ilu(j,l))
                  pu(i,j,k+1) = pu(i,j,k) + dpu(i,j,kn)
               enddo
            enddo
            do l = 1, isv(j)
               do i = max(1, ifv(j,l)), min(ii, ilv(j,l))
                  pv(i,j,k+1) = pv(i,j,k) + dpv(i,j,kn)
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do

      call xctilr(dp(1-nbdy,1-nbdy,k1n), 1, kk, 3, 3, halo_ps)

      !$omp parallel do private(k, kn, l, i)
      do j = -2, jj+3
          do l = 1, isp(j)
             do i = max(-1, ifp(j,l)), min(ii, ilp(j,l))
                util1(i,j) = p(i,j,kk+1)
             enddo
          enddo
         do k = 1, kk
            kn = k + nn
            do l = 1, isp(j)
               do i = max(-2, ifp(j,l)), min(ii+3, ilp(j,l))
                  p(i,j,k+1) = p(i,j,k) + dp(i,j,kn)
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do

      !$omp parallel do private(k,kn,l,i,q)
      do j = -1, jj+2
         do k = 1, kk
            kn = k + nn
            do l = 1, isu(j)
               do i = max(-1, ifu(j,l)), min(ii+2, ilu(j,l))
                  q = min(p(i,j,kk+1), p(i-1,j,kk+1))
                  dpu(i,j,kn) = &
                       .5_r8*( (min(q, p(i-1,j,k+1)) - min(q, p(i-1,j,k))) &
                             + (min(q, p(i  ,j,k+1)) - min(q, p(i  ,j,k))))
               enddo
            enddo
            do l = 1, isv(j)
               do i = max(-1, ifv(j,l)), min(ii+2, ilv(j,l))
                  q = min(p(i,j,kk+1), p(i,j-1,kk+1))
                  dpv(i,j,kn) = &
                       .5_r8*( (min(q, p(i,j-1,k+1)) - min(q, p(i,j-1,k))) &
                             + (min(q, p(i,j  ,k+1)) - min(q, p(i,j  ,k))))
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do

      do j = 1, jj

         do l = 1, isu(j)
            do i = max(1, ifu(j,l)), min(ii, ilu(j,l))

               ! Copy variables into 1D arrays. Rescale source interfaces so the
               ! pressure range of source and destination columns match.
               p_dst(1) = pu(i,j,1)
               do k = 1, kk
                  kn = k + nn
                  v_src(k) = u(i,j,kn)
                  p_dst(k+1) = p_dst(k) + dpu(i,j,kn)
               enddo
               q = min(util1(i-1,j), util1(i,j))/pu(i,j,kk+1)
               do k = 1, kk+1
                  p_src(k) = pu(i,j,k)*q
               enddo

               ! Prepare reconstruction with current interface pressures.
               errstat = prepare_reconstruction(rcgs, p_src, i, 1)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(ale_regrid_remap)')
                  stop '(ale_regrid_remap)'
               endif

               ! Prepare remapping to layer structure with regridded interface
               ! pressures.
               errstat = prepare_remapping(rcgs, rms, p_dst, i, 1)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(ale_regrid_remap)')
                  stop '(ale_regrid_remap)'
               endif

               ! Reconstruct u-component of velocity.
               errstat = reconstruct(rcgs, v_rcss, v_src, i, 1)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(ale_regrid_remap)')
                  stop '(ale_regrid_remap)'
               endif

               ! Remap u-component of velocity to regridded layers.
               errstat = remap(v_rcss, rms, v_rm, i, 1)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(ale_regrid_remap)')
                  stop '(ale_regrid_remap)'
               endif

               ! Update 3D array.
               do k = 1, kk
                  kn = k + nn
                  u(i,j,kn) = v_rm(k)
               enddo

               if (do_acc_uvellvl) then

                  ! Prepare remapping to destination z-levels.
                  p_dst_diazlv(1) = p_src(1)
                  q = 1._r8/ubath(i,j)
                  do kd = 2, ddm
                     p_dst_diazlv(kd) = min(depthslev_bnds(1,kd)*q, 1._r8) &
                                        *(p_src(kk+1) - p_src(1)) + p_src(1)
                  enddo
                  p_dst_diazlv(ddm+1) = p_src(kk+1)
                  errstat = prepare_remapping(rcgs, rms_diazlv, p_dst_diazlv, &
                                              i, 1)
                  if (errstat /= hor3map_noerr) then
                     write(lp,*) trim(hor3map_errstr(errstat))
                     call xchalt('(ale_regrid_remap)')
                     stop '(ale_regrid_remap)'
                  endif

                  ! Remap u-component of velocity to diagnostic z-levels and
                  ! accumulate.
                  errstat = remap(v_rcss, rms_diazlv, v_rm_diazlv, i, 1)
                  if (errstat /= hor3map_noerr) then
                     write(lp,*) trim(hor3map_errstr(errstat))
                     call xchalt('(ale_regrid_remap)')
                     stop '(ale_regrid_remap)'
                  endif
                  do iogrp = 1, nphy
                     if (acc_uvellvl(iogrp) /= 0) then
                        do kd = 1, ddm
                           phylvl(i,j,kd,acc_uvellvl(iogrp)) = &
                                phylvl(i,j,kd,acc_uvellvl(iogrp)) &
                              + v_rm_diazlv(kd) + ub(i,j,n)
                        enddo
                     endif
                  enddo

               endif

            enddo
         enddo

         do l = 1, isv(j)
            do i = max(1, ifv(j,l)), min(ii, ilv(j,l))

               ! Copy variables into 1D arrays. Rescale source interfaces so the
               ! pressure range of source and destination columns match.
               p_dst(1) = pv(i,j,1)
               do k = 1, kk
                  kn = k + nn
                  v_src(k) = v(i,j,kn)
                  p_dst(k+1) = p_dst(k) + dpv(i,j,kn)
               enddo
               q = min(util1(i,j-1), util1(i,j))/pv(i,j,kk+1)
               do k = 1, kk+1
                  p_src(k) = pv(i,j,k)*q
               enddo

               ! Prepare reconstruction with current interface pressures.
               errstat = prepare_reconstruction(rcgs, p_src, i, 1)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(ale_regrid_remap)')
                  stop '(ale_regrid_remap)'
               endif

               ! Prepare remapping to layer structure with regridded interface
               ! pressures.
               errstat = prepare_remapping(rcgs, rms, p_dst, i, 1)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(ale_regrid_remap)')
                  stop '(ale_regrid_remap)'
               endif

               ! Reconstruct v-component of velocity.
               errstat = reconstruct(rcgs, v_rcss, v_src, i, 1)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(ale_regrid_remap)')
                  stop '(ale_regrid_remap)'
               endif

               ! Remap v-component of velocity to regridded layers.
               errstat = remap(v_rcss, rms, v_rm, i, 1)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(ale_regrid_remap)')
                  stop '(ale_regrid_remap)'
               endif

               ! Update 3D array.
               do k = 1, kk
                  kn = k + nn
                  v(i,j,kn) = v_rm(k)
               enddo

               if (do_acc_vvellvl) then

                  ! Prepare remapping to destination z-levels.
                  p_dst_diazlv(1) = p_src(1)
                  q = 1._r8/vbath(i,j)
                  do kd = 2, ddm
                     p_dst_diazlv(kd) = min(depthslev_bnds(1,kd)*q, 1._r8) &
                                        *(p_src(kk+1) - p_src(1)) + p_src(1)
                  enddo
                  p_dst_diazlv(ddm+1) = p_src(kk+1)
                  errstat = prepare_remapping(rcgs, rms_diazlv, p_dst_diazlv, &
                                              i, 1)
                  if (errstat /= hor3map_noerr) then
                     write(lp,*) trim(hor3map_errstr(errstat))
                     call xchalt('(ale_regrid_remap)')
                     stop '(ale_regrid_remap)'
                  endif

                  ! Remap v-component of velocity to diagnostic z-levels and
                  ! accumulate.
                  errstat = remap(v_rcss, rms_diazlv, v_rm_diazlv, i, 1)
                  if (errstat /= hor3map_noerr) then
                     write(lp,*) trim(hor3map_errstr(errstat))
                     call xchalt('(ale_regrid_remap)')
                     stop '(ale_regrid_remap)'
                  endif
                  do iogrp = 1, nphy
                     if (acc_vvellvl(iogrp) /= 0) then
                        do kd = 1, ddm
                           phylvl(i,j,kd,acc_vvellvl(iogrp)) = &
                                phylvl(i,j,kd,acc_vvellvl(iogrp)) &
                              + v_rm_diazlv(kd) + vb(i,j,n)
                        enddo
                     endif
                  enddo

               endif

            enddo
         enddo

      enddo

      if (csdiag) then
         if (mnproc == 1) then
            write (lp,*) 'ale_regrid_remap:'
         endif
         call chksum(dp   (1-nbdy,1-nbdy,k1n), kk, halo_ps, 'dp'   )
         call chksum(temp (1-nbdy,1-nbdy,k1n), kk, halo_ps, 'temp' )
         call chksum(saln (1-nbdy,1-nbdy,k1n), kk, halo_ps, 'saln' )
         call chksum(sigma(1-nbdy,1-nbdy,k1n), kk, halo_ps, 'sigma')
         do nt = 1, ntr
            write(cnt, '(i2.2)') nt
            call chksum(trc(1-nbdy,1-nbdy,k1n,nt), kk, halo_ps, 'trc'//cnt)
         enddo
         if (ltedtp_opt == ltedtp_neutral) then
            call chksum(utflld(1-nbdy,1-nbdy,k1m), kk, halo_uv, 'utflld')
            call chksum(vtflld(1-nbdy,1-nbdy,k1m), kk, halo_vv, 'vtflld')
            call chksum(usflld(1-nbdy,1-nbdy,k1m), kk, halo_uv, 'usflld')
            call chksum(vsflld(1-nbdy,1-nbdy,k1m), kk, halo_vv, 'vsflld')
            call chksum(utflx (1-nbdy,1-nbdy,k1m), kk, halo_uv, 'utflx')
            call chksum(vtflx (1-nbdy,1-nbdy,k1m), kk, halo_vv, 'vtflx')
            call chksum(usflx (1-nbdy,1-nbdy,k1m), kk, halo_uv, 'usflx')
            call chksum(vsflx (1-nbdy,1-nbdy,k1m), kk, halo_vv, 'vsflx')
         endif
         call chksum(dpu(1-nbdy,1-nbdy,k1n), kk, halo_us, 'dpu')
         call chksum(dpv(1-nbdy,1-nbdy,k1n), kk, halo_vs, 'dpv')
         call chksum(u  (1-nbdy,1-nbdy,k1n), kk, halo_uv, 'u'  )
         call chksum(v  (1-nbdy,1-nbdy,k1n), kk, halo_vv, 'v'  )
      endif

   end subroutine ale_regrid_remap

   subroutine ale_remap_diazlv(m, n, mm, nn, k1m, k1n)
   ! ---------------------------------------------------------------------------
   ! Remap mid time-level variables to diagnostic z-levels.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      integer, parameter :: p_ord = 4

      real(r8), dimension(kdm+1,1-nbdy:idm+nbdy) :: p_src_js
      real(r8), dimension(p_ord+1,kdm,ntr_loc,1-nbdy:idm+nbdy) :: tpc_src_js
      real(r8), dimension(2,kdm,ntr_loc,1-nbdy:idm+nbdy) :: t_srcdi_js
      real(r8), dimension(kdm+1) :: pu_tmp, pv_tmp, p_src
      real(r8), dimension(kdm) :: v_src
      real(r8), dimension(ddm+1) :: p_dst_diazlv
      real(r8), dimension(ddm) :: v_rm_diazlv
      real(r8) :: q
      integer, dimension(1-nbdy:idm+nbdy) :: ksmx_js
      integer :: j, i, k, l, km, kd, iogrp, errstat
      logical :: do_acc_templvl  , do_acc_salnlvl, &
                 do_acc_uvellvl  , do_acc_vvellvl, &
                 do_acc_idlagelvl

      ! ------------------------------------------------------------------------
      ! Check if diagnostic variables should be remapped, either to be
      ! accumulated or instantaneously recorded.
      ! ------------------------------------------------------------------------

      do_acc_templvl = sum(acc_templvl  (1:nphy)*alarm_phy(1:nphy)) /= 0
      do_acc_salnlvl = sum(acc_salnlvl  (1:nphy)*alarm_phy(1:nphy)) /= 0
      do_acc_uvellvl = sum(acc_uvellvl  (1:nphy)*alarm_phy(1:nphy)) /= 0
      do_acc_vvellvl = sum(acc_vvellvl  (1:nphy)*alarm_phy(1:nphy)) /= 0
      do_acc_idlagelvl = &
         sum(acc_idlagelvl(1:nphy)*alarm_phy(1:nphy)) /= 0 .and. &
         use_TRC .and. use_IDLAGE

      ! ------------------------------------------------------------------------
      ! Remap tracers.
      ! ------------------------------------------------------------------------

      if (do_acc_templvl .or. do_acc_salnlvl .or. do_acc_idlagelvl) then

         do j = 1, jj

            ! Vertically reconstruct tracers.
            call reconstruct_trc_jslice(p_src_js, ksmx_js, &
                                        tpc_src_js, t_srcdi_js, &
                                        1, ii, j, 1, mm)

            ! Remap tracers to diagnostic z-levels.
            call remap_trc_diazlv_jslice(p_src_js, &
                                         1, ii, j, 1, &
                                         do_acc_templvl, do_acc_salnlvl, &
                                         do_acc_idlagelvl)

         enddo

      endif

      ! ------------------------------------------------------------------------
      ! Remap velocity.
      ! ------------------------------------------------------------------------

      if (.not. (do_acc_uvellvl .or. do_acc_vvellvl)) return

      do j = 1, jj

         if (do_acc_uvellvl) then

            do l = 1, isu(j)
               do i = max(1, ifu(j,l)), min(ii, ilu(j,l))

                  ! Copy variables into 1D arrays. Rescale source interfaces so
                  ! the pressure range of source and destination columns match.
                  pu_tmp(1) = pu(i,j,1)
                  do k = 1, kk
                     km = k + mm
                     v_src(k) = u(i,j,km)
                     pu_tmp(k+1) = pu_tmp(k) + dpu(i,j,km)
                  enddo
                  q = min(p(i-1,j,kk+1), p(i,j,kk+1))/pu_tmp(kk+1)
                  do k = 1, kk+1
                     p_src(k) = pu_tmp(k)*q
                  enddo

                  ! Prepare reconstruction with current interface pressures.
                  errstat = prepare_reconstruction(rcgs, p_src, i, 1)
                  if (errstat /= hor3map_noerr) then
                     write(lp,*) trim(hor3map_errstr(errstat))
                     call xchalt('(ale_regrid_remap)')
                     stop '(ale_regrid_remap)'
                  endif

                  ! Reconstruct u-component of velocity.
                  errstat = reconstruct(rcgs, v_rcss, v_src, i, 1)
                  if (errstat /= hor3map_noerr) then
                     write(lp,*) trim(hor3map_errstr(errstat))
                     call xchalt('(ale_regrid_remap)')
                     stop '(ale_regrid_remap)'
                  endif

                  ! Prepare remapping to destination z-levels.
                  p_dst_diazlv(1) = p_src(1)
                  q = 1._r8/ubath(i,j)
                  do kd = 2, ddm
                     p_dst_diazlv(kd) = min(depthslev_bnds(1,kd)*q, 1._r8) &
                                        *(p_src(kk+1) - p_src(1)) + p_src(1)
                  enddo
                  p_dst_diazlv(ddm+1) = p_src(kk+1)
                  errstat = prepare_remapping(rcgs, rms_diazlv, p_dst_diazlv, &
                                              i, 1)
                  if (errstat /= hor3map_noerr) then
                     write(lp,*) trim(hor3map_errstr(errstat))
                     call xchalt('(ale_regrid_remap)')
                     stop '(ale_regrid_remap)'
                  endif

                  ! Remap u-component of velocity to diagnostic z-levels and
                  ! accumulate.
                  errstat = remap(v_rcss, rms_diazlv, v_rm_diazlv, i, 1)
                  if (errstat /= hor3map_noerr) then
                     write(lp,*) trim(hor3map_errstr(errstat))
                     call xchalt('(ale_regrid_remap)')
                     stop '(ale_regrid_remap)'
                  endif
                  do iogrp = 1, nphy
                     if (acc_uvellvl(iogrp) /= 0) then
                        do kd = 1, ddm
                           phylvl(i,j,kd,acc_uvellvl(iogrp)) = &
                                phylvl(i,j,kd,acc_uvellvl(iogrp)) &
                              + v_rm_diazlv(kd) + ub(i,j,m)
                        enddo
                     endif
                  enddo

               enddo
            enddo

         endif

         if (do_acc_vvellvl) then

            do l = 1, isv(j)
               do i = max(1, ifv(j,l)), min(ii, ilv(j,l))

                  ! Copy variables into 1D arrays. Rescale source interfaces so
                  ! the pressure range of source and destination columns match.
                  pv_tmp(1) = pv(i,j,1)
                  do k = 1, kk
                     km = k + mm
                     v_src(k) = v(i,j,km)
                     pv_tmp(k+1) = pv_tmp(k) + dpv(i,j,km)
                  enddo
                  q = min(p(i,j-1,kk+1), p(i,j,kk+1))/pv_tmp(kk+1)
                  do k = 1, kk+1
                     p_src(k) = pv_tmp(k)*q
                  enddo

                  ! Prepare reconstruction with current interface pressures.
                  errstat = prepare_reconstruction(rcgs, p_src, i, 1)
                  if (errstat /= hor3map_noerr) then
                     write(lp,*) trim(hor3map_errstr(errstat))
                     call xchalt('(ale_regrid_remap)')
                     stop '(ale_regrid_remap)'
                  endif

                  ! Reconstruct v-component of velocity.
                  errstat = reconstruct(rcgs, v_rcss, v_src, i, 1)
                  if (errstat /= hor3map_noerr) then
                     write(lp,*) trim(hor3map_errstr(errstat))
                     call xchalt('(ale_regrid_remap)')
                     stop '(ale_regrid_remap)'
                  endif

                  ! Prepare remapping to destination z-levels.
                  p_dst_diazlv(1) = p_src(1)
                  q = 1._r8/vbath(i,j)
                  do kd = 2, ddm
                     p_dst_diazlv(kd) = min(depthslev_bnds(1,kd)*q, 1._r8) &
                                        *(p_src(kk+1) - p_src(1)) + p_src(1)
                  enddo
                  p_dst_diazlv(ddm+1) = p_src(kk+1)
                  errstat = prepare_remapping(rcgs, rms_diazlv, p_dst_diazlv, &
                                              i, 1)
                  if (errstat /= hor3map_noerr) then
                     write(lp,*) trim(hor3map_errstr(errstat))
                     call xchalt('(ale_regrid_remap)')
                     stop '(ale_regrid_remap)'
                  endif

                  ! Remap v-component of velocity to diagnostic z-levels and
                  ! accumulate.
                  errstat = remap(v_rcss, rms_diazlv, v_rm_diazlv, i, 1)
                  if (errstat /= hor3map_noerr) then
                     write(lp,*) trim(hor3map_errstr(errstat))
                     call xchalt('(ale_regrid_remap)')
                     stop '(ale_regrid_remap)'
                  endif
                  do iogrp = 1, nphy
                     if (acc_vvellvl(iogrp) /= 0) then
                        do kd = 1, ddm
                           phylvl(i,j,kd,acc_vvellvl(iogrp)) = &
                                phylvl(i,j,kd,acc_vvellvl(iogrp)) &
                              + v_rm_diazlv(kd) + vb(i,j,m)
                        enddo
                     endif
                  enddo

               enddo
            enddo

         endif

      enddo

   end subroutine ale_remap_diazlv

end module mod_ale_regrid_remap
