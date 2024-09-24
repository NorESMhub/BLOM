! ------------------------------------------------------------------------------
! Copyright (C) 2021-2024 Mats Bentsen, Mehmet Ilicak
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
   use mod_constants, only: g, epsilp, onem
   use mod_time,      only: delt1
   use mod_xc
   use mod_grid,      only: scuy, scvx, scp2, scuxi, scvyi, scp2i
   use mod_vcoord,    only: vcoord_tag, vcoord_isopyc_bulkml, &
                            vcoord_cntiso_hybrid, vcoord_plevel, &
                            sigmar, plevel
   use mod_eos,       only: sig, dsigdt, dsigds
   use mod_state,     only: u, v, dp, dpu, dpv, temp, saln, sigma, p, pu, pv
   use mod_hor3map,   only: recon_grd_struct, recon_src_struct, remap_struct, &
                            hor3map_plm, hor3map_ppm, hor3map_pqm, &
                            hor3map_monotonic, hor3map_non_oscillatory, &
                            hor3map_non_oscillatory_posdef, &
                            initialize_rcgs, initialize_rcss, initialize_rms, &
                            prepare_reconstruction, reconstruct, &
                            extract_polycoeff, regrid, &
                            prepare_remapping, remap, &
                            hor3map_noerr, hor3map_errstr
   use mod_diffusion, only: ltedtp_opt, ltedtp_neutral, difiso, difmxp
   use mod_ndiff,     only: ndiff_prep_jslice, ndiff_uflx_jslice, &
                            ndiff_vflx_jslice, ndiff_update_trc_jslice
   use mod_checksum,  only: csdiag, chksummsk
   use mod_tracers,   only: ntr, trc

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
        smooth_diff_max        = 50000._r8
   integer :: &
        upper_bndr_ord = 6, &
        lower_bndr_ord = 4, &
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
   type(remap_struct) :: rms

   public :: regrid_method_tag, regrid_method_direct, &
             readnml_ale_regrid_remap, init_ale_regrid_remap, &
             ale_regrid_remap

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
      real(r8) :: beta, sdpsum, smean, dpmin_int, pku, pku_test, pmin, dpt, &
                  pt, ptu1, ptl1, ptu2, ptl2, w1, x
      integer :: l, i, k, kn, ks, ke, kl, ku, errstat
      logical :: thin_layers, layer_added

      ! Minimum potential density difference with respect to pressure for
      ! potential density to be used in regridding.
      beta = bfsq_min/(g*g)

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
            dpmin_int = min(plevel(2) - plevel(1), dpmin_interior)
            ks = max(2, ks)
            ke = min(kk, ke)
            k = ks
            do while (k <= ke)
               if (p_dst(k+1,i) - p_dst(k,i) < dpmin_int) then
                  if (k == ke) then
                     p_dst(k,i) = p_dst(ke+1,i)
                  else
                     ku = k
                     kl = k + 1
                     pku = .5_r8*(p_dst(kl,i) + p_dst(ku,i) - dpmin_int)
                     do
                        layer_added = .false.
                        kl = kl + 1
                        pku_test = ((pku - dpmin_int)*(kl - ku) + p_dst(kl,i)) &
                                   /(kl - ku + 1)
                        if (pku_test + (kl - ku)*dpmin_int > p_dst(kl,i)) then
                           if (kl == ke + 1) exit
                           pku = pku_test
                           layer_added = .true.
                        else
                           kl = kl - 1
                        endif
                        ku = ku - 1
                        pku_test = ((pku - dpmin_int)*(kl - ku) + p_dst(ku,i)) &
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
                                            p_dst(k -1,i) + dpmin_int)
                        enddo
                        do k = kl+1, ke
                           p_dst(k,i) = &
                              min(p_dst(ke+1,i), &
                              max(p_dst(k,i), p_dst(1,i) + dpmin_int*(k - 1)))
                        enddo
                     elseif (kl == ke + 1) then
                        do k = ku, kl
                           p_dst(k,i) = p_dst(ke+1,i)
                        enddo
                     else
                        p_dst(ku,i) = pku
                        do k = ku+1, kl
                           p_dst(k,i) = p_dst(k-1,i) + dpmin_int
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
            dpt = plevel(2) - plevel(1)
            do k = 2, ke
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
               p_dst(k,i) = min(p_dst(ke+1,i), &
                                max(p_dst(k-1,i) + dpmin_int, pt))
            enddo

         enddo
      enddo

   end subroutine regrid_cntiso_hybrid_direct_jslice

   subroutine regrid_cntiso_hybrid_nudge_jslice( &
                 p_src, ksmx, tpc_src, t_srcdi, p_dst, stab_fac, ilb, iub, j)
   ! ---------------------------------------------------------------------------
   ! For vcoord == 'cntiso_hybrid' and regrid_method = 'nudge', nudge the
   ! interface pressures to reduce the deviation from interface target potential
   ! densities, while maintaining minimum layer thicknesses towards the surface.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:,1-nbdy:), intent(in) :: p_src
      integer, dimension(1-nbdy:), intent(in) :: ksmx
      real(r8), dimension(:,:,:,1-nbdy:), intent(in) :: tpc_src, t_srcdi
      real(r8), dimension(:,1-nbdy:), intent(out) :: p_dst, stab_fac
      integer, intent(in) :: ilb, iub, j

      integer, parameter :: &
           it = 1, &
           is = 2

      real(r8), dimension(2,kdm) :: sig_srcdi
      integer, dimension(1-nbdy:idm+nbdy) :: kdmx

      real(r8), dimension(kdm+1) :: sig_trg, sig_pmin
      real(r8), dimension(kdm) :: dsig_trg, pmin
      real(r8) :: sig_max, ckt, sig_up, sig_lo, dk, dki, &
                  dsigdx_up, dsigdx_lo, x, xi, si, t, nudge_fac, &
                  dsig, dsigdx, dp_up, dp_lo, sig_intrp
      integer :: l, i, k, kt, kl, ktzmin, ktzmax

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
            p_dst(1,i) = pmin(1)

            stab_fac(1,i) = 0._r8
            nudge_fac = delt1/regrid_nudge_ts

            ! Find the index of the first interface with potential density at
            ! minimum interface pressure smaller than the reference potential
            ! density of this transition interface. A layer range above and
            ! below the transition interface may be specified, making a
            ! transition zone where interface reference potential densities are
            ! adjusted to achieve a more gradual change from pressure level to
            ! isopycnic interfaces.
            sig_pmin(1) = sig_srcdi(1,1)
            kt = 2
            kl = 1
            do while (kt <= kdmx(i))
               do while (p_src(kl+1,i) < pmin(kt))
                  kl = kl + 1
               enddo
               sig_pmin(kt) = ( (p_src(kl+1,i) - pmin(kt))*sig_srcdi(1,kl) &
                              + (pmin(kt) - p_src(kl,i))*sig_srcdi(2,kl)) &
                              /(p_src(kl+1,i) - p_src(kl,i))
               if (sig_trg(kt) > sig_pmin(kt)) then
                  ktzmin = max(3, kt - dktzu)
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
               stab_fac(kt,i) = 0._r8
               kt = kt + 1
            enddo

            ! Starting at the transition interface, nudge the interface
            ! pressures to reduce the deviation from the interface reference
            ! potential density.

            do k = kt, kk+1
               p_dst(k,i) = p_src(kk+1,i)
               stab_fac(k,i) = 1._r8
            enddo

            do k = kt, min(ksmx(i), kdmx(i))
               if      (sig_trg(k) < sig_srcdi(2,k-1) .and. &
                        sig_trg(k) < sig_srcdi(1,k  )) then
                  dsig = sig_trg(k) - sig_srcdi(2,k-1)
                  dsigdx = dsigdt(t_srcdi(2,k-1,it,i), t_srcdi(2,k-1,is,i)) &
                           *dpeval1(tpc_src(:,k-1,it,i)) &
                         + dsigds(t_srcdi(2,k-1,it,i), t_srcdi(2,k-1,is,i)) &
                           *dpeval1(tpc_src(:,k-1,is,i))
                  stab_fac(k,i) = dsigdx/dsig_trg(k-1)
                  dsigdx = dsig_trg(k-1)*max(stab_fac(k,i), stab_fac_limit)
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
                  stab_fac(k,i) = dsigdx/dsig_trg(k)
                  dsigdx = dsig_trg(k)*max(stab_fac(k,i), stab_fac_limit)
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
                     stab_fac(k,i) = dsigdx/dsig_trg(k-1)
                     dsigdx = dsig_trg(k-1)*max(stab_fac(k,i), stab_fac_limit)
                     p_dst(k,i) = p_src(k,i) &
                                + max(- .5_r8, dsig*nudge_fac/dsigdx) &
                                  *(p_src(k,i) - p_src(k-1,i))
                  else
                     dsigdx = dsigdx_lo + 2._r8*(sig_srcdi(1,k  ) - sig_intrp)
                     stab_fac(k,i) = dsigdx/dsig_trg(k)
                     dsigdx = dsig_trg(k)*max(stab_fac(k,i), stab_fac_limit)
                     p_dst(k,i) = p_src(k,i) &
                                + min(.5_r8, dsig*nudge_fac/dsigdx) &
                                  *(p_src(k+1,i) - p_src(k,i))
                  endif
               endif
               p_dst(k,i) = min(max(p_dst(k,i), pmin(k), &
                                    p_dst(k-1,i) + dpmin_interior), &
                                p_src(kk+1,i))
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
                  stab_fac(k,i) = dsigdx/dsig_trg(ksmx(i)-1)
                  dsigdx = dsig_trg(ksmx(i)-1) &
                           *max(stab_fac(k,i), stab_fac_limit)
                  p_dst(k,i) = p_src(kk+1,i) &
                             + max(- .5_r8, dsig*nudge_fac/dsigdx) &
                               *(p_src(kk+1,i) - p_src(ksmx(i),i))
                  p_dst(k,i) = min(max(p_dst(k,i), pmin(k), &
                                       p_dst(k-1,i) + dpmin_interior), &
                                   p_src(kk+1,i))
               endif
            enddo

         enddo
      enddo

   end subroutine regrid_cntiso_hybrid_nudge_jslice

   subroutine regrid_jslice(p_src, ksmx, tpc_src, t_srcdi, p_dst, stab_fac, &
                            ilb, iub, j, js, nn)
   ! ---------------------------------------------------------------------------
   ! Carry out regridding layer interfaces.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:,1-nbdy:), intent(in) :: p_src
      integer, dimension(1-nbdy:), intent(in) :: ksmx
      real(r8), dimension(:,:,:,1-nbdy:), intent(in) :: tpc_src, t_srcdi
      real(r8), dimension(:,1-nbdy:), intent(out) :: p_dst, stab_fac
      integer, intent(in) :: ilb, iub, j, js, nn

      if (vcoord_tag == vcoord_plevel) then
         call regrid_plevel_jslice(p_src, p_dst, ilb, iub, j)
      else
         if (regrid_method_tag == regrid_method_direct) then
            call regrid_cntiso_hybrid_direct_jslice(p_src, p_dst, &
                                                    ilb, iub, j, js, nn)
         else
            call regrid_cntiso_hybrid_nudge_jslice(p_src, ksmx, tpc_src, &
                                                   t_srcdi, p_dst, stab_fac, &
                                                   ilb, iub, j)
         endif
      endif

   end subroutine regrid_jslice

   subroutine regrid_smooth_jslice(p_dst_js, stab_fac_js, smtflxconv_js, &
                                   ilb, iub, j, js2, js3)
   ! ---------------------------------------------------------------------------
   ! For vcoord == 'cntiso_hybrid' and regrid_method == 'nudge', apply lateral
   ! smoothing of the regridded interfaces when a vertical stability factor is
   ! below a specified threshold.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:,1-nbdy:,:), intent(inout) :: &
         p_dst_js, stab_fac_js, smtflxconv_js
      integer, intent(in) :: ilb, iub, j, js2, js3

      real(r8) :: cdiff, difmx, flxhi, flxlo, flx, q, sdiff
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
               q = .5_r8*( max(0._r8, min(stab_fac_limit, &
                                          stab_fac_js(k,i-1,js3))) &
                         + max(0._r8, min(stab_fac_limit, &
                                          stab_fac_js(k,i  ,js3))))
               sdiff = min((stab_fac_limit - q)*smooth_diff_max &
                           /stab_fac_limit, difmx)
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
               q = .5_r8*( max(0._r8, min(stab_fac_limit, &
                                          stab_fac_js(k,i,js2))) &
                         + max(0._r8, min(stab_fac_limit, &
                                          stab_fac_js(k,i,js3))))
               sdiff = min((stab_fac_limit - q)*smooth_diff_max &
                           /stab_fac_limit, difmx)
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
         regrid_method, regrid_nudge_ts, stab_fac_limit, smooth_diff_max, &
         dktzu, dktzl

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
         call xcbcst(regrid_nudge_ts)
         call xcbcst(stab_fac_limit)
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
         write (lp,*) '  regrid_nudge_ts =        ', regrid_nudge_ts
         write (lp,*) '  stab_fac_limit =         ', stab_fac_limit
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

      ! Change units from [m] to [g cm-1 s-2] of depth interval variables.
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

   end subroutine init_ale_regrid_remap

   subroutine ale_regrid_remap(m, n, mm, nn, k1m, k1n)
   ! ---------------------------------------------------------------------------
   ! Regrid, remap and carry out additional operations that makes use of the
   ! reconstructed vertical profiles, such as neutral diffusion.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      integer, parameter :: p_ord = 4

      real(r8), dimension(kdm+1,1-nbdy:idm+nbdy,3) :: &
         p_src_js, p_dst_js, stab_fac_js, smtflxconv_js
      real(r8), dimension(p_ord+1,kdm,ntr_loc,1-nbdy:idm+nbdy,3) :: tpc_src_js
      real(r8), dimension(2,kdm,ntr_loc,1-nbdy:idm+nbdy,3) :: t_srcdi_js
      real(r8), dimension(2,kdm,1-nbdy:idm+nbdy,3) :: &
         p_srcdi_js, drhodt_srcdi_js, drhods_srcdi_js
      real(r8), dimension(kdm,ntr_loc,1-nbdy:idm+nbdy,3) :: flxconv_js
      real(r8), dimension(kdm,ntr_loc,1-nbdy:idm+nbdy) :: trc_rm
      real(r8), dimension(kdm+1) :: p_1d, p_dst_1d
      real(r8), dimension(kdm) :: u_1d, v_1d
      real(r8) :: q
      integer, dimension(1-nbdy:idm+nbdy,3) :: ksmx_js, kdmx_js
      integer :: ilb1, ilb2, ilb3, iub1, iub2, iub3, jofs2, jofs3, &
                 jlb_regrid_smooth, jlb_ndiff_prep, jlb_ndiff_uflx, &
                 jlb_ndiff_vflx, jlb_ndiff_update_trc, &
                 js1, js2, js3, j, nt, i, k, l, kn, errstat
      logical :: do_regrid_smooth, do_ndiff

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
         call xctilr(dp   (1-nbdy,1-nbdy,k1n), 1, kk, jofs3, jofs3, halo_ps)
         call xctilr(temp (1-nbdy,1-nbdy,k1n), 1, kk, jofs3, jofs3, halo_ps)
         call xctilr(saln (1-nbdy,1-nbdy,k1n), 1, kk, jofs3, jofs3, halo_ps)
         call xctilr(sigma(1-nbdy,1-nbdy,k1n), 1, kk, jofs3, jofs3, halo_ps)
      endif
      if (do_ndiff) then
         do nt = 1, ntr
            call xctilr(trc(1-nbdy,1-nbdy,k1n,nt), 1, kk, 1, 1, halo_ps)
         enddo
         call xctilr(difiso, 1,kk, 1,1, halo_ps)
      end if

      ! Inital j-slice indices.
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
         call regrid_jslice         (p_src_js(:,:,js3), ksmx_js(:,js3), &
                                     tpc_src_js(:,:,:,:,js3), &
                                     t_srcdi_js(:,:,:,:,js3), &
                                     p_dst_js(:,:,js3), stab_fac_js(:,:,js3), &
                                     ilb3, iub3, j+jofs3, js3, nn)

         ! If requested, apply lateral smoothing of the interfaces after
         ! regridding.
         if (j >= jlb_regrid_smooth) &
            call regrid_smooth_jslice   (p_dst_js, stab_fac_js, smtflxconv_js, &
                                         ilb2, iub2, j+jofs2, js2, js3)

         ! If requested, prepare neutral diffusion.
         if (j >= jlb_ndiff_prep) &
            call ndiff_prep_jslice      (p_src_js, ksmx_js, &
                                         tpc_src_js, t_srcdi_js, &
                                         p_dst_js, kdmx_js, p_srcdi_js, &
                                         drhodt_srcdi_js, drhods_srcdi_js, &
                                         flxconv_js, &
                                         ilb2, iub2, j+jofs2, js2, mm)

         ! If requested, compute the contribution of u-component fluxes to the
         ! flux convergence of neutral diffusion.
         if (j >= jlb_ndiff_uflx) &
            call ndiff_uflx_jslice      (ksmx_js, tpc_src_js, t_srcdi_js, &
                                         p_dst_js, kdmx_js, p_srcdi_js, &
                                         drhodt_srcdi_js, drhods_srcdi_js, &
                                         flxconv_js, &
                                         ntr_loc, ilb1, iub2, j, js1, mm, nn)

         ! If requested, compute the contribution of v-component fluxes to the
         ! flux convergence of neutral diffusion.
         if (j >= jlb_ndiff_vflx) &
            call ndiff_vflx_jslice      (ksmx_js, tpc_src_js, t_srcdi_js, &
                                         p_dst_js, kdmx_js, p_srcdi_js, &
                                         drhodt_srcdi_js, drhods_srcdi_js, &
                                         flxconv_js, &
                                         ntr_loc, ilb1, iub1, j+jofs2, &
                                         js1, js2, mm, nn)

         ! Remap tracers to the regridded layers.
         if (j >= 1) &
            call remap_trc_jslice       (p_dst_js(:,:,js1), trc_rm, &
                                         ilb1, iub1, j, js1)

         ! If requested, update the tracers by applying the neutral diffusion
         ! flux convergence.
         if (j >= jlb_ndiff_update_trc) &
            call ndiff_update_trc_jslice(p_dst_js, flxconv_js, trc_rm, &
                                         ntr_loc, ilb1, iub1, j, js1)

         ! Copy from the j-slice array to the full tracer array.
         if (j >= 1) &
            call copy_jslice_to_3d      (p_dst_js(:,:,js1), trc_rm, &
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
               p_dst_1d(1) = pu(i,j,1)
               do k = 1, kk
                  kn = k + nn
                  u_1d(k) = u(i,j,kn)
                  p_dst_1d(k+1) = p_dst_1d(k) + dpu(i,j,kn)
               enddo
               q = p_dst_1d(kk+1)/pu(i,j,kk+1)
               do k = 1, kk+1
                  p_1d(k) = pu(i,j,k)*q
               enddo

               ! Prepare reconstruction with current interface pressures.
               errstat = prepare_reconstruction(rcgs, p_1d, i, 1)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(ale_regrid_remap)')
                  stop '(ale_regrid_remap)'
               endif

               ! Prepare remapping to layer structure with regridded interface
               ! pressures.
               errstat = prepare_remapping(rcgs, rms, p_dst_1d, i, 1)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(ale_regrid_remap)')
                  stop '(ale_regrid_remap)'
               endif

               ! Reconstruct and remap u-component of velocity.
               errstat = reconstruct(rcgs, v_rcss, u_1d, i, 1)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(ale_regrid_remap)')
                  stop '(ale_regrid_remap)'
               endif
               errstat = remap(v_rcss, rms, u_1d, i, 1)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(ale_regrid_remap)')
                  stop '(ale_regrid_remap)'
               endif

               ! Update 3D arrays
               do k = 1, kk
                  kn = k + nn
                  u(i,j,kn) = u_1d(k)
               enddo

            enddo
         enddo

         do l = 1, isv(j)
            do i = max(1, ifv(j,l)), min(ii, ilv(j,l))

               ! Copy variables into 1D arrays. Rescale source interfaces so the
               ! pressure range of source and destination columns match.
               p_dst_1d(1) = pv(i,j,1)
               do k = 1, kk
                  kn = k + nn
                  v_1d(k) = v(i,j,kn)
                  p_dst_1d(k+1) = p_dst_1d(k) + dpv(i,j,kn)
               enddo
               q = p_dst_1d(kk+1)/pv(i,j,kk+1)
               do k = 1, kk+1
                  p_1d(k) = pv(i,j,k)*q
               enddo

               ! Prepare reconstruction with current interface pressures.
               errstat = prepare_reconstruction(rcgs, p_1d, i, 1)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(ale_regrid_remap)')
                  stop '(ale_regrid_remap)'
               endif

               ! Prepare remapping to layer structure with regridded interface
               ! pressures.
               errstat = prepare_remapping(rcgs, rms, p_dst_1d, i, 1)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(ale_regrid_remap)')
                  stop '(ale_regrid_remap)'
               endif

               ! Reconstruct and remap v-component of velocity.
               errstat = reconstruct(rcgs, v_rcss, v_1d, i, 1)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(ale_regrid_remap)')
                  stop '(ale_regrid_remap)'
               endif
               errstat = remap(v_rcss, rms, v_1d, i, 1)
               if (errstat /= hor3map_noerr) then
                  write(lp,*) trim(hor3map_errstr(errstat))
                  call xchalt('(ale_regrid_remap)')
                  stop '(ale_regrid_remap)'
               endif

               ! Update 3D arrays
               do k = 1, kk
                  kn = k + nn
                  v(i,j,kn) = v_1d(k)
               enddo

            enddo
         enddo

      enddo

      if (csdiag) then
         if (mnproc == 1) then
            write (lp,*) 'ale_regrid_remap:'
         endif
         call chksummsk(dp   (1-nbdy,1-nbdy,k1n), ip, kk, 'dp')
         call chksummsk(temp (1-nbdy,1-nbdy,k1n), ip, kk, 'temp')
         call chksummsk(saln (1-nbdy,1-nbdy,k1n), ip, kk, 'saln')
         call chksummsk(sigma(1-nbdy,1-nbdy,k1n), ip, kk, 'sigma')
         do nt = 1, ntr
            call chksummsk(trc(1-nbdy,1-nbdy,k1n,nt), ip, kk, 'trc')
         enddo
         call chksummsk(dpu(1-nbdy,1-nbdy,k1n), iu, kk, 'dpu')
         call chksummsk(dpv(1-nbdy,1-nbdy,k1n), iv, kk, 'dpv')
         call chksummsk(u  (1-nbdy,1-nbdy,k1n), iu, kk, 'u')
         call chksummsk(v  (1-nbdy,1-nbdy,k1n), iv, kk, 'v')
      endif

   end subroutine ale_regrid_remap

end module mod_ale_regrid_remap
