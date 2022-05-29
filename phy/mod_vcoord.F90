! ------------------------------------------------------------------------------
! Copyright (C) 2021-2022 Mats Bentsen
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

module mod_vcoord
! ------------------------------------------------------------------------------
! This module contains parameter, variables and procedures related to the
! vertical coordinate.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_config, only: inst_suffix
   use mod_constants, only: g, epsil, spval, onem
   use mod_xc
   use mod_eos, only: sig, dsigdt, dsigds
   use mod_state, only: u, v, dp, dpu, dpv, temp, saln, sigma, p, pu, pv
   use mod_hor3map, only: recon_grd_struct, recon_src_struct, remap_struct, &
                          hor3map_plm, hor3map_ppm, hor3map_pqm, &
                          hor3map_monotonic, hor3map_non_oscillatory, &
                          hor3map_non_oscillatory_posdef, &
                          initialize_rcgs, initialize_rcss, initialize_rms, &
                          prepare_reconstruction, reconstruct, &
                          extract_polycoeff, regrid2, &
                          prepare_remapping, remap, &
                          hor3map_noerr, hor3map_errstr
   use mod_diffusion, only : ltedtp_opt, ltedtp_neutral, difiso
   use mod_ndiff, only: ndiff_prep_jslice, ndiff_uflx_jslice, &
                        ndiff_vflx_jslice, ndiff_update_trc_jslice
   use mod_checksum, only: csdiag, chksummsk
#ifdef TRC
   use mod_tracers, only: ntr, trc
#endif

   implicit none

   private

   ! Options with default values, modifiable by namelist.
   character(len = 80) :: &
      vcoord_type = 'isopyc_bulkml', &
      reconstruction_method = 'ppm', &
      density_limiting = 'monotonic', &
      tracer_limiting = 'monotonic', &
      velocity_limiting = 'monotonic'
   logical :: &
      density_pc_upper_bndr = .false., &
      density_pc_lower_bndr = .false., &
      tracer_pc_upper_bndr = .true., &
      tracer_pc_lower_bndr = .false., &
      velocity_pc_upper_bndr = .true., &
      velocity_pc_lower_bndr = .false.
   real(r8) :: &
      dpmin_surface          = 1.5_r8, &
      dpmin_inflation_factor = 1._r8, &
      dpmin_interior         = .1_r8, &
      regrid_nudge_factor    = .1_r8

   ! Options derived from string options.
   integer :: &
      vcoord_type_tag, &
      reconstruction_method_tag, &
      density_limiting_tag, &
      tracer_limiting_tag, &
      velocity_limiting_tag

   ! Parameters:
   integer, parameter :: &
      isopyc_bulkml = 1, &        ! Vertical coordinate type: bulk surface mixed
                                  ! layer with isopycnic layers below.
      cntiso_hybrid = 2, &        ! Vertical coordinate type: Hybrid coordinate
                                  ! with pressure coordinates towards the
                                  ! surface and continuous isopycnal below.
#ifdef TRC
      ntr_loc = ntr + 2           ! Local number of tracers where temperature
                                  ! and salinity is added to the ntr parameter.
#else
      ntr_loc = 2                 ! Local number of tracers consisting of
                                  ! temperature and salinity.
#endif
   real(r8), parameter :: &
      bfsq_min      = 1.e-7_r8, & ! Minimum buoyancy frequency squared in
                                  ! monotonized potential density to be used in
                                  ! regridding [s-2].
      regrid_mval   = - 1.e33_r8  ! Missing value for regridding.


   type(recon_grd_struct) :: rcgs
   type(recon_src_struct) :: d_rcss, v_rcss
   type(recon_src_struct) , dimension(ntr_loc) :: trc_rcss
   type(remap_struct) :: rms

   real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: &
      sigmar     ! Reference potential density [g cm-3].

   public :: vcoord_type_tag, isopyc_bulkml, cntiso_hybrid, sigmar, &
             readnml_vcoord, inivar_vcoord, cntiso_hybrid_regrid_direct_remap, &
             cntiso_hybrid_regrid_remap, remap_velocity

contains

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

   subroutine prep_recon_jslice(p_src, i_lb, i_ub, j, j_rs, nn)
   ! ---------------------------------------------------------------------------
   ! Prepare vertical layer reconstruction along a j-slice of the model data.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:,1-nbdy:), intent(out) :: p_src
      integer, intent(in) :: i_lb, i_ub, j, j_rs, nn

      integer :: l, i, k, errstat

      do l = 1, isp(j)
      do i = max(i_lb, ifp(j,l)), min(i_ub, ilp(j,l))

         p_src(1,i) = p(i,j,1)
         do k = 1, kk
            p_src(k+1,i) = p_src(k,i) + dp(i,j,k+nn)
         enddo

         errstat = prepare_reconstruction(rcgs, p_src(:,i), i, j_rs)
         if (errstat /= hor3map_noerr) then
            write(lp,*) trim(hor3map_errstr(errstat))
            call xchalt('(prep_recon_jslice)')
                   stop '(prep_recon_jslice)'
         endif

      enddo
      enddo

   end subroutine prep_recon_jslice

   subroutine recon_trc_jslice(i_lb, i_ub, j, j_rs, nn)
   ! ---------------------------------------------------------------------------
   ! Vertically reconstruct temperature, salinity and additional tracers along a
   ! j-slice of the model data.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: i_lb, i_ub, j, j_rs, nn

      real(r8), dimension(kdm,ntr_loc) :: trc_1d
      integer :: l, i, k, kn, nt, errstat

      do l = 1, isp(j)
      do i = max(i_lb, ifp(j,l)), min(i_ub, ilp(j,l))

         ! Copy variables into 1D arrays.
         do k = 1, kk
            kn = k + nn
            trc_1d(k,1) = temp(i,j,kn)
            trc_1d(k,2) = saln(i,j,kn)
#ifdef TRC
            do nt = 1, ntr
              trc_1d(k,nt+2) = trc(i,j,kn,nt)
            enddo
#endif
         enddo

         ! Reconstruct tracers.
         do nt = 1, ntr_loc
            errstat = reconstruct(rcgs, trc_rcss(nt), trc_1d(:,nt), i, j_rs)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(recon_trc_jslice)')
                      stop '(recon_trc_jslice)'
            endif
         enddo

      enddo
      enddo

   end subroutine recon_trc_jslice

   subroutine remap_trc_jslice(p_dst, trc_rm, i_lb, i_ub, j, j_rs)

      real(r8), dimension(:,1-nbdy:), intent(in) :: p_dst
      real(r8), dimension(:,:,1-nbdy:), intent(out) :: trc_rm
      integer, intent(in) :: i_lb, i_ub, j, j_rs

      integer :: l, i, nt, errstat

      do l = 1, isp(j)
      do i = max(i_lb, ifp(j,l)), min(i_ub, ilp(j,l))

         ! Prepare remapping to target layers.
         errstat = prepare_remapping(rcgs, rms, p_dst(:,i), i, j_rs)
         if (errstat /= hor3map_noerr) then
            write(lp,*) trim(hor3map_errstr(errstat))
            write(lp,*) 'i, j:', i + i0, j + j0
            do nt = 1,kk
               write(lp,*) nt, p_dst(nt+1,i), p_dst(nt+1,i) - p_dst(nt,i)
            enddo
            call xchalt('(remap_trc_jslice)')
                   stop '(remap_trc_jslice)'
         endif

         ! Remap tracers.
         do nt = 1, ntr_loc
            errstat = remap(trc_rcss(nt), rms, trc_rm(:,nt,i), i, j_rs)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(remap_trc_jslice)')
                      stop '(remap_trc_jslice)'
            endif
         enddo

      enddo
      enddo

   end subroutine remap_trc_jslice

   subroutine cntiso_regrid_direct_jslice(p_src, p_dst, i_lb, i_ub, j, j_rs, nn)

      real(r8), dimension(:,1-nbdy:), intent(in) :: p_src
      real(r8), dimension(:,1-nbdy:), intent(out) :: p_dst
      integer, intent(in) :: i_lb, i_ub, j, j_rs, nn

      real(r8), dimension(kdm+1) :: sigmar_1d
      real(r8), dimension(kdm) :: sigma_1d
      real(r8) :: beta, sdpsum, smean, dpmin_max, dpmin_int, dpmin_sfc, &
                  pku, pku_test, pmin, dpt, pt, ptu1, ptl1, ptu2, ptl2, w1, x
      integer :: l, i, k, kn, ks, ke, kl, ku, errstat
      logical :: thin_layers, layer_added

      ! Minimum potential density difference with respect to pressure for
      ! potential density to be used in regridding.
      beta = bfsq_min/(g*g)

      do l = 1, isp(j)
      do i = max(i_lb, ifp(j,l)), min(i_ub, ilp(j,l))

         ! Copy variables into 1D arrays.
         do k = 1, kk
            kn = k + nn
            sigma_1d(k) = sigma(i,j,kn)
            sigmar_1d(k) = sigmar(i,j,k)
         enddo
         sigmar_1d(kk+1) = sigmar_1d(kk)

         ! Make sure potential density to be used in regridding is
         ! monotonically increasing with depth.
         kl = kk
         ku = kl - 1
         do while (ku > 0)
            thin_layers = p_src(kl+1,i) - p_src(ku,i) < epsil
            if (thin_layers .or. &
                sigma_1d(kl) - sigma_1d(ku) &
                < .5_r8*beta*(p_src(kl+1,i) - p_src(ku,i))) then
               sdpsum = sigma_1d(ku)*(p_src(ku+1,i) - p_src(ku,i)) &
                      + sigma_1d(kl)*(p_src(kl+1,i) - p_src(kl,i))
               if (.not. thin_layers) &
                  smean = sdpsum/(p_src(kl+1,i) - p_src(ku,i))
               do
                  layer_added = .false.
                  if (ku > 1) then
                     if (thin_layers) then
                        ku = ku - 1
                        sdpsum = sdpsum &
                               + sigma_1d(ku)*(p_src(ku+1,i) - p_src(ku,i))
                        thin_layers = p_src(kl+1,i) - p_src(ku,i) < epsil
                        if (.not. thin_layers) &
                           smean = sdpsum/(p_src(kl+1,i) - p_src(ku,i))
                        layer_added = .true.
                     else
                        if (smean - sigma_1d(ku-1) &
                            < .5_r8*beta*(p_src(kl+1,i) - p_src(ku-1,i))) then
                           ku = ku - 1
                           sdpsum = sdpsum &
                                  + sigma_1d(ku)*(p_src(ku+1,i) - p_src(ku,i))
                           smean = sdpsum/(p_src(kl+1,i) - p_src(ku,i))
                           layer_added = .true.
                        endif
                     endif
                  endif
                  if (kl < kk) then
                     if (thin_layers) then
                        kl = kl + 1
                        sdpsum = sdpsum &
                               + sigma_1d(kl)*(p_src(kl+1,i) - p_src(kl,i))
                        thin_layers = p_src(kl+1,i) - p_src(ku,i) < epsil
                        if (.not. thin_layers) &
                           smean = sdpsum/(p_src(kl+1,i) - p_src(ku,i))
                        layer_added = .true.
                     else
                        if (sigma_1d(kl+1) - smean &
                            < .5_r8*beta*(p_src(kl+2,i) - p_src(ku,i))) then
                           kl = kl + 1
                           sdpsum = sdpsum &
                                  + sigma_1d(kl)*(p_src(kl+1,i) - p_src(kl,i))
                           smean = sdpsum/(p_src(kl+1,i) - p_src(ku,i))
                           layer_added = .true.
                        endif
                     endif
                  endif
                  if (.not. layer_added) exit
               enddo
               do k = ku, kl
                  sigma_1d(k) = smean &
                              + .5_r8*beta*( p_src(k ,i) + p_src(k +1,i) &
                                           - p_src(ku,i) - p_src(kl+1,i))
               enddo
            endif
            kl = ku
            ku = kl - 1
         enddo

         ! Monotonically reconstruct potential density.
         errstat = reconstruct(rcgs, d_rcss, sigma_1d, i, j_rs)
         if (errstat /= hor3map_noerr) then
            write(lp,*) trim(hor3map_errstr(errstat))
            call xchalt('(cntiso_regrid_direct_jslice)')
                   stop '(cntiso_regrid_direct_jslice)'
         endif

         ! On the basis of the reconstructed potential density, regrid
         ! interface pressures so interface potential densities match target
         ! values.
         errstat = regrid2(d_rcss, sigmar_1d, p_dst(:,i), regrid_mval, &
                           i, j_rs)
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
               sdpsum = sdpsum + sigma_1d(k)*(p_src(k+1,i) - p_src(k,i))
            enddo
            smean = sdpsum/(p_src(kk+1,i) - p_src(1,i))
            ks = 2
            do while (ks <= kk)
               if (smean < sigmar_1d(ks)) exit
               ks = ks + 1
            enddo
            do k = ks, kk
               p_dst(k,i) = p_src(kk+1,i)
            enddo
            ke = ks - 1
         endif

         ! Modify interface pressures so that layer thicknesses are
         ! above a specified threshold.
         dpmin_max = (p_src(kk+1,i) - p_src(1,i))/kk
         dpmin_max = dpmin_surface
         dpmin_int = min(dpmin_max, dpmin_surface, dpmin_interior)
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
                                         p_dst(k-1,i) + dpmin_int)
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
         dpmin_sfc = min(dpmin_max, dpmin_surface)
         pmin = p_src(1,i) + dpmin_sfc
         dpt = dpmin_sfc
         do k = 2, ke
            dpmin_sfc = dpmin_sfc*dpmin_inflation_factor
            dpt = max(p_dst(k+1,i) - p_dst(k,i), dpt, dpmin_sfc)
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
            p_dst(k,i) = min(p_dst(ke+1,i), max(p_dst(k-1,i) + dpmin_int, pt))
            pmin = pmin + dpmin_sfc
         enddo

      enddo
      enddo

   end subroutine cntiso_regrid_direct_jslice

   subroutine cntiso_regrid_nudge_jslice(p_src, p_dst, i_lb, i_ub, j, j_rs, nn)

      real(r8), dimension(:,1-nbdy:), intent(in) :: p_src
      real(r8), dimension(:,1-nbdy:), intent(out) :: p_dst
      integer, intent(in) :: i_lb, i_ub, j, j_rs, nn

      integer, parameter :: &
         p_ord = 4, &
         it = 1, &
         is = 2

      real(r8), dimension(p_ord+1,kdm,2,1-nbdy:idm+nbdy) :: tpc_src
      real(r8), dimension(2,kdm,2,1-nbdy:idm+nbdy) :: t_srcdi
      real(r8), dimension(2,kdm) :: sig_srcdi
      integer, dimension(1-nbdy:idm+nbdy) :: ksmx, kdmx

      real(r8), dimension(kdm+1) :: sigmar_1d, pmin, sig_pmin
      real(r8) :: sig_max, dpmin_sfc, dsig, dsigdx, q
      integer :: l, i, nt, k, kr, kl, klastok, kt, errstat
      logical :: ok

      do l = 1, isp(j)
      do i = max(i_lb, ifp(j, l)), min(i_ub, ilp(j, l))

         ! Extract polynomial coefficients of the reconstructions.
         do nt = 1, 2
            errstat = extract_polycoeff(trc_rcss(nt), &
                                        tpc_src(:,:,nt,i), i, j_rs)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(ndiff_prep_jslice)')
                      stop '(ndiff_prep_jslice)'
            endif
         enddo

         ! Find index of deepest source layer with non-zero thickness.
         ksmx(i) = kk
         do k = kk, 1, -1
            if (p_src(k,i) == p_src(kk+1,i)) ksmx(i) = k - 1
         enddo

         ! Store variables in dual interface arrays with with values
         ! corresponding to upper and lower interface of each layer. Also find
         ! the maximum lower interface potential density of the reconstructed
         ! column.
         sig_max = 0._r8
         do k = 1, ksmx(i)
            do nt = 1, 2
               t_srcdi(1,k,nt,i) = peval0(tpc_src(:,k,nt,i))
               t_srcdi(2,k,nt,i) = peval1(tpc_src(:,k,nt,i))
            enddo
            sig_srcdi(1,k) = sig(t_srcdi(1,k,it,i), t_srcdi(1,k,is,i))
            sig_srcdi(2,k) = sig(t_srcdi(2,k,it,i), t_srcdi(2,k,is,i))
            sig_max = max(sig_max, sig_srcdi(2,k))
         enddo

         ! Copy variables into 1D arrays.
         do k = 1, kk
            sigmar_1d(k) = sigmar(i,j,k)
         enddo
         sigmar_1d(kk+1) = sigmar_1d(kk)

         ! Find the index of the first layer which lower interface reference
         ! potential density is denser than the maximum lower interface
         ! potential density of the reconstructed column.
         do k = kk, 1, -1
            if (sigmar_1d(k) < sig_max) exit
         enddo
         kdmx(i) = max(1, k)

         do k = kdmx(i)+1, kk+1
            p_dst(k,i) = p_src(kk+1,i)
         enddo

         dpmin_sfc = dpmin_surface
         pmin(1) = p_src(1,i)
         do k = 1, kk
            pmin(k+1) = min(pmin(k) + dpmin_sfc, p_src(kk+1,i))
            dpmin_sfc = dpmin_sfc*dpmin_inflation_factor
         enddo
         p_dst(1,i) = pmin(1)

         sig_pmin(1) = sig_srcdi(1,1)
         kr = 2
         kl = 1
         do while (kr <= kdmx(i))
            do while (p_src(kl+1,i) < pmin(kr))
               kl = kl + 1
            enddo
            sig_pmin(kr) = ( (p_src(kl+1,i) - pmin(kr))*sig_srcdi(1,kl) &
                           + (pmin(kr) - p_src(kl,i))*sig_srcdi(2,kl)) &
                           /(p_src(kl+1,i) - p_src(kl,i))
            if (sigmar_1d(kr) > sig_pmin(kr)) exit
            p_dst(kr,i) = pmin(kr)
            kr = kr + 1
         enddo

         klastok = kr - 1
         do k = kr, min(ksmx(i), kdmx(i))
            ok = .true.
            if      (sigmar_1d(k) < sig_srcdi(2,k-1) .and. &
                     sigmar_1d(k) < sig_srcdi(1,k  )) then
               dsig = (sigmar_1d(k) - sig_srcdi(2,k-1))*regrid_nudge_factor
               dsigdx = dsigdt(t_srcdi(2,k-1,it,i), t_srcdi(2,k-1,is,i)) &
                        *dpeval1(tpc_src(:,k-1,it,i)) &
                      + dsigds(t_srcdi(2,k-1,it,i), t_srcdi(2,k-1,is,i)) &
                        *dpeval1(tpc_src(:,k-1,is,i))
               if (- dsig > .5_r8*dsigdx) then
                  dsigdx = sig_srcdi(2,k-1) - sig_srcdi(1,k-1)
                  if (- dsig > .5_r8*dsigdx) ok = .false.
               endif
               if (ok) p_dst(k,i) = p_src(k,i) &
                                  + dsig*(p_src(k,i) - p_src(k-1,i))/dsigdx
            elseif  (sigmar_1d(k) > sig_srcdi(2,k-1) .and. &
                     sigmar_1d(k) > sig_srcdi(1,k  )) then
               dsig = (sigmar_1d(k) - sig_srcdi(1,k))*regrid_nudge_factor
               dsigdx = dsigdt(t_srcdi(1,k,it,i), t_srcdi(1,k,is,i)) &
                        *dpeval0(tpc_src(:,k,it,i)) &
                      + dsigds(t_srcdi(1,k,it,i), t_srcdi(1,k,is,i)) &
                        *dpeval0(tpc_src(:,k,is,i))
               if (dsig > .5_r8*dsigdx) then
                  dsigdx = sig_srcdi(2,k) - sig_srcdi(1,k)
                  if (dsig > .5_r8*dsigdx) ok = .false.
               endif
               if (ok) p_dst(k,i) = p_src(k,i) &
                                  + dsig*(p_src(k+1,i) - p_src(k,i))/dsigdx
            else
               p_dst(k,i) = p_src(k,i)
            endif
            if (ok) then
               p_dst(k,i) = &
                  min(max(p_dst(k,i), pmin(k), &
                          p_dst(klastok,i) + (k - klastok)*dpmin_interior), &
                      p_src(kk+1,i))
               if (k - klastok > 1) then
                  q = (p_dst(k,i) - p_dst(klastok,i))/(k - klastok)
                  do kt = klastok+1, k-1
                     p_dst(kt,i) = min(max(p_dst(kt-1,i) + q, pmin(kt)), &
                                       p_src(kk+1,i))
                  enddo
               endif
               klastok = k
            endif
         enddo

         do k = max(kr, min(ksmx(i), kdmx(i))) + 1, kdmx(i)
            ok = .true.
            if (sigmar_1d(k) < sig_srcdi(2,ksmx(i))) then
               dsig = (sigmar_1d(k) - sig_srcdi(2,ksmx(i)))*regrid_nudge_factor
               dsigdx = dsigdt(t_srcdi(2,ksmx(i),it,i), &
                               t_srcdi(2,ksmx(i),is,i)) &
                        *dpeval1(tpc_src(:,ksmx(i),it,i)) &
                      + dsigds(t_srcdi(2,ksmx(i),it,i), &
                               t_srcdi(2,ksmx(i),is,i)) &
                        *dpeval1(tpc_src(:,ksmx(i),is,i))
               if (- dsig > .5_r8*dsigdx) then
                  dsigdx = sig_srcdi(2,ksmx(i)) - sig_srcdi(1,ksmx(i))
                  if (- dsig > .5_r8*dsigdx) ok = .false.
               endif
               if (ok) p_dst(k,i) = p_src(kk+1,i) &
                                  + dsig*(p_src(kk+1,i) - p_src(ksmx(i),i)) &
                                    /dsigdx
            else
               p_dst(k,i) = p_src(kk+1,i)
            endif
            if (ok) then
               p_dst(k,i) = &
                  min(max(p_dst(k,i), pmin(k), &
                          p_dst(klastok,i) + (k - klastok)*dpmin_interior), &
                      p_src(kk+1,i))
               if (k - klastok > 1) then
                  q = (p_dst(k,i) - p_dst(klastok,i))/(k - klastok)
                  do kt = klastok+1, k-1
                     p_dst(kt,i) = min(max(p_dst(kt-1,i) + q, pmin(kt)), &
                                       p_src(kk+1,i))
                  enddo
               endif
               klastok = k
            endif
         enddo

         if (kdmx(i) - klastok > 0) then
            q = (p_dst(kdmx(i)+1,i) - p_dst(klastok,i))/(kdmx(i) + 1 - klastok)
            do kt = klastok+1, kdmx(i)
               p_dst(kt,i) = min(max(p_dst(kt-1,i) + q, pmin(kt)), &
                                 p_src(kk+1,i))
            enddo
         endif

      enddo
      enddo

   end subroutine cntiso_regrid_nudge_jslice

   subroutine copy_jslice_to_3d(p_dst, trc_rm, i_lb, i_ub, j, nn)

      real(r8), dimension(:,1-nbdy:), intent(in) :: p_dst
      real(r8), dimension(:,:,1-nbdy:), intent(in) :: trc_rm

      integer, intent(in) :: i_lb, i_ub, j, nn

      integer :: l, i, k, kn, nt

      do l = 1, isp(j)
      do i = max(i_lb, ifp(j,l)), min(i_ub, ilp(j,l))

         do k = 1, kk
            kn = k + nn
            temp(i,j,kn) = trc_rm(k,1,i)
            saln(i,j,kn) = trc_rm(k,2,i)
            dp(i,j,kn) = p_dst(k+1,i) - p_dst(k,i)
            sigma(i,j,kn) = sig(trc_rm(k,1,i), trc_rm(k,2,i))
#ifdef TRC
            do nt = 1, ntr
               trc(i,j,kn,nt) = trc_rm(k,nt+2,i)
            enddo
#endif
         enddo

      enddo
      enddo

   end subroutine copy_jslice_to_3d

   subroutine readnml_vcoord
   ! ---------------------------------------------------------------------------
   ! Read variables in the namelist group 'vcoord' and resolve options.
   ! ---------------------------------------------------------------------------

      character(len = 80) :: nml_fname
      integer :: ios
      logical :: fexist

      namelist /vcoord/ &
         vcoord_type, reconstruction_method, &
         density_limiting, tracer_limiting, velocity_limiting, &
         density_pc_upper_bndr, density_pc_lower_bndr, &
         tracer_pc_upper_bndr, tracer_pc_lower_bndr, &
         velocity_pc_upper_bndr, velocity_pc_lower_bndr, &
         dpmin_surface, dpmin_inflation_factor, dpmin_interior, &
         regrid_nudge_factor

      ! Read variables in the namelist group 'vcoord'.
      if (mnproc == 1) then
         nml_fname = 'ocn_in'//trim(inst_suffix)
         inquire(file = nml_fname, exist = fexist)
         if (fexist) then
            open (unit = nfu, file = nml_fname, status = 'old', action = 'read')
         else
            nml_fname = 'limits'//trim(inst_suffix)
            inquire(file = nml_fname, exist = fexist)
            if (fexist) then
               open (unit = nfu, file = nml_fname, status = 'old', &
                     action = 'read')
            else
               write (lp,*) 'readnml_vcoord: could not find namelist file!'
               call xchalt('(readnml_vcoord)')
                      stop '(readnml_vcoord)'
            endif
         endif
         read (unit = nfu, nml = vcoord, iostat = ios)
         close (unit = nfu)
      endif
      call xcbcst(ios)
      if (ios /= 0) then
         if (mnproc == 1) &
            write (lp,*) 'readnml_vcoord: No vertical coordinate variable '//  &
                         'group found in namelist. Using defaults.'
      else
        call xcbcst(vcoord_type)
        call xcbcst(reconstruction_method)
        call xcbcst(density_limiting)
        call xcbcst(tracer_limiting)
        call xcbcst(velocity_limiting)
        call xcbcst(density_pc_upper_bndr)
        call xcbcst(density_pc_lower_bndr)
        call xcbcst(tracer_pc_upper_bndr)
        call xcbcst(tracer_pc_lower_bndr)
        call xcbcst(velocity_pc_upper_bndr)
        call xcbcst(velocity_pc_lower_bndr)
        call xcbcst(dpmin_surface)
        call xcbcst(dpmin_inflation_factor)
        call xcbcst(dpmin_interior)
        call xcbcst(regrid_nudge_factor)
      endif
      if (mnproc == 1) then
         write (lp,*) 'readnml_vcoord: vertical coordinate variables:'
         write (lp,*) '  vcoord_type =            ', &
                      trim(vcoord_type)
         write (lp,*) '  reconstruction_method =  ', &
                      trim(reconstruction_method)
         write (lp,*) '  density_limiting =       ', &
                      trim(density_limiting)
         write (lp,*) '  tracer_limiting =        ', &
                      trim(tracer_limiting)
         write (lp,*) '  velocity_limiting =      ', &
                      trim(velocity_limiting)
         write (lp,*) '  density_pc_upper_bndr =  ', density_pc_upper_bndr
         write (lp,*) '  density_pc_lower_bndr =  ', density_pc_lower_bndr
         write (lp,*) '  tracer_pc_upper_bndr =   ', tracer_pc_upper_bndr
         write (lp,*) '  tracer_pc_lower_bndr =   ', tracer_pc_lower_bndr
         write (lp,*) '  velocity_pc_upper_bndr = ', velocity_pc_upper_bndr
         write (lp,*) '  velocity_pc_lower_bndr = ', velocity_pc_lower_bndr
         write (lp,*) '  dpmin_surface =          ', dpmin_surface
         write (lp,*) '  dpmin_inflation_factor = ', dpmin_inflation_factor
         write (lp,*) '  dpmin_interior =         ', dpmin_interior
         write (lp,*) '  regrid_nudge_factor =    ', regrid_nudge_factor
      endif

      ! Resolve options.
      select case (trim(vcoord_type))
         case ('isopyc_bulkml')
            vcoord_type_tag = isopyc_bulkml
         case ('cntiso_hybrid')
            vcoord_type_tag = cntiso_hybrid
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') &
                  ' readnml_vcoord: vcoord_type = ', &
                  trim(vcoord_type), ' is unsupported!'
            call xcstop('(readnml_vcoord)')
                   stop '(readnml_vcoord)'
      end select
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
                  ' readnml_vcoord: reconstruction_method = ', &
                  trim(reconstruction_method), ' is unsupported!'
            call xcstop('(readnml_vcoord)')
                   stop '(readnml_vcoord)'
      end select
      select case (trim(density_limiting))
         case ('monotonic')
            density_limiting_tag = hor3map_monotonic
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') &
                  ' readnml_vcoord: density_limiting = ', &
                  trim(density_limiting), ' is unsupported!'
            call xcstop('(readnml_vcoord)')
                   stop '(readnml_vcoord)'
      end select
      select case (trim(tracer_limiting))
         case ('monotonic')
            tracer_limiting_tag = hor3map_monotonic
         case ('non_oscillatory')
            tracer_limiting_tag = hor3map_non_oscillatory
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') &
                  ' readnml_vcoord: tracer_limiting = ', &
                  trim(tracer_limiting), ' is unsupported!'
            call xcstop('(readnml_vcoord)')
                   stop '(readnml_vcoord)'
      end select
      select case (trim(velocity_limiting))
         case ('monotonic')
            velocity_limiting_tag = hor3map_monotonic
         case ('non_oscillatory')
            velocity_limiting_tag = hor3map_non_oscillatory
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') &
                  ' readnml_vcoord: velocity_limiting = ', &
                  trim(velocity_limiting), ' is unsupported!'
            call xcstop('(readnml_vcoord)')
                   stop '(readnml_vcoord)'
      end select

      ! Change units from [m] to [g cm-1 s-2] of depth interval variables.
      dpmin_surface = dpmin_surface*onem
      dpmin_interior = dpmin_interior*onem

   end subroutine readnml_vcoord

   subroutine inivar_vcoord
   ! ---------------------------------------------------------------------------
   ! Initialize arrays and data structures.
   ! ---------------------------------------------------------------------------

      integer :: i, j, k, nt, errstat

   !$omp parallel do private(i, k)
      do j = 1-nbdy, jj+nbdy
         do k = 1, kk
            do i = 1-nbdy, ii+nbdy
               sigmar(i,j,k) = spval
            enddo
         enddo
      enddo
   !$omp end parallel do

      ! Configuration of the reconstruction data structure that only depends on
      ! the source grid.
      rcgs%n_src = kk
      if (ltedtp_opt == ltedtp_neutral) then
         rcgs%i_lbound = 0
         rcgs%i_ubound = ii + 1
      else
         rcgs%i_ubound = ii
      endif
      rcgs%j_ubound = 2
      rcgs%method = reconstruction_method_tag

      ! Configuration of reconstruction data structures that is specific to
      ! various source data.

      d_rcss%limiting = density_limiting_tag
      d_rcss%pc_left_bndr = density_pc_upper_bndr
      d_rcss%pc_right_bndr = density_pc_lower_bndr

      trc_rcss(1)%limiting = tracer_limiting_tag
      trc_rcss(1)%pc_left_bndr = tracer_pc_upper_bndr
      trc_rcss(1)%pc_right_bndr = tracer_pc_lower_bndr
      if (tracer_limiting_tag == hor3map_non_oscillatory) then
#ifdef TRC
         do nt = 2, ntr_loc
            trc_rcss(nt)%limiting = hor3map_non_oscillatory_posdef
            trc_rcss(nt)%pc_left_bndr = tracer_pc_upper_bndr
            trc_rcss(nt)%pc_right_bndr = tracer_pc_lower_bndr
         enddo
#endif
      else
#ifdef TRC
         do nt = 2, ntr_loc
            trc_rcss(nt)%limiting = tracer_limiting_tag
            trc_rcss(nt)%pc_left_bndr = tracer_pc_upper_bndr
            trc_rcss(nt)%pc_right_bndr = tracer_pc_lower_bndr
         enddo
#endif
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
         call xchalt('(inivar_vcoord)')
                stop '(inivar_vcoord)'
      endif

      errstat = initialize_rcss(rcgs, d_rcss)
      if (errstat /= hor3map_noerr) then
         write(lp,*) trim(hor3map_errstr(errstat))
         call xchalt('(inivar_vcoord)')
                stop '(inivar_vcoord)'
      endif

      do nt = 1, ntr_loc
         errstat = initialize_rcss(rcgs, trc_rcss(nt))
         if (errstat /= hor3map_noerr) then
            write(lp,*) trim(hor3map_errstr(errstat))
            call xchalt('(inivar_vcoord)')
                   stop '(inivar_vcoord)'
         endif
      enddo

      errstat = initialize_rcss(rcgs, v_rcss)
      if (errstat /= hor3map_noerr) then
         write(lp,*) trim(hor3map_errstr(errstat))
         call xchalt('(inivar_vcoord)')
                stop '(inivar_vcoord)'
      endif

      errstat = initialize_rms(rcgs, rms)
      if (errstat /= hor3map_noerr) then
         write(lp,*) trim(hor3map_errstr(errstat))
         call xchalt('(inivar_vcoord)')
                stop '(inivar_vcoord)'
      endif

   end subroutine inivar_vcoord

   subroutine cntiso_hybrid_regrid_direct_remap(m, n, mm, nn, k1m, k1n)

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      real(r8), dimension(kdm+1,1-nbdy:idm+nbdy,2) :: p_src_rs, p_dst_rs
      real(r8), dimension(kdm,ntr_loc,1-nbdy:idm+nbdy) :: trc_rm
      integer :: j_rs, jm_rs, jp_rs, j, nt

      if (ltedtp_opt /= ltedtp_neutral) then

         j_rs = 1

         do j = 1, jj
            call prep_recon_jslice(p_src_rs(:,:,j_rs), 1, ii, j, j_rs, nn)
            call recon_trc_jslice(1, ii, j, j_rs, nn)
            call cntiso_regrid_direct_jslice(p_src_rs(:,:,j_rs), p_dst_rs(:,:,j_rs), &
                                      1, ii, j, j_rs, nn)
            call remap_trc_jslice(p_dst_rs(:,:,j_rs), trc_rm, &
                                  1, ii, j, j_rs)
            call copy_jslice_to_3d(p_dst_rs(:,:,j_rs), trc_rm, 1, ii, j, nn)
         enddo

      else

         call xctilr(dp   (1-nbdy,1-nbdy,k1n), 1, kk, 1, 1, halo_ps)
         call xctilr(temp (1-nbdy,1-nbdy,k1n), 1, kk, 1, 1, halo_ps)
         call xctilr(saln (1-nbdy,1-nbdy,k1n), 1, kk, 1, 1, halo_ps)
         call xctilr(sigma(1-nbdy,1-nbdy,k1n), 1, kk, 1, 1, halo_ps)
#ifdef TRC
         do nt = 1, ntr
!#  if defined(TKE) && !defined(TKEIDF)
!            if (nt == itrtke .or. nt == itrgls) cycle
!#  endif
            call xctilr(trc(1-nbdy,1-nbdy,k1n,nt), 1, kk, 1, 1, halo_ps)
         enddo
#endif
         call xctilr(difiso, 1,kk, 1,1, halo_ps)

         jm_rs = 1
         jp_rs = 2

         do j = -1, 0
            jm_rs = 3 - jm_rs
            jp_rs = 3 - jp_rs
            call prep_recon_jslice(p_src_rs(:,:,jp_rs), &
                                   0, ii+1, j+1, jp_rs, nn)
            call recon_trc_jslice(0, ii+1, j+1, jp_rs, nn)
            call cntiso_regrid_direct_jslice(p_src_rs(:,:,jp_rs), &
                                      p_dst_rs(:,:,jp_rs), &
                                      0, ii+1, j+1, jp_rs, nn)
            call ndiff_prep_jslice(p_src_rs, p_dst_rs, trc_rcss, &
                                   0, ii+1, j+1, jp_rs, mm)
         enddo

         j = 0
         call ndiff_vflx_jslice(p_dst_rs, 1, ii, j+1, jp_rs, mm, nn)

         do j = 1, jj
            jm_rs = 3 - jm_rs
            jp_rs = 3 - jp_rs
            call prep_recon_jslice(p_src_rs(:,:,jp_rs), &
                                   0, ii+1, j+1, jp_rs, nn)
            call recon_trc_jslice(0, ii+1, j+1, jp_rs, nn)
            call cntiso_regrid_direct_jslice(p_src_rs(:,:,jp_rs), &
                                      p_dst_rs(:,:,jp_rs), &
                                      0, ii+1, j+1, jp_rs, nn)
            call ndiff_prep_jslice(p_src_rs, p_dst_rs, trc_rcss, &
                                   0, ii+1, j+1, jp_rs, mm)
            call ndiff_uflx_jslice(p_dst_rs, 1, ii+1, j, jm_rs, mm, nn)
            call ndiff_vflx_jslice(p_dst_rs, 1, ii, j+1, jp_rs, mm, nn)
            call remap_trc_jslice(p_dst_rs(:,:,jm_rs), trc_rm, &
                                  1, ii, j, jm_rs)
            call ndiff_update_trc_jslice(p_dst_rs, trc_rm, 1, ii, j, jm_rs)
            call copy_jslice_to_3d(p_dst_rs(:,:,jm_rs), trc_rm, 1, ii, j, nn)
         enddo

      endif

      if (csdiag) then
         if (mnproc == 1) then
            write (lp,*) 'cntiso_hybrid_regrid_direct_remap:'
         endif
         call chksummsk(dp   (1-nbdy,1-nbdy,k1n), ip, kk, 'dp')
         call chksummsk(temp (1-nbdy,1-nbdy,k1n), ip, kk, 'temp')
         call chksummsk(saln (1-nbdy,1-nbdy,k1n), ip, kk, 'saln')
         call chksummsk(sigma(1-nbdy,1-nbdy,k1n), ip, kk, 'sigma')
         call chksummsk(sigmar, ip, kk, 'sigmar')
#ifdef TRC
         do nt = 1, ntr
            call chksummsk(trc(1-nbdy,1-nbdy,k1n,nt), ip, kk, 'trc')
         enddo
#endif
      endif

   end subroutine cntiso_hybrid_regrid_direct_remap

   subroutine cntiso_hybrid_regrid_remap(m, n, mm, nn, k1m, k1n)

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      real(r8), dimension(kdm+1,1-nbdy:idm+nbdy,2) :: p_src_rs, p_dst_rs
      real(r8), dimension(kdm,ntr_loc,1-nbdy:idm+nbdy) :: trc_rm
      integer :: j_rs, jm_rs, jp_rs, j, nt

      if (ltedtp_opt /= ltedtp_neutral) then

         j_rs = 1

         do j = 1, jj
            call prep_recon_jslice(p_src_rs(:,:,j_rs), 1, ii, j, j_rs, nn)
            call recon_trc_jslice(1, ii, j, j_rs, nn)
            call cntiso_regrid_nudge_jslice(p_src_rs(:,:,j_rs), p_dst_rs(:,:,j_rs), &
                                      1, ii, j, j_rs, nn)
            call remap_trc_jslice(p_dst_rs(:,:,j_rs), trc_rm, &
                                  1, ii, j, j_rs)
            call copy_jslice_to_3d(p_dst_rs(:,:,j_rs), trc_rm, 1, ii, j, nn)
         enddo

      else

         call xctilr(dp   (1-nbdy,1-nbdy,k1n), 1, kk, 1, 1, halo_ps)
         call xctilr(temp (1-nbdy,1-nbdy,k1n), 1, kk, 1, 1, halo_ps)
         call xctilr(saln (1-nbdy,1-nbdy,k1n), 1, kk, 1, 1, halo_ps)
         call xctilr(sigma(1-nbdy,1-nbdy,k1n), 1, kk, 1, 1, halo_ps)
#ifdef TRC
         do nt = 1, ntr
!#  if defined(TKE) && !defined(TKEIDF)
!            if (nt == itrtke .or. nt == itrgls) cycle
!#  endif
            call xctilr(trc(1-nbdy,1-nbdy,k1n,nt), 1, kk, 1, 1, halo_ps)
         enddo
#endif
         call xctilr(difiso, 1,kk, 1,1, halo_ps)

         jm_rs = 1
         jp_rs = 2

         do j = -1, 0
            jm_rs = 3 - jm_rs
            jp_rs = 3 - jp_rs
            call prep_recon_jslice(p_src_rs(:,:,jp_rs), &
                                   0, ii+1, j+1, jp_rs, nn)
            call recon_trc_jslice(0, ii+1, j+1, jp_rs, nn)
            call cntiso_regrid_nudge_jslice(p_src_rs(:,:,jp_rs), &
                                      p_dst_rs(:,:,jp_rs), &
                                      0, ii+1, j+1, jp_rs, nn)
            call ndiff_prep_jslice(p_src_rs, p_dst_rs, trc_rcss, &
                                   0, ii+1, j+1, jp_rs, mm)
         enddo

         j = 0
         call ndiff_vflx_jslice(p_dst_rs, 1, ii, j+1, jp_rs, mm, nn)

         do j = 1, jj
            jm_rs = 3 - jm_rs
            jp_rs = 3 - jp_rs
            call prep_recon_jslice(p_src_rs(:,:,jp_rs), &
                                   0, ii+1, j+1, jp_rs, nn)
            call recon_trc_jslice(0, ii+1, j+1, jp_rs, nn)
            call cntiso_regrid_nudge_jslice(p_src_rs(:,:,jp_rs), &
                                      p_dst_rs(:,:,jp_rs), &
                                      0, ii+1, j+1, jp_rs, nn)
            call ndiff_prep_jslice(p_src_rs, p_dst_rs, trc_rcss, &
                                   0, ii+1, j+1, jp_rs, mm)
            call ndiff_uflx_jslice(p_dst_rs, 1, ii+1, j, jm_rs, mm, nn)
            call ndiff_vflx_jslice(p_dst_rs, 1, ii, j+1, jp_rs, mm, nn)
            call remap_trc_jslice(p_dst_rs(:,:,jm_rs), trc_rm, &
                                  1, ii, j, jm_rs)
            call ndiff_update_trc_jslice(p_dst_rs, trc_rm, 1, ii, j, jm_rs)
            call copy_jslice_to_3d(p_dst_rs(:,:,jm_rs), trc_rm, 1, ii, j, nn)
         enddo

      endif

      if (csdiag) then
         if (mnproc == 1) then
            write (lp,*) 'cntiso_hybrid_regrid_remap:'
         endif
         call chksummsk(dp   (1-nbdy,1-nbdy,k1n), ip, kk, 'dp')
         call chksummsk(temp (1-nbdy,1-nbdy,k1n), ip, kk, 'temp')
         call chksummsk(saln (1-nbdy,1-nbdy,k1n), ip, kk, 'saln')
         call chksummsk(sigma(1-nbdy,1-nbdy,k1n), ip, kk, 'sigma')
         call chksummsk(sigmar, ip, kk, 'sigmar')
#ifdef TRC
         do nt = 1, ntr
            call chksummsk(trc(1-nbdy,1-nbdy,k1n,nt), ip, kk, 'trc')
         enddo
#endif
      endif

   end subroutine cntiso_hybrid_regrid_remap

   subroutine remap_velocity(m, n, mm, nn, k1m, k1n)

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      real(r8), dimension(kdm+1) :: p_1d, p_dst_1d
      real(r8), dimension(kdm) :: u_1d, v_1d
      real(r8) :: q
      integer :: i, j, k, l, kn, errstat

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
      do j = -2, jj+2
         do k = 1, kk
            kn = k + nn
            do l = 1, isp(j)
            do i = max(-2, ifp(j,l)), min(ii+2, ilp(j,l))
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
               call xchalt('(remap_velocity)')
                      stop '(remap_velocity)'
            endif

            ! Prepare remapping to layer structure with regridded interface
            ! pressures.
            errstat = prepare_remapping(rcgs, rms, p_dst_1d, i, 1)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(remap_velocity)')
                      stop '(remap_velocity)'
            endif

            ! Reconstruct and remap u-component of velocity.
            errstat = reconstruct(rcgs, v_rcss, u_1d, i, 1)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(remap_velocity)')
                      stop '(remap_velocity)'
            endif
            errstat = remap(v_rcss, rms, u_1d, i, 1)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(remap_velocity)')
                      stop '(remap_velocity)'
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
               call xchalt('(remap_velocity)')
                      stop '(remap_velocity)'
            endif

            ! Prepare remapping to layer structure with regridded interface
            ! pressures.
            errstat = prepare_remapping(rcgs, rms, p_dst_1d, i, 1)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(remap_velocity)')
                      stop '(remap_velocity)'
            endif

            ! Reconstruct and remap v-component of velocity.
            errstat = reconstruct(rcgs, v_rcss, v_1d, i, 1)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(remap_velocity)')
                      stop '(remap_velocity)'
            endif
            errstat = remap(v_rcss, rms, v_1d, i, 1)
            if (errstat /= hor3map_noerr) then
               write(lp,*) trim(hor3map_errstr(errstat))
               call xchalt('(remap_velocity)')
                      stop '(remap_velocity)'
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
            write (lp,*) 'remap_velocity:'
         endif
         call chksummsk(dpu(1-nbdy,1-nbdy,k1n), iu, kk, 'dpu')
         call chksummsk(dpv(1-nbdy,1-nbdy,k1n), iv, kk, 'dpv')
         call chksummsk(u  (1-nbdy,1-nbdy,k1n), iu, kk, 'u')
         call chksummsk(v  (1-nbdy,1-nbdy,k1n), iv, kk, 'v')
      endif

   end subroutine remap_velocity

end module mod_vcoord
