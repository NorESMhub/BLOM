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

module mod_vcoord
! ------------------------------------------------------------------------------
! This module contains parameter, variables and procedures related to the
! vertical coordinate.
! ------------------------------------------------------------------------------

   use mod_types,     only: r8
   use mod_config,    only: inst_suffix
   use mod_constants, only: grav, spval, onem, onecm
   use mod_time,      only: nstep, nstep_in_day, nday_of_year, nday_in_year, &
                            baclin
   use mod_xc
   use mod_grid,      only: scp2, area
   use mod_state,     only: dp, sigma, p
   use mod_cmnfld,    only: dpml
   use mod_utility,   only: util1
   use mod_checksum,  only: csdiag, chksummsk

   implicit none
   private

   ! Derived data types.

   type :: sigref_fun_spec_type

     real(r8) :: &
        dsdz_bot, & ! Derivative of sigma with respect to z at Bezier point 4
                    ! (z = 1) [kg m-3]. This defines the tangent line of the
                    ! Bezier curve at Bezier point 4 and Bezier points 2 and 3
                    ! lie on this line.
        sp1, &      ! Sigma value at Bezier point 1 (z = 0) [kg m-3].
        zp2, &      ! z value at Bezier point 2 [].
        zp3, &      ! z value at Bezier point 3 [].
        sp4, &      ! Sigma value at Bezier point 4 (z = 1) [kg m-3].
        z_top, &    ! z value at transition from Bezier curve to parabola
                    ! covering the z-range [0, z_top] [].
        s_top, &    ! Sigma value of parabola at z = 0 [kg m-3].
        z_bot, &    ! z value at transition from Bezier curve to parabola
                    ! covering the z-range [z_bot, 1] [].
        s_bot       ! Sigma value of parabola at z = 1 [kg m-3].
      
   end type sigref_fun_spec_type

   ! Parameters:
   real(r8), parameter :: &
      z_eps = 1.e-12_r8, &
      t_tol = 1.e-12_r8
   integer, parameter :: &
      vcoord_isopyc_bulkml = 1, &  ! Vertical coordinate type: bulk surface
                                   ! mixed layer with isopycnic layers below.
      vcoord_cntiso_hybrid = 2, &  ! Vertical coordinate type: Hybrid
                                   ! coordinate with pressure coordinates
                                   ! towards the surface and continuous
                                   ! isopycnal below.
      vcoord_plevel        = 3, &  ! Vertical coordinate type: pressure
                                   ! coordinate.
      sra_tlev_num         = 12, & ! Number of time levels of yearly statistics
                                   ! for sigref adaptation.
      kdm_max              = 1000  ! Maximum anticipated vertical dimension.

   ! Options with default values, modifiable by namelist.
   character(len = 80) :: &
      vcoord_type            = 'isopyc_bulkml', &
      sigref_spec            = 'inicon', &
      plevel_spec            = 'inflation'
   real(r8) :: &
      dpmin_surface          = 1.5_r8, &
      dpmin_inflation_factor = 1._r8, &
      sra_clim_ts            = 5._r8, &
      sra_param_ts           = 5._r8, &
      sra_massfrac_bot       = .01, &
      sra_massfrac_eps       = .0001
   type(sigref_fun_spec_type) :: &
      sigref_fun_spec
   real(r8), dimension(kdm_max) :: &
      sigref                 = spval, &
      plevel                 = spval
   logical :: &
      sigref_adaption        = .false.

   ! Options derived from string options.
   integer :: &
      vcoord_tag

   real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: &
      sigmar, &                   ! Reference potential density [kg m-3].
      sra_massdc_colsum, &        ! Column sum of mass within density classes
      sra_sigmassdc_colsum        ! Column sum of potential density times mass
                                  ! within density classes.
   real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,sra_tlev_num) :: &
      sra_dpml_sum, &             ! Sum of pressure at mixed layer base.
      sra_sigmlb_sum, &           ! Sum of potential density at mixed layer
                                  ! base.
      sra_dpml_clim, &            ! Climatology of pressure at mixed layer base.
      sra_sigmlb_clim             ! Climatology of potential density at mixed
                                  ! layer base.
   real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
      sra_dpml_dmax, &            ! Daily maximum mixed layer pressure
                                  ! thickness.
      sra_sigmlb_dmax, &          ! Potential density at daily maximum mixed
                                  ! layer depth.
      sra_massgs_colsum, &        ! Column sum of total mass.
      sra_cost_wgt                ! Grid-cell weight in cost estimation.
   real(r8), dimension(kdm) :: &
      sra_sigref_sum              ! Sum of reference potential densities.
   real(r8) :: &
      sra_s_bot_sum               ! Sum of sigma value of parabola at z = 1.
   integer, dimension(sra_tlev_num) :: &
      sra_tlev_accnum             ! Number of accumulated fields in time levels
                                  ! of yearly statistics for sigref adaptation.
   integer :: &
      sra_accnum                  ! Number of accumulated reference potential
                                  ! density parameters.
   type(sigref_fun_spec_type) :: &
      sigref_fun_spec_old, &      ! 
      sigref_fun_spec_new

   public :: vcoord_type, vcoord_tag, vcoord_isopyc_bulkml, &
             vcoord_cntiso_hybrid, vcoord_plevel, sra_tlev_num, sigref_spec, &
             sigmar, sigref_fun_spec, sigref, plevel, sigref_adaption, &
             sra_massdc_colsum, sra_sigmassdc_colsum, sra_massgs_colsum, &
             sra_dpml_sum, sra_sigmlb_sum, sra_dpml_clim, sra_sigmlb_clim, &
             sra_sigref_sum, sra_s_bot_sum, sra_tlev_accnum, sra_accnum, &
             sigref_fun_spec_old, sigref_fun_spec_new, &
             readnml_vcoord, inivar_vcoord, extract_sigref, sigref_adapt

contains

   pure function cubic_root(a, b, c, d, x_ini, x_tol) result(x_new)
   ! ---------------------------------------------------------------------------
   ! Find cubic root by Newton-Raphson iterative method.
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: a, b, c, d, x_ini, x_tol
      real(r8) :: x_new

      real(r8) :: x

      x = x_ini
      do
         x_new = x - (((a*x + b)*x + c)*x + d)/((3._r8*a*x + 2._r8*b)*x + c)
         if (abs(x_new - x) < x_tol) return
         x = x_new
      enddo

   end function cubic_root

   pure function sigref_fun(fun_spec, kmax) result(sigref)
   ! ---------------------------------------------------------------------------
   ! Return reference potential densities based on a functional expression
   ! specified by the parameters in fun_spec of type sigref_fun_spec_type. The
   ! functional expression consist of a cubic Bezier curve matched with
   ! parabolas at the top and bottom of the index range.
   ! ---------------------------------------------------------------------------

      type(sigref_fun_spec_type), intent(in) :: fun_spec
      integer, intent(in) :: kmax

      real(r8), dimension(kmax) :: sigref

      real(r8) :: sp2, sp3, zp1, zp4, az, bz, cz, dz, as, bs, cs, ds, &
                  t, z, f0, ft, dft, q1, q2, a, b, c
      integer :: ktt, ktb, k

      ! Obtain end and start k-indicies for the top and bottom parabolas
      ! respectively.
      if (fun_spec%z_top > z_eps) then
         ktt = int(fun_spec%z_top*(kmax - 1)) + 1
      else
         ktt = 0
      endif
      if (fun_spec%z_bot < 1._r8 - z_eps) then
         ktb = int(fun_spec%z_bot*(kmax - 1)) + 2
      else
         ktb = kmax + 1
      endif

      ! Obtain the missing Bezier points from the input parameters. 
      zp1 = 0._r8
      sp2 = fun_spec%sp4 - fun_spec%dsdz_bot*(1._r8 - fun_spec%zp2)
      sp3 = fun_spec%sp4 - fun_spec%dsdz_bot*(1._r8 - fun_spec%zp3)
      zp4 = 1._r8

      ! Compute reference sigma values from the cubic Bezier curve. The Bezier
      ! curve function argument t corresponding to layer index k is found by
      ! solving a cubic polynomial by Newton-Raphson iterative method.

      az = -       zp1 + 3._r8*fun_spec%zp2 - 3._r8*fun_spec%zp3 + zp4
      bz =   3._r8*zp1 - 6._r8*fun_spec%zp2 + 3._r8*fun_spec%zp3
      cz = - 3._r8*zp1 + 3._r8*fun_spec%zp2

      as = -       fun_spec%sp1 + 3._r8*sp2 - 3._r8*sp3 + fun_spec%sp4
      bs =   3._r8*fun_spec%sp1 - 6._r8*sp2 + 3._r8*sp3
      cs = - 3._r8*fun_spec%sp1 + 3._r8*sp2
      ds =         fun_spec%sp1

      t = 0._r8
      do k = ktt + 1, ktb - 1
         z = real(k - 1, r8)/real(kmax - 1, r8)
         dz = zp1 - z
         t = cubic_root(az, bz, cz, dz, t, t_tol)
         sigref(k) = ((as*t + bs)*t + cs)*t + ds
      enddo

      ! Compute reference sigma values in the range [0, z_top] from a parabola
      ! defined by specified sigma value s_top at z = 0 and matching the Bezier
      ! curve and its derivative with respect to z at z = z_top.
      if (ktt > 0) then
         dz = zp1 - fun_spec%z_top
         t = cubic_root(az, bz, cz, dz, 0._r8, t_tol)
         f0 = fun_spec%s_top
         ft = ((as*t + bs)*t + cs)*t + ds
         dft = ((3._r8*as*t + 2._r8*bs)*t + cs)/((3._r8*az*t + 2._r8*bz)*t + cz)
         q1 = 1._r8/fun_spec%z_top
         q2 = (f0 - ft)*q1
         a = (dft + q2)*q1
         b = - (dft + 2._r8*q2)
         c = f0
         do k = 1, ktt
            z = real(k - 1, r8)/real(kmax - 1, r8)
            sigref(k) = (a*z + b)*z + c
         enddo
      endif

      ! Compute reference sigma values in the range [z_bot, 1] from a parabola
      ! defined by specified sigma value s_bot at z = 1 and matching the Bezier
      ! curve and its derivative with respect to z at z = z_bot.
      if (ktb <= kmax) then
         dz = zp1 - fun_spec%z_bot
         t = cubic_root(az, bz, cz, dz, 1._r8, t_tol)
         f0 = fun_spec%s_bot
         ft = ((as*t + bs)*t + cs)*t + ds
         dft = ((3._r8*as*t + 2._r8*bs)*t + cs)/((3._r8*az*t + 2._r8*bz)*t + cz)
         q1 = 1._r8/(1._r8 - fun_spec%z_bot)
         q1 = q1*q1
         a = ((fun_spec%z_bot - 1._r8)*dft + f0 - ft)*q1
         b = (- (dft*fun_spec%z_bot + 2._r8*(f0 - ft))*fun_spec%z_bot + dft)*q1
         c = (((f0 + dft)*fun_spec%z_bot - 2._r8*ft - dft)*fun_spec%z_bot + ft)*q1
         do k = ktb, kmax
            z = real(k - 1, r8)/real(kmax - 1, r8)
            sigref(k) = (a*z + b)*z + c
         enddo
      endif

   end function sigref_fun

   function sra_cost(plevel_test, sigref_test) result(cost)
   ! ---------------------------------------------------------------------------
   ! Return the cost that measures the deviation from mixed layer pressure
   ! thickness (dpml) and the pressure pressure thickness occupied by constant
   ! pressure levels (dpml_plev). The norm used is:
   !   log(dpml_plev/dpml)**2*sra_cost_wgt
   ! ---------------------------------------------------------------------------

      real(r8), dimension(kdm), intent(in) :: plevel_test, sigref_test

      real(r8) :: cost

      real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: cost_2d
      real(r8) :: dpml, sigmlb, w, dpml_plev, logdiff
      integer :: i, j, k, tlev, l

      cost_2d(:,:) = 0._r8

      do tlev = 1, sra_tlev_num
         do j = 1, jj
            do l = 1, isp(j)
            do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
               dpml = sra_dpml_clim(i,j,tlev)
               if (dpml /= spval) then
                  sigmlb = sra_sigmlb_clim(i,j,tlev)
                  k = 2
                  do while (k <= kdm)
                     if (sigmlb < sigref_test(k)) then
                        w =  (sigmlb         - sigref_test(k-1)) &
                            /(sigref_test(k) - sigref_test(k-1))
                        dpml_plev = (1._r8 - w)*plevel_test(k-1) &
                                  +          w *plevel_test(k  )
                        logdiff = log(dpml_plev/dpml)
                        cost_2d(i,j) = cost_2d(i,j) &
                                     + logdiff*logdiff*sra_cost_wgt(i,j)
                        exit
                     endif
                     k = k + 1
                  enddo
               endif
            enddo
            enddo
         enddo
      enddo

      call xcsum(cost, cost_2d, ips)
      
   end function sra_cost

   function sra_cost_grad(plevel_test, sigref_fun_spec_base, x, dx) &
      result(cost_grad)
   ! ---------------------------------------------------------------------------
   ! Estimate the cost gradient with respect to sp1 and zp2 of the reference
   ! potential density function specification.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(kdm), intent(in) :: plevel_test
      type(sigref_fun_spec_type), intent(in) :: sigref_fun_spec_base
      real(r8), dimension(2), intent(in) :: x, dx

      real(r8), dimension(2) :: cost_grad

      type(sigref_fun_spec_type) :: sigref_fun_spec_test
      real(r8) :: cost_m, cost_p

      sigref_fun_spec_test = sigref_fun_spec_base

      sigref_fun_spec_test%zp2 = x(2)
      sigref_fun_spec_test%sp1 = x(1) - .5_r8*dx(1)
      cost_m  = sra_cost(plevel_test, sigref_fun(sigref_fun_spec_test, kdm))
      sigref_fun_spec_test%sp1 = x(1) + .5_r8*dx(1)
      cost_p  = sra_cost(plevel_test, sigref_fun(sigref_fun_spec_test, kdm))
      cost_grad(1) = (cost_p - cost_m)/dx(1)

      sigref_fun_spec_test%sp1 = x(1)
      sigref_fun_spec_test%zp2 = x(2) - .5_r8*dx(2)
      cost_m  = sra_cost(plevel_test, sigref_fun(sigref_fun_spec_test, kdm))
      sigref_fun_spec_test%zp2 = x(2) + .5_r8*dx(2)
      cost_p  = sra_cost(plevel_test, sigref_fun(sigref_fun_spec_test, kdm))
      cost_grad(2) = (cost_p - cost_m)/dx(2)

   end function sra_cost_grad

   subroutine sra_update(m, n, mm, nn, k1m, k1n)
   ! ---------------------------------------------------------------------------
   ! Update reference potential densities.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      real(r8) :: wgt_tf1, wgt_tf2, sp1_tf1, zp2_tf1, sp4_tf1, s_bot_tf1
      integer :: i, j, k

      ! Time filter weights.
      wgt_tf1 = ( real(nday_of_year - 1, r8) &
                + real(mod(nstep, nstep_in_day), r8)/real(nstep_in_day, r8)) &
                /real(nday_in_year, r8)
      wgt_tf2 = baclin/(86400._r8*real(nday_in_year, r8)*sra_param_ts + baclin)

      ! Create signals, which vary linearly from old to new parameter values
      ! over a year, for the final time filter.
      sp1_tf1   = (1._r8 - wgt_tf1)*sigref_fun_spec_old%sp1 &
                +          wgt_tf1 *sigref_fun_spec_new%sp1
      zp2_tf1   = (1._r8 - wgt_tf1)*sigref_fun_spec_old%zp2 &
                +          wgt_tf1 *sigref_fun_spec_new%zp2
      sp4_tf1   = (1._r8 - wgt_tf1)*sigref_fun_spec_old%sp4 &
                +          wgt_tf1 *sigref_fun_spec_new%sp4
      s_bot_tf1 = (1._r8 - wgt_tf1)*sigref_fun_spec_old%s_bot &
                +          wgt_tf1 *sigref_fun_spec_new%s_bot

      ! Apply final time filter.
      sigref_fun_spec%sp1   = (1._r8 - wgt_tf2)*sigref_fun_spec%sp1 &
                            +          wgt_tf2 *sp1_tf1
      sigref_fun_spec%zp2   = (1._r8 - wgt_tf2)*sigref_fun_spec%zp2 &
                            +          wgt_tf2 *zp2_tf1
      sigref_fun_spec%sp4   = (1._r8 - wgt_tf2)*sigref_fun_spec%sp4 &
                            +          wgt_tf2 *sp4_tf1
      sigref_fun_spec%s_bot = (1._r8 - wgt_tf2)*sigref_fun_spec%s_bot &
                            +          wgt_tf2 *s_bot_tf1

      ! Obtain updated reference potential densities.
      sigref(1:kdm) = sigref_fun(sigref_fun_spec, kdm)
      !$omp parallel do private(i, k)
      do j = 1-nbdy, jj+nbdy
         do k = 1, kk
            do i = 1-nbdy, ii+nbdy
               sigmar(i,j,k) = sigref(k)
            enddo
         enddo
      enddo
      !$omp end parallel do

      if (csdiag) then
         if (mnproc == 1) then
            write (lp,*) 'sra_update:'
         endif
         call chksummsk(sigmar, ip, kk, 'sigmar')
      endif

   end subroutine sra_update

   subroutine sra_find_ml_dmax(m, n, mm, nn, k1m, k1n)
   ! ---------------------------------------------------------------------------
   ! Find pressure thickness and potential density at daily maximum mixed layer
   ! depth.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      real(r8) :: pml, pup, plo, sup, slo
      integer :: i, j, k, l, km

      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
            if (sra_dpml_dmax(i,j) /= spval .and. &
                dpml(i,j) > sra_dpml_dmax(i,j)) then
               sra_dpml_dmax(i,j) = dpml(i,j)
               pml = p(i,j,1) + dpml(i,j)
               pup = p(i,j,1) + .5_r8*dp(i,j,1)
               sup = sigma(i,j,k1m)
               k = 2
               km = k + mm
               do
                  if (dp(i,j,km) > onecm) then
                     plo = p(i,j,k) + .5_r8*dp(i,j,km)
                     slo = sigma(i,j,km)
                     if (pml <= plo) then
                        sra_sigmlb_dmax(i,j) = ( slo*(pml - pup) &
                                               + sup*(plo - pml)) &
                                               /(plo - pup)
                        exit
                     endif
                     pup = plo
                     sup = slo
                  endif
                  k = k + 1
                  if (k > kk) then
                     sra_dpml_dmax(i,j) = spval
                     sra_sigmlb_dmax(i,j) = spval
                     exit
                  endif
                  km = k + mm
               enddo
            endif
         enddo
         enddo
      enddo

      if (csdiag) then
         if (mnproc == 1) then
            write (lp,*) 'sra_find_ml_dmax:'
         endif
         call chksummsk(sra_dpml_dmax, ip, 1, 'sra_dpml_dmax')
         call chksummsk(sra_sigmlb_dmax, ip, 1, 'sra_sigmlb_dmax')
      endif

   end subroutine sra_find_ml_dmax

   subroutine sra_accumulate(m, n, mm, nn, k1m, k1n)
   ! ---------------------------------------------------------------------------
   ! Daily accumulation of mixed layer properties and mass distribution in
   ! density classes.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      real(r8) :: mass
      integer :: i, j, k, l, tlev, km, kdc

      ! ------------------------------------------------------------------------
      ! Accumulate daily maximum mixed layer depth and associated potential
      ! density at mixed layer base.
      ! ------------------------------------------------------------------------

      ! Find time level to accumulate in.
      if (nday_of_year == 1) then
         tlev = sra_tlev_num
      else
         tlev = int(real((nday_of_year - 2)*sra_tlev_num, r8) &
                    /real(nday_in_year, r8)) + 1
      endif

      sra_tlev_accnum(tlev) = sra_tlev_accnum(tlev) + 1
      do j = 1, jj
         do l = 1, isp(j)
         do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
            if (sra_dpml_sum(i,j,tlev) /= spval) then
               if (sra_dpml_dmax(i,j) /= spval) then
                  sra_dpml_sum(i,j,tlev) = sra_dpml_sum(i,j,tlev) &
                                         + sra_dpml_dmax(i,j)
                  sra_sigmlb_sum(i,j,tlev) = sra_sigmlb_sum(i,j,tlev) &
                                           + sra_sigmlb_dmax(i,j)
               else
                  sra_dpml_sum(i,j,tlev) = spval
                  sra_sigmlb_sum(i,j,tlev) = spval
               endif
            endif
         enddo
         enddo
      enddo

      ! Reset daily maximum variables.
      sra_dpml_dmax(:,:) = 0._r8
      sra_sigmlb_dmax(:,:) = spval

      ! ------------------------------------------------------------------------
      ! Accumulate mass distribution in density classes.
      ! ------------------------------------------------------------------------

      do k = 1, kk
         km = k + mm
         do j = 1, jj
            do l = 1, isp(j)
            do i = max(1, ifp(j,l)), min(ii, ilp(j,l))

               ! Find k-index of density class.
               kdc = k
               do
                 if (sigma(i,j,km) >= sigref(kdc)) then
                    if (kdc == kdm) exit
                    if (sigma(i,j,km) < sigref(kdc+1)) exit
                    kdc = kdc + 1
                 else
                    if (kdc == 1) exit
                    kdc = kdc - 1
                 endif
               enddo

               ! Accumulate column sums.
               mass = dp(i,j,km)*scp2(i,j)/grav
               sra_massgs_colsum(i,j) = sra_massgs_colsum(i,j) + mass
               sra_massdc_colsum(i,j,kdc) = sra_massdc_colsum(i,j,kdc) + mass
               sra_sigmassdc_colsum(i,j,kdc) = sra_sigmassdc_colsum(i,j,kdc) &
                                             + sigma(i,j,km)*mass

            enddo
            enddo
         enddo
      enddo

      ! ------------------------------------------------------------------------
      ! Accumulate potential density parameters.
      ! ------------------------------------------------------------------------

      sra_accnum = sra_accnum + 1
      sra_sigref_sum(:) = sra_sigref_sum(:) + sigref(1:kdm)
      sra_s_bot_sum = sra_s_bot_sum + sigref_fun_spec%s_bot

      if (csdiag) then
         if (mnproc == 1) then
            write (lp,*) 'sra_accumulate:'
         endif
         call chksummsk(sra_dpml_sum, ip, sra_tlev_num, 'sra_dpml_sum')
         call chksummsk(sra_sigmlb_sum, ip, sra_tlev_num, 'sra_sigmlb_sum')
         call chksummsk(sra_massgs_colsum, ip, 1, 'sra_massgs_colsum')
         call chksummsk(sra_massdc_colsum, ip, kk, 'sra_massdc_colsum')
         call chksummsk(sra_sigmassdc_colsum, ip, kk, 'sra_sigmassdc_colsum')
      endif

   end subroutine sra_accumulate

   subroutine sra_optimize
   ! ---------------------------------------------------------------------------
   ! Optimize parameters for the potential density function based on accumulated
   ! properties.
   ! ---------------------------------------------------------------------------

      real(r8), parameter :: &
         adam_alpha = .01_r8, &
         adam_beta1 = .9_r8, &
         adam_beta2 = .999_r8, &
         adam_eps   = 1.e-8_r8
      integer, parameter :: &
         adam_maxiter = 500

      real(r8), dimension(kdm) :: massfracdc, sigdc, sigref_mean
      real(r8), dimension(2) :: adam_m, adam_v, adam_mhat, adam_vhat, &
                                x, dx, cost_grad
      real(r8) :: wgt_tf, q, massgs, massdc, sigmassdc, s_bot_mean, rktb, &
                  massfrac_bot, sp4_new, s_bot_new, cost, adam_beta1pt, &
                  adam_beta2pt
      integer :: i, j, l, tlev, kdc, ktb
      type(sigref_fun_spec_type) :: sigref_fun_spec_test

      ! Copy reference potential density function specifications from new to
      ! old.
      sigref_fun_spec_old = sigref_fun_spec_new

      ! Obtain time-level climatology of mixed layer depth and associated
      ! potential density at mixed layer base.
      wgt_tf = 1._r8/(sra_clim_ts + 1._r8)
      do tlev = 1, sra_tlev_num
         q = 1./real(sra_tlev_accnum(tlev), r8)
         do j = 1, jj
            do l = 1, isp(j)
            do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
               if (sra_dpml_sum(i,j,tlev) /= spval) then
                  if (sra_dpml_clim(i,j,tlev) /= spval) then
                     sra_dpml_clim(i,j,tlev) = &
                        (1._r8 - wgt_tf)*sra_dpml_clim(i,j,tlev) &
                      +          wgt_tf *sra_dpml_sum(i,j,tlev)*q
                     sra_sigmlb_clim(i,j,tlev) = &
                        (1._r8 - wgt_tf)*sra_sigmlb_clim(i,j,tlev) &
                      +          wgt_tf *sra_sigmlb_sum(i,j,tlev)*q
                  else
                     sra_dpml_clim(i,j,tlev) = sra_dpml_sum(i,j,tlev)*q
                     sra_sigmlb_clim(i,j,tlev) = sra_sigmlb_sum(i,j,tlev)*q
                  endif
               endif
            enddo
            enddo
         enddo
      enddo

      ! Reset sums.
      sra_tlev_accnum(:) = 0
      sra_dpml_sum(:,:,:) = 0._r8
      sra_sigmlb_sum(:,:,:) = 0._r8

      ! Obtain mass distribution in density classes by summing up column sums
      ! globally and making global means.
      call xcsum(massgs, sra_massgs_colsum, ips)
      do kdc = 1, kdm
         call xcsum(massdc, sra_massdc_colsum(1-nbdy,1-nbdy,kdc), ips)
         call xcsum(sigmassdc, sra_sigmassdc_colsum(1-nbdy,1-nbdy,kdc), ips)
         massfracdc(kdc) = massdc/massgs
         if (massdc > 0._r8) then
            sigdc(kdc) = sigmassdc/massdc
         else
            sigdc(kdc) = spval
         endif
      enddo

      ! Reset sums.
      sra_massgs_colsum(:,:) = 0._r8
      sra_massdc_colsum(:,:,:) = 0._r8
      sra_sigmassdc_colsum(:,:,:) = 0._r8

      ! Obtain time mean of potential density parameters.
      q = 1._r8/real(sra_accnum, r8)
      sigref_mean(:) = sra_sigref_sum(:)*q
      s_bot_mean = sra_s_bot_sum*q

      ! Reset sums.
      sra_accnum = 0
      sra_sigref_sum(:) = 0._r8
      sra_s_bot_sum = 0._r8

      ! ------------------------------------------------------------------------
      ! Adjust parameters related to the maximum reference potential density and
      ! distribution of the densest reference potential densities.
      ! ------------------------------------------------------------------------

      if (sigref_fun_spec%z_bot < 1._r8 - z_eps) then

         ! Diagnose the mass fraction in the range [z_bot 1]
         rktb = sigref_fun_spec%z_bot*real(kdm - 1, r8) + 1._r8
         ktb = int(rktb) + 1
         massfrac_bot = massfracdc(ktb-1)*(real(ktb, r8) - rktb)
         do kdc = ktb, kdm
            massfrac_bot = massfrac_bot + massfracdc(kdc)
         enddo
         if (mnproc == 1) &
            write(lp,*) 'sra_optimize: massfrac_bot:', massfrac_bot

         ! Adjust sp4 so that the mass fraction in the range [z_bot 1]
         ! approaches sra_massfrac_bot.

         massfrac_bot = 0._r8
         kdc = kdm + 1
         do
            kdc = kdc - 1
            if (massfrac_bot + massfracdc(kdc) > sra_massfrac_bot) then
               if     (kdc == kdm) then
                  sp4_new = &
                     sigref_mean(kdc) &
                   + (1._r8 - sigref_fun_spec%z_bot)*sigref_fun_spec%dsdz_bot
               elseif (massfracdc(kdc) < sra_massfrac_eps) then
                  sp4_new = &
                     .5_r8*(sigref_mean(kdc) + sigref_mean(kdc+1)) &
                   + (1._r8 - sigref_fun_spec%z_bot)*sigref_fun_spec%dsdz_bot
               else
                  q = (sra_massfrac_bot - massfrac_bot)/massfracdc(kdc)
                  sp4_new = &
                     sigref_mean(kdc)*q + sigref_mean(kdc+1)*(1._r8 - q) &
                   + (1._r8 - sigref_fun_spec%z_bot)*sigref_fun_spec%dsdz_bot
               endif
               exit
            else
               massfrac_bot = massfrac_bot + massfracdc(kdc)
            endif
         enddo

         ! Adjust s_bot.

         if (massfracdc(kdm) < sra_massfrac_eps) then
            ! If the mass fraction of the densest layer is below
            ! sra_massfrac_eps, set s_bot to the mean reference density of the
            ! densest layer with mass fraction above sra_massfrac_eps.
            kdc = kdm - 1
            do while (massfracdc(kdc) < sra_massfrac_eps)
               kdc = kdc - 1
            enddo
            s_bot_new = sigref_mean(kdc)
         else
            ! Adjust s_bot to balance the mass fraction of the two densest
            ! layers.
            if (massfracdc(kdm-1) > massfracdc(kdm)) then
               s_bot_new = &
                  s_bot_mean &
               - .5_r8*(massfracdc(kdm-1) - massfracdc(kdm)) &
                      *(s_bot_mean - sigref_mean(kdm-1))/massfracdc(kdm-1)
            else
               s_bot_new = &
                  s_bot_mean &
                + (massfracdc(kdm) - massfracdc(kdm-1)) &
                  *(sigdc(kdm) - s_bot_mean)/massfracdc(kdm)
            endif
         endif
         s_bot_new = max(s_bot_new, sp4_new)

         if (mnproc == 1) then
            write(lp,*) 'sra_optimize: sp4   old/new:', &
                        sigref_fun_spec_old%sp4, sp4_new
            write(lp,*) 'sra_optimize: s_bot old/new:', &
                        sigref_fun_spec_old%s_bot, s_bot_new
         endif

      endif

      ! ------------------------------------------------------------------------
      ! Optimize parameters so that the difference between the range of constant
      ! pressure levels and simulated mixed layer depth is minimized.
      ! ------------------------------------------------------------------------

      cost = sra_cost(plevel, sigref_fun(sigref_fun_spec_old, kdm))
      if (mnproc == 1) &
         write(lp,'(a,f15.7)') ' sra_optimize: cost prev. optim. sigref:', cost
      cost = sra_cost(plevel, sigref)
      if (mnproc == 1) &
         write(lp,'(a,f15.7)') ' sra_optimize: cost current sigref:     ', cost

      ! Use ADAM optimizer to find new sp1 and zp2.

      adam_m(:) = 0._r8
      adam_v(:) = 0._r8
      adam_beta1pt = 1._r8
      adam_beta2pt = 1._r8

      x = [sigref_fun_spec%sp1, sigref_fun_spec%zp2]
      dx = [1.e-6_r8, 1.e-6_r8]
      sigref_fun_spec_test = sigref_fun_spec

      do i = 1, adam_maxiter

         cost_grad(:) = sra_cost_grad(plevel, sigref_fun_spec, x, dx)
         adam_m(:) = adam_beta1*adam_m(:) &
                   + (1._r8 - adam_beta1)*cost_grad(:)
         adam_v(:) = adam_beta2*adam_v(:) &
                   + (1._r8 - adam_beta2)*(cost_grad(:)*cost_grad(:))
         adam_beta1pt = adam_beta1pt*adam_beta1
         adam_beta2pt = adam_beta2pt*adam_beta2
         adam_mhat(:) = adam_m(:)/(1._r8 - adam_beta1pt)
         adam_vhat(:) = adam_v(:)/(1._r8 - adam_beta2pt)
         x(:) = x(:) - adam_alpha*adam_mhat(:)/(sqrt(adam_vhat(:)) + adam_eps)

         if ( mod(i, 100) == 0) then
            sigref_fun_spec_test%sp1 = x(1)
            sigref_fun_spec_test%zp2 = x(2)
            cost = sra_cost(plevel, sigref_fun(sigref_fun_spec_test, kdm))
            if (mnproc == 1) &
               write(lp,'(a,3e15.7)') ' sra_optimize: sp1, zp2, cost:', &
                                      x(1), x(2), cost
         endif

      enddo

      sigref_fun_spec_new       = sigref_fun_spec
      sigref_fun_spec_new%sp1   = x(1)
      sigref_fun_spec_new%zp2   = x(2)
      sigref_fun_spec_new%sp4   = sp4_new
      sigref_fun_spec_new%s_bot = s_bot_new

      if (mnproc == 1) then
        write(lp,*) 'sra_optimize: optimized sigref_fun_spec:'
        write(lp,*) '  dsdz_bot:',sigref_fun_spec_new%dsdz_bot
        write(lp,*) '  sp1:     ',sigref_fun_spec_new%sp1
        write(lp,*) '  zp2:     ',sigref_fun_spec_new%zp2
        write(lp,*) '  zp3:     ',sigref_fun_spec_new%zp3
        write(lp,*) '  sp4:     ',sigref_fun_spec_new%sp4
        write(lp,*) '  z_top:   ',sigref_fun_spec_new%z_top
        write(lp,*) '  s_top:   ',sigref_fun_spec_new%s_top
        write(lp,*) '  z_bot:   ',sigref_fun_spec_new%z_bot
        write(lp,*) '  s_bot:   ',sigref_fun_spec_new%s_bot
      endif

   end subroutine sra_optimize

   ! ---------------------------------------------------------------------------
   ! Public procedures.
   ! ---------------------------------------------------------------------------

   subroutine readnml_vcoord
   ! ---------------------------------------------------------------------------
   ! Read variables in the namelist group 'vcoord' and resolve options.
   ! ---------------------------------------------------------------------------

      character(len = 80) :: nml_fname
      real(r8) :: dpmin
      integer :: nfu, ios, k
      logical :: fexist

      namelist /vcoord/ &
         vcoord_type, dpmin_surface, dpmin_inflation_factor, &
         sigref_spec, plevel_spec, sigref_fun_spec, sigref, plevel, &
         sigref_adaption, sra_clim_ts, sra_param_ts, sra_massfrac_bot, &
         sra_massfrac_eps

      ! Read variables in the namelist group 'vcoord'.
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
         call xcbcst(dpmin_surface)
         call xcbcst(dpmin_inflation_factor)
         call xcbcst(sigref_spec)
         call xcbcst(plevel_spec)
         call xcbcst(sigref_fun_spec%dsdz_bot)
         call xcbcst(sigref_fun_spec%sp1)
         call xcbcst(sigref_fun_spec%zp2)
         call xcbcst(sigref_fun_spec%zp3)
         call xcbcst(sigref_fun_spec%sp4)
         call xcbcst(sigref_fun_spec%z_top)
         call xcbcst(sigref_fun_spec%s_top)
         call xcbcst(sigref_fun_spec%z_bot)
         call xcbcst(sigref_fun_spec%s_bot)
         call xcbcst(sigref)
         call xcbcst(plevel)
         call xcbcst(sigref_adaption)
         call xcbcst(sra_clim_ts)
         call xcbcst(sra_param_ts)
         call xcbcst(sra_massfrac_bot)
         call xcbcst(sra_massfrac_eps)
      endif
      if (mnproc == 1) then
         write (lp,*) 'readnml_vcoord: vertical coordinate variables:'
         write (lp,*) '  vcoord_type =            ', trim(vcoord_type)
         write (lp,*) '  dpmin_surface =          ', dpmin_surface
         write (lp,*) '  dpmin_inflation_factor = ', dpmin_inflation_factor
         write (lp,*) '  sigref_fun_spec =        ', sigref_fun_spec
         write (lp,*) '  sigref_spec =            ', trim(sigref_spec)
         write (lp,*) '  plevel_spec =            ', trim(plevel_spec)
         write (lp,*) '  sigref_adaption =        ', sigref_adaption
         write (lp,*) '  sra_clim_ts =            ', sra_clim_ts
         write (lp,*) '  sra_param_ts =           ', sra_param_ts
         write (lp,*) '  sra_massfrac_bot =       ', sra_massfrac_bot
         write (lp,*) '  sra_massfrac_eps =       ', sra_massfrac_eps
      endif

      ! Change units from [m] to [kg m-1 s-2] of depth interval variables.
      dpmin_surface = dpmin_surface*onem

      ! Resolve options.
      select case (trim(vcoord_type))
         case ('isopyc_bulkml')
            vcoord_tag = vcoord_isopyc_bulkml
         case ('cntiso_hybrid')
            vcoord_tag = vcoord_cntiso_hybrid
         case ('plevel')
            vcoord_tag = vcoord_plevel
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') ' readnml_vcoord: vcoord_type = ', &
                                 trim(vcoord_type), ' is unsupported!'
            call xcstop('(readnml_vcoord)')
                   stop '(readnml_vcoord)'
      end select
      if (vcoord_tag /= vcoord_isopyc_bulkml) then
         select case (trim(sigref_spec))
            case ('inicon')
            case ('function')
               sigref(1:kdm) = sigref_fun(sigref_fun_spec,kdm)
            case ('namelist')
               k = 1
               do while (sigref(k) /= spval)
                  k = k + 1
                  if (k > kdm_max) exit
               enddo
               if (k /= kdm + 1) then
                  if (mnproc == 1) &
                     write (lp,'(3a)') &
                        ' readnml_vcoord: number of sigref values does not match vertical dimension!'
                  call xcstop('(readnml_vcoord)')
                         stop '(readnml_vcoord)'
               endif
            case default
               if (mnproc == 1) &
                  write (lp,'(3a)') ' readnml_vcoord: sigref_spec = ', &
                                    trim(sigref_spec), ' is unsupported!'
               call xcstop('(readnml_vcoord)')
                      stop '(readnml_vcoord)'
         end select
         select case (trim(plevel_spec))
            case ('inflation')
               dpmin = dpmin_surface
               plevel(1) = 0._r8
               do k = 1, kk - 1
                  plevel(k+1) = plevel(k) + dpmin
                  dpmin = dpmin*dpmin_inflation_factor
               enddo
            case ('namelist')
               k = 1
               do while (plevel(k) /= spval)
                  k = k + 1
                  if (k > kdm_max) exit
               enddo
               if (k /= kdm + 1) then
                  if (mnproc == 1) &
                     write (lp,'(3a)') &
                        ' readnml_vcoord: number of plevel values does not match vertical dimension!'
                  call xcstop('(readnml_vcoord)')
                  stop '(readnml_vcoord)'
               endif
               ! Change units from [m] to [kg m-1 s-2].
               plevel(:) = plevel(:)*onem
            case default
               if (mnproc == 1) &
                  write (lp,'(3a)') ' readnml_vcoord: plevel_spec = ', &
                                    trim(plevel_spec), ' is unsupported!'
               call xcstop('(readnml_vcoord)')
               stop '(readnml_vcoord)'
         end select
      endif

   end subroutine readnml_vcoord

   subroutine inivar_vcoord
   ! ---------------------------------------------------------------------------
   ! Initialize arrays and data structures.
   ! ---------------------------------------------------------------------------

      integer :: i, j, k, l

      if (vcoord_tag == vcoord_isopyc_bulkml .or. &
          trim(sigref_spec) == 'inicon') then
         sigmar(:,:,:) = spval
      else
         !$omp parallel do private(i, k)
         do j = 1-nbdy, jj+nbdy
            do k = 1, kk
               do i = 1-nbdy, ii+nbdy
                  sigmar(i,j,k) = sigref(k)
               enddo
            enddo
         enddo
         !$omp end parallel do
         if (sigref_adaption) then
            sigref_fun_spec_old = sigref_fun_spec
            sigref_fun_spec_new = sigref_fun_spec
            sra_dpml_dmax(:,:) = 0._r8
            sra_sigmlb_dmax(:,:) = spval
            sra_tlev_accnum(:) = 0
            sra_dpml_sum(:,:,:) = 0._r8
            sra_sigmlb_sum(:,:,:) = 0._r8
            sra_dpml_clim(:,:,:) = spval
            sra_sigmlb_clim(:,:,:) = spval
            sra_massgs_colsum(:,:) = 0._r8
            sra_massdc_colsum(:,:,:) = 0._r8
            sra_sigmassdc_colsum(:,:,:) = 0._r8
            sra_accnum = 0
            sra_sigref_sum(:) = 0._r8
            sra_s_bot_sum = 0._r8
            do j = 1, jj
               do l = 1, isp(j)
               do i = max(1, ifp(j,l)), min(ii, ilp(j,l))
                  sra_cost_wgt(i,j) = scp2(i,j)/(area*real(sra_tlev_num, r8))
               enddo
               enddo
            enddo
         endif
      endif

   end subroutine inivar_vcoord

   subroutine extract_sigref
   ! ---------------------------------------------------------------------------
   ! In case it is not otherwise specified, extract reference potential density
   ! vector representative of the dominating ocean domain.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(itdm,jtdm) :: tmp2d
      integer :: i, j, k, i1, j1
      logical :: sigref_found

      if (vcoord_tag == vcoord_isopyc_bulkml .or. &
          trim(sigref_spec) == 'inicon') then
         !$omp parallel do private(i)
         do j = 1, jj
           do i = 1, ii
             util1(i,j) = real(ipwocn(i,j), r8)
           enddo
         enddo
         !$omp end parallel do
         call xcaget(tmp2d, util1, 1)
         if (mnproc == 1) then
            sigref_found = .false.
            do j = 1, jtdm
               do i = 1, itdm
                 if (tmp2d(i,j) > 0._r8) then
                    i1 = i
                    j1 = j
                    sigref_found = .true.
                    exit
                 endif
              enddo
              if (sigref_found) exit
            enddo
         endif
         call xcbcst(i1)
         call xcbcst(j1)
         do k = 1, kk
            call xceget(sigref(k),sigmar(1-nbdy,1-nbdy,k),i1,j1)
         enddo
      endif
      if (mnproc == 1) then
         write(lp,*) 'sigma layers = ',sigref(1:kk)
      endif

   end subroutine extract_sigref

   subroutine sigref_adapt(m, n, mm, nn, k1m, k1n)

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      if (.not. sigref_adaption) return

      ! Update reference potential densities.
      call sra_update(m, n, mm, nn, k1m, k1n)

      ! Find pressure thickness and potential density at daily maximum mixed
      ! layer depth.
      call sra_find_ml_dmax(m, n, mm, nn, k1m, k1n)

      if (mod(nstep, nstep_in_day) == 0) then

         ! Daily accumulation of mixed layer properties and mass distribution in
         ! density classes.
         call sra_accumulate(m, n, mm, nn, k1m, k1n)

         if (nday_of_year == 1) then

            ! At the end of a model year, optimize parameters for the
            ! potential density function based on accumulated properties.
            call sra_optimize

         endif

      endif

   end subroutine sigref_adapt

end module mod_vcoord
