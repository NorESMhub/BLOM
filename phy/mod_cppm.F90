! ------------------------------------------------------------------------------
! Copyright (C) 2025 Mats Bentsen
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

module mod_cppm
! ------------------------------------------------------------------------------
! This module contains routines for advection of layer pressure thickness and
! tracers using Compatible Piecewise Parabolic Method (CPPM). If monotonic
! limiting is selected, the method ensures monotonicity of layer thickness for
! non-divergent flow and monotonicity of tracers for all flows.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_time, only: nstep
   use mod_xc, only: idm, jdm, kdm, nbdy, ii, jj, ip, &
                     halo_ps, halo_us, halo_vs, halo_uv, halo_vv, &
                     lp, mnproc, nproc, jpr, nreg, itdm, i0, &
                     xctilr, xcstop
   use mod_grid, only: scpx, scpy, scp2i
   use mod_state, only: dp, temp, saln, uflx, vflx, utflx, vtflx, &
                        usflx, vsflx, p, cau, cav, pbu, pbv
   use mod_tracers, only: trc, ntr

   implicit none

   private

   ! Options with default values, modifiable by namelist.
   character(len = 80) :: &
      cppm_compatibility = 'full', &    ! CPPM compatibility. Valid options:
                                        ! 'full', 'partial'.
      cppm_limiting = 'non_oscillatory' ! CPPM limiting. Valid methods:
                                        ! 'monotonic', 'non_oscillatory'.

   ! Options derived from string options.
   integer :: &
      cppm_compatibility_tag, &
      cppm_limiting_tag

   integer, parameter :: &
      cppm_compatibility_full       = 1, & ! Full compatibility.
      cppm_compatibility_partial    = 2, & ! Partial compatibility.
      cppm_limiting_monotonic       = 1, & ! Monotonic limiting.
      cppm_limiting_non_oscillatory = 2, & ! Non-oscillatory limiting.
      stencil_0000 = 0, &
      stencil_1111 = 1, &
      stencil_1110 = 2, &
      stencil_0111 = 3, &
      stencil_1100 = 4, &
      stencil_0110 = 5, &
      stencil_0011 = 6, &
      stencil_0100 = 7, &
      stencil_0010 = 8

   real(r8), parameter :: &
      c0 = 0._r8, c1 = 1._r8, c2 = 2._r8, c3 = 3._r8, c4 = 4._r8, c5 = 5._r8, &
      c6 = 6._r8, c12 = 12._r8, c18 = 18._r8, c42 = 42._r8, c60 = 60._r8, &
      c1_2 =1._r8/2._r8, c1_3 = 1._r8/3._r8, c2_3 = 2._r8/3._r8, &
      c1_4 = 1._r8/4._r8, c3_4 = 3._r8/4._r8, c1_5 = 1._r8/5._r8, &
      c1_6 = 1._r8/6._r8, c1_10 = 1._r8/10._r8, c1_12 = 1._r8/12._r8, & 
      c1_15 = 1._r8/15._r8, c1_20 = 1._r8/20._r8, &
      dpeps = 1.e-12_r8

   integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: stencili
   integer, dimension(1-nbdy:jdm+nbdy,1-nbdy:idm+nbdy) :: stencilj
   real(r8), dimension(12,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
      tmc0i, tmcli, tmcri
   real(r8), dimension(12,1-nbdy:jdm+nbdy,1-nbdy:idm+nbdy) :: &
      tmc0j, tmclj, tmcrj
   real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
      hevc1i, hevc2i, hevc3i, hevc4i, ssci, scci, d2mi
   real(r8), dimension(1-nbdy:jdm+nbdy,1-nbdy:idm+nbdy) :: &
      hevc1j, hevc2j, hevc3j, hevc4j, sscj, sccj, d2mj
   real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: hel_3d, her_3d

   integer :: ntr_loc

   public :: cppm_compatibility, cppm_limiting, init_cppm, cppm

contains

   ! ---------------------------------------------------------------------------
   ! Private procedures.
   ! ---------------------------------------------------------------------------

   pure subroutine set_stencil_coeffs(stencil_mask, dx, stencil, &
                                      hevc1, hevc2, hevc3, hevc4, &
                                      tmc0, tmcl, tmcr)
   ! ---------------------------------------------------------------------------
   ! Set stencil tag and, accounting for the presence of eligible cell mean
   ! values in the 4-cell stencil, set coefficients for thickness edge value
   ! reconstruction and cell width dependent coefficients for the matrix
   ! elements used for solving for tracer edge value coefficients.
   ! ---------------------------------------------------------------------------

      integer, dimension(4), intent(in) :: stencil_mask
      real(r8), dimension(4), intent(in) :: dx
      integer, intent(out) :: stencil
      real(r8), intent(out) :: hevc1, hevc2, hevc3, hevc4
      real(r8), dimension(12), intent(out) :: tmc0, tmcl, tmcr

      real(r8) :: a12, a22, a32, a42, a13, a23, a33, a43, a14, a24, a34, a44

      ! Define elements of the matrix A for the linear system of equations that
      ! relates thickness edge values to cell means for the 1111 stencil.
      a12 = - dx(2) - c1_2*dx(1)
      a22 = - c1_2*dx(2)
      a32 =   c1_2*dx(3)
      a42 =   dx(3) + c1_2*dx(4)
      a13 =   a12*a12 + c1_12*dx(1)*dx(1)
      a23 = - c2_3*a22*dx(2)
      a33 =   c2_3*a32*dx(3)
      a43 =   a42*a42 + c1_12*dx(4)*dx(4)
      a14 =   (a13 + c1_6*dx(1)*dx(1))*a12
      a24 = - c3_4*a23*dx(2)
      a34 =   c3_4*a33*dx(3)
      a44 =   (a43 + c1_6*dx(4)*dx(4))*a42

      ! Cell width dependent coefficients for the matrix elements used for
      ! solving for tracer edge value coefficients.

      tmcl( 1) = - c1_12*dx(1)
      tmcl( 2) =   (c1_10*dx(1) + c1_6*dx(2))*dx(1)
      tmcl( 3) = - (c1_10*(dx(1) + c3*dx(2))*dx(1) + c1_4*dx(2)**2)*dx(1)
      tmcl( 4) = - c1_12*dx(2)
      tmcl( 5) =   c1_10*dx(2)**2
      tmcl( 6) = - c1_10*dx(2)**3
      tmcl( 7) = - c1_12*dx(3)
      tmcl( 8) = - c1_15*dx(3)**2
      tmcl( 9) = - c1_20*dx(3)**3
      tmcl(10) = - c1_12*dx(4)
      tmcl(11) = - (c1_15*dx(4) + c1_6*dx(3))*dx(4)
      tmcl(12) = - (c1_5*(c1_4*dx(4) + dx(3))*dx(4) + c1_4*dx(3)**2)*dx(4)

      tmcr( 1) =   c1_12*dx(1)
      tmcr( 2) = - (c1_15*dx(1) + c1_6*dx(2))*dx(1)
      tmcr( 3) =   (c1_5*(c1_4*dx(1) + dx(2))*dx(1) + c1_4*dx(2)**2)*dx(1)
      tmcr( 4) =   c1_12*dx(2)
      tmcr( 5) = - c1_15*dx(2)**2
      tmcr( 6) =   c1_20*dx(2)**3
      tmcr( 7) =   c1_12*dx(3)
      tmcr( 8) =   c1_10*dx(3)**2
      tmcr( 9) =   c1_10*dx(3)**3
      tmcr(10) =   c1_12*dx(4)
      tmcr(11) =   (c1_10*dx(4) + c1_6*dx(3))*dx(4)
      tmcr(12) =   (c1_10*(dx(4) + c3*dx(3))*dx(4) + c1_4*dx(3)**2)*dx(4)

      tmc0( 1) = a12
      tmc0( 2) = a13 - tmcl( 2) - tmcr( 2)
      tmc0( 3) = a14 - tmcl( 3) - tmcr( 3)
      tmc0( 4) = a22
      tmc0( 5) = a23 - tmcl( 5) - tmcr( 5)
      tmc0( 6) = a24 - tmcl( 6) - tmcr( 6)
      tmc0( 7) = a32
      tmc0( 8) = a33 - tmcl( 8) - tmcr( 8)
      tmc0( 9) = a34 - tmcl( 9) - tmcr( 9)
      tmc0(10) = a42
      tmc0(11) = a43 - tmcl(11) - tmcr(11)
      tmc0(12) = a44 - tmcl(12) - tmcr(12)

      if     (all(stencil_mask == [1, 1, 1, 1])) then
         
         ! Set stencil tag for the 1111 stencil.
         stencil = stencil_1111

         ! LU decomposition of matrix A.
         a22 = a22 - a12
         a32 = a32 - a12
         a42 = a42 - a12
         a23 = (a23 - a13)/a22
         a33 = a33 - a13 - a23*a32
         a43 = a43 - a13 - a23*a42
         a24 = (a24 - a14)/a22
         a34 = a34 - a14 - a24*a32
         a44 = a44 - a14 - a24*a42
         a34 = a34/a33
         a44 = a44 - a34*a43

         ! Forward/backward substitution to obtain coefficients for thickness
         ! edge value reconstruction.
         hevc2 = - a12
         hevc3 = - a13 - a23*hevc2
         hevc4 = - a14 - a24*hevc2 - a34*hevc3
         hevc4 = hevc4/a44
         hevc3 = (hevc3 - a43*hevc4)/a33
         hevc2 = (hevc2 - a32*hevc3 - a42*hevc4)/a22
         hevc1 = c1 - hevc2 - hevc3 - hevc4

      elseif (all(stencil_mask == [1, 1, 1, 0])) then

         ! Set stencil tag for the 1110 stencil.
         stencil = stencil_1110

         ! LU decomposition of the submatrix of A to be used for the 1110
         ! stencil.
         a22 = a22 - a12
         a32 = a32 - a12
         a23 = (a23 - a13)/a22
         a33 = a33 - a13 - a23*a32

         ! Forward/backward substitution to obtain coefficients for thickness
         ! edge value reconstruction.
         hevc2 = - a12
         hevc3 = - a13 - a23*hevc2
         hevc3 = hevc3/a33
         hevc2 = (hevc2 - a32*hevc3)/a22
         hevc1 = c1 - hevc2 - hevc3
         hevc4 = c0

      elseif (all(stencil_mask == [0, 1, 1, 1])) then

         ! Set stencil tag for the 0111 stencil.
         stencil = stencil_0111

         ! LU decomposition of the submatrix of A to be used for the 0111
         ! stencil.
         a32 = a32 - a22
         a42 = a42 - a22
         a33 = (a33 - a23)/a32
         a43 = a43 - a23 - a33*a42

         ! Forward/backward substitution to obtain coefficients for thickness
         ! edge value reconstruction.
         hevc3 = - a22
         hevc4 = - a23 - a33*hevc3
         hevc4 = hevc4/a43
         hevc3 = (hevc3 - a42*hevc4)/a32
         hevc2 = c1 - hevc3 - hevc4
         hevc1 = c0

      elseif (all(stencil_mask == [0, 1, 1, 0])) then

         ! Set stencil tag for the 0110 stencil.
         stencil = stencil_0110

         ! Compute coefficients for thickness edge value reconstruction using a
         ! submatrix of A.
         a32 = a32 - a22
         hevc3 = - a22/a32
         hevc2 = c1 - hevc3
         hevc1 = c0
         hevc4 = c0

      elseif (all(stencil_mask(1:2) == [1, 1])) then

         ! Set stencil tag for the 1100 stencil.
         stencil = stencil_1100

         ! Compute coefficients for thickness edge value reconstruction using a
         ! submatrix of A.
         a22 = a22 - a12
         hevc2 = - a12/a22
         hevc1 = c1 - hevc2
         hevc3 = c0
         hevc4 = c0

      elseif (all(stencil_mask(3:4) == [1, 1])) then

         ! Set stencil tag for the 0011 stencil.
         stencil = stencil_0011

         ! Compute coefficients for thickness edge value reconstruction using a
         ! submatrix of A.
         a42 = a42 - a32
         hevc4 = - a32/a42
         hevc3 = c1 - hevc4
         hevc1 = c0
         hevc2 = c0

      elseif (stencil_mask(2) == 1) then

         ! Set stencil tag for the 0100 stencil.
         stencil = stencil_0100

         ! Coefficients for thickness edge value reconstruction.
         hevc1 = c0
         hevc2 = c1
         hevc3 = c0
         hevc4 = c0

      elseif (stencil_mask(3) == 1) then

         ! Set stencil tag for the 0010 stencil.
         stencil = stencil_0010

         ! Coefficients for thickness edge value reconstruction.
         hevc1 = c0
         hevc2 = c0
         hevc3 = c1
         hevc4 = c0

      else

         ! Set stencil tag for the 0000 stencil.
         stencil = stencil_0000

         ! Coefficients for thickness edge value reconstruction.
         hevc1 = c0
         hevc2 = c0
         hevc3 = c0
         hevc4 = c0

      endif

   end subroutine set_stencil_coeffs

   pure subroutine set_slope_coeffs(stencil_mask, dx, ssc, scc)
   ! ---------------------------------------------------------------------------
   ! Set coefficients for left and right one-sided slopes (ssc) and centered
   ! slopes (scc), accounting for the presence of eligible cell mean values in
   ! the 3-cell stencil.
   ! ---------------------------------------------------------------------------

      integer, dimension(3), intent(in) :: stencil_mask
      real(r8), dimension(3), intent(in) :: dx
      real(r8), intent(out) :: ssc, scc

      if     (any(stencil_mask == [0, 0, 0])) then
         ssc = c0
         scc = c0
      else
         ssc = c2
         scc = c2*dx(2)/(dx(1) + c2*dx(2) + dx(3))
      endif

   end subroutine set_slope_coeffs

   pure subroutine set_d2_mask(stencil_mask, d2m)
   ! ---------------------------------------------------------------------------
   ! Set mask for the value proportional to the second derivative of unlimited
   ! parabolas, accounting for the presence of eligible cell mean values in the
   ! 3-cell stencil.
   ! ---------------------------------------------------------------------------

      integer, dimension(3), intent(in) :: stencil_mask
      real(r8), intent(out) :: d2m

      if     (any(stencil_mask == [0, 0, 0])) then
         d2m = c0
      else
         d2m = c1
      endif

   end subroutine set_d2_mask

   pure subroutine h_edges_nosc(ijdm, ijs, ije, hevc1, hevc2, hevc3, hevc4, &
                                ssc, scc, d2m, hm, hel, her)
   ! ---------------------------------------------------------------------------
   ! Reconstruct thickness edge values using a 4th order scheme and applying
   ! non-oscillatory limiting.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: ijdm, ijs, ije
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(in) :: &
         hevc1, hevc2, hevc3, hevc4, ssc, scc, d2m, hm
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(out) :: hel, her

      real(r8), dimension(1-nbdy:ijdm+nbdy) :: d2h
      real(r8) :: he, sl, sr, sc, d, q, r, a2
      integer :: i

      do i = ijs-1, ije+2
         he = hevc1(i)*hm(i-2) + hevc2(i)*hm(i-1) &
            + hevc3(i)*hm(i  ) + hevc4(i)*hm(i+1)
         hel(i  ) = he
         her(i-1) = he      
      enddo

      do i = ijs-1, ije+1
         d2h(i) = d2m(i)*(hel(i) - c2*hm(i) + her(i))
      enddo

      do i = ijs, ije

         if (d2h(i-1)*d2h(i  ) <= c0 .or. &
             d2h(i  )*d2h(i+1) <= c0) then
            sl = ssc(i)*(hm(i  ) - hm(i-1))
            sr = ssc(i)*(hm(i+1) - hm(i  ))
            if (sl*sr > c0) then
               sc = scc(i)*(hm(i+1) - hm(i-1))
               sc = sign(min(abs(sl), abs(sr), abs(sc)), sc)
               if ((hm(i-1) - hel(i))*(hm(i) - hel(i)) > c0) &
                  hel(i) = hm(i) &
                         - sign(min(c1_2*abs(sc), &
                                    abs(hel(i) - hm(i))), sc)
               if ((hm(i+1) - her(i))*(hm(i) - her(i)) > c0) &
                  her(i) = hm(i) &
                         + sign(min(c1_2*abs(sc), &
                                    abs(her(i) - hm(i))), sc)
               d = her(i) - hel(i)
               q = d*(c2*hm(i) - hel(i) - her(i))
               r = c1_3*d*d
               if     (  q > r) then
                  hel(i) = c3*hm(i) - c2*her(i)
               elseif (- r > q) then
                  her(i) = c3*hm(i) - c2*hel(i)
               endif
            else
               hel(i) = hm(i)
               her(i) = hm(i)
            endif
         endif

         hel(i) = max(hel(i), dpeps)
         her(i) = max(her(i), dpeps)
         sl = c2*(c3*hm(i) - c2*hel(i) - her(i))
         a2 = c3*(hel(i) - c2*hm(i) + her(i))
         sr = sl + c2*a2
         if (sl < c0 .and. sr > c0) then
            if (a2*hel(i) - c1_4*sl*sl < a2*dpeps) then
               q = c3*hm(i)/(c3*sl*sr + c4*a2*a2)
               hel(i) = sl*sl*q
               her(i) = sr*sr*q
            endif
         endif

      enddo

   end subroutine h_edges_nosc

   pure subroutine h_edges_mono(ijdm, ijs, ije, hevc1, hevc2, hevc3, hevc4, &
                                ssc, scc, hm, hel, her)
   ! ---------------------------------------------------------------------------
   ! Reconstruct thickness edge values using a 4th order scheme and applying
   ! monotonic limiting.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: ijdm, ijs, ije
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(in) :: &
         hevc1, hevc2, hevc3, hevc4, ssc, scc, hm
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(out) :: hel, her

      real(r8) :: he, sl, sr, sc, d, q, r
      integer :: i

      do i = ijs, ije+1
         he = hevc1(i)*hm(i-2) + hevc2(i)*hm(i-1) &
            + hevc3(i)*hm(i  ) + hevc4(i)*hm(i+1)
         hel(i  ) = he
         her(i-1) = he      
      enddo

      do i = ijs, ije

         sl = ssc(i)*(hm(i  ) - hm(i-1))
         sr = ssc(i)*(hm(i+1) - hm(i  ))
         if (sl*sr > c0) then
            sc = scc(i)*(hm(i+1) - hm(i-1))
            sc = sign(min(abs(sl), abs(sr), abs(sc)), sc)
            if ((hm(i-1) - hel(i))*(hm(i) - hel(i)) > c0) &
               hel(i) = hm(i) &
                      - sign(min(c1_2*abs(sc), &
                                 abs(hel(i) - hm(i))), sc)
            if ((hm(i+1) - her(i))*(hm(i) - her(i)) > c0) &
               her(i) = hm(i) &
                      + sign(min(c1_2*abs(sc), &
                                 abs(her(i) - hm(i))), sc)
            d = her(i) - hel(i)
            q = d*(c2*hm(i) - hel(i) - her(i))
            r = c1_3*d*d
            if     (  q > r) then
               hel(i) = c3*hm(i) - c2*her(i)
            elseif (- r > q) then
               her(i) = c3*hm(i) - c2*hel(i)
            endif
         else
            hel(i) = hm(i)
            her(i) = hm(i)
         endif

      enddo

   end subroutine h_edges_mono

   pure subroutine parabola_coeffs_fc_nosc(ijdm, ijs, ije, &
                                           stencil, tmc0, tmcl, tmcr, &
                                           ssc, scc, d2m, &
                                           hm, tm, hel, her, &
                                           hpc0, hpc1, hpc2, tpc0, tpc1, tpc2)
   ! ---------------------------------------------------------------------------
   ! Compute the coefficients of the piecewise parabolas. Tracer edge values are
   ! estimated with a 4th order scheme that is compatible with the thickness
   ! reconstruction. Non-oscillatory limiting is applied.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: ijdm, ijs, ije
      integer, dimension(1-nbdy:ijdm+nbdy), intent(in) :: stencil
      real(r8), dimension(12,1-nbdy:ijdm+nbdy), intent(in) :: tmc0, tmcl, tmcr
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(in) :: &
         ssc, scc, d2m, hm, hel, her
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(in) :: tm
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(out) :: hpc0, hpc1, hpc2
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(out) :: &
         tpc0, tpc1, tpc2

      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy) :: d2t, tel, ter
      real(r8), dimension(1-nbdy:ijdm+nbdy) :: &
         hf1m, hf1l, hf1r, hf2m, hf2l, hf2r
      real(r8) :: h1i, h2i, h3i, h4i, a12, a22, a32, a42, a13, a23, a33, a43, &
                  a14, a24, a34, a44, q, tevc1, tevc2, tevc3, tevc4, &
                  te, sl, sr, sc, a2
      integer :: i, nt

      ! Compute compatible tracer edge value coefficients.

      do i = ijs-1, ije+2

         select case (stencil(i))
            case (stencil_1111)

               ! Define elements and complete the LU decomposition of the matrix
               ! for the linear system of equations that relates edge values to
               ! cell means for the 1111 stencil.
               h1i = c1/hm(i-2)
               h2i = c1/hm(i-1)
               h3i = c1/hm(i  )
               h4i = c1/hm(i+1)
               a12 =  tmc0( 1,i) &
                   + (tmcl( 1,i)*hel(i-2) + tmcr( 1,i)*her(i-2))*h1i
               a13 =  tmc0( 2,i) &
                   + (tmcl( 2,i)*hel(i-2) + tmcr( 2,i)*her(i-2))*h1i
               a14 =  tmc0( 3,i) &
                   + (tmcl( 3,i)*hel(i-2) + tmcr( 3,i)*her(i-2))*h1i
               a22 =  tmc0( 4,i) &
                   + (tmcl( 4,i)*hel(i-1) + tmcr( 4,i)*her(i-1))*h2i - a12
               a23 =  tmc0( 5,i) &
                   + (tmcl( 5,i)*hel(i-1) + tmcr( 5,i)*her(i-1))*h2i - a13
               a24 =  tmc0( 6,i) &
                   + (tmcl( 6,i)*hel(i-1) + tmcr( 6,i)*her(i-1))*h2i - a14
               a32 =  tmc0( 7,i) &
                   + (tmcl( 7,i)*hel(i  ) + tmcr( 7,i)*her(i  ))*h3i - a12
               a33 =  tmc0( 8,i) &
                   + (tmcl( 8,i)*hel(i  ) + tmcr( 8,i)*her(i  ))*h3i - a13
               a34 =  tmc0( 9,i) &
                   + (tmcl( 9,i)*hel(i  ) + tmcr( 9,i)*her(i  ))*h3i - a14
               a42 =  tmc0(10,i) &
                   + (tmcl(10,i)*hel(i+1) + tmcr(10,i)*her(i+1))*h4i - a12
               a43 =  tmc0(11,i) &
                   + (tmcl(11,i)*hel(i+1) + tmcr(11,i)*her(i+1))*h4i - a13
               a44 =  tmc0(12,i) &
                   + (tmcl(12,i)*hel(i+1) + tmcr(12,i)*her(i+1))*h4i - a14
               q = c1/a22
               a23 = a23*q
               a33 = a33 - a23*a32
               a43 = a43 - a23*a42
               a24 = a24*q
               a34 = a34 - a24*a32
               a44 = a44 - a24*a42
               a34 = a34/a33
               a44 = a44 - a34*a43

               ! Forward/backward substitution to obtain coefficients for tracer
               ! edge value reconstruction.
               tevc2 = - a12
               tevc3 = - a13 - a23*tevc2
               tevc4 = - a14 - a24*tevc2 - a34*tevc3
               tevc4 = tevc4/a44
               tevc3 = (tevc3 - a43*tevc4)/a33
               tevc2 = (tevc2 - a32*tevc3 - a42*tevc4)/a22
               tevc1 = c1 - tevc2 - tevc3 - tevc4

            case (stencil_0000)

               ! Coefficients for the 0000 stencil.
               tevc1 = c0
               tevc2 = c0
               tevc3 = c0
               tevc4 = c0

            case (stencil_1110)

               ! Define elements and complete the LU decomposition of the matrix
               ! for the linear system of equations that relates edge values to
               ! cell means for the 1110 stencil.
               h1i = c1/hm(i-2)
               h2i = c1/hm(i-1)
               h3i = c1/hm(i  )
               a12 =  tmc0( 1,i) &
                   + (tmcl( 1,i)*hel(i-2) + tmcr( 1,i)*her(i-2))*h1i
               a13 =  tmc0( 2,i) &
                   + (tmcl( 2,i)*hel(i-2) + tmcr( 2,i)*her(i-2))*h1i
               a22 =  tmc0( 4,i) &
                   + (tmcl( 4,i)*hel(i-1) + tmcr( 4,i)*her(i-1))*h2i - a12
               a23 =  tmc0( 5,i) &
                   + (tmcl( 5,i)*hel(i-1) + tmcr( 5,i)*her(i-1))*h2i - a13
               a32 =  tmc0( 7,i) &
                   + (tmcl( 7,i)*hel(i  ) + tmcr( 7,i)*her(i  ))*h3i - a12
               a33 =  tmc0( 8,i) &
                   + (tmcl( 8,i)*hel(i  ) + tmcr( 8,i)*her(i  ))*h3i - a13
               a23 = a23/a22
               a33 = a33 - a23*a32

               ! Forward/backward substitution to obtain coefficients for tracer
               ! edge value reconstruction.
               tevc2 = - a12
               tevc3 = - a13 - a23*tevc2
               tevc3 = tevc3/a33
               tevc2 = (tevc2 - a32*tevc3)/a22
               tevc1 = c1 - tevc2 - tevc3
               tevc4 = c0

            case (stencil_0111)

               ! Define elements and complete the LU decomposition of the matrix
               ! for the linear system of equations that relates edge values to
               ! cell means for the 0111 stencil.
               h2i = c1/hm(i-1)
               h3i = c1/hm(i  )
               h4i = c1/hm(i+1)
               a22 =  tmc0( 4,i) &
                   + (tmcl( 4,i)*hel(i-1) + tmcr( 4,i)*her(i-1))*h2i
               a23 =  tmc0( 5,i) &
                   + (tmcl( 5,i)*hel(i-1) + tmcr( 5,i)*her(i-1))*h2i
               a32 =  tmc0( 7,i) &
                   + (tmcl( 7,i)*hel(i  ) + tmcr( 7,i)*her(i  ))*h3i - a22
               a33 =  tmc0( 8,i) &
                   + (tmcl( 8,i)*hel(i  ) + tmcr( 8,i)*her(i  ))*h3i - a23
               a42 =  tmc0(10,i) &
                   + (tmcl(10,i)*hel(i+1) + tmcr(10,i)*her(i+1))*h4i - a22
               a43 =  tmc0(11,i) &
                   + (tmcl(11,i)*hel(i+1) + tmcr(11,i)*her(i+1))*h4i - a23
               a33 = a33/a32
               a43 = a43 - a33*a42

               ! Forward/backward substitution to obtain coefficients for tracer
               ! edge value reconstruction.
               tevc3 = - a22
               tevc4 = - a23 - a33*tevc3
               tevc4 = tevc4/a43
               tevc3 = (tevc3 - a42*tevc4)/a32
               tevc2 = c1 - tevc3 - tevc4
               tevc1 = c0

            case (stencil_1100)

               ! Define elements and complete the LU decomposition of the matrix
               ! for the linear system of equations that relates edge values to
               ! cell means for the 1100 stencil.
               h1i = c1/hm(i-2)
               h2i = c1/hm(i-1)
               a12 =  tmc0( 1,i) &
                   + (tmcl( 1,i)*hel(i-2) + tmcr( 1,i)*her(i-2))*h1i
               a22 =  tmc0( 4,i) &
                   + (tmcl( 4,i)*hel(i-1) + tmcr( 4,i)*her(i-1))*h2i - a12

               ! Forward/backward substitution to obtain coefficients for tracer
               ! edge value reconstruction.
               tevc2 = - a12/a22
               tevc1 = c1 - tevc2
               tevc3 = c0
               tevc4 = c0

            case (stencil_0110)

               ! Define elements and complete the LU decomposition of the matrix
               ! for the linear system of equations that relates edge values to
               ! cell means for the 0110 stencil.
               h2i = c1/hm(i-1)
               h3i = c1/hm(i  )
               a22 =  tmc0( 4,i) &
                   + (tmcl( 4,i)*hel(i-1) + tmcr( 4,i)*her(i-1))*h2i
               a32 =  tmc0( 7,i) &
                   + (tmcl( 7,i)*hel(i  ) + tmcr( 7,i)*her(i  ))*h3i - a22

               ! Forward/backward substitution to obtain coefficients for tracer
               ! edge value reconstruction.
               tevc3 = - a22/a32
               tevc2 = c1 - tevc3
               tevc1 = c0
               tevc4 = c0

            case (stencil_0011)

               ! Define elements and complete the LU decomposition of the matrix
               ! for the linear system of equations that relates edge values to
               ! cell means for the 0011 stencil.
               h3i = c1/hm(i  )
               h4i = c1/hm(i+1)
               a32 =  tmc0( 7,i) &
                   + (tmcl( 7,i)*hel(i  ) + tmcr( 7,i)*her(i  ))*h3i
               a42 =  tmc0(10,i) &
                   + (tmcl(10,i)*hel(i+1) + tmcr(10,i)*her(i+1))*h4i - a32

               ! Forward/backward substitution to obtain coefficients for tracer
               ! edge value reconstruction.
               tevc4 = - a32/a42
               tevc3 = c1 - tevc4
               tevc1 = c0
               tevc2 = c0

            case (stencil_0100)

               ! Coefficients for the 0100 stencil.
               tevc1 = c0
               tevc2 = c1
               tevc3 = c0
               tevc4 = c0

            case (stencil_0010)

               ! Coefficients for the 0010 stencil.
               tevc1 = c0
               tevc2 = c0
               tevc3 = c1
               tevc4 = c0

         end select

         do nt = 1, ntr_loc
            te = tevc1*tm(nt,i-2) + tevc2*tm(nt,i-1) &
               + tevc3*tm(nt,i  ) + tevc4*tm(nt,i+1)
            tel(nt,i  ) = te
            ter(nt,i-1) = te
         enddo

      enddo

      ! Compute thickness dependent factors for the computation of tracer
      ! parabola coefficients.
      do i = ijs-1, ije+1
         q = c1/(c12*hm(i) - hel(i) - her(i))
         hf1m(i) = c60*hm(i)*q
         hf1l(i) = - (c42*hm(i) + c4*hel(i) - c6*her(i))*q
         hf1r(i) = - (c18*hm(i) - c4*hel(i) + c6*her(i))*q
         hf2m(i) = - hf1m(i)
         hf2l(i) =   c5*(c6*hm(i) + hel(i) - her(i))*q
         hf2r(i) =   c5*(c6*hm(i) - hel(i) + her(i))*q
         do nt = 1, ntr_loc
            d2t(nt,i) = d2m(i)*( hf2m(i)*tm(nt,i) &
                               + hf2l(i)*tel(nt,i) &
                               + hf2r(i)*ter(nt,i))
         enddo
      enddo

      do i = ijs, ije

         do nt = 1, ntr_loc
            if (d2t(nt,i-1)*d2t(nt,i  ) <= c0 .or. &
                d2t(nt,i  )*d2t(nt,i+1) <= c0) then
               sl = ssc(i)*(tm(nt,i  ) - tm(nt,i-1))
               sr = ssc(i)*(tm(nt,i+1) - tm(nt,i  ))
               if (sl*sr > c0) then
                  sc = scc(i)*(tm(nt,i+1) - tm(nt,i-1))
                  sc = sign(min(abs(sl), abs(sr), abs(sc)), sc)
                  if ((tm(nt,i-1) - tel(nt,i))*(tm(nt,i) - tel(nt,i)) > c0) &
                     tel(nt,i) = tm(nt,i) &
                               - sign(min(c1_2*abs(sc), &
                                          abs(tel(nt,i) - tm(nt,i))), sc)
                  if ((tm(nt,i+1) - ter(nt,i))*(tm(nt,i) - ter(nt,i)) > c0) &
                     ter(nt,i) = tm(nt,i) &
                               + sign(min(c1_2*abs(sc), &
                                          abs(ter(nt,i) - tm(nt,i))), sc)
                  sl = hf1m(i)*tm(nt,i) + hf1l(i)*tel(nt,i) + hf1r(i)*ter(nt,i)
                  a2 = hf2m(i)*tm(nt,i) + hf2l(i)*tel(nt,i) + hf2r(i)*ter(nt,i)
                  sr = sl + c2*a2
                  if (sl*sr < c0) then
                     if ((ter(nt,i) - tel(nt,i))*a2 < c0) then
                        tel(nt,i) = - ( (hf1m(i) + c2*hf2m(i))*tm(nt,i) &
                                      + (hf1r(i) + c2*hf2r(i))*ter(nt,i)) &
                                      /(hf1l(i) + c2*hf2l(i))
                     else
                        ter(nt,i) = - (hf1m(i)*tm(nt,i) + hf1l(i)*tel(nt,i)) &
                                      /hf1r(i)
                     endif
                  endif
               else
                  tel(nt,i) = tm(nt,i)
                  ter(nt,i) = tm(nt,i)
               endif
            endif
         enddo

         do nt = 2, ntr_loc
            tel(nt,i) = max(tel(nt,i), c0)
            ter(nt,i) = max(ter(nt,i), c0)
            sl = hf1m(i)*tm(nt,i) + hf1l(i)*tel(nt,i) + hf1r(i)*ter(nt,i)
            a2 = hf2m(i)*tm(nt,i) + hf2l(i)*tel(nt,i) + hf2r(i)*ter(nt,i)
            sr = sl + c2*a2
            if (sl < c0 .and. sr > c0) then
               if (a2*tel(nt,i) - c1_4*sl*sl < c0) then
                  q = c3*tm(nt,i)/(c3*sl*sr + c4*a2*a2)
                  tel(nt,i) = sl*sl*q
                  ter(nt,i) = sr*sr*q
               endif
            endif
         enddo

         hpc0(i) = hel(i)
         hpc1(i) = c6*hm(i) - c4*hel(i) - c2*her(i)
         hpc2(i) = c3*(hel(i) - c2*hm(i) + her(i))
         do nt = 1, ntr_loc
            tpc0(nt,i) = tel(nt,i)
            tpc1(nt,i) = hf1m(i)*tm(nt,i) &
                       + hf1l(i)*tel(nt,i) &
                       + hf1r(i)*ter(nt,i)
            tpc2(nt,i) = hf2m(i)*tm(nt,i) &
                       + hf2l(i)*tel(nt,i) &
                       + hf2r(i)*ter(nt,i)
         enddo

      enddo

   end subroutine parabola_coeffs_fc_nosc

   pure subroutine parabola_coeffs_fc_mono(ijdm, ijs, ije, &
                                           stencil, tmc0, tmcl, tmcr, &
                                           ssc, scc, &
                                           hm, tm, hel, her, &
                                           hpc0, hpc1, hpc2, tpc0, tpc1, tpc2)
   ! ---------------------------------------------------------------------------
   ! Compute the coefficients of the piecewise parabolas. Tracer edge values are
   ! estimated with a 4th order scheme that is compatible with the thickness
   ! reconstruction. Monotonic limiting is applied.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: ijdm, ijs, ije
      integer, dimension(1-nbdy:ijdm+nbdy), intent(in) :: stencil
      real(r8), dimension(12,1-nbdy:ijdm+nbdy), intent(in) :: tmc0, tmcl, tmcr
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(in) :: &
         ssc, scc, hm, hel, her
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(in) :: tm
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(out) :: hpc0, hpc1, hpc2
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(out) :: &
         tpc0, tpc1, tpc2

      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy) :: tel, ter
      real(r8) :: h1i, h2i, h3i, h4i, a12, a22, a32, a42, a13, a23, a33, a43, &
                  a14, a24, a34, a44, q, tevc1, tevc2, tevc3, tevc4, &
                  te, hf1m, hf1l, hf1r, hf2m, hf2l, hf2r, sl, sr, sc, a2
      integer :: i, nt

      ! Compute compatible tracer edge value coefficients.

      do i = ijs, ije+1

         select case (stencil(i))
            case (stencil_1111)

               ! Define elements and complete the LU decomposition of the matrix
               ! for the linear system of equations that relates edge values to
               ! cell means for the 1111 stencil.
               h1i = c1/hm(i-2)
               h2i = c1/hm(i-1)
               h3i = c1/hm(i  )
               h4i = c1/hm(i+1)
               a12 =  tmc0( 1,i) &
                   + (tmcl( 1,i)*hel(i-2) + tmcr( 1,i)*her(i-2))*h1i
               a13 =  tmc0( 2,i) &
                   + (tmcl( 2,i)*hel(i-2) + tmcr( 2,i)*her(i-2))*h1i
               a14 =  tmc0( 3,i) &
                   + (tmcl( 3,i)*hel(i-2) + tmcr( 3,i)*her(i-2))*h1i
               a22 =  tmc0( 4,i) &
                   + (tmcl( 4,i)*hel(i-1) + tmcr( 4,i)*her(i-1))*h2i - a12
               a23 =  tmc0( 5,i) &
                   + (tmcl( 5,i)*hel(i-1) + tmcr( 5,i)*her(i-1))*h2i - a13
               a24 =  tmc0( 6,i) &
                   + (tmcl( 6,i)*hel(i-1) + tmcr( 6,i)*her(i-1))*h2i - a14
               a32 =  tmc0( 7,i) &
                   + (tmcl( 7,i)*hel(i  ) + tmcr( 7,i)*her(i  ))*h3i - a12
               a33 =  tmc0( 8,i) &
                   + (tmcl( 8,i)*hel(i  ) + tmcr( 8,i)*her(i  ))*h3i - a13
               a34 =  tmc0( 9,i) &
                   + (tmcl( 9,i)*hel(i  ) + tmcr( 9,i)*her(i  ))*h3i - a14
               a42 =  tmc0(10,i) &
                   + (tmcl(10,i)*hel(i+1) + tmcr(10,i)*her(i+1))*h4i - a12
               a43 =  tmc0(11,i) &
                   + (tmcl(11,i)*hel(i+1) + tmcr(11,i)*her(i+1))*h4i - a13
               a44 =  tmc0(12,i) &
                   + (tmcl(12,i)*hel(i+1) + tmcr(12,i)*her(i+1))*h4i - a14
               q = c1/a22
               a23 = a23*q
               a33 = a33 - a23*a32
               a43 = a43 - a23*a42
               a24 = a24*q
               a34 = a34 - a24*a32
               a44 = a44 - a24*a42
               a34 = a34/a33
               a44 = a44 - a34*a43

               ! Forward/backward substitution to obtain coefficients for tracer
               ! edge value reconstruction.
               tevc2 = - a12
               tevc3 = - a13 - a23*tevc2
               tevc4 = - a14 - a24*tevc2 - a34*tevc3
               tevc4 = tevc4/a44
               tevc3 = (tevc3 - a43*tevc4)/a33
               tevc2 = (tevc2 - a32*tevc3 - a42*tevc4)/a22
               tevc1 = c1 - tevc2 - tevc3 - tevc4

            case (stencil_0000)

               ! Coefficients for the 0000 stencil.
               tevc1 = c0
               tevc2 = c0
               tevc3 = c0
               tevc4 = c0

            case (stencil_1110)

               ! Define elements and complete the LU decomposition of the matrix
               ! for the linear system of equations that relates edge values to
               ! cell means for the 1110 stencil.
               h1i = c1/hm(i-2)
               h2i = c1/hm(i-1)
               h3i = c1/hm(i  )
               a12 =  tmc0( 1,i) &
                   + (tmcl( 1,i)*hel(i-2) + tmcr( 1,i)*her(i-2))*h1i
               a13 =  tmc0( 2,i) &
                   + (tmcl( 2,i)*hel(i-2) + tmcr( 2,i)*her(i-2))*h1i
               a22 =  tmc0( 4,i) &
                   + (tmcl( 4,i)*hel(i-1) + tmcr( 4,i)*her(i-1))*h2i - a12
               a23 =  tmc0( 5,i) &
                   + (tmcl( 5,i)*hel(i-1) + tmcr( 5,i)*her(i-1))*h2i - a13
               a32 =  tmc0( 7,i) &
                   + (tmcl( 7,i)*hel(i  ) + tmcr( 7,i)*her(i  ))*h3i - a12
               a33 =  tmc0( 8,i) &
                   + (tmcl( 8,i)*hel(i  ) + tmcr( 8,i)*her(i  ))*h3i - a13
               a23 = a23/a22
               a33 = a33 - a23*a32

               ! Forward/backward substitution to obtain coefficients for tracer
               ! edge value reconstruction.
               tevc2 = - a12
               tevc3 = - a13 - a23*tevc2
               tevc3 = tevc3/a33
               tevc2 = (tevc2 - a32*tevc3)/a22
               tevc1 = c1 - tevc2 - tevc3
               tevc4 = c0

            case (stencil_0111)

               ! Define elements and complete the LU decomposition of the matrix
               ! for the linear system of equations that relates edge values to
               ! cell means for the 0111 stencil.
               h2i = c1/hm(i-1)
               h3i = c1/hm(i  )
               h4i = c1/hm(i+1)
               a22 =  tmc0( 4,i) &
                   + (tmcl( 4,i)*hel(i-1) + tmcr( 4,i)*her(i-1))*h2i
               a23 =  tmc0( 5,i) &
                   + (tmcl( 5,i)*hel(i-1) + tmcr( 5,i)*her(i-1))*h2i
               a32 =  tmc0( 7,i) &
                   + (tmcl( 7,i)*hel(i  ) + tmcr( 7,i)*her(i  ))*h3i - a22
               a33 =  tmc0( 8,i) &
                   + (tmcl( 8,i)*hel(i  ) + tmcr( 8,i)*her(i  ))*h3i - a23
               a42 =  tmc0(10,i) &
                   + (tmcl(10,i)*hel(i+1) + tmcr(10,i)*her(i+1))*h4i - a22
               a43 =  tmc0(11,i) &
                   + (tmcl(11,i)*hel(i+1) + tmcr(11,i)*her(i+1))*h4i - a23
               a33 = a33/a32
               a43 = a43 - a33*a42

               ! Forward/backward substitution to obtain coefficients for tracer
               ! edge value reconstruction.
               tevc3 = - a22
               tevc4 = - a23 - a33*tevc3
               tevc4 = tevc4/a43
               tevc3 = (tevc3 - a42*tevc4)/a32
               tevc2 = c1 - tevc3 - tevc4
               tevc1 = c0

            case (stencil_1100)

               ! Define elements and complete the LU decomposition of the matrix
               ! for the linear system of equations that relates edge values to
               ! cell means for the 1100 stencil.
               h1i = c1/hm(i-2)
               h2i = c1/hm(i-1)
               a12 =  tmc0( 1,i) &
                   + (tmcl( 1,i)*hel(i-2) + tmcr( 1,i)*her(i-2))*h1i
               a22 =  tmc0( 4,i) &
                   + (tmcl( 4,i)*hel(i-1) + tmcr( 4,i)*her(i-1))*h2i - a12

               ! Forward/backward substitution to obtain coefficients for tracer
               ! edge value reconstruction.
               tevc2 = - a12/a22
               tevc1 = c1 - tevc2
               tevc3 = c0
               tevc4 = c0

            case (stencil_0110)

               ! Define elements and complete the LU decomposition of the matrix
               ! for the linear system of equations that relates edge values to
               ! cell means for the 0110 stencil.
               h2i = c1/hm(i-1)
               h3i = c1/hm(i  )
               a22 =  tmc0( 4,i) &
                   + (tmcl( 4,i)*hel(i-1) + tmcr( 4,i)*her(i-1))*h2i
               a32 =  tmc0( 7,i) &
                   + (tmcl( 7,i)*hel(i  ) + tmcr( 7,i)*her(i  ))*h3i - a22

               ! Forward/backward substitution to obtain coefficients for tracer
               ! edge value reconstruction.
               tevc3 = - a22/a32
               tevc2 = c1 - tevc3
               tevc1 = c0
               tevc4 = c0

            case (stencil_0011)

               ! Define elements and complete the LU decomposition of the matrix
               ! for the linear system of equations that relates edge values to
               ! cell means for the 0011 stencil.
               h3i = c1/hm(i  )
               h4i = c1/hm(i+1)
               a32 =  tmc0( 7,i) &
                   + (tmcl( 7,i)*hel(i  ) + tmcr( 7,i)*her(i  ))*h3i
               a42 =  tmc0(10,i) &
                   + (tmcl(10,i)*hel(i+1) + tmcr(10,i)*her(i+1))*h4i - a32

               ! Forward/backward substitution to obtain coefficients for tracer
               ! edge value reconstruction.
               tevc4 = - a32/a42
               tevc3 = c1 - tevc4
               tevc1 = c0
               tevc2 = c0

            case (stencil_0100)

               ! Coefficients for the 0100 stencil.
               tevc1 = c0
               tevc2 = c1
               tevc3 = c0
               tevc4 = c0

            case (stencil_0010)

               ! Coefficients for the 0010 stencil.
               tevc1 = c0
               tevc2 = c0
               tevc3 = c1
               tevc4 = c0

         end select

         do nt = 1, ntr_loc
            te = tevc1*tm(nt,i-2) + tevc2*tm(nt,i-1) &
               + tevc3*tm(nt,i  ) + tevc4*tm(nt,i+1)
            tel(nt,i  ) = te
            ter(nt,i-1) = te
         enddo

      enddo

      do i = ijs, ije

         ! Compute thickness dependent factors for the computation of tracer
         ! parabola coefficients.
         q = c1/(c12*hm(i) - hel(i) - her(i))
         hf1m = c60*hm(i)*q
         hf1l = - (c42*hm(i) + c4*hel(i) - c6*her(i))*q
         hf1r = - (c18*hm(i) - c4*hel(i) + c6*her(i))*q
         hf2m = - hf1m
         hf2l =   c5*(c6*hm(i) + hel(i) - her(i))*q
         hf2r =   c5*(c6*hm(i) - hel(i) + her(i))*q

         do nt = 1, ntr_loc
            sl = ssc(i)*(tm(nt,i  ) - tm(nt,i-1))
            sr = ssc(i)*(tm(nt,i+1) - tm(nt,i  ))
            if (sl*sr > c0) then
               sc = scc(i)*(tm(nt,i+1) - tm(nt,i-1))
               sc = sign(min(abs(sl), abs(sr), abs(sc)), sc)
               if ((tm(nt,i-1) - tel(nt,i))*(tm(nt,i) - tel(nt,i)) > c0) &
                  tel(nt,i) = tm(nt,i) &
                            - sign(min(c1_2*abs(sc), &
                                       abs(tel(nt,i) - tm(nt,i))), sc)
               if ((tm(nt,i+1) - ter(nt,i))*(tm(nt,i) - ter(nt,i)) > c0) &
                  ter(nt,i) = tm(nt,i) &
                            + sign(min(c1_2*abs(sc), &
                                       abs(ter(nt,i) - tm(nt,i))), sc)
               sl = hf1m*tm(nt,i) + hf1l*tel(nt,i) + hf1r*ter(nt,i)
               a2 = hf2m*tm(nt,i) + hf2l*tel(nt,i) + hf2r*ter(nt,i)
               sr = sl + c2*a2
               if (sl*sr < c0) then
                  if ((ter(nt,i) - tel(nt,i))*a2 < c0) then
                     tel(nt,i) = - ( (hf1m + c2*hf2m)*tm(nt,i) &
                                   + (hf1r + c2*hf2r)*ter(nt,i)) &
                                   /(hf1l + c2*hf2l)
                  else
                     ter(nt,i) = - (hf1m*tm(nt,i) + hf1l*tel(nt,i))/hf1r
                  endif
               endif
            else
               tel(nt,i) = tm(nt,i)
               ter(nt,i) = tm(nt,i)
            endif
         enddo

         hpc0(i) = hel(i)
         hpc1(i) = c6*hm(i) - c4*hel(i) - c2*her(i)
         hpc2(i) = c3*(hel(i) - c2*hm(i) + her(i))
         do nt = 1, ntr_loc
            tpc0(nt,i) = tel(nt,i)
            tpc1(nt,i) = hf1m*tm(nt,i) + hf1l*tel(nt,i) + hf1r*ter(nt,i)
            tpc2(nt,i) = hf2m*tm(nt,i) + hf2l*tel(nt,i) + hf2r*ter(nt,i)
         enddo

      enddo

   end subroutine parabola_coeffs_fc_mono

   pure subroutine parabola_coeffs_pc_nosc(ijdm, ijs, ije, &
                                           hevc1, hevc2, hevc3, hevc4, &
                                           ssc, scc, d2m, &
                                           hm, tm, &
                                           hpc0, hpc1, hpc2, tpc0, tpc1, tpc2)
   ! ---------------------------------------------------------------------------
   ! Compute the coefficients of the piecewise parabolas. Edge values are
   ! estimated with a 4th order scheme. The tracer reconstruction is not
   ! compatible with the thickness reconstruction. Non-oscillatory limiting is
   ! applied.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: ijdm, ijs, ije
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(in) :: &
         hevc1, hevc2, hevc3, hevc4, ssc, scc, d2m, hm
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(in) :: tm
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(out) :: hpc0, hpc1, hpc2
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(out) :: &
         tpc0, tpc1, tpc2

      real(r8), dimension(1-nbdy:ijdm+nbdy) :: hel, her, d2h
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy) :: tel, ter, d2t
      real(r8) :: he, te, sl, sr, sc, d, q, r, a2
      integer :: i, nt

      do i = ijs-1, ije+2
         he = hevc1(i)*hm(i-2) + hevc2(i)*hm(i-1) &
            + hevc3(i)*hm(i  ) + hevc4(i)*hm(i+1)
         hel(i  ) = he
         her(i-1) = he      
         do nt = 1, ntr_loc
            te = hevc1(i)*tm(nt,i-2) + hevc2(i)*tm(nt,i-1) &
               + hevc3(i)*tm(nt,i  ) + hevc4(i)*tm(nt,i+1)
            tel(nt,i  ) = te
            ter(nt,i-1) = te
         enddo
      enddo

      do i = ijs-1, ije+1
         d2h(i) = d2m(i)*(hel(i) - c2*hm(i) + her(i))
         do nt = 1, ntr_loc
            d2t(nt,i) = d2m(i)*(tel(nt,i) - c2*tm(nt,i) + ter(nt,i))
         enddo
      enddo

      do i = ijs, ije

         if (d2h(i-1)*d2h(i  ) <= c0 .or. &
             d2h(i  )*d2h(i+1) <= c0) then
            sl = ssc(i)*(hm(i  ) - hm(i-1))
            sr = ssc(i)*(hm(i+1) - hm(i  ))
            if (sl*sr > c0) then
               sc = scc(i)*(hm(i+1) - hm(i-1))
               sc = sign(min(abs(sl), abs(sr), abs(sc)), sc)
               if ((hm(i-1) - hel(i))*(hm(i) - hel(i)) > c0) &
                  hel(i) = hm(i) &
                         - sign(min(c1_2*abs(sc), &
                                    abs(hel(i) - hm(i))), sc)
               if ((hm(i+1) - her(i))*(hm(i) - her(i)) > c0) &
                  her(i) = hm(i) &
                         + sign(min(c1_2*abs(sc), &
                                    abs(her(i) - hm(i))), sc)
               d = her(i) - hel(i)
               q = d*(c2*hm(i) - hel(i) - her(i))
               r = c1_3*d*d
               if     (  q > r) then
                  hel(i) = c3*hm(i) - c2*her(i)
               elseif (- r > q) then
                  her(i) = c3*hm(i) - c2*hel(i)
               endif
            else
               hel(i) = hm(i)
               her(i) = hm(i)
            endif
         endif

         hel(i) = max(hel(i), dpeps)
         her(i) = max(her(i), dpeps)
         sl = c2*(c3*hm(i) - c2*hel(i) - her(i))
         a2 = c3*(hel(i) - c2*hm(i) + her(i))
         sr = sl + c2*a2
         if (sl < c0 .and. sr > c0) then
            if (a2*hel(i) - c1_4*sl*sl < a2*dpeps) then
               q = c3*hm(i)/(c3*sl*sr + c4*a2*a2)
               hel(i) = sl*sl*q
               her(i) = sr*sr*q
            endif
         endif

         do nt = 1, ntr_loc
            if (d2t(nt,i-1)*d2t(nt,i  ) <= c0 .or. &
                d2t(nt,i  )*d2t(nt,i+1) <= c0) then
               sl = ssc(i)*(tm(nt,i  ) - tm(nt,i-1))
               sr = ssc(i)*(tm(nt,i+1) - tm(nt,i  ))
               if (sl*sr > c0) then
                  sc = scc(i)*(tm(nt,i+1) - tm(nt,i-1))
                  sc = sign(min(abs(sl), abs(sr), abs(sc)), sc)
                  if ((tm(nt,i-1) - tel(nt,i))*(tm(nt,i) - tel(nt,i)) > c0) &
                     tel(nt,i) = tm(nt,i) &
                               - sign(min(c1_2*abs(sc), &
                                          abs(tel(nt,i) - tm(nt,i))), sc)
                  if ((tm(nt,i+1) - ter(nt,i))*(tm(nt,i) - ter(nt,i)) > c0) &
                     ter(nt,i) = tm(nt,i) &
                               + sign(min(c1_2*abs(sc), &
                                          abs(ter(nt,i) - tm(nt,i))), sc)
                  d = ter(nt,i) - tel(nt,i)
                  q = d*(c2*tm(nt,i) - tel(nt,i) - ter(nt,i))
                  r = c1_3*d*d
                  if     (  q > r) then
                     tel(nt,i) = c3*tm(nt,i) - c2*ter(nt,i)
                  elseif (- r > q) then
                     ter(nt,i) = c3*tm(nt,i) - c2*tel(nt,i)
                  endif
               else
                  tel(nt,i) = tm(nt,i)
                  ter(nt,i) = tm(nt,i)
               endif
            endif
         enddo

         do nt = 2, ntr_loc
            tel(nt,i) = max(tel(nt,i), c0)
            ter(nt,i) = max(ter(nt,i), c0)
            sl = c2*(c3*tm(nt,i) - c2*tel(nt,i) - ter(nt,i))
            a2 = c3*(tel(nt,i) - c2*tm(nt,i) + ter(nt,i))
            sr = sl + c2*a2
            if (sl < c0 .and. sr > c0) then
               if (a2*tel(nt,i) - c1_4*sl*sl < c0) then
                  q = c3*tm(nt,i)/(c3*sl*sr + c4*a2*a2)
                  tel(nt,i) = sl*sl*q
                  ter(nt,i) = sr*sr*q
               endif
            endif
         enddo

         hpc0(i) = hel(i)
         hpc1(i) = c6*hm(i) - c4*hel(i) - c2*her(i)
         hpc2(i) = c3*(hel(i) - c2*hm(i) + her(i))
         do nt = 1, ntr_loc
            tpc0(nt,i) = tel(nt,i)
            tpc1(nt,i) = c6*tm(nt,i) - c4*tel(nt,i) - c2*ter(nt,i)
            tpc2(nt,i) = c3*(tel(nt,i) - c2*tm(nt,i) + ter(nt,i))
         enddo

      enddo

   end subroutine parabola_coeffs_pc_nosc

   pure subroutine parabola_coeffs_pc_mono(ijdm, ijs, ije, &
                                           hevc1, hevc2, hevc3, hevc4, &
                                           ssc, scc, &
                                           hm, tm, &
                                           hpc0, hpc1, hpc2, tpc0, tpc1, tpc2)
   ! ---------------------------------------------------------------------------
   ! Compute the coefficients of the piecewise parabolas. Edge values are
   ! estimated with a 4th order scheme. The tracer reconstruction is not
   ! compatible with the thickness reconstruction. Non-oscillatory limiting is
   ! applied.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: ijdm, ijs, ije
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(in) :: &
         hevc1, hevc2, hevc3, hevc4, ssc, scc, hm
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(in) :: tm
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(out) :: hpc0, hpc1, hpc2
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(out) :: &
         tpc0, tpc1, tpc2

      real(r8), dimension(1-nbdy:ijdm+nbdy) :: hel, her
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy) :: tel, ter
      real(r8) :: he, te, sl, sr, sc, d, q, r
      integer :: i, nt

      do i = ijs, ije+1
         he = hevc1(i)*hm(i-2) + hevc2(i)*hm(i-1) &
            + hevc3(i)*hm(i  ) + hevc4(i)*hm(i+1)
         hel(i  ) = he
         her(i-1) = he      
         do nt = 1, ntr_loc
            te = hevc1(i)*tm(nt,i-2) + hevc2(i)*tm(nt,i-1) &
               + hevc3(i)*tm(nt,i  ) + hevc4(i)*tm(nt,i+1)
            tel(nt,i  ) = te
            ter(nt,i-1) = te
         enddo
      enddo

      do i = ijs, ije

         sl = ssc(i)*(hm(i  ) - hm(i-1))
         sr = ssc(i)*(hm(i+1) - hm(i  ))
         if (sl*sr > c0) then
            sc = scc(i)*(hm(i+1) - hm(i-1))
            sc = sign(min(abs(sl), abs(sr), abs(sc)), sc)
            if ((hm(i-1) - hel(i))*(hm(i) - hel(i)) > c0) &
               hel(i) = hm(i) &
                      - sign(min(c1_2*abs(sc), &
                                 abs(hel(i) - hm(i))), sc)
            if ((hm(i+1) - her(i))*(hm(i) - her(i)) > c0) &
               her(i) = hm(i) &
                      + sign(min(c1_2*abs(sc), &
                                 abs(her(i) - hm(i))), sc)
            d = her(i) - hel(i)
            q = d*(c2*hm(i) - hel(i) - her(i))
            r = c1_3*d*d
            if     (  q > r) then
               hel(i) = c3*hm(i) - c2*her(i)
            elseif (- r > q) then
               her(i) = c3*hm(i) - c2*hel(i)
            endif
         else
            hel(i) = hm(i)
            her(i) = hm(i)
         endif

         do nt = 1, ntr_loc
            sl = ssc(i)*(tm(nt,i  ) - tm(nt,i-1))
            sr = ssc(i)*(tm(nt,i+1) - tm(nt,i  ))
            if (sl*sr > c0) then
               sc = scc(i)*(tm(nt,i+1) - tm(nt,i-1))
               sc = sign(min(abs(sl), abs(sr), abs(sc)), sc)
               if ((tm(nt,i-1) - tel(nt,i))*(tm(nt,i) - tel(nt,i)) > c0) &
                  tel(nt,i) = tm(nt,i) &
                            - sign(min(c1_2*abs(sc), &
                                       abs(tel(nt,i) - tm(nt,i))), sc)
               if ((tm(nt,i+1) - ter(nt,i))*(tm(nt,i) - ter(nt,i)) > c0) &
                  ter(nt,i) = tm(nt,i) &
                            + sign(min(c1_2*abs(sc), &
                                       abs(ter(nt,i) - tm(nt,i))), sc)
               d = ter(nt,i) - tel(nt,i)
               q = d*(c2*tm(nt,i) - tel(nt,i) - ter(nt,i))
               r = c1_3*d*d
               if     (  q > r) then
                  tel(nt,i) = c3*tm(nt,i) - c2*ter(nt,i)
               elseif (- r > q) then
                  ter(nt,i) = c3*tm(nt,i) - c2*tel(nt,i)
               endif
            else
               tel(nt,i) = tm(nt,i)
               ter(nt,i) = tm(nt,i)
            endif
         enddo

         hpc0(i) = hel(i)
         hpc1(i) = c6*hm(i) - c4*hel(i) - c2*her(i)
         hpc2(i) = c3*(hel(i) - c2*hm(i) + her(i))
         do nt = 1, ntr_loc
            tpc0(nt,i) = tel(nt,i)
            tpc1(nt,i) = c6*tm(nt,i) - c4*tel(nt,i) - c2*ter(nt,i)
            tpc2(nt,i) = c3*(tel(nt,i) - c2*tm(nt,i) + ter(nt,i))
         enddo

      enddo

   end subroutine parabola_coeffs_pc_mono

   pure subroutine flux_integration(ijdm, ijs, ije, &
                                    ca, ai, db, du, dl, hpc0, hpc1, hpc2, &
                                    tpc0, tpc1, tpc2, &
                                    hf, htf)
   ! ---------------------------------------------------------------------------
   ! To obtain the fluxes, integrate the upstream piecewise parabolas over the
   ! flux area. For tracer fluxes, integrate the product of thickness and tracer
   ! parabolas. If the lower interface of the upstream cell (dl) is below the
   ! lowest edge interface (db), use a piecewise constant reconstruction of the
   ! upstream thickness (hb) set to the thickness of the upstream cell above db.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: ijdm, ijs, ije
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(in) :: &
         ca, ai, db, du, dl, hpc0, hpc1, hpc2
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(in) :: &
         tpc0, tpc1, tpc2
      real(r8), dimension(1-nbdy:ijdm+nbdy), intent(out) :: hf
      real(r8), dimension(ntr_loc,1-nbdy:ijdm+nbdy), intent(out) :: htf

      real(r8) :: c, hb, p0, p1, p2, q1, q2, q3, q4
      integer :: i, nt

      do i = ijs, ije

         if (ca(i) < c0) then

           c = ca(i)*ai(i)

           if (dl(i) > db(i)) then

              hb = max(c0, db(i) - du(i))

              hf(i) = hb*ca(i)

              p0 =        hb
              p1 = - c1_2*hb*c
              p2 =   c1_3*hb*c*c

           else

              hf(i) = (hpc0(i) - (c1_2*hpc1(i) - c1_3*hpc2(i)*c)*c)*ca(i)

              p0 =         hpc0(i) - (c1_2*hpc1(i) - c1_3*hpc2(i)*c)*c 
              p1 = - (c1_2*hpc0(i) - (c1_3*hpc1(i) - c1_4*hpc2(i)*c)*c)*c
              p2 =   (c1_3*hpc0(i) - (c1_4*hpc1(i) - c1_5*hpc2(i)*c)*c)*c*c

           endif

           do nt = 1, ntr_loc
              htf(nt,i) = ( p0*tpc0(nt,i) &
                          + p1*tpc1(nt,i) &
                          + p2*tpc2(nt,i))*ca(i)
           enddo


         else

           c = ca(i)*ai(i-1)

           q1 = c1 - c1_2*c
           q2 = c1 - (c1 - c1_3*c)*c

           if (dl(i-1) > db(i)) then

              hb = max(c0, db(i) - du(i-1))

              hf(i) = hb*ca(i)

              p0 =    hb
              p1 = q1*hb
              p2 = q2*hb

           else

              hf(i) = (hpc0(i-1) + q1*hpc1(i-1) + q2*hpc2(i-1))*ca(i)

              q3 = c1_4*(c1 + c3*(c1 - c)*q2)
              q4 = c1_5*(c1 + c4*(c1 - c)*q3)
              p0 =    hpc0(i-1) + q1*hpc1(i-1) + q2*hpc2(i-1)
              p1 = q1*hpc0(i-1) + q2*hpc1(i-1) + q3*hpc2(i-1)
              p2 = q2*hpc0(i-1) + q3*hpc1(i-1) + q4*hpc2(i-1)

           endif

           do nt = 1, ntr_loc
              htf(nt,i) = ( p0*tpc0(nt,i-1) &
                          + p1*tpc1(nt,i-1) &
                          + p2*tpc2(nt,i-1))*ca(i)
           enddo

         endif

      enddo

   end subroutine flux_integration

   subroutine cppm_fc_nosc_i(m, n, mm, nn, k1m, k1n, second_pass)
   ! ---------------------------------------------------------------------------
   ! Do CPPM transport in the i-direction, applying non-oscillatory limiting.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n
      logical, intent(in) :: second_pass

      real(r8), dimension(1-nbdy:idm+nbdy) :: &
         db, dl, du, ca, ai, ho, hm, hel, her, hpc0, hpc1, hpc2, hf
      real(r8), dimension(ntr_loc,1-nbdy:idm+nbdy) :: tm, tpc0, tpc1, tpc2, htf

      real(r8) :: he_tmp, hn, hni
      integer :: i, j, k, km, kn, nt

      call xctilr(dp(1-nbdy,1-nbdy,k1n), 1, kdm, 4, 0, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,k1n), 1, kdm, 4, 0, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1n), 1, kdm, 4, 0, halo_ps)
      do nt = 3, ntr_loc
         call xctilr(trc(1-nbdy,1-nbdy,k1n,nt-2), 1, kdm, 4, 0, halo_ps)
      enddo

      do k = 1, kdm
         kn = k + nn
         do j = 1, jdm

            ! Extract 1D arrays from multidimensional arrays.
            do i = -2, idm + 3
               ai(i) = scp2i(i,j)
               hm(i) = max(c0, dp(i,j,kn)) + dpeps
            enddo

            if (second_pass) then
               ! For the second pass of the split scheme, apply a divergence
               ! correction to the thickness used for subsequent flux
               ! estimations in order to preserve an initially uniform thickness
               ! in non-divergent flow.
               do i = -2, idm + 3
                  hm(i) = hm(i)/(c1 - (cav(i,j+1,k) - cav(i,j,k))*ai(i))
               enddo
            endif

            ! Reconstruct thickness edge values and apply limiting.
            call h_edges_nosc(idm, 1, idm, &
                              hevc1i(:,j), hevc2i(:,j), &
                              hevc3i(:,j), hevc4i(:,j), &
                              ssci(:,j), scci(:,j), d2mi(:,j), &
                              hm, hel, her)

            ! Copy thickness edge values to 3D arrays.
            do i = 1, idm
               hel_3d(i,j,k) = hel(i)
               her_3d(i,j,k) = her(i)
            enddo

         enddo
      enddo

      call xctilr(hel_3d, 1, kdm, 4, 0, halo_ps)
      call xctilr(her_3d, 1, kdm, 4, 0, halo_ps)

      ! With arctic patch, swap the order of thickness edge values.
      if (nreg == 2 .and. nproc == jpr) then
         j = jj
         do k = 1, kdm
            do i = - 3, ii + 4
               he_tmp = hel_3d(i,j,k)
               hel_3d(i,j,k) = her_3d(i,j,k)
               her_3d(i,j,k) = he_tmp
            enddo
         enddo
      endif

      do k = 1, kdm
         km = k + mm
         kn = k + nn
         do j = 1, jdm

            ! Extract 1D arrays from multidimensional arrays.
            do i = 1, idm + 1
               ca(i) = cau(i,j,k)
               db(i) = pbu(i,j,n)
            enddo
            do i = 0, idm + 1
               du(i) = p(i,j,k  )
               dl(i) = p(i,j,k+1)
            enddo
            do i = -3, idm + 4
               ai(i) = scp2i(i,j)
               ho(i) = max(c0, dp(i,j,kn)) + dpeps
               hm(i) = ho(i)
               hel(i) = hel_3d(i,j,k)
               her(i) = her_3d(i,j,k)
               tm(1,i) = temp(i,j,kn)
               tm(2,i) = saln(i,j,kn)
               do nt = 3, ntr_loc
                  tm(nt,i) = trc(i,j,kn,nt-2)
               enddo
            enddo

            if (second_pass) then
               ! For the second pass of the split scheme, apply a divergence
               ! correction to the thickness used for subsequent flux
               ! estimations in order to preserve an initially uniform thickness
               ! in non-divergent flow.
               do i = -3, idm + 4
                  hm(i) = hm(i)/(c1 - (cav(i,j+1,k) - cav(i,j,k))*ai(i))
               enddo
            endif

            ! Reconstruct tracer edge values, apply limiting, define parabolas
            ! and obtain fluxes by integrating over the flux area.

            call parabola_coeffs_fc_nosc(idm, 0, idm + 1, &
                                         stencili(:,j), &
                                         tmc0i(:,:,j), &
                                         tmcli(:,:,j), &
                                         tmcri(:,:,j), &
                                         ssci(:,j), scci(:,j), &
                                         d2mi(:,j), &
                                         hm, tm, hel, her, &
                                         hpc0, hpc1, hpc2, tpc0, tpc1, tpc2)

            call flux_integration(idm, 1, idm + 1, &
                                  ca, ai, db, du, dl, hpc0, hpc1, hpc2, &
                                  tpc0, tpc1, tpc2, &
                                  hf, htf)

            ! Update variables with flux divergences.
            do i = 1, idm
               hn = ho(i) - (hf(i+1) - hf(i))*ai(i)
               hni = c1/hn
               temp(i,j,kn) = ( ho(i)*tm(1,i) &
                              - (htf(1,i+1) - htf(1,i))*ai(i))*hni
               saln(i,j,kn) = ( ho(i)*tm(2,i) &
                              - (htf(2,i+1) - htf(2,i))*ai(i))*hni
               do nt = 3, ntr_loc
                  trc(i,j,kn,nt-2) = ( ho(i)*tm(nt,i) &
                                     - (htf(nt,i+1) - htf(nt,i))*ai(i))*hni
               enddo
               dp(i,j,kn) = max(c0, hn - dpeps)
            enddo

            ! Accumulate fluxes.
            do i = 1, idm + 1
               uflx (i,j,km) = uflx (i,j,km) + hf(i)
               utflx(i,j,km) = utflx(i,j,km) + htf(1,i)
               usflx(i,j,km) = usflx(i,j,km) + htf(2,i)
            enddo

         enddo
      enddo

   end subroutine cppm_fc_nosc_i

   subroutine cppm_fc_nosc_j(m, n, mm, nn, k1m, k1n, second_pass)
   ! ---------------------------------------------------------------------------
   ! Do CPPM transport in the j-direction, applying non-oscillatory limiting.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n
      logical, intent(in) :: second_pass

      real(r8), dimension(1-nbdy:jdm+nbdy) :: &
         db, dl, du, ca, ai, ho, hm, hel, her, hpc0, hpc1, hpc2, hf
      real(r8), dimension(ntr_loc,1-nbdy:jdm+nbdy) :: tm, tpc0, tpc1, tpc2, htf

      real(r8) :: he_tmp, hn, hni
      integer :: i, j, k, km, kn, nt

      call xctilr(dp(1-nbdy,1-nbdy,k1n), 1, kdm, 0, 4, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,k1n), 1, kdm, 0, 4, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1n), 1, kdm, 0, 4, halo_ps)
      do nt = 3, ntr_loc
         call xctilr(trc(1-nbdy,1-nbdy,k1n,nt-2), 1, kdm, 0, 4, halo_ps)
      enddo

      do k = 1, kdm
         kn = k + nn
         do i = 1, idm

            ! Extract 1D arrays from multidimensional arrays.
            do j = -2, jdm + 3
               ai(j) = scp2i(i,j)
               hm(j) = max(c0, dp(i,j,kn)) + dpeps
            enddo

            if (second_pass) then
               ! For the second pass of the split scheme, apply a divergence
               ! correction to the thickness used for subsequent flux
               ! estimations in order to preserve an initially uniform thickness
               ! in non-divergent flow.
               do j = -2, jdm + 3
                  hm(j) = hm(j)/(c1 - (cau(i+1,j,k) - cau(i,j,k))*ai(j))
               enddo
            endif

            ! Reconstruct thickness edge values and apply limiting.
            call h_edges_nosc(jdm, 1, jdm, &
                              hevc1j(:,i), hevc2j(:,i), &
                              hevc3j(:,i), hevc4j(:,i), &
                              sscj(:,i), sccj(:,i), d2mj(:,i), &
                              hm, hel, her)

            ! Copy edge values to 3D arrays.
            do j = 1, jdm
               hel_3d(i,j,k) = hel(j)
               her_3d(i,j,k) = her(j)
            enddo

         enddo
      enddo

      call xctilr(hel_3d, 1, kdm, 0, 4, halo_ps)
      call xctilr(her_3d, 1, kdm, 0, 4, halo_ps)

      ! With arctic patch, swap the order of thickness edge values.
      if (nreg == 2 .and. nproc == jpr) then
         do k = 1, kdm
            j = jj
            do i = max(1, itdm/2 - i0 + 1), ii
               he_tmp = hel_3d(i,j,k)
               hel_3d(i,j,k) = her_3d(i,j,k)
               her_3d(i,j,k) = he_tmp
            enddo
            do j = jj + 1, jj + 4
               do i = 1, ii
                  he_tmp = hel_3d(i,j,k)
                  hel_3d(i,j,k) = her_3d(i,j,k)
                  her_3d(i,j,k) = he_tmp
               enddo
            enddo
         enddo
      endif

      do k = 1, kdm
         km = k + mm
         kn = k + nn
         do i = 1, idm

            ! Extract 1D arrays from multidimensional arrays.
            do j = 1, jdm + 1
               ca(j) = cav(i,j,k)
               db(j) = pbv(i,j,n)
            enddo
            do j = 0, jdm + 1
               du(j) = p(i,j,k  )
               dl(j) = p(i,j,k+1)
            enddo
            do j = -3, jdm + 4
               ai(j) = scp2i(i,j)
               ho(j) = max(c0, dp(i,j,kn)) + dpeps
               hm(j) = ho(j)
               hel(j) = hel_3d(i,j,k)
               her(j) = her_3d(i,j,k)
               tm(1,j) = temp(i,j,kn)
               tm(2,j) = saln(i,j,kn)
               do nt = 3, ntr_loc
                  tm(nt,j) = trc(i,j,kn,nt-2)
               enddo
            enddo

            if (second_pass) then
               ! For the second pass of the split scheme, apply a divergence
               ! correction to the thickness used for subsequent flux
               ! estimations in order to preserve an initially uniform thickness
               ! in non-divergent flow.
               do j = -3, jdm + 4
                  hm(j) = hm(j)/(c1 - (cau(i+1,j,k) - cau(i,j,k))*ai(j))
               enddo
            endif

            ! Reconstruct tracer edge values, apply limiting, define parabolas
            ! and obtain fluxes by integrating over the flux area.

            call parabola_coeffs_fc_nosc(jdm, 0, jdm + 1, &
                                         stencilj(:,i), &
                                         tmc0j(:,:,i), &
                                         tmclj(:,:,i), &
                                         tmcrj(:,:,i), &
                                         sscj(:,i), sccj(:,i), &
                                         d2mj(:,i), &
                                         hm, tm, hel, her, &
                                         hpc0, hpc1, hpc2, tpc0, tpc1, tpc2)

            call flux_integration(jdm, 1, jdm + 1, &
                                  ca, ai, db, du, dl, hpc0, hpc1, hpc2, &
                                  tpc0, tpc1, tpc2, &
                                  hf, htf)

            ! Update variables with flux divergences.
            do j = 1, jdm
               hn = ho(j) - (hf(j+1) - hf(j))*ai(j)
               hni = c1/hn
               temp(i,j,kn) = ( ho(j)*tm(1,j) &
                              - (htf(1,j+1) - htf(1,j))*ai(j))*hni
               saln(i,j,kn) = ( ho(j)*tm(2,j) &
                              - (htf(2,j+1) - htf(2,j))*ai(j))*hni
               do nt = 3, ntr_loc
                  trc(i,j,kn,nt-2) = ( ho(j)*tm(nt,j) &
                                     - (htf(nt,j+1) - htf(nt,j))*ai(j))*hni
               enddo
               dp(i,j,kn) = max(c0, hn - dpeps)
            enddo

            ! Accumulate fluxes.
            do j = 1, jdm + 1
               vflx (i,j,km) = vflx (i,j,km) + hf(j)
               vtflx(i,j,km) = vtflx(i,j,km) + htf(1,j)
               vsflx(i,j,km) = vsflx(i,j,km) + htf(2,j)
            enddo

         enddo
      enddo

   end subroutine cppm_fc_nosc_j

   subroutine cppm_fc_mono_i(m, n, mm, nn, k1m, k1n, second_pass)
   ! ---------------------------------------------------------------------------
   ! Do CPPM transport in the i-direction, applying monotonic limiting.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n
      logical, intent(in) :: second_pass

      real(r8), dimension(1-nbdy:idm+nbdy) :: &
         db, dl, du, ca, ai, ho, hm, hel, her, hpc0, hpc1, hpc2, hf
      real(r8), dimension(ntr_loc,1-nbdy:idm+nbdy) :: tm, tpc0, tpc1, tpc2, htf

      real(r8) :: he_tmp, hn, hni
      integer :: i, j, k, km, kn, nt

      call xctilr(dp(1-nbdy,1-nbdy,k1n), 1, kdm, 3, 0, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,k1n), 1, kdm, 3, 0, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1n), 1, kdm, 3, 0, halo_ps)
      do nt = 3, ntr_loc
         call xctilr(trc(1-nbdy,1-nbdy,k1n,nt-2), 1, kdm, 3, 0, halo_ps)
      enddo

      do k = 1, kdm
         kn = k + nn
         do j = 1, jdm

            ! Extract 1D arrays from multidimensional arrays.
            do i = -2, idm + 3
               ai(i) = scp2i(i,j)
               hm(i) = max(c0, dp(i,j,kn)) + dpeps
            enddo

            if (second_pass) then
               ! For the second pass of the split scheme, apply a divergence
               ! correction to the thickness used for subsequent flux
               ! estimations in order to preserve an initially uniform thickness
               ! in non-divergent flow.
               do i = -2, idm + 3
                  hm(i) = hm(i)/(c1 - (cav(i,j+1,k) - cav(i,j,k))*ai(i))
               enddo
            endif

            ! Reconstruct thickness edge values and apply limiting.
            call h_edges_mono(idm, 1, idm, &
                              hevc1i(:,j), hevc2i(:,j), &
                              hevc3i(:,j), hevc4i(:,j), &
                              ssci(:,j), scci(:,j), &
                              hm, hel, her)

            ! Copy thickness edge values to 3D arrays.
            do i = 1, idm
               hel_3d(i,j,k) = hel(i)
               her_3d(i,j,k) = her(i)
            enddo

         enddo
      enddo

      call xctilr(hel_3d, 1, kdm, 3, 0, halo_ps)
      call xctilr(her_3d, 1, kdm, 3, 0, halo_ps)

      ! With arctic patch, swap the order of thickness edge values.
      if (nreg == 2 .and. nproc == jpr) then
         j = jj
         do k = 1, kdm
            do i = - 2, ii + 3
               he_tmp = hel_3d(i,j,k)
               hel_3d(i,j,k) = her_3d(i,j,k)
               her_3d(i,j,k) = he_tmp
            enddo
         enddo
      endif

      do k = 1, kdm
         km = k + mm
         kn = k + nn
         do j = 1, jdm

            ! Extract 1D arrays from multidimensional arrays.
            do i = 1, idm + 1
               ca(i) = cau(i,j,k)
               db(i) = pbu(i,j,n)
            enddo
            do i = 0, idm + 1
               du(i) = p(i,j,k  )
               dl(i) = p(i,j,k+1)
            enddo
            do i = -2, idm + 3
               ai(i) = scp2i(i,j)
               ho(i) = max(c0, dp(i,j,kn)) + dpeps
               hm(i) = ho(i)
               hel(i) = hel_3d(i,j,k)
               her(i) = her_3d(i,j,k)
               tm(1,i) = temp(i,j,kn)
               tm(2,i) = saln(i,j,kn)
               do nt = 3, ntr_loc
                  tm(nt,i) = trc(i,j,kn,nt-2)
               enddo
            enddo

            if (second_pass) then
               ! For the second pass of the split scheme, apply a divergence
               ! correction to the thickness used for subsequent flux
               ! estimations in order to preserve an initially uniform thickness
               ! in non-divergent flow.
               do i = -2, idm + 3
                  hm(i) = hm(i)/(c1 - (cav(i,j+1,k) - cav(i,j,k))*ai(i))
               enddo
            endif

            ! Reconstruct tracer edge values, apply limiting, define parabolas
            ! and obtain fluxes by integrating over the flux area.

            call parabola_coeffs_fc_mono(idm, 0, idm + 1, &
                                         stencili(:,j), &
                                         tmc0i(:,:,j), &
                                         tmcli(:,:,j), &
                                         tmcri(:,:,j), &
                                         ssci(:,j), scci(:,j), &
                                         hm, tm, hel, her, &
                                         hpc0, hpc1, hpc2, tpc0, tpc1, tpc2)

            call flux_integration(idm, 1, idm + 1, &
                                  ca, ai, db, du, dl, hpc0, hpc1, hpc2, &
                                  tpc0, tpc1, tpc2, &
                                  hf, htf)

            ! Update variables with flux divergences.
            do i = 1, idm
               hn = ho(i) - (hf(i+1) - hf(i))*ai(i)
               hni = c1/hn
               temp(i,j,kn) = ( ho(i)*tm(1,i) &
                              - (htf(1,i+1) - htf(1,i))*ai(i))*hni
               saln(i,j,kn) = ( ho(i)*tm(2,i) &
                              - (htf(2,i+1) - htf(2,i))*ai(i))*hni
               do nt = 3, ntr_loc
                  trc(i,j,kn,nt-2) = ( ho(i)*tm(nt,i) &
                                     - (htf(nt,i+1) - htf(nt,i))*ai(i))*hni
               enddo
               dp(i,j,kn) = max(c0, hn - dpeps)
            enddo

            ! Accumulate fluxes.
            do i = 1, idm + 1
               uflx (i,j,km) = uflx (i,j,km) + hf(i)
               utflx(i,j,km) = utflx(i,j,km) + htf(1,i)
               usflx(i,j,km) = usflx(i,j,km) + htf(2,i)
            enddo

         enddo
      enddo

   end subroutine cppm_fc_mono_i

   subroutine cppm_fc_mono_j(m, n, mm, nn, k1m, k1n, second_pass)
   ! ---------------------------------------------------------------------------
   ! Do CPPM transport in the j-direction, applying monotonic limiting.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n
      logical, intent(in) :: second_pass

      real(r8), dimension(1-nbdy:jdm+nbdy) :: &
         db, dl, du, ca, ai, ho, hm, hel, her, hpc0, hpc1, hpc2, hf
      real(r8), dimension(ntr_loc,1-nbdy:jdm+nbdy) :: tm, tpc0, tpc1, tpc2, htf

      real(r8) :: he_tmp, hn, hni
      integer :: i, j, k, km, kn, nt

      call xctilr(dp(1-nbdy,1-nbdy,k1n), 1, kdm, 0, 3, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,k1n), 1, kdm, 0, 3, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1n), 1, kdm, 0, 3, halo_ps)
      do nt = 3, ntr_loc
         call xctilr(trc(1-nbdy,1-nbdy,k1n,nt-2), 1, kdm, 0, 3, halo_ps)
      enddo

      do k = 1, kdm
         kn = k + nn
         do i = 1, idm

            ! Extract 1D arrays from multidimensional arrays.
            do j = -2, jdm + 3
               ai(j) = scp2i(i,j)
               hm(j) = max(c0, dp(i,j,kn)) + dpeps
            enddo

            if (second_pass) then
               ! For the second pass of the split scheme, apply a divergence
               ! correction to the thickness used for subsequent flux
               ! estimations in order to preserve an initially uniform thickness
               ! in non-divergent flow.
               do j = -2, jdm + 3
                  hm(j) = hm(j)/(c1 - (cau(i+1,j,k) - cau(i,j,k))*ai(j))
               enddo
            endif

            ! Reconstruct thickness edge values and apply limiting.
            call h_edges_mono(jdm, 1, jdm, &
                              hevc1j(:,i), hevc2j(:,i), &
                              hevc3j(:,i), hevc4j(:,i), &
                              sscj(:,i), sccj(:,i), &
                              hm, hel, her)

            ! Copy edge values to 3D arrays.
            do j = 1, jdm
               hel_3d(i,j,k) = hel(j)
               her_3d(i,j,k) = her(j)
            enddo

         enddo
      enddo

      call xctilr(hel_3d, 1, kdm, 0, 3, halo_ps)
      call xctilr(her_3d, 1, kdm, 0, 3, halo_ps)

      ! With arctic patch, swap the order of thickness edge values.
      if (nreg == 2 .and. nproc == jpr) then
         do k = 1, kdm
            j = jj
            do i = max(1, itdm/2 - i0 + 1), ii
               he_tmp = hel_3d(i,j,k)
               hel_3d(i,j,k) = her_3d(i,j,k)
               her_3d(i,j,k) = he_tmp
            enddo
            do j = jj + 1, jj + 3
               do i = 1, ii
                  he_tmp = hel_3d(i,j,k)
                  hel_3d(i,j,k) = her_3d(i,j,k)
                  her_3d(i,j,k) = he_tmp
               enddo
            enddo
         enddo
      endif

      do k = 1, kdm
         km = k + mm
         kn = k + nn
         do i = 1, idm

            ! Extract 1D arrays from multidimensional arrays.
            do j = 1, jdm + 1
               ca(j) = cav(i,j,k)
               db(j) = pbv(i,j,n)
            enddo
            do j = 0, jdm + 1
               du(j) = p(i,j,k  )
               dl(j) = p(i,j,k+1)
            enddo
            do j = -2, jdm + 3
               ai(j) = scp2i(i,j)
               ho(j) = max(c0, dp(i,j,kn)) + dpeps
               hm(j) = ho(j)
               hel(j) = hel_3d(i,j,k)
               her(j) = her_3d(i,j,k)
               tm(1,j) = temp(i,j,kn)
               tm(2,j) = saln(i,j,kn)
               do nt = 3, ntr_loc
                  tm(nt,j) = trc(i,j,kn,nt-2)
               enddo
            enddo

            if (second_pass) then
               ! For the second pass of the split scheme, apply a divergence
               ! correction to the thickness used for subsequent flux
               ! estimations in order to preserve an initially uniform thickness
               ! in non-divergent flow.
               do j = -2, jdm + 3
                  hm(j) = hm(j)/(c1 - (cau(i+1,j,k) - cau(i,j,k))*ai(j))
               enddo
            endif

            ! Reconstruct tracer edge values, apply limiting, define parabolas
            ! and obtain fluxes by integrating over the flux area.

            call parabola_coeffs_fc_mono(jdm, 0, jdm + 1, &
                                         stencilj(:,i), &
                                         tmc0j(:,:,i), &
                                         tmclj(:,:,i), &
                                         tmcrj(:,:,i), &
                                         sscj(:,i), sccj(:,i), &
                                         hm, tm, hel, her, &
                                         hpc0, hpc1, hpc2, tpc0, tpc1, tpc2)

            call flux_integration(jdm, 1, jdm + 1, &
                                  ca, ai, db, du, dl, hpc0, hpc1, hpc2, &
                                  tpc0, tpc1, tpc2, &
                                  hf, htf)

            ! Update variables with flux divergences.
            do j = 1, jdm
               hn = ho(j) - (hf(j+1) - hf(j))*ai(j)
               hni = c1/hn
               temp(i,j,kn) = ( ho(j)*tm(1,j) &
                              - (htf(1,j+1) - htf(1,j))*ai(j))*hni
               saln(i,j,kn) = ( ho(j)*tm(2,j) &
                              - (htf(2,j+1) - htf(2,j))*ai(j))*hni
               do nt = 3, ntr_loc
                  trc(i,j,kn,nt-2) = ( ho(j)*tm(nt,j) &
                                     - (htf(nt,j+1) - htf(nt,j))*ai(j))*hni
               enddo
               dp(i,j,kn) = max(c0, hn - dpeps)
            enddo

            ! Accumulate fluxes.
            do j = 1, jdm + 1
               vflx (i,j,km) = vflx (i,j,km) + hf(j)
               vtflx(i,j,km) = vtflx(i,j,km) + htf(1,j)
               vsflx(i,j,km) = vsflx(i,j,km) + htf(2,j)
            enddo

         enddo
      enddo

   end subroutine cppm_fc_mono_j

   subroutine cppm_pc_nosc_i(m, n, mm, nn, k1m, k1n, second_pass)
   ! ---------------------------------------------------------------------------
   ! Do CPPM transport in the i-direction, applying non-oscillatory limiting.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n
      logical, intent(in) :: second_pass

      real(r8), dimension(1-nbdy:idm+nbdy) :: &
         db, dl, du, ca, ai, ho, hm, hel, her, hpc0, hpc1, hpc2, hf
      real(r8), dimension(ntr_loc,1-nbdy:idm+nbdy) :: tm, tpc0, tpc1, tpc2, htf

      real(r8) :: hn, hni
      integer :: i, j, k, km, kn, nt

      call xctilr(dp(1-nbdy,1-nbdy,k1n), 1, kdm, 4, 0, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,k1n), 1, kdm, 4, 0, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1n), 1, kdm, 4, 0, halo_ps)
      do nt = 3, ntr_loc
         call xctilr(trc(1-nbdy,1-nbdy,k1n,nt-2), 1, kdm, 4, 0, halo_ps)
      enddo

      do k = 1, kdm
         km = k + mm
         kn = k + nn
         do j = 1, jdm

            ! Extract 1D arrays from multidimensional arrays.
            do i = 1, idm + 1
               ca(i) = cau(i,j,k)
               db(i) = pbu(i,j,n)
            enddo
            do i = 0, idm + 1
               du(i) = p(i,j,k  )
               dl(i) = p(i,j,k+1)
            enddo
            do i = -3, idm + 4
               ai(i) = scp2i(i,j)
               ho(i) = max(c0, dp(i,j,kn)) + dpeps
               hm(i) = ho(i)
               tm(1,i) = temp(i,j,kn)
               tm(2,i) = saln(i,j,kn)
               do nt = 3, ntr_loc
                  tm(nt,i) = trc(i,j,kn,nt-2)
               enddo
            enddo

            if (second_pass) then
               ! For the second pass of the split scheme, apply a divergence
               ! correction to the thickness used for subsequent flux
               ! estimations in order to preserve an initially uniform thickness
               ! in non-divergent flow.
               do i = -3, idm + 4
                  hm(i) = hm(i)/(c1 - (cav(i,j+1,k) - cav(i,j,k))*ai(i))
               enddo
            endif

            ! Reconstruct tracer edge values, apply limiting, define parabolas
            ! and obtain fluxes by integrating over the flux area.

            call parabola_coeffs_pc_nosc(idm, 0, idm + 1, &
                                         hevc1i(:,j), hevc2i(:,j), &
                                         hevc3i(:,j), hevc4i(:,j), &
                                         ssci(:,j), scci(:,j), &
                                         d2mi(:,j), &
                                         hm, tm, &
                                         hpc0, hpc1, hpc2, tpc0, tpc1, tpc2)

            call flux_integration(idm, 1, idm + 1, &
                                  ca, ai, db, du, dl, hpc0, hpc1, hpc2, &
                                  tpc0, tpc1, tpc2, &
                                  hf, htf)

            ! Update variables with flux divergences.
            do i = 1, idm
               hn = ho(i) - (hf(i+1) - hf(i))*ai(i)
               hni = c1/hn
               temp(i,j,kn) = ( ho(i)*tm(1,i) &
                              - (htf(1,i+1) - htf(1,i))*ai(i))*hni
               saln(i,j,kn) = ( ho(i)*tm(2,i) &
                              - (htf(2,i+1) - htf(2,i))*ai(i))*hni
               do nt = 3, ntr_loc
                  trc(i,j,kn,nt-2) = ( ho(i)*tm(nt,i) &
                                     - (htf(nt,i+1) - htf(nt,i))*ai(i))*hni
               enddo
               dp(i,j,kn) = max(c0, hn - dpeps)
            enddo

            ! Accumulate fluxes.
            do i = 1, idm + 1
               uflx (i,j,km) = uflx (i,j,km) + hf(i)
               utflx(i,j,km) = utflx(i,j,km) + htf(1,i)
               usflx(i,j,km) = usflx(i,j,km) + htf(2,i)
            enddo

         enddo
      enddo

   end subroutine cppm_pc_nosc_i

   subroutine cppm_pc_nosc_j(m, n, mm, nn, k1m, k1n, second_pass)
   ! ---------------------------------------------------------------------------
   ! Do CPPM transport in the j-direction, applying non-oscillatory limiting.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n
      logical, intent(in) :: second_pass

      real(r8), dimension(1-nbdy:jdm+nbdy) :: &
         db, dl, du, ca, ai, ho, hm, hel, her, hpc0, hpc1, hpc2, hf
      real(r8), dimension(ntr_loc,1-nbdy:jdm+nbdy) :: tm, tpc0, tpc1, tpc2, htf

      real(r8) :: hn, hni
      integer :: i, j, k, km, kn, nt

      call xctilr(dp(1-nbdy,1-nbdy,k1n), 1, kdm, 0, 4, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,k1n), 1, kdm, 0, 4, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1n), 1, kdm, 0, 4, halo_ps)
      do nt = 3, ntr_loc
         call xctilr(trc(1-nbdy,1-nbdy,k1n,nt-2), 1, kdm, 0, 4, halo_ps)
      enddo

      do k = 1, kdm
         km = k + mm
         kn = k + nn
         do i = 1, idm

            ! Extract 1D arrays from multidimensional arrays.
            do j = 1, jdm + 1
               ca(j) = cav(i,j,k)
               db(j) = pbv(i,j,n)
            enddo
            do j = 0, jdm + 1
               du(j) = p(i,j,k  )
               dl(j) = p(i,j,k+1)
            enddo
            do j = -3, jdm + 4
               ai(j) = scp2i(i,j)
               ho(j) = max(c0, dp(i,j,kn)) + dpeps
               hm(j) = ho(j)
               tm(1,j) = temp(i,j,kn)
               tm(2,j) = saln(i,j,kn)
               do nt = 3, ntr_loc
                  tm(nt,j) = trc(i,j,kn,nt-2)
               enddo
            enddo

            if (second_pass) then
               ! For the second pass of the split scheme, apply a divergence
               ! correction to the thickness used for subsequent flux
               ! estimations in order to preserve an initially uniform thickness
               ! in non-divergent flow.
               do j = -3, jdm + 4
                  hm(j) = hm(j)/(c1 - (cau(i+1,j,k) - cau(i,j,k))*ai(j))
               enddo
            endif

            ! Reconstruct tracer edge values, apply limiting, define parabolas
            ! and obtain fluxes by integrating over the flux area.

            call parabola_coeffs_pc_nosc(jdm, 0, jdm + 1, &
                                         hevc1j(:,i), hevc2j(:,i), &
                                         hevc3j(:,i), hevc4j(:,i), &
                                         sscj(:,i), sccj(:,i), &
                                         d2mj(:,i), &
                                         hm, tm, &
                                         hpc0, hpc1, hpc2, tpc0, tpc1, tpc2)

            call flux_integration(jdm, 1, jdm + 1, &
                                  ca, ai, db, du, dl, hpc0, hpc1, hpc2, &
                                  tpc0, tpc1, tpc2, &
                                  hf, htf)

            ! Update variables with flux divergences.
            do j = 1, jdm
               hn = ho(j) - (hf(j+1) - hf(j))*ai(j)
               hni = c1/hn
               temp(i,j,kn) = ( ho(j)*tm(1,j) &
                              - (htf(1,j+1) - htf(1,j))*ai(j))*hni
               saln(i,j,kn) = ( ho(j)*tm(2,j) &
                              - (htf(2,j+1) - htf(2,j))*ai(j))*hni
               do nt = 3, ntr_loc
                  trc(i,j,kn,nt-2) = ( ho(j)*tm(nt,j) &
                                     - (htf(nt,j+1) - htf(nt,j))*ai(j))*hni
               enddo
               dp(i,j,kn) = max(c0, hn - dpeps)
            enddo

            ! Accumulate fluxes.
            do j = 1, jdm + 1
               vflx (i,j,km) = vflx (i,j,km) + hf(j)
               vtflx(i,j,km) = vtflx(i,j,km) + htf(1,j)
               vsflx(i,j,km) = vsflx(i,j,km) + htf(2,j)
            enddo

         enddo
      enddo

   end subroutine cppm_pc_nosc_j

   subroutine cppm_pc_mono_i(m, n, mm, nn, k1m, k1n, second_pass)
   ! ---------------------------------------------------------------------------
   ! Do CPPM transport in the i-direction, applying non-oscillatory limiting.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n
      logical, intent(in) :: second_pass

      real(r8), dimension(1-nbdy:idm+nbdy) :: &
         db, dl, du, ca, ai, ho, hm, hel, her, hpc0, hpc1, hpc2, hf
      real(r8), dimension(ntr_loc,1-nbdy:idm+nbdy) :: tm, tpc0, tpc1, tpc2, htf

      real(r8) :: hn, hni
      integer :: i, j, k, km, kn, nt

      call xctilr(dp(1-nbdy,1-nbdy,k1n), 1, kdm, 3, 0, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,k1n), 1, kdm, 3, 0, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1n), 1, kdm, 3, 0, halo_ps)
      do nt = 3, ntr_loc
         call xctilr(trc(1-nbdy,1-nbdy,k1n,nt-2), 1, kdm, 3, 0, halo_ps)
      enddo

      do k = 1, kdm
         km = k + mm
         kn = k + nn
         do j = 1, jdm

            ! Extract 1D arrays from multidimensional arrays.
            do i = 1, idm + 1
               ca(i) = cau(i,j,k)
               db(i) = pbu(i,j,n)
            enddo
            do i = 0, idm + 1
               du(i) = p(i,j,k  )
               dl(i) = p(i,j,k+1)
            enddo
            do i = -2, idm + 3
               ai(i) = scp2i(i,j)
               ho(i) = max(c0, dp(i,j,kn)) + dpeps
               hm(i) = ho(i)
               tm(1,i) = temp(i,j,kn)
               tm(2,i) = saln(i,j,kn)
               do nt = 3, ntr_loc
                  tm(nt,i) = trc(i,j,kn,nt-2)
               enddo
            enddo

            if (second_pass) then
               ! For the second pass of the split scheme, apply a divergence
               ! correction to the thickness used for subsequent flux
               ! estimations in order to preserve an initially uniform thickness
               ! in non-divergent flow.
               do i = -2, idm + 3
                  hm(i) = hm(i)/(c1 - (cav(i,j+1,k) - cav(i,j,k))*ai(i))
               enddo
            endif

            ! Reconstruct tracer edge values, apply limiting, define parabolas
            ! and obtain fluxes by integrating over the flux area.

            call parabola_coeffs_pc_mono(idm, 0, idm + 1, &
                                         hevc1i(:,j), hevc2i(:,j), &
                                         hevc3i(:,j), hevc4i(:,j), &
                                         ssci(:,j), scci(:,j), &
                                         hm, tm, &
                                         hpc0, hpc1, hpc2, tpc0, tpc1, tpc2)

            call flux_integration(idm, 1, idm + 1, &
                                  ca, ai, db, du, dl, hpc0, hpc1, hpc2, &
                                  tpc0, tpc1, tpc2, &
                                  hf, htf)

            ! Update variables with flux divergences.
            do i = 1, idm
               hn = ho(i) - (hf(i+1) - hf(i))*ai(i)
               hni = c1/hn
               temp(i,j,kn) = ( ho(i)*tm(1,i) &
                              - (htf(1,i+1) - htf(1,i))*ai(i))*hni
               saln(i,j,kn) = ( ho(i)*tm(2,i) &
                              - (htf(2,i+1) - htf(2,i))*ai(i))*hni
               do nt = 3, ntr_loc
                  trc(i,j,kn,nt-2) = ( ho(i)*tm(nt,i) &
                                     - (htf(nt,i+1) - htf(nt,i))*ai(i))*hni
               enddo
               dp(i,j,kn) = max(c0, hn - dpeps)
            enddo

            ! Accumulate fluxes.
            do i = 1, idm + 1
               uflx (i,j,km) = uflx (i,j,km) + hf(i)
               utflx(i,j,km) = utflx(i,j,km) + htf(1,i)
               usflx(i,j,km) = usflx(i,j,km) + htf(2,i)
            enddo

         enddo
      enddo

   end subroutine cppm_pc_mono_i

   subroutine cppm_pc_mono_j(m, n, mm, nn, k1m, k1n, second_pass)
   ! ---------------------------------------------------------------------------
   ! Do CPPM transport in the j-direction, applying non-oscillatory limiting.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n
      logical, intent(in) :: second_pass

      real(r8), dimension(1-nbdy:jdm+nbdy) :: &
         db, dl, du, ca, ai, ho, hm, hel, her, hpc0, hpc1, hpc2, hf
      real(r8), dimension(ntr_loc,1-nbdy:jdm+nbdy) :: tm, tpc0, tpc1, tpc2, htf

      real(r8) :: hn, hni
      integer :: i, j, k, km, kn, nt

      call xctilr(dp(1-nbdy,1-nbdy,k1n), 1, kdm, 0, 3, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,k1n), 1, kdm, 0, 3, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k1n), 1, kdm, 0, 3, halo_ps)
      do nt = 3, ntr_loc
         call xctilr(trc(1-nbdy,1-nbdy,k1n,nt-2), 1, kdm, 0, 3, halo_ps)
      enddo

      do k = 1, kdm
         km = k + mm
         kn = k + nn
         do i = 1, idm

            ! Extract 1D arrays from multidimensional arrays.
            do j = 1, jdm + 1
               ca(j) = cav(i,j,k)
               db(j) = pbv(i,j,n)
            enddo
            do j = 0, jdm + 1
               du(j) = p(i,j,k  )
               dl(j) = p(i,j,k+1)
            enddo
            do j = -2, jdm + 3
               ai(j) = scp2i(i,j)
               ho(j) = max(c0, dp(i,j,kn)) + dpeps
               hm(j) = ho(j)
               tm(1,j) = temp(i,j,kn)
               tm(2,j) = saln(i,j,kn)
               do nt = 3, ntr_loc
                  tm(nt,j) = trc(i,j,kn,nt-2)
               enddo
            enddo

            if (second_pass) then
               ! For the second pass of the split scheme, apply a divergence
               ! correction to the thickness used for subsequent flux
               ! estimations in order to preserve an initially uniform thickness
               ! in non-divergent flow.
               do j = -2, jdm + 3
                  hm(j) = hm(j)/(c1 - (cau(i+1,j,k) - cau(i,j,k))*ai(j))
               enddo
            endif

            ! Reconstruct tracer edge values, apply limiting, define parabolas
            ! and obtain fluxes by integrating over the flux area.

            call parabola_coeffs_pc_mono(jdm, 0, jdm + 1, &
                                         hevc1j(:,i), hevc2j(:,i), &
                                         hevc3j(:,i), hevc4j(:,i), &
                                         sscj(:,i), sccj(:,i), &
                                         hm, tm, &
                                         hpc0, hpc1, hpc2, tpc0, tpc1, tpc2)

            call flux_integration(jdm, 1, jdm + 1, &
                                  ca, ai, db, du, dl, hpc0, hpc1, hpc2, &
                                  tpc0, tpc1, tpc2, &
                                  hf, htf)

            ! Update variables with flux divergences.
            do j = 1, jdm
               hn = ho(j) - (hf(j+1) - hf(j))*ai(j)
               hni = c1/hn
               temp(i,j,kn) = ( ho(j)*tm(1,j) &
                              - (htf(1,j+1) - htf(1,j))*ai(j))*hni
               saln(i,j,kn) = ( ho(j)*tm(2,j) &
                              - (htf(2,j+1) - htf(2,j))*ai(j))*hni
               do nt = 3, ntr_loc
                  trc(i,j,kn,nt-2) = ( ho(j)*tm(nt,j) &
                                     - (htf(nt,j+1) - htf(nt,j))*ai(j))*hni
               enddo
               dp(i,j,kn) = max(c0, hn - dpeps)
            enddo

            ! Accumulate fluxes.
            do j = 1, jdm + 1
               vflx (i,j,km) = vflx (i,j,km) + hf(j)
               vtflx(i,j,km) = vtflx(i,j,km) + htf(1,j)
               vsflx(i,j,km) = vsflx(i,j,km) + htf(2,j)
            enddo

         enddo
      enddo

   end subroutine cppm_pc_mono_j

   ! ---------------------------------------------------------------------------
   ! Public procedures.
   ! ---------------------------------------------------------------------------

   subroutine init_cppm
   ! ---------------------------------------------------------------------------
   ! Resolve namelist options and set various coefficients and masks for CPPM.
   ! ---------------------------------------------------------------------------

      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: stencilj_perm
      real(r8), dimension(12,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
         tmc0j_perm, tmclj_perm, tmcrj_perm
      real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
         hevc1j_perm, hevc2j_perm, hevc3j_perm, hevc4j_perm, &
         sscj_perm, sccj_perm, d2mj_perm, tmp2d

      real(r8), dimension(4) :: dx4
      real(r8), dimension(3) :: dx3
      integer, dimension(4) :: sm4
      integer, dimension(3) :: sm3
      real(r8) :: hevc_tmp
      integer :: i, j, k

      ! Resolve namelist options.
      select case (trim(cppm_compatibility))
         case ('full')
            cppm_compatibility_tag = cppm_compatibility_full
         case ('partial')
            cppm_compatibility_tag = cppm_compatibility_partial
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') &
                  ' init_cppm: cppm_compatibility = ', &
                  trim(cppm_compatibility), ' is unsupported!'
            call xcstop('(init_cppm)')
            stop '(init_cppm)'
      end select
      select case (trim(cppm_limiting))
         case ('monotonic')
            cppm_limiting_tag = cppm_limiting_monotonic
         case ('non_oscillatory')
            cppm_limiting_tag = cppm_limiting_non_oscillatory
         case default
            if (mnproc == 1) &
               write (lp,'(3a)') &
                  ' init_cppm: cppm_limiting = ', trim(cppm_limiting), &
                  ' is unsupported!'
            call xcstop('(init_cppm)')
            stop '(init_cppm)'
      end select

      ! Set stencil tag, coefficients for thickness edge value, coefficients for
      ! slope reconstructions and mask for the value proportional to the second
      ! derivative of unlimited parabolas.

      stencili(:,:) = stencil_0000
      hevc1i(:,:) = c0
      hevc2i(:,:) = c0
      hevc3i(:,:) = c0
      hevc4i(:,:) = c0
      tmc0i(:,:,:) = c0
      tmcli(:,:,:) = c0
      tmcri(:,:,:) = c0
      ssci(:,:) = c0
      scci(:,:) = c0
      d2mi(:,:) = c0
      stencilj_perm(:,:) = stencil_0000
      hevc1j_perm(:,:) = c0
      hevc2j_perm(:,:) = c0
      hevc3j_perm(:,:) = c0
      hevc4j_perm(:,:) = c0
      tmc0j_perm(:,:,:) = c0
      tmclj_perm(:,:,:) = c0
      tmcrj_perm(:,:,:) = c0
      sscj_perm(:,:) = c0
      sccj_perm(:,:) = c0
      d2mj_perm(:,:) = c0

      do j = 1, jj
         do i = 1, ii
            sm4 = ip(i-2:i+1,j)
            dx4 = scpx(i-2:i+1,j)
            call set_stencil_coeffs(sm4, dx4, stencili(i,j), &
                                    hevc1i(i,j), hevc2i(i,j), &
                                    hevc3i(i,j), hevc4i(i,j), &
                                    tmc0i(:,i,j), tmcli(:,i,j), tmcri(:,i,j))
            sm3 = ip(i-1:i+1,j)
            dx3 = scpx(i-1:i+1,j)
            call set_slope_coeffs(sm3, dx3, ssci(i,j), scci(i,j))
            call set_d2_mask(sm3, d2mi(i,j))
            sm4 = ip(i,j-2:j+1)
            dx4 = scpy(i,j-2:j+1)
            call set_stencil_coeffs(sm4, dx4, stencilj_perm(i,j), &
                                    hevc1j_perm(i,j), hevc2j_perm(i,j), &
                                    hevc3j_perm(i,j), hevc4j_perm(i,j), &
                                    tmc0j_perm(:,i,j), &
                                    tmclj_perm(:,i,j), &
                                    tmcrj_perm(:,i,j))
            sm3 = ip(i,j-1:j+1)
            dx3 = scpy(i,j-1:j+1)
            call set_slope_coeffs(sm3, dx3, sscj_perm(i,j), sccj_perm(i,j))
            call set_d2_mask(sm3, d2mj_perm(i,j))
         enddo
      enddo

      tmp2d(:,:) = real(stencili(:,:), r8)
      call xctilr(tmp2d, 1, 1, nbdy, 0, halo_us)
      stencili(:,:) = nint(tmp2d(:,:))
      call xctilr(hevc1i, 1, 1, nbdy, 0, halo_us)
      call xctilr(hevc2i, 1, 1, nbdy, 0, halo_us)
      call xctilr(hevc3i, 1, 1, nbdy, 0, halo_us)
      call xctilr(hevc4i, 1, 1, nbdy, 0, halo_us)
      do k = 1, 12
         tmp2d(:,:) = tmc0i(k,:,:)
         call xctilr(tmp2d, 1, 1, nbdy, 0, halo_us)
         tmc0i(k,:,:) = tmp2d(:,:)
         tmp2d(:,:) = tmcli(k,:,:)
         call xctilr(tmp2d, 1, 1, nbdy, 0, halo_us)
         tmcli(k,:,:) = tmp2d(:,:)
         tmp2d(:,:) = tmcri(k,:,:)
         call xctilr(tmp2d, 1, 1, nbdy, 0, halo_us)
         tmcri(k,:,:) = tmp2d(:,:)
      enddo
      call xctilr(ssci, 1, 1, nbdy, 0, halo_ps)
      call xctilr(scci, 1, 1, nbdy, 0, halo_ps)
      call xctilr(d2mi, 1, 1, nbdy, 0, halo_ps)
      tmp2d(:,:) = real(stencilj_perm(:,:), r8)
      call xctilr(tmp2d, 1, 1, 0, nbdy, halo_vs)
      stencilj_perm(:,:) = nint(tmp2d(:,:))
      call xctilr(hevc1j_perm, 1, 1, 0, nbdy, halo_vs)
      call xctilr(hevc2j_perm, 1, 1, 0, nbdy, halo_vs)
      call xctilr(hevc3j_perm, 1, 1, 0, nbdy, halo_vs)
      call xctilr(hevc4j_perm, 1, 1, 0, nbdy, halo_vs)
      do k = 1, 12
         tmp2d(:,:) = tmc0j_perm(k,:,:)
         call xctilr(tmp2d, 1, 1, 0, nbdy, halo_vs)
         tmc0j_perm(k,:,:) = tmp2d(:,:)
         tmp2d(:,:) = tmclj_perm(k,:,:)
         call xctilr(tmp2d, 1, 1, 0, nbdy, halo_vs)
         tmclj_perm(k,:,:) = tmp2d(:,:)
         tmp2d(:,:) = tmcrj_perm(k,:,:)
         call xctilr(tmp2d, 1, 1, 0, nbdy, halo_vs)
         tmcrj_perm(k,:,:) = tmp2d(:,:)
      enddo
      call xctilr(sscj_perm, 1, 1, 0, nbdy, halo_ps)
      call xctilr(sccj_perm, 1, 1, 0, nbdy, halo_ps)
      call xctilr(d2mj_perm, 1, 1, 0, nbdy, halo_ps)

      ! With arctic patch, swap the order of some stencil tags and edge value
      ! coefficients.
      if (nreg == 2 .and. nproc == jpr) then
         j = jj
         do i = 1 - nbdy, ii + nbdy
            select case(stencili(i,j))
               case (stencil_1110)
                  stencili(i,j) = stencil_0111
               case (stencil_0111)
                  stencili(i,j) = stencil_1110
               case (stencil_1100)
                  stencili(i,j) = stencil_0011
               case (stencil_0011)
                  stencili(i,j) = stencil_1100
               case (stencil_0100)
                  stencili(i,j) = stencil_0010
               case (stencil_0010)
                  stencili(i,j) = stencil_0100
            end select
            hevc_tmp = hevc1i(i,j)
            hevc1i(i,j) = hevc4i(i,j)
            hevc4i(i,j) = hevc_tmp
            hevc_tmp = hevc2i(i,j)
            hevc2i(i,j) = hevc3i(i,j)
            hevc3i(i,j) = hevc_tmp
         enddo
         do i = max(1, itdm/2 - i0 + 1), ii
            select case(stencilj_perm(i,j))
               case (stencil_1110)
                  stencilj_perm(i,j) = stencil_0111
               case (stencil_0111)
                  stencilj_perm(i,j) = stencil_1110
               case (stencil_1100)
                  stencilj_perm(i,j) = stencil_0011
               case (stencil_0011)
                  stencilj_perm(i,j) = stencil_1100
               case (stencil_0100)
                  stencilj_perm(i,j) = stencil_0010
               case (stencil_0010)
                  stencilj_perm(i,j) = stencil_0100
            end select
            hevc_tmp = hevc1j_perm(i,j)
            hevc1j_perm(i,j) = hevc4j_perm(i,j)
            hevc4j_perm(i,j) = hevc_tmp
            hevc_tmp = hevc2j_perm(i,j)
            hevc2j_perm(i,j) = hevc3j_perm(i,j)
            hevc3j_perm(i,j) = hevc_tmp
         enddo
         do j = jj + 1, jj + nbdy
            do i = 1, ii
               select case(stencilj_perm(i,j))
                  case (stencil_1110)
                     stencilj_perm(i,j) = stencil_0111
                  case (stencil_0111)
                     stencilj_perm(i,j) = stencil_1110
                  case (stencil_1100)
                     stencilj_perm(i,j) = stencil_0011
                  case (stencil_0011)
                     stencilj_perm(i,j) = stencil_1100
                  case (stencil_0100)
                     stencilj_perm(i,j) = stencil_0010
                  case (stencil_0010)
                     stencilj_perm(i,j) = stencil_0100
               end select
               hevc_tmp = hevc1j_perm(i,j)
               hevc1j_perm(i,j) = hevc4j_perm(i,j)
               hevc4j_perm(i,j) = hevc_tmp
               hevc_tmp = hevc2j_perm(i,j)
               hevc2j_perm(i,j) = hevc3j_perm(i,j)
               hevc3j_perm(i,j) = hevc_tmp
            enddo
         enddo
      endif

      do j = 1 - nbdy, jdm + nbdy
         do i = 1 - nbdy, idm + nbdy
            stencilj(j,i) = stencilj_perm(i,j)
            hevc1j(j,i) = hevc1j_perm(i,j)
            hevc2j(j,i) = hevc2j_perm(i,j)
            hevc3j(j,i) = hevc3j_perm(i,j)
            hevc4j(j,i) = hevc4j_perm(i,j)
            tmc0j(:,j,i) = tmc0j_perm(:,i,j)
            tmclj(:,j,i) = tmclj_perm(:,i,j)
            tmcrj(:,j,i) = tmcrj_perm(:,i,j)
            sscj(j,i) = sscj_perm(i,j)
            sccj(j,i) = sccj_perm(i,j)
            d2mj(j,i) = d2mj_perm(i,j)
         enddo
      enddo

      ! Initialize thickness edge value arrays.
      hel_3d(:,:,:) = c0
      her_3d(:,:,:) = c0

      ! Local number of tracers where temperature and salinity is added to the
      ! ntr parameter.
      ntr_loc = 2 + ntr

   end subroutine init_cppm

   subroutine cppm(m, n, mm, nn, k1m, k1n)
   ! ---------------------------------------------------------------------------
   ! To complete a 2D transport step using CPPM, call 1D transport operators in
   ! alternate directional order (Strang splitting).
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: m, n, mm, nn, k1m, k1n

      if (cppm_compatibility_tag == cppm_compatibility_full) then

         if (cppm_limiting_tag == cppm_limiting_non_oscillatory) then

            call xctilr(cau, 1, kdm, 4, 4, halo_uv)
            call xctilr(cav, 1, kdm, 4, 4, halo_vv)

            if (mod(nstep,2) == 1) then
               call cppm_fc_nosc_i(m, n, mm, nn, k1m, k1n, &
                                   second_pass = .false.)
               call cppm_fc_nosc_j(m, n, mm, nn, k1m, k1n, &
                                   second_pass = .true. )
            else
               call cppm_fc_nosc_j(m, n, mm, nn, k1m, k1n, &
                                   second_pass = .false.)
               call cppm_fc_nosc_i(m, n, mm, nn, k1m, k1n, &
                                   second_pass = .true. )
            endif

         else

            call xctilr(cau, 1, kdm, 3, 3, halo_uv)
            call xctilr(cav, 1, kdm, 3, 3, halo_vv)

            if (mod(nstep,2) == 1) then
               call cppm_fc_mono_i(m, n, mm, nn, k1m, k1n, &
                                   second_pass = .false.)
               call cppm_fc_mono_j(m, n, mm, nn, k1m, k1n, &
                                   second_pass = .true. )
            else
               call cppm_fc_mono_j(m, n, mm, nn, k1m, k1n, &
                                   second_pass = .false.)
               call cppm_fc_mono_i(m, n, mm, nn, k1m, k1n, &
                                   second_pass = .true. )
            endif

         endif

      else

         if (cppm_limiting_tag == cppm_limiting_non_oscillatory) then

            call xctilr(cau, 1, kdm, 4, 4, halo_uv)
            call xctilr(cav, 1, kdm, 4, 4, halo_vv)

            if (mod(nstep,2) == 1) then
               call cppm_pc_nosc_i(m, n, mm, nn, k1m, k1n, &
                                   second_pass = .false.)
               call cppm_pc_nosc_j(m, n, mm, nn, k1m, k1n, &
                                   second_pass = .true. )
            else
               call cppm_pc_nosc_j(m, n, mm, nn, k1m, k1n, &
                                   second_pass = .false.)
               call cppm_pc_nosc_i(m, n, mm, nn, k1m, k1n, &
                                   second_pass = .true. )
            endif

         else

            call xctilr(cau, 1, kdm, 3, 3, halo_uv)
            call xctilr(cav, 1, kdm, 3, 3, halo_vv)

            if (mod(nstep,2) == 1) then
               call cppm_pc_mono_i(m, n, mm, nn, k1m, k1n, &
                                   second_pass = .false.)
               call cppm_pc_mono_j(m, n, mm, nn, k1m, k1n, &
                                   second_pass = .true. )
            else
               call cppm_pc_mono_j(m, n, mm, nn, k1m, k1n, &
                                   second_pass = .false.)
               call cppm_pc_mono_i(m, n, mm, nn, k1m, k1n, &
                                   second_pass = .true. )
            endif

         endif

      endif

   end subroutine cppm

end module mod_cppm
