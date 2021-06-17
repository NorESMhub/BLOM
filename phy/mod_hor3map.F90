! ------------------------------------------------------------------------------
! Copyright (C) 2021 Mats Bentsen
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

module mod_hor3map
! ------------------------------------------------------------------------------
! This module contains routines and data structures for High-order
! One-dimensional Reconstruction, Regridding and ReMAPping (HOR3MAP) on
! nonuniform grids. The remapping is done by integrating a piecewise polynomial
! reconstruction of the source data in segments to establish finite volume
! averages on the destination grid. For computational efficiency when remapping
! multiple data with same source and destination grids, stages of the
! reconstruction and remapping are separated into distinct routines. The module
! also have a routine for regridding, where grid locations are found as
! intersections between a polynomial reconstruction of the source data and
! desired grid cell edge data values.
! ------------------------------------------------------------------------------

   use, intrinsic :: iso_fortran_env, only: real64

   implicit none

   private

   ! Option parameters.
   integer, parameter :: &
      hor3map_pcm                    = 100, & ! Reconstruction methods
      hor3map_plm                    = 101, &
      hor3map_ppm                    = 102, &
      hor3map_no_limiting            = 200, & ! Limiting methods
      hor3map_monotonic              = 201, &
      hor3map_non_oscillatory        = 203, &
      hor3map_non_oscillatory_posdef = 204

   ! Error parameters.
   integer, parameter :: &
      hor3map_noerr                   =  0, &
      hor3map_invalid_method          =  1, &
      hor3map_nonmonotonic_src_edges  =  2, &
      hor3map_src_extent_too_small    =  3, &
      hor3map_failed_to_allocate_rcs  =  4, &
      hor3map_recon_not_prepared      =  5, &
      hor3map_inconsistent_grid_range =  6, &
      hor3map_nonmonotonic_dst_edges  =  7, &
      hor3map_failed_to_allocate_rms  =  8, &
      hor3map_src_size_mismatch       =  9, &
      hor3map_invalid_plm_limiting    = 10, &
      hor3map_invalid_ppm_limiting    = 11, &
      hor3map_recon_not_available     = 12, &
      hor3map_grd_size_mismatch       = 13, &
      hor3map_remap_not_prepared      = 14, &
      hor3map_dst_size_mismatch       = 15, &
      hor3map_errmsg_num = 15
   character(len = 80), dimension(hor3map_errmsg_num), parameter :: errmsg = &
      ["Invalid reconstruction method!                                   ", &
       "Source grid edges do not monotonically increase or decrease!     ", &
       "Source grid extent too small!                                    ", &
       "Failed to allocate reconstruction data structure!                ", &
       "Call 'prepare_reconstruction' first!                             ", &
       "Inconsistent source and destination grid range!                  ", &
       "Destination grid edges do not monotonically increase or decrease!", &
       "Failed to allocate remapping data structure!                     ", &
       "Size mismatch between source grid edges and data array!          ", &
       "Invalid limiting method for PLM!                                 ", &
       "Invalid limiting method for PPM!                                 ", &
       "Call 'reconstruct' first!                                        ", &
       "Size mismatch between grid edge values and locations!            ", &
       "Call 'prepare_remapping' first!                                  ", &
       "Size mismatch between destination grid edges and data array!     "]

   ! Numeric data types.
   integer, parameter :: &
      r8 = real64

   ! Numeric constants.
   real(r8), parameter :: &
      eps  = 1.e-14_r8, &    ! Small non-dimensional value.
      c0   = 0._r8, &
      c1   = 1._r8, &
      c2   = 2._r8, &
      c3   = 3._r8, &
      c4   = 4._r8, &
      c1_2 = 1._r8/2._r8, &
      c1_3 = 1._r8/3._r8, &
      c1_4 = 1._r8/4._r8

   type reconstruction_struct
      real(r8), allocatable, dimension(:, :) :: polycoeff
      real(r8), allocatable, dimension(:) :: &
         x_edge_src, h_src, hi_src, hci_src, alpha, beta, b, c, u_src, uel, uer
      real(r8) :: x_eps, al, bl, cl, ar, br, cr
      integer, allocatable, dimension(:) :: src_index_map
      integer :: n_src_all, n_src, method, limiting
      logical :: prepared = .false., reconstructed = .false.
   end type reconstruction_struct

   type remap_struct
      real(r8), allocatable, dimension(:) :: h_dst, hi_dst, seg_int_lim
      integer, allocatable, dimension(:) :: n_src_seg, seg_dst_index
      integer :: n_dst
      logical :: prepared = .false.
   end type remap_struct

   public :: reconstruction_struct, remap_struct, &
             prepare_reconstruction, prepare_remapping, &
             reconstruct, regrid, remap, free_rcs, free_rms, &
             hor3map_pcm, hor3map_plm, hor3map_ppm, &
             hor3map_no_limiting, hor3map_monotonic, hor3map_non_oscillatory, &
             hor3map_non_oscillatory_posdef, &
             hor3map_noerr, hor3map_errstr

contains

   ! ---------------------------------------------------------------------------
   ! Private procedures.
   ! ---------------------------------------------------------------------------

   function allocate_rcs(rcs) result(errstat)
   ! ---------------------------------------------------------------------------
   ! Allocate arrays in the reconstruction data structure.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs

      integer :: errstat

      integer :: allocstat

      allocate(rcs%x_edge_src(rcs%n_src + 1), rcs%h_src(rcs%n_src), &
               rcs%hi_src(rcs%n_src), rcs%hci_src(2:rcs%n_src - 1), &
               rcs%alpha(2:rcs%n_src + 1), rcs%beta(rcs%n_src), &
               rcs%b(rcs%n_src), rcs%c(rcs%n_src), rcs%u_src(rcs%n_src), &
               rcs%uel(rcs%n_src + 1), rcs%uer(rcs%n_src), &
               rcs%src_index_map(rcs%n_src), rcs%polycoeff(3, rcs%n_src), &
               stat = allocstat)

      if (allocstat == 0) then
         errstat = hor3map_noerr
      else
         errstat = hor3map_failed_to_allocate_rcs
         return
      endif

   end function allocate_rcs

   function allocate_rms(rcs, rms) result(errstat)
   ! ---------------------------------------------------------------------------
   ! Allocate arrays in the remapping data structure.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(in) :: rcs
      type(remap_struct), intent(inout) :: rms

      integer :: errstat

      integer :: allocstat

      allocate(rms%h_dst(rms%n_dst), rms%hi_dst(rms%n_dst), &
               rms%seg_int_lim(rcs%n_src + rms%n_dst), &
               rms%n_src_seg(rcs%n_src), &
               rms%seg_dst_index(rcs%n_src + rms%n_dst), &
               stat = allocstat)

      if (allocstat == 0) then
         errstat = hor3map_noerr
      else
         errstat = hor3map_failed_to_allocate_rms
         return
      endif

   end function allocate_rms

   subroutine edge_coeff_ih4(rcs)
   ! ---------------------------------------------------------------------------
   ! Compute coefficients for implicit 4th order accurate edge estimation.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs

      integer :: ns, j

      ns = rcs%n_src

      ! Compute matrix entries and right-hand side coefficients for the left
      ! boundary with an asymmetric stencil.
      call edge_coeff_ih4_asym(rcs%h_src(1), rcs%h_src(2), rcs%h_src(3), &
                               rcs%beta(1), rcs%al, rcs%bl, rcs%cl)

      ! Compute matrix entries and right-hand side coefficients where a
      ! symmetric stencil is applicable.
      do j = 2, ns
         call edge_coeff_ih4_sym(rcs%h_src(j - 1), rcs%h_src(j), &
                                 rcs%alpha(j), rcs%beta(j), rcs%b(j), rcs%c(j))
      enddo

      ! Compute matrix entries and right-hand side coefficients for the right
      ! boundary with an asymmetric stencil.
      call edge_coeff_ih4_asym(rcs%h_src(ns), rcs%h_src(ns - 1), &
                               rcs%h_src(ns - 2), &
                               rcs%alpha(ns + 1), rcs%ar, rcs%br, rcs%cr)

   end subroutine edge_coeff_ih4

   subroutine edge_coeff_ih4_asym(h0, h1, h2, beta, a, b, c)
   ! ---------------------------------------------------------------------------
   ! Compute coefficients for implicit 4th order accurate edge estimation where
   ! edges are asymmetrically related to cell means by:
   ! 
   !    u(j-1/2) + beta*u(j+1/2) = a*u(j) + b*u(j+1) + c*u(j+2)
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: h0, h1, h2
      real(r8), intent(out) :: beta, a, b, c

      real(r8) :: h01, h012, h02, h12, q1, q2, q3, q4, r, s

      h01 = h0 + h1
      h012 = h0 + h1 + h2
      h02 = h0 + h2
      h12 = h1 + h2
      q1 = c1/h01
      q2 = c1/h1
      q3 = c1/h012
      q4 = c1/h12
      r = h0*h0*q3*q4*q4
      s = c2*h1
      beta = h01*h012*q2*q4
      a = ((h01 + c3*h012 + s)*h0 + s*h12)*q1*q3
      b = (h01 + h12)*(s*h012 + h2*h02)*q1*q2*r
      c = - h01*r

   end subroutine edge_coeff_ih4_asym

   subroutine edge_coeff_ih4_sym(h0, h1, alpha, beta, b, c)
   ! ---------------------------------------------------------------------------
   ! Compute coefficients for implicit 4th order accurate edge estimation where
   ! edges are symmetrically related to cell means by:
   ! 
   !    alpha*u(j-3/2) + u(j-1/2) + beta*u(j+1/2) = b*u(j-1) + c*u(j)
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: h0, h1
      real(r8), intent(out) :: alpha, beta, b, c

      real(r8) :: q

      q = c1/(h0 + h1)
      alpha = h1*h1*q*q
      beta  = h0*h0*q*q
      b = c2*alpha*(h1 + c2*h0)*q
      c = c2*beta *(h0 + c2*h1)*q

   end subroutine edge_coeff_ih4_sym

   subroutine reconstruct_plm_no_limiting(rcs)
   ! ---------------------------------------------------------------------------
   ! Carry out a reconstruction with piecewise lines
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs

      real(r8) :: sc
      integer :: ns, j

      ns = rcs%n_src

      sc = c2*(rcs%u_src(2) - rcs%u_src(1)) &
             /(rcs%h_src(2) + rcs%h_src(1))
      rcs%polycoeff(2, 1) = sc*rcs%h_src(1)
      rcs%polycoeff(1, 1) = rcs%u_src(1) - c1_2*rcs%polycoeff(2, 1)
      do j = 2, ns - 1
         sc = (rcs%u_src(j + 1) - rcs%u_src(j - 1))*rcs%hci_src(j)
         rcs%polycoeff(2, j) = sc*rcs%h_src(j)
         rcs%polycoeff(1, j) = rcs%u_src(j) - c1_2*rcs%polycoeff(2, j)
      enddo
      sc = c2*(rcs%u_src(ns) - rcs%u_src(ns - 1)) &
             /(rcs%h_src(ns) + rcs%h_src(ns - 1))
      rcs%polycoeff(2, ns) = sc*rcs%h_src(ns)
      rcs%polycoeff(1, ns) = rcs%u_src(ns) - c1_2*rcs%polycoeff(2, ns)
      rcs%polycoeff(3, :) = c0

   end subroutine reconstruct_plm_no_limiting

   subroutine reconstruct_plm_monotonic(rcs, pc_boundary_cells)
   ! ---------------------------------------------------------------------------
   ! Carry out a reconstruction with piecewise lines and apply limiting to
   ! ensure a monotonic reconstruction.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs
      logical, intent(in) :: pc_boundary_cells

      real(r8) :: sl, sr, sc
      integer :: ns, j

      ns = rcs%n_src

      ! Use monotonized central-difference limiter.

      do j = 2, ns - 1
         sl = c2*(rcs%u_src(j) - rcs%u_src(j - 1))*rcs%hi_src(j)
         sr = c2*(rcs%u_src(j + 1) - rcs%u_src(j))*rcs%hi_src(j)
         if (sl*sr > c0) then
            sc = (rcs%u_src(j + 1) - rcs%u_src(j - 1))*rcs%hci_src(j)
            sc = sign(min(abs(sl), abs(sr), abs(sc)), sc)
         else
            sc = c0
         endif
         rcs%polycoeff(2, j) = sc*rcs%h_src(j)
         rcs%polycoeff(1, j) = rcs%u_src(j) - c1_2*rcs%polycoeff(2, j)
      enddo

      if (pc_boundary_cells) then
         rcs%polycoeff(1, 1) = rcs%u_src(1)
         rcs%polycoeff(2, 1) = c0
         rcs%polycoeff(1, ns) = rcs%u_src(ns)
         rcs%polycoeff(2, ns) = c0
      else
         sc = c2*(rcs%u_src(2) - rcs%u_src(1)) &
                /(rcs%h_src(2) + rcs%h_src(1))
         rcs%polycoeff(2, 1) = sc*rcs%h_src(1)
         rcs%polycoeff(1, 1) = rcs%u_src(1) - c1_2*rcs%polycoeff(2, 1)
         sc = c2*(rcs%u_src(ns) - rcs%u_src(ns - 1)) &
                /(rcs%h_src(ns) + rcs%h_src(ns - 1))
         rcs%polycoeff(2, ns) = sc*rcs%h_src(ns)
         rcs%polycoeff(1, ns) = rcs%u_src(ns) &
                              - c1_2*rcs%polycoeff(2, ns)
      endif
      rcs%polycoeff(3, :) = c0

   end subroutine reconstruct_plm_monotonic

   subroutine reconstruct_ppm_edge_values(rcs)
   ! ---------------------------------------------------------------------------
   ! Reconstruct edge values using an implicit 4th order scheme.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs

      real(r8), dimension(rcs%n_src + 1) :: rhs, gam
      real(r8) :: bei
      integer :: ns, j

      ns = rcs%n_src

      ! Obtain right hand side of tridiagonal system of equations.
      rhs(1) = rcs%al*rcs%u_src(1) &
             + rcs%bl*rcs%u_src(2) &
             + rcs%cl*rcs%u_src(3)
      do j = 2, ns
         rhs(j) = rcs%b(j)*rcs%u_src(j - 1) + rcs%c(j)*rcs%u_src(j)
      enddo
      rhs(ns + 1) = rcs%ar*rcs%u_src(ns) &
                  + rcs%br*rcs%u_src(ns - 1) &
                  + rcs%cr*rcs%u_src(ns - 2)

      ! Solve tridiagonal system of equations to obtain edge values of
      ! source data.
      bei = c1
      rcs%uel(1) = rhs(1)
      do j = 2, ns + 1
         gam(j) = rcs%beta(j - 1)*bei
         bei = c1/(c1 - rcs%alpha(j)*gam(j))
         rcs%uel(j) = (rhs(j) - rcs%alpha(j)*rcs%uel(j - 1))*bei
      enddo
      do j = ns, 1, - 1
         rcs%uel(j) = rcs%uel(j) - gam(j + 1)*rcs%uel(j + 1)
         rcs%uer(j) = rcs%uel(j + 1)
      enddo

   end subroutine reconstruct_ppm_edge_values

   subroutine reconstruct_ppm_interior_monotonic(rcs)
   ! ---------------------------------------------------------------------------
   ! Apply limiting to ensure a monotonic reconstruction of piecewise parabolas
   ! for interior grid cells.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs

      real(r8) :: sl, sr, sc, d, q, r
      integer :: j

      do j = 2, rcs%n_src - 1

         sl = c2*(rcs%u_src(j) - rcs%u_src(j - 1))*rcs%hi_src(j)
         sr = c2*(rcs%u_src(j + 1) - rcs%u_src(j))*rcs%hi_src(j)

         if (sl*sr > c0) then

            sc = (rcs%u_src(j + 1) - rcs%u_src(j - 1))*rcs%hci_src(j)
            sc = sign(min(abs(sl), abs(sr), abs(sc)), sc)
            if ( (rcs%u_src(j - 1) - rcs%uel(j)) &
                *(rcs%u_src(j    ) - rcs%uel(j)) > c0) &
               rcs%uel(j) = rcs%u_src(j) &
                          - sign(min(c1_2*rcs%h_src(j)*abs(sc), &
                                     abs(rcs%uel(j) - rcs%u_src(j))), sc)
            if ( (rcs%u_src(j + 1) - rcs%uer(j)) &
                *(rcs%u_src(j    ) - rcs%uer(j)) > c0) &
               rcs%uer(j) = rcs%u_src(j) &
                          + sign(min(c1_2*rcs%h_src(j)*abs(sc), &
                                     abs(rcs%uer(j) - rcs%u_src(j))), sc)

            d = rcs%uer(j) - rcs%uel(j)
            q = d*(c2*rcs%u_src(j) - rcs%uel(j) - rcs%uer(j))
            r = c1_3*d*d
            if     (  q > r) then
               rcs%uel(j) = c3*rcs%u_src(j) - c2*rcs%uer(j)
            elseif (- r > q) then
               rcs%uer(j) = c3*rcs%u_src(j) - c2*rcs%uel(j)
            endif

         else
            rcs%uel(j) = rcs%u_src(j)
            rcs%uer(j) = rcs%u_src(j)
         endif

      enddo

   end subroutine reconstruct_ppm_interior_monotonic

   subroutine reconstruct_ppm_interior_non_oscillatory(rcs)
   ! ---------------------------------------------------------------------------
   ! Apply limiting to prevent a oscillatory reconstruction of piecewise
   ! parabolas for interior grid cells.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs

      real(r8), dimension(rcs%n_src) :: d2
      real(r8) :: sl, sr, sc, d, q, r
      integer :: j

      ! Obtain values proportional to the second derivative of the unlimited
      ! parabolas. 
      do j = 1, rcs%n_src
         d2(j) = rcs%uel(j) - c2*rcs%u_src(j) + rcs%uer(j)
      enddo

      do j = 2, rcs%n_src - 1

         ! Only apply limiting if the sign of the second 
         ! derivative differs from the sign of second derivatives
         ! of any of the neighbouring parabolas.
         if (d2(j - 1)*d2(j) < c0 .or. d2(j)*d2(j + 1) < c0) then

            sl = c2*(rcs%u_src(j) - rcs%u_src(j - 1))*rcs%hi_src(j)
            sr = c2*(rcs%u_src(j + 1) - rcs%u_src(j))*rcs%hi_src(j)

            if (sl*sr > c0) then

               sc = (rcs%u_src(j + 1) - rcs%u_src(j - 1))*rcs%hci_src(j)
               sc = sign(min(abs(sl), abs(sr), abs(sc)), sc)
               if ( (rcs%u_src(j - 1) - rcs%uel(j)) &
                   *(rcs%u_src(j    ) - rcs%uel(j)) > c0) &
                  rcs%uel(j) = rcs%u_src(j) &
                             - sign(min(c1_2*rcs%h_src(j)*abs(sc), &
                                        abs(rcs%uel(j) - rcs%u_src(j))), sc)
               if ( (rcs%u_src(j + 1) - rcs%uer(j)) &
                   *(rcs%u_src(j    ) - rcs%uer(j)) > c0) &
                  rcs%uer(j) = rcs%u_src(j) &
                             + sign(min(c1_2*rcs%h_src(j)*abs(sc), &
                                        abs(rcs%uer(j) - rcs%u_src(j))), sc)

               d = rcs%uer(j) - rcs%uel(j)
               q = d*(c2*rcs%u_src(j) - rcs%uel(j) - rcs%uer(j))
               r = c1_3*d*d
               if     (  q > r) then
                  rcs%uel(j) = c3*rcs%u_src(j) - c2*rcs%uer(j)
               elseif (- r > q) then
                  rcs%uer(j) = c3*rcs%u_src(j) - c2*rcs%uel(j)
               endif

            else
               rcs%uel(j) = rcs%u_src(j)
               rcs%uer(j) = rcs%u_src(j)
            endif

         endif

      enddo

   end subroutine reconstruct_ppm_interior_non_oscillatory

   subroutine reconstruct_ppm_boundary_limited(rcs, pc_boundary_cells)
   ! ---------------------------------------------------------------------------
   ! Handle piecewise parabola limiting of boundary cells.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs
      logical, intent(in) :: pc_boundary_cells

      real(r8) :: dl, dr, d, q, r
      integer :: ns

      ns = rcs%n_src

      if (pc_boundary_cells) then

         ! Piecewise constant reconstruction of boundary cells.
         rcs%uel(1) = rcs%u_src(1)
         rcs%uer(1) = rcs%u_src(1)
         rcs%uel(ns) = rcs%u_src(ns)
         rcs%uer(ns) = rcs%u_src(ns)

      else

         ! Do not treat boundary cells as local extrema, but ensure
         ! that the piecewise parabolas are monotonic within the
         ! boundary cells.

         if ((rcs%u_src(2) - rcs%uer(1))*(rcs%u_src(1) - rcs%uer(1)) > c0) then
            rcs%uel(1) = rcs%u_src(1)
            rcs%uer(1) = rcs%u_src(1)
         elseif (( (rcs%u_src(3) - rcs%u_src(2)) &
                  *(rcs%h_src(1) + rcs%h_src(2)) &
                 - (rcs%u_src(2) - rcs%u_src(1)) &
                  *(rcs%h_src(2) + rcs%h_src(3))) &
                 *(rcs%uel(1) - c2*rcs%u_src(1) + rcs%uer(1)) < c0) then
            rcs%uel(1) = c1_2*(c3*rcs%u_src(1) - rcs%uer(1))
         else
            dl = rcs%u_src(1) - rcs%uel(1)
            dr = rcs%uer(1) - rcs%u_src(1)
            if (dl*dr < c0) then
               rcs%uel(1) = c1_2*(c3*rcs%u_src(1) - rcs%uer(1))
            else
               d = dl + dr
               q = d*(dl - dr)
               r = c1_3*d*d
               if     (  q > r) then
                  rcs%uel(1) = c3*rcs%u_src(1) - c2*rcs%uer(1)
               elseif (- r > q) then
                  rcs%uer(1) = c3*rcs%u_src(1) - c2*rcs%uel(1)
               endif
            endif
         endif

         if ( (rcs%u_src(ns    ) - rcs%uel(ns)) &
             *(rcs%u_src(ns - 1) - rcs%uel(ns)) > c0) then
            rcs%uel(ns) = rcs%u_src(ns)
            rcs%uer(ns) = rcs%u_src(ns)
         elseif (( (rcs%u_src(ns    ) - rcs%u_src(ns - 1)) &
                  *(rcs%h_src(ns - 2) + rcs%h_src(ns - 1)) &
                 - (rcs%u_src(ns - 1) - rcs%u_src(ns - 2)) &
                  *(rcs%h_src(ns - 1) + rcs%h_src(ns    ))) &
                 *(rcs%uel(ns) - c2*rcs%u_src(ns) + rcs%uer(ns)) < c0) then
            rcs%uer(ns) = c1_2*(c3*rcs%u_src(ns) - rcs%uel(ns))
         else
            dl = rcs%u_src(ns) - rcs%uel(ns)
            dr = rcs%uer(ns) - rcs%u_src(ns)
            if (dl*dr < c0) then
               rcs%uer(ns) = c1_2*(c3*rcs%u_src(ns) - rcs%uel(ns))
            else
               d = dl + dr
               q = d*(dl - dr)
               r = c1_3*d*d
               if     (  q > r) then
                  rcs%uel(ns) = c3*rcs%u_src(ns) - c2*rcs%uer(ns)
               elseif (- r > q) then
                  rcs%uer(ns) = c3*rcs%u_src(ns) - c2*rcs%uel(ns)
               endif
            endif
         endif

      endif

   end subroutine reconstruct_ppm_boundary_limited

   subroutine reconstruct_ppm_posdef(rcs)
   ! ---------------------------------------------------------------------------
   ! Modify piecewise parabolas so they are never negative within the grid cell.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs

      real(r8) :: min_u_0, sl, a2, sr, q
      integer :: j

      do j = 1, rcs%n_src
         min_u_0 = min(rcs%u_src(j), c0)
         rcs%uel(j) = max(rcs%uel(j), min_u_0)
         rcs%uer(j) = max(rcs%uer(j), min_u_0)
         sl = c2*(c3*rcs%u_src(j) - c2*rcs%uel(j) - rcs%uer(j))
         a2 = c3*(rcs%uel(j) - c2*rcs%u_src(j) + rcs%uer(j))
         sr = sl + c2*a2
         if (sl < c0 .and. sr > c0) then
            if (a2*rcs%uel(j) - c1_4*sl*sl < a2*min_u_0) then
               q = c3*rcs%u_src(j)/(c3*sl*sr + c4*a2*a2)
               rcs%uel(j) = sl*sl*q
               rcs%uer(j) = sr*sr*q
            endif
         endif
      enddo

   end subroutine reconstruct_ppm_posdef

   subroutine reconstruct_ppm_polycoeff(rcs)
   ! ---------------------------------------------------------------------------
   ! Obtain coefficients for piecewise parabolas from grid cell means and left
   ! and right edge values.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs

      integer :: j

      do j = 1, rcs%n_src
         rcs%polycoeff(1, j) = rcs%uel(j)
         rcs%polycoeff(2, j) = c2*(c3*rcs%u_src(j) - c2*rcs%uel(j) - rcs%uer(j))
         rcs%polycoeff(3, j) = c3*(rcs%uel(j) - c2*rcs%u_src(j) + rcs%uer(j))
      enddo

   end subroutine reconstruct_ppm_polycoeff

   pure function parabola_intersection(pc, u, u_eps, xil, xir) result(xi)

      real(r8), dimension(3), intent(in) :: pc
      real(r8), intent(in) :: u, u_eps, xil, xir

      real(r8) :: xi

      real(r8) :: q, s, xi1, xi2, xim

      if (abs(pc(3)) < u_eps) then
         if (abs(pc(2)) < u_eps) then
            xi = xil
         else
            xi = (u - pc(1))/pc(2)
         endif
      else
         q = c1_2/pc(3)
         s = sqrt(max(c0, pc(2)*pc(2) - c4*pc(3)*(pc(1) - u)))
         xi1 = - (pc(2) + s)*q
         xi2 = - (pc(2) - s)*q
         xim = c1_2*(xil + xir)
         if (abs(xi1 - xim) < abs(xi2 - xim)) then
            xi = xi1
         else
            xi = xi2
         endif
      endif

      xi = max(xil, min(xir, xi))

   end function parabola_intersection

   ! ---------------------------------------------------------------------------
   ! Public procedures.
   ! ---------------------------------------------------------------------------

   function prepare_reconstruction(x_edge_src, method, rcs) result(errstat)
   ! ---------------------------------------------------------------------------
   ! Prepare reconstruction based on edge locations of source grid cells and
   ! requested reconstruction method. Reconstruction data is stored in a
   ! reconstruction data structure.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:), intent(in) :: x_edge_src
      integer, intent(in) :: method
      type(reconstruction_struct), intent(inout) :: rcs

      integer :: errstat

      integer :: n_src, i, j

      errstat = hor3map_noerr

      ! Check reconstruction method.
      if (method /= hor3map_pcm .and. method /= hor3map_plm .and. &
          method /= hor3map_ppm) then
         errstat = hor3map_invalid_method
         return
      endif

      ! Number of source grid cells.
      rcs%n_src_all = size(x_edge_src) - 1

      ! Check that source grid edges are monotonically increasing or decreasing.
      if (x_edge_src(rcs%n_src_all + 1) - x_edge_src(1) > c0) then
         do j = 1, rcs%n_src_all
            if (x_edge_src(j + 1) < x_edge_src(j)) then
               errstat = hor3map_nonmonotonic_src_edges
               return
            endif
         enddo
      else
         do j = 1, rcs%n_src_all
            if (x_edge_src(j + 1) > x_edge_src(j)) then
               errstat = hor3map_nonmonotonic_src_edges
               return
            endif
         enddo
      endif

      ! Set small value with same dimensions as edge locations.
      rcs%x_eps = max(abs( x_edge_src(rcs%n_src_all + 1) &
                         - x_edge_src(1)), eps)*eps

      ! Number of source grid cells excluding empty grid cells.
      n_src = 0
      do j = 1, rcs%n_src_all
         if (abs(x_edge_src(j + 1) - x_edge_src(j)) > rcs%x_eps) &
            n_src = n_src + 1
      enddo

      ! If no non-empty grid cells are present, return with error code.
      if (n_src == 0) then
         errstat = hor3map_src_extent_too_small
         return
      endif

      ! Select actual reconstruction method to be used, considering requested
      ! method and availability of non-empty grid cells.
      if     (method == hor3map_pcm .or. n_src == 1) then
         rcs%method = hor3map_pcm
      elseif (method == hor3map_plm .or. n_src < 4) then
         rcs%method = hor3map_plm
      else
         rcs%method = method
      endif

      ! If needed, allocate arrays in reconstruction data structure.
      if (.not. rcs%prepared) then
         rcs%n_src = n_src
         errstat = allocate_rcs(rcs)
         if (errstat /= hor3map_noerr) return
      elseif (rcs%n_src /= n_src) then
         call free_rcs(rcs)
         rcs%n_src = n_src
         errstat = allocate_rcs(rcs)
         if (errstat /= hor3map_noerr) return
      endif

      ! Set index map from an array with only non-empty source grid cells to a
      ! full source array.
      j = 0
      do i = 1, rcs%n_src_all
         if (abs(x_edge_src(i + 1) - x_edge_src(i)) > rcs%x_eps) then
            j = j + 1
            rcs%src_index_map(j) = i
         endif
      enddo
      ! Store source edge locations in reconstruction data structure.
      rcs%x_edge_src(1) = x_edge_src(1)
      do j = 1, rcs%n_src
         rcs%x_edge_src(j + 1) = x_edge_src(rcs%src_index_map(j) + 1)
      enddo

      ! From edge locations, obtain source grid cell widths and their
      ! multiplicative inverse.
      do j = 1, rcs%n_src
         rcs%h_src(j) = abs(rcs%x_edge_src(j + 1) - rcs%x_edge_src(j))
         rcs%hi_src(j) = c1/rcs%h_src(j)
      enddo
      
      ! Store the multiplicative inverse of cell width used for estimating
      ! centered linear slope.
      do j = 2, rcs%n_src - 1
         rcs%hci_src(j) = c2/( rcs%h_src(j - 1) + c2*rcs%h_src(j) &
                             + rcs%h_src(j + 1))
      enddo

      if (rcs%method == hor3map_ppm) call edge_coeff_ih4(rcs)

      rcs%prepared = .true.

   end function prepare_reconstruction

   function prepare_remapping(rcs, x_edge_dst, rms) result(errstat)
   ! ---------------------------------------------------------------------------
   ! Prepare remapping based on a reconstruction data structure and edge
   ! locations of destination grid cells. Remapping data is stored in a remap
   ! data structure.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(in) :: rcs
      real(r8), dimension(:), intent(in) :: x_edge_dst
      type(remap_struct), intent(out) :: rms

      integer :: errstat

      real(r8) :: xil
      integer :: n_dst, j, js, jd, iseg

      errstat = hor3map_noerr

      ! Check that the reconstruction has been prepared.
      if (.not. rcs%prepared) then
         errstat = hor3map_recon_not_prepared
         return
      endif

      ! Number of destination grid cells.
      n_dst = size(x_edge_dst) - 1

      ! Check for consistency between the source and destination grid range.
      if (abs( rcs%x_edge_src(1) - x_edge_dst(1)) > rcs%x_eps .or. &
          abs( rcs%x_edge_src(rcs%n_src + 1) &
             - x_edge_dst(n_dst + 1)) > rcs%x_eps) then
         errstat = hor3map_inconsistent_grid_range
         return
      endif

      ! Check that destination grid edges are monotonically increasing or
      ! decreasing.
      if (x_edge_dst(n_dst + 1) - x_edge_dst(1) > c0) then
         do j = 1, n_dst
            if (x_edge_dst(j + 1) < x_edge_dst(j)) then
               errstat = hor3map_nonmonotonic_dst_edges
               return
            endif
         enddo
      else
         do j = 1, n_dst
            if (x_edge_dst(j + 1) > x_edge_dst(j)) then
               errstat = hor3map_nonmonotonic_dst_edges
               return
            endif
         enddo
      endif

      ! If needed, allocate arrays in remap data structure.
      if (.not. rms%prepared) then
         rms%n_dst = n_dst
         errstat = allocate_rms(rcs, rms)
         if (errstat /= hor3map_noerr) return
      elseif (rms%n_dst /= n_dst) then
         call free_rms(rms)
         rms%n_dst = n_dst
         errstat = allocate_rms(rcs, rms)
         if (errstat /= hor3map_noerr) return
      endif

      ! From edge locations, obtain destination grid cell widths and their
      ! multiplicative inverse.
      do j = 1, rms%n_dst
         rms%h_dst(j) = abs(x_edge_dst(j + 1) - x_edge_dst(j))
         if (rms%h_dst(j) > rcs%x_eps) then
            rms%hi_dst(j) = c1/rms%h_dst(j)
         else
            rms%hi_dst(j) = c0
         endif
      enddo

      ! Locate all segments that require integration of the reconstructed source
      ! data and obtain integration limits and destination index for the
      ! accumulation of integrals.

      js = 1
      jd = 1
      iseg = 0
      rms%n_src_seg(js) = 0
      xil = c0

      if (x_edge_dst(n_dst + 1) - x_edge_dst(1) > c0) then

         do
            iseg = iseg + 1
            rms%n_src_seg(js) = rms%n_src_seg(js) + 1
            rms%seg_dst_index(iseg) = jd
            if     (  abs(rcs%x_edge_src(js + 1) - x_edge_dst(jd + 1)) &
                   <= rcs%x_eps) then
               if (rms%hi_dst(jd) == c0) then
                  rms%seg_int_lim(iseg) = xil
               else
                  rms%seg_int_lim(iseg) = c1
               endif
               if (js == rcs%n_src) exit
               xil = c0
               js = js + 1
               jd = jd + 1
               rms%n_src_seg(js) = 0
            elseif (rcs%x_edge_src(js + 1) < x_edge_dst(jd + 1)) then
               rms%seg_int_lim(iseg) = c1
               xil = c0
               js = js + 1
               rms%n_src_seg(js) = 0
            else
               if (rms%hi_dst(jd) == c0) then
                  rms%seg_int_lim(iseg) = xil
               else
                  rms%seg_int_lim(iseg) = &
                     (x_edge_dst(jd + 1) - rcs%x_edge_src(js))*rcs%hi_src(js)
                  xil = rms%seg_int_lim(iseg)
               endif
               jd = jd + 1
            endif
         enddo

      else

         do
            iseg = iseg + 1
            rms%n_src_seg(js) = rms%n_src_seg(js) + 1
            rms%seg_dst_index(iseg) = jd
            if     (  abs(rcs%x_edge_src(js + 1) - x_edge_dst(jd + 1)) &
                   <= rcs%x_eps) then
               if (rms%hi_dst(jd) == c0) then
                  rms%seg_int_lim(iseg) = xil
               else
                  rms%seg_int_lim(iseg) = c1
               endif
               if (js == rcs%n_src) exit
               xil = c0
               js = js + 1
               jd = jd + 1
               rms%n_src_seg(js) = 0
            elseif (rcs%x_edge_src(js + 1) > x_edge_dst(jd + 1)) then
               rms%seg_int_lim(iseg) = c1
               xil = c0
               js = js + 1
               rms%n_src_seg(js) = 0
            else
               if (rms%hi_dst(jd) == c0) then
                  rms%seg_int_lim(iseg) = xil
               else
                  rms%seg_int_lim(iseg) = &
                     (rcs%x_edge_src(js) - x_edge_dst(jd + 1))*rcs%hi_src(js)
                  xil = rms%seg_int_lim(iseg)
               endif
               jd = jd + 1
            endif
         enddo

      endif

      rms%prepared = .true.

   end function prepare_remapping

   function reconstruct(rcs, u_src, limiting, pc_boundary_cells) result(errstat)
   ! ---------------------------------------------------------------------------
   ! Carry out the piecewise polynomial reconstruction of the source data with
   ! desired limiting method and handling of boundaries.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs
      real(r8), dimension(:), intent(in) :: u_src
      integer, intent(in) :: limiting
      logical, intent(in) :: pc_boundary_cells

      integer :: errstat

      integer :: j

      errstat = hor3map_noerr

      ! Check that the reconstruction has been prepared.
      if (.not. rcs%prepared) then
         errstat = hor3map_recon_not_prepared
         return
      endif

      if (size(u_src) /= rcs%n_src_all) then
         errstat = hor3map_src_size_mismatch
         return
      endif

      rcs%limiting = limiting

      ! Copy source data array to array that excludes empty grid cells.
      do j = 1, rcs%n_src
        rcs%u_src(j) = u_src(rcs%src_index_map(j))
      enddo

      select case (rcs%method)
         case (hor3map_plm)
            select case (limiting)
               case (hor3map_no_limiting)
                  call reconstruct_plm_no_limiting(rcs)
               case (hor3map_monotonic)
                  call reconstruct_plm_monotonic(rcs, pc_boundary_cells)
               case default
                  errstat = hor3map_invalid_plm_limiting
                  return
            end select
         case (hor3map_ppm)
            call reconstruct_ppm_edge_values(rcs)
            select case (limiting)
               case (hor3map_no_limiting)
               case (hor3map_monotonic)
                  call reconstruct_ppm_interior_monotonic(rcs)
                  call reconstruct_ppm_boundary_limited(rcs, pc_boundary_cells)
               case (hor3map_non_oscillatory)
                  call reconstruct_ppm_interior_non_oscillatory(rcs)
                  call reconstruct_ppm_boundary_limited(rcs, pc_boundary_cells)
               case (hor3map_non_oscillatory_posdef)
                  call reconstruct_ppm_interior_non_oscillatory(rcs)
                  call reconstruct_ppm_boundary_limited(rcs, pc_boundary_cells)
                  call reconstruct_ppm_posdef(rcs)
               case default
                  errstat = hor3map_invalid_ppm_limiting
                  return
            end select
            call reconstruct_ppm_polycoeff(rcs)
      end select

      rcs%reconstructed = .true.

   end function reconstruct

   function regrid(rcs, u_edge_grd, x_edge_grd, missing_value) result(errstat)
   ! ---------------------------------------------------------------------------
   ! Find grid locations where desired grid cell edge data values intersect with
   ! a reconstruction of the source data.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(in) :: rcs
      real(r8), dimension(:), intent(in) :: u_edge_grd
      real(r8), dimension(:), intent(out) :: x_edge_grd
      real(r8), intent(in) :: missing_value

      integer :: errstat

      real(r8), dimension(3) :: pcl, pcr
      real(r8) :: u_min, u_max, u_eps, u_sgn, umr, uml, xi, &
                  duml, dumr, uerl, uelr
      integer :: n_grd, jg, js

      errstat = hor3map_noerr

      ! Check that the reconstruction is available.
      if (.not. rcs%reconstructed) then
         errstat = hor3map_recon_not_available
         return
      endif

      ! Number of grid edges.
      n_grd = size(u_edge_grd)

      ! Check grid array size consistency.
      if (size(x_edge_grd) /= n_grd) then
         errstat = hor3map_grd_size_mismatch
         return
      endif

      ! Initialize grid intersections as missing value.
      x_edge_grd(:) = missing_value

      ! Return in case PCM method is used.
      if (rcs%method == hor3map_pcm) return

      ! Set small value with same dimensions as source data.
      u_min = minval(rcs%u_src(:))
      u_max = maxval(rcs%u_src(:))
      if (abs(u_max - u_min) < eps) then
         return
      endif
      u_eps = abs(u_max - u_min)*eps

      ! To indicate monotonically increasing or decreasing source values, use
      ! the sign of the difference of the source boundary values.
      u_sgn = sign(c1, rcs%u_src(rcs%n_src) - rcs%u_src(1))

      ! Find possible intersections in the first half of the first source grid
      ! cell. 
      jg = 1
      do
         if ((u_edge_grd(jg) - rcs%uel(1))*u_sgn >= c0) exit
         jg = jg + 1
         if (jg > n_grd) return
      enddo
      js = 1
      umr =      rcs%polycoeff(1, js) &
          + c1_2*rcs%polycoeff(2, js) &
          + c1_4*rcs%polycoeff(3, js)
      do
         if ((u_edge_grd(jg) - umr)*u_sgn > c0) exit
         xi = parabola_intersection(rcs%polycoeff(:, js), u_edge_grd(jg), &
                                    u_eps, c0, c1_2)
         x_edge_grd(jg) = rcs%x_edge_src(js) &
                        + (rcs%x_edge_src(js + 1) - rcs%x_edge_src(js))*xi
         jg = jg + 1
         if (jg > n_grd) return
      enddo

      outer: do 

         ! For the current grid edge index, find the index of the first source
         ! grid cell with mid point reconstructed value larger than the grid
         ! edge value.
         do
            uml = umr
            umr =      rcs%polycoeff(1, js) &
                + c1_2*rcs%polycoeff(2, js) &
                + c1_4*rcs%polycoeff(3, js)
            if ((u_edge_grd(jg) - umr)*u_sgn <= c0) exit
            js = js + 1
            if (js > rcs%n_src) exit outer
         enddo

         ! Construct new parabolas left and right of the edge that are
         ! continuous and smooth across the edge and with the original piecewise
         ! parabolas left and right of the edge at the mid points of their
         ! respective grid cells.
         duml = rcs%polycoeff(2, js - 1) + rcs%polycoeff(3, js - 1)
         dumr = rcs%polycoeff(2, js    ) + rcs%polycoeff(3, js    )
         pcr(2) = (c4*(umr - uml) - duml - dumr) &
                  *rcs%h_src(js)/(rcs%h_src(js - 1) + rcs%h_src(js))
         pcr(1) = umr - c1_4*(dumr + pcr(2))
         if (pcr(2)*(rcs%u_src(js) - rcs%u_src(js - 1)) < c0) then
            ! If the slope of the new parabolas are non-monotonic at the
            ! edge, set the edge slope to zero and enforce that the new
            ! parabolas cross the edge within the interval spanned by the
            ! edge values of the original piecewise parabolas. Smoothness
            ! with the original piecewise parabolas at grid cell mid points
            ! is then not guaranteed.
            pcr(2) = c0
            uerl = rcs%uer(js - 1)
            uelr = rcs%uel(js)
            pcr(1) = min(max(pcr(1), min(uerl, uelr)), max(uerl, uelr))
            pcr(3) = c4*(umr - pcr(1))
            pcl(1) = c4*uml - c3*pcr(1)
            pcl(2) = c2*(pcr(1) - pcl(1))
            pcl(3) = - c1_2*pcl(2)
         else
            pcr(3) = dumr - pcr(2)
            pcl(1) = pcr(1) - duml
            pcl(2) = c4*(uml - pcl(1)) - duml
            pcl(3) = duml - pcl(2)
         endif

         ! Find all intersections with piecewise parabola in the last half of
         ! the source grid cell left of the edge.
         do
            if ((u_edge_grd(jg) - pcr(1))*u_sgn > c0) exit
            xi = parabola_intersection(pcl, u_edge_grd(jg), u_eps, c1_2, c1)
            x_edge_grd(jg) = rcs%x_edge_src(js - 1) &
                           + (rcs%x_edge_src(js) - rcs%x_edge_src(js - 1))*xi
            jg = jg + 1
            if (jg > n_grd) return
         enddo

         ! Find all intersections with piecewise parabola in the first half of
         ! the source grid cell right of the edge.
         do
            if ((u_edge_grd(jg) - umr)*u_sgn > c0) exit
            xi = parabola_intersection(pcr, u_edge_grd(jg), u_eps, c0, c1_2)
            x_edge_grd(jg) = rcs%x_edge_src(js) &
                           + (rcs%x_edge_src(js + 1) - rcs%x_edge_src(js))*xi
            jg = jg + 1
            if (jg > n_grd) return
         enddo

      enddo outer

      ! Find possible intersections in the last half of the last source grid
      ! cell. 
      js = rcs%n_src
      do
         if ((u_edge_grd(jg) - rcs%uer(js))*u_sgn > c0) return
         xi = parabola_intersection(rcs%polycoeff(:, js), u_edge_grd(jg), &
                                    u_eps, c1_2, c1)
         x_edge_grd(jg) = rcs%x_edge_src(js) &
                        + (rcs%x_edge_src(js + 1) - rcs%x_edge_src(js))*xi
         jg = jg + 1
         if (jg > n_grd) return
      enddo

   end function regrid

   function remap(rcs, rms, u_dst) result(errstat)
   ! ---------------------------------------------------------------------------
   ! Carry out the remapping of a piecewise polynomial reconstruction of the
   ! source data to a destination grid.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(in) :: rcs
      type(remap_struct), intent(in) :: rms
      real(r8), dimension(:), intent(out) :: u_dst

      integer :: errstat

      real(r8) :: xil, xir, adl, adr
      integer :: iseg, js, jd, i_src_seg

      errstat = hor3map_noerr

      ! Check that the remapping has been prepared
      if (.not. rms%prepared) then
         errstat = hor3map_remap_not_prepared
         return
      endif

      ! Check that the reconstruction is available.
      if (.not. rcs%reconstructed) then
         errstat = hor3map_recon_not_available
         return
      endif

      if (size(u_dst) /= rms%n_dst) then
         errstat = hor3map_dst_size_mismatch
         return
      endif

      u_dst(:) = 0._r8
      iseg = 0

      select case (rcs%method)

         case (hor3map_pcm)

            do js = 1, rcs%n_src
               if (rms%n_src_seg(js) == 1) then
                  ! No integration needed
                  iseg = iseg + 1
                  jd = rms%seg_dst_index(iseg)
                  u_dst(jd) = u_dst(jd) &
                            + rcs%u_src(js)*rcs%h_src(js)*rms%hi_dst(jd)
               else
                  ! Integrate the required segments of each source grid cell in
                  ! succession, adding the integrals to the appropriate
                  ! destination grid cells.
                  xil = c0
                  do i_src_seg = 1, rms%n_src_seg(js)
                     iseg = iseg + 1
                     xir = rms%seg_int_lim(iseg)
                     jd = rms%seg_dst_index(iseg)
                     if (xil == xir) then
                        u_dst(jd) = rcs%u_src(js)
                     else
                        u_dst(jd) = u_dst(jd) &
                                  + rcs%u_src(js)*(xir - xil)*rcs%h_src(js) &
                                    *rms%hi_dst(jd)
                        xil = xir
                     endif
                  enddo
               endif
            enddo

         case (hor3map_plm)

            do js = 1, rcs%n_src
               if (rms%n_src_seg(js) == 1) then
                  ! No integration needed
                  iseg = iseg + 1
                  jd = rms%seg_dst_index(iseg)
                  u_dst(jd) = u_dst(jd) &
                            + rcs%u_src(js)*rcs%h_src(js)*rms%hi_dst(jd)
               else
                  ! Integrate the required segments of each source grid cell in
                  ! succession, adding the integrals to the appropriate
                  ! destination grid cells.
                  xil = c0
                  adl = c0
                  do i_src_seg = 1, rms%n_src_seg(js)
                     iseg = iseg + 1
                     xir = rms%seg_int_lim(iseg)
                     jd = rms%seg_dst_index(iseg)
                     if (xil == xir) then
                        u_dst(jd) = rcs%polycoeff(1, js) &
                                  + rcs%polycoeff(2, js)*xir
                     else
                        adr = (      rcs%polycoeff(1, js) &
                              + c1_2*rcs%polycoeff(2, js)*xir)*xir
                        u_dst(jd) = u_dst(jd) &
                                  + (adr - adl)*rcs%h_src(js)*rms%hi_dst(jd)
                        xil = xir
                        adl = adr
                     endif
                  enddo
               endif
            enddo

         case (hor3map_ppm)

            do js = 1, rcs%n_src
               if (rms%n_src_seg(js) == 1) then
                  ! No integration needed
                  iseg = iseg + 1
                  jd = rms%seg_dst_index(iseg)
                  u_dst(jd) = u_dst(jd) &
                            + rcs%u_src(js)*rcs%h_src(js)*rms%hi_dst(jd)
               else
                  ! Integrate the required segments of each source grid cell in
                  ! succession, adding the integrals to the appropriate
                  ! destination grid cells.
                  xil = c0
                  adl = c0
                  do i_src_seg = 1, rms%n_src_seg(js)
                     iseg = iseg + 1
                     xir = rms%seg_int_lim(iseg)
                     jd = rms%seg_dst_index(iseg)
                     if (xil == xir) then
                        u_dst(jd) =   rcs%polycoeff(1, js) &
                                  + ( rcs%polycoeff(2, js) &
                                    + rcs%polycoeff(3, js)*xir)*xir
                     else
                        adr = (        rcs%polycoeff(1, js) &
                              + ( c1_2*rcs%polycoeff(2, js) &
                                + c1_3*rcs%polycoeff(3, js)*xir)*xir)*xir
                        u_dst(jd) = u_dst(jd) &
                                  + (adr - adl)*rcs%h_src(js)*rms%hi_dst(jd)
                        xil = xir
                        adl = adr
                     endif
                  enddo
               endif
            enddo

      end select

   end function remap

   subroutine free_rcs(rcs)
   ! ---------------------------------------------------------------------------
   ! Deallocate arrays and reset flags.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs

      deallocate(rcs%x_edge_src, rcs%h_src, rcs%hi_src, rcs%hci_src, &
                 rcs%alpha, rcs%beta, rcs%b, rcs%c, rcs%u_src, rcs%uel, &
                 rcs%uer, rcs%src_index_map, rcs%polycoeff)

      rcs%prepared = .false.
      rcs%reconstructed = .false.

   end subroutine free_rcs

   subroutine free_rms(rms)
   ! ---------------------------------------------------------------------------
   ! Deallocate arrays and reset flags.
   ! ---------------------------------------------------------------------------

      type(remap_struct), intent(inout) :: rms

      deallocate(rms%h_dst, rms%hi_dst, rms%seg_int_lim, rms%n_src_seg, &
                 rms%seg_dst_index)

      rms%prepared = .false.

   end subroutine free_rms

   pure function hor3map_errstr(errstat) result(errstr)
   ! ---------------------------------------------------------------------------
   ! Returns static reference to an error message string corresponding to a
   ! HOR3MAP error status.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: errstat ! Error status.

      character(len = 80) :: errstr  ! Error message string.

      if (errstat > 0 .and. errstat <= hor3map_errmsg_num) then
         errstr = errmsg(errstat)
      else
         errstr = 'Unknown error status!'
      endif

   end function hor3map_errstr

end module mod_hor3map
