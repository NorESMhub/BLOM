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

#undef DEBUG
#define DEBUG

   use, intrinsic :: iso_fortran_env, only: real64
#ifdef DEBUG
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_signaling_nan
#endif

   implicit none

   private

   ! Option parameters.
   integer, parameter :: &
      hor3map_pcm                    = 100, & ! Reconstruction methods
      hor3map_plm                    = 101, &
      hor3map_ppm                    = 102, &
      hor3map_pqm                    = 103, &
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
      hor3map_invalid_pqm_limiting    = 12, &
      hor3map_recon_not_available     = 13, &
      hor3map_grd_size_mismatch       = 14, &
      hor3map_remap_not_prepared      = 15, &
      hor3map_dst_size_mismatch       = 16, &
      hor3map_errmsg_num = 16
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
       "Invalid limiting method for PQM!                                 ", &
       "Call 'reconstruct' first!                                        ", &
       "Size mismatch between grid edge values and locations!            ", &
       "Call 'prepare_remapping' first!                                  ", &
       "Size mismatch between destination grid edges and data array!     "]

   ! Numeric data types.
   integer, parameter :: &
      r8 = real64

   ! Maximum order of accuracy in edge estimation.
   integer, parameter :: maxord = 6

   ! Small non-dimensional value.
   real(r8), parameter :: eps = 1.e-14_r8

   ! Lower bounds of prod(h(i))/max(h(i))^n, where h(i) are cell widths
   ! belonging to a stencil of n grid cells. Respecting these bounds will ensure
   ! condition numbers below 10^8 of matrices involved in various linear
   ! equation systems.
   real(r8), parameter :: &
      hplim_ih4 = 5.e-7_r8, &
      hplim_ih6 = 5.e-5_r8, &
      hplim_eh4 = 3.e-10_r8, &
      hplim_eh6 = 3.e-7_r8

   ! Numeric constants.
   real(r8), parameter :: &
      c0 = 0._r8, c1 = 1._r8, c2 = 2._r8, c3 = 3._r8, c4 = 4._r8, c5 = 5._r8, &
      c6 = 6._r8, c12 = 12._r8, c15 = 15._r8, c18 = 18._r8, c28 = 28._r8, &
      c30 = 30._r8, c32 = 32._r8, c60 = 60._r8, &
      c1_2 = 1._r8/2._r8, c1_3 = 1._r8/3._r8, c1_4 = 1._r8/4._r8, &
      c1_5 = 1._r8/5._r8, c1_6 = 1._r8/6._r8, c1_12 = 1._r8/12._r8, &
      c1_24 = 1._r8/24._r8, c1_80 = 1._r8/80._r8, c1_120 = 1._r8/120._r8, &
      c3_4 = 3._r8/4._r8, c3_2 = 3._r8/2._r8, c5_2 = 5._r8/2._r8, &
      c9_2 = 9._r8/2._r8

   type reconstruction_struct
      real(r8), allocatable, dimension(:, :) :: &
         tdecoeff, tdscoeff, lblu, rblu, polycoeff
      real(r8), allocatable, dimension(:) :: &
         x_edge_src, h_src, hi_src, hci_src, src_dst_weight, &
         u_src, uel, uer, usl, usr
      real(r8) :: x_eps
      integer, allocatable, dimension(:) :: src_dst_index
      integer :: n_src_all, n_src, method, limiting
      logical :: &
         alloced       = .false., &
         prepared      = .false., &
         reconstructed = .false.
   end type reconstruction_struct

   type remap_struct
      real(r8), allocatable, dimension(:) :: h_dst, hi_dst, seg_int_lim
      integer, allocatable, dimension(:) :: n_src_seg, seg_dst_index
      integer :: n_dst
      logical :: &
         alloced  = .false., &
         prepared = .false.
   end type remap_struct

   public :: reconstruction_struct, remap_struct, &
             prepare_reconstruction, prepare_remapping, &
             reconstruct, regrid, remap, free_rcs, free_rms, &
             hor3map_pcm, hor3map_plm, hor3map_ppm, hor3map_pqm, &
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

      allocate(rcs%x_edge_src(rcs%n_src_all + 1), rcs%h_src(rcs%n_src_all), &
               rcs%hi_src(rcs%n_src_all), rcs%hci_src(2:rcs%n_src_all - 1), &
               rcs%src_dst_weight(rcs%n_src_all), &
               rcs%tdecoeff(maxord, rcs%n_src_all), &
               rcs%tdscoeff(maxord, rcs%n_src_all), &
               rcs%lblu(maxord, maxord), rcs%rblu(maxord, maxord), &
               rcs%u_src(rcs%n_src_all), &
               rcs%uel(rcs%n_src_all), rcs%uer(rcs%n_src_all), &
               rcs%usl(rcs%n_src_all), rcs%usr(rcs%n_src_all), &
               rcs%src_dst_index(rcs%n_src_all), &
               rcs%polycoeff(maxord - 1, rcs%n_src_all), &
               stat = allocstat)

      if (allocstat == 0) then
         errstat = hor3map_noerr
      else
         errstat = hor3map_failed_to_allocate_rcs
         return
      endif

      rcs%alloced = .true.

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
               rms%seg_int_lim(rcs%n_src_all + rms%n_dst), &
               rms%n_src_seg(rcs%n_src_all), &
               rms%seg_dst_index(rcs%n_src_all + rms%n_dst), &
               stat = allocstat)

      if (allocstat == 0) then
         errstat = hor3map_noerr
      else
         errstat = hor3map_failed_to_allocate_rms
         return
      endif

      rms%alloced = .true.

   end function allocate_rms

   subroutine lu_decompose(n, a)
   ! ---------------------------------------------------------------------------
   ! Replace the n x n input matrix A with its LU decomposition.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: n
      real(r8), dimension(:, :), intent(inout) :: a

      integer :: i, j, k

      do k = 1, n - 1
         do i = k + 1, n
            a(i, k) = a(i, k)/a(k, k)
            do j = k + 1, n
               a(i, j) = a(i, j) - a(i, k)*a(k, j)
            enddo
         enddo
      enddo

   end subroutine lu_decompose

   subroutine lu_solve(n, lu, x)
   ! ---------------------------------------------------------------------------
   ! Solve the linear system of equations A*x = b using the LU decomposition of
   ! the n x n matrix A. The argument x has b as input and is replaced with the
   ! solution upon return.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: n
      real(r8), dimension(:, :), intent(in) :: lu
      real(r8), dimension(:), intent(inout) :: x

      integer :: i, j

      ! Forward substitution.
      do i = 2, n
         do j = 1, i - 1
            x(i) = x(i) - lu(i, j)*x(j)
         enddo
      enddo

      ! Back substitution.
      x(n) = x(n)/lu(n, n)
      do i = n - 1, 1, -1
         do j = i + 1, n
            x(i) = x(i) - lu(i, j)*x(j)
         enddo
         x(i) = x(i)/lu(i, i)
      enddo

   end subroutine lu_solve

   subroutine edge_ih4_coeff(h, tdecoeff)
   ! ---------------------------------------------------------------------------
   ! Compute row coefficients for the tridiagonal system of equations to be
   ! solved for 4th order accurate edge estimates.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:), intent(in) :: h
      real(r8), dimension(:), intent(inout) :: tdecoeff

      real(r8) :: q

      q = c1/(h(1) + h(2))
      tdecoeff(1) = h(2)*h(2)*q*q
      tdecoeff(2) = h(1)*h(1)*q*q
      tdecoeff(3) = c2*tdecoeff(1)*(h(2) + c2*h(1))*q
      tdecoeff(4) = c2*tdecoeff(2)*(h(1) + c2*h(2))*q

   end subroutine edge_ih4_coeff

   subroutine edge_ih6_slope_ih5_coeff_common(a, tdecoeff, tdscoeff)
   ! ---------------------------------------------------------------------------
   ! Common procedure for the various stencils for the computation of row
   ! coefficients for the tridiagonal system of equations to be solved for 6th
   ! and 5th order accurate edge and slope estimates, respectively.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:, :), intent(inout) :: a
      real(r8), dimension(:), intent(inout) :: tdecoeff, tdscoeff

      real(r8), dimension(6, 6) :: b

      ! Define matrix for linear system to be solved for slope coefficients.

      b(1:5, 3:6) = a(2:6, 3:6)

      b(1, 1) = c1
      b(2, 1) = c2*a(2, 1)
      b(3, 1) = c3*a(3, 1)
      b(4, 1) = c4*a(4, 1)
      b(5, 1) = c5*a(5, 1)
      b(6, 1) = c0

      b(1, 2) = c1
      b(2, 2) = c2*a(2, 2)
      b(3, 2) = c3*a(3, 2)
      b(4, 2) = c4*a(4, 2)
      b(5, 2) = c5*a(5, 2)
      b(6, 2) = c0

      b(6, 3:6) = c1

      ! Solve linear system for edge coefficients.
      tdecoeff(:) = [ - c1, c0, c0, c0, c0, c0]
      call lu_decompose(6, a)
      call lu_solve(6, a, tdecoeff)

      ! Solve linear system for slope coefficients.
      tdscoeff(:) = [ - c1, c0, c0, c0, c0, c0]
      call lu_decompose(6, b)
      call lu_solve(6, b, tdscoeff)

   end subroutine edge_ih6_slope_ih5_coeff_common

   subroutine edge_ih6_slope_ih5_coeff_asymleft(h, tdecoeff, tdscoeff)
   ! ---------------------------------------------------------------------------
   ! With an asymmetrical stencil, where edge values are shifted left compared
   ! to cell mean values, compute row coefficients for the tridiagonal system of
   ! equations to be solved for 6th and 5th order accurate edge and slope
   ! estimates, respectively.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:), intent(in) :: h
      real(r8), dimension(:), intent(inout) :: tdecoeff, tdscoeff

      real(r8), dimension(6, 6) :: a
      real(r8) :: a25sq, a26sq, h3sq, h4sq

      ! Define matrix for linear system to be solved for edge coefficients.

      a(1, 1) = c1
      a(2, 1) = - h(1)
      a(3, 1) = - a(2, 1)*h(1)
      a(4, 1) = - a(3, 1)*h(1)
      a(5, 1) = - a(4, 1)*h(1)
      a(6, 1) = - a(5, 1)*h(1)

      a(1, 2) = c1
      a(2, 2) = h(2)
      a(3, 2) = a(2, 2)*h(2)
      a(4, 2) = a(3, 2)*h(2)
      a(5, 2) = a(4, 2)*h(2)
      a(6, 2) = a(5, 2)*h(2)

      a(1, 3) = - c1
      a(2, 3) = - c1_2*a(2, 1)
      a(3, 3) = - c1_3*a(3, 1)
      a(4, 3) = - c1_4*a(4, 1)
      a(5, 3) = - c1_5*a(5, 1)
      a(6, 3) = - c1_6*a(6, 1)

      a(1, 4) = - c1
      a(2, 4) = - c1_2*a(2, 2)
      a(3, 4) = - c1_3*a(3, 2)
      a(4, 4) = - c1_4*a(4, 2)
      a(5, 4) = - c1_5*a(5, 2)
      a(6, 4) = - c1_6*a(6, 2)

      a(1, 5) = - c1
      a(2, 5) = - h(2) - c1_2*h(3)
      a25sq = a(2, 5)*a(2, 5)
      h3sq = h(3)*h(3)
      a(3, 5) = - a25sq - c1_12*h3sq
      a(4, 5) = a(2, 5)*(a25sq + c1_4*h3sq)
      a(5, 5) = - a25sq*(a25sq + c1_2*h3sq) - c1_80*h3sq*h3sq
      a(6, 5) = a(2, 5)*(a25sq + c3_4*h3sq)*(a25sq + c1_12*h3sq)

      a(1, 6) = - c1
      a(2, 6) = - h(2) - h(3) - c1_2*h(4)
      a26sq = a(2, 6)*a(2, 6)
      h4sq = h(4)*h(4)
      a(3, 6) = - a26sq - c1_12*h4sq
      a(4, 6) = a(2, 6)*(a26sq + c1_4*h4sq)
      a(5, 6) = - a26sq*(a26sq + c1_2*h4sq) - c1_80*h4sq*h4sq
      a(6, 6) = a(2, 6)*(a26sq + c3_4*h4sq)*(a26sq + c1_12*h4sq)

      call edge_ih6_slope_ih5_coeff_common(a, tdecoeff, tdscoeff)

   end subroutine edge_ih6_slope_ih5_coeff_asymleft

   subroutine edge_ih6_slope_ih5_coeff_sym(h, tdecoeff, tdscoeff)
   ! ---------------------------------------------------------------------------
   ! With a symmetrical stencil, compute row coefficients for the tridiagonal
   ! system of equations to be solved for 6th and 5th order accurate edge and
   ! slope estimates, respectively.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:), intent(in) :: h
      real(r8), dimension(:), intent(inout) :: tdecoeff, tdscoeff

      real(r8), dimension(6, 6) :: a
      real(r8) :: a23sq, a26sq, h1sq, h4sq

      ! Define matrix for linear system to be solved for edge coefficients.

      a(1, 1) = c1
      a(2, 1) = - h(2)
      a(3, 1) = - a(2, 1)*h(2)
      a(4, 1) = - a(3, 1)*h(2)
      a(5, 1) = - a(4, 1)*h(2)
      a(6, 1) = - a(5, 1)*h(2)

      a(1, 2) = c1
      a(2, 2) = h(3)
      a(3, 2) = a(2, 2)*h(3)
      a(4, 2) = a(3, 2)*h(3)
      a(5, 2) = a(4, 2)*h(3)
      a(6, 2) = a(5, 2)*h(3)

      a(1, 3) = - c1
      a(2, 3) = c1_2*h(1) + h(2)
      a23sq = a(2, 3)*a(2, 3)
      h1sq = h(1)*h(1)
      a(3, 3) = - a23sq - c1_12*h1sq
      a(4, 3) = a(2, 3)*(a23sq + c1_4*h1sq)
      a(5, 3) = - a23sq*(a23sq + c1_2*h1sq) - c1_80*h1sq*h1sq
      a(6, 3) = a(2, 3)*(a23sq + c3_4*h1sq)*(a23sq + c1_12*h1sq)

      a(1, 4) = - c1
      a(2, 4) = - c1_2*a(2, 1)
      a(3, 4) = - c1_3*a(3, 1)
      a(4, 4) = - c1_4*a(4, 1)
      a(5, 4) = - c1_5*a(5, 1)
      a(6, 4) = - c1_6*a(6, 1)

      a(1, 5) = - c1
      a(2, 5) = - c1_2*a(2, 2)
      a(3, 5) = - c1_3*a(3, 2)
      a(4, 5) = - c1_4*a(4, 2)
      a(5, 5) = - c1_5*a(5, 2)
      a(6, 5) = - c1_6*a(6, 2)

      a(1, 6) = - c1
      a(2, 6) = - h(3) - c1_2*h(4)
      a26sq = a(2, 6)*a(2, 6)
      h4sq = h(4)*h(4)
      a(3, 6) = - a26sq - c1_12*h4sq
      a(4, 6) = a(2, 6)*(a26sq + c1_4*h4sq)
      a(5, 6) = - a26sq*(a26sq + c1_2*h4sq) - c1_80*h4sq*h4sq
      a(6, 6) = a(2, 6)*(a26sq + c3_4*h4sq)*(a26sq + c1_12*h4sq)

      call edge_ih6_slope_ih5_coeff_common(a, tdecoeff, tdscoeff)

   end subroutine edge_ih6_slope_ih5_coeff_sym

   subroutine edge_ih6_slope_ih5_coeff_asymright(h, tdecoeff, tdscoeff)
   ! ---------------------------------------------------------------------------
   ! With an asymmetrical stencil, where edge values are shifted left compared
   ! to cell mean values, compute row coefficients for the tridiagonal system of
   ! equations to be solved for 6th and 5th order accurate edge and slope
   ! estimates, respectively.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:), intent(in) :: h
      real(r8), dimension(:), intent(inout) :: tdecoeff, tdscoeff

      real(r8), dimension(6, 6) :: a
      real(r8) :: a23sq, a24sq, h1sq, h2sq

      ! Define matrix for linear system to be solved for edge coefficients.

      a(1, 1) = c1
      a(2, 1) = - h(3)
      a(3, 1) = - a(2, 1)*h(3)
      a(4, 1) = - a(3, 1)*h(3)
      a(5, 1) = - a(4, 1)*h(3)
      a(6, 1) = - a(5, 1)*h(3)

      a(1, 2) = c1
      a(2, 2) = h(4)
      a(3, 2) = a(2, 2)*h(4)
      a(4, 2) = a(3, 2)*h(4)
      a(5, 2) = a(4, 2)*h(4)
      a(6, 2) = a(5, 2)*h(4)

      a(1, 3) = - c1
      a(2, 3) = c1_2*h(1) + h(2) + h(3)
      a23sq = a(2, 3)*a(2, 3)
      h1sq = h(1)*h(1)
      a(3, 3) = - a23sq - c1_12*h1sq
      a(4, 3) = a(2, 3)*(a23sq + c1_4*h1sq)
      a(5, 3) = - a23sq*(a23sq + c1_2*h1sq) - c1_80*h1sq*h1sq
      a(6, 3) = a(2, 3)*(a23sq + c3_4*h1sq)*(a23sq + c1_12*h1sq)

      a(1, 4) = - c1
      a(2, 4) = c1_2*h(2) + h(3)
      a24sq = a(2, 4)*a(2, 4)
      h2sq = h(2)*h(2)
      a(3, 4) = - a24sq - c1_12*h2sq
      a(4, 4) = a(2, 4)*(a24sq + c1_4*h2sq)
      a(5, 4) = - a24sq*(a24sq + c1_2*h2sq) - c1_80*h2sq*h2sq
      a(6, 4) = a(2, 4)*(a24sq + c3_4*h2sq)*(a24sq + c1_12*h2sq)

      a(1, 5) = - c1
      a(2, 5) = - c1_2*a(2, 1)
      a(3, 5) = - c1_3*a(3, 1)
      a(4, 5) = - c1_4*a(4, 1)
      a(5, 5) = - c1_5*a(5, 1)
      a(6, 5) = - c1_6*a(6, 1)

      a(1, 6) = - c1
      a(2, 6) = - c1_2*a(2, 2)
      a(3, 6) = - c1_3*a(3, 2)
      a(4, 6) = - c1_4*a(4, 2)
      a(5, 6) = - c1_5*a(5, 2)
      a(6, 6) = - c1_6*a(6, 2)

      call edge_ih6_slope_ih5_coeff_common(a, tdecoeff, tdscoeff)

   end subroutine edge_ih6_slope_ih5_coeff_asymright

   subroutine edge_eh4_lblu(h, a)
   ! ---------------------------------------------------------------------------
   ! Compute LU matrix for explicitly estimating 4th order accurate left
   ! boundary edge value.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:), intent(in) :: h
      real(r8), dimension(:, :), intent(inout) :: a

      real(r8) :: a22sq, a32sq, a42sq, h2sq, h3sq, h4sq

      ! Define matrix for linear system to be solved for edge value.

      a(1:4, 1) = c1

      a(1, 2) = c1_2*h(1)
      a(2, 2) = a(1, 2) + c1_2*(h(1) + h(2))
      a(3, 2) = a(2, 2) + c1_2*(h(2) + h(3))
      a(4, 2) = a(3, 2) + c1_2*(h(3) + h(4))

      a22sq = a(2, 2)*a(2, 2)
      a32sq = a(3, 2)*a(3, 2)
      a42sq = a(4, 2)*a(4, 2)
      h2sq = h(2)*h(2)
      h3sq = h(3)*h(3)
      h4sq = h(4)*h(4)

      a(1, 3) = c1_3*a(1, 2)*h(1)
      a(2, 3) = c1_2*(a22sq + c1_12*h2sq)
      a(3, 3) = c1_2*(a32sq + c1_12*h3sq)
      a(4, 3) = c1_2*(a42sq + c1_12*h4sq)

      a(1, 4) = c1_4*a(1, 3)*h(1)
      a(2, 4) = c1_6*a(2, 2)*(a22sq + c1_4*h2sq)
      a(3, 4) = c1_6*a(3, 2)*(a32sq + c1_4*h3sq)
      a(4, 4) = c1_6*a(4, 2)*(a42sq + c1_4*h4sq)

      ! LU decomposition.
      call lu_decompose(4, a)

   end subroutine edge_eh4_lblu

   subroutine edge_eh4_rblu(h, a)
   ! ---------------------------------------------------------------------------
   ! Compute LU matrix for explicitly estimating 4th order accurate right
   ! boundary edge value.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:), intent(in) :: h
      real(r8), dimension(:, :), intent(inout) :: a

      real(r8) :: a12sq, a22sq, a32sq, h1sq, h2sq, h3sq

      ! Define matrix for linear system to be solved for edge value.

      a(1:4, 1) = c1

      a(4, 2) = - c1_2*h(4)
      a(3, 2) = a(4, 2) - c1_2*(h(4) + h(3))
      a(2, 2) = a(3, 2) - c1_2*(h(3) + h(2))
      a(1, 2) = a(2, 2) - c1_2*(h(2) + h(1))

      a12sq = a(1, 2)*a(1, 2)
      a22sq = a(2, 2)*a(2, 2)
      a32sq = a(3, 2)*a(3, 2)
      h1sq = h(1)*h(1)
      h2sq = h(2)*h(2)
      h3sq = h(3)*h(3)

      a(1, 3) = c1_2*(a12sq + c1_12*h1sq)
      a(2, 3) = c1_2*(a22sq + c1_12*h2sq)
      a(3, 3) = c1_2*(a32sq + c1_12*h3sq)
      a(4, 3) = - c1_3*a(4, 2)*h(4)

      a(1, 4) = c1_6*a(1, 2)*(a12sq + c1_4*h1sq)
      a(2, 4) = c1_6*a(2, 2)*(a22sq + c1_4*h2sq)
      a(3, 4) = c1_6*a(3, 2)*(a32sq + c1_4*h3sq)
      a(4, 4) = - c1_4*a(4, 3)*h(4)

      ! LU decomposition.
      call lu_decompose(4, a)

   end subroutine edge_eh4_rblu

   subroutine edge_eh6_slope_eh5_lblu(h, a)
   ! ---------------------------------------------------------------------------
   ! Compute LU matrix for explicitly estimating 6th and 5th order accurate left
   ! edge and slope values, respectively.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:), intent(in) :: h
      real(r8), dimension(:, :), intent(inout) :: a

      real(r8) :: a22sq, a32sq, a42sq, a52sq, a62sq, &
                  h2sq, h3sq, h4sq, h5sq, h6sq

      ! Define matrix for linear system to be solved for edge and slope values.

      a(1:6, 1) = c1

      a(1, 2) = c1_2*h(1)
      a(2, 2) = a(1, 2) + c1_2*(h(1) + h(2))
      a(3, 2) = a(2, 2) + c1_2*(h(2) + h(3))
      a(4, 2) = a(3, 2) + c1_2*(h(3) + h(4))
      a(5, 2) = a(4, 2) + c1_2*(h(4) + h(5))
      a(6, 2) = a(5, 2) + c1_2*(h(5) + h(6))

      a22sq = a(2, 2)*a(2, 2)
      a32sq = a(3, 2)*a(3, 2)
      a42sq = a(4, 2)*a(4, 2)
      a52sq = a(5, 2)*a(5, 2)
      a62sq = a(6, 2)*a(6, 2)
      h2sq = h(2)*h(2)
      h3sq = h(3)*h(3)
      h4sq = h(4)*h(4)
      h5sq = h(5)*h(5)
      h6sq = h(6)*h(6)

      a(1, 3) = c1_3*a(1, 2)*h(1)
      a(2, 3) = c1_2*(a22sq + c1_12*h2sq)
      a(3, 3) = c1_2*(a32sq + c1_12*h3sq)
      a(4, 3) = c1_2*(a42sq + c1_12*h4sq)
      a(5, 3) = c1_2*(a52sq + c1_12*h5sq)
      a(6, 3) = c1_2*(a62sq + c1_12*h6sq)

      a(1, 4) = c1_4*a(1, 3)*h(1)
      a(2, 4) = c1_6*a(2, 2)*(a22sq + c1_4*h2sq)
      a(3, 4) = c1_6*a(3, 2)*(a32sq + c1_4*h3sq)
      a(4, 4) = c1_6*a(4, 2)*(a42sq + c1_4*h4sq)
      a(5, 4) = c1_6*a(5, 2)*(a52sq + c1_4*h5sq)
      a(6, 4) = c1_6*a(6, 2)*(a62sq + c1_4*h6sq)

      a(1, 5) = c1_5*a(1, 4)*h(1)
      a(2, 5) = c1_24*(a22sq*(a22sq + c1_2*h2sq) + c1_80*h2sq*h2sq)
      a(3, 5) = c1_24*(a32sq*(a32sq + c1_2*h3sq) + c1_80*h3sq*h3sq)
      a(4, 5) = c1_24*(a42sq*(a42sq + c1_2*h4sq) + c1_80*h4sq*h4sq)
      a(5, 5) = c1_24*(a52sq*(a52sq + c1_2*h5sq) + c1_80*h5sq*h5sq)
      a(6, 5) = c1_24*(a62sq*(a62sq + c1_2*h6sq) + c1_80*h6sq*h6sq)

      a(1, 6) = c1_6*a(1, 5)*h(1)
      a(2, 6) = c1_120*a(2, 2)*(a22sq + c3_4*h2sq)*(a22sq + c1_12*h2sq)
      a(3, 6) = c1_120*a(3, 2)*(a32sq + c3_4*h3sq)*(a32sq + c1_12*h3sq)
      a(4, 6) = c1_120*a(4, 2)*(a42sq + c3_4*h4sq)*(a42sq + c1_12*h4sq)
      a(5, 6) = c1_120*a(5, 2)*(a52sq + c3_4*h5sq)*(a52sq + c1_12*h5sq)
      a(6, 6) = c1_120*a(6, 2)*(a62sq + c3_4*h6sq)*(a62sq + c1_12*h6sq)

      ! LU decomposition.
      call lu_decompose(6, a)

   end subroutine edge_eh6_slope_eh5_lblu

   subroutine edge_eh6_slope_eh5_rblu(h, a)
   ! ---------------------------------------------------------------------------
   ! Compute LU matrix for explicitly estimating 6th and 5th order accurate
   ! right edge and slope values, respectively.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:), intent(in) :: h
      real(r8), dimension(:, :), intent(inout) :: a

      real(r8) :: a12sq, a22sq, a32sq, a42sq, a52sq, &
                  h1sq, h2sq, h3sq, h4sq, h5sq

      ! Define matrix for linear system to be solved for edge and slope values.

      a(1:6, 1) = c1

      a(6, 2) = - c1_2*h(6)
      a(5, 2) = a(6, 2) - c1_2*(h(6) + h(5))
      a(4, 2) = a(5, 2) - c1_2*(h(5) + h(4))
      a(3, 2) = a(4, 2) - c1_2*(h(4) + h(3))
      a(2, 2) = a(3, 2) - c1_2*(h(3) + h(2))
      a(1, 2) = a(2, 2) - c1_2*(h(2) + h(1))

      a12sq = a(1, 2)*a(1, 2)
      a22sq = a(2, 2)*a(2, 2)
      a32sq = a(3, 2)*a(3, 2)
      a42sq = a(4, 2)*a(4, 2)
      a52sq = a(5, 2)*a(5, 2)
      h1sq = h(1)*h(1)
      h2sq = h(2)*h(2)
      h3sq = h(3)*h(3)
      h4sq = h(4)*h(4)
      h5sq = h(5)*h(5)

      a(1, 3) = c1_2*(a12sq + c1_12*h1sq)
      a(2, 3) = c1_2*(a22sq + c1_12*h2sq)
      a(3, 3) = c1_2*(a32sq + c1_12*h3sq)
      a(4, 3) = c1_2*(a42sq + c1_12*h4sq)
      a(5, 3) = c1_2*(a52sq + c1_12*h5sq)
      a(6, 3) = - c1_3*a(6, 2)*h(6)

      a(1, 4) = c1_6*a(1, 2)*(a12sq + c1_4*h1sq)
      a(2, 4) = c1_6*a(2, 2)*(a22sq + c1_4*h2sq)
      a(3, 4) = c1_6*a(3, 2)*(a32sq + c1_4*h3sq)
      a(4, 4) = c1_6*a(4, 2)*(a42sq + c1_4*h4sq)
      a(5, 4) = c1_6*a(5, 2)*(a52sq + c1_4*h5sq)
      a(6, 4) = - c1_4*a(6, 3)*h(6)

      a(1, 5) = c1_24*(a12sq*(a12sq + c1_2*h1sq) + c1_80*h1sq*h1sq)
      a(2, 5) = c1_24*(a22sq*(a22sq + c1_2*h2sq) + c1_80*h2sq*h2sq)
      a(3, 5) = c1_24*(a32sq*(a32sq + c1_2*h3sq) + c1_80*h3sq*h3sq)
      a(4, 5) = c1_24*(a42sq*(a42sq + c1_2*h4sq) + c1_80*h4sq*h4sq)
      a(5, 5) = c1_24*(a52sq*(a52sq + c1_2*h5sq) + c1_80*h5sq*h5sq)
      a(6, 5) = - c1_5*a(6, 4)*h(6)

      a(1, 6) = c1_120*a(1, 2)*(a12sq + c3_4*h1sq)*(a12sq + c1_12*h1sq)
      a(2, 6) = c1_120*a(2, 2)*(a22sq + c3_4*h2sq)*(a22sq + c1_12*h2sq)
      a(3, 6) = c1_120*a(3, 2)*(a32sq + c3_4*h3sq)*(a32sq + c1_12*h3sq)
      a(4, 6) = c1_120*a(4, 2)*(a42sq + c3_4*h4sq)*(a42sq + c1_12*h4sq)
      a(5, 6) = c1_120*a(5, 2)*(a52sq + c3_4*h5sq)*(a52sq + c1_12*h5sq)
      a(6, 6) = - c1_6*a(6, 5)*h(6)

      ! LU decomposition.
      call lu_decompose(6, a)

   end subroutine edge_eh6_slope_eh5_rblu

   subroutine prepare_pqm(rcs, x_edge_src)
   ! ---------------------------------------------------------------------------
   ! Prepare reconstruction with piecewise quartics using implicit 6th order
   ! accurate edge and 5th order accurate slope estimation.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs
      real(r8), dimension(:), intent(in) :: x_edge_src

      integer, dimension(rcs%n_src_all) :: prev_index, next_index
      real(r8) :: hp, h_max, h_min
      integer :: jp, j, last_index, jf, jl, n, j_min, jn, jd, js
      integer :: first_index = 0 ! Initialized to avoid compiler warning.

      ! Exclude near-empty grid cells and establish a doubly linked list that
      ! connects the remaining grid cells.
      rcs%n_src = 0
      jp = 0
      do j = 1, rcs%n_src_all
         rcs%h_src(j) = abs(x_edge_src(j + 1) - x_edge_src(j))
         if (rcs%h_src(j) > c2*rcs%x_eps) then
            rcs%n_src = rcs%n_src + 1
            rcs%src_dst_index(j) = 1
            prev_index(j) = jp
            if (jp == 0) then
               first_index = j
            else
               next_index(jp) = j
            endif
            jp = j
         else
            rcs%src_dst_index(j) = 0
         endif
      enddo
      last_index = jp
      next_index(jp) = 0
      if (rcs%n_src < 6) return

      ! Exclude grid cells that may lead to large condition numbers for the
      ! linear systems to be solved in edge_ih6_slope_ih5_coeff_asymleft,
      ! edge_ih6_slope_ih5_coeff_sym and edge_ih6_slope_ih5_coeff_asymright.
      ! Excluded grid cells are merged with the non-excluded neighbour grid cell
      ! having the smallest grid cell width.
      jf = first_index
      outer: do
         j = jf
         hp = rcs%h_src(j)
         h_max = rcs%h_src(j)
         do n = 1, 3
            j = next_index(j)
            if (j == 0) exit outer
            hp = hp*rcs%h_src(j)
            h_max = max(h_max, rcs%h_src(j))
         enddo
         if (hp > hplim_ih6*h_max**4) then
            jf = next_index(jf)
         else
            rcs%n_src = rcs%n_src - 1
            if (rcs%n_src < 6) return
            j = jf
            h_min = rcs%h_src(j)
            j_min = j
            do n = 1, 3
               j = next_index(j)
               if (rcs%h_src(j) < h_min) then
                  h_min = rcs%h_src(j)
                  j_min = j
               endif
            enddo
            jp = prev_index(j_min)
            jn = next_index(j_min)
            if     (jp == 0) then
               rcs%src_dst_index(j_min) = - jn
               rcs%h_src(jn) = rcs%h_src(jn) + rcs%h_src(j_min)
               first_index = jn
               prev_index(jn) = 0
               jf = jn
            elseif (jn == 0) then
               rcs%src_dst_index(j_min) = - jp
               rcs%h_src(jp) = rcs%h_src(jp) + rcs%h_src(j_min)
               next_index(jp) = 0
               last_index = jp
               exit
            else
               if (rcs%h_src(jn) < rcs%h_src(jp)) then
                  rcs%src_dst_index(j_min) = - jn
                  rcs%h_src(jn) = rcs%h_src(jn) + rcs%h_src(j_min)
               else
                  rcs%src_dst_index(j_min) = - jp
                  rcs%h_src(jp) = rcs%h_src(jp) + rcs%h_src(j_min)
               endif
               next_index(jp) = jn
               prev_index(jn) = jp
               jf = jp
               if (jf /= first_index) then
                  jf = prev_index(jf)
                  if (jf /= first_index) jf = prev_index(jf)
               endif
            endif
         endif
      enddo outer

      ! Exclude grid cells that may lead to a large condition number for the
      ! linear system to be solved in edge_eh6_slope_eh5_lblu. Excluded grid
      ! cells are merged with the non-excluded neighbour grid cell having the
      ! smallest grid cell width.
      jf = first_index
      do
         j = jf
         hp = rcs%h_src(j)
         h_max = rcs%h_src(j)
         do n = 1, 5
            j = next_index(j)
            hp = hp*rcs%h_src(j)
            h_max = max(h_max, rcs%h_src(j))
         enddo
         if (hp > hplim_eh6*h_max**6) then
            exit
         else
            rcs%n_src = rcs%n_src - 1
            if (rcs%n_src < 6) return
            j = jf
            h_min = rcs%h_src(j)
            j_min = j
            do n = 1, 5
               j = next_index(j)
               if (rcs%h_src(j) < h_min) then
                  h_min = rcs%h_src(j)
                  j_min = j
               endif
            enddo
            jp = prev_index(j_min)
            jn = next_index(j_min)
            if (jp == 0) then
               rcs%src_dst_index(j_min) = - jn
               rcs%h_src(jn) = rcs%h_src(jn) + rcs%h_src(j_min)
               first_index = jn
               prev_index(jn) = 0
               jf = jn
            else
               if (rcs%h_src(jn) < rcs%h_src(jp)) then
                  rcs%src_dst_index(j_min) = - jn
                  rcs%h_src(jn) = rcs%h_src(jn) + rcs%h_src(j_min)
               else
                  rcs%src_dst_index(j_min) = - jp
                  rcs%h_src(jp) = rcs%h_src(jp) + rcs%h_src(j_min)
               endif
               next_index(jp) = jn
               prev_index(jn) = jp
            endif
         endif
      enddo

      ! Exclude grid cells that may lead to a large condition number for the
      ! linear system to be solved in edge_eh6_slope_eh5_rblu. Excluded grid
      ! cells are merged with the non-excluded neighbour grid cell having the
      ! smallest grid cell width.
      jl = last_index
      do
         j = jl
         hp = rcs%h_src(j)
         h_max = rcs%h_src(j)
         do n = 1, 5
            j = prev_index(j)
            hp = hp*rcs%h_src(j)
            h_max = max(h_max, rcs%h_src(j))
         enddo
         if (hp > hplim_eh6*h_max**6) then
            exit
         else
            rcs%n_src = rcs%n_src - 1
            if (rcs%n_src < 6) return
            j = jl
            h_min = rcs%h_src(j)
            j_min = j
            do n = 1, 5
               j = prev_index(j)
               if (rcs%h_src(j) < h_min) then
                  h_min = rcs%h_src(j)
                  j_min = j
               endif
            enddo
            jp = prev_index(j_min)
            jn = next_index(j_min)
            if (jn == 0) then
               rcs%src_dst_index(j_min) = - jp
               rcs%h_src(jp) = rcs%h_src(jp) + rcs%h_src(j_min)
               next_index(jp) = 0
               jl = jp
            else
               if (rcs%h_src(jn) < rcs%h_src(jp)) then
                  rcs%src_dst_index(j_min) = - jn
                  rcs%h_src(jn) = rcs%h_src(jn) + rcs%h_src(j_min)
               else
                  rcs%src_dst_index(j_min) = - jp
                  rcs%h_src(jp) = rcs%h_src(jp) + rcs%h_src(j_min)
               endif
               next_index(jp) = jn
               prev_index(jn) = jp
            endif
         endif
      enddo

      ! For the non-excluded grid cells, assign the destination index in the
      ! continuous array of grid cells to be used in the reconstruction. Also
      ! set the grid cell widths of the continuous array.
      jd = 0
      do js = 1, rcs%n_src_all
         if (rcs%src_dst_index(js) > 0) then
            jd = jd + 1
            rcs%src_dst_index(js) = jd
            rcs%h_src(jd) = rcs%h_src(js)
            rcs%hi_src(jd) = c1/rcs%h_src(jd)
         endif
      enddo

      ! Find the destination index of excluded grid cells to be merged and
      ! compute the mapping weights.
      do js = 1, rcs%n_src_all
         jd = rcs%src_dst_index(js)
         do while (jd < 0)
            jd = rcs%src_dst_index(- jd)
         enddo
         rcs%src_dst_index(js) = jd
         if (jd > 0) &
            rcs%src_dst_weight(js) = ( x_edge_src(js + 1) &
                                     - x_edge_src(js    ))*rcs%hi_src(jd)
      enddo

      ! Set source edge values in the continuous reconstruction array.
      rcs%x_edge_src(1) = x_edge_src(1)
      js = 1
      do j = 1, rcs%n_src - 1
         do
           js = js + 1
           if (rcs%src_dst_index(js) /= j .and. rcs%src_dst_index(js) /= 0) exit
         enddo
         rcs%x_edge_src(j + 1) = x_edge_src(js)
      enddo
      rcs%x_edge_src(rcs%n_src + 1) = x_edge_src(rcs%n_src_all + 1)

      ! Compute the multiplicative inverse of cell width used for estimating
      ! centered linear slope.
      do j = 2, rcs%n_src - 1
         rcs%hci_src(j) = c2/( rcs%h_src(j - 1) + c2*rcs%h_src(j) &
                             + rcs%h_src(j + 1))
      enddo


      ! Compute coefficients for the tridiagonal system of equations for the
      ! estimation of interior edge and slope values.
      call edge_ih6_slope_ih5_coeff_asymleft( &
              rcs%h_src(1:4), &
              rcs%tdecoeff(:, 2), rcs%tdscoeff(:, 2))
      do j = 3, rcs%n_src - 1
         call edge_ih6_slope_ih5_coeff_sym( &
                 rcs%h_src((j - 2):(j + 1)), &
                 rcs%tdecoeff(:, j), rcs%tdscoeff(:, j))
      enddo
      call edge_ih6_slope_ih5_coeff_asymright( &
              rcs%h_src((rcs%n_src - 3):rcs%n_src), &
              rcs%tdecoeff(:, rcs%n_src), rcs%tdscoeff(:, rcs%n_src))

      ! Compute LU matrices for the explicit estimation of boundary edge and
      ! slope values.
      call edge_eh6_slope_eh5_lblu(rcs%h_src(1:6), &
                                   rcs%lblu)
      call edge_eh6_slope_eh5_rblu(rcs%h_src((rcs%n_src - 5):rcs%n_src), &
                                   rcs%rblu)

   end subroutine prepare_pqm

   subroutine prepare_ppm(rcs, x_edge_src)
   ! ---------------------------------------------------------------------------
   ! Prepare reconstruction with piecewise parabolas using implicit 4th order
   ! accurate edge estimation.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs
      real(r8), dimension(:), intent(in) :: x_edge_src

      integer, dimension(rcs%n_src_all) :: prev_index, next_index
      real(r8) :: hp, h_max, h_min
      integer :: jp, j, last_index, jf, jl, n, j_min, jn, jd, js
      integer :: first_index = 0 ! Initialized to avoid compiler warning.

      ! Exclude near-empty grid cells and establish a doubly linked list that
      ! connects the remaining grid cells.
      rcs%n_src = 0
      jp = 0
      do j = 1, rcs%n_src_all
         rcs%h_src(j) = abs(x_edge_src(j + 1) - x_edge_src(j))
         if (rcs%h_src(j) > c2*rcs%x_eps) then
            rcs%n_src = rcs%n_src + 1
            rcs%src_dst_index(j) = 1
            prev_index(j) = jp
            if (jp == 0) then
               first_index = j
            else
               next_index(jp) = j
            endif
            jp = j
         else
            rcs%src_dst_index(j) = 0
         endif
      enddo
      last_index = jp
      next_index(jp) = 0
      if (rcs%n_src < 4) return

      ! Exclude grid cells that may lead to large condition numbers for the
      ! linear systems to be solved in edge_ih4_coeff. Excluded grid cells are
      ! merged with the non-excluded neighbour grid cell having the smallest
      ! grid cell width.
      jf = first_index
      jl = next_index(jf)
      do
         if (rcs%h_src(jf)*rcs%h_src(jl) > &
             hplim_ih4*max(rcs%h_src(jf), rcs%h_src(jl))**2) then
            jf = jl
            jl = next_index(jf)
            if (jl == 0) exit
         else
            rcs%n_src = rcs%n_src - 1
            if (rcs%n_src < 4) return
            if (rcs%h_src(jf) < rcs%h_src(jl)) then
               j = jf
               jf = prev_index(jf)
               prev_index(jl) = jf
               if (jf == 0) then
                  rcs%src_dst_index(j) = - jl
                  rcs%h_src(jl) = rcs%h_src(jl) + rcs%h_src(j)
                  first_index = jl
                  jf = jl
                  jl = next_index(jf)
                  if (jl == 0) exit
               else
                  if (rcs%h_src(jf) < rcs%h_src(jl)) then
                     rcs%src_dst_index(j) = - jf
                     rcs%h_src(jf) = rcs%h_src(jf) + rcs%h_src(j)
                  else
                     rcs%src_dst_index(j) = - jl
                     rcs%h_src(jl) = rcs%h_src(jl) + rcs%h_src(j)
                  endif
                  next_index(jf) = jl
               endif
            else
               j = jl
               jl = next_index(jl)
               next_index(jf) = jl
               if (jl == 0) then
                  rcs%src_dst_index(j) = - jf
                  rcs%h_src(jf) = rcs%h_src(jf) + rcs%h_src(j)
                  last_index = jf
                  exit
               endif
               if (rcs%h_src(jf) < rcs%h_src(jl)) then
                  rcs%src_dst_index(j) = - jf
                  rcs%h_src(jf) = rcs%h_src(jf) + rcs%h_src(j)
               else
                  rcs%src_dst_index(j) = - jl
                  rcs%h_src(jl) = rcs%h_src(jl) + rcs%h_src(j)
               endif
               prev_index(jl) = jf
            endif
         endif
      enddo

      ! Exclude grid cells that may lead to a large condition number for the
      ! linear system to be solved in edge_eh4_lblu. Excluded grid cells are
      ! merged with the non-excluded neighbour grid cell having the smallest
      ! grid cell width.
      jf = first_index
      do
         j = jf
         hp = rcs%h_src(j)
         h_max = rcs%h_src(j)
         do n = 1, 3
            j = next_index(j)
            hp = hp*rcs%h_src(j)
            h_max = max(h_max, rcs%h_src(j))
         enddo
         if (hp > hplim_eh4*h_max**4) then
            exit
         else
            rcs%n_src = rcs%n_src - 1
            if (rcs%n_src < 4) return
            j = jf
            h_min = rcs%h_src(j)
            j_min = j
            do n = 1, 3
               j = next_index(j)
               if (rcs%h_src(j) < h_min) then
                  h_min = rcs%h_src(j)
                  j_min = j
               endif
            enddo
            jp = prev_index(j_min)
            jn = next_index(j_min)
            if (jp == 0) then
               rcs%src_dst_index(j_min) = - jn
               rcs%h_src(jn) = rcs%h_src(jn) + rcs%h_src(j_min)
               first_index = jn
               prev_index(jn) = 0
               jf = jn
            else
               if (rcs%h_src(jn) < rcs%h_src(jp)) then
                  rcs%src_dst_index(j_min) = - jn
                  rcs%h_src(jn) = rcs%h_src(jn) + rcs%h_src(j_min)
               else
                  rcs%src_dst_index(j_min) = - jp
                  rcs%h_src(jp) = rcs%h_src(jp) + rcs%h_src(j_min)
               endif
               next_index(jp) = jn
               prev_index(jn) = jp
            endif
         endif
      enddo

      ! Exclude grid cells that may lead to a large condition number for the
      ! linear system to be solved in edge_eh4_rblu. Excluded grid cells are
      ! merged with the non-excluded neighbour grid cell having the smallest
      ! grid cell width.
      jl = last_index
      do
         j = jl
         hp = rcs%h_src(j)
         h_max = rcs%h_src(j)
         do n = 1, 3
            j = prev_index(j)
            hp = hp*rcs%h_src(j)
            h_max = max(h_max, rcs%h_src(j))
         enddo
         if (hp > hplim_eh4*h_max**4) then
            exit
         else
            rcs%n_src = rcs%n_src - 1
            if (rcs%n_src < 4) return
            j = jl
            h_min = rcs%h_src(j)
            j_min = j
            do n = 1, 3
               j = prev_index(j)
               if (rcs%h_src(j) < h_min) then
                  h_min = rcs%h_src(j)
                  j_min = j
               endif
            enddo
            jp = prev_index(j_min)
            jn = next_index(j_min)
            if (jn == 0) then
               rcs%src_dst_index(j_min) = - jp
               rcs%h_src(jp) = rcs%h_src(jp) + rcs%h_src(j_min)
               next_index(jp) = 0
               jl = jp
            else
               if (rcs%h_src(jn) < rcs%h_src(jp)) then
                  rcs%src_dst_index(j_min) = - jn
                  rcs%h_src(jn) = rcs%h_src(jn) + rcs%h_src(j_min)
               else
                  rcs%src_dst_index(j_min) = - jp
                  rcs%h_src(jp) = rcs%h_src(jp) + rcs%h_src(j_min)
               endif
               next_index(jp) = jn
               prev_index(jn) = jp
            endif
         endif
      enddo

      ! For the non-excluded grid cells, assign the destination index in the
      ! continuous array of grid cells to be used in the reconstruction. Also
      ! set the grid cell widths of the continuous array.
      jd = 0
      do js = 1, rcs%n_src_all
         if (rcs%src_dst_index(js) > 0) then
            jd = jd + 1
            rcs%src_dst_index(js) = jd
            rcs%h_src(jd) = rcs%h_src(js)
            rcs%hi_src(jd) = c1/rcs%h_src(jd)
         endif
      enddo

      ! Find the destination index of excluded grid cells to be merged and
      ! compute the mapping weights.
      do js = 1, rcs%n_src_all
         jd = rcs%src_dst_index(js)
         do while (jd < 0)
            jd = rcs%src_dst_index(- jd)
         enddo
         rcs%src_dst_index(js) = jd
         if (jd > 0) &
            rcs%src_dst_weight(js) = ( x_edge_src(js + 1) &
                                     - x_edge_src(js    ))*rcs%hi_src(jd)
      enddo

      ! Set source edge values in the continuous reconstruction array.
      rcs%x_edge_src(1) = x_edge_src(1)
      js = 1
      do j = 1, rcs%n_src - 1
         do
           js = js + 1
           if (rcs%src_dst_index(js) /= j .and. rcs%src_dst_index(js) /= 0) exit
         enddo
         rcs%x_edge_src(j + 1) = x_edge_src(js)
      enddo
      rcs%x_edge_src(rcs%n_src + 1) = x_edge_src(rcs%n_src_all + 1)

      ! Compute the multiplicative inverse of cell width used for estimating
      ! centered linear slope.
      do j = 2, rcs%n_src - 1
         rcs%hci_src(j) = c2/( rcs%h_src(j - 1) + c2*rcs%h_src(j) &
                             + rcs%h_src(j + 1))
      enddo

      ! Compute coefficients for the tridiagonal system of equations for the
      ! estimation of interior edge values.
      do j = 2, rcs%n_src
         call edge_ih4_coeff(rcs%h_src((j - 1):j), rcs%tdecoeff(:, j))
      enddo

      ! Compute LU matrices for the explicit estimation of boundary edge values.
      call edge_eh4_lblu(rcs%h_src(1:4), rcs%lblu)
      call edge_eh4_rblu(rcs%h_src((rcs%n_src - 3):rcs%n_src), rcs%rblu)

   end subroutine prepare_ppm

   subroutine prepare_plm(rcs, x_edge_src)
   ! ---------------------------------------------------------------------------
   ! Prepare reconstruction with piecewise lines.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs
      real(r8), dimension(:), intent(in) :: x_edge_src

      integer :: j, js

      ! Exclude near-empty grid cells and assign the destination index in the
      ! continuous array of grid cells to be used in the reconstruction.
      rcs%n_src = 0
      do j = 1, rcs%n_src_all
         if (abs(x_edge_src(j + 1) - x_edge_src(j)) > c2*rcs%x_eps) then
            rcs%n_src = rcs%n_src + 1
            rcs%src_dst_index(j) = rcs%n_src
         else
            rcs%src_dst_index(j) = 0
         endif
      enddo
      if (rcs%n_src < 2) return

      ! Set source edge values in the continuous reconstruction array.
      rcs%x_edge_src(1) = x_edge_src(1)
      js = 1
      do j = 1, rcs%n_src - 1
         do
            js = js + 1
            if (rcs%src_dst_index(js) /= 0) exit
         enddo
         rcs%x_edge_src(j + 1) = x_edge_src(js)
      enddo
      rcs%x_edge_src(rcs%n_src + 1) = x_edge_src(rcs%n_src_all + 1)

      ! From edge locations, obtain source grid cell widths and their
      ! multiplicative inverse.
      do j = 1, rcs%n_src
         rcs%h_src(j) = abs(rcs%x_edge_src(j + 1) - rcs%x_edge_src(j))
         rcs%hi_src(j) = c1/rcs%h_src(j)
      enddo

      ! Compute the multiplicative inverse of cell width used for estimating
      ! centered linear slope.
      do j = 2, rcs%n_src - 1
         rcs%hci_src(j) = c2/( rcs%h_src(j - 1) + c2*rcs%h_src(j) &
                             + rcs%h_src(j + 1))
      enddo

   end subroutine prepare_plm

   subroutine prepare_pcm(rcs, x_edge_src)
   ! ---------------------------------------------------------------------------
   ! Prepare piecewise constant reconstruction.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs
      real(r8), dimension(:), intent(in) :: x_edge_src

      integer :: j, js

      ! Exclude near-empty grid cells and assign the destination index in the
      ! continuous array of grid cells to be used in the reconstruction.
      rcs%n_src = 0
      do j = 1, rcs%n_src_all
         if (abs(x_edge_src(j + 1) - x_edge_src(j)) > c2*rcs%x_eps) then
            rcs%n_src = rcs%n_src + 1
            rcs%src_dst_index(j) = rcs%n_src
         else
            rcs%src_dst_index(j) = 0
         endif
      enddo
      if (rcs%n_src == 0) return

      ! Set source edge values in the continuous reconstruction array.
      rcs%x_edge_src(1) = x_edge_src(1)
      js = 1
      do j = 1, rcs%n_src - 1
         do
            js = js + 1
            if (rcs%src_dst_index(js) /= 0) exit
         enddo
         rcs%x_edge_src(j + 1) = x_edge_src(js)
      enddo
      rcs%x_edge_src(rcs%n_src + 1) = x_edge_src(rcs%n_src_all + 1)

      ! From edge locations, obtain source grid cell widths and their
      ! multiplicative inverse.
      do j = 1, rcs%n_src
         rcs%h_src(j) = abs(rcs%x_edge_src(j + 1) - rcs%x_edge_src(j))
         rcs%hi_src(j) = c1/rcs%h_src(j)
      enddo

   end subroutine prepare_pcm

   subroutine reconstruct_plm_no_limiting(rcs)
   ! ---------------------------------------------------------------------------
   ! Carry out a reconstruction with piecewise lines.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs

      real(r8) :: sc
      integer :: ns, j

      ns = rcs%n_src

      sc = c2*(rcs%u_src(2) - rcs%u_src(1)) &
             /(rcs%h_src(2) + rcs%h_src(1))
      rcs%polycoeff(2, 1) = sc*rcs%h_src(1)
      rcs%polycoeff(1, 1) = rcs%u_src(1) - c1_2*rcs%polycoeff(2, 1)
      rcs%uel(1) = rcs%polycoeff(1, 1)
      rcs%uer(1) = rcs%polycoeff(1, 1) + rcs%polycoeff(2, 1)
      do j = 2, ns - 1
         sc = (rcs%u_src(j + 1) - rcs%u_src(j - 1))*rcs%hci_src(j)
         rcs%polycoeff(2, j) = sc*rcs%h_src(j)
         rcs%polycoeff(1, j) = rcs%u_src(j) - c1_2*rcs%polycoeff(2, j)
         rcs%uel(j) = rcs%polycoeff(1, j)
         rcs%uer(j) = rcs%polycoeff(1, j) + rcs%polycoeff(2, j)
      enddo
      sc = c2*(rcs%u_src(ns) - rcs%u_src(ns - 1)) &
             /(rcs%h_src(ns) + rcs%h_src(ns - 1))
      rcs%polycoeff(2, ns) = sc*rcs%h_src(ns)
      rcs%polycoeff(1, ns) = rcs%u_src(ns) - c1_2*rcs%polycoeff(2, ns)
      rcs%uel(ns) = rcs%polycoeff(1, ns)
      rcs%uer(ns) = rcs%polycoeff(1, ns) + rcs%polycoeff(2, ns)

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
         rcs%uel(j) = rcs%polycoeff(1, j)
         rcs%uer(j) = rcs%polycoeff(1, j) + rcs%polycoeff(2, j)
      enddo

      if (pc_boundary_cells) then
         rcs%polycoeff(1, 1) = rcs%u_src(1)
         rcs%polycoeff(2, 1) = c0
         rcs%uel(1) = rcs%u_src(1)
         rcs%uer(1) = rcs%u_src(1)
         rcs%polycoeff(1, ns) = rcs%u_src(ns)
         rcs%polycoeff(2, ns) = c0
         rcs%uel(ns) = rcs%u_src(ns)
         rcs%uer(ns) = rcs%u_src(ns)
      else
         sc = c2*(rcs%u_src(2) - rcs%u_src(1)) &
                /(rcs%h_src(2) + rcs%h_src(1))
         rcs%polycoeff(2, 1) = sc*rcs%h_src(1)
         rcs%polycoeff(1, 1) = rcs%u_src(1) - c1_2*rcs%polycoeff(2, 1)
         rcs%uel(1) = rcs%polycoeff(1, 1)
         rcs%uer(1) = rcs%polycoeff(1, 1) + rcs%polycoeff(2, 1)
         sc = c2*(rcs%u_src(ns) - rcs%u_src(ns - 1)) &
                /(rcs%h_src(ns) + rcs%h_src(ns - 1))
         rcs%polycoeff(2, ns) = sc*rcs%h_src(ns)
         rcs%polycoeff(1, ns) = rcs%u_src(ns) &
                              - c1_2*rcs%polycoeff(2, ns)
         rcs%uel(ns) = rcs%polycoeff(1, ns)
         rcs%uer(ns) = rcs%polycoeff(1, ns) + rcs%polycoeff(2, ns)
      endif

      rcs%polycoeff(3, :) = c0

   end subroutine reconstruct_plm_monotonic

   subroutine reconstruct_ppm_edge_values(rcs)
   ! ---------------------------------------------------------------------------
   ! Reconstruct edge values using an implicit 4th order scheme.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs

      real(r8), dimension(4) :: x
      real(r8), dimension(rcs%n_src + 1) :: uedge
      real(r8), dimension(rcs%n_src) :: rhs, gam
      real(r8) :: bei
      integer :: ns, j

      ns = rcs%n_src

      ! Obtain the left boundary edge value.
      x(:) = rcs%u_src(1:4)
      call lu_solve(4, rcs%lblu, x)
      uedge(1) = x(1)

      ! Obtain the right boundary edge value.
      x(:) = rcs%u_src((ns - 3):ns)
      call lu_solve(4, rcs%rblu, x)
      uedge(ns + 1) = x(1)

      ! Obtain right hand side of tridiagonal system of equations.
      do j = 2, ns
         rhs(j) = rcs%tdecoeff(3, j)*rcs%u_src(j - 1) &
                + rcs%tdecoeff(4, j)*rcs%u_src(j)
      enddo

      ! Solve tridiagonal system of equations to obtain interior edge values.
      gam(1) = c0
      do j = 2, ns
         bei = c1/(c1 - rcs%tdecoeff(1, j)*gam(j - 1))
         uedge(j) = (rhs(j) - rcs%tdecoeff(1, j)*uedge(j - 1))*bei
         gam(j) = rcs%tdecoeff(2, j)*bei
      enddo
      do j = ns, 2, - 1
         uedge(j) = uedge(j) - gam(j)*uedge(j + 1)
      enddo

      ! Set left and right edge values for each grid cell.
      rcs%uel(1:ns) = uedge(1:ns)
      rcs%uer(1:ns) = uedge(2:(ns + 1))

   end subroutine reconstruct_ppm_edge_values

   subroutine reconstruct_pqm_edge_slope_values(rcs)
   ! ---------------------------------------------------------------------------
   ! Reconstruct edge and slope values using implicit 6th and 5th order schemes,
   ! respectively.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs

      real(r8), dimension(6) :: x
      real(r8), dimension(rcs%n_src + 1) :: uedge, uslope
      real(r8), dimension(rcs%n_src) :: rhs, gam
      real(r8) :: bei
      integer :: ns, j

      ns = rcs%n_src

      ! Obtain the left boundary edge and slope values.
      x(:) = rcs%u_src(1:6)
      call lu_solve(6, rcs%lblu, x)
      uedge(1) = x(1)
      uslope(1) = x(2)

      ! Obtain the right boundary edge and slope values.
      x(:) = rcs%u_src((ns - 5):ns)
      call lu_solve(6, rcs%rblu, x)
      uedge(ns + 1) = x(1)
      uslope(ns + 1) = x(2)

      ! Obtain right hand side of tridiagonal system of equations for edge
      ! values.
      rhs(2) = rcs%tdecoeff(3, 2)*rcs%u_src(1) &
             + rcs%tdecoeff(4, 2)*rcs%u_src(2) &
             + rcs%tdecoeff(5, 2)*rcs%u_src(3) &
             + rcs%tdecoeff(6, 2)*rcs%u_src(4)
      do j = 3, ns - 1
         rhs(j) = rcs%tdecoeff(3, j)*rcs%u_src(j - 2) &
                + rcs%tdecoeff(4, j)*rcs%u_src(j - 1) &
                + rcs%tdecoeff(5, j)*rcs%u_src(j) &
                + rcs%tdecoeff(6, j)*rcs%u_src(j + 1)
      enddo
      rhs(ns) = rcs%tdecoeff(3, ns)*rcs%u_src(ns - 3) &
              + rcs%tdecoeff(4, ns)*rcs%u_src(ns - 2) &
              + rcs%tdecoeff(5, ns)*rcs%u_src(ns - 1) &
              + rcs%tdecoeff(6, ns)*rcs%u_src(ns)

      ! Solve tridiagonal system of equations to obtain interior edge values.
      gam(1) = c0
      do j = 2, ns
         bei = c1/(c1 - rcs%tdecoeff(1, j)*gam(j - 1))
         uedge(j) = (rhs(j) - rcs%tdecoeff(1, j)*uedge(j - 1))*bei
         gam(j) = rcs%tdecoeff(2, j)*bei
      enddo
      do j = ns, 2, - 1
         uedge(j) = uedge(j) - gam(j)*uedge(j + 1)
      enddo

      ! Obtain right hand side of tridiagonal system of equations for slope
      ! values.
      rhs(2) = rcs%tdscoeff(3, 2)*rcs%u_src(1) &
             + rcs%tdscoeff(4, 2)*rcs%u_src(2) &
             + rcs%tdscoeff(5, 2)*rcs%u_src(3) &
             + rcs%tdscoeff(6, 2)*rcs%u_src(4)
      do j = 3, ns - 1
         rhs(j) = rcs%tdscoeff(3, j)*rcs%u_src(j - 2) &
                + rcs%tdscoeff(4, j)*rcs%u_src(j - 1) &
                + rcs%tdscoeff(5, j)*rcs%u_src(j) &
                + rcs%tdscoeff(6, j)*rcs%u_src(j + 1)
      enddo
      rhs(ns) = rcs%tdscoeff(3, ns)*rcs%u_src(ns - 3) &
              + rcs%tdscoeff(4, ns)*rcs%u_src(ns - 2) &
              + rcs%tdscoeff(5, ns)*rcs%u_src(ns - 1) &
              + rcs%tdscoeff(6, ns)*rcs%u_src(ns)

      ! Solve tridiagonal system of equations to obtain interior slope values.
      gam(1) = c0
      do j = 2, ns
         bei = c1/(c1 - rcs%tdscoeff(1, j)*gam(j - 1))
         uslope(j) = (rhs(j) - rcs%tdscoeff(1, j)*uslope(j - 1))*bei
         gam(j) = rcs%tdscoeff(2, j)*bei
      enddo
      do j = ns, 2, - 1
         uslope(j) = uslope(j) - gam(j)*uslope(j + 1)
      enddo

      ! Set left and right edge values for each grid cell.
      rcs%uel(1:ns) = uedge(1:ns)
      rcs%uer(1:ns) = uedge(2:(ns + 1))

      ! Set left and right slope values for each grid cell and scale the slope
      ! values with the grid cell widths.
      rcs%usl(1:ns) = uslope(1:ns)*rcs%h_src(1:ns)
      rcs%usr(1:ns) = uslope(2:(ns + 1))*rcs%h_src(1:ns)

   end subroutine reconstruct_pqm_edge_slope_values

   subroutine limit_ppm_interior_monotonic(rcs)
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

   end subroutine limit_ppm_interior_monotonic

   subroutine limit_ppm_interior_non_oscillatory(rcs)
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

   end subroutine limit_ppm_interior_non_oscillatory

   subroutine limit_ppm_boundary(rcs, pc_boundary_cells)
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

   end subroutine limit_ppm_boundary

   subroutine limit_ppm_posdef(rcs)
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

   end subroutine limit_ppm_posdef

   subroutine polycoeff_ppm(rcs)
   ! ---------------------------------------------------------------------------
   ! Obtain coefficients for piecewise parabolas from grid cell means and left
   ! and right edge values.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs

      integer :: j

      do j = 1, rcs%n_src
         rcs%polycoeff(1, j) = rcs%uel(j)
         rcs%polycoeff(2, j) = c6*rcs%u_src(j) - c4*rcs%uel(j) - c2*rcs%uer(j)
         rcs%polycoeff(3, j) = c3*(rcs%uel(j) - c2*rcs%u_src(j) + rcs%uer(j))
      enddo

   end subroutine polycoeff_ppm

   subroutine polycoeff_pqm(rcs)
   ! ---------------------------------------------------------------------------
   ! Obtain coefficients for piecewise quartics from grid cell means and left
   ! and right edge and slope values.
   ! ---------------------------------------------------------------------------

      type(reconstruction_struct), intent(inout) :: rcs

      integer :: j

      do j = 1, rcs%n_src
         rcs%polycoeff(1, j) = rcs%uel(j)
         rcs%polycoeff(2, j) = rcs%usl(j)
         rcs%polycoeff(3, j) = &
              c30*rcs%u_src(j) - c18*rcs%uel(j) - c12*rcs%uer(j) &
            - c9_2*rcs%usl(j) + c3_2*rcs%usr(j) 
         rcs%polycoeff(4, j) = &
            - c60*rcs%u_src(j) + c32*rcs%uel(j) + c28*rcs%uer(j) &
            + c6*rcs%usl(j) - c4*rcs%usr(j)
         rcs%polycoeff(5, j) = &
              c30*rcs%u_src(j) - c15*(rcs%uel(j) + rcs%uer(j)) &
            - c5_2*(rcs%usl(j) - rcs%usr(j))
      enddo

   end subroutine polycoeff_pqm

   pure function parabola_intersection(pc, u, u_eps, xil, xir) result(xi)

      real(r8), dimension(:), intent(in) :: pc
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

      integer :: n_src_all, j

      rcs%prepared = .false.
      errstat = hor3map_noerr

      ! Check reconstruction method.
      if (method /= hor3map_pcm .and. method /= hor3map_plm .and. &
          method /= hor3map_ppm .and. method /= hor3map_pqm) then
         errstat = hor3map_invalid_method
         return
      endif

      ! Number of source grid cells.
      n_src_all = size(x_edge_src) - 1

      ! Check that source grid edges are monotonically increasing or decreasing.
      if (x_edge_src(n_src_all + 1) - x_edge_src(1) > c0) then
         do j = 1, n_src_all
            if (x_edge_src(j + 1) < x_edge_src(j)) then
               errstat = hor3map_nonmonotonic_src_edges
               return
            endif
         enddo
      else
         do j = 1, n_src_all
            if (x_edge_src(j + 1) > x_edge_src(j)) then
               errstat = hor3map_nonmonotonic_src_edges
               return
            endif
         enddo
      endif

      ! Set small value with same dimensions as edge locations.
      rcs%x_eps = max(abs( x_edge_src(n_src_all + 1) &
                         - x_edge_src(1)), eps)*eps

      ! If needed, allocate arrays in reconstruction data structure.
      if (.not. rcs%alloced) then
         rcs%n_src_all = n_src_all
         errstat = allocate_rcs(rcs)
         if (errstat /= hor3map_noerr) return
      elseif (rcs%n_src_all /= n_src_all) then
         call free_rcs(rcs)
         rcs%n_src_all = n_src_all
         errstat = allocate_rcs(rcs)
         if (errstat /= hor3map_noerr) return
      endif
#ifdef DEBUG
      rcs%x_edge_src(:) = ieee_value(1._r8, ieee_signaling_nan)
      rcs%h_src(:) = ieee_value(1._r8, ieee_signaling_nan)
      rcs%hi_src(:) = ieee_value(1._r8, ieee_signaling_nan)
      rcs%hci_src(:) = ieee_value(1._r8, ieee_signaling_nan)
      rcs%src_dst_weight(:) = ieee_value(1._r8, ieee_signaling_nan)
      rcs%tdecoeff(:, :) = ieee_value(1._r8, ieee_signaling_nan)
      rcs%tdscoeff(:, :) = ieee_value(1._r8, ieee_signaling_nan)
      rcs%lblu(:, :) = ieee_value(1._r8, ieee_signaling_nan)
      rcs%rblu(:, :) = ieee_value(1._r8, ieee_signaling_nan)
      rcs%u_src(:) = ieee_value(1._r8, ieee_signaling_nan)
      rcs%uel(:) = ieee_value(1._r8, ieee_signaling_nan)
      rcs%uer(:) = ieee_value(1._r8, ieee_signaling_nan)
      rcs%src_dst_index = -9999
      rcs%polycoeff(:, :) = ieee_value(1._r8, ieee_signaling_nan)
#endif

      ! Based on the requested reconstruction method, prepare the data structure
      ! for the various methods. Arrays with indices and weights are
      ! constructed that will map the source data to a continuous array of
      ! grid cells that are non-empty and with widths that will ensure
      ! condition numbers below a specified threshold of matrices in linear
      ! equation systems to be solved. If insufficient grid cells are available
      ! for the requested method, lower order methods are tried.

      rcs%method = method

      if (rcs%method == hor3map_pqm) then
         call prepare_pqm(rcs, x_edge_src)
         if (rcs%n_src < 6) then
            rcs%method = hor3map_ppm
         else
            rcs%method = hor3map_pqm
         endif
      endif

      if (rcs%method == hor3map_ppm) then
         call prepare_ppm(rcs, x_edge_src)
         if (rcs%n_src < 4) then
            rcs%method = hor3map_plm
         else
            rcs%method = hor3map_ppm
         endif
      endif

      if (rcs%method == hor3map_plm) then
         call prepare_plm(rcs, x_edge_src)
         if (rcs%n_src < 2) then
            rcs%method = hor3map_pcm
         else
            rcs%method = hor3map_plm
         endif
      endif

      if (rcs%method == hor3map_pcm) then
         call prepare_pcm(rcs, x_edge_src)
         if (rcs%n_src == 0) then
            errstat = hor3map_src_extent_too_small
            return
         else
            rcs%method = hor3map_pcm
         endif
      endif

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

      rms%prepared = .false.
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
      if (.not. rms%alloced) then
         rms%n_dst = n_dst
         errstat = allocate_rms(rcs, rms)
         if (errstat /= hor3map_noerr) return
      elseif (rms%n_dst /= n_dst) then
         call free_rms(rms)
         rms%n_dst = n_dst
         errstat = allocate_rms(rcs, rms)
         if (errstat /= hor3map_noerr) return
      endif
#ifdef DEBUG
      rms%h_dst(:) = ieee_value(1._r8, ieee_signaling_nan)
      rms%hi_dst(:) = ieee_value(1._r8, ieee_signaling_nan)
      rms%seg_int_lim(:) = ieee_value(1._r8, ieee_signaling_nan)
      rms%n_src_seg(:) = -9999
      rms%seg_dst_index(:) = -9999
#endif

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

      integer :: js, jd

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

      ! Copy source data array to continuous array of grid cells to be used in
      ! the reconstruction.
      if (rcs%method == hor3map_pcm .or. rcs%method == hor3map_plm) then
         do js = 1, rcs%n_src_all
            jd = rcs%src_dst_index(js)
            if (jd /= 0) rcs%u_src(jd) = u_src(js)
         enddo
      else
         rcs%u_src(1:rcs%n_src) = c0
         do js = 1, rcs%n_src_all
            jd = rcs%src_dst_index(js)
            if (jd /= 0) rcs%u_src(jd) = rcs%u_src(jd) &
                                       + rcs%src_dst_weight(js)*u_src(js)
         enddo
      endif

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
                  call limit_ppm_interior_monotonic(rcs)
                  call limit_ppm_boundary(rcs, pc_boundary_cells)
               case (hor3map_non_oscillatory)
                  call limit_ppm_interior_non_oscillatory(rcs)
                  call limit_ppm_boundary(rcs, pc_boundary_cells)
               case (hor3map_non_oscillatory_posdef)
                  call limit_ppm_interior_non_oscillatory(rcs)
                  call limit_ppm_boundary(rcs, pc_boundary_cells)
                  call limit_ppm_posdef(rcs)
               case default
                  errstat = hor3map_invalid_ppm_limiting
                  return
            end select
            call polycoeff_ppm(rcs)
         case (hor3map_pqm)
            call reconstruct_pqm_edge_slope_values(rcs)
            select case (limiting)
               case (hor3map_no_limiting)
               case default
                  errstat = hor3map_invalid_pqm_limiting
                  return
            end select
            call polycoeff_pqm(rcs)
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
      u_min = minval(rcs%u_src(1:rcs%n_src))
      u_max = maxval(rcs%u_src(1:rcs%n_src))
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

         case (hor3map_pqm)

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
                        u_dst(jd) =       rcs%polycoeff(1, js) &
                                  + (     rcs%polycoeff(2, js) &
                                    + (   rcs%polycoeff(3, js) &
                                      + ( rcs%polycoeff(4, js) &
                                        + rcs%polycoeff(5, js) &
                                    *xir)*xir)*xir)*xir
                     else
                        adr = (            rcs%polycoeff(1, js) &
                              + (     c1_2*rcs%polycoeff(2, js) &
                                + (   c1_3*rcs%polycoeff(3, js) &
                                  + ( c1_4*rcs%polycoeff(4, js) &
                                    + c1_5*rcs%polycoeff(5, js) &
                              *xir)*xir)*xir)*xir)*xir
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
                 rcs%src_dst_weight, rcs%tdecoeff, rcs%tdscoeff, rcs%lblu, &
                 rcs%rblu, rcs%u_src, rcs%uel, rcs%uer, rcs%usl, rcs%usr, &
                 rcs%src_dst_index, rcs%polycoeff)

      rcs%alloced = .false.
      rcs%reconstructed = .false.

   end subroutine free_rcs

   subroutine free_rms(rms)
   ! ---------------------------------------------------------------------------
   ! Deallocate arrays and reset flags.
   ! ---------------------------------------------------------------------------

      type(remap_struct), intent(inout) :: rms

      deallocate(rms%h_dst, rms%hi_dst, rms%seg_int_lim, rms%n_src_seg, &
                 rms%seg_dst_index)

      rms%alloced = .false.

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
