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

#define DEBUG
#undef DEBUG

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
      hor3map_non_oscillatory_posdef = 204, &
      hor3map_regrid_method_1        = 301, & ! Regrid methods
      hor3map_regrid_method_2        = 302

   ! Error parameters.
   integer, parameter :: &
      hor3map_noerr                     =  0, &
      hor3map_invalid_recon_method      =  1, &
      hor3map_resizing_initialized_rcgs =  2, &
      hor3map_nonmonotonic_src_edges    =  3, &
      hor3map_src_extent_too_small      =  4, &
      hor3map_failed_to_allocate_rcgs   =  5, &
      hor3map_recon_not_prepared        =  6, &
      hor3map_resizing_initialized_rms  =  7, &
      hor3map_inconsistent_grid_range   =  8, &
      hor3map_nonmonotonic_dst_edges    =  9, &
      hor3map_failed_to_allocate_rms    = 10, &
      hor3map_src_size_mismatch         = 11, &
      hor3map_failed_to_allocate_rcss   = 12, &
      hor3map_invalid_plm_limiting      = 13, &
      hor3map_invalid_ppm_limiting      = 14, &
      hor3map_invalid_pqm_limiting      = 15, &
      hor3map_recon_not_available       = 16, &
      hor3map_invalid_regrid_method     = 17, &
      hor3map_grd_size_mismatch         = 18, &
      hor3map_remap_not_prepared        = 19, &
      hor3map_dst_size_mismatch         = 20, &
      hor3map_index_out_of_bounds       = 21, &
      hor3map_inconsistent_rcgs         = 22, &
      hor3map_errmsg_num = 22
   character(len = 80), dimension(hor3map_errmsg_num), parameter :: errmsg = &
      ["Invalid reconstruction method!                                   ", &
       "Cannot resize initialized reconstruction grid data structure!    ", &
       "Source grid edges do not monotonically increase or decrease!     ", &
       "Source grid extent too small!                                    ", &
       "Failed to allocate reconstruction grid data structure!           ", &
       "Call 'prepare_reconstruction' first!                             ", &
       "Cannot resize initialized remapping data structure!              ", &
       "Inconsistent source and destination grid range!                  ", &
       "Destination grid edges do not monotonically increase or decrease!", &
       "Failed to allocate remapping data structure!                     ", &
       "Size mismatch between source grid edges and data array!          ", &
       "Failed to allocate reconstruction source data structure!         ", &
       "Invalid limiting method for PLM!                                 ", &
       "Invalid limiting method for PPM!                                 ", &
       "Invalid limiting method for PQM!                                 ", &
       "Call 'reconstruct' first!                                        ", &
       "Invalid regrid method!                                           ", &
       "Size mismatch between grid edge values and locations!            ", &
       "Call 'prepare_remapping' first!                                  ", &
       "Size mismatch between destination grid edges and data array!     ", &
       "Array index of data structure is out of bounds!                  ", &
       "Inconsistent data structure for reconstruction and remapping!    "]

   ! Numeric data types.
   integer, parameter :: &
      r8 = real64

   ! Small non-dimensional value.
   real(r8), parameter :: eps = 1.e-14_r8

   ! Minimum number of source cells for reconstructions.
   integer, parameter :: &
      n_src_min_plm = 2, &
      n_src_min_ppm = 3, &
      n_src_min_pqm = 4

   ! Maximum order of explicit schemes for boundary edge and slope estimates.
   integer, parameter :: &
      eb_ord_max_ppm = 4, &
      eb_ord_max_pqm = 6

   ! Lower bounds of prod(h(i))/max(h(i))^n, where h(i) are cell widths
   ! belonging to a stencil of n grid cells. Respecting these bounds will ensure
   ! acceptable condition numbers of matrices involved in various linear
   ! equation systems.
   real(r8), parameter :: &
      hplim_ih4 = 5.e-7_r8, &
      hplim_ih6 = 1.e-7_r8
   real(r8), parameter, dimension(6) :: &
      hplim_eb = [1.e-10_r8, 1.e-10_r8, 1.e-10_r8, &
                  1.e-10_r8,  1.e-8_r8,  1.e-7_r8]

   ! Numeric constants.
   real(r8), parameter :: &
      c0 = 0._r8, c1 = 1._r8, c2 = 2._r8, c3 = 3._r8, c4 = 4._r8, c5 = 5._r8, &
      c6 = 6._r8, c10 = 10._r8, c12 = 12._r8, c15 = 15._r8, c18 = 18._r8, &
      c20 = 20._r8, c28 = 28._r8, c30 = 30._r8, c32 = 32._r8, c60 = 60._r8, &
      c1_2 = 1._r8/2._r8, c1_3 = 1._r8/3._r8, c1_4 = 1._r8/4._r8, &
      c1_5 = 1._r8/5._r8, c1_6 = 1._r8/6._r8, c1_8 = 1._r8/8._r8, &
      c1_12 = 1._r8/12._r8, c1_16 = 1._r8/16._r8, c1_24 = 1._r8/24._r8, &
      c1_80 = 1._r8/80._r8, c1_120 = 1._r8/120._r8, &
      c2_3 = 2._r8/3._r8, c3_4 = 3._r8/4._r8, c3_2 = 3._r8/2._r8, &
      c5_2 = 5._r8/2._r8, c8_3 = 8._r8/3._r8, c10_3 = 10._r8/3._r8, &
      c9_2 = 9._r8/2._r8

   ! Derived data types.

   type :: recon_grd_struct
      ! Reconstruction grid data structure.

      integer :: &
         i_lbound       = 1, &
         i_ubound       = 1, &
         j_lbound       = 1, &
         j_ubound       = 1, &
         i_index        = 1, &
         j_index        = 1, &
         i_index_curr   = 0, &
         j_index_curr   = 0, &
         method         = hor3map_ppm, &
         left_bndr_ord  = 0, &
         right_bndr_ord = 0
      logical :: &
         initialized  = .false.
      integer :: n_src, p_ord

      real(r8), allocatable, dimension(:,:,:) :: &
         tdecoeff_data, tdscoeff_data, lblu_data, rblu_data
      real(r8), allocatable, dimension(:,:) :: &
         x_edge_src_data, h_src_data, hi_src_data, hci_src_data, &
         src_dst_weight_data
      real(r8), allocatable, dimension(:) :: &
         x_eps_data
      integer, allocatable, dimension(:,:) :: &
         src_dst_index_data
      integer, allocatable, dimension(:) :: &
         n_src_actual_data, method_actual_data, &
         left_bndr_ord_actual_data, right_bndr_ord_actual_data
      logical, allocatable, dimension(:) :: &
         prepared_data

      real(r8), dimension(:,:), pointer :: &
         tdecoeff, tdscoeff, lblu, rblu
      real(r8), dimension(:), pointer :: &
         x_edge_src, h_src, hi_src, hci_src, &
         src_dst_weight
      real(r8), pointer :: &
         x_eps
      integer, dimension(:), pointer :: &
         src_dst_index
      integer, pointer :: &
         n_src_actual, method_actual, left_bndr_ord_actual, &
         right_bndr_ord_actual
      logical, pointer :: &
         prepared

      type(recon_src_struct), pointer :: rcss_dep_head
      type(remap_struct), pointer :: rms_dep_head

   end type recon_grd_struct

   type :: recon_src_struct
      ! Reconstruction source data structure.

      integer :: &
         limiting      = hor3map_monotonic, &
         i_index_curr  = 0, &
         j_index_curr  = 0
      logical :: &
         pc_left_bndr  = .true., &
         pc_right_bndr = .true., &
         initialized   = .false.

      real(r8), allocatable, dimension(:,:,:) :: &
         polycoeff_data
      real(r8), allocatable, dimension(:,:) :: &
         u_src_data, uel_data, uer_data, usl_data, usr_data
      real(r8), allocatable, dimension(:) :: &
         u_range_data, u_eps_data, uu_eps_data
      logical, allocatable, dimension(:) :: &
         reconstructed_data

      real(r8), dimension(:,:), pointer :: &
         polycoeff
      real(r8), dimension(:), pointer :: &
         u_src, uel, uer, usl, usr
      real(r8), pointer :: &
         u_range, u_eps, uu_eps
      logical, pointer :: &
         reconstructed

      type(recon_grd_struct), pointer :: rcgs
      type(recon_src_struct), pointer :: rcss_dep_next

   end type recon_src_struct

   type :: remap_struct
      ! Remapping data structure.

      integer :: &
         i_index_curr = 0, &
         j_index_curr = 0
      logical :: &
         initialized  = .false.
      integer :: n_dst

      real(r8), allocatable, dimension(:,:) :: &
         seg_int_lim_data, seg_weight_data
      integer, allocatable, dimension(:,:) :: &
         n_src_seg_data, seg_dst_index_data
      logical, allocatable, dimension(:) :: &
         prepared_data

      real(r8), dimension(:), pointer :: seg_int_lim, seg_weight
      integer, dimension(:), pointer :: n_src_seg, seg_dst_index
      logical, pointer :: prepared

      type(recon_grd_struct), pointer :: rcgs
      type(remap_struct), pointer :: rms_dep_next

   end type remap_struct

   public :: recon_grd_struct, recon_src_struct, remap_struct, &
             initialize_rcgs, initialize_rcss, initialize_rms, &
             prepare_reconstruction, prepare_remapping, &
             reconstruct, extract_polycoeff, regrid, remap, &
             free_rcgs, free_rcss, free_rms, &
             hor3map_pcm, hor3map_plm, hor3map_ppm, hor3map_pqm, &
             hor3map_no_limiting, hor3map_monotonic, hor3map_non_oscillatory, &
             hor3map_non_oscillatory_posdef, &
             hor3map_noerr, hor3map_errstr

contains

   ! ---------------------------------------------------------------------------
   ! Private procedures.
   ! ---------------------------------------------------------------------------

   function assign_ptr_rcgs(rcgs) result(errstat)
   ! ---------------------------------------------------------------------------
   ! Assign array pointers within reconstruction grid data structure.
   ! ---------------------------------------------------------------------------

      type(recon_grd_struct), target, intent(inout) :: rcgs

      integer :: errstat

      integer :: ij_index

      errstat = hor3map_noerr

      ! Check if new pointer assignments are needed.
      if (rcgs%i_index == rcgs%i_index_curr .and. &
          rcgs%j_index == rcgs%j_index_curr) return

      ! Check index bounds.
      if (rcgs%i_index < rcgs%i_lbound .or. rcgs%i_index > rcgs%i_ubound .or. &
          rcgs%j_index < rcgs%j_lbound .or. rcgs%j_index > rcgs%j_ubound) then
         errstat = hor3map_index_out_of_bounds
         return
      endif

      ! Assign array pointers within the reconstruction grid data structure.

      ij_index = rcgs%i_index - rcgs%i_lbound + 1 &
               + (rcgs%j_index - rcgs%j_lbound) &
                 *(rcgs%i_ubound - rcgs%i_lbound + 1)

      rcgs%x_eps => rcgs%x_eps_data(ij_index)
      rcgs%x_edge_src => rcgs%x_edge_src_data(:,ij_index)
      rcgs%h_src => rcgs%h_src_data(:,ij_index)
      rcgs%hi_src => rcgs%hi_src_data(:,ij_index)
      rcgs%src_dst_index => rcgs%src_dst_index_data(:,ij_index)
      rcgs%n_src_actual => rcgs%n_src_actual_data(ij_index)
      rcgs%method_actual => rcgs%method_actual_data(ij_index)
      rcgs%prepared => rcgs%prepared_data(ij_index)
      if (rcgs%method /= hor3map_pcm) then
         rcgs%hci_src => rcgs%hci_src_data(:,ij_index)
      endif
      if (rcgs%method == hor3map_ppm .or. rcgs%method == hor3map_pqm) then
         rcgs%src_dst_weight => rcgs%src_dst_weight_data(:,ij_index)
         rcgs%tdecoeff => rcgs%tdecoeff_data(:,:,ij_index)
         rcgs%tdscoeff => rcgs%tdscoeff_data(:,:,ij_index)
         rcgs%lblu => rcgs%lblu_data(:,:,ij_index)
         rcgs%rblu => rcgs%rblu_data(:,:,ij_index)
         rcgs%left_bndr_ord_actual => rcgs%left_bndr_ord_actual_data(ij_index)
         rcgs%right_bndr_ord_actual => rcgs%right_bndr_ord_actual_data(ij_index)
      endif

      rcgs%i_index_curr = rcgs%i_index
      rcgs%j_index_curr = rcgs%j_index

   end function assign_ptr_rcgs

   function assign_ptr_rcss(rcss) result(errstat)
   ! ---------------------------------------------------------------------------
   ! Assign array pointers within reconstruction grid and source data
   ! structures.
   ! ---------------------------------------------------------------------------

      type(recon_src_struct), target, intent(inout) :: rcss

      integer :: errstat

      integer :: ij_index

      errstat = hor3map_noerr

      ! Check if new pointer assignments are needed.
      if (rcss%rcgs%i_index == rcss%i_index_curr .and. &
          rcss%rcgs%j_index == rcss%j_index_curr) return

      ! Check index bounds.
      if (rcss%rcgs%i_index < rcss%rcgs%i_lbound .or. &
          rcss%rcgs%i_index > rcss%rcgs%i_ubound .or. &
          rcss%rcgs%j_index < rcss%rcgs%j_lbound .or. &
          rcss%rcgs%j_index > rcss%rcgs%j_ubound) then
         errstat = hor3map_index_out_of_bounds
         return
      endif

      ij_index = rcss%rcgs%i_index - rcss%rcgs%i_lbound + 1 &
               + (rcss%rcgs%j_index - rcss%rcgs%j_lbound) &
                *(rcss%rcgs%i_ubound - rcss%rcgs%i_lbound + 1)

      rcss%u_src => rcss%u_src_data(:,ij_index)
      rcss%u_range => rcss%u_range_data(ij_index)
      rcss%u_eps => rcss%u_eps_data(ij_index)
      rcss%uu_eps => rcss%uu_eps_data(ij_index)
      rcss%uel => rcss%uel_data(:,ij_index)
      rcss%uer => rcss%uer_data(:,ij_index)
      rcss%polycoeff => rcss%polycoeff_data(:,:,ij_index)
      rcss%reconstructed => rcss%reconstructed_data(ij_index)
      if (rcss%rcgs%method == hor3map_pqm) then
         rcss%usl => rcss%usl_data(:,ij_index)
         rcss%usr => rcss%usr_data(:,ij_index)
      endif

      rcss%i_index_curr = rcss%rcgs%i_index
      rcss%j_index_curr = rcss%rcgs%j_index

   end function assign_ptr_rcss

   function assign_ptr_rms(rms) result(errstat)
   ! ---------------------------------------------------------------------------
   ! Assign array pointers within the remapping data structure.
   ! ---------------------------------------------------------------------------

      type(remap_struct), target, intent(inout) :: rms

      integer :: errstat

      integer :: ij_index

      errstat = hor3map_noerr

      ! Check if new pointer assignments are needed.
      if (rms%rcgs%i_index == rms%i_index_curr .and. &
          rms%rcgs%j_index == rms%j_index_curr) return

      ! Check index bounds.
      if (rms%rcgs%i_index < rms%rcgs%i_lbound .or. &
          rms%rcgs%i_index > rms%rcgs%i_ubound .or. &
          rms%rcgs%j_index < rms%rcgs%j_lbound .or. &
          rms%rcgs%j_index > rms%rcgs%j_ubound) then
         errstat = hor3map_index_out_of_bounds
         return
      endif

      ij_index = rms%rcgs%i_index - rms%rcgs%i_lbound + 1 &
               + (rms%rcgs%j_index - rms%rcgs%j_lbound) &
                *(rms%rcgs%i_ubound - rms%rcgs%i_lbound + 1)

      rms%seg_int_lim => rms%seg_int_lim_data(:,ij_index)
      rms%seg_weight => rms%seg_weight_data(:,ij_index)
      rms%n_src_seg => rms%n_src_seg_data(:,ij_index)
      rms%seg_dst_index => rms%seg_dst_index_data(:,ij_index)
      rms%prepared => rms%prepared_data(ij_index)

      rms%i_index_curr = rms%rcgs%i_index
      rms%j_index_curr = rms%rcgs%j_index

   end function assign_ptr_rms

   pure subroutine left_bndr_cond(rcgs, prev_index, next_index, first_index, &
                                  last_index, lb_ord, ns, ns_min)
   ! ---------------------------------------------------------------------------
   ! Exclude grid cells that may lead to a large condition number for the linear
   ! system to be solved for edge and slope values at the left boundary.
   ! Excluded grid cells are merged with the non-excluded neighbour grid cell
   ! having the smallest grid cell width. Depending on the actual number of grid
   ! cells left, the order of reconstruction will be automatically adjusted.
   ! ---------------------------------------------------------------------------

      type(recon_grd_struct), intent(inout) :: rcgs
      integer, dimension(:), intent(inout) :: prev_index, next_index
      integer, intent(inout) :: first_index, last_index, lb_ord, ns
      integer, intent(in) :: ns_min

      real(r8) :: hp, h_max, h_min
      integer :: jf, j, n, j_min, jp, jn

      jf = first_index

      do
         j = jf
         hp = rcgs%h_src(j)
         h_max = rcgs%h_src(j)
         do n = 1, lb_ord-1
            j = next_index(j)
            hp = hp*rcgs%h_src(j)
            h_max = max(h_max, rcgs%h_src(j))
         enddo
         if (hp > hplim_eb(lb_ord)*h_max**lb_ord) then
            return
         else
            ns = ns - 1
            if (ns < ns_min) return
            j = jf
            h_min = rcgs%h_src(j)
            j_min = j
            do n = 1, lb_ord-1
               j = next_index(j)
               if (rcgs%h_src(j) < h_min) then
                  h_min = rcgs%h_src(j)
                  j_min = j
               endif
            enddo
            jp = prev_index(j_min)
            jn = next_index(j_min)
            if     (jp == 0) then
               rcgs%src_dst_index(j_min) = - jn
               rcgs%h_src(jn) = rcgs%h_src(jn) + rcgs%h_src(j_min)
               first_index = jn
               prev_index(jn) = 0
               jf = jn
            elseif (jn == 0) then
               rcgs%src_dst_index(j_min) = - jp
               rcgs%h_src(jp) = rcgs%h_src(jp) + rcgs%h_src(j_min)
               next_index(jp) = 0
               last_index = jp
            else
               if (rcgs%h_src(jn) < rcgs%h_src(jp)) then
                  rcgs%src_dst_index(j_min) = - jn
                  rcgs%h_src(jn) = rcgs%h_src(jn) + rcgs%h_src(j_min)
               else
                  rcgs%src_dst_index(j_min) = - jp
                  rcgs%h_src(jp) = rcgs%h_src(jp) + rcgs%h_src(j_min)
               endif
               next_index(jp) = jn
               prev_index(jn) = jp
            endif
            lb_ord = min(ns, lb_ord)
         endif
      enddo

   end subroutine left_bndr_cond

   pure subroutine right_bndr_cond(rcgs, prev_index, next_index, last_index, &
                                   rb_ord, ns, ns_min)
   ! ---------------------------------------------------------------------------
   ! Exclude grid cells that may lead to a large condition number for the linear
   ! system to be solved for edge and slope values at the right boundary.
   ! Excluded grid cells are merged with the non-excluded neighbour grid cell
   ! having the smallest grid cell width. Depending on the actual number of grid
   ! cells left, the order of reconstruction will be automatically adjusted.
   ! ---------------------------------------------------------------------------

      type(recon_grd_struct), intent(inout) :: rcgs
      integer, dimension(:), intent(inout) :: prev_index, next_index
      integer, intent(inout) :: rb_ord, ns
      integer, intent(in) :: last_index, ns_min

      real(r8) :: hp, h_max, h_min
      integer :: jl, j, n, j_min, jp, jn

      jl = last_index

      do
         j = jl
         hp = rcgs%h_src(j)
         h_max = rcgs%h_src(j)
         do n = 1, rb_ord-1
            j = prev_index(j)
            hp = hp*rcgs%h_src(j)
            h_max = max(h_max, rcgs%h_src(j))
         enddo
         if (hp > hplim_eb(rb_ord)*h_max**rb_ord) then
            return
         else
            ns = ns - 1
            if (ns < ns_min) return
            j = jl
            h_min = rcgs%h_src(j)
            j_min = j
            do n = 1, rb_ord-1
               j = prev_index(j)
               if (rcgs%h_src(j) < h_min) then
                  h_min = rcgs%h_src(j)
                  j_min = j
               endif
            enddo
            jp = prev_index(j_min)
            jn = next_index(j_min)
            if     (jp == 0) then
               rcgs%src_dst_index(j_min) = - jn
               rcgs%h_src(jn) = rcgs%h_src(jn) + rcgs%h_src(j_min)
               prev_index(jn) = 0
            elseif (jn == 0) then
               rcgs%src_dst_index(j_min) = - jp
               rcgs%h_src(jp) = rcgs%h_src(jp) + rcgs%h_src(j_min)
               next_index(jp) = 0
               jl = jp
            else
               if (rcgs%h_src(jn) < rcgs%h_src(jp)) then
                  rcgs%src_dst_index(j_min) = - jn
                  rcgs%h_src(jn) = rcgs%h_src(jn) + rcgs%h_src(j_min)
               else
                  rcgs%src_dst_index(j_min) = - jp
                  rcgs%h_src(jp) = rcgs%h_src(jp) + rcgs%h_src(j_min)
               endif
               next_index(jp) = jn
               prev_index(jn) = jp
            endif
            rb_ord = min(ns, rb_ord)
         endif
      enddo

   end subroutine right_bndr_cond

   pure subroutine lu_decompose(n, a)
   ! ---------------------------------------------------------------------------
   ! Replace the n x n input matrix A with its LU decomposition.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: n
      real(r8), dimension(:,:), intent(inout) :: a

      real(r8) :: q
      integer :: i, j, k

      do k = 1, n-1
         q = c1/a(k,k)
         do i = k+1, n
            a(i,k) = a(i,k)*q
            do j = k+1, n
               a(i,j) = a(i,j) - a(i,k)*a(k,j)
            enddo
         enddo
      enddo

   end subroutine lu_decompose

   pure subroutine lu_solve(n, lu, x)
   ! ---------------------------------------------------------------------------
   ! Solve the linear system of equations A*x = b using the LU decomposition of
   ! the n x n matrix A. The argument x has b as input and is replaced with the
   ! solution upon return.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: n
      real(r8), dimension(:,:), intent(in) :: lu
      real(r8), dimension(:), intent(inout) :: x

      integer :: i, j

      ! Forward substitution.
      do i = 2, n
         do j = 1, i-1
            x(i) = x(i) - lu(i,j)*x(j)
         enddo
      enddo

      ! Back substitution.
      x(n) = x(n)/lu(n, n)
      do i = n-1, 1, -1
         do j = i+1, n
            x(i) = x(i) - lu(i,j)*x(j)
         enddo
         x(i) = x(i)/lu(i,i)
      enddo

   end subroutine lu_solve

   pure subroutine edge_ih4_coeff(h, tdecoeff)
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

   pure subroutine slope_ih3_coeff(h, tdscoeff)
   ! ---------------------------------------------------------------------------
   ! Compute row coefficients for the tridiagonal system of equations to be
   ! solved for 3rd order accurate slope estimates.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:), intent(in) :: h
      real(r8), dimension(:), intent(inout) :: tdscoeff

      real(r8) :: h11, h22, h12, q

      h11 = h(1)*h(1)
      h22 = h(2)*h(2)
      h12 = h(1)*h(2)
      q = c1/((h(1) + h(2))*(h11 + c3*h12 + h22))
      tdscoeff(1) = h(2)*(h11 + h(2)*(h(1) - h(2)))*q
      tdscoeff(2) = h(1)*(h22 + h(1)*(h(2) - h(1)))*q
      tdscoeff(3) = - c12*h12*q
      tdscoeff(4) = - tdscoeff(3)

   end subroutine slope_ih3_coeff

   pure subroutine edge_ih6_slope_ih5_coeff_common(a, tdecoeff, tdscoeff)
   ! ---------------------------------------------------------------------------
   ! Common procedure for the various stencils for the computation of row
   ! coefficients for the tridiagonal system of equations to be solved for 6th
   ! and 5th order accurate edge and slope estimates, respectively.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:,:), intent(inout) :: a
      real(r8), dimension(:), intent(inout) :: tdecoeff, tdscoeff

      real(r8), dimension(6,6) :: b

      ! Define matrix for linear system to be solved for slope coefficients.

      b(1:5,3:6) = a(2:6,3:6)

      b(1,1) = c1
      b(2,1) = c2*a(2,1)
      b(3,1) = c3*a(3,1)
      b(4,1) = c4*a(4,1)
      b(5,1) = c5*a(5,1)
      b(6,1) = c0

      b(1,2) = c1
      b(2,2) = c2*a(2,2)
      b(3,2) = c3*a(3,2)
      b(4,2) = c4*a(4,2)
      b(5,2) = c5*a(5,2)
      b(6,2) = c0

      b(6,3:6) = c1

      ! Solve linear system for edge coefficients.
      tdecoeff(:) = [ - c1, c0, c0, c0, c0, c0]
      call lu_decompose(6, a)
      call lu_solve(6, a, tdecoeff)

      ! Solve linear system for slope coefficients.
      tdscoeff(:) = [ - c1, c0, c0, c0, c0, c0]
      call lu_decompose(6, b)
      call lu_solve(6, b, tdscoeff)

   end subroutine edge_ih6_slope_ih5_coeff_common

   pure subroutine edge_ih6_slope_ih5_coeff_asymleft(h, tdecoeff, tdscoeff)
   ! ---------------------------------------------------------------------------
   ! With an asymmetrical stencil, where edge values are shifted left compared
   ! to cell mean values, compute row coefficients for the tridiagonal system of
   ! equations to be solved for 6th and 5th order accurate edge and slope
   ! estimates, respectively.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:), intent(in) :: h
      real(r8), dimension(:), intent(inout) :: tdecoeff, tdscoeff

      real(r8), dimension(6,6) :: a
      real(r8) :: a25sq, a26sq, h3sq, h4sq

      ! Define matrix for linear system to be solved for edge coefficients.

      a(1,1) = c1
      a(2,1) = - h(1)
      a(3,1) = - a(2,1)*h(1)
      a(4,1) = - a(3,1)*h(1)
      a(5,1) = - a(4,1)*h(1)
      a(6,1) = - a(5,1)*h(1)

      a(1,2) = c1
      a(2,2) = h(2)
      a(3,2) = a(2,2)*h(2)
      a(4,2) = a(3,2)*h(2)
      a(5,2) = a(4,2)*h(2)
      a(6,2) = a(5,2)*h(2)

      a(1,3) = - c1
      a(2,3) = - c1_2*a(2,1)
      a(3,3) = - c1_3*a(3,1)
      a(4,3) = - c1_4*a(4,1)
      a(5,3) = - c1_5*a(5,1)
      a(6,3) = - c1_6*a(6,1)

      a(1,4) = - c1
      a(2,4) = - c1_2*a(2,2)
      a(3,4) = - c1_3*a(3,2)
      a(4,4) = - c1_4*a(4,2)
      a(5,4) = - c1_5*a(5,2)
      a(6,4) = - c1_6*a(6,2)

      a(1,5) = - c1
      a(2,5) = - h(2) - c1_2*h(3)
      a25sq = a(2,5)*a(2,5)
      h3sq = h(3)*h(3)
      a(3,5) = - a25sq - c1_12*h3sq
      a(4,5) = a(2,5)*(a25sq + c1_4*h3sq)
      a(5,5) = - a25sq*(a25sq + c1_2*h3sq) - c1_80*h3sq*h3sq
      a(6,5) = a(2,5)*(a25sq + c3_4*h3sq)*(a25sq + c1_12*h3sq)

      a(1,6) = - c1
      a(2,6) = - h(2) - h(3) - c1_2*h(4)
      a26sq = a(2,6)*a(2,6)
      h4sq = h(4)*h(4)
      a(3,6) = - a26sq - c1_12*h4sq
      a(4,6) = a(2,6)*(a26sq + c1_4*h4sq)
      a(5,6) = - a26sq*(a26sq + c1_2*h4sq) - c1_80*h4sq*h4sq
      a(6,6) = a(2,6)*(a26sq + c3_4*h4sq)*(a26sq + c1_12*h4sq)

      call edge_ih6_slope_ih5_coeff_common(a, tdecoeff, tdscoeff)

   end subroutine edge_ih6_slope_ih5_coeff_asymleft

   pure subroutine edge_ih6_slope_ih5_coeff_sym(h, tdecoeff, tdscoeff)
   ! ---------------------------------------------------------------------------
   ! With a symmetrical stencil, compute row coefficients for the tridiagonal
   ! system of equations to be solved for 6th and 5th order accurate edge and
   ! slope estimates, respectively.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:), intent(in) :: h
      real(r8), dimension(:), intent(inout) :: tdecoeff, tdscoeff

      real(r8), dimension(6,6) :: a
      real(r8) :: a23sq, a26sq, h1sq, h4sq

      ! Define matrix for linear system to be solved for edge coefficients.

      a(1,1) = c1
      a(2,1) = - h(2)
      a(3,1) = - a(2,1)*h(2)
      a(4,1) = - a(3,1)*h(2)
      a(5,1) = - a(4,1)*h(2)
      a(6,1) = - a(5,1)*h(2)

      a(1,2) = c1
      a(2,2) = h(3)
      a(3,2) = a(2,2)*h(3)
      a(4,2) = a(3,2)*h(3)
      a(5,2) = a(4,2)*h(3)
      a(6,2) = a(5,2)*h(3)

      a(1,3) = - c1
      a(2,3) = c1_2*h(1) + h(2)
      a23sq = a(2,3)*a(2,3)
      h1sq = h(1)*h(1)
      a(3,3) = - a23sq - c1_12*h1sq
      a(4,3) = a(2,3)*(a23sq + c1_4*h1sq)
      a(5,3) = - a23sq*(a23sq + c1_2*h1sq) - c1_80*h1sq*h1sq
      a(6,3) = a(2,3)*(a23sq + c3_4*h1sq)*(a23sq + c1_12*h1sq)

      a(1,4) = - c1
      a(2,4) = - c1_2*a(2,1)
      a(3,4) = - c1_3*a(3,1)
      a(4,4) = - c1_4*a(4,1)
      a(5,4) = - c1_5*a(5,1)
      a(6,4) = - c1_6*a(6,1)

      a(1,5) = - c1
      a(2,5) = - c1_2*a(2,2)
      a(3,5) = - c1_3*a(3,2)
      a(4,5) = - c1_4*a(4,2)
      a(5,5) = - c1_5*a(5,2)
      a(6,5) = - c1_6*a(6,2)

      a(1,6) = - c1
      a(2,6) = - h(3) - c1_2*h(4)
      a26sq = a(2,6)*a(2,6)
      h4sq = h(4)*h(4)
      a(3,6) = - a26sq - c1_12*h4sq
      a(4,6) = a(2,6)*(a26sq + c1_4*h4sq)
      a(5,6) = - a26sq*(a26sq + c1_2*h4sq) - c1_80*h4sq*h4sq
      a(6,6) = a(2,6)*(a26sq + c3_4*h4sq)*(a26sq + c1_12*h4sq)

      call edge_ih6_slope_ih5_coeff_common(a, tdecoeff, tdscoeff)

   end subroutine edge_ih6_slope_ih5_coeff_sym

   pure subroutine edge_ih6_slope_ih5_coeff_asymright(h, tdecoeff, tdscoeff)
   ! ---------------------------------------------------------------------------
   ! With an asymmetrical stencil, where edge values are shifted left compared
   ! to cell mean values, compute row coefficients for the tridiagonal system of
   ! equations to be solved for 6th and 5th order accurate edge and slope
   ! estimates, respectively.
   ! ---------------------------------------------------------------------------

      real(r8), dimension(:), intent(in) :: h
      real(r8), dimension(:), intent(inout) :: tdecoeff, tdscoeff

      real(r8), dimension(6,6) :: a
      real(r8) :: a23sq, a24sq, h1sq, h2sq

      ! Define matrix for linear system to be solved for edge coefficients.

      a(1,1) = c1
      a(2,1) = - h(3)
      a(3,1) = - a(2,1)*h(3)
      a(4,1) = - a(3,1)*h(3)
      a(5,1) = - a(4,1)*h(3)
      a(6,1) = - a(5,1)*h(3)

      a(1,2) = c1
      a(2,2) = h(4)
      a(3,2) = a(2,2)*h(4)
      a(4,2) = a(3,2)*h(4)
      a(5,2) = a(4,2)*h(4)
      a(6,2) = a(5,2)*h(4)

      a(1,3) = - c1
      a(2,3) = c1_2*h(1) + h(2) + h(3)
      a23sq = a(2,3)*a(2,3)
      h1sq = h(1)*h(1)
      a(3,3) = - a23sq - c1_12*h1sq
      a(4,3) = a(2,3)*(a23sq + c1_4*h1sq)
      a(5,3) = - a23sq*(a23sq + c1_2*h1sq) - c1_80*h1sq*h1sq
      a(6,3) = a(2,3)*(a23sq + c3_4*h1sq)*(a23sq + c1_12*h1sq)

      a(1,4) = - c1
      a(2,4) = c1_2*h(2) + h(3)
      a24sq = a(2,4)*a(2,4)
      h2sq = h(2)*h(2)
      a(3,4) = - a24sq - c1_12*h2sq
      a(4,4) = a(2,4)*(a24sq + c1_4*h2sq)
      a(5,4) = - a24sq*(a24sq + c1_2*h2sq) - c1_80*h2sq*h2sq
      a(6,4) = a(2,4)*(a24sq + c3_4*h2sq)*(a24sq + c1_12*h2sq)

      a(1,5) = - c1
      a(2,5) = - c1_2*a(2,1)
      a(3,5) = - c1_3*a(3,1)
      a(4,5) = - c1_4*a(4,1)
      a(5,5) = - c1_5*a(5,1)
      a(6,5) = - c1_6*a(6,1)

      a(1,6) = - c1
      a(2,6) = - c1_2*a(2,2)
      a(3,6) = - c1_3*a(3,2)
      a(4,6) = - c1_4*a(4,2)
      a(5,6) = - c1_5*a(5,2)
      a(6,6) = - c1_6*a(6,2)

      call edge_ih6_slope_ih5_coeff_common(a, tdecoeff, tdscoeff)

   end subroutine edge_ih6_slope_ih5_coeff_asymright

   pure subroutine edge_slope_lblu(lb_ord, h, a)
   ! ---------------------------------------------------------------------------
   ! Compute LU matrix for explicitly estimating lb_ord and lb_ord - 1 order
   ! accurate left edge and slope values, respectively.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: lb_ord
      real(r8), dimension(:), intent(in) :: h
      real(r8), dimension(:,:), intent(inout) :: a

      real(r8), dimension(lb_ord) :: a2sq, hsq
      integer :: i

      ! Define matrix for linear system to be solved for edge and slope values.

      a(1:lb_ord,1) = c1

      a(1,2) = c1_2*h(1)
      do i = 2, lb_ord
         a(i,2) = a(i-1,2) + c1_2*(h(i-1) + h(i))
      enddo

      if (lb_ord > 2) then

         a(1,3) = c1_3*a(1,2)*h(1)
         do i = 2, lb_ord
            a2sq(i) = a(i,2)*a(i,2)
            hsq(i) = h(i)*h(i)
            a(i,3) = c1_2*(a2sq(i) + c1_12*hsq(i))
         enddo

         if (lb_ord > 3) then

            a(1,4) = c1_4*a(1,3)*h(1)
            do i = 2, lb_ord
               a(i,4) = c1_6*a(i,2)*(a2sq(i) + c1_4*hsq(i))
            enddo

            if (lb_ord > 4) then

               a(1,5) = c1_5*a(1,4)*h(1)
               do i = 2, lb_ord
                  a(i,5) = c1_24*( a2sq(i)*(a2sq(i) + c1_2*hsq(i)) &
                                 + c1_80*hsq(i)*hsq(i))
               enddo

               if (lb_ord > 5) then

                  a(1,6) = c1_6*a(1,5)*h(1)
                  do i = 2, lb_ord
                     a(i,6) = c1_120*a(i,2)*(a2sq(i) + c3_4 *hsq(i)) &
                                           *(a2sq(i) + c1_12*hsq(i))
                  enddo

               endif
            endif
         endif
      endif

      ! LU decomposition.
      call lu_decompose(lb_ord, a)

   end subroutine edge_slope_lblu

   pure subroutine edge_slope_rblu(rb_ord, h, a)
   ! ---------------------------------------------------------------------------
   ! Compute LU matrix for explicitly estimating rb_ord and rb_ord - 1 order
   ! accurate right edge and slope values, respectively.
   ! ---------------------------------------------------------------------------

      integer, intent(in) :: rb_ord
      real(r8), dimension(:), intent(in) :: h
      real(r8), dimension(:,:), intent(inout) :: a

      real(r8), dimension(rb_ord) :: a2sq, hsq
      integer :: i

      ! Define matrix for linear system to be solved for edge and slope values.

      a(1:rb_ord,1) = c1

      a(rb_ord,2) = - c1_2*h(rb_ord)
      do i = rb_ord-1, 1, -1
         a(i,2) = a(i+1,2) - c1_2*(h(i+1) + h(i))
      enddo

      if (rb_ord > 2) then

         do i = 1, rb_ord-1
            a2sq(i) = a(i,2)*a(i,2)
            hsq(i) = h(i)*h(i)
            a(i,3) = c1_2*(a2sq(i) + c1_12*hsq(i))
         enddo
         a(rb_ord,3) = - c1_3*a(rb_ord,2)*h(rb_ord)

         if (rb_ord > 3) then

            do i = 1, rb_ord-1
               a(i,4) = c1_6*a(i,2)*(a2sq(i) + c1_4*hsq(i))
            enddo
            a(rb_ord,4) = - c1_4*a(rb_ord,3)*h(rb_ord)

            if (rb_ord > 4) then

               do i = 1, rb_ord-1
                  a(i,5) = c1_24*( a2sq(i)*(a2sq(i) + c1_2*hsq(i)) &
                                 + c1_80*hsq(i)*hsq(i))
               enddo
               a(rb_ord,5) = - c1_5*a(rb_ord,4)*h(rb_ord)

               if (rb_ord > 5) then

                  do i = 1, rb_ord-1
                     a(i,6) = c1_120*a(i,2)*(a2sq(i) + c3_4 *hsq(i)) &
                                           *(a2sq(i) + c1_12*hsq(i))
                  enddo
                  a(rb_ord,6) = - c1_6*a(rb_ord,5)*h(rb_ord)

               endif
            endif
         endif
      endif

      ! LU decomposition.
      call lu_decompose(rb_ord, a)

   end subroutine edge_slope_rblu

   pure subroutine prepare_pqm(rcgs, x_edge_src)
   ! ---------------------------------------------------------------------------
   ! Prepare reconstruction with piecewise quartics using implicit 6th order
   ! accurate edge and 5th order accurate slope estimation.
   ! ---------------------------------------------------------------------------

      type(recon_grd_struct), intent(inout) :: rcgs
      real(r8), dimension(:), intent(in) :: x_edge_src

      integer, dimension(rcgs%n_src) :: prev_index, next_index
      real(r8) :: hp, h_max, h_min, h
      integer :: ns, jp, j, last_index, jf, n, j_min, jn, lb_ord, rb_ord, jd, js
      integer :: first_index
!     integer :: first_index = 0 ! Initialized to avoid compiler warning.

      ! Exclude near-empty grid cells and establish a doubly linked list that
      ! connects the remaining grid cells.
      ns = 0
      jp = 0
      do j = 1, rcgs%n_src
         rcgs%h_src(j) = abs(x_edge_src(j+1) - x_edge_src(j))
         if (rcgs%h_src(j) > c2*rcgs%x_eps) then
            ns = ns + 1
            rcgs%src_dst_index(j) = 1
            prev_index(j) = jp
            if (jp == 0) then
               first_index = j
            else
               next_index(jp) = j
            endif
            jp = j
         else
            rcgs%src_dst_index(j) = 0
         endif
      enddo
      last_index = jp
      next_index(jp) = 0
      if (ns < n_src_min_pqm) then
         rcgs%n_src_actual = ns
         return
      endif

      ! Exclude grid cells that may lead to large condition numbers for the
      ! linear systems to be solved in edge_ih6_slope_ih5_coeff_asymleft,
      ! edge_ih6_slope_ih5_coeff_sym and edge_ih6_slope_ih5_coeff_asymright.
      ! Excluded grid cells are merged with the non-excluded neighbour grid cell
      ! having the smallest grid cell width.
      jf = first_index
      outer: do
         j = jf
         hp = rcgs%h_src(j)
         h_max = rcgs%h_src(j)
         do n = 1, 3
            j = next_index(j)
            if (j == 0) exit outer
            hp = hp*rcgs%h_src(j)
            h_max = max(h_max, rcgs%h_src(j))
         enddo
         if (hp > hplim_ih6*h_max**4) then
            jf = next_index(jf)
         else
            ns = ns - 1
            if (ns < n_src_min_pqm) then
               rcgs%n_src_actual = ns
               return
            endif
            j = jf
            h_min = rcgs%h_src(j)
            j_min = j
            do n = 1, 3
               j = next_index(j)
               if (rcgs%h_src(j) < h_min) then
                  h_min = rcgs%h_src(j)
                  j_min = j
               endif
            enddo
            jp = prev_index(j_min)
            jn = next_index(j_min)
            if     (jp == 0) then
               rcgs%src_dst_index(j_min) = - jn
               rcgs%h_src(jn) = rcgs%h_src(jn) + rcgs%h_src(j_min)
               first_index = jn
               prev_index(jn) = 0
               jf = jn
            elseif (jn == 0) then
               rcgs%src_dst_index(j_min) = - jp
               rcgs%h_src(jp) = rcgs%h_src(jp) + rcgs%h_src(j_min)
               next_index(jp) = 0
               last_index = jp
               exit
            else
               if (rcgs%h_src(jn) < rcgs%h_src(jp)) then
                  rcgs%src_dst_index(j_min) = - jn
                  rcgs%h_src(jn) = rcgs%h_src(jn) + rcgs%h_src(j_min)
               else
                  rcgs%src_dst_index(j_min) = - jp
                  rcgs%h_src(jp) = rcgs%h_src(jp) + rcgs%h_src(j_min)
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
      ! linear systems to be solved for edge and slope values. Excluded grid
      ! cells are merged with the non-excluded neighbour grid cell having the
      ! smallest grid cell width.
      lb_ord = min(ns, rcgs%left_bndr_ord)
      call left_bndr_cond(rcgs, prev_index, next_index, first_index, &
                          last_index, lb_ord, ns, n_src_min_pqm)
      if (ns < n_src_min_pqm) then
         rcgs%n_src_actual = ns
         return
      endif
      rb_ord = min(ns, rcgs%right_bndr_ord)
      call right_bndr_cond(rcgs, prev_index, next_index, last_index, &
                           rb_ord, ns, n_src_min_pqm)
      if (ns < n_src_min_pqm) then
         rcgs%n_src_actual = ns
         return
      endif

      ! For the non-excluded grid cells, assign the destination index in the
      ! continuous array of grid cells to be used in the reconstruction. Also
      ! set the grid cell widths of the continuous array.
      jd = 0
      do js = 1, rcgs%n_src
         if (rcgs%src_dst_index(js) > 0) then
            jd = jd + 1
            rcgs%src_dst_index(js) = jd
            rcgs%h_src(jd) = rcgs%h_src(js)
            rcgs%hi_src(jd) = c1/rcgs%h_src(jd)
         endif
      enddo

      ! Find the destination index of excluded grid cells to be merged and
      ! compute the mapping weights.
      do js = 1, rcgs%n_src
         jd = rcgs%src_dst_index(js)
         do while (jd < 0)
            jd = rcgs%src_dst_index(- jd)
         enddo
         rcgs%src_dst_index(js) = jd
         if (jd > 0) then
            h = abs(x_edge_src(js+1) - x_edge_src(js))
            if (abs(h - rcgs%h_src(jd)) < rcgs%x_eps) then
               rcgs%src_dst_weight(js) = c1
            else
               rcgs%src_dst_weight(js) = h*rcgs%hi_src(jd)
            endif
         endif
      enddo

      ! Set source edge values in the continuous reconstruction array.
      rcgs%x_edge_src(1) = x_edge_src(1)
      js = 1
      do j = 1, ns-1
         do
            js = js + 1
            if (rcgs%src_dst_index(js) /= j .and. &
                rcgs%src_dst_index(js) /= 0) exit
         enddo
         rcgs%x_edge_src(j+1) = x_edge_src(js)
      enddo
      rcgs%x_edge_src(ns+1) = x_edge_src(rcgs%n_src+1)

      ! Compute the multiplicative inverse of cell width used for estimating
      ! centered linear slope.
      do j = 2, ns-1
         rcgs%hci_src(j) = c2/( rcgs%h_src(j-1) + c2*rcgs%h_src(j) &
                              + rcgs%h_src(j+1))
      enddo

      ! Compute coefficients for the tridiagonal system of equations for the
      ! estimation of interior edge and slope values. If the row coefficients
      ! for 6th order accurate edge estimates or 5th order accurate slope
      ! estimates leads to a tridiagonal matrix that is not diagonally dominant,
      ! use instead row coefficients for 4th and 3rd order accurate edge and
      ! slope estimates, respectively.
      if (lb_ord < 5) then
         call edge_ih4_coeff(rcgs%h_src(1:2), rcgs%tdecoeff(:,2))
         rcgs%tdecoeff(5,2) = c0
         rcgs%tdecoeff(6,2) = c0
         call slope_ih3_coeff(rcgs%h_src(1:2), rcgs%tdscoeff(:,2))
         rcgs%tdscoeff(5,2) = c0
         rcgs%tdscoeff(6,2) = c0
      else
         call edge_ih6_slope_ih5_coeff_asymleft(rcgs%h_src(1:4), &
                                                rcgs%tdecoeff(:,2), &
                                                rcgs%tdscoeff(:,2))
         if (abs(rcgs%tdecoeff(1,2)) + abs(rcgs%tdecoeff(2,2)) > c1 .or. &
             abs(rcgs%tdscoeff(1,2)) + abs(rcgs%tdscoeff(2,2)) > c1) then
            call edge_ih4_coeff(rcgs%h_src(1:2), rcgs%tdecoeff(:,2))
            rcgs%tdecoeff(5,2) = c0
            rcgs%tdecoeff(6,2) = c0
            call slope_ih3_coeff(rcgs%h_src(1:2), rcgs%tdscoeff(:,2))
            rcgs%tdscoeff(5,2) = c0
            rcgs%tdscoeff(6,2) = c0
         endif
      endif
      do j = 3, ns-1
         call edge_ih6_slope_ih5_coeff_sym(rcgs%h_src((j-2):(j+1)), &
                                           rcgs%tdecoeff(:,j), &
                                           rcgs%tdscoeff(:,j))
         if (abs(rcgs%tdecoeff(1,j)) + abs(rcgs%tdecoeff(2,j)) > c1 .or. &
             abs(rcgs%tdscoeff(1,j)) + abs(rcgs%tdscoeff(2,j)) > c1) then
            call edge_ih4_coeff(rcgs%h_src((j-1):j), rcgs%tdecoeff(:,j))
            rcgs%tdecoeff(5,j) = rcgs%tdecoeff(4,j)
            rcgs%tdecoeff(4,j) = rcgs%tdecoeff(3,j)
            rcgs%tdecoeff(3,j) = c0
            rcgs%tdecoeff(6,j) = c0
            call slope_ih3_coeff(rcgs%h_src((j-1):j), rcgs%tdscoeff(:,j))
            rcgs%tdscoeff(5,j) = rcgs%tdscoeff(4,j)
            rcgs%tdscoeff(4,j) = rcgs%tdscoeff(3,j)
            rcgs%tdscoeff(3,j) = c0
            rcgs%tdscoeff(6,j) = c0
         endif
      enddo
      if (rb_ord < 5) then
         call edge_ih4_coeff(rcgs%h_src((ns-1):ns), rcgs%tdecoeff(:,ns))
         rcgs%tdecoeff(5,ns) = rcgs%tdecoeff(3,ns)
         rcgs%tdecoeff(6,ns) = rcgs%tdecoeff(4,ns)
         rcgs%tdecoeff(3,ns) = c0
         rcgs%tdecoeff(4,ns) = c0
         call slope_ih3_coeff(rcgs%h_src((ns-1):ns), rcgs%tdscoeff(:,ns))
         rcgs%tdscoeff(5,ns) = rcgs%tdscoeff(3,ns)
         rcgs%tdscoeff(6,ns) = rcgs%tdscoeff(4,ns)
         rcgs%tdscoeff(3,ns) = c0
         rcgs%tdscoeff(4,ns) = c0
      else
         call edge_ih6_slope_ih5_coeff_asymright(rcgs%h_src((ns-3):ns), &
                                                 rcgs%tdecoeff(:,ns), &
                                                 rcgs%tdscoeff(:,ns))
         if (abs(rcgs%tdecoeff(1,ns)) + abs(rcgs%tdecoeff(2,ns)) > c1 .or. &
             abs(rcgs%tdscoeff(1,ns)) + abs(rcgs%tdscoeff(2,ns)) > c1) then
            call edge_ih4_coeff(rcgs%h_src((ns-1):ns), rcgs%tdecoeff(:,ns))
            rcgs%tdecoeff(5,ns) = rcgs%tdecoeff(3,ns)
            rcgs%tdecoeff(6,ns) = rcgs%tdecoeff(4,ns)
            rcgs%tdecoeff(3,ns) = c0
            rcgs%tdecoeff(4,ns) = c0
            call slope_ih3_coeff(rcgs%h_src((ns-1):ns), rcgs%tdscoeff(:,ns))
            rcgs%tdscoeff(5,ns) = rcgs%tdscoeff(3,ns)
            rcgs%tdscoeff(6,ns) = rcgs%tdscoeff(4,ns)
            rcgs%tdscoeff(3,ns) = c0
            rcgs%tdscoeff(4,ns) = c0
         endif
      endif

      ! Compute LU matrices for the explicit estimation of boundary edge and
      ! slope values.
      if (lb_ord > 1) &
         call edge_slope_lblu(lb_ord, rcgs%h_src(1:lb_ord), rcgs%lblu)
      if (rb_ord > 1) &
         call edge_slope_rblu(rb_ord, rcgs%h_src((ns-rb_ord+1):ns), rcgs%rblu)

      rcgs%n_src_actual = ns
      rcgs%left_bndr_ord_actual = lb_ord
      rcgs%right_bndr_ord_actual = rb_ord

   end subroutine prepare_pqm

   pure subroutine prepare_ppm(rcgs, x_edge_src)
   ! ---------------------------------------------------------------------------
   ! Prepare reconstruction with piecewise parabolas using implicit 4th order
   ! accurate edge estimation.
   ! ---------------------------------------------------------------------------

      type(recon_grd_struct), intent(inout) :: rcgs
      real(r8), dimension(:), intent(in) :: x_edge_src

      integer, dimension(rcgs%n_src) :: prev_index, next_index
      real(r8) :: h
      integer :: ns, jp, j, last_index, jf, jl, lb_ord, rb_ord, jd, js
      integer :: first_index
!     integer :: first_index = 0 ! Initialized to avoid compiler warning.

      ! Exclude near-empty grid cells and establish a doubly linked list that
      ! connects the remaining grid cells.
      ns = 0
      jp = 0
      do j = 1, rcgs%n_src
         rcgs%h_src(j) = abs(x_edge_src(j+1) - x_edge_src(j))
         if (rcgs%h_src(j) > c2*rcgs%x_eps) then
            ns = ns + 1
            rcgs%src_dst_index(j) = 1
            prev_index(j) = jp
            if (jp == 0) then
               first_index = j
            else
               next_index(jp) = j
            endif
            jp = j
         else
            rcgs%src_dst_index(j) = 0
         endif
      enddo
      last_index = jp
      next_index(jp) = 0
      if (ns < n_src_min_ppm) then
         rcgs%n_src_actual = ns
         return
      endif

      ! Exclude grid cells that may lead to large condition numbers for the
      ! linear systems to be solved in edge_ih4_coeff. Excluded grid cells are
      ! merged with the non-excluded neighbour grid cell having the smallest
      ! grid cell width.
      jf = first_index
      jl = next_index(jf)
      do
         if (rcgs%h_src(jf)*rcgs%h_src(jl) > &
             hplim_ih4*max(rcgs%h_src(jf), rcgs%h_src(jl))**2) then
            jf = jl
            jl = next_index(jf)
            if (jl == 0) exit
         else
            ns = ns - 1
            if (ns < n_src_min_ppm) then
               rcgs%n_src_actual = ns
               return
            endif
            if (rcgs%h_src(jf) < rcgs%h_src(jl)) then
               j = jf
               jf = prev_index(jf)
               prev_index(jl) = jf
               if (jf == 0) then
                  rcgs%src_dst_index(j) = - jl
                  rcgs%h_src(jl) = rcgs%h_src(jl) + rcgs%h_src(j)
                  first_index = jl
                  jf = jl
                  jl = next_index(jf)
                  if (jl == 0) exit
               else
                  if (rcgs%h_src(jf) < rcgs%h_src(jl)) then
                     rcgs%src_dst_index(j) = - jf
                     rcgs%h_src(jf) = rcgs%h_src(jf) + rcgs%h_src(j)
                  else
                     rcgs%src_dst_index(j) = - jl
                     rcgs%h_src(jl) = rcgs%h_src(jl) + rcgs%h_src(j)
                  endif
                  next_index(jf) = jl
               endif
            else
               j = jl
               jl = next_index(jl)
               next_index(jf) = jl
               if (jl == 0) then
                  rcgs%src_dst_index(j) = - jf
                  rcgs%h_src(jf) = rcgs%h_src(jf) + rcgs%h_src(j)
                  last_index = jf
                  exit
               endif
               if (rcgs%h_src(jf) < rcgs%h_src(jl)) then
                  rcgs%src_dst_index(j) = - jf
                  rcgs%h_src(jf) = rcgs%h_src(jf) + rcgs%h_src(j)
               else
                  rcgs%src_dst_index(j) = - jl
                  rcgs%h_src(jl) = rcgs%h_src(jl) + rcgs%h_src(j)
               endif
               prev_index(jl) = jf
            endif
         endif
      enddo

      ! Exclude grid cells that may lead to a large condition number for the
      ! linear systems to be solved for edge and slope values. Excluded grid
      ! cells are merged with the non-excluded neighbour grid cell having the
      ! smallest grid cell width.
      lb_ord = min(ns, rcgs%left_bndr_ord, eb_ord_max_ppm)
      call left_bndr_cond(rcgs, prev_index, next_index, first_index, &
                          last_index, lb_ord, ns, n_src_min_ppm)
      if (ns < n_src_min_ppm) then
         rcgs%n_src_actual = ns
         return
      endif
      rb_ord = min(ns, rcgs%right_bndr_ord, eb_ord_max_ppm)
      call right_bndr_cond(rcgs, prev_index, next_index, last_index, &
                           rb_ord, ns, n_src_min_ppm)
      if (ns < n_src_min_ppm) then
         rcgs%n_src_actual = ns
         return
      endif

      ! For the non-excluded grid cells, assign the destination index in the
      ! continuous array of grid cells to be used in the reconstruction. Also
      ! set the grid cell widths of the continuous array.
      jd = 0
      do js = 1, rcgs%n_src
         if (rcgs%src_dst_index(js) > 0) then
            jd = jd + 1
            rcgs%src_dst_index(js) = jd
            rcgs%h_src(jd) = rcgs%h_src(js)
            rcgs%hi_src(jd) = c1/rcgs%h_src(jd)
         endif
      enddo

      ! Find the destination index of excluded grid cells to be merged and
      ! compute the mapping weights.
      do js = 1, rcgs%n_src
         jd = rcgs%src_dst_index(js)
         do while (jd < 0)
            jd = rcgs%src_dst_index(- jd)
         enddo
         rcgs%src_dst_index(js) = jd
         if (jd > 0) then
            h = abs(x_edge_src(js+1) - x_edge_src(js))
            if (abs(h - rcgs%h_src(jd)) < rcgs%x_eps) then
               rcgs%src_dst_weight(js) = c1
            else
               rcgs%src_dst_weight(js) = h*rcgs%hi_src(jd)
            endif
         endif
      enddo

      ! Set source edge values in the continuous reconstruction array.
      rcgs%x_edge_src(1) = x_edge_src(1)
      js = 1
      do j = 1, ns-1
         do
            js = js + 1
            if (rcgs%src_dst_index(js) /= j .and. &
                rcgs%src_dst_index(js) /= 0) exit
         enddo
         rcgs%x_edge_src(j+1) = x_edge_src(js)
      enddo
      rcgs%x_edge_src(ns+1) = x_edge_src(rcgs%n_src+1)

      ! Compute the multiplicative inverse of cell width used for estimating
      ! centered linear slope.
      do j = 2, ns-1
         rcgs%hci_src(j) = c2/( rcgs%h_src(j-1) + c2*rcgs%h_src(j) &
                              + rcgs%h_src(j+1))
      enddo

      ! Compute coefficients for the tridiagonal system of equations for the
      ! estimation of interior edge values.
      do j = 2, ns
         call edge_ih4_coeff(rcgs%h_src((j-1):j), rcgs%tdecoeff(:,j))
      enddo

      ! Compute LU matrices for the explicit estimation of boundary edge values.
      if (lb_ord > 1) &
         call edge_slope_lblu(lb_ord, rcgs%h_src(1:lb_ord), rcgs%lblu)
      if (rb_ord > 1) &
         call edge_slope_rblu(rb_ord, rcgs%h_src((ns-rb_ord+1):ns), rcgs%rblu)

      rcgs%n_src_actual = ns
      rcgs%left_bndr_ord_actual = lb_ord
      rcgs%right_bndr_ord_actual = rb_ord

   end subroutine prepare_ppm

   pure subroutine prepare_plm(rcgs, x_edge_src)
   ! ---------------------------------------------------------------------------
   ! Prepare reconstruction with piecewise lines.
   ! ---------------------------------------------------------------------------

      type(recon_grd_struct), intent(inout) :: rcgs
      real(r8), dimension(:), intent(in) :: x_edge_src

      integer :: ns, j, js

      ! Exclude near-empty grid cells and assign the destination index in the
      ! continuous array of grid cells to be used in the reconstruction.
      ns = 0
      do j = 1, rcgs%n_src
         if (abs(x_edge_src(j+1) - x_edge_src(j)) > c2*rcgs%x_eps) then
            ns = ns + 1
            rcgs%src_dst_index(j) = ns
         else
            rcgs%src_dst_index(j) = 0
         endif
      enddo
      if (ns < n_src_min_plm) then
         rcgs%n_src_actual = ns
         return
      endif

      ! Set source edge values in the continuous reconstruction array.
      rcgs%x_edge_src(1) = x_edge_src(1)
      js = 1
      do j = 1, ns-1
         do
            js = js + 1
            if (rcgs%src_dst_index(js) /= j .and. &
                rcgs%src_dst_index(js) /= 0) exit
         enddo
         rcgs%x_edge_src(j+1) = x_edge_src(js)
      enddo
      rcgs%x_edge_src(ns+1) = x_edge_src(rcgs%n_src+1)

      ! From edge locations, obtain source grid cell widths and their
      ! multiplicative inverse.
      do j = 1, ns
         rcgs%h_src(j) = abs(rcgs%x_edge_src(j+1) - rcgs%x_edge_src(j))
         rcgs%hi_src(j) = c1/rcgs%h_src(j)
      enddo

      ! Compute the multiplicative inverse of cell width used for estimating
      ! centered linear slope.
      do j = 2, ns-1
         rcgs%hci_src(j) = c2/( rcgs%h_src(j-1) + c2*rcgs%h_src(j) &
                              + rcgs%h_src(j+1))
      enddo

      rcgs%n_src_actual = ns

   end subroutine prepare_plm

   pure subroutine prepare_pcm(rcgs, x_edge_src)
   ! ---------------------------------------------------------------------------
   ! Prepare piecewise constant reconstruction.
   ! ---------------------------------------------------------------------------

      type(recon_grd_struct), intent(inout) :: rcgs
      real(r8), dimension(:), intent(in) :: x_edge_src

      integer :: ns, j, js

      ! Exclude near-empty grid cells and assign the destination index in the
      ! continuous array of grid cells to be used in the reconstruction.
      ns = 0
      do j = 1, rcgs%n_src
         if (abs(x_edge_src(j+1) - x_edge_src(j)) > c2*rcgs%x_eps) then
            ns = ns + 1
            rcgs%src_dst_index(j) = ns
         else
            rcgs%src_dst_index(j) = 0
         endif
      enddo
      if (ns == 0) then
         rcgs%n_src_actual = ns
         return
      endif

      ! Set source edge values in the continuous reconstruction array.
      rcgs%x_edge_src(1) = x_edge_src(1)
      js = 1
      do j = 1, ns-1
         do
            js = js + 1
            if (rcgs%src_dst_index(js) /= j .and. &
                rcgs%src_dst_index(js) /= 0) exit
         enddo
         rcgs%x_edge_src(j+1) = x_edge_src(js)
      enddo
      rcgs%x_edge_src(ns+1) = x_edge_src(rcgs%n_src+1)

      ! From edge locations, obtain source grid cell widths and their
      ! multiplicative inverse.
      do j = 1, ns
         rcgs%h_src(j) = abs(rcgs%x_edge_src(j+1) - rcgs%x_edge_src(j))
         rcgs%hi_src(j) = c1/rcgs%h_src(j)
      enddo

      rcgs%n_src_actual = ns

   end subroutine prepare_pcm

   pure subroutine reconstruct_plm_no_limiting(rcss)
   ! ---------------------------------------------------------------------------
   ! Carry out a reconstruction with piecewise lines.
   ! ---------------------------------------------------------------------------

      type(recon_src_struct), intent(inout) :: rcss

      real(r8) :: sc
      integer :: ns, j

      ns = rcss%rcgs%n_src_actual

      sc = c2*(rcss%u_src(2) - rcss%u_src(1)) &
             /(rcss%rcgs%h_src(2) + rcss%rcgs%h_src(1))
      rcss%polycoeff(2,1) = sc*rcss%rcgs%h_src(1)
      rcss%polycoeff(1,1) = rcss%u_src(1) - c1_2*rcss%polycoeff(2,1)
      rcss%uel(1) = rcss%polycoeff(1,1)
      rcss%uer(1) = rcss%polycoeff(1,1) + rcss%polycoeff(2,1)
      do j = 2, ns-1
         sc = (rcss%u_src(j+1) - rcss%u_src(j-1))*rcss%rcgs%hci_src(j)
         rcss%polycoeff(2,j) = sc*rcss%rcgs%h_src(j)
         rcss%polycoeff(1,j) = rcss%u_src(j) - c1_2*rcss%polycoeff(2,j)
         rcss%uel(j) = rcss%polycoeff(1,j)
         rcss%uer(j) = rcss%polycoeff(1,j) + rcss%polycoeff(2,j)
      enddo
      sc = c2*(rcss%u_src(ns) - rcss%u_src(ns-1)) &
             /(rcss%rcgs%h_src(ns) + rcss%rcgs%h_src(ns-1))
      rcss%polycoeff(2,ns) = sc*rcss%rcgs%h_src(ns)
      rcss%polycoeff(1,ns) = rcss%u_src(ns) - c1_2*rcss%polycoeff(2,ns)
      rcss%uel(ns) = rcss%polycoeff(1,ns)
      rcss%uer(ns) = rcss%polycoeff(1,ns) + rcss%polycoeff(2,ns)

      rcss%polycoeff(3:rcss%rcgs%p_ord+1,:) = c0

   end subroutine reconstruct_plm_no_limiting

   pure subroutine reconstruct_plm_monotonic(rcss)
   ! ---------------------------------------------------------------------------
   ! Carry out a reconstruction with piecewise lines and apply limiting to
   ! ensure a monotonic reconstruction.
   ! ---------------------------------------------------------------------------

      type(recon_src_struct), intent(inout) :: rcss

      real(r8) :: sl, sr, sc
      integer :: ns, j

      ns = rcss%rcgs%n_src_actual

      ! Use monotonized central-difference limiter for interior grid cells.
      do j = 2, ns-1
         sl = c2*(rcss%u_src(j) - rcss%u_src(j-1))*rcss%rcgs%hi_src(j)
         sr = c2*(rcss%u_src(j+1) - rcss%u_src(j))*rcss%rcgs%hi_src(j)
         if (sl*sr > c0) then
            sc = (rcss%u_src(j+1) - rcss%u_src(j-1))*rcss%rcgs%hci_src(j)
            sc = sign(min(abs(sl), abs(sr), abs(sc)), sc)
         else
            sc = c0
         endif
         rcss%polycoeff(2,j) = sc*rcss%rcgs%h_src(j)
         rcss%polycoeff(1,j) = rcss%u_src(j) - c1_2*rcss%polycoeff(2,j)
         rcss%uel(j) = rcss%polycoeff(1,j)
         rcss%uer(j) = rcss%polycoeff(1,j) + rcss%polycoeff(2,j)
      enddo

      if (rcss%pc_left_bndr) then
         ! Piecewise constant reconstruction of left boundary cell.
         rcss%polycoeff(1,1) = rcss%u_src(1)
         rcss%polycoeff(2,1) = c0
         rcss%uel(1) = rcss%u_src(1)
         rcss%uer(1) = rcss%u_src(1)
      else
         ! Piecewise linear reconstruction of left boundary cell.
         sc = c2*(rcss%u_src(2) - rcss%u_src(1)) &
                /(rcss%rcgs%h_src(2) + rcss%rcgs%h_src(1))
         rcss%polycoeff(2,1) = sc*rcss%rcgs%h_src(1)
         rcss%polycoeff(1,1) = rcss%u_src(1) - c1_2*rcss%polycoeff(2,1)
         rcss%uel(1) = rcss%polycoeff(1,1)
         rcss%uer(1) = rcss%polycoeff(1,1) + rcss%polycoeff(2,1)
      endif

      if (rcss%pc_right_bndr) then
         ! Piecewise constant reconstruction of right boundary cell.
         rcss%polycoeff(1,ns) = rcss%u_src(ns)
         rcss%polycoeff(2,ns) = c0
         rcss%uel(ns) = rcss%u_src(ns)
         rcss%uer(ns) = rcss%u_src(ns)
      else
         ! Piecewise linear reconstruction of right boundary cell.
         sc = c2*(rcss%u_src(ns) - rcss%u_src(ns-1)) &
                /(rcss%rcgs%h_src(ns) + rcss%rcgs%h_src(ns-1))
         rcss%polycoeff(2,ns) = sc*rcss%rcgs%h_src(ns)
         rcss%polycoeff(1,ns) = rcss%u_src(ns) - c1_2*rcss%polycoeff(2,ns)
         rcss%uel(ns) = rcss%polycoeff(1,ns)
         rcss%uer(ns) = rcss%polycoeff(1,ns) + rcss%polycoeff(2,ns)
      endif

      rcss%polycoeff(3:rcss%rcgs%p_ord+1,:) = c0

   end subroutine reconstruct_plm_monotonic

   pure subroutine reconstruct_ppm_edge_values(rcss)
   ! ---------------------------------------------------------------------------
   ! Reconstruct edge values using an implicit 4th order scheme.
   ! ---------------------------------------------------------------------------

      type(recon_src_struct), intent(inout) :: rcss

      real(r8), dimension(eb_ord_max_ppm) :: x
      real(r8), dimension(rcss%rcgs%n_src_actual+1) :: uedge
      real(r8), dimension(rcss%rcgs%n_src_actual) :: rhs, gam
      real(r8) :: bei
      integer :: ns, lb_ord, rb_ord, j

      ns = rcss%rcgs%n_src_actual
      lb_ord = rcss%rcgs%left_bndr_ord_actual
      rb_ord = rcss%rcgs%right_bndr_ord_actual

      ! Obtain the left boundary edge value.
      if (lb_ord == 1) then
         uedge(1) = rcss%u_src(1)
      else
         x(1:lb_ord) = rcss%u_src(1:lb_ord)
         call lu_solve(lb_ord, rcss%rcgs%lblu, x)
         uedge(1) = x(1)
      endif

      ! Obtain the right boundary edge value.
      if (lb_ord == 1) then
         uedge(ns+1) = rcss%u_src(ns)
      else
         x(1:rb_ord) = rcss%u_src((ns-rb_ord+1):ns)
         call lu_solve(rb_ord, rcss%rcgs%rblu, x)
         uedge(ns+1) = x(1)
      endif

      ! Obtain right hand side of tridiagonal system of equations.
      do j = 2, ns
         rhs(j) = rcss%rcgs%tdecoeff(3,j)*rcss%u_src(j-1) &
                + rcss%rcgs%tdecoeff(4,j)*rcss%u_src(j  )
      enddo

      ! Solve tridiagonal system of equations to obtain interior edge values.
      gam(1) = c0
      do j = 2, ns
         bei = c1/(c1 - rcss%rcgs%tdecoeff(1,j)*gam(j-1))
         uedge(j) = (rhs(j) - rcss%rcgs%tdecoeff(1,j)*uedge(j-1))*bei
         gam(j) = rcss%rcgs%tdecoeff(2,j)*bei
      enddo
      do j = ns, 2, -1
         uedge(j) = uedge(j) - gam(j)*uedge(j+1)
      enddo

      ! Set left and right edge values for each grid cell.
      rcss%uel(1:ns) = uedge(1:ns)
      rcss%uer(1:ns) = uedge(2:(ns+1))

   end subroutine reconstruct_ppm_edge_values

   pure subroutine reconstruct_pqm_edge_slope_values(rcss)
   ! ---------------------------------------------------------------------------
   ! Reconstruct edge and slope values using implicit 6th and 5th order schemes,
   ! respectively.
   ! ---------------------------------------------------------------------------

      type(recon_src_struct), intent(inout) :: rcss

      real(r8), dimension(eb_ord_max_pqm) :: x
      real(r8), dimension(rcss%rcgs%n_src_actual+1) :: uedge, uslope
      real(r8), dimension(rcss%rcgs%n_src_actual) :: rhs, gam
      real(r8) :: bei
      integer :: ns, lb_ord, rb_ord, j

      ns = rcss%rcgs%n_src_actual
      lb_ord = rcss%rcgs%left_bndr_ord_actual
      rb_ord = rcss%rcgs%right_bndr_ord_actual

      ! Obtain the left boundary edge and slope values.
      if (lb_ord == 1) then
         uedge(1) = rcss%u_src(1)
         uslope(1) = c0
      else
         x(1:lb_ord) = rcss%u_src(1:lb_ord)
         call lu_solve(lb_ord, rcss%rcgs%lblu, x)
         uedge(1) = x(1)
         uslope(1) = x(2)
      endif

      ! Obtain the right boundary edge and slope values.
      if (rb_ord == 1) then
         uedge(ns+1) = rcss%u_src(ns)
         uslope(ns+1) = c0
      else
         x(1:rb_ord) = rcss%u_src((ns-rb_ord+1):ns)
         call lu_solve(rb_ord, rcss%rcgs%rblu, x)
         uedge(ns+1) = x(1)
         uslope(ns+1) = x(2)
      endif

      ! Obtain right hand side of tridiagonal system of equations for edge
      ! values.
      rhs(2) = rcss%rcgs%tdecoeff(3,2)*rcss%u_src(1) &
             + rcss%rcgs%tdecoeff(4,2)*rcss%u_src(2) &
             + rcss%rcgs%tdecoeff(5,2)*rcss%u_src(3) &
             + rcss%rcgs%tdecoeff(6,2)*rcss%u_src(4)
      do j = 3, ns-1
         rhs(j) = rcss%rcgs%tdecoeff(3,j)*rcss%u_src(j-2) &
                + rcss%rcgs%tdecoeff(4,j)*rcss%u_src(j-1) &
                + rcss%rcgs%tdecoeff(5,j)*rcss%u_src(j  ) &
                + rcss%rcgs%tdecoeff(6,j)*rcss%u_src(j+1)
      enddo
      rhs(ns) = rcss%rcgs%tdecoeff(3,ns)*rcss%u_src(ns-3) &
              + rcss%rcgs%tdecoeff(4,ns)*rcss%u_src(ns-2) &
              + rcss%rcgs%tdecoeff(5,ns)*rcss%u_src(ns-1) &
              + rcss%rcgs%tdecoeff(6,ns)*rcss%u_src(ns  )

      ! Solve tridiagonal system of equations to obtain interior edge values.
      gam(1) = c0
      do j = 2, ns
         bei = c1/(c1 - rcss%rcgs%tdecoeff(1,j)*gam(j-1))
         uedge(j) = (rhs(j) - rcss%rcgs%tdecoeff(1,j)*uedge(j-1))*bei
         gam(j) = rcss%rcgs%tdecoeff(2,j)*bei
      enddo
      do j = ns, 2, -1
         uedge(j) = uedge(j) - gam(j)*uedge(j+1)
      enddo

      ! Obtain right hand side of tridiagonal system of equations for slope
      ! values.
      rhs(2) = rcss%rcgs%tdscoeff(3,2)*rcss%u_src(1) &
             + rcss%rcgs%tdscoeff(4,2)*rcss%u_src(2) &
             + rcss%rcgs%tdscoeff(5,2)*rcss%u_src(3) &
             + rcss%rcgs%tdscoeff(6,2)*rcss%u_src(4)
      do j = 3, ns-1
         rhs(j) = rcss%rcgs%tdscoeff(3,j)*rcss%u_src(j-2) &
                + rcss%rcgs%tdscoeff(4,j)*rcss%u_src(j-1) &
                + rcss%rcgs%tdscoeff(5,j)*rcss%u_src(j  ) &
                + rcss%rcgs%tdscoeff(6,j)*rcss%u_src(j+1)
      enddo
      rhs(ns) = rcss%rcgs%tdscoeff(3,ns)*rcss%u_src(ns-3) &
              + rcss%rcgs%tdscoeff(4,ns)*rcss%u_src(ns-2) &
              + rcss%rcgs%tdscoeff(5,ns)*rcss%u_src(ns-1) &
              + rcss%rcgs%tdscoeff(6,ns)*rcss%u_src(ns  )

      ! Solve tridiagonal system of equations to obtain interior slope values.
      gam(1) = c0
      do j = 2, ns
         bei = c1/(c1 - rcss%rcgs%tdscoeff(1,j)*gam(j-1))
         uslope(j) = (rhs(j) - rcss%rcgs%tdscoeff(1,j)*uslope(j-1))*bei
         gam(j) = rcss%rcgs%tdscoeff(2,j)*bei
      enddo
      do j = ns, 2, -1
         uslope(j) = uslope(j) - gam(j)*uslope(j+1)
      enddo

      ! Set left and right edge values for each grid cell.
      rcss%uel(1:ns) = uedge(1:ns)
      rcss%uer(1:ns) = uedge(2:(ns+1))

      ! Set left and right slope values for each grid cell and scale the slope
      ! values with the grid cell widths.
      rcss%usl(1:ns) = uslope(1:ns)*rcss%rcgs%h_src(1:ns)
      rcss%usr(1:ns) = uslope(2:(ns+1))*rcss%rcgs%h_src(1:ns)

   end subroutine reconstruct_pqm_edge_slope_values

   pure subroutine limit_ppm_interior_monotonic(rcss)
   ! ---------------------------------------------------------------------------
   ! Apply limiting to ensure a monotonic reconstruction of piecewise parabolas
   ! for interior grid cells.
   ! ---------------------------------------------------------------------------

      type(recon_src_struct), intent(inout) :: rcss

      real(r8) :: sl, sr, sc, d, q, r
      integer :: ns, j

      ns = rcss%rcgs%n_src_actual

      do j = 2, ns-1
         sl = c2*(rcss%u_src(j) - rcss%u_src(j-1))*rcss%rcgs%hi_src(j)
         sr = c2*(rcss%u_src(j+1) - rcss%u_src(j))*rcss%rcgs%hi_src(j)
         if (sl*sr > c0) then
            sc = (rcss%u_src(j+1) - rcss%u_src(j-1))*rcss%rcgs%hci_src(j)
            sc = sign(min(abs(sl), abs(sr), abs(sc)), sc)
            if ( (rcss%u_src(j-1) - rcss%uel(j)) &
                *(rcss%u_src(j  ) - rcss%uel(j)) > c0) &
               rcss%uel(j) = rcss%u_src(j) &
                           - sign(min(c1_2*rcss%rcgs%h_src(j)*abs(sc), &
                                      abs(rcss%uel(j) - rcss%u_src(j))), sc)
            if ( (rcss%u_src(j+1) - rcss%uer(j)) &
                *(rcss%u_src(j  ) - rcss%uer(j)) > c0) &
               rcss%uer(j) = rcss%u_src(j) &
                           + sign(min(c1_2*rcss%rcgs%h_src(j)*abs(sc), &
                                      abs(rcss%uer(j) - rcss%u_src(j))), sc)
         else
            rcss%uel(j) = rcss%u_src(j)
            rcss%uer(j) = rcss%u_src(j)
         endif
      enddo

!     do j = 2, ns-1
      do j = 3, ns-1
         if ( (rcss%uel(j) - rcss%uer(j-1)) &
             *(rcss%u_src(j) - rcss%u_src(j-1)) < c0) then
            rcss%uel(j) = c1_2*(rcss%uer(j-1) + rcss%uel(j))
            rcss%uer(j-1) = rcss%uel(j)
         endif
      enddo

      do j = 2, ns-1
         d = rcss%uer(j) - rcss%uel(j)
         q = d*(c2*rcss%u_src(j) - rcss%uel(j) - rcss%uer(j))
         r = c1_3*d*d
         if     (  q > r) then
            rcss%uel(j) = c3*rcss%u_src(j) - c2*rcss%uer(j)
         elseif (- r > q) then
            rcss%uer(j) = c3*rcss%u_src(j) - c2*rcss%uel(j)
         endif
      enddo

   end subroutine limit_ppm_interior_monotonic

   pure subroutine limit_ppm_interior_non_oscillatory(rcss)
   ! ---------------------------------------------------------------------------
   ! Apply limiting to prevent a oscillatory reconstruction of piecewise
   ! parabolas for interior grid cells.
   ! ---------------------------------------------------------------------------

      type(recon_src_struct), intent(inout) :: rcss

      real(r8), dimension(rcss%rcgs%n_src_actual) :: d2
      real(r8) :: sl, sr, sc, d, q, r
      integer :: ns, j

      ns = rcss%rcgs%n_src_actual

      ! Obtain values proportional to the second derivative of the unlimited
      ! parabolas.
      do j = 1, ns
         d2(j) = rcss%uel(j) - c2*rcss%u_src(j) + rcss%uer(j)
      enddo

      do j = 2, ns-1
         ! Only apply limiting if the sign of the second
         ! derivative differs from the sign of second derivatives
         ! of any of the neighbouring parabolas.
         if (d2(j-1)*d2(j) < c0 .or. d2(j)*d2(j+1) < c0) then
            sl = c2*(rcss%u_src(j) - rcss%u_src(j-1))*rcss%rcgs%hi_src(j)
            sr = c2*(rcss%u_src(j+1) - rcss%u_src(j))*rcss%rcgs%hi_src(j)
            if (sl*sr > c0) then
               sc = (rcss%u_src(j+1) - rcss%u_src(j-1))*rcss%rcgs%hci_src(j)
               sc = sign(min(abs(sl), abs(sr), abs(sc)), sc)
               if ( (rcss%u_src(j-1) - rcss%uel(j)) &
                   *(rcss%u_src(j  ) - rcss%uel(j)) > c0) &
                  rcss%uel(j) = rcss%u_src(j) &
                              - sign(min(c1_2*rcss%rcgs%h_src(j)*abs(sc), &
                                         abs(rcss%uel(j) - rcss%u_src(j))), sc)
               if ( (rcss%u_src(j+1) - rcss%uer(j)) &
                   *(rcss%u_src(j  ) - rcss%uer(j)) > c0) &
                  rcss%uer(j) = rcss%u_src(j) &
                              + sign(min(c1_2*rcss%rcgs%h_src(j)*abs(sc), &
                                         abs(rcss%uer(j) - rcss%u_src(j))), sc)
            else
               rcss%uel(j) = rcss%u_src(j)
               rcss%uer(j) = rcss%u_src(j)
            endif
         endif
      enddo

!     do j = 2, ns-1
      do j = 3, ns-1
         if ( (rcss%uel(j) - rcss%uer(j-1)) &
             *(rcss%u_src(j) - rcss%u_src(j-1)) < c0) then
            rcss%uel(j) = c1_2*(rcss%uer(j-1) + rcss%uel(j))
            rcss%uer(j-1) = rcss%uel(j)
         endif
      enddo

      do j = 2, ns-1
         if (d2(j-1)*d2(j) < c0 .or. d2(j)*d2(j+1) < c0) then
            d = rcss%uer(j) - rcss%uel(j)
            q = d*(c2*rcss%u_src(j) - rcss%uel(j) - rcss%uer(j))
            r = c1_3*d*d
            if     (  q > r) then
               rcss%uel(j) = c3*rcss%u_src(j) - c2*rcss%uer(j)
            elseif (- r > q) then
               rcss%uer(j) = c3*rcss%u_src(j) - c2*rcss%uel(j)
            endif
         endif
      enddo

   end subroutine limit_ppm_interior_non_oscillatory

   pure subroutine limit_ppm_boundary(rcss)
   ! ---------------------------------------------------------------------------
   ! Handle piecewise parabola limiting of boundary cells.
   ! ---------------------------------------------------------------------------

      type(recon_src_struct), intent(inout) :: rcss

      real(r8) :: s
      integer :: ns

      ns = rcss%rcgs%n_src_actual

      if (rcss%pc_left_bndr) then
         ! Piecewise constant reconstruction of the left boundary cell.
         rcss%uel(1) = rcss%u_src(1)
         rcss%uer(1) = rcss%u_src(1)
      else
         ! Do not treat the left boundary cell as a local extrema, but ensure
         ! that the piecewise parabola is monotonic within the cell.
         if ( (rcss%u_src(2) - rcss%uer(1)) &
             *(rcss%u_src(1) - rcss%uer(1)) > c0) then
            rcss%uel(1) = rcss%u_src(1)
            rcss%uer(1) = rcss%u_src(1)
         else
            s = c2*(rcss%u_src(3) - rcss%u_src(2)) &
                  /(rcss%rcgs%h_src(2) + rcss%rcgs%h_src(3))
            if (s > 0) then
               rcss%uer(1) = &
                  max(rcss%u_src(1), &
                      min(rcss%uer(1), &
                          rcss%u_src(1) + c1_3*s*rcss%rcgs%h_src(1)))
            else
               rcss%uer(1) = &
                  min(rcss%u_src(1), &
                      max(rcss%uer(1), &
                          rcss%u_src(1) + c1_3*s*rcss%rcgs%h_src(1)))
            endif
            rcss%uel(1) = c1_2*(c3*rcss%u_src(1) - rcss%uer(1))
         endif
      endif

      if (rcss%pc_right_bndr) then
         ! Piecewise constant reconstruction of the right boundary cell.
         rcss%uel(ns) = rcss%u_src(ns)
         rcss%uer(ns) = rcss%u_src(ns)
      else
         ! Do not treat the right boundary cell as a local extrema, but ensure
         ! that the piecewise parabola is monotonic within the cell.
         if ( (rcss%u_src(ns  ) - rcss%uel(ns)) &
             *(rcss%u_src(ns-1) - rcss%uel(ns)) > c0) then
            rcss%uel(ns) = rcss%u_src(ns)
            rcss%uer(ns) = rcss%u_src(ns)
         else
            s = c2*(rcss%u_src(ns-1) - rcss%u_src(ns-2)) &
                  /(rcss%rcgs%h_src(ns-2) + rcss%rcgs%h_src(ns-1))
            if (s > 0) then
               rcss%uel(ns) = &
                  min(rcss%u_src(ns), &
                      max(rcss%uel(ns), &
                          rcss%u_src(ns) - c1_3*s*rcss%rcgs%h_src(ns)))
            else
               rcss%uel(ns) = &
                  max(rcss%u_src(ns), &
                      min(rcss%uel(ns), &
                          rcss%u_src(ns) - c1_3*s*rcss%rcgs%h_src(ns)))
            endif
            rcss%uer(ns) = c1_2*(c3*rcss%u_src(ns) - rcss%uel(ns))
         endif
      endif

   end subroutine limit_ppm_boundary

   pure subroutine limit_ppm_posdef(rcss)
   ! ---------------------------------------------------------------------------
   ! Modify piecewise parabolas so they are never negative within the grid cell.
   ! ---------------------------------------------------------------------------

      type(recon_src_struct), intent(inout) :: rcss

      real(r8) :: min_u_0, sl, a2, sr, q
      integer :: j

      do j = 1, rcss%rcgs%n_src_actual
         min_u_0 = min(rcss%u_src(j), c0)
         rcss%uel(j) = max(rcss%uel(j), min_u_0)
         rcss%uer(j) = max(rcss%uer(j), min_u_0)
         sl = c2*(c3*rcss%u_src(j) - c2*rcss%uel(j) - rcss%uer(j))
         a2 = c3*(rcss%uel(j) - c2*rcss%u_src(j) + rcss%uer(j))
         sr = sl + c2*a2
         if (sl < c0 .and. sr > c0) then
            if (a2*rcss%uel(j) - c1_4*sl*sl < a2*min_u_0) then
               q = c3*rcss%u_src(j)/(c3*sl*sr + c4*a2*a2)
               rcss%uel(j) = sl*sl*q
               rcss%uer(j) = sr*sr*q
            endif
         endif
      enddo

   end subroutine limit_ppm_posdef

   pure subroutine polycoeff_ppm(rcss)
   ! ---------------------------------------------------------------------------
   ! Obtain coefficients for piecewise parabolas from grid cell means and left
   ! and right edge values.
   ! ---------------------------------------------------------------------------

      type(recon_src_struct), intent(inout) :: rcss

      integer :: j

      do j = 1, rcss%rcgs%n_src_actual
         rcss%polycoeff(1,j) = rcss%uel(j)
         rcss%polycoeff(2,j) = c6*rcss%u_src(j) - c4*rcss%uel(j) &
                             - c2*rcss%uer(j)
         rcss%polycoeff(3,j) = c3*(rcss%uel(j) - c2*rcss%u_src(j) + rcss%uer(j))
      enddo

   end subroutine polycoeff_ppm

   pure subroutine limit_pqm_monotonic(rcss)
   ! ---------------------------------------------------------------------------
   ! Apply limiting to ensure a monotonic reconstruction of piecewise quartics
   ! for interior grid cells.
   ! ---------------------------------------------------------------------------

      type(recon_src_struct), intent(inout) :: rcss

      real(r8), dimension(rcss%rcgs%n_src_actual) :: sl, sr, sc
      real(r8) :: a0, a1, a2, a3, b0, b1, b2, q1, q2, q3, s, xi
      integer :: ns, j
      logical :: incon_inflex

      ns = rcss%rcgs%n_src_actual

      do j = 2, ns-1
         sl(j) = c2*(rcss%u_src(j) - rcss%u_src(j-1))*rcss%rcgs%hi_src(j)
         sr(j) = c2*(rcss%u_src(j+1) - rcss%u_src(j))*rcss%rcgs%hi_src(j)
         sc(j) = (rcss%u_src(j+1) - rcss%u_src(j-1))*rcss%rcgs%hci_src(j)
         sc(j) = sign(min(abs(sl(j)), abs(sr(j)), abs(sc(j))), sc(j))
         if (sl(j)*sr(j) > c0) then
            if ( (rcss%u_src(j-1) - rcss%uel(j)) &
                *(rcss%u_src(j  ) - rcss%uel(j)) > c0) &
               rcss%uel(j) = rcss%u_src(j) &
                           - sign(min(c1_2*rcss%rcgs%h_src(j)*abs(sc(j)), &
                                      abs(rcss%uel(j) - rcss%u_src(j))), sc(j))
            if ( (rcss%u_src(j+1) - rcss%uer(j)) &
                *(rcss%u_src(j  ) - rcss%uer(j)) > c0) &
               rcss%uer(j) = rcss%u_src(j) &
                           + sign(min(c1_2*rcss%rcgs%h_src(j)*abs(sc(j)), &
                                      abs(rcss%uer(j) - rcss%u_src(j))), sc(j))
!           if (rcss%usl(j)*sc(j) < c0) rcss%usl(j) = sc(j)
!           if (rcss%usr(j)*sc(j) < c0) rcss%usr(j) = sc(j)
            if (rcss%usl(j)*sc(j) < c0) rcss%usl(j) = c0
            if (rcss%usr(j)*sc(j) < c0) rcss%usr(j) = c0
         else
            rcss%uel(j) = rcss%u_src(j)
            rcss%uer(j) = rcss%u_src(j)
            rcss%usl(j) = c0
            rcss%usr(j) = c0
         endif
      enddo

      do j = 3, ns-1
         if ( (rcss%uel(j) - rcss%uer(j-1)) &
             *(rcss%u_src(j) - rcss%u_src(j-1)) < c0) then
            rcss%uel(j) = c1_2*(rcss%uer(j-1) + rcss%uel(j))
            rcss%uer(j-1) = rcss%uel(j)
         endif
      enddo

      do j = 2, ns-1

         ! Compute polynomial coefficients for 1. derivative of the
         ! reconstruction.
         a0 = rcss%usl(j)
         a1 = c2*( c30*rcss%u_src(j) - c18*rcss%uel(j) - c12*rcss%uer(j) &
                 - c9_2*rcss%usl(j) + c3_2*rcss%usr(j))
         a2 = c3*(- c60*rcss%u_src(j) + c32*rcss%uel(j) + c28*rcss%uer(j) &
                  + c6*rcss%usl(j) - c4*rcss%usr(j))
         a3 = c4*( c30*rcss%u_src(j) - c15*(rcss%uel(j) + rcss%uer(j)) &
                 - c5_2*(rcss%usl(j) - rcss%usr(j)))

         ! Compute polynomial coefficients for 2. derivative of the
         ! reconstruction.
         b0 = a1
         b1 = c2*a2
         b2 = c3*a3

         ! Check for inconsistent inflextion points.
         incon_inflex = .false.
         q1 = b0*b2
         q2 = b1*b1 - c4*q1
         if (q2 > c0) then
            if (b0*(b0 + b1 + b2) < c0) then
               ! One inflection point.
               if (abs(b2) < rcss%u_eps) then
                  if (abs(b1) > rcss%u_eps) then
                     xi = - b0/b1
                     if ((a0 + xi*(a1 + xi*(a2 + xi*a3)))*sc(j) < c0) &
                        incon_inflex = .true.
                  endif
               else
                  q3 = c1_2/b2
                  s = sqrt(q2)
                  xi = - (b1 + s)*q3
                  if (xi > c0 .and. xi < c1) then
                     if ((a0 + xi*(a1 + xi*(a2 + xi*a3)))*sc(j) < c0) &
                        incon_inflex = .true.
                  else
                     xi = - (b1 - s)*q3
                     if ((a0 + xi*(a1 + xi*(a2 + xi*a3)))*sc(j) < c0) &
                        incon_inflex = .true.
                  endif
               endif
            elseif (q1 > rcss%uu_eps) then ! Should imply b2 != 0
               ! Two inflection points.
               q3 = c1_2/b2
               s = sqrt(q2)
               xi = - (b1 + s)*q3
               if ((a0 + xi*(a1 + xi*(a2 + xi*a3)))*sc(j) < c0) then
                  incon_inflex = .true.
               else
                  xi = - (b1 - s)*q3
                  if ((a0 + xi*(a1 + xi*(a2 + xi*a3)))*sc(j) < c0) &
                     incon_inflex = .true.
               endif
            endif
         endif

         if (incon_inflex) then
            if (abs(sl(j)) < abs(sr(j))) then
               rcss%usl(j) = c10_3*rcss%u_src(j) - c8_3*rcss%uel(j) &
                           - c2_3*rcss%uer(j)
               if (rcss%usl(j)*sc(j) < c0) then
                  rcss%usl(j) = c0
                  rcss%uer(j) = c5*rcss%u_src(j) - c4*rcss%uel(j)
                  rcss%usr(j) = c20*(rcss%u_src(j) - rcss%uel(j))
               else
                  rcss%usr(j) = c4*rcss%uel(j) + c6*rcss%uer(j) &
                              - c10*rcss%u_src(j)
                  if (rcss%usr(j)*sc(j) < c0) then
                     rcss%usr(j) = c0
                     rcss%uel(j) = c5_2*rcss%u_src(j) - c3_2*rcss%uer(j)
                     rcss%usl(j) = c10_3*(rcss%uer(j) - rcss%u_src(j))
                  endif
               endif
            else
               rcss%usr(j) = c8_3*rcss%uer(j) + c2_3*rcss%uel(j) &
                           - c10_3*rcss%u_src(j)
               if (rcss%usr(j)*sc(j) < c0) then
                  rcss%usr(j) = c0
                  rcss%uel(j) = c5*rcss%u_src(j) - c4*rcss%uer(j)
                  rcss%usl(j) = c20*(rcss%uer(j) - rcss%u_src(j))
               else
                  rcss%usl(j) = c10*rcss%u_src(j) - c4*rcss%uer(j) &
                              - c6*rcss%uel(j)
                  if (rcss%usl(j)*sc(j) < c0) then
                     rcss%usl(j) = c0
                     rcss%uer(j) = c5_2*rcss%u_src(j) - c3_2*rcss%uel(j)
                     rcss%usr(j) = c10_3*(rcss%u_src(j) - rcss%uel(j))
                  endif
               endif
            endif
         endif

      enddo

      if (rcss%pc_left_bndr) then
         ! Piecewise constant reconstruction of the left boundary cell.
         rcss%uel(1) = rcss%u_src(1)
         rcss%uer(1) = rcss%u_src(1)
         rcss%usl(1) = c0
         rcss%usr(1) = c0
      else
         ! Do not treat the left boundary cell as a local extrema, but ensure
         ! that the piecewise parabola is monotonic within the cell.
         if ( (rcss%u_src(2) - rcss%uer(1)) &
             *(rcss%u_src(1) - rcss%uer(1)) > c0) then
            rcss%uel(1) = rcss%u_src(1)
            rcss%uer(1) = rcss%u_src(1)
            rcss%usl(1) = c0
            rcss%usr(1) = c0
         else
            s = c2*(rcss%u_src(3) - rcss%u_src(2)) &
                  /(rcss%rcgs%h_src(2) + rcss%rcgs%h_src(3))
            if (s > 0) then
               rcss%uer(1) = &
                  max(rcss%u_src(1), &
                      min(rcss%uel(2), &
                          rcss%u_src(1) + c1_3*s*rcss%rcgs%h_src(1)))
            else
               rcss%uer(1) = &
                  min(rcss%u_src(1), &
                      max(rcss%uel(2), &
                          rcss%u_src(1) + c1_3*s*rcss%rcgs%h_src(1)))
            endif
            rcss%uel(1) = c1_2*(c3*rcss%u_src(1) - rcss%uer(1))
            rcss%usl(1) = c6*rcss%u_src(1) - c4*rcss%uel(1) - c2*rcss%uer(1)
            rcss%usr(1) = c2*rcss%uel(1) + c4*rcss%uer(1) - c6*rcss%u_src(1)
         endif
      endif

      if (rcss%pc_right_bndr) then
         ! Piecewise constant reconstruction of the right boundary cell.
         rcss%uel(ns) = rcss%u_src(ns)
         rcss%uer(ns) = rcss%u_src(ns)
         rcss%usl(ns) = c0
         rcss%usr(ns) = c0
      else
         ! Do not treat the right boundary cell as a local extrema, but ensure
         ! that the piecewise parabola is monotonic within the cell.
         if ( (rcss%u_src(ns  ) - rcss%uel(ns)) &
             *(rcss%u_src(ns-1) - rcss%uel(ns)) > c0) then
            rcss%uel(ns) = rcss%u_src(ns)
            rcss%uer(ns) = rcss%u_src(ns)
            rcss%usl(ns) = c0
            rcss%usr(ns) = c0
         else
            s = c2*(rcss%u_src(ns-1) - rcss%u_src(ns-2)) &
                  /(rcss%rcgs%h_src(ns-2) + rcss%rcgs%h_src(ns-1))
            if (s > 0) then
               rcss%uel(ns) = &
                  min(rcss%u_src(ns), &
                      max(rcss%uer(ns-1), &
                          rcss%u_src(ns) - c1_3*s*rcss%rcgs%h_src(ns)))
            else
               rcss%uel(ns) = &
                  max(rcss%u_src(ns), &
                      min(rcss%uer(ns-1), &
                          rcss%u_src(ns) - c1_3*s*rcss%rcgs%h_src(ns)))
            endif
            rcss%uer(ns) = c1_2*(c3*rcss%u_src(ns) - rcss%uel(ns))
            rcss%usl(ns) = c6*rcss%u_src(ns) - c4*rcss%uel(ns) - c2*rcss%uer(ns)
            rcss%usr(ns) = c2*rcss%uel(ns) + c4*rcss%uer(ns) - c6*rcss%u_src(ns)
         endif
      endif

   end subroutine limit_pqm_monotonic

   pure subroutine limit_pqm_non_oscillatory(rcss)
   ! ---------------------------------------------------------------------------
   ! Apply limiting to ensure a monotonic reconstruction of piecewise quartics
   ! for interior grid cells.
   ! ---------------------------------------------------------------------------

      type(recon_src_struct), intent(inout) :: rcss

      real(r8), dimension(rcss%rcgs%n_src_actual) :: d2, sl, sr, sc
      logical, dimension(rcss%rcgs%n_src_actual) :: smooth
      real(r8) :: a0, a1, a2, a3, b0, b1, b2, q1, q2, q3, s, xi
      integer :: ns, j
      logical :: inflex, incon_inflex

      ns = rcss%rcgs%n_src_actual

      ! Obtain values proportional to the second derivative of the unlimited
      ! parabolas.
      do j = 1, ns
         d2(j) = rcss%uel(j) - c2*rcss%u_src(j) + rcss%uer(j)
      enddo

      do j = 2, ns-1
         ! Set flag if the reconstruction is considered smooth, that is
         ! the sign of the second derivative equals the sign of the second
         ! derivatives of both the neighbouring parabolas.
         smooth(j) = d2(j-1)*d2(j) >= c0 .and. d2(j)*d2(j+1) >= c0

         if (smooth(j)) then

            ! Slopes of a parabolic reconstruction.
            sl(j) = c6*rcss%u_src(j) - c4*rcss%uel(j) - c2*rcss%uer(j)
            sr(j) = c2*rcss%uel(j) + c4*rcss%uer(j) - c6*rcss%u_src(j)

            if (sl(j) < c0 .and. sr(j) > c0) then

               ! If the slopes of a parabolic reconstruction has different
               ! signs, the parabolic reconstruction is chosen.
               rcss%usl(j) = sl(j)
               rcss%usr(j) = sr(j)

            else

               ! If the quartic reconstruction has one or more inflextion
               ! points, a parabolic reconstruction is chosen.
               b0 = c2*( c30*rcss%u_src(j) - c18*rcss%uel(j) - c12*rcss%uer(j) &
                       - c9_2*rcss%usl(j) + c3_2*rcss%usr(j))
               b1 = c6*(- c60*rcss%u_src(j) + c32*rcss%uel(j) &
                        + c28*rcss%uer(j) + c6*rcss%usl(j) - c4*rcss%usr(j))
               b2 = c12*( c30*rcss%u_src(j) - c15*(rcss%uel(j) + rcss%uer(j)) &
                        - c5_2*(rcss%usl(j) - rcss%usr(j)))
               q1 = b0*b2
               q2 = b1*b1 - c4*q1
               if (q2 > c0 .and. &
                   (b0*(b0 + b1 + b2) < c0 .or. q1 > rcss%uu_eps)) then
                  inflex = .true.
               else
                  inflex = .false.
               endif
               if (inflex) then
                  rcss%usl(j) = sl(j)
                  rcss%usr(j) = sr(j)
               endif

            endif

         else

            ! Apply limiting for unsmooth reconstruction.
            sl(j) = c2*(rcss%u_src(j) - rcss%u_src(j-1))*rcss%rcgs%hi_src(j)
            sr(j) = c2*(rcss%u_src(j+1) - rcss%u_src(j))*rcss%rcgs%hi_src(j)
            sc(j) = (rcss%u_src(j+1) - rcss%u_src(j-1))*rcss%rcgs%hci_src(j)
            sc(j) = sign(min(abs(sl(j)), abs(sr(j)), abs(sc(j))), sc(j))
            if (sl(j)*sr(j) > c0) then
               if ( (rcss%u_src(j-1) - rcss%uel(j)) &
                   *(rcss%u_src(j  ) - rcss%uel(j)) > c0) &
                  rcss%uel(j) = rcss%u_src(j) &
                              - sign(min(c1_2*rcss%rcgs%h_src(j)*abs(sc(j)), &
                                         abs(rcss%uel(j) - rcss%u_src(j))), &
                                     sc(j))
               if ( (rcss%u_src(j+1) - rcss%uer(j)) &
                   *(rcss%u_src(j  ) - rcss%uer(j)) > c0) &
                  rcss%uer(j) = rcss%u_src(j) &
                              + sign(min(c1_2*rcss%rcgs%h_src(j)*abs(sc(j)), &
                                         abs(rcss%uer(j) - rcss%u_src(j))), &
                                     sc(j))
!              if (rcss%usl(j)*sc(j) < c0) rcss%usl(j) = sc(j)
!              if (rcss%usr(j)*sc(j) < c0) rcss%usr(j) = sc(j)
               if (rcss%usl(j)*sc(j) < c0) rcss%usl(j) = c0
               if (rcss%usr(j)*sc(j) < c0) rcss%usr(j) = c0
            else
               rcss%uel(j) = rcss%u_src(j)
               rcss%uer(j) = rcss%u_src(j)
               rcss%usl(j) = c0
               rcss%usr(j) = c0
            endif

         endif
      enddo

      do j = 3, ns-1
         if ( (rcss%uel(j) - rcss%uer(j-1)) &
             *(rcss%u_src(j) - rcss%u_src(j-1)) < c0) then
            if     (smooth(j-1)) then
               rcss%uel(j) = rcss%uer(j-1)
            elseif (smooth(j  )) then
               rcss%uer(j-1) = rcss%uel(j)
            else
               rcss%uel(j) = c1_2*(rcss%uer(j-1) + rcss%uel(j))
               rcss%uer(j-1) = rcss%uel(j)
            endif
         endif
      enddo

      do j = 2, ns-1

         if (.not.smooth(j)) then

            ! Compute polynomial coefficients for 1. derivative of the
            ! reconstruction.
            a0 = rcss%usl(j)
            a1 = c2*( c30*rcss%u_src(j) - c18*rcss%uel(j) - c12*rcss%uer(j) &
                    - c9_2*rcss%usl(j) + c3_2*rcss%usr(j))
            a2 = c3*(- c60*rcss%u_src(j) + c32*rcss%uel(j) + c28*rcss%uer(j) &
                     + c6*rcss%usl(j) - c4*rcss%usr(j))
            a3 = c4*( c30*rcss%u_src(j) - c15*(rcss%uel(j) + rcss%uer(j)) &
                    - c5_2*(rcss%usl(j) - rcss%usr(j)))

            ! Compute polynomial coefficients for 2. derivative of the
            ! reconstruction.
            b0 = a1
            b1 = c2*a2
            b2 = c3*a3

            ! Check for inconsistent inflextion points.
            incon_inflex = .false.
            q1 = b0*b2
            q2 = b1*b1 - c4*q1
            if (q2 > c0) then
               if (b0*(b0 + b1 + b2) < c0) then
                  ! One inflection point.
                  if (abs(b2) < rcss%u_eps) then
                     if (abs(b1) > rcss%u_eps) then
                        xi = - b0/b1
                        if ((a0 + xi*(a1 + xi*(a2 + xi*a3)))*sc(j) < c0) &
                           incon_inflex = .true.
                     endif
                  else
                     q3 = c1_2/b2
                     s = sqrt(q2)
                     xi = - (b1 + s)*q3
                     if (xi > c0 .and. xi < c1) then
                        if ((a0 + xi*(a1 + xi*(a2 + xi*a3)))*sc(j) < c0) &
                           incon_inflex = .true.
                     else
                        xi = - (b1 - s)*q3
                        if ((a0 + xi*(a1 + xi*(a2 + xi*a3)))*sc(j) < c0) &
                           incon_inflex = .true.
                     endif
                  endif
               elseif (q1 > rcss%uu_eps) then ! Should imply b2 != 0
                  ! Two inflection points.
                  q3 = c1_2/b2
                  s = sqrt(q2)
                  xi = - (b1 + s)*q3
                  if ((a0 + xi*(a1 + xi*(a2 + xi*a3)))*sc(j) < c0) then
                     incon_inflex = .true.
                  else
                     xi = - (b1 - s)*q3
                     if ((a0 + xi*(a1 + xi*(a2 + xi*a3)))*sc(j) < c0) &
                        incon_inflex = .true.
                  endif
               endif
            endif

            if (incon_inflex) then
               if (abs(sl(j)) < abs(sr(j))) then
                  rcss%usl(j) = c10_3*rcss%u_src(j) - c8_3*rcss%uel(j) &
                              - c2_3*rcss%uer(j)
                  if (rcss%usl(j)*sc(j) < c0) then
                     rcss%usl(j) = c0
                     rcss%uer(j) = c5*rcss%u_src(j) - c4*rcss%uel(j)
                     rcss%usr(j) = c20*(rcss%u_src(j) - rcss%uel(j))
                  else
                     rcss%usr(j) = c4*rcss%uel(j) + c6*rcss%uer(j) &
                                 - c10*rcss%u_src(j)
                     if (rcss%usr(j)*sc(j) < c0) then
                        rcss%usr(j) = c0
                        rcss%uel(j) = c5_2*rcss%u_src(j) - c3_2*rcss%uer(j)
                        rcss%usl(j) = c10_3*(rcss%uer(j) - rcss%u_src(j))
                     endif
                  endif
               else
                  rcss%usr(j) = c8_3*rcss%uer(j) + c2_3*rcss%uel(j) &
                              - c10_3*rcss%u_src(j)
                  if (rcss%usr(j)*sc(j) < c0) then
                     rcss%usr(j) = c0
                     rcss%uel(j) = c5*rcss%u_src(j) - c4*rcss%uer(j)
                     rcss%usl(j) = c20*(rcss%uer(j) - rcss%u_src(j))
                  else
                     rcss%usl(j) = c10*rcss%u_src(j) - c4*rcss%uer(j) &
                                 - c6*rcss%uel(j)
                     if (rcss%usl(j)*sc(j) < c0) then
                        rcss%usl(j) = c0
                        rcss%uer(j) = c5_2*rcss%u_src(j) - c3_2*rcss%uel(j)
                        rcss%usr(j) = c10_3*(rcss%u_src(j) - rcss%uel(j))
                     endif
                  endif
               endif
            endif

         endif

      enddo

      if (rcss%pc_left_bndr) then
         ! Piecewise constant reconstruction of the left boundary cell.
         rcss%uel(1) = rcss%u_src(1)
         rcss%uer(1) = rcss%u_src(1)
         rcss%usl(1) = c0
         rcss%usr(1) = c0
      else
         ! Do not treat the left boundary cell as a local extrema, but ensure
         ! that the piecewise parabola is monotonic within the cell.
         if ( (rcss%u_src(2) - rcss%uer(1)) &
             *(rcss%u_src(1) - rcss%uer(1)) > c0) then
            rcss%uel(1) = rcss%u_src(1)
            rcss%uer(1) = rcss%u_src(1)
            rcss%usl(1) = c0
            rcss%usr(1) = c0
         else
            s = c2*(rcss%u_src(3) - rcss%u_src(2)) &
                  /(rcss%rcgs%h_src(2) + rcss%rcgs%h_src(3))
            if (s > 0) then
               rcss%uer(1) = &
                  max(rcss%u_src(1), &
                      min(rcss%uel(2), &
                          rcss%u_src(1) + c1_3*s*rcss%rcgs%h_src(1)))
            else
               rcss%uer(1) = &
                  min(rcss%u_src(1), &
                      max(rcss%uel(2), &
                          rcss%u_src(1) + c1_3*s*rcss%rcgs%h_src(1)))
            endif
            rcss%uel(1) = c1_2*(c3*rcss%u_src(1) - rcss%uer(1))
            rcss%usl(1) = c6*rcss%u_src(1) - c4*rcss%uel(1) - c2*rcss%uer(1)
            rcss%usr(1) = c2*rcss%uel(1) + c4*rcss%uer(1) - c6*rcss%u_src(1)
         endif
      endif

      if (rcss%pc_right_bndr) then
         ! Piecewise constant reconstruction of the right boundary cell.
         rcss%uel(ns) = rcss%u_src(ns)
         rcss%uer(ns) = rcss%u_src(ns)
         rcss%usl(ns) = c0
         rcss%usr(ns) = c0
      else
         ! Do not treat the right boundary cell as a local extrema, but ensure
         ! that the piecewise parabola is monotonic within the cell.
         if ( (rcss%u_src(ns  ) - rcss%uel(ns)) &
             *(rcss%u_src(ns-1) - rcss%uel(ns)) > c0) then
            rcss%uel(ns) = rcss%u_src(ns)
            rcss%uer(ns) = rcss%u_src(ns)
            rcss%usl(ns) = c0
            rcss%usr(ns) = c0
         else
            s = c2*(rcss%u_src(ns-1) - rcss%u_src(ns-2)) &
                  /(rcss%rcgs%h_src(ns-2) + rcss%rcgs%h_src(ns-1))
            if (s > 0) then
               rcss%uel(ns) = &
                  min(rcss%u_src(ns), &
                      max(rcss%uer(ns-1), &
                          rcss%u_src(ns) - c1_3*s*rcss%rcgs%h_src(ns)))
            else
               rcss%uel(ns) = &
                  max(rcss%u_src(ns), &
                      min(rcss%uer(ns-1), &
                          rcss%u_src(ns) - c1_3*s*rcss%rcgs%h_src(ns)))
            endif
            rcss%uer(ns) = c1_2*(c3*rcss%u_src(ns) - rcss%uel(ns))
            rcss%usl(ns) = c6*rcss%u_src(ns) - c4*rcss%uel(ns) - c2*rcss%uer(ns)
            rcss%usr(ns) = c2*rcss%uel(ns) + c4*rcss%uer(ns) - c6*rcss%u_src(ns)
         endif
      endif

   end subroutine limit_pqm_non_oscillatory

   pure subroutine limit_pqm_non_oscillatory_posdef(rcss)
   ! ---------------------------------------------------------------------------
   ! Apply limiting to ensure a monotonic reconstruction of piecewise quartics
   ! for interior grid cells.
   ! ---------------------------------------------------------------------------

      type(recon_src_struct), intent(inout) :: rcss

      real(r8), dimension(rcss%rcgs%n_src_actual) :: d2, sl, sr, sc
      logical, dimension(rcss%rcgs%n_src_actual) :: smooth
      real(r8) :: min_u_0, a0, a1, a2, a3, b0, b1, b2, q1, q2, q3, s, xi
      integer :: ns, j
      logical :: inflex, incon_inflex

      ns = rcss%rcgs%n_src_actual

      ! Obtain values proportional to the second derivative of the unlimited
      ! parabolas.
      do j = 1, ns
         d2(j) = rcss%uel(j) - c2*rcss%u_src(j) + rcss%uer(j)
      enddo

      do j = 2, ns-1
         ! Set flag if the reconstruction is considered smooth, that is
         ! the sign of the second derivative equals the sign of the second
         ! derivatives of both the neighbouring parabolas.
         smooth(j) = d2(j-1)*d2(j) >= c0 .and. d2(j)*d2(j+1) >= c0

         if (smooth(j)) then

            ! Ensure edge values of smooth reconstruction is positive definite.
            min_u_0 = min(rcss%u_src(j), c0)
            rcss%uel(j) = max(rcss%uel(j), min_u_0)
            rcss%uer(j) = max(rcss%uer(j), min_u_0)

            ! Slopes of a parabolic reconstruction.
            sl(j) = c6*rcss%u_src(j) - c4*rcss%uel(j) - c2*rcss%uer(j)
            sr(j) = c2*rcss%uel(j) + c4*rcss%uer(j) - c6*rcss%u_src(j)

            if (sl(j) < c0 .and. sr(j) > c0) then

               ! If the slopes of a parabolic reconstruction has different
               ! signs, the parabolic reconstruction is chosen. If needed,
               ! modify the parabola it is positive definite.
               a2 = c1_2*(sr(j) - sl(j))
               if (a2*rcss%uel(j) - c1_4*sl(j)*sl(j) < a2*min_u_0) then
                  q1 = c3*rcss%u_src(j)/(c3*sl(j)*sr(j) + c4*a2*a2)
                  rcss%uel(j) = sl(j)*sl(j)*q1
                  rcss%uer(j) = sr(j)*sr(j)*q1
                  rcss%usl(j) = c6*rcss%u_src(j) - c4*rcss%uel(j) &
                              - c2*rcss%uer(j)
                  rcss%usr(j) = c2*rcss%uel(j) + c4*rcss%uer(j) &
                              - c6*rcss%u_src(j)
               else
                  rcss%usl(j) = sl(j)
                  rcss%usr(j) = sr(j)
               endif

            else

               ! If the quartic reconstruction has one or more inflextion
               ! points, a parabolic reconstruction is chosen.
               b0 = c2*( c30*rcss%u_src(j) - c18*rcss%uel(j) - c12*rcss%uer(j) &
                       - c9_2*rcss%usl(j) + c3_2*rcss%usr(j))
               b1 = c6*(- c60*rcss%u_src(j) + c32*rcss%uel(j) &
                        + c28*rcss%uer(j) + c6*rcss%usl(j) - c4*rcss%usr(j))
               b2 = c12*( c30*rcss%u_src(j) - c15*(rcss%uel(j) + rcss%uer(j)) &
                        - c5_2*(rcss%usl(j) - rcss%usr(j)))
               q1 = b0*b2
               q2 = b1*b1 - c4*q1
               if (q2 > c0 .and. &
                   (b0*(b0 + b1 + b2) < c0 .or. q1 > rcss%uu_eps)) then
                  inflex = .true.
               else
                  inflex = .false.
               endif
               if (inflex) then
                  rcss%usl(j) = sl(j)
                  rcss%usr(j) = sr(j)
               endif

            endif

         else

            ! Apply limiting for unsmooth reconstruction.
            sl(j) = c2*(rcss%u_src(j) - rcss%u_src(j-1))*rcss%rcgs%hi_src(j)
            sr(j) = c2*(rcss%u_src(j+1) - rcss%u_src(j))*rcss%rcgs%hi_src(j)
            sc(j) = (rcss%u_src(j+1) - rcss%u_src(j-1))*rcss%rcgs%hci_src(j)
            sc(j) = sign(min(abs(sl(j)), abs(sr(j)), abs(sc(j))), sc(j))
            if (sl(j)*sr(j) > c0) then
               if ( (rcss%u_src(j-1) - rcss%uel(j)) &
                   *(rcss%u_src(j  ) - rcss%uel(j)) > c0) &
                  rcss%uel(j) = rcss%u_src(j) &
                              - sign(min(c1_2*rcss%rcgs%h_src(j)*abs(sc(j)), &
                                         abs(rcss%uel(j) - rcss%u_src(j))), &
                                     sc(j))
               if ( (rcss%u_src(j+1) - rcss%uer(j)) &
                   *(rcss%u_src(j  ) - rcss%uer(j)) > c0) &
                  rcss%uer(j) = rcss%u_src(j) &
                              + sign(min(c1_2*rcss%rcgs%h_src(j)*abs(sc(j)), &
                                         abs(rcss%uer(j) - rcss%u_src(j))), &
                                     sc(j))
!              if (rcss%usl(j)*sc(j) < c0) rcss%usl(j) = sc(j)
!              if (rcss%usr(j)*sc(j) < c0) rcss%usr(j) = sc(j)
               if (rcss%usl(j)*sc(j) < c0) rcss%usl(j) = c0
               if (rcss%usr(j)*sc(j) < c0) rcss%usr(j) = c0
            else
               rcss%uel(j) = rcss%u_src(j)
               rcss%uer(j) = rcss%u_src(j)
               rcss%usl(j) = c0
               rcss%usr(j) = c0
            endif

         endif
      enddo

      do j = 3, ns-1
         if ( (rcss%uel(j) - rcss%uer(j-1)) &
             *(rcss%u_src(j) - rcss%u_src(j-1)) < c0) then
            if     (smooth(j-1)) then
               rcss%uel(j) = rcss%uer(j-1)
            elseif (smooth(j  )) then
               rcss%uer(j-1) = rcss%uel(j)
            else
               rcss%uel(j) = c1_2*(rcss%uer(j-1) + rcss%uel(j))
               rcss%uer(j-1) = rcss%uel(j)
            endif
         endif
      enddo

      do j = 2, ns-1

         if (.not.smooth(j)) then

            ! Compute polynomial coefficients for 1. derivative of the
            ! reconstruction.
            a0 = rcss%usl(j)
            a1 = c2*( c30*rcss%u_src(j) - c18*rcss%uel(j) - c12*rcss%uer(j) &
                    - c9_2*rcss%usl(j) + c3_2*rcss%usr(j))
            a2 = c3*(- c60*rcss%u_src(j) + c32*rcss%uel(j) + c28*rcss%uer(j) &
                     + c6*rcss%usl(j) - c4*rcss%usr(j))
            a3 = c4*( c30*rcss%u_src(j) - c15*(rcss%uel(j) + rcss%uer(j)) &
                    - c5_2*(rcss%usl(j) - rcss%usr(j)))

            ! Compute polynomial coefficients for 2. derivative of the
            ! reconstruction.
            b0 = a1
            b1 = c2*a2
            b2 = c3*a3

            ! Check for inconsistent inflextion points.
            incon_inflex = .false.
            q1 = b0*b2
            q2 = b1*b1 - c4*q1
            if (q2 > c0) then
               if (b0*(b0 + b1 + b2) < c0) then
                  ! One inflection point.
                  if (abs(b2) < rcss%u_eps) then
                     if (abs(b1) > rcss%u_eps) then
                        xi = - b0/b1
                        if ((a0 + xi*(a1 + xi*(a2 + xi*a3)))*sc(j) < c0) &
                           incon_inflex = .true.
                     endif
                  else
                     q3 = c1_2/b2
                     s = sqrt(q2)
                     xi = - (b1 + s)*q3
                     if (xi > c0 .and. xi < c1) then
                        if ((a0 + xi*(a1 + xi*(a2 + xi*a3)))*sc(j) < c0) &
                           incon_inflex = .true.
                     else
                        xi = - (b1 - s)*q3
                        if ((a0 + xi*(a1 + xi*(a2 + xi*a3)))*sc(j) < c0) &
                           incon_inflex = .true.
                     endif
                  endif
               elseif (q1 > rcss%uu_eps) then ! Should imply b2 != 0
                  ! Two inflection points.
                  q3 = c1_2/b2
                  s = sqrt(q2)
                  xi = - (b1 + s)*q3
                  if ((a0 + xi*(a1 + xi*(a2 + xi*a3)))*sc(j) < c0) then
                     incon_inflex = .true.
                  else
                     xi = - (b1 - s)*q3
                     if ((a0 + xi*(a1 + xi*(a2 + xi*a3)))*sc(j) < c0) &
                        incon_inflex = .true.
                  endif
               endif
            endif

            if (incon_inflex) then
               if (abs(sl(j)) < abs(sr(j))) then
                  rcss%usl(j) = c10_3*rcss%u_src(j) - c8_3*rcss%uel(j) &
                              - c2_3*rcss%uer(j)
                  if (rcss%usl(j)*sc(j) < c0) then
                     rcss%usl(j) = c0
                     rcss%uer(j) = c5*rcss%u_src(j) - c4*rcss%uel(j)
                     rcss%usr(j) = c20*(rcss%u_src(j) - rcss%uel(j))
                  else
                     rcss%usr(j) = c4*rcss%uel(j) + c6*rcss%uer(j) &
                                 - c10*rcss%u_src(j)
                     if (rcss%usr(j)*sc(j) < c0) then
                        rcss%usr(j) = c0
                        rcss%uel(j) = c5_2*rcss%u_src(j) - c3_2*rcss%uer(j)
                        rcss%usl(j) = c10_3*(rcss%uer(j) - rcss%u_src(j))
                     endif
                  endif
               else
                  rcss%usr(j) = c8_3*rcss%uer(j) + c2_3*rcss%uel(j) &
                              - c10_3*rcss%u_src(j)
                  if (rcss%usr(j)*sc(j) < c0) then
                     rcss%usr(j) = c0
                     rcss%uel(j) = c5*rcss%u_src(j) - c4*rcss%uer(j)
                     rcss%usl(j) = c20*(rcss%uer(j) - rcss%u_src(j))
                  else
                     rcss%usl(j) = c10*rcss%u_src(j) - c4*rcss%uer(j) &
                                 - c6*rcss%uel(j)
                     if (rcss%usl(j)*sc(j) < c0) then
                        rcss%usl(j) = c0
                        rcss%uer(j) = c5_2*rcss%u_src(j) - c3_2*rcss%uel(j)
                        rcss%usr(j) = c10_3*(rcss%u_src(j) - rcss%uel(j))
                     endif
                  endif
               endif
            endif

         endif

      enddo

      if (rcss%pc_left_bndr) then
         ! Piecewise constant reconstruction of the left boundary cell.
         rcss%uel(1) = rcss%u_src(1)
         rcss%uer(1) = rcss%u_src(1)
         rcss%usl(1) = c0
         rcss%usr(1) = c0
      else
         ! Do not treat the left boundary cell as a local extrema, but ensure
         ! that the piecewise parabola is monotonic within the cell.
         if ( (rcss%u_src(2) - rcss%uer(1)) &
             *(rcss%u_src(1) - rcss%uer(1)) > c0) then
            rcss%uel(1) = rcss%u_src(1)
            rcss%uer(1) = rcss%u_src(1)
            rcss%usl(1) = c0
            rcss%usr(1) = c0
         else
            s = c2*(rcss%u_src(3) - rcss%u_src(2)) &
                  /(rcss%rcgs%h_src(2) + rcss%rcgs%h_src(3))
            if (s > 0) then
               rcss%uer(1) = &
                  max(rcss%u_src(1), &
                      min(rcss%uel(2), &
                          rcss%u_src(1) + c1_3*s*rcss%rcgs%h_src(1)))
               rcss%uel(1) = max(min(rcss%u_src(1), c0), &
                                 c1_2*(c3*rcss%u_src(1) - rcss%uer(1)))
               rcss%uer(1) = c3*rcss%u_src(1) - c2*rcss%uel(1)
            else
               rcss%uer(1) = &
                  min(rcss%u_src(1), &
                      max(rcss%uel(2), &
                          rcss%u_src(1) + c1_3*s*rcss%rcgs%h_src(1)))
               rcss%uel(1) = c1_2*(c3*rcss%u_src(1) - rcss%uer(1))
            endif
            rcss%usl(1) = c6*rcss%u_src(1) - c4*rcss%uel(1) - c2*rcss%uer(1)
            rcss%usr(1) = c2*rcss%uel(1) + c4*rcss%uer(1) - c6*rcss%u_src(1)
         endif
      endif

      if (rcss%pc_right_bndr) then
         ! Piecewise constant reconstruction of the right boundary cell.
         rcss%uel(ns) = rcss%u_src(ns)
         rcss%uer(ns) = rcss%u_src(ns)
         rcss%usl(ns) = c0
         rcss%usr(ns) = c0
      else
         ! Do not treat the right boundary cell as a local extrema, but ensure
         ! that the piecewise parabola is monotonic within the cell.
         if ( (rcss%u_src(ns  ) - rcss%uel(ns)) &
             *(rcss%u_src(ns-1) - rcss%uel(ns)) > c0) then
            rcss%uel(ns) = rcss%u_src(ns)
            rcss%uer(ns) = rcss%u_src(ns)
            rcss%usl(ns) = c0
            rcss%usr(ns) = c0
         else
            s = c2*(rcss%u_src(ns-1) - rcss%u_src(ns-2)) &
                  /(rcss%rcgs%h_src(ns-2) + rcss%rcgs%h_src(ns-1))
            if (s > 0) then
               rcss%uel(ns) = &
                  min(rcss%u_src(ns), &
                      max(rcss%uer(ns-1), &
                          rcss%u_src(ns) - c1_3*s*rcss%rcgs%h_src(ns)))
               rcss%uer(ns) = c1_2*(c3*rcss%u_src(ns) - rcss%uel(ns))
            else
               rcss%uel(ns) = &
                  max(rcss%u_src(ns), &
                      min(rcss%uer(ns-1), &
                          rcss%u_src(ns) - c1_3*s*rcss%rcgs%h_src(ns)))
               rcss%uer(ns) = max(min(rcss%u_src(ns), c0), &
                                  c1_2*(c3*rcss%u_src(ns) - rcss%uel(ns)))
               rcss%uel(ns) = c3*rcss%u_src(ns) - c2*rcss%uer(ns)
            endif
            rcss%usl(ns) = c6*rcss%u_src(ns) - c4*rcss%uel(ns) - c2*rcss%uer(ns)
            rcss%usr(ns) = c2*rcss%uel(ns) + c4*rcss%uer(ns) - c6*rcss%u_src(ns)
         endif
      endif

   end subroutine limit_pqm_non_oscillatory_posdef

   pure subroutine polycoeff_pqm(rcss)
   ! ---------------------------------------------------------------------------
   ! Obtain coefficients for piecewise quartics from grid cell means and left
   ! and right edge and slope values.
   ! ---------------------------------------------------------------------------

      type(recon_src_struct), intent(inout) :: rcss

      integer :: j

      do j = 1, rcss%rcgs%n_src_actual
         rcss%polycoeff(1,j) = rcss%uel(j)
         rcss%polycoeff(2,j) = rcss%usl(j)
         rcss%polycoeff(3,j) = &
              c30*rcss%u_src(j) - c18*rcss%uel(j) - c12*rcss%uer(j) &
            - c9_2*rcss%usl(j) + c3_2*rcss%usr(j)
         rcss%polycoeff(4,j) = &
            - c60*rcss%u_src(j) + c32*rcss%uel(j) + c28*rcss%uer(j) &
            + c6*rcss%usl(j) - c4*rcss%usr(j)
         rcss%polycoeff(5,j) = &
              c30*rcss%u_src(j) - c15*(rcss%uel(j) + rcss%uer(j)) &
            - c5_2*(rcss%usl(j) - rcss%usr(j))
      enddo

   end subroutine polycoeff_pqm

   pure function line_intersection(pc, u, u_eps, xil, xir) result(xi)

      real(r8), dimension(2), intent(in) :: pc
      real(r8), intent(in) :: u, u_eps, xil, xir

      real(r8) :: xi

      if (abs(pc(2)) < u_eps) then
         xi = xil
      else
         xi = max(xil, min(xir, (u - pc(1))/pc(2)))
      endif

   end function line_intersection

   pure function parabola_intersection(pc, u, u_eps, xil, xir) result(xi)

      real(r8), dimension(3), intent(in) :: pc
      real(r8), intent(in) :: u, u_eps, xil, xir

      real(r8) :: xi

      real(r8) :: q, s, xi1, xi2, xim

      if (abs(pc(3)) < u_eps) then
         xi = line_intersection(pc(1:2), u, u_eps, xil, xir)
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
         xi = max(xil, min(xir, xi))
      endif

   end function parabola_intersection

   pure function quartic_intersection(pc, u, u_eps, xil, xir) result(xi)

      real(r8), dimension(5), intent(in) :: pc
      real(r8), intent(in) :: u, u_eps, xil, xir

      real(r8) :: xi

      real(r8) :: r, drdx, xi_old
      integer :: n

      if (abs(pc(4)) < u_eps .and. abs(pc(5)) < u_eps) then
         xi = parabola_intersection(pc(1:3), u, u_eps, xil, xir)
      else
         xi = c1_2*(xil + xir)
         do n = 1, 10
            r = pc(1) + (pc(2) + (pc(3) + (pc(4) + pc(5)*xi)*xi)*xi)*xi - u
            drdx = pc(2) + (c2*pc(3) + (c3*pc(4) + c4*pc(5)*xi)*xi)*xi
            xi_old = xi
            xi = max(xil, min(xir, xi_old - r/sign(max(eps, abs(drdx)), drdx)))
            if (abs(xi - xi_old) < 1.e-9_r8) return
         enddo
      endif

   end function quartic_intersection

   pure subroutine regrid_plm_method_1(rcss, u_sgn, u_edge_grd, x_edge_grd)

      type(recon_src_struct), intent(in) :: rcss
      real(r8), dimension(:), intent(in) :: u_edge_grd
      real(r8), dimension(:), intent(inout) :: x_edge_grd
      real(r8), intent(in) :: u_sgn

      real(r8) :: ue_min, ue_max, xi
      integer :: ns, ng, jg, js

      ! Number of source grid cells.
      ns = rcss%rcgs%n_src_actual

      ! Number of grid edges.
      ng = size(u_edge_grd)

      jg = 1
      do
         if ((u_edge_grd(jg) - rcss%uel(1))*u_sgn >= c0) exit
         jg = jg + 1
         if (jg > ng) return
      enddo

      js = 1
      do
         if (js + 1 > ns) exit
         ue_min = min(rcss%uer(js)*u_sgn, rcss%uel(js+1)*u_sgn)
         do
            if (u_edge_grd(jg)*u_sgn >= ue_min) exit
            xi = line_intersection(rcss%polycoeff(:,js), u_edge_grd(jg), &
                                   rcss%u_eps, c0, c1)
            x_edge_grd(jg) = rcss%rcgs%x_edge_src(js) &
                           + ( rcss%rcgs%x_edge_src(js+1) &
                             - rcss%rcgs%x_edge_src(js  ))*xi
            jg = jg + 1
            if (jg > ng) return
         enddo
         ue_max = max(rcss%uer(js)*u_sgn, rcss%uel(js+1)*u_sgn)
         do
            if (u_edge_grd(jg)*u_sgn > ue_max) exit
            x_edge_grd(jg) = rcss%rcgs%x_edge_src(js+1)
            jg = jg + 1
            if (jg > ng) return
         enddo
         js = js + 1
      enddo

      do
         if ((u_edge_grd(jg) - rcss%uer(js))*u_sgn > c0) return
         xi = line_intersection(rcss%polycoeff(:,js), u_edge_grd(jg), &
                                rcss%u_eps, c0, c1)
         x_edge_grd(jg) = rcss%rcgs%x_edge_src(js) &
                        + ( rcss%rcgs%x_edge_src(js+1) &
                          - rcss%rcgs%x_edge_src(js  ))*xi
         jg = jg + 1
         if (jg > ng) return
      enddo

   end subroutine regrid_plm_method_1

   pure subroutine regrid_ppm_method_1(rcss, u_sgn, u_edge_grd, x_edge_grd)

      type(recon_src_struct), intent(in) :: rcss
      real(r8), dimension(:), intent(in) :: u_edge_grd
      real(r8), dimension(:), intent(inout) :: x_edge_grd
      real(r8), intent(in) :: u_sgn

      real(r8) :: ue_min, ue_max, xi
      integer :: ns, ng, jg, js

      ! Number of source grid cells.
      ns = rcss%rcgs%n_src_actual

      ! Number of grid edges.
      ng = size(u_edge_grd)

      jg = 1
      do
         if ((u_edge_grd(jg) - rcss%uel(1))*u_sgn >= c0) exit
         jg = jg + 1
         if (jg > ng) return
      enddo

      js = 1
      do
         if (js + 1 > ns) exit
         ue_min = min(rcss%uer(js)*u_sgn, rcss%uel(js+1)*u_sgn)
         do
            if (u_edge_grd(jg)*u_sgn >= ue_min) exit
            xi = parabola_intersection(rcss%polycoeff(:,js), u_edge_grd(jg), &
                                       rcss%u_eps, c0, c1)
            x_edge_grd(jg) = rcss%rcgs%x_edge_src(js) &
                           + ( rcss%rcgs%x_edge_src(js+1) &
                             - rcss%rcgs%x_edge_src(js  ))*xi
            jg = jg + 1
            if (jg > ng) return
         enddo
         ue_max = max(rcss%uer(js)*u_sgn, rcss%uel(js+1)*u_sgn)
         do
            if (u_edge_grd(jg)*u_sgn > ue_max) exit
            x_edge_grd(jg) = rcss%rcgs%x_edge_src(js+1)
            jg = jg + 1
            if (jg > ng) return
         enddo
         js = js + 1
      enddo

      do
         if ((u_edge_grd(jg) - rcss%uer(js))*u_sgn > c0) return
         xi = parabola_intersection(rcss%polycoeff(:,js), u_edge_grd(jg), &
                                    rcss%u_eps, c0, c1)
         x_edge_grd(jg) = rcss%rcgs%x_edge_src(js) &
                        + ( rcss%rcgs%x_edge_src(js+1) &
                          - rcss%rcgs%x_edge_src(js  ))*xi
         jg = jg + 1
         if (jg > ng) return
      enddo

   end subroutine regrid_ppm_method_1

   pure subroutine regrid_pqm_method_1(rcss, u_sgn, u_edge_grd, x_edge_grd)

      type(recon_src_struct), intent(in) :: rcss
      real(r8), dimension(:), intent(in) :: u_edge_grd
      real(r8), dimension(:), intent(inout) :: x_edge_grd
      real(r8), intent(in) :: u_sgn

      real(r8) :: ue_min, ue_max, xi
      integer :: ns, ng, jg, js

      ! Number of source grid cells.
      ns = rcss%rcgs%n_src_actual

      ! Number of grid edges.
      ng = size(u_edge_grd)

      jg = 1
      do
         if ((u_edge_grd(jg) - rcss%uel(1))*u_sgn >= c0) exit
         jg = jg + 1
         if (jg > ng) return
      enddo

      js = 1
      do
         if (js + 1 > ns) exit
         ue_min = min(rcss%uer(js)*u_sgn, rcss%uel(js+1)*u_sgn)
         do
            if (u_edge_grd(jg)*u_sgn >= ue_min) exit
            xi = quartic_intersection(rcss%polycoeff(:,js), u_edge_grd(jg), &
                                      rcss%u_eps, c0, c1)
            x_edge_grd(jg) = rcss%rcgs%x_edge_src(js) &
                           + ( rcss%rcgs%x_edge_src(js+1) &
                             - rcss%rcgs%x_edge_src(js  ))*xi
            jg = jg + 1
            if (jg > ng) return
         enddo
         ue_max = max(rcss%uer(js)*u_sgn, rcss%uel(js+1)*u_sgn)
         do
            if (u_edge_grd(jg)*u_sgn > ue_max) exit
            x_edge_grd(jg) = rcss%rcgs%x_edge_src(js+1)
            jg = jg + 1
            if (jg > ng) return
         enddo
         js = js + 1
      enddo

      do
         if ((u_edge_grd(jg) - rcss%uer(js))*u_sgn > c0) return
         xi = quartic_intersection(rcss%polycoeff(:,js), u_edge_grd(jg), &
                                   rcss%u_eps, c0, c1)
         x_edge_grd(jg) = rcss%rcgs%x_edge_src(js) &
                        + ( rcss%rcgs%x_edge_src(js+1) &
                          - rcss%rcgs%x_edge_src(js  ))*xi
         jg = jg + 1
         if (jg > ng) return
      enddo

   end subroutine regrid_pqm_method_1

   pure subroutine regrid_plm_method_2(rcss, u_sgn, u_edge_grd, x_edge_grd)

      type(recon_src_struct), intent(in) :: rcss
      real(r8), dimension(:), intent(in) :: u_edge_grd
      real(r8), dimension(:), intent(inout) :: x_edge_grd
      real(r8), intent(in) :: u_sgn

      real(r8), dimension(3) :: pcl, pcr
      real(r8) :: umr, uml, xi, duml, dumr, uerl, uelr
      integer :: ns, ng, jg, js

      ! Number of source grid cells.
      ns = rcss%rcgs%n_src_actual

      ! Number of grid edges.
      ng = size(u_edge_grd)

      ! Find possible intersections in the first half of the first source grid
      ! cell.
      jg = 1
      do
         if ((u_edge_grd(jg) - rcss%uel(1))*u_sgn >= c0) exit
         jg = jg + 1
         if (jg > ng) return
      enddo
      js = 1
      umr =      rcss%polycoeff(1,js) &
          + c1_2*rcss%polycoeff(2,js)
      do
         if ((u_edge_grd(jg) - umr)*u_sgn > c0) exit
         xi = line_intersection(rcss%polycoeff(:,js), u_edge_grd(jg), &
                                rcss%u_eps, c0, c1_2)
         x_edge_grd(jg) = rcss%rcgs%x_edge_src(js) &
                        + ( rcss%rcgs%x_edge_src(js+1) &
                          - rcss%rcgs%x_edge_src(js  ))*xi
         jg = jg + 1
         if (jg > ng) return
      enddo

      outer: do

         ! For the current grid edge index, find the index of the first source
         ! grid cell with mid point reconstructed value larger than the grid
         ! edge value.
         do
            uml = umr
            umr =      rcss%polycoeff(1,js) &
                + c1_2*rcss%polycoeff(2,js)
            if ((u_edge_grd(jg) - umr)*u_sgn <= c0) exit
            js = js + 1
            if (js > ns) exit outer
         enddo

         ! Construct new parabolas left and right of the edge that are
         ! continuous and smooth across the edge and with the original piecewise
         ! lines left and right of the edge at the mid points of their
         ! respective grid cells.
         duml = rcss%polycoeff(2,js-1)
         dumr = rcss%polycoeff(2,js  )
         pcr(2) = (c4*(umr - uml) - duml - dumr)*rcss%rcgs%h_src(js) &
                  /(rcss%rcgs%h_src(js-1) + rcss%rcgs%h_src(js))
         pcr(1) = umr - c1_4*(dumr + pcr(2))
         if (pcr(2)*(rcss%u_src(js) - rcss%u_src(js-1)) < c0) then
            ! If the slope of the new parabolas are non-monotonic at the
            ! edge, set the edge slope to zero and enforce that the new
            ! parabolas cross the edge within the interval spanned by the
            ! edge values of the original piecewise lines. Smoothness
            ! with the original piecewise lines at grid cell mid points
            ! is then not guaranteed.
            pcr(2) = c0
            uerl = rcss%uer(js-1)
            uelr = rcss%uel(js)
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
            xi = parabola_intersection(pcl, u_edge_grd(jg), &
                                       rcss%u_eps, c1_2, c1)
            x_edge_grd(jg) = rcss%rcgs%x_edge_src(js-1) &
                           + ( rcss%rcgs%x_edge_src(js  ) &
                             - rcss%rcgs%x_edge_src(js-1))*xi
            jg = jg + 1
            if (jg > ng) return
         enddo

         ! Find all intersections with piecewise parabola in the first half of
         ! the source grid cell right of the edge.
         do
            if ((u_edge_grd(jg) - umr)*u_sgn > c0) exit
            xi = parabola_intersection(pcr, u_edge_grd(jg), &
                                       rcss%u_eps, c0, c1_2)
            x_edge_grd(jg) = rcss%rcgs%x_edge_src(js) &
                           + ( rcss%rcgs%x_edge_src(js+1) &
                             - rcss%rcgs%x_edge_src(js  ))*xi
            jg = jg + 1
            if (jg > ng) return
         enddo

      enddo outer

      ! Find possible intersections in the last half of the last source grid
      ! cell.
      js = ns
      do
         if ((u_edge_grd(jg) - rcss%uer(js))*u_sgn > c0) return
         xi = line_intersection(rcss%polycoeff(:,js), u_edge_grd(jg), &
                                rcss%u_eps, c1_2, c1)
         x_edge_grd(jg) = rcss%rcgs%x_edge_src(js) &
                        + ( rcss%rcgs%x_edge_src(js+1) &
                          - rcss%rcgs%x_edge_src(js  ))*xi
         jg = jg + 1
         if (jg > ng) return
      enddo

   end subroutine regrid_plm_method_2

   pure subroutine regrid_ppm_method_2(rcss, u_sgn, u_edge_grd, x_edge_grd)

      type(recon_src_struct), intent(in) :: rcss
      real(r8), dimension(:), intent(in) :: u_edge_grd
      real(r8), dimension(:), intent(inout) :: x_edge_grd
      real(r8), intent(in) :: u_sgn

      real(r8), dimension(3) :: pcl, pcr
      real(r8) :: umr, uml, xi, duml, dumr, uerl, uelr
      integer :: ns, ng, jg, js

      ! Number of source grid cells.
      ns = rcss%rcgs%n_src_actual

      ! Number of grid edges.
      ng = size(u_edge_grd)

      ! Find possible intersections in the first half of the first source grid
      ! cell.
      jg = 1
      do
         if ((u_edge_grd(jg) - rcss%uel(1))*u_sgn >= c0) exit
         jg = jg + 1
         if (jg > ng) return
      enddo
      js = 1
      umr =      rcss%polycoeff(1,js) &
          + c1_2*rcss%polycoeff(2,js) &
          + c1_4*rcss%polycoeff(3,js)
      do
         if ((u_edge_grd(jg) - umr)*u_sgn > c0) exit
         xi = parabola_intersection(rcss%polycoeff(:,js), u_edge_grd(jg), &
                                    rcss%u_eps, c0, c1_2)
         x_edge_grd(jg) = rcss%rcgs%x_edge_src(js) &
                        + ( rcss%rcgs%x_edge_src(js+1) &
                          - rcss%rcgs%x_edge_src(js  ))*xi
         jg = jg + 1
         if (jg > ng) return
      enddo

      outer: do

         ! For the current grid edge index, find the index of the first source
         ! grid cell with mid point reconstructed value larger than the grid
         ! edge value.
         do
            uml = umr
            umr =      rcss%polycoeff(1,js) &
                + c1_2*rcss%polycoeff(2,js) &
                + c1_4*rcss%polycoeff(3,js)
            if ((u_edge_grd(jg) - umr)*u_sgn <= c0) exit
            js = js + 1
            if (js > ns) exit outer
         enddo

         ! Construct new parabolas left and right of the edge that are
         ! continuous and smooth across the edge and with the original piecewise
         ! parabolas left and right of the edge at the mid points of their
         ! respective grid cells.
         duml = rcss%polycoeff(2,js-1) + rcss%polycoeff(3,js-1)
         dumr = rcss%polycoeff(2,js  ) + rcss%polycoeff(3,js  )
         pcr(2) = (c4*(umr - uml) - duml - dumr)*rcss%rcgs%h_src(js) &
                  /(rcss%rcgs%h_src(js-1) + rcss%rcgs%h_src(js))
         pcr(1) = umr - c1_4*(dumr + pcr(2))
         if (pcr(2)*(rcss%u_src(js) - rcss%u_src(js-1)) < c0) then
            ! If the slope of the new parabolas are non-monotonic at the
            ! edge, set the edge slope to zero and enforce that the new
            ! parabolas cross the edge within the interval spanned by the
            ! edge values of the original piecewise parabolas. Smoothness
            ! with the original piecewise parabolas at grid cell mid points
            ! is then not guaranteed.
            pcr(2) = c0
            uerl = rcss%uer(js-1)
            uelr = rcss%uel(js)
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
            xi = parabola_intersection(pcl, u_edge_grd(jg), &
                                       rcss%u_eps, c1_2, c1)
            x_edge_grd(jg) = rcss%rcgs%x_edge_src(js-1) &
                           + ( rcss%rcgs%x_edge_src(js  ) &
                             - rcss%rcgs%x_edge_src(js-1))*xi
            jg = jg + 1
            if (jg > ng) return
         enddo

         ! Find all intersections with piecewise parabola in the first half of
         ! the source grid cell right of the edge.
         do
            if ((u_edge_grd(jg) - umr)*u_sgn > c0) exit
            xi = parabola_intersection(pcr, u_edge_grd(jg), &
                                       rcss%u_eps, c0, c1_2)
            x_edge_grd(jg) = rcss%rcgs%x_edge_src(js) &
                           + ( rcss%rcgs%x_edge_src(js+1) &
                             - rcss%rcgs%x_edge_src(js  ))*xi
            jg = jg + 1
            if (jg > ng) return
         enddo

      enddo outer

      ! Find possible intersections in the last half of the last source grid
      ! cell.
      js = ns
      do
         if ((u_edge_grd(jg) - rcss%uer(js))*u_sgn > c0) return
         xi = parabola_intersection(rcss%polycoeff(:,js), u_edge_grd(jg), &
                                    rcss%u_eps, c1_2, c1)
         x_edge_grd(jg) = rcss%rcgs%x_edge_src(js) &
                        + ( rcss%rcgs%x_edge_src(js+1) &
                          - rcss%rcgs%x_edge_src(js  ))*xi
         jg = jg + 1
         if (jg > ng) return
      enddo

   end subroutine regrid_ppm_method_2

   pure subroutine regrid_pqm_method_2(rcss, u_sgn, u_edge_grd, x_edge_grd)

      type(recon_src_struct), intent(in) :: rcss
      real(r8), dimension(:), intent(in) :: u_edge_grd
      real(r8), dimension(:), intent(inout) :: x_edge_grd
      real(r8), intent(in) :: u_sgn

      real(r8), dimension(3) :: pcl, pcr
      real(r8) :: umr, uml, xi, duml, dumr, uerl, uelr
      integer :: ns, ng, jg, js

      ! Number of source grid cells.
      ns = rcss%rcgs%n_src_actual

      ! Number of grid edges.
      ng = size(u_edge_grd)

      ! Find possible intersections in the first half of the first source grid
      ! cell.
      jg = 1
      do
         if ((u_edge_grd(jg) - rcss%uel(1))*u_sgn >= c0) exit
         jg = jg + 1
         if (jg > ng) return
      enddo
      js = 1
      umr =       rcss%polycoeff(1,js) &
          + c1_2 *rcss%polycoeff(2,js) &
          + c1_4 *rcss%polycoeff(3,js) &
          + c1_8 *rcss%polycoeff(4,js) &
          + c1_16*rcss%polycoeff(5,js)
      do
         if ((u_edge_grd(jg) - umr)*u_sgn > c0) exit
         xi = quartic_intersection(rcss%polycoeff(:,js), u_edge_grd(jg), &
                                   rcss%u_eps, c0, c1_2)
         x_edge_grd(jg) = rcss%rcgs%x_edge_src(js) &
                        + ( rcss%rcgs%x_edge_src(js+1) &
                          - rcss%rcgs%x_edge_src(js  ))*xi
         jg = jg + 1
         if (jg > ng) return
      enddo

      outer: do

         ! For the current grid edge index, find the index of the first source
         ! grid cell with mid point reconstructed value larger than the grid
         ! edge value.
         do
            uml = umr
            umr =       rcss%polycoeff(1,js) &
                + c1_2 *rcss%polycoeff(2,js) &
                + c1_4 *rcss%polycoeff(3,js) &
                + c1_8 *rcss%polycoeff(4,js) &
                + c1_16*rcss%polycoeff(5,js)
            if ((u_edge_grd(jg) - umr)*u_sgn <= c0) exit
            js = js + 1
            if (js > ns) exit outer
         enddo

         ! Construct new parabolas left and right of the edge that are
         ! continuous and smooth across the edge and with the original piecewise
         ! quartics left and right of the edge at the mid points of their
         ! respective grid cells.
         duml =      rcss%polycoeff(2,js-1) +      rcss%polycoeff(3,js-1) &
              + c3_4*rcss%polycoeff(4,js-1) + c1_2*rcss%polycoeff(5,js-1)
         dumr =      rcss%polycoeff(2,js  ) +      rcss%polycoeff(3,js  ) &
              + c3_4*rcss%polycoeff(4,js  ) + c1_2*rcss%polycoeff(5,js  )
         pcr(2) = (c4*(umr - uml) - duml - dumr)*rcss%rcgs%h_src(js) &
                  /(rcss%rcgs%h_src(js-1) + rcss%rcgs%h_src(js))
         pcr(1) = umr - c1_4*(dumr + pcr(2))
         if (pcr(2)*(rcss%u_src(js) - rcss%u_src(js-1)) < c0) then
            ! If the slope of the new parabolas are non-monotonic at the
            ! edge, set the edge slope to zero and enforce that the new
            ! parabolas cross the edge within the interval spanned by the
            ! edge values of the original piecewise quartics. Smoothness
            ! with the original piecewise quartics at grid cell mid points
            ! is then not guaranteed.
            pcr(2) = c0
            uerl = rcss%uer(js-1)
            uelr = rcss%uel(js)
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
            xi = parabola_intersection(pcl, u_edge_grd(jg), &
                                       rcss%u_eps, c1_2, c1)
            x_edge_grd(jg) = rcss%rcgs%x_edge_src(js-1) &
                           + ( rcss%rcgs%x_edge_src(js  ) &
                             - rcss%rcgs%x_edge_src(js-1))*xi
            jg = jg + 1
            if (jg > ng) return
         enddo

         ! Find all intersections with piecewise parabola in the first half of
         ! the source grid cell right of the edge.
         do
            if ((u_edge_grd(jg) - umr)*u_sgn > c0) exit
            xi = parabola_intersection(pcr, u_edge_grd(jg), &
                                       rcss%u_eps, c0, c1_2)
            x_edge_grd(jg) = rcss%rcgs%x_edge_src(js) &
                           + ( rcss%rcgs%x_edge_src(js+1) &
                             - rcss%rcgs%x_edge_src(js  ))*xi
            jg = jg + 1
            if (jg > ng) return
         enddo

      enddo outer

      ! Find possible intersections in the last half of the last source grid
      ! cell.
      js = ns
      do
         if ((u_edge_grd(jg) - rcss%uer(js))*u_sgn > c0) return
         xi = quartic_intersection(rcss%polycoeff(:,js), u_edge_grd(jg), &
                                   rcss%u_eps, c1_2, c1)
         x_edge_grd(jg) = rcss%rcgs%x_edge_src(js) &
                        + ( rcss%rcgs%x_edge_src(js+1) &
                          - rcss%rcgs%x_edge_src(js  ))*xi
         jg = jg + 1
         if (jg > ng) return
      enddo

   end subroutine regrid_pqm_method_2

   ! ---------------------------------------------------------------------------
   ! Public procedures.
   ! ---------------------------------------------------------------------------

   function initialize_rcgs(rcgs) result(errstat)
   ! ---------------------------------------------------------------------------
   ! Initialize reconstruction grid data structure.
   ! ---------------------------------------------------------------------------

      type(recon_grd_struct), intent(inout) :: rcgs

      integer :: errstat

      integer :: ij_size, allocstat

      ! Check requested reconstruction method and set the required order of
      ! polynomials for the piecewise reconstruction.
      select case (rcgs%method)
         case (hor3map_pcm)
            rcgs%p_ord = 0
         case (hor3map_plm)
            rcgs%p_ord = 1
         case (hor3map_ppm)
            rcgs%p_ord = 2
            if (rcgs%left_bndr_ord == 0) then
               rcgs%left_bndr_ord = eb_ord_max_ppm
            else
               rcgs%left_bndr_ord = &
                  max(1, min(eb_ord_max_ppm, rcgs%left_bndr_ord))
            endif
            if (rcgs%right_bndr_ord == 0) then
               rcgs%right_bndr_ord = eb_ord_max_ppm
            else
               rcgs%right_bndr_ord = &
                  max(1, min(eb_ord_max_ppm, rcgs%right_bndr_ord))
            endif
         case (hor3map_pqm)
            rcgs%p_ord = 4
            if (rcgs%left_bndr_ord == 0) then
               rcgs%left_bndr_ord = eb_ord_max_pqm
            else
               rcgs%left_bndr_ord = &
                  max(1, min(eb_ord_max_pqm, rcgs%left_bndr_ord))
            endif
            if (rcgs%right_bndr_ord == 0) then
               rcgs%right_bndr_ord = eb_ord_max_pqm
            else
               rcgs%right_bndr_ord = &
                  max(1, min(eb_ord_max_pqm, rcgs%right_bndr_ord))
            endif
         case default
            errstat = hor3map_invalid_recon_method
            return
      end select

      ! Allocate data arrays.

      ij_size = (rcgs%i_ubound - rcgs%i_lbound + 1) &
               *(rcgs%j_ubound - rcgs%j_lbound + 1)

      allocate(rcgs%x_eps_data(ij_size), &
               rcgs%x_edge_src_data(rcgs%n_src+1,ij_size), &
               rcgs%h_src_data(rcgs%n_src,ij_size), &
               rcgs%hi_src_data(rcgs%n_src,ij_size), &
               rcgs%src_dst_index_data(rcgs%n_src,ij_size), &
               rcgs%n_src_actual_data(ij_size), &
               rcgs%method_actual_data(ij_size), &
               rcgs%prepared_data(ij_size), &
               stat = allocstat)
      if (allocstat /= 0) then
         errstat = hor3map_failed_to_allocate_rcgs
         return
      endif
#ifdef DEBUG
      rcgs%x_eps_data(:) = ieee_value(1._r8, ieee_signaling_nan)
      rcgs%x_edge_src_data(:,:) = ieee_value(1._r8, ieee_signaling_nan)
      rcgs%h_src_data(:,:) = ieee_value(1._r8, ieee_signaling_nan)
      rcgs%hi_src_data(:,:) = ieee_value(1._r8, ieee_signaling_nan)
      rcgs%src_dst_index_data(:,:) = - 9999
      rcgs%n_src_actual_data(:) = - 9999
      rcgs%method_actual_data(:) = - 9999
#endif

      if (rcgs%method /= hor3map_pcm) then
         allocate(rcgs%hci_src_data(rcgs%n_src,ij_size), &
                  stat = allocstat)
         if (allocstat /= 0) then
            errstat = hor3map_failed_to_allocate_rcgs
            return
         endif
#ifdef DEBUG
         rcgs%hci_src_data(:,:) = ieee_value(1._r8, ieee_signaling_nan)
#endif
      endif

      if (rcgs%method == hor3map_ppm .or. rcgs%method == hor3map_pqm) then
         allocate(rcgs%src_dst_weight_data(rcgs%n_src,ij_size), &
                  rcgs%tdecoeff_data(rcgs%p_ord+2,rcgs%n_src,ij_size), &
                  rcgs%tdscoeff_data(rcgs%p_ord+2,rcgs%n_src,ij_size), &
                  rcgs%lblu_data(rcgs%left_bndr_ord,rcgs%left_bndr_ord,&
                                 ij_size), &
                  rcgs%rblu_data(rcgs%right_bndr_ord,rcgs%right_bndr_ord, &
                                 ij_size), &
                  rcgs%left_bndr_ord_actual_data(ij_size), &
                  rcgs%right_bndr_ord_actual_data(ij_size), &
                  stat = allocstat)
         if (allocstat /= 0) then
            errstat = hor3map_failed_to_allocate_rcgs
            return
         endif
#ifdef DEBUG
         rcgs%src_dst_weight_data(:,:) = ieee_value(1._r8, ieee_signaling_nan)
         rcgs%tdecoeff_data(:,:,:) = ieee_value(1._r8, ieee_signaling_nan)
         rcgs%tdscoeff_data(:,:,:) = ieee_value(1._r8, ieee_signaling_nan)
         rcgs%lblu_data(:,:,:) = ieee_value(1._r8, ieee_signaling_nan)
         rcgs%rblu_data(:,:,:) = ieee_value(1._r8, ieee_signaling_nan)
         rcgs%left_bndr_ord_actual_data(:) = - 9999
         rcgs%right_bndr_ord_actual_data(:) = - 9999
#endif
      endif

      rcgs%prepared_data(:) = .false.
      rcgs%initialized = .true.

      errstat = hor3map_noerr

   end function initialize_rcgs

   function initialize_rcss(rcgs, rcss) result(errstat)
   ! ---------------------------------------------------------------------------
   ! Initialize reconstruction source data structure.
   ! ---------------------------------------------------------------------------

      type(recon_grd_struct), target, intent(inout) :: rcgs
      type(recon_src_struct), target, intent(inout) :: rcss

      integer :: errstat

      integer :: ij_size, allocstat

      ij_size = (rcgs%i_ubound - rcgs%i_lbound + 1) &
               *(rcgs%j_ubound - rcgs%j_lbound + 1)

      allocate(rcss%u_src_data(rcgs%n_src,ij_size), &
               rcss%u_range_data(ij_size), &
               rcss%u_eps_data(ij_size), &
               rcss%uu_eps_data(ij_size), &
               rcss%uel_data(rcgs%n_src,ij_size), &
               rcss%uer_data(rcgs%n_src,ij_size), &
               rcss%polycoeff_data(rcgs%p_ord+1,rcgs%n_src,ij_size), &
               rcss%reconstructed_data(ij_size), &
               stat = allocstat)
      if (allocstat /= 0) then
         errstat = hor3map_failed_to_allocate_rcss
         return
      endif
#ifdef DEBUG
      rcss%u_src_data(:,:) = ieee_value(1._r8, ieee_signaling_nan)
      rcss%u_range_data(:) = ieee_value(1._r8, ieee_signaling_nan)
      rcss%u_eps_data(:) = ieee_value(1._r8, ieee_signaling_nan)
      rcss%uu_eps_data(:) = ieee_value(1._r8, ieee_signaling_nan)
      rcss%uel_data(:,:) = ieee_value(1._r8, ieee_signaling_nan)
      rcss%uer_data(:,:) = ieee_value(1._r8, ieee_signaling_nan)
      rcss%polycoeff_data(:,:,:) = ieee_value(1._r8, ieee_signaling_nan)
#endif
      if (rcgs%method == hor3map_pqm) then
         allocate(rcss%usl_data(rcgs%n_src,ij_size), &
                  rcss%usr_data(rcgs%n_src,ij_size), &
                  stat = allocstat)
         if (allocstat /= 0) then
            errstat = hor3map_failed_to_allocate_rcss
            return
         endif
#ifdef DEBUG
         rcss%usl_data(:,:) = ieee_value(1._r8, ieee_signaling_nan)
         rcss%usr_data(:,:) = ieee_value(1._r8, ieee_signaling_nan)
#endif
      endif

      rcss%rcgs => rcgs
      rcss%rcss_dep_next => rcgs%rcss_dep_head
      rcgs%rcss_dep_head => rcss
      rcss%reconstructed_data(:) = .false.
      rcss%initialized = .true.

      errstat = hor3map_noerr

   end function initialize_rcss

   function initialize_rms(rcgs, rms) result(errstat)
   ! ---------------------------------------------------------------------------
   ! Initialize remapping data structure.
   ! ---------------------------------------------------------------------------

      type(recon_grd_struct), target, intent(inout) :: rcgs
      type(remap_struct), target, intent(inout) :: rms

      integer :: errstat

      integer :: ij_size, allocstat

      ij_size = (rcgs%i_ubound - rcgs%i_lbound + 1) &
               *(rcgs%j_ubound - rcgs%j_lbound + 1)

      allocate(rms%seg_int_lim_data(rcgs%n_src+rms%n_dst,ij_size), &
               rms%seg_weight_data(rcgs%n_src+rms%n_dst,ij_size), &
               rms%n_src_seg_data(rcgs%n_src,ij_size), &
               rms%seg_dst_index_data(rcgs%n_src+rms%n_dst,ij_size), &
               rms%prepared_data(ij_size), &
               stat = allocstat)
      if (allocstat /= 0) then
         errstat = hor3map_failed_to_allocate_rms
         return
      endif
#ifdef DEBUG
      rms%seg_int_lim_data(:,:) = ieee_value(1._r8, ieee_signaling_nan)
      rms%seg_weight_data(:,:) = ieee_value(1._r8, ieee_signaling_nan)
      rms%n_src_seg_data(:,:) = - 9999
      rms%seg_dst_index_data(:,:) = - 9999
#endif

      rms%rcgs => rcgs
      rms%rms_dep_next => rcgs%rms_dep_head
      rcgs%rms_dep_head => rms
      rms%prepared_data(:) = .false.
      rms%initialized = .true.

      errstat = hor3map_noerr

   end function initialize_rms

   function prepare_reconstruction(rcgs, x_edge_src, i_index, j_index) &
      result(errstat)
   ! ---------------------------------------------------------------------------
   ! Prepare reconstruction based on edge locations of source grid cells and
   ! requested reconstruction method. Reconstruction data is stored in a
   ! reconstruction grid data structure.
   ! ---------------------------------------------------------------------------

      type(recon_grd_struct), intent(inout) :: rcgs
      real(r8), dimension(:), intent(in) :: x_edge_src
      integer, optional, intent(in) :: i_index, j_index

      integer :: errstat

      integer :: n_src, j

      errstat = hor3map_noerr

      ! Check optional arguments.
      if (present(i_index)) rcgs%i_index = i_index
      if (present(j_index)) rcgs%j_index = j_index

      ! Number of source grid cells.
      n_src = size(x_edge_src) - 1

      ! Check that source grid edges are monotonically increasing or decreasing.
      if (x_edge_src(n_src+1) - x_edge_src(1) > c0) then
         do j = 1, n_src
            if (x_edge_src(j+1) < x_edge_src(j)) then
               errstat = hor3map_nonmonotonic_src_edges
               return
            endif
         enddo
      else
         do j = 1, n_src
            if (x_edge_src(j+1) > x_edge_src(j)) then
               errstat = hor3map_nonmonotonic_src_edges
               return
            endif
         enddo
      endif

      ! If needed, initialize reconstruction grid data structure.
      if (.not. rcgs%initialized) then
         rcgs%n_src = n_src
         errstat = initialize_rcgs(rcgs)
         if (errstat /= hor3map_noerr) return
      elseif (rcgs%n_src /= n_src) then
         if (rcgs%i_lbound == 1 .and. rcgs%i_ubound == 1 .and. &
             rcgs%j_lbound == 1 .and. rcgs%j_ubound == 1) then
            call free_rcgs(rcgs)
            rcgs%n_src = n_src
            errstat = initialize_rcgs(rcgs)
            if (errstat /= hor3map_noerr) return
         else
            errstat = hor3map_resizing_initialized_rcgs
            return
         endif
      endif

      ! Assign array pointers within reconstruction grid data structure.
      errstat = assign_ptr_rcgs(rcgs)
      if (errstat /= hor3map_noerr) return

      rcgs%prepared = .false.

      ! Set small value with same dimensions as edge locations.
      rcgs%x_eps = max(abs(x_edge_src(n_src+1) - x_edge_src(1)), eps)*eps

      ! Based on the requested reconstruction method, prepare the data structure
      ! for the various methods. Arrays with indices and weights are
      ! constructed that will map the source data to a continuous array of
      ! grid cells that are non-empty and with widths that will ensure
      ! condition numbers below a specified threshold of matrices in linear
      ! equation systems to be solved. If insufficient grid cells are available
      ! for the requested method, lower order methods are tried.

      rcgs%method_actual = rcgs%method

      if (rcgs%method_actual == hor3map_pqm) then
         call prepare_pqm(rcgs, x_edge_src)
         if (rcgs%n_src_actual < n_src_min_pqm) then
            rcgs%method_actual = hor3map_ppm
         endif
      endif

      if (rcgs%method_actual == hor3map_ppm) then
         call prepare_ppm(rcgs, x_edge_src)
         if (rcgs%n_src_actual < n_src_min_ppm) then
            rcgs%method_actual = hor3map_plm
         endif
      endif

      if (rcgs%method_actual == hor3map_plm) then
         call prepare_plm(rcgs, x_edge_src)
         if (rcgs%n_src_actual < n_src_min_plm) then
            rcgs%method_actual = hor3map_pcm
         endif
      endif

      if (rcgs%method_actual == hor3map_pcm) then
         call prepare_pcm(rcgs, x_edge_src)
         if (rcgs%n_src_actual == 0) then
            errstat = hor3map_src_extent_too_small
            return
         endif
      endif

      ! Set flag to indicate the reconstruction has been prepared.
      rcgs%prepared = .true.

   end function prepare_reconstruction

   function prepare_remapping(rcgs, rms, x_edge_dst, i_index, j_index) &
      result(errstat)
   ! ---------------------------------------------------------------------------
   ! Prepare remapping based on a reconstruction data structure and edge
   ! locations of destination grid cells. Remapping data is stored in a remap
   ! data structure.
   ! ---------------------------------------------------------------------------

      type(recon_grd_struct), target, intent(inout) :: rcgs
      real(r8), dimension(:), intent(in) :: x_edge_dst
      type(remap_struct), intent(inout) :: rms
      integer, optional, intent(in) :: i_index, j_index

      integer :: errstat

      real(r8), dimension(size(x_edge_dst)-1) :: h_dst
      real(r8) :: xil
      integer :: n_dst, j, js, jd, iseg

      ! Check optional arguments.
      if (present(i_index)) rcgs%i_index = i_index
      if (present(j_index)) rcgs%j_index = j_index

      ! Check that the reconstruction grid data structure has been initialized.
      if (.not. rcgs%initialized) then
         errstat = hor3map_recon_not_prepared
         return
      endif

      ! Number of destination grid cells.
      n_dst = size(x_edge_dst) - 1

      ! If needed, initialize remapping data structure.
      if (.not. rms%initialized) then
         rms%n_dst = n_dst
         errstat = initialize_rms(rcgs, rms)
         if (errstat /= hor3map_noerr) return
      elseif (.not. associated(rms%rcgs, rcgs)) then
         call free_rms(rms)
         rms%n_dst = n_dst
         errstat = initialize_rms(rcgs, rms)
         if (errstat /= hor3map_noerr) return
      elseif (rms%n_dst /= n_dst) then
         if (rcgs%i_lbound == 1 .and. rcgs%i_ubound == 1 .and. &
             rcgs%j_lbound == 1 .and. rcgs%j_ubound == 1) then
            call free_rms(rms)
            rms%n_dst = n_dst
            errstat = initialize_rms(rcgs, rms)
            if (errstat /= hor3map_noerr) return
         else
            errstat = hor3map_resizing_initialized_rms
            return
         endif
      endif

      ! Assign array pointers within reconstruction grid and source
      ! data structures.
      errstat = assign_ptr_rcgs(rms%rcgs)
      if (errstat /= hor3map_noerr) return
      errstat = assign_ptr_rms(rms)
      if (errstat /= hor3map_noerr) return

      ! Check that the reconstruction has been prepared.
      if (.not. rcgs%prepared) then
         errstat = hor3map_recon_not_prepared
         return
      endif

      rms%prepared = .false.

      ! Check for consistency between the source and destination grid range.
      if (abs(rcgs%x_edge_src(1) - x_edge_dst(1)) > rcgs%x_eps .or. &
          abs(rcgs%x_edge_src(rcgs%n_src_actual+1) - x_edge_dst(n_dst+1)) &
          > rcgs%x_eps) then
         errstat = hor3map_inconsistent_grid_range
         return
      endif

      ! Check that destination grid edges are monotonically increasing or
      ! decreasing.
      if (x_edge_dst(n_dst+1) - x_edge_dst(1) > c0) then
         do j = 1, n_dst
            if (x_edge_dst(j+1) < x_edge_dst(j)) then
               errstat = hor3map_nonmonotonic_dst_edges
               return
            endif
         enddo
      else
         do j = 1, n_dst
            if (x_edge_dst(j+1) > x_edge_dst(j)) then
               errstat = hor3map_nonmonotonic_dst_edges
               return
            endif
         enddo
      endif

      ! From edge locations, obtain destination grid cell widths.
      do j = 1, n_dst
         h_dst(j) = abs(x_edge_dst(j+1) - x_edge_dst(j))
      enddo

      ! Locate all segments that require integration of the reconstructed source
      ! data and obtain integration limits, segment weights and destination
      ! index for the accumulation of integrals.

      js = 1
      jd = 1
      do while (h_dst(jd) <= rcgs%x_eps)
         jd = jd + 1
      enddo
      iseg = 0
      rms%n_src_seg(js) = 0
      xil = c0

      if (x_edge_dst(n_dst+1) - x_edge_dst(1) > c0) then

         do
            iseg = iseg + 1
            rms%n_src_seg(js) = rms%n_src_seg(js) + 1
            rms%seg_dst_index(iseg) = jd
            if     (  abs(rcgs%x_edge_src(js+1) - x_edge_dst(jd+1)) &
                   <= rcgs%x_eps) then
               if (h_dst(jd) > rcgs%x_eps) then
                  rms%seg_int_lim(iseg) = c1
                  rms%seg_weight(iseg) = (c1 - xil)*rcgs%h_src(js)/h_dst(jd)
               else
                  rms%seg_int_lim(iseg) = xil
               endif
               if (js == rcgs%n_src_actual) exit
               xil = c0
               js = js + 1
               jd = jd + 1
               rms%n_src_seg(js) = 0
            elseif (rcgs%x_edge_src(js+1) < x_edge_dst(jd+1)) then
               rms%seg_int_lim(iseg) = c1
               rms%seg_weight(iseg) = (c1 - xil)*rcgs%h_src(js)/h_dst(jd)
               xil = c0
               js = js + 1
               rms%n_src_seg(js) = 0
            else
               if (h_dst(jd) > rcgs%x_eps) then
                  rms%seg_int_lim(iseg) = &
                     (x_edge_dst(jd+1) - rcgs%x_edge_src(js))*rcgs%hi_src(js)
                  rms%seg_weight(iseg) = (rms%seg_int_lim(iseg) - xil) &
                                         *rcgs%h_src(js)/h_dst(jd)
                  xil = rms%seg_int_lim(iseg)
               else
                  rms%seg_int_lim(iseg) = xil
               endif
               jd = jd + 1
            endif
         enddo

      else

         do
            iseg = iseg + 1
            rms%n_src_seg(js) = rms%n_src_seg(js) + 1
            rms%seg_dst_index(iseg) = jd
            if     (  abs(rcgs%x_edge_src(js+1) - x_edge_dst(jd+1)) &
                   <= rcgs%x_eps) then
               if (h_dst(jd) > rcgs%x_eps) then
                  rms%seg_int_lim(iseg) = c1
                  rms%seg_weight(iseg) = (c1 - xil)*rcgs%h_src(js)/h_dst(jd)
               else
                  rms%seg_int_lim(iseg) = xil
               endif
               if (js == rcgs%n_src_actual) exit
               xil = c0
               js = js + 1
               jd = jd + 1
               rms%n_src_seg(js) = 0
            elseif (rcgs%x_edge_src(js+1) > x_edge_dst(jd+1)) then
               rms%seg_int_lim(iseg) = c1
               rms%seg_weight(iseg) = (c1 - xil)*rcgs%h_src(js)/h_dst(jd)
               xil = c0
               js = js + 1
               rms%n_src_seg(js) = 0
            else
               if (h_dst(jd) > rcgs%x_eps) then
                  rms%seg_int_lim(iseg) = &
                     (rcgs%x_edge_src(js) - x_edge_dst(jd+1))*rcgs%hi_src(js)
                  rms%seg_weight(iseg) = (rms%seg_int_lim(iseg) - xil) &
                                         *rcgs%h_src(js)/h_dst(jd)
                  xil = rms%seg_int_lim(iseg)
               else
                  rms%seg_int_lim(iseg) = xil
               endif
               jd = jd + 1
            endif
         enddo

      endif

      rms%prepared = .true.

   end function prepare_remapping

   function reconstruct(rcgs, rcss, u_src, i_index, j_index) result(errstat)
   ! ---------------------------------------------------------------------------
   ! Carry out the piecewise polynomial reconstruction of the source data with
   ! desired limiting method and handling of boundaries.
   ! ---------------------------------------------------------------------------

      type(recon_grd_struct), target, intent(inout) :: rcgs
      real(r8), dimension(:), intent(in) :: u_src
      type(recon_src_struct), intent(inout) :: rcss
      integer, optional, intent(in) :: i_index, j_index

      integer :: errstat

      integer :: js, jd

      ! Check optional arguments.
      if (present(i_index)) rcgs%i_index = i_index
      if (present(j_index)) rcgs%j_index = j_index

      ! Check that the reconstruction grid data structure has been initialized.
      if (.not. rcgs%initialized) then
         errstat = hor3map_recon_not_prepared
         return
      endif

      ! Check consistency of number of source grid cells.
      if (size(u_src) /= rcgs%n_src) then
         errstat = hor3map_src_size_mismatch
         return
      endif

      ! If needed, initialize reconstruction source data structure.
      if (.not. rcss%initialized) then
         errstat = initialize_rcss(rcgs, rcss)
         if (errstat /= hor3map_noerr) return
      elseif (.not. associated(rcss%rcgs, rcgs)) then
         call free_rcss(rcss)
         errstat = initialize_rcss(rcgs, rcss)
         if (errstat /= hor3map_noerr) return
      endif

      ! Assign array pointers within reconstruction grid and source data
      ! structures.
      errstat = assign_ptr_rcgs(rcgs)
      if (errstat /= hor3map_noerr) return
      errstat = assign_ptr_rcss(rcss)
      if (errstat /= hor3map_noerr) return

      ! Check that the reconstruction has been prepared.
      if (.not. rcgs%prepared) then
         errstat = hor3map_recon_not_prepared
         return
      endif

      ! Copy source data array to continuous array of grid cells to be used in
      ! the reconstruction.
      if (rcgs%method_actual == hor3map_pcm .or. &
          rcgs%method_actual == hor3map_plm) then
         do js = 1, rcgs%n_src
            jd = rcgs%src_dst_index(js)
            if (jd /= 0) rcss%u_src(jd) = u_src(js)
         enddo
      else
         rcss%u_src(1:rcgs%n_src_actual) = c0
         do js = 1, rcgs%n_src
            jd = rcgs%src_dst_index(js)
            if (jd /= 0) rcss%u_src(jd) = rcss%u_src(jd) &
                                        + rcgs%src_dst_weight(js)*u_src(js)
         enddo
      endif

      ! Set source data range and associated small values.
      rcss%u_range = abs( minval(rcss%u_src(1:rcgs%n_src_actual)) &
                        - maxval(rcss%u_src(1:rcgs%n_src_actual)))
      rcss%u_eps = max(rcss%u_range, eps*eps)*eps
      rcss%uu_eps = max(rcss%u_range, eps*eps)*rcss%u_eps

      select case (rcgs%method_actual)
         case (hor3map_plm)
            select case (rcss%limiting)
               case (hor3map_no_limiting)
                  call reconstruct_plm_no_limiting(rcss)
               case (hor3map_monotonic, hor3map_non_oscillatory, &
                     hor3map_non_oscillatory_posdef)
                  call reconstruct_plm_monotonic(rcss)
               case default
                  errstat = hor3map_invalid_plm_limiting
                  return
            end select
         case (hor3map_ppm)
            call reconstruct_ppm_edge_values(rcss)
            select case (rcss%limiting)
               case (hor3map_no_limiting)
               case (hor3map_monotonic)
                  call limit_ppm_interior_monotonic(rcss)
                  call limit_ppm_boundary(rcss)
               case (hor3map_non_oscillatory)
                  call limit_ppm_interior_non_oscillatory(rcss)
                  call limit_ppm_boundary(rcss)
               case (hor3map_non_oscillatory_posdef)
                  call limit_ppm_interior_non_oscillatory(rcss)
                  call limit_ppm_boundary(rcss)
                  call limit_ppm_posdef(rcss)
               case default
                  errstat = hor3map_invalid_ppm_limiting
                  return
            end select
            call polycoeff_ppm(rcss)
         case (hor3map_pqm)
            call reconstruct_pqm_edge_slope_values(rcss)
            select case (rcss%limiting)
               case (hor3map_no_limiting)
               case (hor3map_monotonic)
                  call limit_pqm_monotonic(rcss)
               case (hor3map_non_oscillatory)
                  call limit_pqm_non_oscillatory(rcss)
               case (hor3map_non_oscillatory_posdef)
                  call limit_pqm_non_oscillatory_posdef(rcss)
               case default
                  errstat = hor3map_invalid_pqm_limiting
                  return
            end select
            call polycoeff_pqm(rcss)
      end select

      rcss%reconstructed = .true.

   end function reconstruct

   function extract_polycoeff(rcss, polycoeff, i_index, j_index) result(errstat)
   ! ---------------------------------------------------------------------------
   ! Extract polynomial coefficients of a reconstruction. For grid cells that
   ! have been merged due to potential for ill-conditioned linear systems,
   ! polynomial coefficients will be constructed that are consistent with the
   ! reconstruction of the merged cells. Near-empty grid cells are set to a
   ! constant reconstruction.
   ! ---------------------------------------------------------------------------

      type(recon_src_struct), intent(inout) :: rcss
      real(r8), dimension(:,:), intent(out) :: polycoeff
      integer, optional, intent(in) :: i_index, j_index

      integer :: errstat

      real(r8) :: xi0, q
      integer :: js0, js, jd, jd_prev

      ! Check that reconstruction source data structure has been initialized.
      if (.not. rcss%initialized) then
         errstat = hor3map_recon_not_available
         return
      endif

      ! Check optional arguments.
      if (present(i_index)) rcss%rcgs%i_index = i_index
      if (present(j_index)) rcss%rcgs%j_index = j_index

      ! Assign array pointers within data structures.
      errstat = assign_ptr_rcgs(rcss%rcgs)
      if (errstat /= hor3map_noerr) return
      errstat = assign_ptr_rcss(rcss)
      if (errstat /= hor3map_noerr) return

      ! Check that the reconstruction is available.
      if (.not. rcss%reconstructed) then
         errstat = hor3map_recon_not_available
         return
      endif

      ! Extract polynomial coefficients.

      polycoeff(:,:) = c0

      select case (rcss%rcgs%method_actual)

         case (hor3map_pcm)

            js0 = 1
            do
               jd = rcss%rcgs%src_dst_index(js0)
               if (jd == 0) then
                  polycoeff(1,js0) = rcss%u_src(1)
               else
                  polycoeff(1,js0) = rcss%u_src(jd)
                  exit
               endif
               js0 = js0 + 1
               if (js0 > rcss%rcgs%n_src) exit
            enddo
            do js = js0+1, rcss%rcgs%n_src
               jd = rcss%rcgs%src_dst_index(js)
               if (jd == 0) then
                  polycoeff(1,js) = polycoeff(1,js-1)
               else
                  polycoeff(1,js) = rcss%u_src(jd)
               endif
            enddo

         case (hor3map_plm)

            js0 = 1
            do
               jd = rcss%rcgs%src_dst_index(js0)
               if (jd == 0) then
                  polycoeff(1,js0) = rcss%polycoeff(1,1)
               else
                  polycoeff(1:2,js0) = rcss%polycoeff(1:2,1)
                  exit
               endif
               js0 = js0 + 1
               if (js0 > rcss%rcgs%n_src) exit
            enddo
            do js = js0+1, rcss%rcgs%n_src
               jd = rcss%rcgs%src_dst_index(js)
               if (jd == 0) then
                  polycoeff(1,js) = polycoeff(1,js-1) + polycoeff(2,js-1)
               else
                  polycoeff(1:2,js) = rcss%polycoeff(1:2,jd)
               endif
            enddo

         case (hor3map_ppm)

            js0 = 1
            do
               jd = rcss%rcgs%src_dst_index(js0)
               if (jd == 0) then
                  polycoeff(1,js0) = rcss%polycoeff(1,1)
               else
                  exit
               endif
               js0 = js0 + 1
               if (js0 > rcss%rcgs%n_src) exit
            enddo
            jd_prev = -1
            do js = js0, rcss%rcgs%n_src
               jd = rcss%rcgs%src_dst_index(js)
               if (jd == 0) then
                  polycoeff(1,js) = polycoeff(1,js-1) &
                                  + polycoeff(2,js-1) &
                                  + polycoeff(3,js-1)
               else
                  if (rcss%rcgs%src_dst_weight(js) == c1) then
                     polycoeff(1:3,js) = rcss%polycoeff(1:3,jd)
                  else
                     if (jd /= jd_prev) xi0 = c0
                     jd_prev = jd
                     polycoeff(1,js) = (    rcss%polycoeff(3,jd) *xi0 &
                                       +    rcss%polycoeff(2,jd))*xi0 &
                                       +    rcss%polycoeff(1,jd)
                     q =   rcss%rcgs%src_dst_weight(js)
                     polycoeff(2,js) = ( c2*rcss%polycoeff(3,jd) *xi0 &
                                       +    rcss%polycoeff(2,jd))*q
                     q = q*rcss%rcgs%src_dst_weight(js)
                     polycoeff(3,js) =      rcss%polycoeff(3,jd) *q
                     xi0 = xi0 + rcss%rcgs%src_dst_weight(js)
                  endif
               endif
            enddo

         case (hor3map_pqm)

            js0 = 1
            do
               jd = rcss%rcgs%src_dst_index(js0)
               if (jd == 0) then
                  polycoeff(1,js0) = rcss%polycoeff(1,1)
               else
                  exit
               endif
               js0 = js0 + 1
               if (js0 > rcss%rcgs%n_src) exit
            enddo
            jd_prev = -1
            do js = js0, rcss%rcgs%n_src
               jd = rcss%rcgs%src_dst_index(js)
               if (jd == 0) then
                  polycoeff(1,js) = polycoeff(1,js-1) &
                                  + polycoeff(2,js-1) &
                                  + polycoeff(3,js-1) &
                                  + polycoeff(4,js-1) &
                                  + polycoeff(5,js-1)
               else
                  if (rcss%rcgs%src_dst_weight(js) == c1) then
                     polycoeff(1:5,js) = rcss%polycoeff(1:5,jd)
                  else
                     if (jd /= jd_prev) xi0 = c0
                     jd_prev = jd
                     polycoeff(1,js) = (((    rcss%polycoeff(5,jd) *xi0 &
                                         +    rcss%polycoeff(4,jd))*xi0 &
                                         +    rcss%polycoeff(3,jd))*xi0 &
                                         +    rcss%polycoeff(2,jd))*xi0 &
                                         +    rcss%polycoeff(1,jd)
                     q =   rcss%rcgs%src_dst_weight(js)
                     polycoeff(2,js) = ((( c4*rcss%polycoeff(5,jd) *xi0 &
                                         + c3*rcss%polycoeff(4,jd))*xi0 &
                                         + c2*rcss%polycoeff(3,jd))*xi0 &
                                         +    rcss%polycoeff(2,jd))*q
                     q = q*rcss%rcgs%src_dst_weight(js)
                     polycoeff(3,js) =  (( c6*rcss%polycoeff(5,jd) *xi0 &
                                         + c3*rcss%polycoeff(4,jd))*xi0 &
                                         +    rcss%polycoeff(3,jd))*q
                     q = q*rcss%rcgs%src_dst_weight(js)
                     polycoeff(4,js) =   ( c4*rcss%polycoeff(5,jd) *xi0 &
                                         +    rcss%polycoeff(4,jd))*q
                     q = q*rcss%rcgs%src_dst_weight(js)
                     polycoeff(5,js) =        rcss%polycoeff(5,jd) *q
                     xi0 = xi0 + rcss%rcgs%src_dst_weight(js)
                  endif
               endif
            enddo

      end select

   end function extract_polycoeff

   function regrid(rcss, u_edge_grd, x_edge_grd, missing_value, &
                   i_index, j_index, regrid_method) result(errstat)
   ! ---------------------------------------------------------------------------
   ! Find grid locations where desired grid cell edge data values intersect with
   ! a reconstruction of the source data. The following methods are supported:
   ! - Method 1: Find intersections with the piecewise polynomial reconstruction
   !             provided in rcss. 
   ! - Method 2: Find intersections with parabolas constructed to be continuous
   !             and smooth across a grid cell edge and with the original
   !             piecewise polynomials adjacent of the edge at the mid points of
   !             their respective grid cells.
   ! ---------------------------------------------------------------------------

      type(recon_src_struct), intent(inout) :: rcss
      real(r8), dimension(:), intent(in) :: u_edge_grd
      real(r8), dimension(:), intent(out) :: x_edge_grd
      real(r8), intent(in) :: missing_value
      integer, optional, intent(in) :: i_index, j_index, regrid_method

      integer :: errstat, regrid_method_resolved

      real(r8) :: u_sgn

      errstat = hor3map_noerr

      ! Check that reconstruction source data structure has been initialized.
      if (.not. rcss%initialized) then
         errstat = hor3map_recon_not_available
         return
      endif

      ! Check optional arguments.
      if (present(i_index)) rcss%rcgs%i_index = i_index
      if (present(j_index)) rcss%rcgs%j_index = j_index
      if (present(regrid_method)) then
         if (regrid_method /= hor3map_regrid_method_1 .and. &
             regrid_method /= hor3map_regrid_method_2) then
            errstat = hor3map_invalid_regrid_method
            return
         endif
         regrid_method_resolved = regrid_method
      else
         regrid_method_resolved = hor3map_regrid_method_1
      endif

      ! Assign array pointers within data structures.
      errstat = assign_ptr_rcgs(rcss%rcgs)
      if (errstat /= hor3map_noerr) return
      errstat = assign_ptr_rcss(rcss)
      if (errstat /= hor3map_noerr) return

      ! Check that the reconstruction is available.
      if (.not. rcss%reconstructed) then
         errstat = hor3map_recon_not_available
         return
      endif

      ! Check grid array size consistency.
      if (size(x_edge_grd) /= size(u_edge_grd)) then
         errstat = hor3map_grd_size_mismatch
         return
      endif

      ! Initialize grid intersections as missing value.
      x_edge_grd(:) = missing_value

      ! Return in case PCM method is used.
      if (rcss%rcgs%method_actual == hor3map_pcm) return

      ! Return in case the source data range is small.
      if (rcss%u_range < eps) return

      ! To indicate monotonically increasing or decreasing source values, use
      ! the sign of the difference of the source boundary values.
      u_sgn = sign(c1, rcss%u_src(rcss%rcgs%n_src_actual) - rcss%u_src(1))

      if (regrid_method_resolved == hor3map_regrid_method_1) then
         select case (rcss%rcgs%method_actual)
            case (hor3map_plm)
               call regrid_plm_method_1(rcss, u_sgn, u_edge_grd, x_edge_grd)
            case (hor3map_ppm)
               call regrid_ppm_method_1(rcss, u_sgn, u_edge_grd, x_edge_grd)
            case (hor3map_pqm)
               call regrid_pqm_method_1(rcss, u_sgn, u_edge_grd, x_edge_grd)
         end select
      else
         select case (rcss%rcgs%method_actual)
            case (hor3map_plm)
               call regrid_plm_method_2(rcss, u_sgn, u_edge_grd, x_edge_grd)
            case (hor3map_ppm)
               call regrid_ppm_method_2(rcss, u_sgn, u_edge_grd, x_edge_grd)
            case (hor3map_pqm)
               call regrid_pqm_method_2(rcss, u_sgn, u_edge_grd, x_edge_grd)
         end select
      endif

   end function regrid

   function remap(rcss, rms, u_dst, i_index, j_index) result(errstat)
   ! ---------------------------------------------------------------------------
   ! Carry out the remapping of a piecewise polynomial reconstruction of the
   ! source data to a destination grid.
   ! ---------------------------------------------------------------------------

      type(recon_src_struct), intent(inout) :: rcss
      type(remap_struct), intent(inout) :: rms
      real(r8), dimension(:), intent(out) :: u_dst
      integer, optional, intent(in) :: i_index, j_index

      integer :: errstat

      real(r8) :: xil, xir, b1, b2, b3, b4, b5
      integer :: ns, iseg, js, jd, i_src_seg

      errstat = hor3map_noerr

      ! Check that reconstruction source data structure has been initialized.
      if (.not. rcss%initialized) then
         errstat = hor3map_recon_not_available
         return
      endif

      ! Check optional arguments.
      if (present(i_index)) rcss%rcgs%i_index = i_index
      if (present(j_index)) rcss%rcgs%j_index = j_index

      ! Check that remapping data structure has been initialized.
      if (.not. rms%initialized) then
         errstat = hor3map_remap_not_prepared
         return
      endif

      ! Check that data structures have consistent associations.
      if (.not. associated(rcss%rcgs, rms%rcgs)) then
         errstat = hor3map_inconsistent_rcgs
         return
      endif

      ! Assign array pointers within data structures.
      errstat = assign_ptr_rcgs(rcss%rcgs)
      if (errstat /= hor3map_noerr) return
      errstat = assign_ptr_rcss(rcss)
      if (errstat /= hor3map_noerr) return
      errstat = assign_ptr_rms(rms)
      if (errstat /= hor3map_noerr) return

      ! Check that the remapping has been prepared
      if (.not. rms%prepared) then
         errstat = hor3map_remap_not_prepared
         return
      endif

      ! Check that the reconstruction is available.
      if (.not. rcss%reconstructed) then
         errstat = hor3map_recon_not_available
         return
      endif

      if (size(u_dst) /= rms%n_dst) then
         errstat = hor3map_dst_size_mismatch
         return
      endif

      u_dst(:) = 0._r8
      ns = rcss%rcgs%n_src_actual
      iseg = 0

      select case (rcss%rcgs%method_actual)

         case (hor3map_pcm)

            ! Integrate the required segments of each source grid cell in
            ! succession, adding the integrals to the appropriate destination
            ! grid cells.
            do js = 1, ns
               if (rms%n_src_seg(js) == 1) then
                  iseg = iseg + 1
                  jd = rms%seg_dst_index(iseg)
                  u_dst(jd) = u_dst(jd) + rcss%u_src(js)*rms%seg_weight(iseg)
               else
                  xil = c0
                  do i_src_seg = 1, rms%n_src_seg(js)
                     iseg = iseg + 1
                     xir = rms%seg_int_lim(iseg)
                     jd = rms%seg_dst_index(iseg)
                     if (xil == xir) then
                        u_dst(jd) = rcss%u_src(js)
                     else
                        u_dst(jd) = u_dst(jd) &
                                  + rcss%u_src(js)*rms%seg_weight(iseg)
                        xil = xir
                     endif
                  enddo
               endif
            enddo

            ! Set values for any near-empty destination grid cells at the start
            ! and the end of the array.
            do jd = 1, rms%seg_dst_index(1)-1
               u_dst(jd) = rcss%u_src(1)
            enddo
            do jd = rms%seg_dst_index(iseg)+1, rms%n_dst
               u_dst(jd) = rcss%u_src(ns)
            enddo

         case (hor3map_plm)

            ! Integrate the required segments of each source grid cell in
            ! succession, adding the integrals to the appropriate destination
            ! grid cells.
            do js = 1, ns
               if (rms%n_src_seg(js) == 1) then
                  iseg = iseg + 1
                  jd = rms%seg_dst_index(iseg)
                  u_dst(jd) = u_dst(jd) + rcss%u_src(js)*rms%seg_weight(iseg)
               else
                  xil = c0
                  do i_src_seg = 1, rms%n_src_seg(js)
                     iseg = iseg + 1
                     xir = rms%seg_int_lim(iseg)
                     jd = rms%seg_dst_index(iseg)
                     if (xil == xir) then
                        u_dst(jd) = rcss%polycoeff(2,js)*xir &
                                  + rcss%polycoeff(1,js)
                     elseif (xil == c0) then
                        u_dst(jd) = u_dst(jd) &
                                  + ( c1_2*rcss%polycoeff(2,js)*xir &
                                    +      rcss%polycoeff(1,js)) &
                                    *rms%seg_weight(iseg)
                        xil = xir
                     elseif (xir == c1) then
                        b2 = c1_2*rcss%polycoeff(2,js)
                        b1 =      rcss%polycoeff(1,js) + b2
                        u_dst(jd) = u_dst(jd) &
                                  + (b2*xil + b1)*rms%seg_weight(iseg)
                        xil = xir
                     else
                        b2 = c1_2*rcss%polycoeff(2,js)
                        b1 =      rcss%polycoeff(1,js) + b2*xir
                        u_dst(jd) = u_dst(jd) &
                                  + (b2*xil + b1)*rms%seg_weight(iseg)
                        xil = xir
                     endif
                  enddo
               endif
            enddo

            ! Set values for any near-empty destination grid cells at the start
            ! and the end of the array.
            do jd = 1, rms%seg_dst_index(1)-1
               u_dst(jd) = rcss%polycoeff(1,1)
            enddo
            if (rms%seg_dst_index(iseg) < rms%n_dst) then
               jd = rms%seg_dst_index(iseg) + 1
               u_dst(jd) = rcss%polycoeff(1,ns) + rcss%polycoeff(2,ns)
               do jd = rms%seg_dst_index(iseg)+2, rms%n_dst
                  u_dst(jd) = u_dst(rms%seg_dst_index(iseg)+1)
               enddo
            endif

         case (hor3map_ppm)

            ! Integrate the required segments of each source grid cell in
            ! succession, adding the integrals to the appropriate destination
            ! grid cells.
            do js = 1, ns
               if (rms%n_src_seg(js) == 1) then
                  iseg = iseg + 1
                  jd = rms%seg_dst_index(iseg)
                  u_dst(jd) = u_dst(jd) + rcss%u_src(js)*rms%seg_weight(iseg)
               else
                  xil = c0
                  do i_src_seg = 1, rms%n_src_seg(js)
                     iseg = iseg + 1
                     xir = rms%seg_int_lim(iseg)
                     jd = rms%seg_dst_index(iseg)
                     if (xil == xir) then
                        u_dst(jd) = ( rcss%polycoeff(3,js) *xir &
                                    + rcss%polycoeff(2,js))*xir &
                                    + rcss%polycoeff(1,js)
                     elseif (xil == c0) then
                        u_dst(jd) = u_dst(jd) &
                                  + (( c1_3*rcss%polycoeff(3,js) *xir &
                                     + c1_2*rcss%polycoeff(2,js))*xir &
                                     +      rcss%polycoeff(1,js)) &
                                    *rms%seg_weight(iseg)
                        xil = xir
                     elseif (xir == c1) then
                        b3 = c1_3*rcss%polycoeff(3,js)
                        b2 = c1_2*rcss%polycoeff(2,js) + b3
                        b1 =      rcss%polycoeff(1,js) + b2
                        u_dst(jd) = u_dst(jd) &
                                  + ((b3*xil + b2)*xil + b1) &
                                    *rms%seg_weight(iseg)
                        xil = xir
                     else
                        b3 = c1_3*rcss%polycoeff(3,js)
                        b2 = c1_2*rcss%polycoeff(2,js) + b3*xir
                        b1 =      rcss%polycoeff(1,js) + b2*xir
                        u_dst(jd) = u_dst(jd) &
                                  + ((b3*xil + b2)*xil + b1) &
                                    *rms%seg_weight(iseg)
                        xil = xir
                     endif
                  enddo
               endif
            enddo

            ! Set values for any near-empty destination grid cells at the start
            ! and the end of the array.
            do jd = 1, rms%seg_dst_index(1)-1
               u_dst(jd) = rcss%polycoeff(1,1)
            enddo
            if (rms%seg_dst_index(iseg) < rms%n_dst) then
               jd = rms%seg_dst_index(iseg) + 1
               u_dst(jd) = rcss%polycoeff(1,ns) + rcss%polycoeff(2,ns) &
                         + rcss%polycoeff(3,ns)
               do jd = rms%seg_dst_index(iseg)+2, rms%n_dst
                  u_dst(jd) = u_dst(rms%seg_dst_index(iseg)+1)
               enddo
            endif

         case (hor3map_pqm)

            ! Integrate the required segments of each source grid cell in
            ! succession, adding the integrals to the appropriate destination
            ! grid cells.
            do js = 1, ns
               if (rms%n_src_seg(js) == 1) then
                  iseg = iseg + 1
                  jd = rms%seg_dst_index(iseg)
                  u_dst(jd) = u_dst(jd) + rcss%u_src(js)*rms%seg_weight(iseg)
               else
                  xil = c0
                  do i_src_seg = 1, rms%n_src_seg(js)
                     iseg = iseg + 1
                     xir = rms%seg_int_lim(iseg)
                     jd = rms%seg_dst_index(iseg)
                     if (xil == xir) then
                        u_dst(jd) = ((( rcss%polycoeff(5,js) *xir &
                                      + rcss%polycoeff(4,js))*xir &
                                      + rcss%polycoeff(3,js))*xir &
                                      + rcss%polycoeff(2,js))*xir &
                                      + rcss%polycoeff(1,js)
                     elseif (xil == c0) then
                        u_dst(jd) = u_dst(jd) &
                                  + (((( c1_5*rcss%polycoeff(5,js) *xir &
                                       + c1_4*rcss%polycoeff(4,js))*xir &
                                       + c1_3*rcss%polycoeff(3,js))*xir &
                                       + c1_2*rcss%polycoeff(2,js))*xir &
                                       +      rcss%polycoeff(1,js)) &
                                    *rms%seg_weight(iseg)
                        xil = xir
                     elseif (xir == c1) then
                        b5 = c1_5*rcss%polycoeff(5,js)
                        b4 = c1_4*rcss%polycoeff(4,js) + b5
                        b3 = c1_3*rcss%polycoeff(3,js) + b4
                        b2 = c1_2*rcss%polycoeff(2,js) + b3
                        b1 =      rcss%polycoeff(1,js) + b2
                        u_dst(jd) = u_dst(jd) &
                           + ((((b5*xil + b4)*xil + b3)*xil + b2)*xil + b1) &
                             *rms%seg_weight(iseg)
                        xil = xir
                     else
                        b5 = c1_5*rcss%polycoeff(5,js)
                        b4 = c1_4*rcss%polycoeff(4,js) + b5*xir
                        b3 = c1_3*rcss%polycoeff(3,js) + b4*xir
                        b2 = c1_2*rcss%polycoeff(2,js) + b3*xir
                        b1 =      rcss%polycoeff(1,js) + b2*xir
                        u_dst(jd) = u_dst(jd) &
                           + ((((b5*xil + b4)*xil + b3)*xil + b2)*xil + b1) &
                             *rms%seg_weight(iseg)
                        xil = xir
                     endif
                  enddo
               endif
            enddo

            ! Set values for any near-empty destination grid cells at the start
            ! and the end of the array.
            do jd = 1, rms%seg_dst_index(1)-1
               u_dst(jd) = rcss%polycoeff(1,1)
            enddo
            if (rms%seg_dst_index(iseg) < rms%n_dst) then
               jd = rms%seg_dst_index(iseg) + 1
               u_dst(jd) = rcss%polycoeff(1,ns) + rcss%polycoeff(2,ns) &
                         + rcss%polycoeff(3,ns) + rcss%polycoeff(4,ns) &
                         + rcss%polycoeff(5,ns)
               do jd = rms%seg_dst_index(iseg)+2, rms%n_dst
                  u_dst(jd) = u_dst(rms%seg_dst_index(iseg)+1)
               enddo
            endif

      end select

   end function remap

   subroutine free_rcgs(rcgs)
   ! ---------------------------------------------------------------------------
   ! Nullify pointers, deallocate arrays and reset flags.
   ! ---------------------------------------------------------------------------

      type(recon_grd_struct), intent(inout) :: rcgs

      type(recon_src_struct), pointer :: rcss_dep, rcss_dep_next
      type(remap_struct), pointer :: rms_dep, rms_dep_next

      ! Free data structures that depends on this
      ! reconstruction grid data structure.
      rcss_dep => rcgs%rcss_dep_head
      do while (associated(rcss_dep))
         rcss_dep_next => rcss_dep%rcss_dep_next
         call free_rcss(rcss_dep)
         rcss_dep => rcss_dep_next
      enddo
      rms_dep => rcgs%rms_dep_head
      do while (associated(rms_dep))
         rms_dep_next => rms_dep%rms_dep_next
         call free_rms(rms_dep)
         rms_dep => rms_dep_next
      enddo

      nullify(rcgs%x_eps, rcgs%x_edge_src, rcgs%h_src, rcgs%hi_src, &
              rcgs%src_dst_index, rcgs%n_src_actual, rcgs%method_actual, &
              rcgs%prepared, rcgs%rcss_dep_head, rcgs%rms_dep_head)
      deallocate(rcgs%x_eps_data, rcgs%x_edge_src_data, rcgs%h_src_data, &
                 rcgs%hi_src_data, rcgs%src_dst_index_data, &
                 rcgs%n_src_actual_data, rcgs%method_actual_data, &
                 rcgs%prepared_data)

      if (rcgs%method /= hor3map_pcm) then
         nullify(rcgs%hci_src)
         deallocate(rcgs%hci_src_data)
      endif

      if (rcgs%method == hor3map_ppm .or. rcgs%method == hor3map_pqm) then
         nullify(rcgs%src_dst_weight, rcgs%tdecoeff, rcgs%tdscoeff, &
                 rcgs%lblu, rcgs%rblu, rcgs%left_bndr_ord_actual, &
                 rcgs%right_bndr_ord_actual)
         deallocate(rcgs%src_dst_weight_data, rcgs%tdecoeff_data, &
                    rcgs%tdscoeff_data, rcgs%lblu_data, rcgs%rblu_data, &
                    rcgs%left_bndr_ord_actual_data, &
                    rcgs%right_bndr_ord_actual_data)
      endif

      rcgs%i_index_curr = 0
      rcgs%j_index_curr = 0
      rcgs%initialized = .false.

   end subroutine free_rcgs

   subroutine free_rcss(rcss)
   ! ---------------------------------------------------------------------------
   ! Nullify pointers, deallocate arrays and reset flags.
   ! ---------------------------------------------------------------------------

      type(recon_src_struct), intent(inout) :: rcss

      nullify(rcss%u_src, rcss%u_range, rcss%u_eps, rcss%uu_eps, &
              rcss%uel, rcss%uer, rcss%polycoeff, &
              rcss%reconstructed, rcss%rcss_dep_next)
      deallocate(rcss%u_src_data, rcss%u_range_data, rcss%u_eps_data, &
                 rcss%uu_eps_data, rcss%uel_data, rcss%uer_data, &
                 rcss%polycoeff_data, rcss%reconstructed_data)

      if (rcss%rcgs%method == hor3map_pqm) then
         nullify(rcss%usl, rcss%usr)
         deallocate(rcss%usl_data, rcss%usr_data)
      endif

      rcss%i_index_curr = 0
      rcss%j_index_curr = 0
      rcss%initialized = .false.

   end subroutine free_rcss

   subroutine free_rms(rms)
   ! ---------------------------------------------------------------------------
   ! Nullify pointers, deallocate arrays and reset flags.
   ! ---------------------------------------------------------------------------

      type(remap_struct), intent(inout) :: rms

      nullify(rms%seg_int_lim, rms%seg_weight, rms%n_src_seg, &
              rms%seg_dst_index, rms%prepared, rms%rms_dep_next)
      deallocate(rms%seg_int_lim_data, rms%seg_weight_data, &
                 rms%n_src_seg_data, rms%seg_dst_index_data, rms%prepared_data)

      rms%i_index_curr = 0
      rms%j_index_curr = 0
      rms%initialized = .false.

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
