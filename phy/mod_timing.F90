! ------------------------------------------------------------------------------
! Copyright (C) 2006-2022 Mats Bentsen, Alok Kumar Gupta
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

module mod_timing
! ------------------------------------------------------------------------------
! This module contains variables and procedures related to timing.
! ------------------------------------------------------------------------------

   use mod_types, only: r8
   use mod_xc, only: mnproc

   implicit none

   private

   real(r8) :: &
      total_time, &
      total_xio_time, &
      auxil_total_time, &
      getfrc_total_time, &
      tmsmt1_total_time, &
      advdif_total_time, &
      sfcstr_total_time, &
      momtum_total_time, &
      pgforc_total_time, &
      barotp_total_time, &
      pbcor2_total_time, &
      convec_total_time, &
      diapfl_total_time, &
      thermf_total_time, &
      mxlayr_total_time, &
      tmsmt2_total_time, &
      diaacc_total_time, &
      io_total_time, &
      wtimeold             ! time at initialisation 

   real(r8), external :: &
      wtime                ! external timing function 

   public :: total_time, total_xio_time, auxil_total_time, getfrc_total_time, &
             tmsmt1_total_time, advdif_total_time, sfcstr_total_time, &
             momtum_total_time, pgforc_total_time, barotp_total_time, &
             pbcor2_total_time, convec_total_time, diapfl_total_time, &
             thermf_total_time, mxlayr_total_time, tmsmt2_total_time, &
             diaacc_total_time, io_total_time, init_timing, get_time

contains

   subroutine init_timing
   ! ---------------------------------------------------------------------------
   ! Initializes variables used for timing. Must be called before call to
   ! 'get_time'.
   ! ---------------------------------------------------------------------------

      if (mnproc == 1) wtimeold = wtime()

      total_time        = 0._r8
      total_xio_time    = 0._r8
      auxil_total_time  = 0._r8
      getfrc_total_time = 0._r8
      tmsmt1_total_time = 0._r8
      advdif_total_time = 0._r8
      sfcstr_total_time = 0._r8
      momtum_total_time = 0._r8
      pgforc_total_time = 0._r8
      barotp_total_time = 0._r8
      pbcor2_total_time = 0._r8
      convec_total_time = 0._r8
      diapfl_total_time = 0._r8
      thermf_total_time = 0._r8
      mxlayr_total_time = 0._r8
      tmsmt2_total_time = 0._r8
      diaacc_total_time = 0._r8
      io_total_time     = 0._r8

   end subroutine init_timing

   real(r8) function get_time()
   ! ---------------------------------------------------------------------------
   ! Return time in seconds since last call to either init_timing or get_time.
   ! ---------------------------------------------------------------------------

      if (mnproc == 1) then
         get_time = wtime() - wtimeold
         wtimeold = get_time + wtimeold
      else
         get_time = 1._r8
      endif

   end function get_time

end module mod_timing
