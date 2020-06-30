! ------------------------------------------------------------------------------
! Copyright (C) 2006-2020 Mats Bentsen, Alok Kumar Gupta
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

      subroutine initime
c
c --- ------------------------------------------------------------------
c --- Initializes variables used for timing. Must be called before use
c --- of 'gettime'
c --- ------------------------------------------------------------------
c
      use mod_xc, only: mnproc
      implicit none
c
#include "timing.h"
c
      if (mnproc.eq.1) wtimeold=wtime()
c
      total_time       =0.
      total_xio_time   =0.
      auxil_total_time =0.
      getfrc_total_time=0.
      tmsmt1_total_time=0.
      advdif_total_time=0.
      sfcstr_total_time=0.
      momtum_total_time=0.
      pgforc_total_time=0.
      barotp_total_time=0.
      pbtcor_total_time=0.
      convec_total_time=0.
      diapfl_total_time=0.
      thermf_total_time=0.
      mxlayr_total_time=0.
      tmsmt2_total_time=0.
      diaacc_total_time=0.
      io_total_time    =0.
c
      return
      end
