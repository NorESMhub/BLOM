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

      real
     .  total_time,
     .  total_xio_time,
     .  auxil_total_time,
     .  getfrc_total_time,
     .  tmsmt1_total_time,
     .  advdif_total_time,
     .  sfcstr_total_time,
     .  momtum_total_time,
     .  pgforc_total_time,
     .  barotp_total_time,
     .  pbtcor_total_time,
     .  convec_total_time,
     .  diapfl_total_time,
     .  thermf_total_time,
     .  mxlayr_total_time,
     .  tmsmt2_total_time,
     .  diaacc_total_time,
     .  io_total_time
c
      real*8 
     .  wtime, ! external timing function 
     .  wtimeold ! time at initialisation 
c
      common /timing/
     .  total_time,
     .  total_xio_time,
     .  auxil_total_time,
     .  getfrc_total_time,
     .  tmsmt1_total_time,
     .  advdif_total_time,
     .  sfcstr_total_time,
     .  momtum_total_time,
     .  pgforc_total_time,
     .  barotp_total_time,
     .  pbtcor_total_time,
     .  convec_total_time,
     .  diapfl_total_time,
     .  thermf_total_time,
     .  mxlayr_total_time,
     .  tmsmt2_total_time,
     .  diaacc_total_time,
     .  io_total_time,
     .  wtimeold 

