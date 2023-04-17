!>
!! @par (c) Copyright
!! This software is provided under:
!!
!! The 3-Clause BSD License
!! SPDX short identifier: BSD-3-Clause
!! See https://opensource.org/licenses/BSD-3-Clause
!! 
!! (c) Copyright 2016-2021 MPI-M, Joeran Maerz, Irene Stemmler; 
!!     first published 2020
!!
!! Redistribution and use in source and binary forms, with or without 
!! modification, are permitted provided that the following conditions are met:
!!
!! 1. Redistributions of source code must retain the above copyright notice, 
!!    this list of conditions and the following disclaimer.
!! 2. Redistributions in binary form must reproduce the above copyright notice, 
!!    this list of conditions and the following disclaimer in the documentation 
!!    and/or other materials provided with the distribution.
!! 3. Neither the name of the copyright holder nor the names of its contributors
!!    may be used to endorse or promote products derived from this software 
!!    without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
!! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
!! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
!! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
!! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
!! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.[7]
!! 
!!
!! -----------------------------------------------------------------------------
!! -----------------------------------------------------------------------------
!! @file mo_m4ago.F90
!! @brief Module for Marine Aggregates:
!!        The Microstructure, Multiscale, Mechanistic, Marine Aggregates
!!        in the Global Ocean (M4AGO) sinking scheme
!!
!! The mo_aggregates module contains routines to calculate:
!!      - aggregate properties
!!      - mean sinking velocity of aggregates
!!
!! See:
!! Maerz et al. 2020: Microstructure and composition of marine aggregates 
!!                    as co-determinants for vertical particulate organic 
!!                    carbon transfer in the global ocean. 
!!                    Biogeosciences, 17, 1765-1803, 
!!                    https://doi.org/10.5194/bg-17-1765-2020
!!
!! This module is written within the project: 
!! Multiscale Approach on the Role of Marine Aggregates (MARMA)
!! funded by the Max Planck Society (MPG)
!!
!! @author: joeran maerz (joeran.maerz@mpimet.mpg.de), MPI-M, HH
!! 2019, June, revised by Irene Stemmler (refactoring, cleaning), MPI-M, HH
!!
!! 2023 adopted to iHAMOCC by joeran maerz, UiB, Bergen 
!!
!! -----------------------------------------------------------------------------
!! -----------------------------------------------------------------------------
!! 
!!

       kav_rho_p            =  2, &
       kav_d_C              =  3, &
       kws_agg              =  4, &
       kdf_agg              =  5, &
       kstickiness_agg      =  6, &
       kb_agg               =  7, &
       kstickiness_frustule =  8, &
       kLmax_agg            =  9, &
       kdynvis              = 10, &
       kav_rhof_V           = 11, &
       kav_por_V            = 12, &
       naggdiag             = 12

     REAL, DIMENSION (:,:,:,:), ALLOCATABLE, TARGET :: aggregate_diagnostics

   CONTAINS

  !===================================================================================== m4ago_init_params
  SUBROUTINE init_m4ago_nml_params
      POM_remin_q10     = 2.1 ! Bidle et al. 2002: Regulation of Oceanic Silicon...
      POM_remin_Tref    = 10.
  END SUBROUTINE init_m4ago_nml_params

  SUBROUTINE init_m4ago_params

  END SUBROUTINE init_m4ago_params

  SUBROUTINE alloc_mem_m4ago(kpie, kpje, kpke)
     IMPLICIT NONE

     INTEGER, INTENT(in)  :: kpie                  !< 1st REAL of model grid.
     INTEGER, INTENT(in)  :: kpje                  !< 2nd REAL of model grid.
     INTEGER, INTENT(in)  :: kpke                  !< 3rd (vertical) REAL of model grid.

     ! allocate memory space for aggregate properties 
     ALLOCATE(aggregate_diagnostics(kpie, kpje, kpke, naggdiag))
     
     ! mean sinking velocity
     ALLOCATE(ws_agg(kpie,kpje,kpke))

     aggregate_diagnostics = 0.

  END SUBROUTINE alloc_mem_m4ago

  SUBROUTINE mean_aggregate_sinking_speed(kpie, kpje, kpke, kbnd, pddpo, omask, ptho, psao, ppao, prho)
     IMPLICIT NONE
     
     INTEGER, INTENT(in)  :: kpie                  !< 1st REAL of model grid.
     INTEGER, INTENT(in)  :: kpje                  !< 2nd REAL of model grid.
     INTEGER, INTENT(in)  :: kpke                  !< 3rd (vertical) REAL of model grid.
     INTEGER, INTENT(in)  :: kbnd
     REAL, INTENT(in) :: pddpo(kpie,kpje,kpke) !< size of scalar grid cell (3rd dimension) [m]
     REAL, INTENT(in) :: omask(kpie,kpje) 
     REAL, INTENT(in) :: ptho (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) !< potential temperature [deg C]
     REAL, INTENT(in) :: psao (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) !< salinity [psu.].
     REAL, INTENT(in) :: ppao (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd) !< pressure at sea level [Pa].
     REAL, INTENT(in) :: prho (kpie,kpje,kpke) !< density [g/cm3]
  END SUBROUTINE mean_aggregate_sinking_speed

END MODULE mo_m4ago
