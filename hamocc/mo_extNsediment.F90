! Copyright (C) 2022  j. maerz
!
! This file is part of BLOM/iHAMOCC.
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
! along with BLOM. If not, see https://www.gnu.org/licenses/.

MODULE mo_extNsediment
!**********************************************************************
! 
! MODULE mo_extNsediment - extended nitrogen cycle processes 
!                          in the sediment
!
! j.maerz 13.09.2022
! 
! Pupose:
! -------
!   - initialization of sediment related parameters of the 
!     extended nitrogen cycle 
!   - representation of microbial processes
!
! Description:
! ------------
! The module holds the sequentially operated processes of:
!   - nitrification 
!   - denitrification/dissimilatory nitrate reduction from NO3 to NO2  
!   - anammox
!   - denitrification processes from NO2 -> N2O -> N2 and DNRA 
!     (dissimilatory nitrite reduction to ammonium)
!
! The process of ammonification in the sediment for the extended 
! nitrogen cycle is handled inside powach.F90. 
!
!**********************************************************************
  implicit none
  private

  ! public functions
  public :: extNsediment_param_init,sed_nitrification,sed_denit_NO3_to_NO2,sed_anammox,sed_denit_DNRA

  ! public parameters
  !public ::

  ! extended nitrogen cycle sediment parameters 
  real :: sn

  contains
  ! ================================================================================================================================
  subroutine extNsediment_param_init()

  end subroutine extNsediment_param_init

  ! ================================================================================================================================  
  subroutine sed_nitrification()

  end subroutine sed_nitrification

  ! ================================================================================================================================
  subroutine sed_denit_NO3_to_NO2()

  end subroutine sed_denit_NO3_to_NO2

  ! ================================================================================================================================
  subroutine sed_anammox()

  end subroutine sed_anammox

  ! ================================================================================================================================
  subroutine sed_denit_DNRA()

  end subroutine sed_denit_DNRA

END MODULE mo_extNsediment
