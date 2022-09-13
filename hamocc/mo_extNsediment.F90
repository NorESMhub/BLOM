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
  use mo_param1_bgc, only: issso12,ipowaic,ipowaal,ipowaph,ipowaox,ipown2,ipowno3,ipownh4,ipown2o,ipowno2,ks 
  use mo_vgrid,      only: kbo
 
  implicit none

  private

  ! public functions
  public :: extNsediment_param_init,sed_nitrification,sed_denit_NO3_to_NO2,sed_anammox,sed_denit_DNRA

  ! public parameters
  !public ::

  ! extended nitrogen cycle sediment parameters 
  !real :: sn

  contains
  ! ================================================================================================================================
  subroutine extNsediment_param_init()

  end subroutine extNsediment_param_init

  ! ================================================================================================================================  
  subroutine sed_nitrification(j,kpie,kpje,kpke,kbnd,ptho,omask,aerob)
    integer, intent(in) :: j,kpie,kpje,kpke,kbnd
    real,    intent(in) :: omask(kpie,kpje)
    real,    intent(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)
    ! for calculation of pore water DIC and alkalinity changes [P-units]! 
    real,    intent(inout) :: aerob(kpie,ks)
   
    ! local variables
    integer :: i,k
 
  end subroutine sed_nitrification

  ! ================================================================================================================================
  subroutine sed_denit_NO3_to_NO2(j,kpie,kpje,kpke,kbnd,ptho,omask,anaerob)
    integer, intent(in) :: j,kpie,kpje,kpke,kbnd
    real,    intent(in) :: omask(kpie,kpje)
    real,    intent(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)
    ! for calculation of pore water DIC and alkalinity changes [P-units]! 
    real,    intent(inout) :: anaerob(kpie,ks)

    ! local variables
    integer :: i,k

  end subroutine sed_denit_NO3_to_NO2

  ! ================================================================================================================================
  subroutine sed_anammox(j,kpie,kpje,kpke,kbnd,ptho,omask,anaerob)
    integer, intent(in) :: j,kpie,kpje,kpke,kbnd
    real,    intent(in) :: omask(kpie,kpje)
    real,    intent(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)
    ! for calculation of pore water DIC and alkalinity changes [P-units]! 
    real,    intent(inout) :: anaerob(kpie,ks)

    ! local variables
    integer :: i,k

  end subroutine sed_anammox

  ! ================================================================================================================================
  subroutine sed_denit_DNRA(j,kpie,kpje,kpke,kbnd,ptho,omask,anaerob)
    integer, intent(in) :: j,kpie,kpje,kpke,kbnd
    real,    intent(in) :: omask(kpie,kpje)
    real,    intent(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)
    ! for calculation of pore water DIC and alkalinity changes [P-units]! 
    real,    intent(inout) :: anaerob(kpie,ks)
    
    ! local variables
    integer :: i,k

  end subroutine sed_denit_DNRA

END MODULE mo_extNsediment
