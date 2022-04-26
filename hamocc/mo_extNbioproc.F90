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

      MODULE mo_extNbioproc
      !****************************************************************
      !
      ! MODULE mo_extNbioproc - (microbial) biological processes of the
      !                         extended nitrogen cycle 
      !
      ! j.maerz 25.04.2022
      !
      ! Purpose:
      ! --------
      ! - initialization of parameters related to the extended nitrogen cycle
      ! - representing major biological parts of the extended nitrogen cycle
      !
      ! Description:
      ! ------------
      ! The module holds the sequentially operated processes of
      ! - nitrification 
      ! - denitrification/dissimilatory nitrate reduction from NO3 to NO2  
      ! - anammox
      ! - denitrification processes from NO2 -> N2O -> N2 and DNRA 
      !   (dissimilatory nitrite reduction to ammonium)
      !
      ! The process of ammonium and nitrate uptake by phytoplankton 
      ! is handled in ocprod.
      !
      ! Ammonification (PON -> NH4) is also handled in ocprod.
      !
      ! Explicit cyanobacteria?
      !
      ! Sediment processes?
      !
      !****************************************************************
      use mo_vgrid,       only: dp_min
      use mod_xc,         only: mnproc
      use mo_control_bgc, only: io_stdo_bgc
      use mo_biomod,      only: 

      implicit none

      private

      public :: extNbioparam_init,nitrification,denit_NO3_to_NO2,&
              & anammox,denit_dnra,extN_inv_check


      CONTAINS

!==================================================================================================================================      
      subroutine extNbioparam_init()
      ! Initialization of model parameters for the extended nitrogen cycle

      end subroutine extNbioparam_init
     
!==================================================================================================================================      
      subroutine nitrification(kpie,kpje,kpke,pddpo,omask)
      ! Nitrification processes (NH4 -> NO2, NO2 -> NO3)

      integer, intent(in) :: kpie,kpje,kpke
      real,    intent(in) :: omask(kpie,kpje)       
      real,    intent(in) :: pddpo(kpie,kpje,kpke)

      !local variables
      integer :: i,j,k


      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do j = 1,kpje
       do i = 1,kpie
        do k = 1,kpke
         if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then
         
         endif
        enddo
       enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine nitrification

!==================================================================================================================================      
      subroutine denit_NO3_to_NO2(kpie,kpje,kpke,pddpo,omask)
      ! Denitrification / dissimilatory nitrate reduction (NO3 -> NO2)

      integer, intent(in) :: kpie,kpje,kpke
      real,    intent(in) :: omask(kpie,kpje)       
      real,    intent(in) :: pddpo(kpie,kpje,kpke)

      !local variables
      integer :: i,j,k


      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do j = 1,kpje
       do i = 1,kpie
        do k = 1,kpke
         if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then
         
         endif
        enddo
       enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine denit_NO3_to_NO2

!==================================================================================================================================      
      subroutine anammox(kpie,kpje,kpke,pddpo,omask)
      ! Aanammox

      integer, intent(in) :: kpie,kpje,kpke
      real,    intent(in) :: omask(kpie,kpje)       
      real,    intent(in) :: pddpo(kpie,kpje,kpke)

      !local variables
      integer :: i,j,k


      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do j = 1,kpje
       do i = 1,kpie
        do k = 1,kpke
         if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then
         
         endif
        enddo
       enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine anammox

!==================================================================================================================================      
      subroutine denit_dnra(kpie,kpje,kpke,pddpo,omask)
      ! Denitrification processes (NO2 -> N2O -> N2) and dissmilatory nitrite reduction (NO2 -> NH4)

      integer, intent(in) :: kpie,kpje,kpke
      real,    intent(in) :: omask(kpie,kpje)       
      real,    intent(in) :: pddpo(kpie,kpje,kpke)

      !local variables
      integer :: i,j,k


      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do j = 1,kpje
       do i = 1,kpie
        do k = 1,kpke
         if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then
         
         endif
        enddo
       enddo
      enddo
      !$OMP END PARALLEL DO
      end subroutine
!==================================================================================================================================      
      subroutine extN_inv_check(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,inv_message)
      integer, intent(in) :: kpie,kpje,kpke
      real,    intent(in) :: omask(kpie,kpje)       
      real,    intent(in) :: pdlxp(kpie,kpje),pdlyp(kpie,kpje),pddpo(kpie,kpje,kpke)
 
      character (len=*),intent(in) :: inv_message

#ifdef PBGC_OCNP_TIMESTEP
      if (mnproc == 1) then
        write(io_stdo_bgc,*)' '
        write(io_stdo_bgc,*)inv_message
      endif
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif
      end subroutine extN_inv_check

      END MODULE

      
