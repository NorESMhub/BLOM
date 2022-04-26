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
      use mo_control_bgc, only: io_stdo_bgc,dtb
      use mo_param1_bgc,  only: ialkali,ianh4,iano2,ian2o,iano3,idet,igasnit,iiron,ioxygen,iphosph,isco212
      use mo_carbch,      only: ocetra
      use mo_biomod,      only: riron

      implicit none

      private

      public :: extNbioparam_init,nitrification,denit_NO3_to_NO2,&
              & anammox,denit_dnra,extN_inv_check

      real   :: q10ano3denit,sc_ano3denit,Trefano3denit,rano3denit,bkano3denit, &
              & rano2anmx,q10anmx,Trefanmx,alphaanmx,bkoxanmx,bkano2anmx,bkanh4anmx 


      CONTAINS

!==================================================================================================================================      
      subroutine extNbioparam_init()
      !===========================================================================
      ! Initialization of model parameters for the extended nitrogen cycle
      
      ! === Denitrification step NO3 -> NO2:
      rano3denit    = 0.15*dtb ! Maximum growth rate (1/d)
      q10ano3denit  = 2.       ! Q10 factor for denitrification on NO3 (-)
      Trefano3denit = 10.      ! Reference temperature for denitrification on NO3 (degr C) 
      sc_ano3denit  = 0.05e6   ! Shape factor for NO3 denitrification oxygen inhibition function (m3/kmol)
      bkano3denit   = 5e-6     ! Half-saturation constant for NO3 denitrification (kmol/m3)

      ! === Anammox
      rano2anmx     = 0.05*dtb ! Maximum growth rate (1/d)
      q10anmx       = 1.6      ! Q10 factor for anammox (-)
      Trefanmx      = 10.      ! Reference temperature for anammox (degr C)
      alphaanmx     = 0.45e6   ! Shape factor for anammox oxygen inhibition function (m3/kmol)
      bkoxanmx      = 11.3e-6  ! Half-saturation constant for oxygen inhibition function (kmol/m3)
      bkano2anmx    = 5.       ! Half-saturation constant for NO2 limitation (kmol/m3)
      bkanh4anmx    = bkano2anmx * 880./1144. !Half constant function for NH4 limitation

      !===========================================================================
      end subroutine extNbioparam_init
     
!==================================================================================================================================      
      subroutine nitrification(kpie,kpje,kpke,pddpo,omask)
      ! Nitrification processes (NH4 -> NO2, NO2 -> NO3)

      integer, intent(in) :: kpie,kpje,kpke
      real,    intent(in) :: omask(kpie,kpje)       
      real,    intent(in) :: pddpo(kpie,kpje,kpke)

      !local variables
      integer :: i,j,k
      real    :: anh4Tdep,ano2Tdep     

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
      subroutine denit_NO3_to_NO2(kpie,kpje,kpke,pddpo,omask,ptho)
      ! Denitrification / dissimilatory nitrate reduction (NO3 -> NO2)

      integer, intent(in) :: kpie,kpje,kpke
      real,    intent(in) :: omask(kpie,kpje)       
      real,    intent(in) :: pddpo(kpie,kpje,kpke)
      real,    intent(in) :: ptho(kpie,kpje,kpke)

      !local variables
      integer :: i,j,k
      real    :: Tdep,O2inhib,nutlim,ano3new,ano3denit

      !$OMP PARALLEL DO PRIVATE(i,j,k,Tdep,O2inhib,nutlim,ano3new,ano3denit)
      do j = 1,kpje
       do i = 1,kpie
        do k = 1,kpke
         if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then
         
            Tdep      = q10ano3denit**((ptho(i,j,k)-Trefano3denit)/10.) 
            O2inhib   = 1. - tanh(sc_ano3denit*ocetra(i,j,k,ioxygen)) 
            nutlim    = ocetra(i,j,k,iano3)/(ocetra(i,j,k,iano3) + bkano3denit)
 
            ano3new   = ocetra(i,j,k,iano3)/(1. + rano3denit*Tdep*O2inhib*nutlim) 

            ano3denit = min(ocetra(i,j,k,iano3) - ano3new, ocetra(i,j,k,idet)*280.)

            ocetra(i,j,k,iano3)   = ocetra(i,j,k,iano3)   - ano3denit
            ocetra(i,j,k,iano2)   = ocetra(i,j,k,iano2)   + ano3denit
            ocetra(i,j,k,idet)    = ocetra(i,j,k,idet)    - ano3denit/280.
            ocetra(i,j,k,ianh4)   = ocetra(i,j,k,ianh4)   + ano3denit*16./280.
            ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212) + ano3denit*122./280.
            ocetra(i,j,k,iphosph) = ocetra(i,j,k,iphosph) + ano3denit/280.
            ocetra(i,j,k,iiron)   = ocetra(i,j,k,iiron)   + ano3denit*riron/280.
            ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali) + ano3denit*15./280.    
         endif
        enddo
       enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine denit_NO3_to_NO2

!==================================================================================================================================      
      subroutine anammox(kpie,kpje,kpke,pddpo,omask,ptho)
      ! Aanammox

      integer, intent(in) :: kpie,kpje,kpke
      real,    intent(in) :: omask(kpie,kpje)       
      real,    intent(in) :: pddpo(kpie,kpje,kpke)
      real,    intent(in) :: ptho(kpie,kpje,kpke)

      !local variables
      integer :: i,j,k
      real    :: Tdep,O2inhib,nut1lim,nut2lim,ano2new,ano2anmx


      !$OMP PARALLEL DO PRIVATE(i,j,k,Tdep,O2inhib,nut1lim,nut2lim,ano2new,ano2anmx)
      do j = 1,kpje
       do i = 1,kpie
        do k = 1,kpke
         if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then
    
           Tdep     = q10anmx**((ptho(i,j,k)-Trefanmx)/10.)         
           O2inhib  = 1. - exp(alphaanmx*(ocetra(i,j,k,ioxygen)-bkoxanmx))/(1.+ exp(alphaanmx*(ocetra(i,j,k,ioxygen)-bkoxanmx))) 
           nut1lim  = ocetra(i,j,k,iano2)/(ocetra(i,j,k,iano2)+bkano2anmx)
           nut2lim  = ocetra(i,j,k,ianh4)/(ocetra(i,j,k,ianh4)+bkanh4anmx)

           ano2new  = ocetra(i,j,k,iano2)/(1. + rano2anmx*Tdep*O2inhib*nut1lim*nut2lim) 

           ano2anmx = min(ocetra(i,j,k,iano2) - ano2new, ocetra(i,j,k,ianh4)*1144./880., ocetra(i,j,k,isco212)*1144./122., &
                         & ocetra(i,j,k,iphosph)*1144.,  ocetra(i,j,k,iiron)*1144./(riron*16.), ocetra(i,j,k,ialkali)*1144./15.)

           ocetra(i,j,k,iano2)   = ocetra(i,j,k,iano2)   - ano2anmx
           ocetra(i,j,k,ianh4)   = ocetra(i,j,k,ianh4)   - ano2anmx*880./1144.
           ocetra(i,j,k,igasnit) = ocetra(i,j,k,igasnit) + ano2anmx*864./1144.
           ocetra(i,j,k,iano3)   = ocetra(i,j,k,iano3)   + ano2anmx*280./1144.
           ocetra(i,j,k,idet)    = ocetra(i,j,k,idet)    + ano2anmx/1144.
           ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212) - ano2anmx*122/1144.
           ocetra(i,j,k,iphosph) = ocetra(i,j,k,iphosph) - ano2anmx/1144.
           ocetra(i,j,k,iiron)   = ocetra(i,j,k,iiron)   - ano2anmx*riron*16./1144.
           ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali) - ano2anmx*15./1144.

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
      ! provide inventory calculation for extended nitrogen cycle

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

!==================================================================================================================================      
      END MODULE

      
