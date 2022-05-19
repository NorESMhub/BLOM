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

      ! public functions
      public :: extNbioparam_init,nitrification,denit_NO3_to_NO2,&
              & anammox,denit_dnra,extN_inv_check
      ! public parameters
      public :: bkphyanh4,bkphyano3,bkphosph,bkiron


      real   :: q10ano3denit,sc_ano3denit,Trefano3denit,rano3denit,bkano3denit,      &
              & rano2anmx,q10anmx,Trefanmx,alphaanmx,bkoxanmx,bkano2anmx,bkanh4anmx, &
              & rano2denit,q10ano2denit,Trefano2denit,bkoxano2denit,bkano2denit,     &
              & ran2odenit,q10an2odenit,Trefan2odenit,bkoxan2odenit,bkan2odenit,     &
              & rdnra,q10dnra,Trefdnra,bkoxdnra,bkdnra,ranh4nitr,q10anh4nitr,        &
              & Trefanh4nitr,bkoxamox,bkanh4nitr,bkamoxn2o,bkamoxno2,bkyamox,        &
              & rano2nitr,q10ano2nitr,Trefano2nitr,bkoxnitr,bkano2nitr,n2omaxy,      &
              & n2oybeta,bkphyanh4,bkphyano3,bkphosph,bkiron

      real :: eps

      CONTAINS

!==================================================================================================================================      
      subroutine extNbioparam_init()
      !===========================================================================
      ! Initialization of model parameters for the extended nitrogen cycle
      
      ! Phytoplankton growth     
      bkphyanh4     = 0.1e-6    ! Half-saturation constant for NH4 uptake by bulk phytoplankton (kmol/m3)
      bkphyano3     = 0.16e-6   ! Half-saturation constant for NO3 uptake by bulk phytoplankton (kmol/m3)
      bkphosph      = 0.01e-6   ! Half-saturation constant for PO4 uptake by bulk phytoplankton (kmol/m3)
      bkiron        = bkphosph*riron ! Half-saturation constant for Fe uptake by bulk phytoplankton (kmol/m3)

      ! === Denitrification step NO3 -> NO2:
      rano3denit    = 0.15*dtb ! Maximum growth rate denitrification on NO3 at reference T (1/d -> 1/dt)
      q10ano3denit  = 2.       ! Q10 factor for denitrification on NO3 (-)
      Trefano3denit = 10.      ! Reference temperature for denitrification on NO3 (degr C) 
      sc_ano3denit  = 0.05e6   ! Shape factor for NO3 denitrification oxygen inhibition function (m3/kmol)
      bkano3denit   = 5e-6     ! Half-saturation constant for NO3 denitrification (kmol/m3)

      ! === Anammox
      rano2anmx     = 0.05*dtb ! Maximum growth rate for anammox at reference T (1/d -> 1/dt)
      q10anmx       = 1.6      ! Q10 factor for anammox (-)
      Trefanmx      = 10.      ! Reference temperature for anammox (degr C)
      alphaanmx     = 0.45e6   ! Shape factor for anammox oxygen inhibition function (m3/kmol)
      bkoxanmx      = 11.3e-6  ! Half-saturation constant for oxygen inhibition function (kmol/m3)
      bkano2anmx    = 5.e-6    ! Half-saturation constant for NO2 limitation (kmol/m3)
      bkanh4anmx    = bkano2anmx * 880./1144. !Half-saturation constant for NH4 limitation of anammox (kmol/m3)

      ! === Denitrification step NO2 -> N2O
      rano2denit    = 0.12*dtb ! Maximum growth rate denitrification on NO2 at reference T (1/d -> 1/dt) 
      q10ano2denit  = 2.0      ! Q10 factor for denitrification on NO2 (-)
      Trefano2denit = 10.      ! Reference temperature for denitrification on NO2 (degr C)
      bkoxano2denit = 2.e-6    ! Half-saturation constant for (quadratic) oxygen inhibition function of denitrification on NO2 (kmol/m3)
      bkano2denit   = 5.6e-6   ! Half-saturation constant for denitrification on NO2 (kmol/m3)

      ! === Denitrification step N2O -> N2
      ran2odenit    = 0.16*dtb ! Maximum growth rate denitrification on N2O at reference T (1/d -> 1/dt)
      q10an2odenit  = 3.       ! Q1- factor for denitrificationj on N2O (-)
      Trefan2odenit = 10.      ! Reference temperature for denitrification on N2O (degr C)
      bkoxan2odenit = 5e-6     ! Half-saturation constant for (quadratic) oxygen inhibition function of denitrification on N2O (kmol/m3)
      bkan2odenit   = 1e-6     ! Half-saturation constant for denitrification on N2O (kmol/m3)

      ! === DNRA NO2 -> NH4
      rdnra         = 0.1*dtb  ! Maximum growth rate DNRA on NO2 at reference T (1/d -> 1/dt)
      q10dnra       = 2.       ! Q10 factor for DNRA on NO2 (-)
      Trefdnra      = 10.      ! Reference temperature for DNRA (degr C)
      bkoxdnra      = 2.5e-6   ! Half saturation constant for (quadratic) oxygen inhibition function of DNRA on NO2 (kmol/m3)
      bkdnra        = 0.05e-6  ! Half-saturation constant for DNRA on NO2 (kmol/m3)

      ! === Nitrification on NH4
      ranh4nitr     = 1.*dtb   ! Maximum growth rate nitrification on NH4 at reference T (1/d -> 1/dt) 
      q10anh4nitr   = 3.3      ! Q10 factor for nitrification on NH4 (-)
      Trefanh4nitr  = 20.      ! Reference temperature for nitrification on NH4 (degr C)
      bkoxamox      = 0.333e-6 ! Half-saturation constant for oxygen limitation of nitrification on NH4 (kmol/m3)
      bkanh4nitr    = 0.133e-6 ! Half-saturation constant for nitrification on NH4 (kmol/m3)
      bkamoxn2o     = 0.453e-6 ! Half saturation constant for pathway splitting function N2O for nitrification on NH4 (kmol/m3)
      bkamoxno2     = 0.479e-6 ! Half saturation constant for pathway splitting function N2O for nitrification on NH4 (kmol/m3)
      n2omaxy       = 0.006    ! Maximum yield of OM on NH4 nitrification (-)
      n2oybeta      = 18.      ! Decay factor for inhibition function for yield during nitrification on NH4 (kmol/m3)
      bkyamox       = 0.333e-6 ! Half saturation constant for pathway splitting function OM-yield for nitrification on NH4 (kmol/m3)
  
      ! === Nitrification on NO2
      rano2nitr     = 1.54*dtb ! Maximum growth rate nitrification on NO2 at reference T (1/d -> 1/dt) 
      q10ano2nitr   = 2.7      ! Q10 factor for nitrification on NO2 (-)
      Trefano2nitr  = 20.      ! Reference temperature for nitrification on NO2 (degr C)
      bkoxnitr      = 0.788e-6 ! Half-saturation constant for oxygen limitation of nitrification on NO2 (kmol/m3)
      bkano2nitr    = 0.287e-6 ! Half-saturation constant for NO2 for nitrification on NO2 (kmol/m3)

      eps = 1e-12
      !===========================================================================
      end subroutine extNbioparam_init
     
!==================================================================================================================================      
      subroutine nitrification(kpie,kpje,kpke,pddpo,omask,ptho)
      ! Nitrification processes (NH4 -> NO2, NO2 -> NO3) accompanied
      ! by dark carbon fixation and O2-dependent N2O production 

      integer, intent(in) :: kpie,kpje,kpke
      real,    intent(in) :: omask(kpie,kpje)       
      real,    intent(in) :: pddpo(kpie,kpje,kpke)
      real,    intent(in) :: ptho(kpie,kpje,kpke)

      !local variables
      integer :: i,j,k
      real    :: Tdepanh4,O2limanh4,nut1lim,anh4new,potdnh4amox,fdetamox,fno2,fn2o,ftotnh4 
      real    :: Tdepano2,O2limano2,nut2lim,ano2new,potdno2nitr,fdetnitr,fno3,ftotno2
      real    :: amoxfrac,nitrfrac,totd,amox,nitr
 

      !$OMP PARALLEL DO PRIVATE(i,j,k,Tdepanh4,O2limanh4,nut1lim,anh4new,potdnh4amox,fdetamox,fno2,fn2o,ftotnh4, & 
      !$OMP                     Tdepano2,O2limano2,nut2lim,ano2new,potdno2nitr,fdetnitr,fno3,ftotno2,amoxfrac,   &
      !$OMP                     nitrfrac,totd,amox,nitr)

      do j = 1,kpje
       do i = 1,kpie
        do k = 1,kpke
         if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then
       
            ! Ammonium oxidation step of nitrification
            Tdepanh4    = q10anh4nitr**((ptho(i,j,k)-Trefanh4nitr)/10.) 
            O2limanh4   = ocetra(i,j,k,ioxygen)/(ocetra(i,j,k,ioxygen) + bkoxamox)
            nut1lim     = ocetra(i,j,k,ianh4)/(ocetra(i,j,k,ianh4) + bkanh4nitr)
            anh4new     = ocetra(i,j,k,ianh4)/(1. + ranh4nitr*Tdepanh4*O2limanh4*nut1lim)
            potdnh4amox = ocetra(i,j,k,ianh4) - anh4new
            
            ! pathway splitting functions according to Goreau 1980
            fn2o     = 1. - ocetra(i,j,k,ioxygen)/(ocetra(i,j,k,ioxygen) + bkamoxn2o)
            fno2     = ocetra(i,j,k,ioxygen)/(ocetra(i,j,k,ioxygen) + bkamoxno2)
            fdetamox = n2omaxy*2.*(1. + n2oybeta)*ocetra(i,j,k,ioxygen)*bkyamox &
                     & /(ocetra(i,j,k,ioxygen)**2 + 2.*ocetra(i,j,k,ioxygen)*bkyamox + bkyamox**2)

            ! normalization of pathway splitting functions to sum=1
            ftotnh4  = fn2o + fno2 + fdetamox + eps
            fn2o     = fn2o/ftotnh4
            fno2     = fno2/ftotnh4
            fdetamox = 1. - (fn2o + fno2)

            ! NO2 oxidizing step of nitrification
            Tdepano2    = q10ano2nitr**((ptho(i,j,k)-Trefano2nitr)/10.) 
            O2limano2   = ocetra(i,j,k,ioxygen)/(ocetra(i,j,k,ioxygen) + bkoxnitr)
            nut2lim     = ocetra(i,j,k,iano2)/(ocetra(i,j,k,iano2) + bkano2nitr)
            ano2new     = ocetra(i,j,k,iano2)/(1. + rano2nitr*Tdepano2*O2limano2*nut2lim)
            potdno2nitr = ocetra(i,j,k,iano2) - ano2new

            ! pathway splitting functions for NO2 nitrification - assuming to be the same as for NH4
            fno3     = fno2 + fn2o! no N2O prod in this step - NO2 enters instead NO3
            fdetnitr = fdetamox

            ! normalization of pathway splitting functions for NO2 nitrification
            ftotno2  = fno2 + fdetamox + eps
            fno3     = fno3/ftotno2
            fdetnitr = 1. - fno3

            ! limitation of the two processes through available nutrients, etc.
            totd     = potdnh4amox + potdno2nitr
            amoxfrac =  potdnh4amox/(totd + eps)
            nitrfrac = 1. - amoxfrac
            totd     = max(0.,                                                                                                     &
                     &   min(totd,                                                                                                 &
                     &       ocetra(i,j,k,ianh4)/(amoxfrac + fdetamox*nitrfrac + eps),                                             & ! ammonium
                     &       ocetra(i,j,k,isco212)/((122./16.)*(fdetamox*amoxfrac + fdetnitr*nitrfrac) + eps),                     & ! CO2
                     &       ocetra(i,j,k,iphosph)/((fdetamox*amoxfrac + fdetnitr*nitrfrac)/16. + eps),                            & ! PO4
                     &       ocetra(i,j,k,iiron)/((fdetamox*amoxfrac + fdetnitr*nitrfrac)/(16.*riron) + eps),                      & ! Fe
                     &       ocetra(i,j,k,ioxygen)                                                                                 &
                     &       /((1.5*fno2 + fn2o - 140./16.*fdetamox)*amoxfrac + (0.5*fno3 + 140./16.*fdetnitr)*nitrfrac +eps),     & ! O2
                     &       ocetra(i,j,k,ialkali)                                                                                 &
                     &       /((2.*fno2 + fn2o + 15./16.*fdetamox)*amoxfrac + (15./16.*fdetnitr)*nitrfrac + eps)))                   ! alkalinity
            amox     = amoxfrac*totd 
            nitr     = nitrfrac*totd

            ocetra(i,j,k,ianh4)   = ocetra(i,j,k,ianh4) - amox - fdetnitr*nitr
            ocetra(i,j,k,ian2o)   = ocetra(i,j,k,ian2o) + 0.5*fn2o*amox
            ocetra(i,j,k,iano2)   = ocetra(i,j,k,iano2) + fno2*amox - nitr
            ocetra(i,j,k,iano3)   = ocetra(i,j,k,iano3) + nitr
            ocetra(i,j,k,idet)    = ocetra(i,j,k,idet)  + fdetamox/16.*amox + fdetnitr/16.*nitr
            ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212) - 122./16.*fdetamox*amox - 122./16.*fdetnitr*nitr
            ocetra(i,j,k,iphosph) = ocetra(i,j,k,iphosph) - fdetamox/16.*amox - fdetnitr/16.*nitr
            ocetra(i,j,k,iiron)   = ocetra(i,j,k,iiron)   - riron/16.*fdetamox*amox - riron/16.*fdetnitr*nitr
            ocetra(i,j,k,ioxygen) = ocetra(i,j,k,ioxygen) - (1.5*fno2 + fn2o - 140./16.*fdetamox)*amox   &
                                  &                       - (0.5*fno3 - 140./16.*fdetnitr)*nitr
            ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali) - (2.*fno2 + fn2o + 15./16.*fdetamox)*amox - 15./16.*fdetnitr*nitr
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

            ano3denit = max(0.,min(ocetra(i,j,k,iano3) - ano3new, ocetra(i,j,k,idet)*280.))

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

           ano2anmx = max(0.,min(ocetra(i,j,k,iano2) - ano2new, ocetra(i,j,k,ianh4)*1144./880., ocetra(i,j,k,isco212)*1144./122., &
                         & ocetra(i,j,k,iphosph)*1144.,  ocetra(i,j,k,iiron)*1144./(riron*16.), ocetra(i,j,k,ialkali)*1144./15.))

           ocetra(i,j,k,iano2)   = ocetra(i,j,k,iano2)   - ano2anmx
           ocetra(i,j,k,ianh4)   = ocetra(i,j,k,ianh4)   - ano2anmx*880./1144.
           ocetra(i,j,k,igasnit) = ocetra(i,j,k,igasnit) + ano2anmx*864./1144.
           ocetra(i,j,k,iano3)   = ocetra(i,j,k,iano3)   + ano2anmx*280./1144.
           ocetra(i,j,k,idet)    = ocetra(i,j,k,idet)    + ano2anmx/1144.
           ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212) - ano2anmx*122./1144.
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
      subroutine denit_dnra(kpie,kpje,kpke,pddpo,omask,ptho)
      ! Denitrification processes (NO2 -> N2O -> N2) and dissmilatory nitrite reduction (NO2 -> NH4)

      integer, intent(in) :: kpie,kpje,kpke
      real,    intent(in) :: omask(kpie,kpje)       
      real,    intent(in) :: pddpo(kpie,kpje,kpke)
      real,    intent(in) :: ptho(kpie,kpje,kpke)

      !local variables
      integer :: i,j,k
      real    :: Tdepano2,O2inhibano2,nutlimano2,detlimano2,rpotano2denit,ano2denit
      real    :: Tdepdnra,O2inhibdnra,nutlimdnra,detlimdnra,rpotano2dnra,ano2dnra
      real    :: fdenit,fdnra,potano2new,potdano2,potddet,fdetano2denit,fdetan2odenit,fdetdnra  
      real    :: Tdepan2o,O2inhiban2o,nutliman2o,detliman2o,an2onew,an2odenit

      !$OMP PARALLEL DO PRIVATE(i,j,k,Tdepano2,O2inhibano2,nutlimano2,detlimano2,ano2denit,    &
      !$OMP                     Tdepan2o,O2inhiban2o,nutliman2o,detliman2o,an2onew,an2odenit,  &
      !$OMP                     rpotano2denit,rpotano2dnra,                                    &
      !$OMP                     fdenit,fdnra,potano2new,potdano2,potddet,fdetano2denit,        &
      !$OMP                     fdetan2odenit,fdetdnra,                                        &
      !$OMP                     Tdepdnra,O2inhibdnra,nutlimdnra,detlimdnra,ano2dnra)

      do j = 1,kpje
       do i = 1,kpie
        do k = 1,kpke
         if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then
           
           ! denitrification on NO2
           Tdepano2    =  q10ano2denit**((ptho(i,j,k)-Trefano2denit)/10.) 
           O2inhibano2 = 1. - ocetra(i,j,k,ioxygen)**2/(ocetra(i,j,k,ioxygen)**2 + bkoxano2denit**2) 
           nutlimano2  = ocetra(i,j,k,iano2)/(ocetra(i,j,k,iano2) + bkano2denit)
           rpotano2denit = max(0.,rano2denit*Tdepano2*O2inhibano2*nutlimano2) ! potential rate of denit
           
           ! DNRA on NO2
           Tdepdnra    = q10dnra**((ptho(i,j,k)-Trefdnra)/10.) 
           O2inhibdnra = 1. - ocetra(i,j,k,ioxygen)**2/(ocetra(i,j,k,ioxygen)**2 + bkoxdnra**2) 
           nutlimdnra  = ocetra(i,j,k,iano2)/(ocetra(i,j,k,iano2) + bkdnra)
           rpotano2dnra = max(0.,rdnra*Tdepdnra*O2inhibdnra*nutlimdnra) ! pot. rate of dnra

           ! === limitation due to NO2:
           ! fraction on potential change of NO2:
           fdenit = rpotano2denit/(rpotano2denit + rpotano2dnra + eps)
           fdnra  = 1. - fdenit

           ! potential new conc of NO2 due to denitrification and DNRA
           potano2new = ocetra(i,j,k,iano2)/(1. + rpotano2denit + rpotano2dnra)
           potdano2   = max(0.,min(ocetra(i,j,k,iano2), ocetra(i,j,k,iano2) - potano2new))
           
           ! potential fractional change
           ano2denit  = fdenit * potdano2   
           ano2dnra   = fdnra  * potdano2

           ! === denitrification on N2O
           Tdepan2o    = q10an2odenit**((ptho(i,j,k)-Trefan2odenit)/10.) 
           O2inhiban2o = 1. - ocetra(i,j,k,ioxygen)**2/(ocetra(i,j,k,ioxygen)**2 + bkoxan2odenit**2) 
           nutliman2o  = ocetra(i,j,k,ian2o)/(ocetra(i,j,k,ian2o) + bkan2odenit)   
           an2onew     = ocetra(i,j,k,ian2o)/(1. + ran2odenit*Tdepan2o*O2inhiban2o*nutliman2o)  
           an2odenit   = max(0.,min(ocetra(i,j,k,ian2o),ocetra(i,j,k,ian2o) - an2onew))

           ! limitation of processes due to detritus
           potddet       = 1./280.*(ano2denit + an2odenit) + 1./(93. + 1./3.)*ano2dnra  ! P units              
           fdetano2denit = 1./280.*ano2denit/(potddet + eps)
           fdetan2odenit = 1./280.*an2odenit/(potddet + eps)
           fdetdnra      = 1. - fdetano2denit - fdetan2odenit 
           potddet       = max(0.,min(potddet,ocetra(i,j,k,idet))) 
       
           ! change of NO2 and N2O in N units
           ano2denit     = fdetano2denit*280.*potddet
           an2odenit     = fdetan2odenit*280.*potddet
           ano2dnra      = fdetdnra * (93. + 1./3.)*potddet

           ! change in tracer concentrations due to denit (NO2->N2O->N2) and DNRA (NO2->NH4)
           ocetra(i,j,k,iano2)   = ocetra(i,j,k,iano2)   - ano2denit - ano2dnra
           ocetra(i,j,k,ian2o)   = ocetra(i,j,k,ian2o)   - an2odenit + 0.5*ano2denit
           ocetra(i,j,k,igasnit) = ocetra(i,j,k,igasnit) + an2odenit
           ocetra(i,j,k,ianh4)   = ocetra(i,j,k,ianh4)   + 16./280. * (ano2denit+an2odenit) + (109.+1./3.)/(93.+1./3.)*ano2dnra
           ocetra(i,j,k,idet)    = ocetra(i,j,k,idet)    - (ano2denit + an2odenit)/280. - ano2dnra/(93.+1./3.)
           ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212) + 122./280.*(ano2denit + an2odenit) + 122./(93.+1./3.) * ano2dnra
           ocetra(i,j,k,iphosph) = ocetra(i,j,k,iphosph) + (ano2denit + an2odenit)/280. + ano2dnra/(93.+1./3.)
           ocetra(i,j,k,iiron)   = ocetra(i,j,k,iiron)   + riron/280.*(ano2denit + an2odenit) + riron/(93.+1./3.) * ano2dnra
           ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali) + (295.*ano2denit + 15.*an2odenit)/280. &
                                 &                       + (201.+1./3.)/(93.+1./3.) * ano2dnra
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

      
