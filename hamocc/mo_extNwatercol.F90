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

      MODULE mo_extNwatercol
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
      ! The respective sediment processes are handled in:
      !     - powach.F90 and
      !     - mo_extNsediment.F90
      !
      !****************************************************************
      use mo_vgrid,       only: dp_min
      use mod_xc,         only: mnproc
      use mo_control_bgc, only: io_stdo_bgc,dtb
      use mo_param1_bgc,  only: ialkali,ianh4,iano2,ian2o,iano3,idet,igasnit,iiron,ioxygen,iphosph,isco212
      use mo_carbch,      only: ocetra
      use mo_biomod,      only: riron,rnit,rcar,rnoi, nitr_NH4,nitr_NO2,nitr_N2O_prod,nitr_NH4_OM,nitr_NO2_OM,denit_NO3,denit_NO2, &
                              & denit_N2O,DNRA_NO2,anmx_N2_prod,anmx_OM_prod

      implicit none

      private

      ! public functions
      public :: extNwatercol_param_init,nitrification,denit_NO3_to_NO2,&
              & anammox,denit_dnra,extN_inv_check,extNwatercol_param_update,extNwatercol_param_write

      ! public parameters for primary production
      public :: bkphyanh4,bkphyano3,bkphosph,bkiron,ro2utammo

      ! Public parameters for extended nitrogen cycle in the sediment.
      ! The basic idea is that we have the same temperature dependence
      ! and same nutrient sensitivities, 
      ! while only the rates vary between sediment and water column
      ! (Thus far, we keep the rates public in order to enable to write them to the log in beleg_parm)
      public ::  q10ano3denit,sc_ano3denit,Trefano3denit,rano3denit,bkano3denit,     &
              & rano2anmx,q10anmx,Trefanmx,alphaanmx,bkoxanmx,bkano2anmx,bkanh4anmx, &
              & rano2denit,q10ano2denit,Trefano2denit,bkoxano2denit,bkano2denit,     &
              & ran2odenit,q10an2odenit,Trefan2odenit,bkoxan2odenit,bkan2odenit,     &
              & rdnra,q10dnra,Trefdnra,bkoxdnra,bkdnra,ranh4nitr,q10anh4nitr,        &
              & Trefanh4nitr,bkoxamox,bkanh4nitr,bkamoxn2o,bkyamox,                  &
              & rano2nitr,q10ano2nitr,Trefano2nitr,bkoxnitr,bkano2nitr,n2omaxy,      &
              & n2oybeta,NOB2AOAy,bn2o,mufn2o,   &
              & rc2n,ro2nnit,rnoxp,rnoxpi,rno2anmx,rno2anmxi,rnh4anmx,     &
              & rnh4anmxi,rno2dnra,rno2dnrai,rnh4dnra,rnh4dnrai,rnm1  


      real   :: q10ano3denit,sc_ano3denit,Trefano3denit,rano3denit,bkano3denit,      &
              & rano2anmx,q10anmx,Trefanmx,alphaanmx,bkoxanmx,bkano2anmx,bkanh4anmx, &
              & rano2denit,q10ano2denit,Trefano2denit,bkoxano2denit,bkano2denit,     &
              & ran2odenit,q10an2odenit,Trefan2odenit,bkoxan2odenit,bkan2odenit,     &
              & rdnra,q10dnra,Trefdnra,bkoxdnra,bkdnra,ranh4nitr,q10anh4nitr,        &
              & Trefanh4nitr,bkoxamox,bkanh4nitr,bkamoxn2o,bkyamox,                  &
              & rano2nitr,q10ano2nitr,Trefano2nitr,bkoxnitr,bkano2nitr,n2omaxy,      &
              & n2oybeta,bkphyanh4,bkphyano3,bkphosph,bkiron,NOB2AOAy,bn2o,mufn2o!,bkamoxno2, 

      real   :: rc2n,ro2utammo,ro2nnit,rnoxp,rnoxpi,rno2anmx,rno2anmxi,rnh4anmx,     &
              & rnh4anmxi,rno2dnra,rno2dnrai,rnh4dnra,rnh4dnrai,rnm1  

      real :: eps,minlim

      CONTAINS

!==================================================================================================================================      
      subroutine extNwatercol_param_init()
      !===========================================================================
      ! Initialization of model parameters for the extended nitrogen cycle
      rc2n          = rcar/rnit       ! iHAMOCC C:N ratio
      ro2utammo     = 140.            ! Oxygen utilization per mol detritus during ammonification
      ro2nnit       = ro2utammo/rnit  !  
      rnoxp         = 280.            ! consumption of NOx per mol detritus during denitrification
      rnoxpi        = 1./rnoxp        ! inverse
      rno2anmx      = 1144.           ! consumption of NO2 per mol organic production by anammox
      rno2anmxi     = 1./rno2anmx     ! inverse
      rnh4anmx      = 880.            ! consumption of NH4 per mol organic production by anammox
      rnh4anmxi     = 1./rnh4anmx     ! inverse
      rno2dnra      = 93. + 1./3.     ! consumption of NO2 per mol OM degradation during DNRA
      rno2dnrai     = 1./rno2dnra     ! inverse
      rnh4dnra      = rno2dnra + rnit ! production of NH4 per mol OM during DNRA
      rnh4dnrai     = 1./rnh4dnra     ! inverse
      rnm1          = rnit - 1.       

      ! Phytoplankton growth     
      bkphyanh4     = 0.12e-6    ! Half-saturation constant for NH4 uptake by bulk phytoplankton (kmol/m3)
      bkphyano3     = 0.16e-6   ! Half-saturation constant for NO3 uptake by bulk phytoplankton (kmol/m3)
      bkphosph      = 0.01e-6   ! Half-saturation constant for PO4 uptake by bulk phytoplankton (kmol/m3)
      bkiron        = bkphosph*riron ! Half-saturation constant for Fe uptake by bulk phytoplankton (kmol/m3)

      ! === Denitrification step NO3 -> NO2:
      !rano3denit    = 0.15*dtb ! Maximum growth rate denitrification on NO3 at reference T (1/d -> 1/dt)
      rano3denit    = 0.05     ! Maximum growth rate denitrification on NO3 at reference T (1/d -> 1/dt)
      q10ano3denit  = 2.       ! Q10 factor for denitrification on NO3 (-)
      Trefano3denit = 10.      ! Reference temperature for denitrification on NO3 (degr C) 
      !sc_ano3denit  = 0.05e6   ! Shape factor for NO3 denitrification oxygen inhibition function (m3/kmol)
      sc_ano3denit  = 0.12e6   ! Shape factor for NO3 denitrification oxygen inhibition function (m3/kmol)
      bkano3denit   = 5.e-6    ! Half-saturation constant for NO3 denitrification (kmol/m3)

      ! === Anammox
      rano2anmx     = 0.05     ! Maximum growth rate for anammox at reference T (1/d -> 1/dt)
      q10anmx       = 1.6      ! Q10 factor for anammox (-)
      Trefanmx      = 10.      ! Reference temperature for anammox (degr C)
      alphaanmx     = 0.45e6   ! Shape factor for anammox oxygen inhibition function (m3/kmol)
      bkoxanmx      = 11.3e-6  ! Half-saturation constant for oxygen inhibition function (kmol/m3)
      bkano2anmx    = 5.e-6    ! Half-saturation constant for NO2 limitation (kmol/m3)
      bkanh4anmx    = bkano2anmx * rnh4anmx/rno2anmx !Half-saturation constant for NH4 limitation of anammox (kmol/m3)

      ! === Denitrification step NO2 -> N2O
      rano2denit    = 0.12     ! Maximum growth rate denitrification on NO2 at reference T (1/d -> 1/dt) 
      q10ano2denit  = 2.0      ! Q10 factor for denitrification on NO2 (-)
      Trefano2denit = 10.      ! Reference temperature for denitrification on NO2 (degr C)
      bkoxano2denit = 2.e-6    ! Half-saturation constant for (quadratic) oxygen inhibition function of denitrification on NO2 (kmol/m3)
      bkano2denit   = 5.6e-6   ! Half-saturation constant for denitrification on NO2 (kmol/m3)

      ! === Denitrification step N2O -> N2
      ran2odenit    = 0.16     ! Maximum growth rate denitrification on N2O at reference T (1/d -> 1/dt)
      q10an2odenit  = 3.       ! Q1- factor for denitrificationj on N2O (-)
      Trefan2odenit = 10.      ! Reference temperature for denitrification on N2O (degr C)
      bkoxan2odenit = 5.e-6    ! Half-saturation constant for (quadratic) oxygen inhibition function of denitrification on N2O (kmol/m3)
      bkan2odenit   = 1.e-6    ! Half-saturation constant for denitrification on N2O (kmol/m3)

      ! === DNRA NO2 -> NH4
      rdnra         = 0.1      ! Maximum growth rate DNRA on NO2 at reference T (1/d -> 1/dt)
      q10dnra       = 2.       ! Q10 factor for DNRA on NO2 (-)
      Trefdnra      = 10.      ! Reference temperature for DNRA (degr C)
      bkoxdnra      = 2.5e-6   ! Half saturation constant for (quadratic) oxygen inhibition function of DNRA on NO2 (kmol/m3)
      bkdnra        = 0.05e-6  ! Half-saturation constant for DNRA on NO2 (kmol/m3)

      ! === Nitrification on NH4
      ranh4nitr     = 1.       ! Maximum growth rate nitrification on NH4 at reference T (1/d -> 1/dt) 
      q10anh4nitr   = 3.3      ! Q10 factor for nitrification on NH4 (-)
      Trefanh4nitr  = 20.      ! Reference temperature for nitrification on NH4 (degr C)
      bkoxamox      = 0.333e-6 ! Half-saturation constant for oxygen limitation of nitrification on NH4 (kmol/m3)
      bkanh4nitr    = 0.133e-6 ! Half-saturation constant for nitrification on NH4 (kmol/m3)
!======
! OLD VERSION OF pathway splitting function
      !bkamoxn2o     = 0.453e-6 ! Half saturation constant for O2 in pathway splitting function N2O for nitrification on NH4 (kmol/m3)
! NEW version similar to Santoros 2021, Ji 2018:      
      bkamoxn2o     = 0.5e-6 ! Half saturation constant for NH4 in pathway splitting function N2O for nitrification on NH4 (kmol/m3)
      mufn2o        = 0.11/(50.*1e6*bkoxamox) !=6.61e-3  0.11/(50*1e6)=2.2e-9 - ~Santoro et al. 2011 with simple MM,  
      bn2o          = 0.077/(50.*mufn2o)  !=0.2331 - before set to 0.3 - base fraction entering N2O 
!======
      !bkamoxno2     = 0.479e-6 ! Half saturation constant for pathway splitting function N2O for nitrification on NH4 (kmol/m3)
!      bkamoxno2     = 0.1e-6 ! Half saturation constant for pathway splitting function N2O for nitrification on NH4 (kmol/m3)
      n2omaxy       = 0.003    ! Maximum yield of OM on NH4 nitrification (-)
      n2oybeta      = 18.      ! Decay factor for inhibition function for yield during nitrification on NH4 (kmol/m3)
      bkyamox       = 0.333e-6 ! Half saturation constant for pathway splitting function OM-yield for nitrification on NH4 (kmol/m3)
  
      ! === Nitrification on NO2
      rano2nitr     = 1.54     ! Maximum growth rate nitrification on NO2 at reference T (1/d -> 1/dt) 
      q10ano2nitr   = 2.7      ! Q10 factor for nitrification on NO2 (-)
      Trefano2nitr  = 20.      ! Reference temperature for nitrification on NO2 (degr C)
      bkoxnitr      = 0.788e-6 ! Half-saturation constant for oxygen limitation of nitrification on NO2 (kmol/m3)
      bkano2nitr    = 0.287e-6 ! Half-saturation constant for NO2 for nitrification on NO2 (kmol/m3)
      NOB2AOAy      = 0.44     ! Ratio of NOB versus AOA yield per energy source ~0.043/0.098 according to Zakem et al. 2022

      eps    = 1.e-25 ! safe division etc. 
      minlim = 1.e-9  ! minimum for limitation functions (e.g. nutlim or oxlim/inh can only decrease to minlim) 
      !===========================================================================

      ! Tweaked parameters:
      rano3denit    = 0.0005 ! Maximum growth rate denitrification on NO3 at reference T (1/d -> 1/dt)
      rano2anmx     = 0.001  ! Maximum growth rate for anammox at reference T (1/d -> 1/dt)
      rano2denit    = 0.001  ! Maximum growth rate denitrification on NO2 at reference T (1/d -> 1/dt) 
      ran2odenit    = 0.0012 ! Maximum growth rate denitrification on N2O at reference T (1/d -> 1/dt)
      rdnra         = 0.001  ! Maximum growth rate DNRA on NO2 at reference T (1/d -> 1/dt)

      end subroutine extNwatercol_param_init
 
!==================================================================================================================================      
      subroutine extNwatercol_param_update()

        rano3denit = rano3denit *dtb ! Maximum growth rate denitrification on NO3 at reference T (1/d -> 1/dt)
        rano2anmx  = rano2anmx  *dtb ! Maximum growth rate for anammox at reference T (1/d -> 1/dt)
        rano2denit = rano2denit *dtb ! Maximum growth rate denitrification on NO2 at reference T (1/d -> 1/dt) 
        ran2odenit = ran2odenit *dtb ! Maximum growth rate denitrification on N2O at reference T (1/d -> 1/dt)
        rdnra      = rdnra      *dtb ! Maximum growth rate DNRA on NO2 at reference T (1/d -> 1/dt)
        ranh4nitr  = ranh4nitr  *dtb ! Maximum growth rate nitrification on NH4 at reference T (1/d -> 1/dt) 
        rano2nitr  = rano2nitr  *dtb ! Maximum growth rate nitrification on NO2 at reference T (1/d -> 1/dt) 
      
      end subroutine extNwatercol_param_update        

!==================================================================================================================================      
      subroutine extNwatercol_param_write()

          REAL :: dtbinv
          dtbinv = 1./dtb
          WRITE(io_stdo_bgc,*) '****************************************************************'
          WRITE(io_stdo_bgc,*) '* HAMOCC extended nitrogen cycle
parameters water column:' 
          WRITE(io_stdo_bgc,*) '*          rc2n          = ',rc2n
          WRITE(io_stdo_bgc,*) '*          ro2utammo     = ',ro2utammo          
          WRITE(io_stdo_bgc,*) '*          ro2nnit       = ',ro2nnit
          WRITE(io_stdo_bgc,*) '*          rnoxp         = ',rnoxp
          WRITE(io_stdo_bgc,*) '*          rnoxpi        = ',rnoxpi
          WRITE(io_stdo_bgc,*) '*          rno2anmx      = ',rno2anmx
          WRITE(io_stdo_bgc,*) '*          rno2anmxi     = ',rno2anmxi
          WRITE(io_stdo_bgc,*) '*          rnh4anmx      = ',rnh4anmx
          WRITE(io_stdo_bgc,*) '*          rnh4anmxi     = ',rnh4anmxi
          WRITE(io_stdo_bgc,*) '*          rno2dnra      = ',rno2dnra
          WRITE(io_stdo_bgc,*) '*          rno2dnrai     = ',rno2dnrai
          WRITE(io_stdo_bgc,*) '*          rnh4dnra      = ',rnh4dnra
          WRITE(io_stdo_bgc,*) '*          rnh4dnrai     = ',rnh4dnrai
          WRITE(io_stdo_bgc,*) '*          rnm1          = ',rnm1
          WRITE(io_stdo_bgc,*) '*          bkphyanh4     = ',bkphyanh4
          WRITE(io_stdo_bgc,*) '*          bkphyano3     = ',bkphyano3
          WRITE(io_stdo_bgc,*) '*          bkphosph      = ',bkphosph
          WRITE(io_stdo_bgc,*) '*          bkiron        = ',bkiron
          WRITE(io_stdo_bgc,*) '*          rano3denit    = ',rano3denit    *dtbinv
          WRITE(io_stdo_bgc,*) '*          q10ano3denit  = ',q10ano3denit
          WRITE(io_stdo_bgc,*) '*          Trefano3denit = ',Trefano3denit
          WRITE(io_stdo_bgc,*) '*          sc_ano3denit  = ',sc_ano3denit
          WRITE(io_stdo_bgc,*) '*          bkano3denit   = ',bkano3denit
          WRITE(io_stdo_bgc,*) '*          rano2anmx     = ',rano2anmx     *dtbinv
          WRITE(io_stdo_bgc,*) '*          q10anmx       = ',q10anmx
          WRITE(io_stdo_bgc,*) '*          Trefanmx      = ',Trefanmx
          WRITE(io_stdo_bgc,*) '*          alphaanmx     = ',alphaanmx
          WRITE(io_stdo_bgc,*) '*          bkoxanmx      = ',bkoxanmx
          WRITE(io_stdo_bgc,*) '*          bkano2anmx    = ',bkano2anmx
          WRITE(io_stdo_bgc,*) '*          bkanh4anmx    = ',bkanh4anmx
          WRITE(io_stdo_bgc,*) '*          rano2denit    = ',rano2denit    *dtbinv
          WRITE(io_stdo_bgc,*) '*          q10ano2denit  = ',q10ano2denit
          WRITE(io_stdo_bgc,*) '*          Trefano2denit = ',Trefano2denit
          WRITE(io_stdo_bgc,*) '*          bkoxano2denit = ',bkoxano2denit
          WRITE(io_stdo_bgc,*) '*          bkano2denit   = ',bkano2denit
          WRITE(io_stdo_bgc,*) '*          ran2odenit    = ',ran2odenit    *dtbinv
          WRITE(io_stdo_bgc,*) '*          q10an2odenit  = ',q10an2odenit
          WRITE(io_stdo_bgc,*) '*          Trefan2odenit = ',Trefan2odenit
          WRITE(io_stdo_bgc,*) '*          bkoxan2odenit = ',bkoxan2odenit
          WRITE(io_stdo_bgc,*) '*          bkan2odenit   = ',bkan2odenit
          WRITE(io_stdo_bgc,*) '*          rdnra         = ',rdnra         *dtbinv
          WRITE(io_stdo_bgc,*) '*          q10dnra       = ',q10dnra
          WRITE(io_stdo_bgc,*) '*          Trefdnra      = ',Trefdnra
          WRITE(io_stdo_bgc,*) '*          bkoxdnra      = ',bkoxdnra
          WRITE(io_stdo_bgc,*) '*          bkdnra        = ',bkdnra
          WRITE(io_stdo_bgc,*) '*          ranh4nitr     = ',ranh4nitr     *dtbinv
          WRITE(io_stdo_bgc,*) '*          q10anh4nitr   = ',q10anh4nitr
          WRITE(io_stdo_bgc,*) '*          Trefanh4nitr  = ',Trefanh4nitr
          WRITE(io_stdo_bgc,*) '*          bkoxamox      = ',bkoxamox
          WRITE(io_stdo_bgc,*) '*          bkanh4nitr    = ',bkanh4nitr
          WRITE(io_stdo_bgc,*) '*          bkamoxn2o     = ',bkamoxn2o
          WRITE(io_stdo_bgc,*) '*          mufn2o        = ',mufn2o
          WRITE(io_stdo_bgc,*) '*          bn2o          = ',bn2o
          WRITE(io_stdo_bgc,*) '*          n2omaxy       = ',n2omaxy
          WRITE(io_stdo_bgc,*) '*          n2oybeta      = ',n2oybeta
          WRITE(io_stdo_bgc,*) '*          bkyamox       = ',bkyamox
          WRITE(io_stdo_bgc,*) '*          rano2nitr     = ',rano2nitr     *dtbinv
          WRITE(io_stdo_bgc,*) '*          q10ano2nitr   = ',q10ano2nitr
          WRITE(io_stdo_bgc,*) '*          Trefano2nitr  = ',Trefano2nitr
          WRITE(io_stdo_bgc,*) '*          bkoxnitr      = ',bkoxnitr
          WRITE(io_stdo_bgc,*) '*          bkano2nitr    = ',bkano2nitr
          WRITE(io_stdo_bgc,*) '*          NOB2AOAy      = ',NOB2AOAy

      end subroutine extNwatercol_param_write

!==================================================================================================================================      
      subroutine nitrification(kpie,kpje,kpke,kbnd,pddpo,omask,ptho)
      ! Nitrification processes (NH4 -> NO2, NO2 -> NO3) accompanied
      ! by dark carbon fixation and O2-dependent N2O production 

      integer, intent(in) :: kpie,kpje,kpke,kbnd
      real,    intent(in) :: omask(kpie,kpje)       
      real,    intent(in) :: pddpo(kpie,kpje,kpke)
      real,    intent(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)

      !local variables
      integer :: i,j,k
      real    :: Tdepanh4,O2limanh4,nut1lim,anh4new,potdnh4amox,fdetamox,fno2,fn2o,ftotnh4 
      real    :: Tdepano2,O2limano2,nut2lim,ano2new,potdno2nitr,fdetnitr,ftotno2,no2fn2o,no2fno2,no2fdetamox
      real    :: amoxfrac,nitrfrac,totd,amox,nitr,temp

      real    :: minlim_oxnh4,minlim_nh4,minlim_oxno2,minlim_no2 ! minimum conc for limitation functions 
 
      minlim_oxnh4  = bkoxamox*minlim/(1. - minlim) 
      minlim_oxno2  = bkoxnitr*minlim/(1. - minlim)  
      minlim_nh4    = bkanh4nitr*minlim/(1. - minlim) 
      minlim_no2    = bkano2nitr*minlim/(1. - minlim)

      ! Set output-related fields to zero
      nitr_NH4      = 0.
      nitr_NO2      = 0.
      nitr_N2O_prod = 0.
      nitr_NH4_OM   = 0.
      nitr_NO2_OM   = 0.

      !$OMP PARALLEL DO PRIVATE(i,k,Tdepanh4,O2limanh4,nut1lim,anh4new,potdnh4amox,fdetamox,fno2,fn2o,ftotnh4,   & 
      !$OMP                     Tdepano2,O2limano2,nut2lim,ano2new,potdno2nitr,fdetnitr,ftotno2,amoxfrac,        &
      !$OMP                     nitrfrac,totd,amox,nitr,temp,no2fn2o,no2fno2,no2fdetamox)

      do j = 1,kpje
      do i = 1,kpie
      do k = 1,kpke
        if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then
          potdnh4amox = 0.
          fn2o        = 0.
          fno2        = 0.
          fdetamox    = 0.
          potdno2nitr = 0.
          fdetnitr    = 0.

!          if(ocetra(i,j,k,ioxygen)>minlim_oxnh4 .and. ocetra(i,j,k,ianh4)>minlim_nh4)then
           temp = merge(ptho(i,j,k),10.,ptho(i,j,k)<40.)
           ! Ammonium oxidation step of nitrification
           Tdepanh4    = q10anh4nitr**((temp-Trefanh4nitr)/10.) 
           O2limanh4   = ocetra(i,j,k,ioxygen)/(ocetra(i,j,k,ioxygen) + bkoxamox)
           nut1lim     = ocetra(i,j,k,ianh4)/(ocetra(i,j,k,ianh4) + bkanh4nitr)
           anh4new     = ocetra(i,j,k,ianh4)/(1. + ranh4nitr*Tdepanh4*O2limanh4*nut1lim)
           potdnh4amox = max(0.,ocetra(i,j,k,ianh4) - anh4new)
            
            ! pathway splitting functions according to Goreau 1980
       !=====
       ! OLD version according to Goreau
            !fn2o     = 1. - ocetra(i,j,k,ioxygen)/(ocetra(i,j,k,ioxygen) + bkamoxn2o)
       ! NEW version similar to Santoros et al. 2021, Ji et al. 2018
           fn2o     = mufn2o * (bn2o + (1.-bn2o)*bkoxamox/(ocetra(i,j,k,ioxygen)+bkoxamox))                                        &
                    &        * ocetra(i,j,k,ianh4)/(ocetra(i,j,k,ianh4)+bkamoxn2o)
       !=====
           fno2     = ocetra(i,j,k,ioxygen)/(ocetra(i,j,k,ioxygen) + bkoxamox)
           fdetamox = n2omaxy*2.*(1. + n2oybeta)*ocetra(i,j,k,ioxygen)*bkyamox &
                     & /(ocetra(i,j,k,ioxygen)**2 + 2.*ocetra(i,j,k,ioxygen)*bkyamox + bkyamox**2)

           ! normalization of pathway splitting functions to sum=1
           ftotnh4  = fn2o + fno2 + fdetamox + eps
           fn2o     = fn2o/ftotnh4
           fno2     = fno2/ftotnh4
           fdetamox = 1. - (fn2o + fno2)
!          endif

!          if(ocetra(i,j,k,ioxygen)>minlim_oxno2 .and. ocetra(i,j,k,iano2)>minlim_no2)then
           temp = merge(ptho(i,j,k),10.,ptho(i,j,k)<40.)
           ! NO2 oxidizing step of nitrification
           Tdepano2    = q10ano2nitr**((temp-Trefano2nitr)/10.) 
           O2limano2   = ocetra(i,j,k,ioxygen)/(ocetra(i,j,k,ioxygen) + bkoxnitr)
           nut2lim     = ocetra(i,j,k,iano2)/(ocetra(i,j,k,iano2) + bkano2nitr)
           ano2new     = ocetra(i,j,k,iano2)/(1. + rano2nitr*Tdepano2*O2limano2*nut2lim)
           potdno2nitr = max(0.,ocetra(i,j,k,iano2) - ano2new)

           ! pathway splitting functions for NO2 nitrification - assuming to be the same as for NH4
           ! but with reduced OM gain per used NO2 as energy source (in amox: NH4)
        !=====
        ! OLD version according to Goreau
           ! no2fn2o     = 1. - ocetra(i,j,k,ioxygen)/(ocetra(i,j,k,ioxygen) + bkamoxn2o)
        ! NEW version
           no2fn2o     = mufn2o * (bn2o + (1.-bn2o)*bkoxamox/(ocetra(i,j,k,ioxygen)+bkoxamox))                                     &
                       &        * ocetra(i,j,k,ianh4)/(ocetra(i,j,k,ianh4)+bkamoxn2o)
        !=====
           no2fno2     = ocetra(i,j,k,ioxygen)/(ocetra(i,j,k,ioxygen) + bkoxamox)
           no2fdetamox = NOB2AOAy*n2omaxy*2.*(1. + n2oybeta)*ocetra(i,j,k,ioxygen)*bkyamox &
                     & /(ocetra(i,j,k,ioxygen)**2 + 2.*ocetra(i,j,k,ioxygen)*bkyamox + bkyamox**2)

           fdetnitr = no2fdetamox/(no2fno2 + no2fn2o)   ! yield to energy usage ratio for NO2 -> ratio equals 16:x
!          endif

          ! limitation of the two processes through available nutrients, etc.
          totd     = potdnh4amox + potdno2nitr
          amoxfrac = potdnh4amox/(totd + eps)
          nitrfrac = 1. - amoxfrac
           
          totd     = max(0.,                                                                                                      &
                   &   min(totd,                                                                                                  &
                   &       ocetra(i,j,k,ianh4)/(amoxfrac + fdetnitr*nitrfrac + eps),                                              & ! ammonium
                   &       ocetra(i,j,k,isco212)/(rc2n*(fdetamox*amoxfrac + fdetnitr*nitrfrac) + eps),                            & ! CO2
                   &       ocetra(i,j,k,iphosph)/(rnoi*(fdetamox*amoxfrac + fdetnitr*nitrfrac) + eps),                            & ! PO4
                   &       ocetra(i,j,k,iiron)/(riron*rnoi*(fdetamox*amoxfrac + fdetnitr*nitrfrac) + eps),                        & ! Fe
                   &       ocetra(i,j,k,ioxygen)                                                                                  &
                   &       /((1.5*fno2 + fn2o - ro2nnit*fdetamox)*amoxfrac + (0.5 - ro2nnit*fdetnitr)*nitrfrac + eps),            & ! O2
                   &       ocetra(i,j,k,ialkali)                                                                                  &
                   &       /((2.*fno2 + fn2o + rnm1*rnoi*fdetamox)*amoxfrac + (rnm1*rnoi*fdetnitr)*nitrfrac + eps)))                ! alkalinity
          amox     = amoxfrac*totd 
          nitr     = nitrfrac*totd

          ocetra(i,j,k,ianh4)   = ocetra(i,j,k,ianh4) - amox - fdetnitr*nitr
          ocetra(i,j,k,ian2o)   = ocetra(i,j,k,ian2o) + 0.5*fn2o*amox
          ocetra(i,j,k,iano2)   = ocetra(i,j,k,iano2) + fno2*amox - nitr
          ocetra(i,j,k,iano3)   = ocetra(i,j,k,iano3) + nitr
          ocetra(i,j,k,idet)    = ocetra(i,j,k,idet)  + rnoi*(fdetamox*amox + fdetnitr*nitr)
          ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212) - rc2n*(fdetamox*amox + fdetnitr*nitr)
          ocetra(i,j,k,iphosph) = ocetra(i,j,k,iphosph) - rnoi*(fdetamox*amox + fdetnitr*nitr)
          ocetra(i,j,k,iiron)   = ocetra(i,j,k,iiron)   - riron*rnoi*(fdetamox*amox + fdetnitr*nitr)
          ocetra(i,j,k,ioxygen) = ocetra(i,j,k,ioxygen) - (1.5*fno2 + fn2o - ro2nnit*fdetamox)*amox   &
                                &                       - (0.5 - ro2nnit*fdetnitr)*nitr
          ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali) - (2.*fno2 + fn2o + rnm1*rnoi*fdetamox)*amox - rnm1*rnoi*fdetnitr*nitr

          ! Output
          nitr_NH4(i,j,k)       = amox               ! kmol N/m3/dtb   - NH4 consumption for nitrification on NH4-incl. usage for biomass
          nitr_NO2(i,j,k)       = nitr               ! kmol N/m3/dtb   - NO2 consumption for nitrification on NO2
          nitr_N2O_prod(i,j,k)  = 0.5*fn2o*amox      ! kmol N2O/m3/dtb - N2O production during aerob ammonium oxidation
          nitr_NH4_OM(i,j,k)    = rnoi*fdetamox*amox ! kmol P/m3/dtb   - organic matter production during aerob NH4 oxidation
          nitr_NO2_OM(i,j,k)    = rnoi*fdetnitr*nitr ! kmol P/m3/dtb   - organic matter production during aerob NO2 oxidation
       
        endif
      enddo
      enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine nitrification

!==================================================================================================================================      
      subroutine denit_NO3_to_NO2(kpie,kpje,kpke,kbnd,pddpo,omask,ptho)
      ! Denitrification / dissimilatory nitrate reduction (NO3 -> NO2)

      integer, intent(in) :: kpie,kpje,kpke,kbnd
      real,    intent(in) :: omask(kpie,kpje)       
      real,    intent(in) :: pddpo(kpie,kpje,kpke)
      real,    intent(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)

      !local variables
      integer :: i,j,k
      real    :: Tdep,O2inhib,nutlim,ano3new,ano3denit,temp

      real    :: minlim_ox,minlim_no3 ! minimum conc for limitation functions 

      minlim_ox  = log(2./minlim-1.)/(2.*sc_ano3denit) 
      minlim_no3 = bkano3denit*minlim/(1.-minlim)

      ! Sett output-related field to zero
      denit_NO3  = 0.

      !$OMP PARALLEL DO PRIVATE(i,k,Tdep,O2inhib,nutlim,ano3new,ano3denit,temp)
      do j = 1,kpje
      do i = 1,kpie
      do k = 1,kpke
        if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then
!          if(ocetra(i,j,k,ioxygen) < minlim_ox .and. ocetra(i,j,k,iano3)>minlim_no3)then
            temp      = merge(ptho(i,j,k),10.,ptho(i,j,k)<40.)
            Tdep      = q10ano3denit**((temp-Trefano3denit)/10.) 
            O2inhib   = 1. - tanh(sc_ano3denit*ocetra(i,j,k,ioxygen)) 
            nutlim    = ocetra(i,j,k,iano3)/(ocetra(i,j,k,iano3) + bkano3denit)
 
            ano3new   = ocetra(i,j,k,iano3)/(1. + rano3denit*Tdep*O2inhib*nutlim) 

            ano3denit = max(0.,min(ocetra(i,j,k,iano3) - ano3new, ocetra(i,j,k,idet)*rnoxp))

            ocetra(i,j,k,iano3)   = ocetra(i,j,k,iano3)   - ano3denit
            ocetra(i,j,k,iano2)   = ocetra(i,j,k,iano2)   + ano3denit
            ocetra(i,j,k,idet)    = ocetra(i,j,k,idet)    - ano3denit*rnoxpi
            ocetra(i,j,k,ianh4)   = ocetra(i,j,k,ianh4)   + ano3denit*rnit*rnoxpi
            ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212) + ano3denit*rcar*rnoxpi
            ocetra(i,j,k,iphosph) = ocetra(i,j,k,iphosph) + ano3denit*rnoxpi
            ocetra(i,j,k,iiron)   = ocetra(i,j,k,iiron)   + ano3denit*riron*rnoxpi
            ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali) + ano3denit*rnm1*rnoxpi

            ! Output
            denit_NO3(i,j,k) = ano3denit ! kmol NO3/m3/dtb   - NO3 usage for denit on NO3    
!          endif
        endif
      enddo
      enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine denit_NO3_to_NO2

!==================================================================================================================================      
      subroutine anammox(kpie,kpje,kpke,kbnd,pddpo,omask,ptho)
      ! Aanammox

      integer, intent(in) :: kpie,kpje,kpke,kbnd
      real,    intent(in) :: omask(kpie,kpje)       
      real,    intent(in) :: pddpo(kpie,kpje,kpke)
      real,    intent(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)

      !local variables
      integer :: i,j,k
      real    :: Tdep,O2inhib,nut1lim,nut2lim,ano2new,ano2anmx,temp

      real    :: minlim_ox,minlim_nh4,minlim_no2 ! minimum conc for limitation functions

      minlim_ox  = log((1.-minlim)/minlim)/alphaanmx + bkoxanmx
      minlim_nh4 = bkanh4anmx*minlim/(1.-minlim)
      minlim_no2 = bkano2anmx*minlim/(1.-minlim) 

      ! Set output-related field to zero
      anmx_N2_prod = 0.
      anmx_OM_prod = 0.

      !$OMP PARALLEL DO PRIVATE(i,k,Tdep,O2inhib,nut1lim,nut2lim,ano2new,ano2anmx,temp)
      do j = 1,kpje
      do i = 1,kpie
      do k = 1,kpke
        if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then
!          if(ocetra(i,j,k,iano2)>minlim_no2 .and. ocetra(i,j,k,ianh4)>minlim_nh4 .and. ocetra(i,j,k,ioxygen)<minlim_ox) then
           temp     = merge(ptho(i,j,k),10.,ptho(i,j,k)<40.)
           Tdep     = q10anmx**((temp-Trefanmx)/10.)         
           O2inhib  = 1. - exp(alphaanmx*(ocetra(i,j,k,ioxygen)-bkoxanmx))/(1.+ exp(alphaanmx*(ocetra(i,j,k,ioxygen)-bkoxanmx))) 
           nut1lim  = ocetra(i,j,k,iano2)/(ocetra(i,j,k,iano2)+bkano2anmx)
           nut2lim  = ocetra(i,j,k,ianh4)/(ocetra(i,j,k,ianh4)+bkanh4anmx)

           ano2new  = ocetra(i,j,k,iano2)/(1. + rano2anmx*Tdep*O2inhib*nut1lim*nut2lim) 

           ano2anmx = max(0.,min(ocetra(i,j,k,iano2) - ano2new, ocetra(i,j,k,ianh4)*rno2anmx*rnh4anmxi,                            &
                             ocetra(i,j,k,isco212)*rno2anmx/rcar, ocetra(i,j,k,iphosph)*rno2anmx,                                  &
                             ocetra(i,j,k,iiron)*rno2anmx/riron, ocetra(i,j,k,ialkali)*rno2anmx/rnm1))

           ocetra(i,j,k,iano2)   = ocetra(i,j,k,iano2)   - ano2anmx
           ocetra(i,j,k,ianh4)   = ocetra(i,j,k,ianh4)   - ano2anmx*rnh4anmx*rno2anmxi
           ocetra(i,j,k,igasnit) = ocetra(i,j,k,igasnit) + ano2anmx*(rnh4anmx-rnit)*rno2anmxi
           ocetra(i,j,k,iano3)   = ocetra(i,j,k,iano3)   + ano2anmx*rnoxp*rno2anmxi
           ocetra(i,j,k,idet)    = ocetra(i,j,k,idet)    + ano2anmx*rno2anmxi
           ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212) - ano2anmx*rcar*rno2anmxi
           ocetra(i,j,k,iphosph) = ocetra(i,j,k,iphosph) - ano2anmx*rno2anmxi
           ocetra(i,j,k,iiron)   = ocetra(i,j,k,iiron)   - ano2anmx*riron*rno2anmxi
           ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali) - ano2anmx*rnm1*rno2anmxi

           ! Output
           anmx_N2_prod(i,j,k) = ano2anmx*(rnh4anmx-rnit)*rno2anmxi  ! kmol N2/m3/dtb - N2 prod through anammox
           anmx_OM_prod(i,j,k) = ano2anmx*rno2anmxi                  ! kmol P/m3/dtb  - OM production by anammox
!          endif
        endif
      enddo
      enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine anammox

!==================================================================================================================================      
      subroutine denit_dnra(kpie,kpje,kpke,kbnd,pddpo,omask,ptho)
      ! Denitrification processes (NO2 -> N2O -> N2) and dissmilatory nitrite reduction (NO2 -> NH4)

      integer, intent(in) :: kpie,kpje,kpke,kbnd
      real,    intent(in) :: omask(kpie,kpje)       
      real,    intent(in) :: pddpo(kpie,kpje,kpke)
      real,    intent(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)

      !local variables
      integer :: i,j,k
      real    :: Tdepano2,O2inhibano2,nutlimano2,detlimano2,rpotano2denit,ano2denit
      real    :: Tdepdnra,O2inhibdnra,nutlimdnra,detlimdnra,rpotano2dnra,ano2dnra
      real    :: fdenit,fdnra,potano2new,potdano2,potddet,fdetano2denit,fdetan2odenit,fdetdnra  
      real    :: Tdepan2o,O2inhiban2o,nutliman2o,detliman2o,an2onew,an2odenit

      real    :: temp


      real    :: minlim_ox,minlim_oxn2o,minlim_no2,minlim_n2o

      minlim_ox     = min(bkoxano2denit,bkoxdnra)/sqrt(minlim)
      minlim_oxn2o  = bkoxan2odenit/sqrt(minlim)
      minlim_no2    = min(bkdnra,bkano2denit)*minlim/(1. - minlim)
      minlim_n2o    = bkan2odenit*minlim/(1. - minlim)
      
      ! Set output-related field to zero
      denit_NO2 = 0.
      denit_N2O = 0.
      DNRA_NO2  = 0.

      !$OMP PARALLEL DO PRIVATE(i,k,Tdepano2,O2inhibano2,nutlimano2,detlimano2,ano2denit,      &
      !$OMP                     Tdepan2o,O2inhiban2o,nutliman2o,detliman2o,an2onew,an2odenit,  &
      !$OMP                     rpotano2denit,rpotano2dnra,                                    &
      !$OMP                     fdenit,fdnra,potano2new,potdano2,potddet,fdetano2denit,        &
      !$OMP                     fdetan2odenit,fdetdnra,                                        &
      !$OMP                     Tdepdnra,O2inhibdnra,nutlimdnra,detlimdnra,ano2dnra,temp)

      do j = 1,kpje
      do i = 1,kpie
      do k = 1,kpke
        if(pddpo(i,j,k) > dp_min .and. omask(i,j) > 0.5) then
          potddet       = 0.
          an2odenit     = 0.
          ano2denit     = 0.
          ano2dnra      = 0.

!          if(0.<=ocetra(i,j,k,ioxygen) .and. ocetra(i,j,k,ioxygen)<minlim_oxn2o .and. ocetra(i,j,k,ian2o)>minlim_n2o)then
           temp = merge(ptho(i,j,k),10.,ptho(i,j,k)<40.)
           ! === denitrification on N2O
           Tdepan2o    = q10an2odenit**((temp-Trefan2odenit)/10.) 
           O2inhiban2o = bkoxan2odenit**2/(ocetra(i,j,k,ioxygen)**2 + bkoxan2odenit**2) 
           nutliman2o  = ocetra(i,j,k,ian2o)/(ocetra(i,j,k,ian2o) + bkan2odenit)   
           an2onew     = ocetra(i,j,k,ian2o)/(1. + ran2odenit*Tdepan2o*O2inhiban2o*nutliman2o)  
           an2odenit   = max(0.,min(ocetra(i,j,k,ian2o),ocetra(i,j,k,ian2o) - an2onew))
!          endif

!          if(0.<=ocetra(i,j,k,ioxygen) .and. ocetra(i,j,k,ioxygen)<minlim_ox .and. ocetra(i,j,k,iano2)>minlim_no2)then 
           temp = merge(ptho(i,j,k),10.,ptho(i,j,k)<40.)
           ! denitrification on NO2
           Tdepano2    =  q10ano2denit**((temp-Trefano2denit)/10.) 
           O2inhibano2 = bkoxano2denit**2/(ocetra(i,j,k,ioxygen)**2 + bkoxano2denit**2) 
           nutlimano2  = ocetra(i,j,k,iano2)/(ocetra(i,j,k,iano2) + bkano2denit)
           rpotano2denit = max(0.,rano2denit*Tdepano2*O2inhibano2*nutlimano2) ! potential rate of denit
           
           ! DNRA on NO2
           Tdepdnra    = q10dnra**((temp-Trefdnra)/10.) 
           O2inhibdnra = bkoxdnra**2/(ocetra(i,j,k,ioxygen)**2 + bkoxdnra**2) 
           nutlimdnra  = ocetra(i,j,k,iano2)/(ocetra(i,j,k,iano2) + bkdnra)
           rpotano2dnra = max(0.,rdnra*Tdepdnra*O2inhibdnra*nutlimdnra) ! pot. rate of dnra

           ! potential new conc of NO2 due to denitrification and DNRA
           potano2new = ocetra(i,j,k,iano2)/(1. + rpotano2denit + rpotano2dnra)
           potdano2   = max(0.,min(ocetra(i,j,k,iano2), ocetra(i,j,k,iano2) - potano2new))
           
           ! === limitation due to NO2:
           ! fraction on potential change of NO2:
           fdenit = rpotano2denit/(rpotano2denit + rpotano2dnra + eps)
           fdnra  = 1. - fdenit
           
           ! potential fractional change
           ano2denit  = fdenit * potdano2   
           ano2dnra   = fdnra  * potdano2
!          endif

          ! limitation of processes due to detritus
          potddet       = rnoxpi*(ano2denit + an2odenit) + rno2dnrai*ano2dnra  ! P units              
          fdetano2denit = rnoxpi*ano2denit/(potddet + eps)
          fdetan2odenit = rnoxpi*an2odenit/(potddet + eps)
          fdetdnra      = 1. - fdetano2denit - fdetan2odenit 
          potddet       = max(0.,min(potddet,ocetra(i,j,k,idet))) 
       
!          if(potddet>0.)then
           ! change of NO2 and N2O in N units
           ano2denit     = fdetano2denit*rnoxp*potddet
           an2odenit     = fdetan2odenit*rnoxp*potddet
           ano2dnra      = fdetdnra*rno2dnra*potddet

           ! change in tracer concentrations due to denit (NO2->N2O->N2) and DNRA (NO2->NH4)
           ocetra(i,j,k,iano2)   = ocetra(i,j,k,iano2)   - ano2denit - ano2dnra
           ocetra(i,j,k,ian2o)   = ocetra(i,j,k,ian2o)   - an2odenit + 0.5*ano2denit
           ocetra(i,j,k,igasnit) = ocetra(i,j,k,igasnit) + an2odenit
           ocetra(i,j,k,ianh4)   = ocetra(i,j,k,ianh4)   + rnit*rnoxpi*(ano2denit+an2odenit) + rnh4dnra*rno2dnrai*ano2dnra
           ocetra(i,j,k,idet)    = ocetra(i,j,k,idet)    - (ano2denit + an2odenit)*rnoxpi - ano2dnra*rno2dnrai
           ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212) + rcar*rnoxpi*(ano2denit + an2odenit) + rcar*rno2dnrai*ano2dnra
           ocetra(i,j,k,iphosph) = ocetra(i,j,k,iphosph) + (ano2denit + an2odenit)*rnoxpi + ano2dnra*rno2dnrai
           ocetra(i,j,k,iiron)   = ocetra(i,j,k,iiron)   + riron*rnoxpi*(ano2denit + an2odenit) + riron*rno2dnrai*ano2dnra
           ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali) + (295.*ano2denit + rnm1*an2odenit)*rnoxpi &
                                 &                       + (rno2dnra + rnh4dnra - 1.)*rno2dnrai * ano2dnra
           ! Output
           denit_NO2(i,j,k) = ano2denit ! kmol NO2/m3/dtb - denitrification on NO2
           denit_N2O(i,j,k) = an2odenit ! kmol N2O/m3/dtb - denitrification on N2O
           DNRA_NO2(i,j,k)  = ano2dnra  ! kmol NO2/m3/dtb - DNRA on NO2
!          endif
        endif
      enddo
      enddo
      enddo
      !$OMP END PARALLEL DO
      end subroutine denit_dnra



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

      
