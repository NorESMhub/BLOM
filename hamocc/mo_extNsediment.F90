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
  use mo_param1_bgc,  only: issso12,ipowaic,ipowaal,ipowaph,ipowaox,ipown2,ipowno3,ipownh4,ipown2o,ipowno2,ks 
  use mo_vgrid,       only: kbo
  use mo_param_bgc,   only: rnit,rcar,rnoi
  use mo_control_bgc, only: io_stdo_bgc,dtb
  use mo_sedmnt,      only: powtra,sedlay,porsol,porwat
  use mo_extNwatercol,only: rc2n,ro2utammo,ro2nnit,rnoxp,rnoxpi,rno2anmx,rno2anmxi,rnh4anmx,                                       &
                         & rnh4anmxi,rno2dnra,rno2dnrai,rnh4dnra,rnh4dnrai,rnm1  

  implicit none

  private

  ! public functions
  public :: extNsediment_param_init,sed_nitrification,sed_denit_NO3_to_NO2,sed_anammox,sed_denit_DNRA,alloc_mem_extNsediment_diag, &
            extNsediment_param_update,extNsediment_param_write

  ! public parameters and fields
  public :: ised_nitr_NH4,ised_nitr_NO2,ised_nitr_N2O_prod,ised_nitr_NH4_OM,ised_nitr_NO2_OM,ised_denit_NO3,ised_denit_NO2,        &
            ised_denit_N2O,ised_DNRA_NO2,ised_anmx_N2_prod,ised_anmx_OM_prod,ised_remin_aerob,ised_remin_sulf,extNsed_diagnostics, &
            POM_remin_q10_sed, POM_remin_Tref_sed,bkox_drempoc_sed,                                                                &
            rano3denit_sed,rano2anmx_sed,rano2denit_sed,ran2odenit_sed,rdnra_sed,ranh4nitr_sed,rano2nitr_sed


  ! extended nitrogen cycle sediment parameters 
  real   ::  q10ano3denit_sed,sc_ano3denit_sed,Trefano3denit_sed,rano3denit_sed,bkano3denit_sed,                                   &
          & rano2anmx_sed,q10anmx_sed,Trefanmx_sed,alphaanmx_sed,bkoxanmx_sed,bkano2anmx_sed,bkanh4anmx_sed,                       &
          & rano2denit_sed,q10ano2denit_sed,Trefano2denit_sed,bkoxano2denit_sed,bkano2denit_sed,                                   &
          & ran2odenit_sed,q10an2odenit_sed,Trefan2odenit_sed,bkoxan2odenit_sed,bkan2odenit_sed,                                   &
          & rdnra_sed,q10dnra_sed,Trefdnra_sed,bkoxdnra_sed,bkdnra_sed,ranh4nitr_sed,q10anh4nitr_sed,                              &
          & Trefanh4nitr_sed,bkoxamox_sed,bkanh4nitr_sed,bkamoxn2o_sed,bkyamox_sed,                                                &
          & rano2nitr_sed,q10ano2nitr_sed,Trefano2nitr_sed,bkoxnitr_sed,bkano2nitr_sed,n2omaxy_sed,                                &
          & n2oybeta_sed,NOB2AOAy_sed,bn2o_sed,mufn2o_sed,POM_remin_q10_sed, POM_remin_Tref_sed,bkox_drempoc_sed 

  ! output
  real, dimension (:,:,:,:), allocatable :: extNsed_diagnostics
  integer, parameter :: &
             ised_nitr_NH4      = 1,  & 
             ised_nitr_NO2      = 2,  &
             ised_nitr_N2O_prod = 3,  &
             ised_nitr_NH4_OM   = 4,  &
             ised_nitr_NO2_OM   = 5,  &
             ised_denit_NO3     = 6,  &
             ised_denit_NO2     = 7,  &
             ised_denit_N2O     = 8,  &
             ised_DNRA_NO2      = 9,  &
             ised_anmx_N2_prod  = 10, &
             ised_anmx_OM_prod  = 11, &
             ised_remin_aerob   = 12, &
             ised_remin_sulf    = 13, &
             n_seddiag          = 13
  
  real :: eps,minlim

  contains

  ! ================================================================================================================================
  subroutine alloc_mem_extNsediment_diag(kpie,kpje,ksed)
  use mod_xc,         only: mnproc
  use mo_control_bgc, only: io_stdo_bgc
  
  implicit none

  INTEGER, intent(in) :: kpie,kpje,ksed ! ksed = ks
  INTEGER             :: errstat


    IF (mnproc.eq.1) THEN
     WRITE(io_stdo_bgc,*)'Memory allocation for sediment output of the extended nitrogen cycle ...'
     WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
     WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
     WRITE(io_stdo_bgc,*)'Third dimension    : ',ksed
     WRITE(io_stdo_bgc,*)'Fourth dimension   : ',n_seddiag
    ENDIF

    ALLOCATE (extNsed_diagnostics(kpie,kpje,ksed,n_seddiag),stat=errstat)

    if(errstat.ne.0) stop 'not enough memory extended nitrogen cycle'

  end subroutine alloc_mem_extNsediment_diag

  ! ================================================================================================================================
  subroutine extNsediment_param_init()
  use mo_extNwatercol,only: q10ano3denit,sc_ano3denit,Trefano3denit,bkano3denit,                                                   &
                         & q10anmx,Trefanmx,alphaanmx,bkoxanmx,bkano2anmx,                                                         &
                         & q10ano2denit,Trefano2denit,bkoxano2denit,bkano2denit,                                                   &
                         & q10an2odenit,Trefan2odenit,bkoxan2odenit,bkan2odenit,                                                   &
                         & q10dnra,Trefdnra,bkoxdnra,bkdnra,                                                                       &
                         & q10anh4nitr,Trefanh4nitr,bkoxamox,bkanh4nitr,bkamoxn2o,bkyamox,n2omaxy,n2oybeta,                        &
                         & q10ano2nitr,Trefano2nitr,bkoxnitr,bkano2nitr,NOB2AOAy,rno2anmx,rnh4anmx 
  use mo_param_bgc,   only: bkox_drempoc,POM_remin_q10,POM_remin_Tref

  implicit none
 
      ! === Ammonification in the sediment
      POM_remin_q10_sed  = POM_remin_q10  ! ammonification Q10 in sediment 
      POM_remin_Tref_sed = POM_remin_Tref ! ammonification Tref in sediment
      bkox_drempoc_sed   = bkox_drempoc   ! half saturation constant for O2 limitatio of ammonification in sediment

      ! === Denitrification step NO3 -> NO2:
      !rano3denit_sed    = 0.15*dtb       ! Maximum growth rate denitrification on NO3 at reference T (1/d -> 1/dt)
      rano3denit_sed    = 0.05           ! Maximum growth rate denitrification on NO3 at reference T (1/d -> 1/dt)
      q10ano3denit_sed  = q10ano3denit   ! Q10 factor for denitrification on NO3 (-)
      Trefano3denit_sed = Trefano3denit  ! Reference temperature for denitrification on NO3 (degr C) 
      !sc_ano3denit_sed  = 0.05e6         ! Shape factor for NO3 denitrification oxygen inhibition function (m3/kmol)
      sc_ano3denit_sed  = sc_ano3denit   ! Shape factor for NO3 denitrification oxygen inhibition function (m3/kmol)
      bkano3denit_sed   = bkano3denit    ! Half-saturation constant for NO3 denitrification (kmol/m3)

      ! === Anammox
      rano2anmx_sed     = 0.05           ! Maximum growth rate for anammox at reference T (1/d -> 1/dt)
      q10anmx_sed       = q10anmx        ! Q10 factor for anammox (-)
      Trefanmx_sed      = Trefanmx       ! Reference temperature for anammox (degr C)
      alphaanmx_sed     = alphaanmx      ! Shape factor for anammox oxygen inhibition function (m3/kmol)
      bkoxanmx_sed      = bkoxanmx       ! Half-saturation constant for oxygen inhibition function (kmol/m3)
      bkano2anmx_sed    = bkano2anmx     ! Half-saturation constant for NO2 limitation (kmol/m3)
      bkanh4anmx_sed    = bkano2anmx_sed * rnh4anmx/rno2anmx !Half-saturation constant for NH4 limitation of anammox (kmol/m3)

      ! === Denitrification step NO2 -> N2O
      rano2denit_sed    = 0.12           ! Maximum growth rate denitrification on NO2 at reference T (1/d -> 1/dt) 
      q10ano2denit_sed  = q10ano2denit   ! Q10 factor for denitrification on NO2 (-)
      Trefano2denit_sed = Trefano2denit  ! Reference temperature for denitrification on NO2 (degr C)
      bkoxano2denit_sed = bkoxano2denit  ! Half-saturation constant for (quadratic) oxygen inhibition function of denitrification on NO2 (kmol/m3)
      bkano2denit_sed   = bkano2denit    ! Half-saturation constant for denitrification on NO2 (kmol/m3)

      ! === Denitrification step N2O -> N2
      ran2odenit_sed    = 0.16           ! Maximum growth rate denitrification on N2O at reference T (1/d -> 1/dt)
      q10an2odenit_sed  = q10an2odenit   ! Q1- factor for denitrificationj on N2O (-)
      Trefan2odenit_sed = Trefan2odenit  ! Reference temperature for denitrification on N2O (degr C)
      bkoxan2odenit_sed = bkoxan2odenit  ! Half-saturation constant for (quadratic) oxygen inhibition function of denitrification on N2O (kmol/m3)
      bkan2odenit_sed   = bkan2odenit    ! Half-saturation constant for denitrification on N2O (kmol/m3)

      ! === DNRA NO2 -> NH4
      rdnra_sed         = 0.1            ! Maximum growth rate DNRA on NO2 at reference T (1/d -> 1/dt)
      q10dnra_sed       = q10dnra        ! Q10 factor for DNRA on NO2 (-)
      Trefdnra_sed      = Trefdnra       ! Reference temperature for DNRA (degr C)
      bkoxdnra_sed      = bkoxdnra       ! Half saturation constant for (quadratic) oxygen inhibition function of DNRA on NO2 (kmol/m3)
      bkdnra_sed        = bkdnra         ! Half-saturation constant for DNRA on NO2 (kmol/m3)

      ! === Nitrification on NH4
      ranh4nitr_sed     = 1.             ! Maximum growth rate nitrification on NH4 at reference T (1/d -> 1/dt) 
      q10anh4nitr_sed   = q10anh4nitr    ! Q10 factor for nitrification on NH4 (-)
      Trefanh4nitr_sed  = Trefanh4nitr   ! Reference temperature for nitrification on NH4 (degr C)
      bkoxamox_sed      = bkoxamox       ! Half-saturation constant for oxygen limitation of nitrification on NH4 (kmol/m3)
      bkanh4nitr_sed    = bkanh4nitr     ! Half-saturation constant for nitrification on NH4 (kmol/m3)
!======
! OLD VERSION OF pathway splitting function
      !bkamoxn2o_sed     = 0.453e-6 ! Half saturation constant for O2 in pathway splitting function N2O for nitrification on NH4 (kmol/m3)
! NEW version similar to Santoros 2021, Ji 2018:      
      bkamoxn2o_sed     = bkamoxn2o      ! Half saturation constant for NH4 in pathway splitting function N2O for nitrification on NH4 (kmol/m3)
      mufn2o_sed        = 0.11/(50.*1e6*bkoxamox_sed)  !=6.61e-3  0.11/(50*1e6)=2.2e-9 - ~Santoro et al. 2011 with simple MM,  
      bn2o_sed          = 0.077/(50.*mufn2o_sed)       !=0.2331 - before set to 0.3 - base fraction entering N2O 
!======
      !bkamoxno2_sed    = 0.479e-6       ! Half saturation constant for pathway splitting function N2O for nitrification on NH4 (kmol/m3)
    !  bkamoxno2_sed     = bkamoxno2      ! Half saturation constant for pathway splitting function N2O for nitrification on NH4 (kmol/m3)
      n2omaxy_sed       = n2omaxy        ! Maximum yield of OM on NH4 nitrification (-)
      n2oybeta_sed      = n2oybeta       ! Decay factor for inhibition function for yield during nitrification on NH4 (kmol/m3)
      bkyamox_sed       = bkyamox        ! Half saturation constant for pathway splitting function OM-yield for nitrification on NH4 (kmol/m3)
  
      ! === Nitrification on NO2
      rano2nitr_sed     = 1.54           ! Maximum growth rate nitrification on NO2 at reference T (1/d -> 1/dt) 
      q10ano2nitr_sed   = q10ano2nitr    ! Q10 factor for nitrification on NO2 (-)
      Trefano2nitr_sed  = Trefano2nitr   ! Reference temperature for nitrification on NO2 (degr C)
      bkoxnitr_sed      = bkoxnitr       ! Half-saturation constant for oxygen limitation of nitrification on NO2 (kmol/m3)
      bkano2nitr_sed    = bkano2nitr     ! Half-saturation constant for NO2 for nitrification on NO2 (kmol/m3)
      NOB2AOAy_sed      = NOB2AOAy       ! Ratio of NOB versus AOA yield per energy source ~0.043/0.098 according to Zakem et al. 2022

      eps    = 1.e-25 ! safe division etc. 
      minlim = 1.e-9  ! minimum for limitation functions (e.g. nutlim or oxlim/inh can only decrease to minlim) 
  end subroutine extNsediment_param_init

  ! ================================================================================================================================  
  subroutine extNsediment_param_update()

        rano3denit_sed = rano3denit_sed *dtb ! Maximum growth rate denitrification on NO3 at reference T (1/d -> 1/dt)
        rano2anmx_sed  = rano2anmx_sed  *dtb ! Maximum growth rate for anammox at reference T (1/d -> 1/dt)
        rano2denit_sed = rano2denit_sed *dtb ! Maximum growth rate denitrification on NO2 at reference T (1/d -> 1/dt) 
        ran2odenit_sed = ran2odenit_sed *dtb ! Maximum growth rate denitrification on N2O at reference T (1/d -> 1/dt)
        rdnra_sed      = rdnra_sed      *dtb ! Maximum growth rate DNRA on NO2 at reference T (1/d -> 1/dt)
        ranh4nitr_sed  = ranh4nitr_sed  *dtb ! Maximum growth rate nitrification on NH4 at reference T (1/d -> 1/dt) 
        rano2nitr_sed  = rano2nitr_sed  *dtb ! Maximum growth rate nitrification on NO2 at reference T (1/d -> 1/dt) 

  end subroutine extNsediment_param_update     
  
  ! ================================================================================================================================  
  subroutine extNsediment_param_write()

          REAL :: dtbinv
          dtbinv = 1./dtb

          WRITE(io_stdo_bgc,*) '****************************************************************'
          WRITE(io_stdo_bgc,*) '* HAMOCC extended nitrogen cycle parameters sediment:' 
          WRITE(io_stdo_bgc,*) '*          POM_remin_q10_sed = ', POM_remin_q10_sed
          WRITE(io_stdo_bgc,*) '*          POM_remin_Tref_sed= ', POM_remin_Tref_sed   
          WRITE(io_stdo_bgc,*) '*          bkox_drempoc_sed  = ', bkox_drempoc_sed   
          WRITE(io_stdo_bgc,*) '*          rano3denit_sed    = ',rano3denit_sed    *dtbinv
          WRITE(io_stdo_bgc,*) '*          q10ano3denit_sed  = ',q10ano3denit_sed
          WRITE(io_stdo_bgc,*) '*          Trefano3denit_sed = ',Trefano3denit_sed
          WRITE(io_stdo_bgc,*) '*          sc_ano3denit_sed  = ',sc_ano3denit_sed
          WRITE(io_stdo_bgc,*) '*          bkano3denit_sed   = ',bkano3denit_sed
          WRITE(io_stdo_bgc,*) '*          rano2anmx_sed     = ',rano2anmx_sed     *dtbinv
          WRITE(io_stdo_bgc,*) '*          q10anmx_sed       = ',q10anmx_sed
          WRITE(io_stdo_bgc,*) '*          Trefanmx_sed      = ',Trefanmx_sed
          WRITE(io_stdo_bgc,*) '*          alphaanmx_sed     = ',alphaanmx_sed
          WRITE(io_stdo_bgc,*) '*          bkoxanmx_sed      = ',bkoxanmx_sed
          WRITE(io_stdo_bgc,*) '*          bkano2anmx_sed    = ',bkano2anmx_sed
          WRITE(io_stdo_bgc,*) '*          bkanh4anmx_sed    = ',bkanh4anmx_sed
          WRITE(io_stdo_bgc,*) '*          rano2denit_sed    = ',rano2denit_sed    *dtbinv
          WRITE(io_stdo_bgc,*) '*          q10ano2denit_sed  = ',q10ano2denit_sed
          WRITE(io_stdo_bgc,*) '*          Trefano2denit_sed = ',Trefano2denit_sed
          WRITE(io_stdo_bgc,*) '*          bkoxano2denit_sed = ',bkoxano2denit_sed
          WRITE(io_stdo_bgc,*) '*          bkano2denit_sed   = ',bkano2denit_sed
          WRITE(io_stdo_bgc,*) '*          ran2odenit_sed    = ',ran2odenit_sed    *dtbinv
          WRITE(io_stdo_bgc,*) '*          q10an2odenit_sed  = ',q10an2odenit_sed
          WRITE(io_stdo_bgc,*) '*          Trefan2odenit_sed = ',Trefan2odenit_sed
          WRITE(io_stdo_bgc,*) '*          bkoxan2odenit_sed = ',bkoxan2odenit_sed
          WRITE(io_stdo_bgc,*) '*          bkan2odenit_sed   = ',bkan2odenit_sed
          WRITE(io_stdo_bgc,*) '*          rdnra_sed         = ',rdnra_sed         *dtbinv
          WRITE(io_stdo_bgc,*) '*          q10dnra_sed       = ',q10dnra_sed
          WRITE(io_stdo_bgc,*) '*          Trefdnra_sed      = ',Trefdnra_sed
          WRITE(io_stdo_bgc,*) '*          bkoxdnra_sed      = ',bkoxdnra_sed
          WRITE(io_stdo_bgc,*) '*          bkdnra_sed        = ',bkdnra_sed
          WRITE(io_stdo_bgc,*) '*          ranh4nitr_sed     = ',ranh4nitr_sed     *dtbinv
          WRITE(io_stdo_bgc,*) '*          q10anh4nitr_sed   = ',q10anh4nitr_sed
          WRITE(io_stdo_bgc,*) '*          Trefanh4nitr_sed  = ',Trefanh4nitr_sed
          WRITE(io_stdo_bgc,*) '*          bkoxamox_sed      = ',bkoxamox_sed
          WRITE(io_stdo_bgc,*) '*          bkanh4nitr_sed    = ',bkanh4nitr_sed
          WRITE(io_stdo_bgc,*) '*          bkamoxn2o_sed     = ',bkamoxn2o_sed
          WRITE(io_stdo_bgc,*) '*          mufn2o_sed        = ',mufn2o_sed
          WRITE(io_stdo_bgc,*) '*          bn2o_sed          = ',bn2o_sed
          WRITE(io_stdo_bgc,*) '*          n2omaxy_sed       = ',n2omaxy_sed
          WRITE(io_stdo_bgc,*) '*          n2oybeta_sed      = ',n2oybeta_sed
          WRITE(io_stdo_bgc,*) '*          bkyamox_sed       = ',bkyamox_sed
          WRITE(io_stdo_bgc,*) '*          rano2nitr_sed     = ',rano2nitr_sed     *dtbinv
          WRITE(io_stdo_bgc,*) '*          q10ano2nitr_sed   = ',q10ano2nitr_sed
          WRITE(io_stdo_bgc,*) '*          Trefano2nitr_sed  = ',Trefano2nitr_sed
          WRITE(io_stdo_bgc,*) '*          bkoxnitr_sed      = ',bkoxnitr_sed
          WRITE(io_stdo_bgc,*) '*          bkano2nitr_sed    = ',bkano2nitr_sed
          WRITE(io_stdo_bgc,*) '*          NOB2AOAy_sed      = ',NOB2AOAy_sed

  end subroutine extNsediment_param_write

  ! ================================================================================================================================  
  subroutine sed_nitrification(j,kpie,kpje,kpke,kbnd,ptho,omask,ex_ddic,ex_dalk)
    integer, intent(in) :: j,kpie,kpje,kpke,kbnd
    real,    intent(in) :: omask(kpie,kpje)
    real,    intent(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)
    ! for calculation of pore water DIC and alkalinity changes [P-units]! 
    real,    intent(inout) :: ex_ddic(kpie,ks)
    real,    intent(inout) :: ex_dalk(kpie,ks)
   
    ! local variables
    integer :: i,k

    real    :: Tdepanh4,O2limanh4,nut1lim,anh4new,potdnh4amox,fdetamox,fno2,fn2o,ftotnh4 
    real    :: Tdepano2,O2limano2,nut2lim,ano2new,potdno2nitr,fdetnitr,no2fn2o,no2fno2,no2fdetamox
    real    :: amoxfrac,nitrfrac,totd,amox,nitr,temp,w2s

    do i = 1,kpie
    do k = 1,ks
       if(omask(i,j)>0.5) then
          potdnh4amox = 0.
          fn2o        = 0.
          fno2        = 0.
          fdetamox    = 0.
          potdno2nitr = 0.
          fdetnitr    = 0.
          w2s         = porwat(i,j,k) / porsol(i,j,k)

!         if(ocetra(i,j,k,ioxygen)>minlim_oxnh4 .and. ocetra(i,j,k,ianh4)>minlim_nh4)then
           temp = merge(ptho(i,j,kbo(i,j)),10.,ptho(i,j,kbo(i,j))<40.)
           ! Ammonium oxidation step of nitrification
           Tdepanh4    = q10anh4nitr_sed**((temp-Trefanh4nitr_sed)/10.) 
           O2limanh4   = powtra(i,j,k,ipowaox)/(powtra(i,j,k,ipowaox) + bkoxamox_sed)
           nut1lim     = powtra(i,j,k,ipownh4)/(powtra(i,j,k,ipownh4) + bkanh4nitr_sed)
           anh4new     = powtra(i,j,k,ipownh4)/(1. + ranh4nitr_sed*Tdepanh4*O2limanh4*nut1lim)
           potdnh4amox = max(0.,powtra(i,j,k,ipownh4) - anh4new)
            
            ! pathway splitting functions according to Goreau 1980
       !=====
       ! OLD version according to Goreau
            !fn2o     = 1. - ocetra(i,j,k,ioxygen)/(ocetra(i,j,k,ioxygen) + bkamoxn2o)
       ! NEW version similar to Santoros et al. 2021, Ji et al. 2018
           fn2o     = mufn2o_sed * (bn2o_sed + (1.-bn2o_sed)*bkoxamox_sed/(powtra(i,j,k,ipowaox)+bkoxamox_sed))                    &
                    &        * powtra(i,j,k,ipownh4)/(powtra(i,j,k,ipownh4)+bkamoxn2o_sed)
       !=====
           fno2     = powtra(i,j,k,ipowaox)/(powtra(i,j,k,ipowaox) + bkoxamox_sed)
           fdetamox = n2omaxy_sed*2.*(1. + n2oybeta_sed)*powtra(i,j,k,ipowaox)*bkyamox_sed                                         &
                     & /(powtra(i,j,k,ipowaox)**2 + 2.*powtra(i,j,k,ipowaox)*bkyamox_sed + bkyamox_sed**2)

           ! normalization of pathway splitting functions to sum=1
           ftotnh4  = fn2o + fno2 + fdetamox + eps
           fn2o     = fn2o/ftotnh4
           fno2     = fno2/ftotnh4
           fdetamox = 1. - (fn2o + fno2)
!          endif

!          if(ocetra(i,j,k,ioxygen)>minlim_oxno2 .and. ocetra(i,j,k,iano2)>minlim_no2)then
!           temp = merge(ptho(i,j,kbo(i,j)),10.,ptho(i,j,kbo(i,j))<40.)
           ! NO2 oxidizing step of nitrification
           Tdepano2    = q10ano2nitr_sed**((temp-Trefano2nitr_sed)/10.) 
           O2limano2   = powtra(i,j,k,ipowaox)/(powtra(i,j,k,ipowaox) + bkoxnitr_sed)
           nut2lim     = powtra(i,j,k,ipowno2)/(powtra(i,j,k,ipowno2) + bkano2nitr_sed)
           ano2new     = powtra(i,j,k,ipowno2)/(1. + rano2nitr_sed*Tdepano2*O2limano2*nut2lim)
           potdno2nitr = max(0.,powtra(i,j,k,ipowno2) - ano2new)

           ! pathway splitting functions for NO2 nitrification - assuming to be the same as for NH4
           ! but with reduced OM gain per used NO2 as energy source (in amox: NH4)
        
           no2fn2o     = mufn2o_sed * (bn2o_sed + (1.-bn2o_sed)*bkoxamox_sed/(powtra(i,j,k,ipowaox)+bkoxamox_sed))                 &
                       &        * powtra(i,j,k,ipownh4)/(powtra(i,j,k,ipownh4)+bkamoxn2o_sed)
           no2fno2     = powtra(i,j,k,ipowaox)/(powtra(i,j,k,ipowaox) + bkoxamox_sed)
           no2fdetamox = NOB2AOAy_sed*n2omaxy_sed*2.*(1. + n2oybeta_sed)*powtra(i,j,k,ipowaox)*bkyamox_sed                         &
                     & /(powtra(i,j,k,ipowaox)**2 + 2.*powtra(i,j,k,ipowaox)*bkyamox_sed + bkyamox_sed**2)

           fdetnitr = no2fdetamox/(no2fno2 + no2fn2o)   ! yield to energy usage ratio for NO2 -> ratio equals 16:x
!          endif

          ! limitation of the two processes through available nutrients, etc.
          totd     = potdnh4amox + potdno2nitr
          amoxfrac = potdnh4amox/(totd + eps)
          nitrfrac = 1. - amoxfrac
           
          ! Account for potential earlier changes in DIC and alkalinity in finiding the minimum 
          totd     = max(0.,                                                                                                       &
                   &   min(totd,                                                                                                   &
                   &       powtra(i,j,k,ipownh4)/(amoxfrac + fdetnitr*nitrfrac + eps),                                             & ! ammonium
                   &       (powtra(i,j,k,ipowaic)+ex_ddic(i,k))/(rc2n*(fdetamox*amoxfrac + fdetnitr*nitrfrac) + eps),              & ! CO2
                   &       powtra(i,j,k,ipowaph)/(rnoi*(fdetamox*amoxfrac + fdetnitr*nitrfrac) + eps),                             & ! PO4
                   &       powtra(i,j,k,ipowaox)                                                                                   &
                   &       /((1.5*fno2 + fn2o - ro2nnit*fdetamox)*amoxfrac + (0.5 - ro2nnit*fdetnitr)*nitrfrac + eps),             & ! O2
                   &       (powtra(i,j,k,ipowaal) + ex_dalk(i,k))                                                                  &
                   &       /((2.*fno2 + fn2o + rnm1*rnoi*fdetamox)*amoxfrac + (rnm1*rnoi*fdetnitr)*nitrfrac + eps)))                ! alkalinity
          amox     = amoxfrac*totd 
          nitr     = nitrfrac*totd

          powtra(i,j,k,ipownh4)   = powtra(i,j,k,ipownh4) - amox - fdetnitr*nitr
          powtra(i,j,k,ipown2o)   = powtra(i,j,k,ipown2o) + 0.5*fn2o*amox
          powtra(i,j,k,ipowno2)   = powtra(i,j,k,ipowno2) + fno2*amox - nitr
          powtra(i,j,k,ipowno3)   = powtra(i,j,k,ipowno3) + nitr
          sedlay(i,j,k,issso12)   = sedlay(i,j,k,issso12)  + rnoi*(fdetamox*amox + fdetnitr*nitr) * w2s
!          ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212) - rc2n*(fdetamox*amox + fdetnitr*nitr)
          powtra(i,j,k,ipowaph)   = powtra(i,j,k,ipowaph) - rnoi*(fdetamox*amox + fdetnitr*nitr)
!          ocetra(i,j,k,iiron)   = ocetra(i,j,k,iiron)   - riron*rnoi*(fdetamox*amox + fdetnitr*nitr)
          powtra(i,j,k,ipowaox)   = powtra(i,j,k,ipowaox) - (1.5*fno2 + fn2o - ro2nnit*fdetamox)*amox   &
                                &                       - (0.5 - ro2nnit*fdetnitr)*nitr
!          ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali) - (2.*fno2 + fn2o + rnm1*rnoi*fdetamox)*amox - rnm1*rnoi*fdetnitr*nitr

          ! update of DIC and alkalinity through ex_ddic and ex_dalk fields 
          ! at later stage, when undersaturation of CaCO3 has been calculted  
          ex_ddic(i,k) = ex_ddic(i,k) - rc2n*(fdetamox*amox + fdetnitr*nitr)
          ex_dalk(i,k) = ex_dalk(i,k) - (2.*fno2 + fn2o + rnm1*rnoi*fdetamox)*amox - rnm1*rnoi*fdetnitr*nitr

          ! output:
          extNsed_diagnostics(i,j,k,ised_nitr_NH4)      = amox
          extNsed_diagnostics(i,j,k,ised_nitr_NO2)      = nitr
          extNsed_diagnostics(i,j,k,ised_nitr_N2O_prod) = 0.5*fn2o*amox
          extNsed_diagnostics(i,j,k,ised_nitr_NH4_OM)   = rnoi*fdetamox*amox * w2s 
          extNsed_diagnostics(i,j,k,ised_nitr_NO2_OM)   = rnoi*fdetnitr*nitr * w2s
       endif
    enddo
    enddo 
  end subroutine sed_nitrification

  ! ================================================================================================================================
  subroutine sed_denit_NO3_to_NO2(j,kpie,kpje,kpke,kbnd,ptho,omask,ex_ddic,ex_dalk)
    integer, intent(in) :: j,kpie,kpje,kpke,kbnd
    real,    intent(in) :: omask(kpie,kpje)
    real,    intent(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)
    ! for calculation of pore water DIC and alkalinity changes [P-units]! 
    real,    intent(inout) :: ex_ddic(kpie,ks)
    real,    intent(inout) :: ex_dalk(kpie,ks)

    ! local variables
    integer :: i,k
    real    :: Tdep,O2inhib,nutlim,ano3new,ano3denit,temp,s2w

    do i = 1,kpie
    do k = 1,ks
       if(omask(i,j)>0.5) then
            s2w       = porsol(i,j,k) / porwat(i,j,k)
            temp      = merge(ptho(i,j,kbo(i,j)),10.,ptho(i,j,kbo(i,j))<40.)
            Tdep      = q10ano3denit_sed**((temp-Trefano3denit_sed)/10.) 
            O2inhib   = 1. - tanh(sc_ano3denit_sed*powtra(i,j,k,ipowaox)) 
            nutlim    = powtra(i,j,k,ipowno3)/(powtra(i,j,k,ipowno3) + bkano3denit_sed)
 
            ano3new   = powtra(i,j,k,ipowno3)/(1. + rano3denit_sed*Tdep*O2inhib*nutlim) 

            ano3denit = max(0.,min(powtra(i,j,k,ipowno3) - ano3new, sedlay(i,j,k,issso12)*rnoxp*s2w))

            powtra(i,j,k,ipowno3) = powtra(i,j,k,ipowno3) - ano3denit
            powtra(i,j,k,ipowno2) = powtra(i,j,k,ipowno2) + ano3denit
            sedlay(i,j,k,issso12) = sedlay(i,j,k,issso12) - ano3denit*rnoxpi/s2w
            powtra(i,j,k,ipownh4) = powtra(i,j,k,ipownh4) + ano3denit*rnit*rnoxpi
            powtra(i,j,k,ipowaph) = powtra(i,j,k,ipowaph) + ano3denit*rnoxpi
            !ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212) + ano3denit*rcar*rnoxpi
            !ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali) + ano3denit*rnm1*rnoxpi
 
            ! update of DIC and alkalinity through ex_ddic and ex_dalk fields 
            ! at later stage, when undersaturation of CaCO3 has been calculted  
            ex_ddic(i,k) = ex_ddic(i,k) + ano3denit*rcar*rnoxpi
            ex_dalk(i,k) = ex_dalk(i,k) + ano3denit*rnm1*rnoxpi

            ! Output:
            extNsed_diagnostics(i,j,k,ised_denit_NO3) = ano3denit
       endif
    enddo
    enddo 

  end subroutine sed_denit_NO3_to_NO2

  ! ================================================================================================================================
  subroutine sed_anammox(j,kpie,kpje,kpke,kbnd,ptho,omask,ex_ddic,ex_dalk)
    integer, intent(in) :: j,kpie,kpje,kpke,kbnd
    real,    intent(in) :: omask(kpie,kpje)
    real,    intent(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)
    ! for calculation of pore water DIC and alkalinity changes [P-units]! 
    real,    intent(inout) :: ex_ddic(kpie,ks)
    real,    intent(inout) :: ex_dalk(kpie,ks)

    ! local variables
    integer :: i,k
    real    :: Tdep,O2inhib,nut1lim,nut2lim,ano2new,ano2anmx,temp,w2s

    do i = 1,kpie
    do k = 1,ks
       if(omask(i,j)>0.5) then
         w2s         = porwat(i,j,k) / porsol(i,j,k)
         temp     = merge(ptho(i,j,kbo(i,j)),10.,ptho(i,j,kbo(i,j))<40.)
         Tdep     = q10anmx_sed**((temp-Trefanmx_sed)/10.)         
         O2inhib  = 1. - exp(alphaanmx_sed*(powtra(i,j,k,ipowaox)-bkoxanmx_sed))                                                   &
                  &       /(1.+ exp(alphaanmx_sed*(powtra(i,j,k,ipowaox)-bkoxanmx_sed))) 
         nut1lim  = powtra(i,j,k,ipowno2)/(powtra(i,j,k,ipowno2)+bkano2anmx_sed)
         nut2lim  = powtra(i,j,k,ipownh4)/(powtra(i,j,k,ipownh4)+bkanh4anmx_sed)

         ano2new  = powtra(i,j,k,ipowno2)/(1. + rano2anmx_sed*Tdep*O2inhib*nut1lim*nut2lim) 

         ! Account for former changes in DIC and alkalinity
         ano2anmx = max(0.,min(powtra(i,j,k,ipowno2) - ano2new, powtra(i,j,k,ipownh4)*rno2anmx*rnh4anmxi,                          &
                               (powtra(i,j,k,ipowaic)+ex_ddic(i,k))*rno2anmx/rcar, powtra(i,j,k,ipowaph)*rno2anmx,                 &
                               (powtra(i,j,k,ipowaal)+ex_dalk(i,k))*rno2anmx/rnm1))

         powtra(i,j,k,ipowno2) = powtra(i,j,k,ipowno2) - ano2anmx
         powtra(i,j,k,ipownh4) = powtra(i,j,k,ipownh4) - ano2anmx*rnh4anmx*rno2anmxi
         powtra(i,j,k,ipown2)  = powtra(i,j,k,ipown2)  + ano2anmx*(rnh4anmx-rnit)*rno2anmxi
         powtra(i,j,k,ipowno3) = powtra(i,j,k,ipowno3) + ano2anmx*rnoxp*rno2anmxi
         sedlay(i,j,k,issso12) = sedlay(i,j,k,issso12) + ano2anmx*rno2anmxi*w2s
         powtra(i,j,k,ipowaph) = powtra(i,j,k,ipowaph) - ano2anmx*rno2anmxi
!         ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212) - ano2anmx*rcar*rno2anmxi
!         ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali) - ano2anmx*rnm1*rno2anmxi
            
         ! update of DIC and alkalinity through ex_ddic and ex_dalk fields 
         ! at later stage, when undersaturation of CaCO3 has been calculted  
         ex_ddic(i,k) = ex_ddic(i,k) - ano2anmx*rcar*rno2anmxi
         ex_dalk(i,k) = ex_dalk(i,k) - ano2anmx*rnm1*rno2anmxi

         ! Output:
         extNsed_diagnostics(i,j,k,ised_anmx_N2_prod) = ano2anmx*(rnh4anmx-rnit)*rno2anmxi  ! kmol N2/m3/dtb - N2 prod through anammox
         extNsed_diagnostics(i,j,k,ised_anmx_OM_prod) = ano2anmx*rno2anmxi*w2s
       endif
    enddo
    enddo 

  end subroutine sed_anammox

  ! ================================================================================================================================
  subroutine sed_denit_DNRA(j,kpie,kpje,kpke,kbnd,ptho,omask,ex_ddic,ex_dalk)
    integer, intent(in) :: j,kpie,kpje,kpke,kbnd
    real,    intent(in) :: omask(kpie,kpje)
    real,    intent(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)
    ! for calculation of pore water DIC and alkalinity changes [P-units]! 
    real,    intent(inout) :: ex_ddic(kpie,ks)
    real,    intent(inout) :: ex_dalk(kpie,ks)
    
    ! local variables
    integer :: i,k
    real    :: Tdepano2,O2inhibano2,nutlimano2,rpotano2denit,ano2denit
    real    :: Tdepdnra,O2inhibdnra,nutlimdnra,rpotano2dnra,ano2dnra
    real    :: fdenit,fdnra,potano2new,potdano2,potddet,fdetano2denit,fdetan2odenit,fdetdnra  
    real    :: Tdepan2o,O2inhiban2o,nutliman2o,an2onew,an2odenit
    real    :: temp,s2w


    do i = 1,kpie
    do k = 1,ks
       if(omask(i,j)>0.5) then
          potddet       = 0.
          an2odenit     = 0.
          ano2denit     = 0.
          ano2dnra      = 0.
          s2w           =  porsol(i,j,k) / porwat(i,j,k)
!          if(0.<=ocetra(i,j,k,ioxygen) .and. ocetra(i,j,k,ioxygen)<minlim_oxn2o .and. ocetra(i,j,k,ian2o)>minlim_n2o)then
           temp = merge(ptho(i,j,kbo(i,j)),10.,ptho(i,j,kbo(i,j))<40.)
           ! === denitrification on N2O
           Tdepan2o    = q10an2odenit_sed**((temp-Trefan2odenit_sed)/10.) 
           O2inhiban2o = bkoxan2odenit_sed**2/(powtra(i,j,k,ipowaox)**2 + bkoxan2odenit_sed**2) 
           nutliman2o  = powtra(i,j,k,ipown2o)/(powtra(i,j,k,ipown2o) + bkan2odenit_sed)   
           an2onew     = powtra(i,j,k,ipown2o)/(1. + ran2odenit_sed*Tdepan2o*O2inhiban2o*nutliman2o)  
           an2odenit   = max(0.,min(powtra(i,j,k,ipown2o),powtra(i,j,k,ipown2o) - an2onew))
!          endif

!          if(0.<=ocetra(i,j,k,ioxygen) .and. ocetra(i,j,k,ioxygen)<minlim_ox .and. ocetra(i,j,k,iano2)>minlim_no2)then 
!           temp = merge(ptho(i,j,kbo(i,j)),10.,ptho(i,j,kbo(i,j))<40.)
           ! denitrification on NO2
           Tdepano2    = q10ano2denit_sed**((temp-Trefano2denit_sed)/10.) 
           O2inhibano2 = bkoxano2denit_sed**2/(powtra(i,j,k,ipowaox)**2 + bkoxano2denit_sed**2) 
           nutlimano2  = powtra(i,j,k,ipowno2)/(powtra(i,j,k,ipowno2) + bkano2denit_sed)
           rpotano2denit = max(0.,rano2denit_sed*Tdepano2*O2inhibano2*nutlimano2) ! potential rate of denit
           
           ! DNRA on NO2
           Tdepdnra    = q10dnra_sed**((temp-Trefdnra_sed)/10.) 
           O2inhibdnra = bkoxdnra_sed**2/(powtra(i,j,k,ipowaox)**2 + bkoxdnra_sed**2) 
           nutlimdnra  = powtra(i,j,k,ipowno2)/(powtra(i,j,k,ipowno2) + bkdnra_sed)
           rpotano2dnra = max(0.,rdnra_sed*Tdepdnra*O2inhibdnra*nutlimdnra) ! pot. rate of dnra

           ! potential new conc of NO2 due to denitrification and DNRA
           potano2new = powtra(i,j,k,ipowno2)/(1. + rpotano2denit + rpotano2dnra)
           potdano2   = max(0.,min(powtra(i,j,k,ipowno2), powtra(i,j,k,ipowno2) - potano2new))
           
           ! === limitation due to NO2:
           ! fraction on potential change of NO2:
           fdenit = rpotano2denit/(rpotano2denit + rpotano2dnra + eps)
           fdnra  = 1. - fdenit
           
           ! potential fractional change
           ano2denit  = fdenit * potdano2   
           ano2dnra   = fdnra  * potdano2
 !         endif

          ! limitation of processes due to detritus (based on pore water volume)
          potddet       = rnoxpi*(ano2denit + an2odenit) + rno2dnrai*ano2dnra  ! P units              
          fdetano2denit = rnoxpi*ano2denit/(potddet + eps)
          fdetan2odenit = rnoxpi*an2odenit/(potddet + eps)
          fdetdnra      = 1. - fdetano2denit - fdetan2odenit 
          potddet       = max(0.,min(potddet,powtra(i,j,k,issso12)*s2w)) 
       
!          if(potddet>0.)then
           ! change of NO2 and N2O in N units
           ano2denit     = fdetano2denit*rnoxp*potddet
           an2odenit     = fdetan2odenit*rnoxp*potddet
           ano2dnra      = fdetdnra*rno2dnra*potddet

           ! change in tracer concentrations due to denit (NO2->N2O->N2) and DNRA (NO2->NH4)
           powtra(i,j,k,ipowno2) = powtra(i,j,k,ipowno2) - ano2denit - ano2dnra
           powtra(i,j,k,ipown2o) = powtra(i,j,k,ipown2o) - an2odenit + 0.5*ano2denit
           powtra(i,j,k,ipown2)  = powtra(i,j,k,ipown2)  + an2odenit
           powtra(i,j,k,ipownh4) = powtra(i,j,k,ipownh4) + rnit*rnoxpi*(ano2denit+an2odenit) + rnh4dnra*rno2dnrai*ano2dnra
           sedlay(i,j,k,issso12) = sedlay(i,j,k,issso12) - ((ano2denit + an2odenit)*rnoxpi + ano2dnra*rno2dnrai)/s2w
           powtra(i,j,k,ipowaph) = powtra(i,j,k,ipowaph) + (ano2denit + an2odenit)*rnoxpi + ano2dnra*rno2dnrai
!           ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212) + rcar*rnoxpi*(ano2denit + an2odenit) + rcar*rno2dnrai*ano2dnra
!           ocetra(i,j,k,iiron)   = ocetra(i,j,k,iiron)   + riron*rnoxpi*(ano2denit + an2odenit) + riron*rno2dnrai*ano2dnra
!           ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali) + (295.*ano2denit + rnm1*an2odenit)*rnoxpi &
!                                 &                       + (rno2dnra + rnh4dnra - 1.)*rno2dnrai * ano2dnra

         ! update of DIC and alkalinity through ex_ddic and ex_dalk fields 
         ! at later stage, when undersaturation of CaCO3 has been calculted  
         ex_ddic(i,k) = ex_ddic(i,k) + rcar*rnoxpi*(ano2denit + an2odenit) + rcar*rno2dnrai*ano2dnra
         ex_dalk(i,k) = ex_dalk(i,k) + (295.*ano2denit + rnm1*an2odenit)*rnoxpi + (rno2dnra + rnh4dnra - 1.)*rno2dnrai * ano2dnra

         extNsed_diagnostics(i,j,k,ised_denit_NO2) = ano2denit
         extNsed_diagnostics(i,j,k,ised_denit_N2O) = an2odenit
         extNsed_diagnostics(i,j,k,ised_DNRA_NO2)  = ano2dnra
       endif
    enddo
    enddo 

  end subroutine sed_denit_DNRA

END MODULE mo_extNsediment
