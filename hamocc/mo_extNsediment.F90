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

module mo_extNsediment
  !**********************************************************************
  !
  ! MODULE mo_extNsediment - extended nitrogen cycle processes
  !                          in the sediment
  !
  ! j.maerz 13.09.2022
  !
  ! Pupose:
  ! -------
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
  use mo_param1_bgc,  only: issso12,ipowaic,ipowaal,ipowaph,ipowaox,ipown2,ipowno3,ipownh4,ipown2o,&
                          & ipowno2,ks
  use mo_vgrid,       only: kbo
  use mo_param_bgc,   only: rnit,rcar,rnoi,                                                        &
                          & rc2n,ro2utammo,ro2nnit,rnoxp,rnoxpi,rno2anmx,rno2anmxi,rnh4anmx,       &
                          & rnh4anmxi,rno2dnra,rno2dnrai,rnh4dnra,rnh4dnrai,rnm1,                  &
                          & q10ano3denit_sed,sc_ano3denit_sed,Trefano3denit_sed,rano3denit_sed,    &
                          & bkano3denit_sed,rano2anmx_sed,q10anmx_sed,Trefanmx_sed,alphaanmx_sed,  &
                          & bkoxanmx_sed,bkano2anmx_sed,bkanh4anmx_sed,rano2denit_sed,             &
                          & q10ano2denit_sed,Trefano2denit_sed,bkoxano2denit_sed,bkano2denit_sed,  &
                          & ran2odenit_sed,q10an2odenit_sed,Trefan2odenit_sed,bkoxan2odenit_sed,   &
                          & bkan2odenit_sed,rdnra_sed,q10dnra_sed,Trefdnra_sed,bkoxdnra_sed,       &
                          & bkdnra_sed,ranh4nitr_sed,q10anh4nitr_sed,Trefanh4nitr_sed,bkoxamox_sed,&
                          & bkanh4nitr_sed,bkamoxn2o_sed,bkyamox_sed,rano2nitr_sed,q10ano2nitr_sed,&
                          & Trefano2nitr_sed,bkoxnitr_sed,bkano2nitr_sed,n2omaxy_sed,n2oybeta_sed, &
                          & NOB2AOAy_sed,bn2o_sed,mufn2o_sed,POM_remin_q10_sed, POM_remin_Tref_sed,&
                          & bkox_drempoc_sed,max_limiter
  use mo_control_bgc, only: io_stdo_bgc
  use mo_sedmnt,      only: powtra,sedlay,porsol,porwat

  implicit none

  private

  ! public functions
  public :: sed_nitrification,sed_denit_NO3_to_NO2,sed_anammox,sed_denit_DNRA,                     &
          & alloc_mem_extNsediment_diag

  ! public parameters and fields
  public :: ised_nitr_NH4,ised_nitr_NO2,ised_nitr_N2O_prod,ised_nitr_NH4_OM,ised_nitr_NO2_OM,      &
          & ised_denit_NO3,ised_denit_NO2,ised_denit_N2O,ised_DNRA_NO2,ised_anmx_N2_prod,          &
          & ised_anmx_OM_prod,ised_remin_aerob,ised_remin_sulf,extNsed_diagnostics

  ! output
  real, dimension (:,:,:,:), allocatable :: extNsed_diagnostics
  integer, parameter ::               &
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

  real :: eps    = 1.e-25
  real :: minlim = 1.e-9

contains

  ! ================================================================================================================================
  subroutine alloc_mem_extNsediment_diag(kpie,kpje,ksed)
    use mod_xc,         only: mnproc
    use mo_control_bgc, only: io_stdo_bgc

    implicit none

    integer, intent(in) :: kpie,kpje,ksed ! ksed = ks

    integer             :: errstat

    if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for sediment output of the extended nitrogen cycle ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',ksed
      write(io_stdo_bgc,*)'Fourth dimension   : ',n_seddiag
    endif

    allocate (extNsed_diagnostics(kpie,kpje,ksed,n_seddiag),stat=errstat)

    if(errstat.ne.0) stop 'not enough memory extended nitrogen cycle'
  end subroutine alloc_mem_extNsediment_diag

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
        if (omask(i,j) > 0.5) then
          potdnh4amox = 0.
          fn2o        = 0.
          fno2        = 0.
          fdetamox    = 0.
          potdno2nitr = 0.
          fdetnitr    = 0.
          w2s         = porwat(i,j,k) / porsol(i,j,k)

          temp = merge(ptho(i,j,kbo(i,j)),10.,ptho(i,j,kbo(i,j))<40.)
          ! Ammonium oxidation step of nitrification
          Tdepanh4    = q10anh4nitr_sed**((temp-Trefanh4nitr_sed)/10.)
          O2limanh4   = powtra(i,j,k,ipowaox)/(powtra(i,j,k,ipowaox) + bkoxamox_sed)
          nut1lim     = powtra(i,j,k,ipownh4)/(powtra(i,j,k,ipownh4) + bkanh4nitr_sed)
          anh4new     = powtra(i,j,k,ipownh4)/(1. + ranh4nitr_sed*Tdepanh4*O2limanh4*nut1lim)
          potdnh4amox = max(0.,powtra(i,j,k,ipownh4) - anh4new)

          ! pathway splitting functions similar to Santoros et al. 2021, Ji et al. 2018
          fn2o     = mufn2o_sed * (bn2o_sed + (1.-bn2o_sed)*bkoxamox_sed                           &
                   &                          /(powtra(i,j,k,ipowaox)+bkoxamox_sed))               &
                   &            * powtra(i,j,k,ipownh4)/(powtra(i,j,k,ipownh4)+bkamoxn2o_sed)

          fno2     = powtra(i,j,k,ipowaox)/(powtra(i,j,k,ipowaox) + bkoxamox_sed)
          fdetamox = n2omaxy_sed*2.*(1. + n2oybeta_sed)*powtra(i,j,k,ipowaox)*bkyamox_sed          &
                   & /(powtra(i,j,k,ipowaox)**2 + 2.*powtra(i,j,k,ipowaox)*bkyamox_sed + bkyamox_sed**2)

          ! normalization of pathway splitting functions to sum=1
          ftotnh4  = fn2o + fno2 + fdetamox + eps
          fn2o     = fn2o/ftotnh4
          fno2     = fno2/ftotnh4
          fdetamox = 1. - (fn2o + fno2)

          ! NO2 oxidizing step of nitrification
          Tdepano2    = q10ano2nitr_sed**((temp-Trefano2nitr_sed)/10.) 
          O2limano2   = powtra(i,j,k,ipowaox)/(powtra(i,j,k,ipowaox) + bkoxnitr_sed)
          nut2lim     = powtra(i,j,k,ipowno2)/(powtra(i,j,k,ipowno2) + bkano2nitr_sed)
          ano2new     = powtra(i,j,k,ipowno2)/(1. + rano2nitr_sed*Tdepano2*O2limano2*nut2lim)
          potdno2nitr = max(0.,powtra(i,j,k,ipowno2) - ano2new)

          ! pathway splitting functions for NO2 nitrification - assuming to be the same as for NH4
          ! but with reduced OM gain per used NO2 as energy source (in amox: NH4)
          no2fn2o     = mufn2o_sed * (bn2o_sed + (1.-bn2o_sed)*bkoxamox_sed                        &
                      &                          /(powtra(i,j,k,ipowaox)+bkoxamox_sed))            &
                      &        * powtra(i,j,k,ipownh4)/(powtra(i,j,k,ipownh4)+bkamoxn2o_sed)
          no2fno2     = powtra(i,j,k,ipowaox)/(powtra(i,j,k,ipowaox) + bkoxamox_sed)
          no2fdetamox = NOB2AOAy_sed*n2omaxy_sed*2.*(1. + n2oybeta_sed)*powtra(i,j,k,ipowaox)*bkyamox_sed  &
                      & /(powtra(i,j,k,ipowaox)**2 + 2.*powtra(i,j,k,ipowaox)*bkyamox_sed + bkyamox_sed**2)

          fdetnitr = no2fdetamox/(no2fno2 + no2fn2o)   ! yield to energy usage ratio for NO2 -> ratio equals 16:x

          ! limitation of the two processes through available nutrients, etc.
          totd     = potdnh4amox + potdno2nitr
          amoxfrac = potdnh4amox/(totd + eps)
          nitrfrac = 1. - amoxfrac

          ! Account for potential earlier changes in DIC and alkalinity in finiding the minimum
          totd  = max(0.,                                                                          &
                &    min(totd,                                                                     &
                &       max_limiter*powtra(i,j,k,ipownh4)/(amoxfrac + fdetnitr*nitrfrac + eps),    & ! ammonium
                &       max_limiter*(powtra(i,j,k,ipowaic) + ex_ddic(i,k))                         &
                &                                /(rc2n*(fdetamox*amoxfrac + fdetnitr*nitrfrac)    &
                &                                  + eps),                                         & ! CO2
                &       max_limiter*powtra(i,j,k,ipowaph)                                          &
                &                            /(rnoi*(fdetamox*amoxfrac+fdetnitr*nitrfrac) + eps),  & ! PO4
                &       max_limiter*powtra(i,j,k,ipowaox)                                          &
                &       /((1.5*fno2 + fn2o - ro2nnit*fdetamox)*amoxfrac                            &
                &                          + (0.5 - ro2nnit*fdetnitr)*nitrfrac + eps),             & ! O2
                &       max_limiter*(powtra(i,j,k,ipowaal) + ex_dalk(i,k))                         &
                &       /((2.*fno2 + fn2o + rnm1*rnoi*fdetamox)*amoxfrac                           &
                &                         + (rnm1*rnoi*fdetnitr)*nitrfrac + eps)))                   ! alkalinity
          amox  = amoxfrac*totd
          nitr  = nitrfrac*totd

          powtra(i,j,k,ipownh4) = powtra(i,j,k,ipownh4) - amox - fdetnitr*nitr
          powtra(i,j,k,ipown2o) = powtra(i,j,k,ipown2o) + 0.5*fn2o*amox
          powtra(i,j,k,ipowno2) = powtra(i,j,k,ipowno2) + fno2*amox - nitr
          powtra(i,j,k,ipowno3) = powtra(i,j,k,ipowno3) + nitr
          sedlay(i,j,k,issso12) = sedlay(i,j,k,issso12) + rnoi*(fdetamox*amox + fdetnitr*nitr)*w2s
          powtra(i,j,k,ipowaph) = powtra(i,j,k,ipowaph) - rnoi*(fdetamox*amox + fdetnitr*nitr)
          powtra(i,j,k,ipowaox) = powtra(i,j,k,ipowaox) - (1.5*fno2 + fn2o - ro2nnit*fdetamox)*amox&
                                &                       - (0.5 - ro2nnit*fdetnitr)*nitr

          ! update of DIC and alkalinity through ex_ddic and ex_dalk fields
          ! at later stage, when undersaturation of CaCO3 has been calculted
          ex_ddic(i,k) = ex_ddic(i,k) - rc2n*(fdetamox*amox + fdetnitr*nitr)
          ex_dalk(i,k) = ex_dalk(i,k) - (2.*fno2 + fn2o + rnm1*rnoi*fdetamox)*amox                 &
                                      - rnm1*rnoi*fdetnitr*nitr

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
        if (omask(i,j) > 0.5) then
          s2w       = porsol(i,j,k) / porwat(i,j,k)
          temp      = merge(ptho(i,j,kbo(i,j)),10.,ptho(i,j,kbo(i,j)) < 40.)
          Tdep      = q10ano3denit_sed**((temp-Trefano3denit_sed)/10.)
          O2inhib   = 1. - tanh(sc_ano3denit_sed*powtra(i,j,k,ipowaox))
          nutlim    = powtra(i,j,k,ipowno3)/(powtra(i,j,k,ipowno3) + bkano3denit_sed)

          ano3new   = powtra(i,j,k,ipowno3)/(1. + rano3denit_sed*Tdep*O2inhib*nutlim)

          ano3denit = max(0.,min(powtra(i,j,k,ipowno3) - ano3new,                                  &
                    &                      max_limiter*sedlay(i,j,k,issso12)*rnoxp*s2w))

          powtra(i,j,k,ipowno3) = powtra(i,j,k,ipowno3) - ano3denit
          powtra(i,j,k,ipowno2) = powtra(i,j,k,ipowno2) + ano3denit
          sedlay(i,j,k,issso12) = sedlay(i,j,k,issso12) - ano3denit*rnoxpi/s2w
          powtra(i,j,k,ipownh4) = powtra(i,j,k,ipownh4) + ano3denit*rnit*rnoxpi
          powtra(i,j,k,ipowaph) = powtra(i,j,k,ipowaph) + ano3denit*rnoxpi

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
          temp     = merge(ptho(i,j,kbo(i,j)),10.,ptho(i,j,kbo(i,j)) < 40.)
          Tdep     = q10anmx_sed**((temp-Trefanmx_sed)/10.)
          O2inhib  = 1. - exp(alphaanmx_sed*(powtra(i,j,k,ipowaox)-bkoxanmx_sed))                  &
                   &      /(1.+ exp(alphaanmx_sed*(powtra(i,j,k,ipowaox)-bkoxanmx_sed)))
          nut1lim  = powtra(i,j,k,ipowno2)/(powtra(i,j,k,ipowno2)+bkano2anmx_sed)
          nut2lim  = powtra(i,j,k,ipownh4)/(powtra(i,j,k,ipownh4)+bkanh4anmx_sed)

          ano2new  = powtra(i,j,k,ipowno2)/(1. + rano2anmx_sed*Tdep*O2inhib*nut1lim*nut2lim)

          ! Account for former changes in DIC and alkalinity
          ano2anmx = max(0.,min(max_limiter*powtra(i,j,k,ipowno2) - ano2new,                       &
                                max_limiter*powtra(i,j,k,ipownh4)*rno2anmx*rnh4anmxi,              &
                                max_limiter*(powtra(i,j,k,ipowaic)+ex_ddic(i,k))*rno2anmx/rcar,    &
                                max_limiter*powtra(i,j,k,ipowaph)*rno2anmx,                        &
                                max_limiter*(powtra(i,j,k,ipowaal)+ex_dalk(i,k))*rno2anmx/rnm1))

          powtra(i,j,k,ipowno2) = powtra(i,j,k,ipowno2) - ano2anmx
          powtra(i,j,k,ipownh4) = powtra(i,j,k,ipownh4) - ano2anmx*rnh4anmx*rno2anmxi
          powtra(i,j,k,ipown2)  = powtra(i,j,k,ipown2)  + ano2anmx*(rnh4anmx-rnit)*rno2anmxi
          powtra(i,j,k,ipowno3) = powtra(i,j,k,ipowno3) + ano2anmx*rnoxp*rno2anmxi
          sedlay(i,j,k,issso12) = sedlay(i,j,k,issso12) + ano2anmx*rno2anmxi*w2s
          powtra(i,j,k,ipowaph) = powtra(i,j,k,ipowaph) - ano2anmx*rno2anmxi

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
        if (omask(i,j) > 0.5) then
          potddet       = 0.
          an2odenit     = 0.
          ano2denit     = 0.
          ano2dnra      = 0.
          s2w           =  porsol(i,j,k) / porwat(i,j,k)
          temp          = merge(ptho(i,j,kbo(i,j)),10.,ptho(i,j,kbo(i,j)) < 40.)

           ! === denitrification on N2O
          Tdepan2o    = q10an2odenit_sed**((temp-Trefan2odenit_sed)/10.)
          O2inhiban2o = bkoxan2odenit_sed**2/(powtra(i,j,k,ipowaox)**2 + bkoxan2odenit_sed**2)
          nutliman2o  = powtra(i,j,k,ipown2o)/(powtra(i,j,k,ipown2o) + bkan2odenit_sed)
          an2onew     = powtra(i,j,k,ipown2o)/(1. + ran2odenit_sed*Tdepan2o*O2inhiban2o*nutliman2o)
          an2odenit   = max(0.,min(powtra(i,j,k,ipown2o),powtra(i,j,k,ipown2o) - an2onew))

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

          ! limitation of processes due to detritus (based on pore water volume)
          potddet       = rnoxpi*(ano2denit + an2odenit) + rno2dnrai*ano2dnra  ! P units
          fdetano2denit = rnoxpi*ano2denit/(potddet + eps)
          fdetan2odenit = rnoxpi*an2odenit/(potddet + eps)
          fdetdnra      = 1. - fdetano2denit - fdetan2odenit
          potddet       = max(0.,min(potddet,max_limiter*sedlay(i,j,k,issso12)*s2w))

          ! change of NO2 and N2O in N units
          ano2denit     = fdetano2denit*rnoxp*potddet
          an2odenit     = fdetan2odenit*rnoxp*potddet
          ano2dnra      = fdetdnra*rno2dnra*potddet

          ! change in tracer concentrations due to denit (NO2->N2O->N2) and DNRA (NO2->NH4)
          powtra(i,j,k,ipowno2) = powtra(i,j,k,ipowno2) - ano2denit - ano2dnra
          powtra(i,j,k,ipown2o) = powtra(i,j,k,ipown2o) - an2odenit + 0.5*ano2denit
          powtra(i,j,k,ipown2)  = powtra(i,j,k,ipown2)  + an2odenit
          powtra(i,j,k,ipownh4) = powtra(i,j,k,ipownh4) + rnit*rnoxpi*(ano2denit+an2odenit)        &
                                &                       + rnh4dnra*rno2dnrai*ano2dnra
          sedlay(i,j,k,issso12) = sedlay(i,j,k,issso12)                                            &
                                &      - ((ano2denit + an2odenit)*rnoxpi + ano2dnra*rno2dnrai)/s2w
          powtra(i,j,k,ipowaph) = powtra(i,j,k,ipowaph) + (ano2denit + an2odenit)*rnoxpi           &
                                &                       + ano2dnra*rno2dnrai

          ! update of DIC and alkalinity through ex_ddic and ex_dalk fields
          ! at later stage, when undersaturation of CaCO3 has been calculted
          ex_ddic(i,k)          = ex_ddic(i,k) + rcar*rnoxpi*(ano2denit + an2odenit)               &
                                &              + rcar*rno2dnrai*ano2dnra
          ex_dalk(i,k)          = ex_dalk(i,k) + (295.*ano2denit + rnm1*an2odenit)*rnoxpi          &
                                &                 + (rno2dnra + rnh4dnra - 1.)*rno2dnrai*ano2dnra

          extNsed_diagnostics(i,j,k,ised_denit_NO2) = ano2denit
          extNsed_diagnostics(i,j,k,ised_denit_N2O) = an2odenit
          extNsed_diagnostics(i,j,k,ised_DNRA_NO2)  = ano2dnra
        endif
      enddo
    enddo
  end subroutine sed_denit_DNRA

end module mo_extNsediment
