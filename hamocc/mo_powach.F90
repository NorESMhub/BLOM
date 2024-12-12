! Copyright (C) 2001  Ernst Maier-Reimer, S. Legutke
! Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger
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

module mo_powach

  implicit none
  private

  public :: powach

contains

  subroutine powach(kpie,kpje,kpke,kbnd,prho,omask,psao,ptho,lspin)

    !***********************************************************************************************
    ! Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
    ! Modified: S.Legutke,   *MPI-MaD, HH*    10.04.01
    !***********************************************************************************************

    use mo_control_bgc, only: dtbgc,use_cisonew,use_extNcycle,lTO2depremin
    use mo_param1_bgc,  only: ioxygen,ipowaal,ipowaic,ipowaox,ipowaph,ipowasi,ipown2,ipowno3,      &
                              isilica,isssc12,issso12,issssil,issster,ks,ipowc13,ipowc14,isssc13,  &
                              isssc14,issso13,issso14,safediv,ipownh4
    use mo_carbch,      only: co3,keqb,ocetra,sedfluxo
    use mo_chemcon,     only: calcon
    use mo_param_bgc,   only: rnit,rcar,rdnit1,rdnit2,ro2ut,disso_sil,silsat,disso_poc,sed_denit,  &
                            & disso_caco3,ro2utammo,                                               &
                            & POM_remin_q10_sed,POM_remin_Tref_sed,bkox_drempoc_sed
    use mo_sedmnt,      only: porwat,porsol,powtra,produs,prcaca,prorca,seddw,sedhpl,sedlay,       &
                              silpro,pror13,pror14,prca13,prca14
    use mo_vgrid,       only: kbo,bolay
    use mo_powadi,      only: powadi
    use mo_carchm,      only: carchm_solve
    use mo_dipowa,      only: dipowa
    use mo_extNsediment,only: sed_nitrification,sed_denit_NO3_to_NO2,sed_anammox,sed_denit_DNRA,   &
                            & extNsed_diagnostics,ised_remin_aerob,ised_remin_sulf

    ! Arguments
    integer, intent(in) :: kpie                                         ! 1st dimension of model grid.
    integer, intent(in) :: kpje                                         ! 2nd dimension of model grid.
    integer, intent(in) :: kpke                                         ! 3rd (vertical) dimension of model grid.
    integer, intent(in) :: kbnd                                         ! nb of halo grid points
    real,    intent(in) :: prho(kpie,kpje,kpke)                         ! seawater density [g/cm^3].
    real,    intent(in) :: omask(kpie,kpje)                             ! land/ocean mask.
    real,    intent(in) :: psao(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) ! salinity [psu].
    real,    intent(in) :: ptho(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke) ! Pot. temperature [deg C].
    logical, intent(in) :: lspin

    ! Local variables
    integer :: i,j,k,l
    real    :: sedb1(kpie,0:ks),sediso(kpie,0:ks)
    real    :: solrat(kpie,ks),powcar(kpie,ks)
    real    :: aerob(kpie,ks),anaerob(kpie,ks),sulf(kpie,ks)
    real    :: ex_ddic(kpie,ks),ex_dalk(kpie,ks) !sum of DIC and alk changes related to extended nitrogen cycle
    real    :: ex_disso_poc
    real    :: aerob13(kpie,ks),anaerob13(kpie,ks),sulf13(kpie,ks) ! cisonew
    real    :: aerob14(kpie,ks),anaerob14(kpie,ks),sulf14(kpie,ks) ! cisonew
    real    :: dissot, undsa, posol
    real    :: umfa, denit, saln, rrho, alk, c, sit, pt
    real    :: K1, K2, Kb, Kw, Ks1, Kf, Ksi, K1p, K2p, K3p
    real    :: ah1, ac, cu, cb, cc, satlev
    real    :: ratc13, ratc14, rato13, rato14, poso13, poso14
    integer, parameter :: niter = 5 ! number of iterations for carchm_solve

    ! Set array for saving diffusive sediment-water-column fluxes to zero
    !********************************************************************

    sedfluxo(:,:,:) = 0.0

    ! A LOOP OVER J
    ! RJ: This loop must go from 1 to kpje in the parallel version,
    !     otherways we had to do a boundary exchange

    !$OMP  PARALLEL DO                                                      &
    !$OMP& PRIVATE(sedb1,sediso,solrat,powcar,aerob,anaerob,                &
    !$OMP&         ex_dalk,ex_ddic,ex_disso_poc,                            &
    !$OMP&         dissot,undsa,posol,                                      &
    !$OMP&         umfa,denit,saln,rrho,alk,c,sit,pt,                       &
    !$OMP&         K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,                      &
    !$OMP&         ah1,ac,cu,cb,cc,satlev,                                  &
    !$OMP&         ratc13,ratc14,rato13,rato14,poso13,poso14,               &
    !$OMP&         k,i)

    j_loop: do j = 1, kpje

      do k = 1, ks
        do i = 1, kpie
          solrat(i,k) = 0.
          powcar(i,k) = 0.
          if  (use_extNcycle) then
            ex_ddic(i,k) = 0.
            ex_dalk(i,k) = 0.
          else
            anaerob(i,k) = 0.
          endif
          aerob(i,k)  = 0.
          sulf(i,k)   = 0.
          if (use_cisonew) then
            anaerob13(i,k)=0.
            aerob13(i,k)  =0.
            sulf13(i,k)   =0.
            anaerob14(i,k)=0.
            aerob14(i,k)  =0.
            sulf14(i,k)   =0.
          endif
        enddo
      enddo

      do k = 0, ks
        do i = 1, kpie
          sedb1(i,k) = 0.
          sediso(i,k) = 0.
        enddo
      enddo

      ! Calculate silicate-opal cycle and simultaneous silicate diffusion
      !******************************************************************

      ! Dissolution rate constant of opal (disso) [1/(kmol Si(OH)4/m3)*1/sec]*dtbgc
      dissot=disso_sil

      ! Evaluate boundary conditions for sediment-water column exchange.
      ! Current undersaturation of bottom water: sedb(i,0) and
      ! Approximation for new solid sediment, as from sedimentation flux: solrat(i,1)

      do i = 1, kpie
        if(omask(i,j) > 0.5) then
          undsa = silsat - powtra(i,j,1,ipowasi)
          sedb1(i,0) = bolay(i,j) * (silsat - ocetra(i,j,kbo(i,j),isilica))
          solrat(i,1) = ( sedlay(i,j,1,issssil)                                                    &
               + silpro(i,j) / (porsol(i,j,1) * seddw(1)) )                                        &
               * dissot / (1. + dissot * undsa) * porsol(i,j,1) / porwat(i,j,1)
        endif
      enddo

      ! Evaluate sediment undersaturation and degradation.
      ! Current undersaturation in pore water: sedb(i,k) and
      ! Approximation for new solid sediment, as from degradation: solrat(i,k)

      do k = 1, ks
        do i = 1, kpie
          if(omask(i,j) > 0.5) then
            undsa = silsat - powtra(i,j,k,ipowasi)
            sedb1(i,k) = seddw(k) * porwat(i,j,k) * (silsat - powtra(i,j,k,ipowasi))
            if ( k > 1 ) solrat(i,k) = sedlay(i,j,k,issssil)                                       &
                 * dissot / (1. + dissot * undsa) * porsol(i,j,k) / porwat(i,j,k)
          endif
        enddo
      enddo

      ! Solve for new undersaturation sediso, from current undersaturation sedb1,
      ! and first guess of new solid sediment solrat.

      call powadi(j,kpie,kpje,solrat,sedb1,sediso,omask)

      ! Update water column silicate, and store the flux for budget.
      ! Add sedimentation to first layer.

      do i = 1, kpie
        if(omask(i,j) > 0.5) then
          if(.not. lspin) then
            sedfluxo(i,j,ipowasi) =                                                                &
                 -(silsat - sediso(i,0) - ocetra(i,j,kbo(i,j),isilica))                            &
                 * bolay(i,j)
            ocetra(i,j,kbo(i,j),isilica) = silsat - sediso(i,0)
          endif
          sedlay(i,j,1,issssil) =                                                                  &
               sedlay(i,j,1,issssil) + silpro(i,j) / (porsol(i,j,1) * seddw(1))
        endif
      enddo


      ! Calculate updated degradation rate from updated undersaturation.
      ! Calculate new solid sediment.
      ! Update pore water concentration from new undersaturation.

      do k = 1, ks
        do i = 1, kpie
          if(omask(i,j) > 0.5) then
            umfa = porsol(i,j,k)/porwat(i,j,k)
            solrat(i,k) = sedlay(i,j,k,issssil) * dissot / (1. + dissot * sediso(i,k))
            posol = sediso(i,k) * solrat(i,k)
            sedlay(i,j,k,issssil) = sedlay(i,j,k,issssil) - posol
            powtra(i,j,k,ipowasi) = silsat - sediso(i,k)
          endif
        enddo
      enddo

      ! Calculate oxygen-POC cycle and simultaneous oxygen diffusion
      !*************************************************************

      ! Degradation rate constant of POP (disso) [1/(kmol O2/m3)*1/sec]*dtbgc
      dissot = disso_poc

      ! This scheme is not based on undersaturation, but on O2 itself

      ! Evaluate boundary conditions for sediment-water column exchange.
      ! Current concentration of bottom water: sedb(i,0) and
      ! Approximation for new solid sediment, as from sedimentation flux: solrat(i,1)

      do i = 1, kpie
        if(omask(i,j) > 0.5) then
          undsa = powtra(i,j,1,ipowaox)
          sedb1(i,0) = bolay(i,j) * ocetra(i,j,kbo(i,j),ioxygen)
          ex_disso_poc=merge(dissot*powtra(i,j,1,ipowaox)/(powtra(i,j,1,ipowaox)+bkox_drempoc_sed) & ! oxygen limitation
                         &       *POM_remin_q10_sed**((ptho(i,j,kbo(i,j))-POM_remin_Tref_sed)/10.) & ! T-dep
                         &   ,dissot,lTO2depremin)
          if ( .not.  use_extNcycle) then
            solrat(i,1) = ( sedlay(i,j,1,issso12) + prorca(i,j)                                    &
                        &  / (porsol(i,j,1) * seddw(1)) )                                          &
                        &   * ro2ut * ex_disso_poc / (1. + ex_disso_poc * undsa)                   &
                        &   * porsol(i,j,1) / porwat(i,j,1)
          else
            ! extended nitrogen cycle - 140mol O2/mol POP O2-consumption
            ! O2 and T-dep
            solrat(i,1) = ( sedlay(i,j,1,issso12) + prorca(i,j)                                    &
                        & / (porsol(i,j,1) * seddw(1)) )                                           &
                        & * ro2utammo * ex_disso_poc / (1. + ex_disso_poc * undsa)                 &
                        & * porsol(i,j,1) / porwat(i,j,1)
          endif
        endif
      enddo

      ! Evaluate sediment concentration and degradation.
      ! Current concentration in pore water: sedb(i,k) and
      ! Approximation for new solid sediment, as from degradation: solrat(i,k)

      do k = 1, ks
        do i = 1, kpie
          if(omask(i,j) > 0.5) then
            undsa = powtra(i,j,k,ipowaox)
            sedb1(i,k) = seddw(k) * porwat(i,j,k) * powtra(i,j,k,ipowaox)
            ex_disso_poc=merge(dissot*powtra(i,j,k,ipowaox)/(powtra(i,j,k,ipowaox)+bkox_drempoc_sed) & ! oxygen limitation
                         &       *POM_remin_q10_sed**((ptho(i,j,kbo(i,j))-POM_remin_Tref_sed)/10.) & ! T-dep
                         &   ,dissot,lTO2depremin)
            if ( .not. use_extNcycle) then
              if (k > 1) solrat(i,k) = sedlay(i,j,k,issso12) * ro2ut * ex_disso_poc                &
                                     & / (1. + ex_disso_poc*undsa) * porsol(i,j,k) / porwat(i,j,k)
            else
              ! extended nitrogen cycle - 140mol O2/mol POP O2-consumption
              if (k > 1) solrat(i,k) = sedlay(i,j,k,issso12) * ro2utammo * ex_disso_poc            &
                                     & /(1. + ex_disso_poc*undsa) * porsol(i,j,k) / porwat(i,j,k)
            endif
          endif
        enddo
      enddo

      ! Solve for new O2 concentration sediso, from current concentration sedb1,
      ! and first guess of new solid sediment solrat.

      call powadi(j,kpie,kpje,solrat,sedb1,sediso,omask)

      ! Update water column oxygen, and store the diffusive flux for budget (sedfluxo,
      ! positive downward). Add sedimentation to first layer.

      do i = 1, kpie
        if(omask(i,j) > 0.5) then
          if(.not. lspin) then
            sedfluxo(i,j,ipowaox) = -(sediso(i,0) - ocetra(i,j,kbo(i,j),ioxygen)) * bolay(i,j)
            ocetra(i,j,kbo(i,j),ioxygen) = sediso(i,0)
          endif
          sedlay(i,j,1,issso12) = sedlay(i,j,1,issso12) + prorca(i,j) / (porsol(i,j,1)*seddw(1))
          if (use_cisonew) then
            sedlay(i,j,1,issso13) = sedlay(i,j,1,issso13) + pror13(i,j) / (porsol(i,j,1)*seddw(1))
            sedlay(i,j,1,issso14) = sedlay(i,j,1,issso14) + pror14(i,j) / (porsol(i,j,1)*seddw(1))
          endif
        endif
      enddo


      ! Calculate updated degradation rate from updated concentration.
      ! Calculate new solid sediment.
      ! Update pore water concentration.
      ! Store flux in array aerob, for later computation of DIC and alkalinity.
      do k = 1, ks
        do i = 1, kpie
          if(omask(i,j) > 0.5) then
            umfa = porsol(i,j,k) / porwat(i,j,k)
            ex_disso_poc=merge(dissot*powtra(i,j,k,ipowaox)/(powtra(i,j,k,ipowaox)+bkox_drempoc_sed) & ! oxygen limitation
                         &       *POM_remin_q10_sed**((ptho(i,j,kbo(i,j))-POM_remin_Tref_sed)/10.) & ! T-dep
                         &   ,dissot,lTO2depremin)
            if (.not. use_extNcycle) then
              solrat(i,k) = sedlay(i,j,k,issso12) * ex_disso_poc/(1. + ex_disso_poc*sediso(i,k))
            else
              solrat(i,k) = sedlay(i,j,k,issso12) * ex_disso_poc/(1. + ex_disso_poc*sediso(i,k))
            endif
            posol = sediso(i,k)*solrat(i,k)
            if (use_cisonew) then
              rato13 = sedlay(i,j,k,issso13) / (sedlay(i,j,k,issso12) + safediv)
              rato14 = sedlay(i,j,k,issso14) / (sedlay(i,j,k,issso12) + safediv)
              poso13 = posol*rato13
              poso14 = posol*rato14
              aerob13(i,k) = poso13*umfa  !this has P units: kmol P/m3 of pore water
              aerob14(i,k) = poso14*umfa  !this has P units: kmol P/m3 of pore water
            endif
            sedlay(i,j,k,issso12) = sedlay(i,j,k,issso12) - posol
            powtra(i,j,k,ipowaph) = powtra(i,j,k,ipowaph) + posol*umfa
            if (.not. use_extNcycle) then
              powtra(i,j,k,ipowno3) = powtra(i,j,k,ipowno3) + posol*rnit*umfa
              aerob(i,k) = posol*umfa     !this has P units: kmol P/m3 of pore water
            else
              powtra(i,j,k,ipownh4) = powtra(i,j,k,ipownh4) + posol*rnit*umfa
              ex_ddic(i,k) = rcar*posol*umfa ! C-units kmol C/m3 of pore water
              ex_dalk(i,k) = (rnit-1.)*posol*umfa ! alkalinity units
              extNsed_diagnostics(i,j,k,ised_remin_aerob) = posol*rnit*umfa ! Output
            endif
            powtra(i,j,k,ipowaox) = sediso(i,k)
            if (use_cisonew) then
              sedlay(i,j,k,issso13) = sedlay(i,j,k,issso13) - poso13
              sedlay(i,j,k,issso14) = sedlay(i,j,k,issso14) - poso14
            endif
          endif
        enddo
      enddo

      ! Calculate nitrate reduction under anaerobic conditions explicitely
      !*******************************************************************
      if (.not. use_extNcycle) then
        ! Denitrification rate constant of POP (disso) [1/sec]*dtbgc
        denit = sed_denit

        ! Store flux in array anaerob, for later computation of DIC and alkalinity.
        do k = 1, ks
          do i = 1, kpie
            if(omask(i,j) > 0.5) then
              if(powtra(i,j,k,ipowaox) < 1.e-6) then
                posol = denit * min(0.25*powtra(i,j,k,ipowno3)/rdnit2, sedlay(i,j,k,issso12))
                umfa = porsol(i,j,k)/porwat(i,j,k)
                anaerob(i,k) = posol*umfa     !this has P units: kmol P/m3 of pore water
                if (use_cisonew) then
                  rato13 = sedlay(i,j,k,issso13) / (sedlay(i,j,k,issso12) + safediv)
                  rato14 = sedlay(i,j,k,issso14) / (sedlay(i,j,k,issso12) + safediv)
                  poso13 = posol * rato13
                  poso14 = posol * rato14
                  anaerob13(i,k) = poso13*umfa  !this has P units: kmol P/m3 of pore water
                  anaerob14(i,k) = poso14*umfa  !this has P units: kmol P/m3 of pore water
                endif
                sedlay(i,j,k,issso12) = sedlay(i,j,k,issso12) - posol
                powtra(i,j,k,ipowaph) = powtra(i,j,k,ipowaph) + posol*umfa
                powtra(i,j,k,ipowno3) = powtra(i,j,k,ipowno3) - rdnit1*posol*umfa
                powtra(i,j,k,ipown2)  = powtra(i,j,k,ipown2)  + rdnit2*posol*umfa
                if (use_cisonew) then
                  sedlay(i,j,k,issso13) = sedlay(i,j,k,issso13) - poso13
                  sedlay(i,j,k,issso14) = sedlay(i,j,k,issso14) - poso14
                endif
              endif
            endif
          enddo
        enddo
      else
        !======>>>> extended nitrogen cycle processes (aerobic and anaerobic) that follow ammonification
        call sed_nitrification(j,kpie,kpje,kpke,kbnd,ptho,omask,ex_ddic,ex_dalk)
        call sed_denit_NO3_to_NO2(j,kpie,kpje,kpke,kbnd,ptho,omask,ex_ddic,ex_dalk)
        call sed_anammox(j,kpie,kpje,kpke,kbnd,ptho,omask,ex_ddic,ex_dalk)
        call sed_denit_dnra(j,kpie,kpje,kpke,kbnd,ptho,omask,ex_ddic,ex_dalk)
      endif

      ! sulphate reduction in sediments
      do k = 1, ks
        do i = 1, kpie
          if(omask(i,j) > 0.5) then
            if(powtra(i,j,k,ipowaox) < 3.e-6 .and. powtra(i,j,k,ipowno3) < 3.e-6) then
              posol = denit * sedlay(i,j,k,issso12)         ! remineralization of poc
              umfa = porsol(i,j,k) / porwat(i,j,k)
              sulf(i,k) = posol*umfa      !this has P units: kmol P/m3 of pore water
              if (use_cisonew) then
                rato13 = sedlay(i,j,k,issso13) / (sedlay(i,j,k,issso12)+safediv)
                rato14 = sedlay(i,j,k,issso14) / (sedlay(i,j,k,issso12)+safediv)
                poso13 = posol * rato13
                poso14 = posol * rato14
                sulf13(i,k) = poso13*umfa !this has P units: kmol P/m3 of pore water
                sulf14(i,k) = poso14*umfa !this has P units: kmol P/m3 of pore water
              endif
              sedlay(i,j,k,issso12) = sedlay(i,j,k,issso12) - posol
              powtra(i,j,k,ipowaph) = powtra(i,j,k,ipowaph) + posol*umfa
              powtra(i,j,k,ipowno3) = powtra(i,j,k,ipowno3) + posol*umfa*rnit
              if (use_cisonew) then
                sedlay(i,j,k,issso13) = sedlay(i,j,k,issso13) - poso13
                sedlay(i,j,k,issso14) = sedlay(i,j,k,issso14) - poso14
              endif
              if (use_extNcycle) then
                extNsed_diagnostics(i,j,k,ised_remin_sulf) = posol*umfa ! Output
              endif
            endif
          endif
        enddo
      enddo   ! end sulphate reduction


      ! Calculate CaCO3-CO3 cycle and simultaneous CO3-undersaturation diffusion
      !*************************************************************************

      ! Compute new powcar, carbonate ion concentration in the sediment
      ! from changed alkalinity (nitrate production during remineralisation)
      ! and DIC gain. Iterate 5 times. This changes pH (sedhpl) of sediment.

      do k = 1, ks
        do i = 1, kpie
          if(omask(i,j) > 0.5) then
            saln= min( 40., max( 0., psao(i,j,kbo(i,j))))
            rrho= prho(i,j,kbo(i,j))
            if (use_extNcycle) then
              alk = (powtra(i,j,k,ipowaal) - (sulf(i,k)+aerob(i,k))*(rnit+1.) + ex_dalk(i,k))  / rrho
              c   = (powtra(i,j,k,ipowaic) + (aerob(i,k)+sulf(i,k))*rcar + ex_ddic(i,k)) / rrho
            else
              alk = (powtra(i,j,k,ipowaal) - (sulf(i,k)+aerob(i,k))*(rnit+1.) + anaerob(i,k)*(rdnit1-1.))  / rrho
              c   = (powtra(i,j,k,ipowaic) + (anaerob(i,k)+aerob(i,k)+sulf(i,k))*rcar) / rrho
            endif
            sit =  powtra(i,j,k,ipowasi) / rrho
            pt  =  powtra(i,j,k,ipowaph) / rrho
            ah1 = sedhpl(i,j,k)
            K1  = keqb( 1,i,j)
            K2  = keqb( 2,i,j)
            Kb  = keqb( 3,i,j)
            Kw  = keqb( 4,i,j)
            Ks1 = keqb( 5,i,j)
            Kf  = keqb( 6,i,j)
            Ksi = keqb( 7,i,j)
            K1p = keqb( 8,i,j)
            K2p = keqb( 9,i,j)
            K3p = keqb(10,i,j)

            call carchm_solve(saln,c,alk,sit,pt,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,ah1,ac,niter)

            cu = ( 2. * c - ac ) / ( 2. + K1 / ah1 )
            cb = K1 * cu / ah1
            cc = K2 * cb / ah1
            sedhpl(i,j,k) = max( 1.e-20, ah1 )
            powcar(i,k)   = cc * rrho
          endif
        enddo
      enddo

      ! Dissolution rate constant of CaCO3 (disso) [1/(kmol CO3--/m3)*1/sec]*dtbgc
      dissot = disso_caco3

      ! Evaluate boundary conditions for sediment-water column exchange.
      ! Current undersaturation of bottom water: sedb(i,0) and
      ! Approximation for new solid sediment, as from sedimentation flux: solrat(i,1)

      ! CO3 saturation concentration is aksp/calcon as in CARCHM
      ! (calcon defined in MO_CHEMCON with 1.028e-2; 1/calcon =~ 97.)

      do i = 1, kpie
        if(omask(i,j) > 0.5) then
          satlev = keqb(11,i,j) / calcon + 2.e-5
          undsa = max( satlev-powcar(i,1), 0. )
          sedb1(i,0) = bolay(i,j) * (satlev-co3(i,j,kbo(i,j)))
          solrat(i,1) = (sedlay(i,j,1,isssc12)                                                     &
               &   + prcaca(i,j) / (porsol(i,j,1)*seddw(1)))                                       &
               &   * dissot / (1.+dissot*undsa) * porsol(i,j,1) / porwat(i,j,1)
        endif
      enddo

      ! Evaluate sediment undersaturation and degradation.
      ! Current undersaturation in pore water: sedb(i,k) and
      ! Approximation for new solid sediment, as from degradation: solrat(i,k)

      do k = 1, ks
        do i = 1, kpie
          if(omask(i,j) > 0.5) then
            undsa = max( keqb(11,i,j) / calcon - powcar(i,k), 0. )
            sedb1(i,k) = seddw(k) * porwat(i,j,k) * undsa
            if (k > 1) then
              solrat(i,k) = sedlay(i,j,k,isssc12) * dissot/(1.+dissot*undsa) * porsol(i,j,k)/porwat(i,j,k)
            end if
            if (undsa <= 0.) then
              solrat(i,k) = 0.
            end if
          endif
        enddo
      enddo

      ! Solve for new undersaturation sediso, from current undersaturation sedb1,
      ! and first guess of new solid sediment solrat.

      call powadi(j,kpie,kpje,solrat,sedb1,sediso,omask)

      ! There is no exchange between water and sediment with respect to co3 so far.
      ! Add sedimentation to first layer.
      do i = 1, kpie
        if(omask(i,j) > 0.5) then
          sedlay(i,j,1,isssc12) =                                                                  &
                 &   sedlay(i,j,1,isssc12) + prcaca(i,j) / (porsol(i,j,1)*seddw(1))
          if (use_cisonew) then
            sedlay(i,j,1,isssc13) =                                                                &
                 &   sedlay(i,j,1,isssc13) + prca13(i,j) / (porsol(i,j,1)*seddw(1))
            sedlay(i,j,1,isssc14) =                                                                &
                 &   sedlay(i,j,1,isssc14) + prca14(i,j) / (porsol(i,j,1)*seddw(1))
          endif
        endif
      enddo

      ! Calculate updated degradation rate from updated undersaturation.
      ! Calculate new solid sediment.
      ! No update of powcar pore water concentration from new undersaturation so far.
      ! Instead, only update DIC, and, of course, alkalinity.
      ! This also includes gains from aerobic and anaerobic decomposition.

      do k = 1, ks
        do i = 1, kpie
          if(omask(i,j) > 0.5) then
            umfa = porsol(i,j,k) / porwat(i,j,k)
            solrat(i,k) = sedlay(i,j,k,isssc12) * dissot / (1. + dissot * sediso(i,k))
            posol = sediso(i,k) * solrat(i,k)
            if (use_cisonew) then
              ratc13 = sedlay(i,j,k,isssc13) / (sedlay(i,j,k,isssc12) + safediv)
              ratc14 = sedlay(i,j,k,isssc14) / (sedlay(i,j,k,isssc12) + safediv)
              poso13 = posol * ratc13
              poso14 = posol * ratc14
            endif
            sedlay(i,j,k,isssc12) = sedlay(i,j,k,isssc12) - posol
            if (use_extNcycle) then
              powtra(i,j,k,ipowaic) = powtra(i,j,k,ipowaic)                                        &
                  &   + posol * umfa + (aerob(i,k) + sulf(i,k)) * rcar + ex_ddic(i,k)
              powtra(i,j,k,ipowaal) = powtra(i,j,k,ipowaal)                                        &
                  &   + 2. * posol * umfa - (rnit+1.)*(aerob(i,k) + sulf(i,k))  + ex_dalk(i,k)
            else
              powtra(i,j,k,ipowaic) = powtra(i,j,k,ipowaic)                                        &
                  &   + posol * umfa + (aerob(i,k) + anaerob(i,k) + sulf(i,k)) * rcar
              powtra(i,j,k,ipowaal) = powtra(i,j,k,ipowaal)                                        &
                  &   + 2. * posol * umfa - (rnit+1.)*(aerob(i,k) + sulf(i,k)) + (rdnit1-1.)*anaerob(i,k)
            endif
            if (use_cisonew) then
              sedlay(i,j,k,isssc13) = sedlay(i,j,k,isssc13) - poso13
              sedlay(i,j,k,isssc14) = sedlay(i,j,k,isssc14) - poso14
              powtra(i,j,k,ipowc13) = powtra(i,j,k,ipowc13) + poso13 * umfa                        &
                 &   + (aerob13(i,k) + anaerob13(i,k) + sulf13(i,k)) * rcar
              powtra(i,j,k,ipowc14) = powtra(i,j,k,ipowc14) + poso14 * umfa                        &
                 &   + (aerob14(i,k) + anaerob14(i,k) + sulf14(i,k)) * rcar
            endif
          endif
        enddo
      enddo

    enddo j_loop

    !$OMP END PARALLEL DO

    call dipowa(kpie,kpje,kpke,omask,lspin)

    !ik add clay sedimentation onto sediment
    !ik this is currently assumed to depend on total and corg sedimentation:
    !ik f(POC) [kg C] / f(total) [kg] = 0.05
    !ik thus it is
    !$OMP PARALLEL DO PRIVATE(i)
    do j = 1, kpje
      do i = 1, kpie
        sedlay(i,j,1,issster) = sedlay(i,j,1,issster) + produs(i,j) / (porsol(i,j,1) * seddw(1))
      enddo
    enddo
    !$OMP END PARALLEL DO

    if(.not. lspin) then
      !$OMP PARALLEL DO PRIVATE(i)
      do j = 1, kpje
        do i = 1, kpie
          silpro(i,j) = 0.
          prorca(i,j) = 0.
          prcaca(i,j) = 0.
          if (use_cisonew) then
            pror13(i,j) = 0.
            pror14(i,j) = 0.
            prca13(i,j) = 0.
            prca14(i,j) = 0.
          endif
          produs(i,j) = 0.
        enddo
      enddo
      !$OMP END PARALLEL DO
    endif

  end subroutine powach

end module mo_powach
