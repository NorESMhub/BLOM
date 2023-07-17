! Copyright (C) 2001  Ernst Maier-Reimer
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


      SUBROUTINE HAMOCC4BCM(kpie,kpje,kpke,kbnd,kplyear,kplmon,kplday,kldtday,&
                            pdlxp,pdlyp,pddpo,prho,pglat,omask,               &
                            dust,rivin,ndep,oafx,pi_ph,                       &
                            pfswr,psicomo,ppao,pfu10,ptho,psao,               &
                            patmco2,pflxco2,pflxdms,patmbromo,pflxbromo,      &
                            patmn2o,pflxn2o,patmnh3,pflxnh3,patmnhxdep,patmnoydep)
!******************************************************************************
!
! HAMOCC4BGC - main routine of iHAMOCC.
!
! Modified
! --------
!  J.Schwinger       *GFI, Bergen*    2013-10-21
!  - added GNEWS2 option for riverine input of carbon and nutrients
!  - code cleanup
!
!  J.Schwinger       *GFI, Bergen*    2014-05-21
!  - moved copying of tracer field to ocetra to micom2hamocc 
!    and hamocc2micom
!
!  J.Schwinger,      *Uni Research, Bergen*   2018-04-12
!  - moved accumulation of all output fields to seperate subroutine,
!    related code-restructuring
!  - added sediment bypass preprocessor option
!
!  J.Schwinger,      *NORCE Climate, Bergen*   2020-05-28
!  - restructuring of iHAMOCC code, cleanup parameter list
!  - boundary conditions (dust, riverinput, N-deposition) are now passed as 
!    an argument
!
! Parameter list:
! ---------------
!
!  *INTEGER* *kpie*       - 1st dimension of model grid.
!  *INTEGER* *kpje*       - 2nd dimension of model grid.
!  *INTEGER* *kpke*       - 3rd (vertical) dimension of model grid.
!  *INTEGER* *kbnd*       - nb of halo grid points.
!  *INTEGER* *kplyear*    - current year.
!  *INTEGER* *kplmon*     - current month.
!  *INTEGER* *kplday*     - current day.
!  *INTEGER* *kldtday*    - number of time step in current day.
!  *REAL*    *pdlxp*      - size of grid cell (longitudinal) [m].
!  *REAL*    *pdlyp*      - size of grid cell (latitudinal) [m].
!  *REAL*    *pddpo*      - size of grid cell (depth) [m].
!  *REAL*    *prho*       - density [kg/m^3].
!  *REAL*    *pglat*      - latitude of grid cells [deg north].
!  *REAL*    *omask*      - land/ocean mask.
!  *REAL*    *dust*       - dust deposition flux [kg/m2/month].
!  *REAL*    *rivin*      - riverine input [kmol m-2 yr-1].
!  *REAL*    *ndep*       - nitrogen deposition [kmol m-2 yr-1].
!  *REAL*    *oaflx*      - alkalinity flux from alkalinization [kmol m-2 yr-1]
!  *REAL*    *pfswr*      - solar radiation [W/m**2].
!  *REAL*    *psicomo*    - sea ice concentration
!  *REAL*    *ppao*       - sea level pressure [Pascal].
!  *REAL*    *pfu10*      - absolute wind speed at 10m height [m/s]
!  *REAL*    *ptho*       - potential temperature [deg C].
!  *REAL*    *psao*       - salinity [psu.].
!  *REAL*    *patmco2*    - atmospheric CO2 concentration [ppm] used in 
!                           fully coupled mode (prognostic/diagnostic CO2).
!  *REAL*    *pflxco2*    - CO2 flux [kg/m^2/s].
!  *REAL*    *pflxdms*    - DMS flux [kg/m^2/s].
!  *REAL*    *patmbromo*  - atmospheric bromoform concentration [ppt] used in 
!                           fully coupled mode.
!  *REAL*    *pflxbromo*  - Bromoform flux [kg/m^2/s].
!  *REAL*    *patmn2o*    - atmospheric nitrous oxide concentration [ppt] used in 
!                           fully coupled mode.
!  *REAL*    *pflxn2o*    - Nitrous oxide flux [kg N2O /m^2/s].
!  *REAL*    *patmnh3*    - atmospheric ammonia concentration [ppt] used in 
!                           fully coupled mode.
!  *REAL*    *pflxnh3*    - Ammonia flux [kg NH3 /m^2/s].
!  *REAL*    *patmnhxdep* - Atmospheric NHx deposition kgN/m2/s
!  *REAL*    *patmnoydep* - Atmospheric NOy deposition kgN/m2/s
!
!******************************************************************************
      use mod_xc,         only: mnproc
      use mo_carbch,      only: atmflx,ocetra,atm
      use mo_biomod,      only: strahl
      use mo_control_bgc, only: ldtrunbgc,dtbgc,ldtbgc,io_stdo_bgc,dtbgc,ndtdaybgc, &
                                do_sedspinup,sedspin_yr_s,sedspin_yr_e,sedspin_ncyc,&
                                do_ndep_coupled
      use mo_param1_bgc,  only: iatmco2,iatmdms,nocetra,nriv
      use mo_vgrid,       only: set_vgrid
      use mo_apply_fedep, only: apply_fedep
      use mo_apply_rivin, only: apply_rivin
      use mo_apply_ndep,  only: apply_ndep
      use mo_apply_oafx,  only: apply_oafx
#if defined(BOXATM)
      use mo_boxatm,      only: update_boxatm
#endif
#ifdef BROMO
      use mo_param1_bgc,  only: iatmbromo
#endif
#ifdef CFC
      use mo_carbch,      only: atm_cfc11_nh,atm_cfc11_sh,atm_cfc12_nh,atm_cfc12_sh,atm_sf6_nh,atm_sf6_sh
#endif
#ifdef extNcycle
      use mo_param1_bgc,  only: iatmn2o,iatmnh3
#endif
      implicit none

      INTEGER, intent(in)  :: kpie,kpje,kpke,kbnd
      INTEGER, intent(in)  :: kplyear,kplmon,kplday,kldtday
      REAL,    intent(in)  :: pdlxp  (kpie,kpje)
      REAL,    intent(in)  :: pdlyp  (kpie,kpje)
      REAL,    intent(in)  :: pddpo  (kpie,kpje,kpke)
      REAL,    intent(in)  :: prho   (kpie,kpje,kpke)
      REAL,    intent(in)  :: pglat  (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in)  :: omask  (kpie,kpje)
      REAL,    intent(in)  :: dust   (kpie,kpje)
      REAL,    intent(in)  :: rivin  (kpie,kpje,nriv)
      REAL,    intent(inout):: ndep   (kpie,kpje,2)
      REAL,    intent(in)  :: oafx   (kpie,kpje)
      REAL,    intent(in)  :: pi_ph  (kpie,kpje)
      REAL,    intent(in)  :: pfswr  (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in)  :: psicomo(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in)  :: ppao   (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in)  :: pfu10  (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in)  :: ptho   (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)
      REAL,    intent(in)  :: psao   (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)
      REAL,    intent(in)  :: patmco2(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(out) :: pflxco2(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(out) :: pflxdms(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in)  :: patmbromo(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(out) :: pflxbromo(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in)  :: patmn2o(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(out) :: pflxn2o(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in)  :: patmnh3(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(out) :: pflxnh3(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in)  :: patmnhxdep(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in)  :: patmnoydep(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)

      INTEGER :: i,j,k,l
      INTEGER :: nspin,it
      LOGICAL :: lspin
      REAL    :: fatmndep

      IF (mnproc.eq.1) THEN
      write(io_stdo_bgc,*) 'iHAMOCC',KLDTDAY,LDTRUNBGC,NDTDAYBGC
      ENDIF


!--------------------------------------------------------------------
! Increment bgc time step counter of run (initialized in HAMOCC_INIT).
!
      ldtrunbgc = ldtrunbgc + 1


!--------------------------------------------------------------------
! Increment bgc time step counter of experiment.
!
      ldtbgc = ldtbgc + 1


!--------------------------------------------------------------------
! Calculate variables related to the vertical grid
!
      call set_vgrid(kpie,kpje,kpke,pddpo)


!--------------------------------------------------------------------
! Pass net solar radiation
!
!$OMP PARALLEL DO PRIVATE(i)
      DO  j=1,kpje
      DO  i=1,kpie
        strahl(i,j)=pfswr(i,j)
      ENDDO
      ENDDO
!$OMP END PARALLEL DO


!--------------------------------------------------------------------
! Pass atmospheric co2 if coupled to an active atmosphere model
!
#if defined(PROGCO2) || defined(DIAGCO2)
!$OMP PARALLEL DO PRIVATE(i)
      DO  j=1,kpje
      DO  i=1,kpie
        atm(i,j,iatmco2)=patmco2(i,j)
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      !if (mnproc.eq.1) write (io_stdo_bgc,*) 'iHAMOCC: getting co2 from atm'
#endif

#ifdef BROMO
!$OMP PARALLEL DO PRIVATE(i)
      DO  j=1,kpje
      DO  i=1,kpie
        IF (patmbromo(i,j).gt.0.) THEN
         atm(i,j,iatmbromo)=patmbromo(i,j)
        ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      if (mnproc.eq.1) write (io_stdo_bgc,*) 'iHAMOCC: getting bromoform from atm'
#endif

#ifdef extNcycle
!$OMP PARALLEL DO PRIVATE(i)
      DO  j=1,kpje
      DO  i=1,kpie
        IF (patmn2o(i,j).gt.0.) THEN
         atm(i,j,iatmn2o)=patmn2o(i,j)
        ENDIF
        IF (patmnh3(i,j).gt.0.) THEN
         atm(i,j,iatmnh3)=patmnh3(i,j)
        ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      if (mnproc.eq.1) write (io_stdo_bgc,*) 'iHAMOCC: getting N2O and NH3 conc. from atm'
      
      IF(do_ndep_coupled) THEN
        fatmndep = 365.*86400./14.00674
        ndep(:,:,:) = 0.
!$OMP PARALLEL DO PRIVATE(i)
        DO  j=1,kpje
        DO  i=1,kpie
          ! convert from kgN/m2/s to climatological input file units: kmolN/m2/yr 
          IF (patmnoydep(i,j).gt.0.) THEN
            ndep(i,j,1) = patmnoydep(i,j)*fatmndep
          ENDIF
          IF (patmnhxdep(i,j).gt.0.) THEN
            ndep(i,j,2) = patmnhxdep(i,j)*fatmndep
          ENDIF
        ENDDO
        ENDDO
!$OMP END PARALLEL DO
      if (mnproc.eq.1) write (io_stdo_bgc,*) 'iHAMOCC: getting NOy and NHx deposition from atm'
      ENDIF
#endif

!--------------------------------------------------------------------
! Read atmospheric cfc concentrations
!
#ifdef CFC
      call get_cfc(kplyear,atm_cfc11_nh,atm_cfc12_nh,atm_sf6_nh,        &
                           atm_cfc11_sh,atm_cfc12_sh,atm_sf6_sh)
#endif

#ifdef PBGC_CK_TIMESTEP
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'before BGC: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif


!---------------------------------------------------------------------
! Biogeochemistry
!
      ! Apply dust (iron) deposition
      ! This routine should be moved to the other routines that handle 
      ! external inputs below for consistency. For now we keep it here
      ! to maintain bit-for-bit reproducibility with the CMIP6 version of 
      ! the model
      CALL apply_fedep(kpie,kpje,kpke,pddpo,omask,dust)

      CALL OCPROD(kpie,kpje,kpke,kbnd,pdlxp,pdlyp,pddpo,omask,ptho,pi_ph,psao,ppao,prho)

#ifdef PBGC_CK_TIMESTEP   
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after OCPROD: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif

 
      do l=1,nocetra
      do K=1,kpke
!$OMP PARALLEL DO PRIVATE(i)
      do J=1,kpje
      do I=1,kpie
        if (OMASK(I,J) .gt. 0.5 ) then
          OCETRA(I,J,K,L)=MAX(0.,OCETRA(I,J,K,L))
        endif
      enddo
      enddo
!$OMP END PARALLEL DO
      enddo
      enddo

#ifdef PBGC_CK_TIMESTEP   
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after LIMIT: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif


      CALL CYANO(kpie,kpje,kpke,kbnd,pddpo,omask,ptho)

#ifdef PBGC_CK_TIMESTEP   
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after CYANO: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif


      CALL CARCHM(kpie,kpje,kpke,kbnd,pdlxp,pdlyp,pddpo,prho,pglat,omask,      &
                  psicomo,ppao,pfu10,ptho,psao)

#ifdef PBGC_CK_TIMESTEP   
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after CARCHM: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif


      ! Apply n-deposition
      CALL apply_ndep(kpie,kpje,kpke,pddpo,omask,ndep)

#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after N deposition: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif	 

      ! Apply riverine input of carbon and nutrients
      call apply_rivin(kpie,kpje,kpke,pddpo,omask,rivin)

#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after river input: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif	 

      ! Apply alkalinity flux due to ocean alkalinization
      call apply_oafx(kpie,kpje,kpke,pddpo,omask,oafx)

#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after ocean alkalinization: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif	 

      ! Apply alkalinity flux due to ocean alkalinization
      call apply_oafx(kpie,kpje,kpke,pddpo,omask,oafx)

#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after ocean alkalinization: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif	 

      ! Update atmospheric pCO2 [ppm]
#if defined(BOXATM)
      CALL update_boxatm(kpie,kpje,pdlxp,pdlyp)
#endif	 

#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after ATMOTR: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif	 

      ! update preformed tracers
      CALL PREFTRC(kpie,kpje,omask)


!--------------------------------------------------------------------
!     Sediment module

#ifndef sedbypass
! jump over sediment if sedbypass is defined

      if(do_sedspinup .and. kplyear>=sedspin_yr_s                              &
                      .and. kplyear<=sedspin_yr_e) then
        nspin = sedspin_ncyc
        if(mnproc == 1) then
          write(io_stdo_bgc,*)
          write(io_stdo_bgc,*) 'iHAMOCC: sediment spinup activated with ',     &
                                nspin, ' subcycles' 
        endif
      else
        nspin = 1
      endif
      
      ! Loop for sediment spinup. If deactivated then nspin=1 and lspin=.false.
      do it=1,nspin

        if( it<nspin ) then
          lspin=.true.
        else
          lspin=.false.      
        endif

        call POWACH(kpie,kpje,kpke,kbnd,prho,omask,psao,ptho,lspin)

      enddo

#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after POWACH: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif	 

!     sediment is shifted once a day (on both time levels!)
      IF(KLDTDAY .EQ. 1 .OR. KLDTDAY .EQ. 2) THEN
         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)                                           &
     &   'Sediment shifting ...'
         ENDIF

         CALL SEDSHI(kpie,kpje,omask)

      ENDIF
#endif


#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after BGC: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif	 

!--------------------------------------------------------------------
! Pass co2 flux. Convert unit from kmol/m^2 to kg/m^2/s.

!$OMP PARALLEL DO PRIVATE(i)
      DO  j=1,kpje
      DO  i=1,kpie
        if(omask(i,j) .gt. 0.5) pflxco2(i,j)=-44.*atmflx(i,j,iatmco2)/dtbgc
      ENDDO
      ENDDO
!$OMP END PARALLEL DO


!--------------------------------------------------------------------
! Pass dms flux. Convert unit from kmol/m^2 to kg/m^2/s.

!$OMP PARALLEL DO PRIVATE(i)
      DO  j=1,kpje
      DO  i=1,kpie
        if(omask(i,j) .gt. 0.5) pflxdms(i,j)=-62.13*atmflx(i,j,iatmdms)/dtbgc
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

!--------------------------------------------------------------------
! Pass bromoform flux. Convert unit from kmol CHBr3/m^2 to kg/m^2/s.
! Negative values to the atmosphere
!$OMP PARALLEL DO PRIVATE(i)
      DO  j=1,kpje
      DO  i=1,kpie
#ifdef BROMO
        if(omask(i,j) .gt. 0.5) pflxbromo(i,j)=-252.7*atmflx(i,j,iatmbromo)/dtbgc
#else
        if(omask(i,j) .gt. 0.5) pflxbromo(i,j)=0.0
#endif
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
!--------------------------------------------------------------------
! Pass nitrous oxide and ammonia fluxes. Convert unit from kmol N2O (NH3)/m2/Delta t to kg/m2/s
! negative values to the atmosphere 
!$OMP PARALLEL DO PRIVATE(i)
      DO  j=1,kpje
      DO  i=1,kpie
#ifdef extNcycle
        if(omask(i,j) .gt. 0.5) pflxn2o(i,j)=-44.012880*atmflx(i,j,iatmn2o)/dtbgc  ! conversion factor checked against CAM 
        if(omask(i,j) .gt. 0.5) pflxnh3(i,j)=-17.028940*atmflx(i,j,iatmnh3)/dtbgc  ! conversion factor checked against CAM
#else
        if(omask(i,j) .gt. 0.5) pflxn2o(i,j)=0.0
        if(omask(i,j) .gt. 0.5) pflxnh3(i,j)=0.0
#endif
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
!--------------------------------------------------------------------
      RETURN
      END
