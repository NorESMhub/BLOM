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
                            dust,rivin,ndep,                                  &
                            pfswr,psicomo,ppao,pfu10,ptho,psao,               &
                            patmco2,pflxco2,pflxdms)
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
!  *REAL*    *rivin*      - riverine input [kmol m-2 yr-2].
!  *REAL*    *ndep*       - nitrogen deposition [kmol m-2 yr-2].
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
!
!******************************************************************************
      use mod_xc
      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      USE mo_bgcmean
      USE mo_control_bgc
      use mo_param1_bgc
      use mo_vgrid,     only: set_vgrid
      use mo_riverinpt, only: riverinpt,nriv
      use mo_ndep,      only: n_deposition
#if defined(BOXATM)
      use mo_boxatm
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
      REAL,    intent(in)  :: dust   (kpie,kpje,nriv)
      REAL,    intent(in)  :: rivin  (kpie,kpje,nriv)
      REAL,    intent(in)  :: ndep   (kpie,kpje)
      REAL,    intent(in)  :: pfswr  (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in)  :: psicomo(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in)  :: ppao   (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in)  :: pfu10  (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in)  :: ptho   (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)
      REAL,    intent(in)  :: psao   (1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)
      REAL,    intent(in)  :: patmco2(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(out) :: pflxco2(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(out) :: pflxdms(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)

      INTEGER :: i,j,k,l

      IF (mnproc.eq.1) THEN
      write(io_stdo_bgc,*) 'iHAMOCC',KLDTDAY,LDTRUNBGC,NDTDAYBGC
      ENDIF


!--------------------------------------------------------------------
! Increment bgc time step counter of run (initialized in INI_BGC).
!
      ldtrunbgc = ldtrunbgc + 1


!--------------------------------------------------------------------
! Increment bgc time step counter of experiment (initialized if IAUFR=0).
!
      ldtbgc = ldtbgc + 1


!--------------------------------------------------------------------
! Calculate variables related to the vertical grid
!
      call set_vgrid(kpie,kpje,kpke,pddpo)


!--------------------------------------------------------------------
! Pass net solar radiation 
!
!$OMP PARALLEL DO
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
!$OMP PARALLEL DO
      DO  j=1,kpje
      DO  i=1,kpie
        atm(i,j,iatmco2)=patmco2(i,j)
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      !if (mnproc.eq.1) write (io_stdo_bgc,*) 'iHAMOCC: getting co2 from atm'
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
      CALL OCPROD(kpie,kpje,kpke,kbnd,pdlxp,pdlyp,pddpo,omask,dust,ptho)

#ifdef PBGC_CK_TIMESTEP   
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after OCPROD: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif

 
      do l=1,nocetra
      do K=1,kpke
!$OMP PARALLEL DO
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
      CALL n_deposition(kpie,kpje,kpke,pddpo,omask,ndep)

      ! Apply riverine input of carbon and nutrients
      call riverinpt(kpie,kpje,kpke,pddpo,omask,rivin)

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

      CALL POWACH(kpie,kpje,kpke,kbnd,prho,omask,psao)


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

!$OMP PARALLEL DO
      DO  j=1,kpje
      DO  i=1,kpie
        if(omask(i,j) .gt. 0.5) pflxco2(i,j)=-44.*atmflx(i,j,iatmco2)/dtbgc
      ENDDO
      ENDDO
!$OMP END PARALLEL DO


!--------------------------------------------------------------------
! Pass dms flux. Convert unit from kmol/m^2 to kg/m^2/s.

!$OMP PARALLEL DO
      DO  j=1,kpje
      DO  i=1,kpie
        if(omask(i,j) .gt. 0.5) pflxdms(i,j)=-62.13*atmflx(i,j,iatmdms)/dtbgc
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

!--------------------------------------------------------------------
      RETURN
      END
