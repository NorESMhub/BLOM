      SUBROUTINE HAMOCC4BCM(kpie,kpje,kpke,pglat,                        &
     &    pfswr,psicomo,ptho,psao,ppao,prho,pddpo,pdlxp,pdlyp,ptiestu,   &
     &    ptiestw,pfu10,patmco2,pflxco2,kplyear,kplmon,kplday,           &
     &    kmonlen,kldtmon,kldtday,omask,days_in_yr,pflxdms)       

!**********************************************************************
!
!**** *BGC* - .
!
!     Modified
!     --------
!     J.Schwinger       *GFI, Bergen*    2013-10-21
!     - added GNEWS2 option for riverine input of carbon and nutrients
!     - code cleanup
!
!     J.Schwinger       *GFI, Bergen*    2014-05-21
!     - moved copying of tracer field to ocetra to micom2hamocc 
!       and hamocc2micom
!
!     J.Schwinger,      *Uni Research, Bergen*   2018-04-12
!     - moved accumulation of all output fields to seperate subroutine,
!       related code-restructuring
!     - added sediment bypass preprocessor option
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*       - 1st dimension of model grid.
!     *INTEGER* *kpje*       - 2nd dimension of model grid.
!     *INTEGER* *kpke*       - 3rd (vertical) dimension of model grid.
!     *REAL*    *pglat*      - latitude og grid cells [deg north].
!     *REAL*    *pfswr*      - solar radiation [W/m**2].
!     *REAL*    *psicomo*    - sea ice concentration
!     *REAL*    *ptho*       - potential temperature [deg C].
!     *REAL*    *psao*       - salinity [psu.].
!     *REAL*    *ppao*       - sea level pressure [Pascal].
!     *REAL*    *prho*       - density [kg/m^3].
!     *REAL*    *pddpo*      - size of scalar grid cell (depth) [m].
!     *REAL*    *pdlxp*      - size of scalar grid cell (longitudinal) [m].
!     *REAL*    *pdlyp*      - size of scalar grid cell (latitudinal) [m].
!     *REAL*    *ptiestu*    - 
!     *REAL*    *ptiestw*    - 
!     *REAL*    *pfu10*      - 
!     *INTEGER* *kmonlen*    - length of current month in days.
!     *INTEGER* *kldtmon*    - monthly time stap in OCE.
!     *INTEGER* *kldtday*    - daily time stap in OCE.
!     *REAL*    *omask*      - land/ocean mask
!     *INTEGER* *days_in_yr* - number of days in year
!
!**********************************************************************

      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      USE mo_bgcmean
      USE mo_control_bgc
      use mo_param1_bgc 
      use mod_xc
#ifdef DIFFAT
      use mo_satm
#endif
#ifdef RIV_GNEWS
      use mo_riverinpt
#endif
      use mo_ndep, only: n_deposition


      implicit none

      INTEGER :: kpie,kpje,kpke
      REAL    :: pglat  (kpie,kpje)
      REAL    :: pfswr  (kpie,kpje)
      REAL    :: psicomo(kpie,kpje)
      REAL    :: pfu10  (kpie,kpje)
      REAL    :: patmco2(kpie,kpje)
      REAL    :: pflxco2(kpie,kpje)
      REAL    :: pflxdms(kpie,kpje)
      REAL    :: ptho   (kpie,kpje,kpke)
      REAL    :: psao   (kpie,kpje,kpke)
      REAL    :: ppao   (kpie,kpje)
      REAL    :: prho   (kpie,kpje,kpke)
      REAL    :: pddpo  (kpie,kpje,kpke)
      REAL    :: pdlxp  (kpie,kpje)
      REAL    :: pdlyp  (kpie,kpje)
      REAL    :: ptiestu(kpie,kpje,kpke+1)
      REAL    :: ptiestw(kpie,kpje,kpke+1)
      REAL    :: omask  (kpie,kpje)
      INTEGER :: kplyear,kplmon,kplday,kmonlen,kldtmon,kldtday
      INTEGER :: days_in_yr

      INTEGER :: i,j,k,l
      REAL    :: emissions      

      IF (mnproc.eq.1) THEN
      write(io_stdo_bgc,*) 'HAMOCC',KLDTDAY,KLDTMON,LDTRUNBGC,NDTDAYBGC
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
! set limits for temp and saln
!
      DO k=1,kpke
!$OMP PARALLEL DO
      DO  j=1,kpje
      DO  i=1,kpie
        ptho(i,j,k)=min(40.,max(-3.,ptho(i,j,k)))
        psao(i,j,k)=min(40.,max( 0.,psao(i,j,k)))
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      ENDDO


!--------------------------------------------------------------------
! Net solar radiation 
!
!$OMP PARALLEL DO
      DO  j=1,kpje
      DO  i=1,kpie
        strahl(i,j)=pfswr(i,j)
      ENDDO
      ENDDO
!$OMP END PARALLEL DO


!--------------------------------------------------------------------
! Pass atmospheric co2
!
#if defined(PROGCO2) || defined(DIAGCO2)
!$OMP PARALLEL DO
      DO  j=1,kpje
      DO  i=1,kpie
        atm(i,j,iatmco2)=patmco2(i,j)
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
         if (mnproc.eq.1) then 
           write (io_stdo_bgc,*) 'jt: getting x2o co2'
         endif

#else
!$OMP PARALLEL DO
      DO  j=1,kpje
      DO  i=1,kpie
        if (patmco2(i,j).lt.0) then
          atm(i,j,iatmco2)=atm_co2
        else
          atm(i,j,iatmco2)=patmco2(i,j)
        endif
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
#endif


!--------------------------------------------------------------------
! Read atmospheric cfc concentrations
!
#ifdef CFC
      call get_cfc(kplyear,atm_cfc11_nh,atm_cfc12_nh,atm_sf6_nh,        &
                           atm_cfc11_sh,atm_cfc12_sh,atm_sf6_sh)
#endif


!---------------------------------------------------------------------
! Read emission data
!
#ifdef EMS_CO2
#ifdef DIFFAT 
      IF (kldtmon.eq.1.and.kldtmon.eq.1) THEN 
             
         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*) 'CO2_EMS gerufen bei kldtmon: ',          &
     &                         kplmon,kplday,kmonlen,kldtmon
         ENDIF
      call co2_ems(kplyear,emissions)
!      emissions= 2000.
      emission = emissions/1000. ! million metric tons --> GigaTons
      ems_per_step=(1.e12/12.)*emission/(float(days_in_yr)*20.)
      write(io_stdo_bgc,*) 'CO2_EMS',emissions,emission,ems_per_step
      write(io_stdo_bgc,*) (1.e12/12.)*emission,float(days_in_yr)*20.

      ENDIF ! First timestep of the month
#endif
#endif


#ifdef PBGC_CK_TIMESTEP
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'before BGC: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif


!---------------------------------------------------------------------
! Recalculate the bottom-most mass containing layer

      call calc_bot(kpie,kpje,kpke,pddpo)


!---------------------------------------------------------------------
!     Biogeochemistry
      CALL OCPROD(kpie,kpje,kpke,ptho,pddpo,pdlxp,pdlyp,ptiestu,        &
     &            ptiestw,kplmon,omask)

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


      CALL CYANO(kpie,kpje,kpke,ptho,pddpo,omask)

#ifdef PBGC_CK_TIMESTEP   
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after CYANO: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif


      CALL CARCHM(kpie,kpje,kpke,pglat,pddpo,pdlxp,pdlyp,psao,ppao,     &
     &            ptho,prho,psicomo,pfu10,ptiestu,omask)

#ifdef PBGC_CK_TIMESTEP   
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after CARCHM: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif


      ! Apply n-deposition
      CALL n_deposition(kpie,kpje,kpke,kplyear,kplmon,pddpo,omask)

#ifdef RIV_GNEWS
      ! Apply riverine input of carbon and nutrients
      call riverinpt(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,omask)
#endif


#ifdef DIFFAT     
      CALL SATM_STEP(atmflx,atm)
#endif	 

#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after ATMOTR: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif	 

!     update preformed tracers
      CALL PREFTRC(kpie,kpje,omask)


!--------------------------------------------------------------------
!     Sediment module

#ifndef sedbypass
! jump over sediment if sedbypass is defined

      CALL POWACH(kpie,kpje,kpke,psao,prho,omask)

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
! Accumulate fields and write output

      CALL ACCFIELDS(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask)       


!--------------------------------------------------------------------
! Pass co2 flux. Convert unit from kmol/m^2 to kg/m^2/s.

!$OMP PARALLEL DO
      DO  j=1,kpje
      DO  i=1,kpie
        pflxco2(i,j)=-44.*atmflx(i,j,iatmco2)/dtbgc
      ENDDO
      ENDDO
!$OMP END PARALLEL DO


!--------------------------------------------------------------------
! Pass dms flux. Convert unit from kmol/m^2 to kg/m^2/s.

!$OMP PARALLEL DO
      DO  j=1,kpje
      DO  i=1,kpie
        pflxdms(i,j)=-62.13*atmflx(i,j,iatmdms)/dtbgc
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

!--------------------------------------------------------------------
      RETURN
      END
