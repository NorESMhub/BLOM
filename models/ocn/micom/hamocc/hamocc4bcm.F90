      SUBROUTINE HAMOCC4BCM(kpie,kpje,kpke,kpbe,                         &
     &    pfswr,psicomo,ptho,psao,ppao,prho,pddpo,pdlxp,pdlyp,ptiestu,   &
     &    ptiestw,pdpio,pfu10,patmco2,pflxco2,kplyear,kplmon,kplday,     &
     &    kmonlen,kldtmon,kldtday,omask,dummy_tr,ntr,ntrbgc,itrbgc,      &
     &    days_in_yr)       

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
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*       - 1st dimension of model grid.
!     *INTEGER* *kpje*       - 2nd dimension of model grid.
!     *INTEGER* *kpke*       - 3rd (vertical) dimension of model grid.
!     *INTEGER* *kpbe*       - nb. of halo rows for tracer field
!     *REAL*    *pfswr*      - solar radiation [W/m**2].
!     *REAL*    *psicomo*    - sea ice concentration
!     *REAL*    *ptho*       - potential temperature [deg C].
!     *REAL*    *psao*       - salinity [psu.].
!     *REAL*    *ppao*       - sea level pressure [Pascal].
!     *REAL*    *prho*       - density [kg/m^3].
!     *REAL*    *pddpo*      - size of scalar grid cell (depth) [m].
!     *REAL*    *pdlxp*      - size of scalar grid cell (longitudinal) [m].
!     *REAL*    *pdlyp*      - size of scalar grid cell (latitudinal) [m].
!     *REAL*    *pdpio*      - inverse size of grid cell (1/depth)[1/m].
!     *INTEGER* *kmonlen*    - length of current month in days.
!     *INTEGER* *kldtmon*    - monthly time stap in OCE.
!     *INTEGER* *kldtday*    - daily time stap in OCE.
!     *REAL*    *omask*      - land/ocean mask
!     *REAL*    *dummy_tr*   - initial/restart tracer field to be passed to the 
!                              ocean model [mol/kg]
!     *INTEGER* *ntr*        - number of tracers in tracer field
!     *INTEGER* *ntrbgc*     - number of biogechemical tracers in tracer field
!     *INTEGER* *itrbgc*     - start index for biogeochemical tracers in tracer field
!     *INTEGER* *days_in_yr* - number of days in year
!
!**********************************************************************

      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      USE mo_bgcmean
      USE mo_control_bgc
!      USE mo_timeser_bgc
      use mo_param1_bgc 
      use mod_xc
#ifdef DIFFAT
      use mo_satm
#endif
#ifdef RIV_GNEWS
      use mo_riverinpt
#endif

      implicit none

      INTEGER :: kpie,kpje,kpke,kpbe,ntr,ntrbgc,itrbgc
      REAL    :: pfswr  (kpie,kpje)
      REAL    :: psicomo(kpie,kpje)
      REAL    :: pfu10  (kpie,kpje)
      REAL    :: patmco2(kpie,kpje)
      REAL    :: pflxco2(kpie,kpje)
      REAL    :: ptho   (kpie,kpje,kpke)
      REAL    :: psao   (kpie,kpje,kpke)
      REAL    :: ppao   (kpie,kpje)
      REAL    :: prho   (kpie,kpje,kpke)
      REAL    :: pddpo  (kpie,kpje,kpke)
      REAL    :: pdlxp  (kpie,kpje)
      REAL    :: pdlyp  (kpie,kpje)
      REAL    :: pdpio  (kpie,kpje,kpke)
      REAL    :: ptiestu(kpie,kpje,kpke+1)
      REAL    :: ptiestw(kpie,kpje,kpke+1)
      REAL    :: omask  (kpie,kpje)
      REAL    :: dummy_tr(1-kpbe:kpie+kpbe,1-kpbe:kpje+kpbe,kpke,ntr)
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
! pass tracer fields in from ocean model; convert mol/kg -> kmol/m^3
!
      ocetra(:,:,:,:)=dummy_tr(1:kpie,1:kpje,:,itrbgc:itrbgc+ntrbgc-1)

      do k=1,kpke
!$OMP PARALLEL DO
      do j=1,kpje
      do i=1,kpie
        if (omask(I,J) .gt. 0.5 ) then
          ocetra(i,j,k,:)=ocetra(i,j,k,:)*prho(i,j,k)
        endif
      enddo
      enddo
!$OMP END PARALLEL DO
      enddo


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
! Net solar radiation: multiply  with sea ice concentration
!
!$OMP PARALLEL DO
      DO  j=1,kpje
      DO  i=1,kpie
        strahl(i,j)=pfswr(i,j)*(1.-MIN(psicomo(i,j),0.9))
      ENDDO
      ENDDO
!$OMP END PARALLEL DO


!--------------------------------------------------------------------
! Pass atmospheric co2
!
#if defined(DIFFAT) || defined(CCSMCOUPLED)
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
#endif 


!--------------------------------------------------------------------
! Read atmospheric cfc concentrations
!
#ifdef CFC
      call get_cfc(kplyear,atm_cfc11,atm_cfc12,atm_sf6)
      IF (mnproc.EQ.1) THEN
      write(*,*)'Carchm getcfc:',kplyear,atm_cfc11,atm_cfc12,atm_sf6
      ENDIF
#endif


!---------------------------------------------------------------------
! Read emission data
!
#ifdef EMS_CO2
#ifdef DIFFAT 
      IF (kldtmon.eq.1.and.kldtmon.eq.1) THEN 
             
         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*) 'CO2_EMS gerufen bei kldtmon: ',       &
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

      CALL OCPROD(kpie,kpje,kpke,ptho,pddpo,pdlxp,pdlyp,pdpio,ptiestu, &
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


      CALL CYANO(kpie,kpje,kpke,pddpo,omask)

#ifdef PBGC_CK_TIMESTEP   
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after CYANO: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif


      CALL CARCHM(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,psao,ppao, &
     &            ptho,prho,psicomo,pfu10,ptiestu,omask)

#ifdef PBGC_CK_TIMESTEP   
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after CARCHM: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif


#ifdef RIV_GNEWS
      ! Apply riverine input of carbon and nutrients
      call riverinpt(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,omask)
#endif

!      
#ifdef DIFFAT     
      CALL SATM_STEP(atmflx,atm)
      call accsrf(jatmco2,atm(1,1,iatmco2),omask,0)
      call accsrf(jatmo2 ,atm(1,1,iatmo2),omask,0)
      call accsrf(jatmn2 ,atm(1,1,iatmn2),omask,0)
#elif defined(CCSMCOUPLED)
      call accsrf(jatmco2,atm(1,1,iatmco2),omask,0)
#endif
!
#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after ATMOTR: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif	 

!--------------------------------------------------------------------
!     Sediment module

      CALL POWACH(kpie,kpje,kpke,pdlxp,pdlyp,psao,omask)
!
#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after POWACH: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif	 

      IF(KLDTDAY .EQ. 1) THEN
         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)                                          &
     &   'Sediment shifting ...'
         ENDIF

         CALL SEDSHI(kpie,kpje,omask)

      ENDIF

!     accumulate sediments
      call accsdm(jpowaic,powtra(1,1,1,ipowaic))
      call accsdm(jpowaal,powtra(1,1,1,ipowaal))
      call accsdm(jpowaph,powtra(1,1,1,ipowaph))
      call accsdm(jpowaox,powtra(1,1,1,ipowaox))
      call accsdm(jpown2 ,powtra(1,1,1,ipown2) )
      call accsdm(jpowno3,powtra(1,1,1,ipowno3))
      call accsdm(jpowasi,powtra(1,1,1,ipowasi))
      call accsdm(jssso12,sedlay(1,1,1,issso12))
      call accsdm(jssssil,sedlay(1,1,1,issssil))
      call accsdm(jsssc12,sedlay(1,1,1,isssc12))
      call accsdm(jssster,sedlay(1,1,1,issster))

!     accumulate sediment burial
      call accbur(jburssso12,burial(1,1,issso12))
      call accbur(jburssssil,burial(1,1,issssil))
      call accbur(jbursssc12,burial(1,1,isssc12))
      call accbur(jburssster,burial(1,1,issster))


#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after BGC: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif	 


!---------------------------------------------------------------------
! write output files

      DO l=1,nbgc 
        nacc_bgc(l)=nacc_bgc(l)+1
        if (bgcwrt(l).gt.0.5) then
          if (GLB_INVENTORY(l).ne.0)                                    & 
     &      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
          call ncwrt_bgc(l)
          nacc_bgc(l)=0 
        endif
      ENDDO

!--------------------------------------------------------------------
! Pass co2 flux. Convert unit from kmol/m^2 to kg/m^2/s.

#if defined(DIFFAT) || defined(CCSMCOUPLED)
!$OMP PARALLEL DO
      DO  j=1,kpje
      DO  i=1,kpie
        pflxco2(i,j)=-44.*atmflx(i,j,iatmco2)/dtbgc
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
#endif	 


!--------------------------------------------------------------------
! pass tracer fields out to ocean model; convert kmol/m^3 -> mol/kg

      do k=1,kpke
!$OMP PARALLEL DO
      do j=1,kpje
      do i=1,kpie
        if (omask(i,j) .gt. 0.5 ) then
          ocetra(i,j,k,:)=ocetra(i,j,k,:)/prho(i,j,k)
        endif
      enddo
      enddo
!$OMP END PARALLEL DO
      enddo

      dummy_tr(1:kpie,1:kpje,:,itrbgc:itrbgc+ntrbgc-1)=ocetra(:,:,:,:)

!--------------------------------------------------------------------
      RETURN
      END
