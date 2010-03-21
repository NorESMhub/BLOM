      SUBROUTINE HAMOCC4BCM(kpie,kpje,kpke,kpbe,                        &
     &    pfswr,psicomo,ptho,psao,pddpo,pdlxp,pdlyp,ptiestu,ptiestw,    &
     &    pdpio,pfu10,patmco2,pflxco2,kplyear,kplmon,kplday,kmonlen,    &
     &    kldtmon,kldtday,omask,dummy_tr,ntr,ntrbgc,itrbgc,             &
     &    days_in_yr)       

!**********************************************************************
!
!**** *BGC* - .
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st REAL of model grid.
!     *INTEGER* *kpje*    - 2nd REAL of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) REAL of model grid.
!     *REAL*    *pfswr*   - solar radiation [W/m**2].
!     *REAL*    *psicomo* - sea ice concentration
!     *REAL*    *ptho*    - potential temperature [deg C].
!     *REAL*    *psao*    - salinity [psu.].
!     *REAL*    *ppao*    - sea level pressure [Pascal].
!     *REAL*    *pddpo*   - size of scalar grid cell (depth) [m].
!     *REAL*    *pdlxp*   - size of scalar grid cell (longitudinal) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (latitudinal) [m].
!     *REAL*    *pdpio*   - inverse size of grid cell (depth)[m].
!     *INTEGER* *kldtmon* - monthly time stap in OCE.
!     *INTEGER* *kldtday* - daily time stap in OCE.
!
!**********************************************************************

      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      USE mo_bgcmean
      USE mo_control_bgc
!      USE mo_timeser_bgc
      use mo_param1_bgc 
#ifdef DIFFAT
      use mo_satm
#endif
#ifdef PDYNAMIC_BGC
      use mo_dynamic
#endif /* PDYNAMIC_BGC */ 
      use mod_xc

      implicit none
      INTEGER :: kpie,kpje,kpke,kpbe,i,j,k,l,ntr,ntrbgc,itrbgc

      REAL pfswr  (kpie,kpje)
      REAL psicomo(kpie,kpje)
      REAL pfu10(kpie,kpje)
      REAL patmco2(kpie,kpje)
      REAL pflxco2(kpie,kpje)
 
      REAL ptho (kpie,kpje,kpke)
      REAL psao (kpie,kpje,kpke)
      REAL pddpo(kpie,kpje,kpke)
      REAL pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL pdpio(kpie,kpje,kpke)
!      REAL ppao(kpie,kpje)            
      REAL ptiestu(kpie,kpje,kpke+1)
      REAL ptiestw(kpie,kpje,kpke+1)
!      REAL zo(kpie,kpje),sicsno(kpie,kpje),sictho(kpie,kpje)
      REAL omask(kpie,kpje)
      REAL dummy_tr(1-kpbe:kpie+kpbe,1-kpbe:kpje+kpbe,kpke,ntr)
      INTEGER :: kplyear,kplmon,kplday,kmonlen,kldtmon,kldtday
      INTEGER :: days_in_yr

      REAL emissions

      IF (mnproc.eq.1) THEN
      write(io_stdo_bgc,*) 'HAMOCC',KLDTDAY,KLDTMON,LDTRUNBGC,NDTDAYBGC
      ENDIF

!
! Increment bgc time step counter of run (initialized in INI_BGC).
!
      ldtrunbgc = ldtrunbgc + 1
!
! Increment bgc time step counter of experiment (initialized if IAUFR=0).
!
      ldtbgc = ldtbgc + 1
!
! Increment timeseries-1 sample counter (initialized in INI_TIMSER_BGC).
!

!--------------------------------------------------------------------
! pass tracer fields in from ocean model
!
      ocetra(:,:,:,:)=dummy_tr(1:kpie,1:kpje,:,itrbgc:itrbgc+ntrbgc-1)

      DO k=1,kpke
!$OMP PARALLEL DO
      DO  j=1,kpje
      DO  i=1,kpie
        ptho(i,j,k)=min(40.,max(-3.,ptho(i,j,k)))
        psao(i,j,k)=min(40.,max(0.,psao(i,j,k)))
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      ENDDO
!
!--------------------------------------------------------------------
! Net solar radiation: multiply  with sea ice concentration

!$OMP PARALLEL DO
      DO  j=1,kpje
      DO  i=1,kpie
        strahl(i,j)=pfswr(i,j)*(1.-MIN(psicomo(i,j),0.9))
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
!
!--------------------------------------------------------------------
! Pass atmospheric co2

#if defined(DIFFAT) || defined(CCSMCOUPLED)
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

#ifdef PBGC_CK_TIMESTEP
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'before BGC: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask &
     &                  ,0)
!      call INVENTORY_BGC(kpie,kpje,kpke)
#endif

!
!--------------------------------------------------------------------
! Get flux or partial pressure of atmospheric CO2.
!
!      IF (MOD(ldtrunbgc-1,ndtdaybgc) .EQ. 0) THEN
!         CALL GET_ATMCO2
!      ENDIF
!_______________________________________________________________
!
! First timestep of the month --> kldtmon eq 1
      
      IF (kldtmon .EQ. 1) THEN 
             
! Calculate the chemical constants of the month.
!         IF (mnproc.eq.1) THEN
!         WRITE(io_stdo_bgc,*) 'CHEMCON gerufen bei kldtmon: ',          &
!     &                         kplmon,kplday,kmonlen,kldtmon
!         END
         CALL CHEMCON(-26,kpie,kpje,kpke,psao,ptho,                     &
     &                 pddpo,pdlxp,pdlyp,ptiestu,kplmon,omask)     

      ENDIF ! First timestep of the month
!      call chksumbgc(satoxy,kpke,'satoxy')
!
!_______________________________________________________________
!
! First timestep of the year
      
#ifdef EMS_CO2
#ifdef DIFFAT 
      IF (kldtmon.eq.1.and.kldtmon.eq.1) THEN 
             
! Calculate the chemical constants of the month.
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
!
!
!_______________________________________________________________
! Recalculate the bottom-most mass containing layer

      call calc_bot(kpie,kpje,kpke,pddpo)
!_______________________________________________________________
!
!
!     Biogeochemistry

      CALL OCPROD(kpie,kpje,kpke,ptho,pddpo,pdlxp,pdlyp,pdpio,ptiestu, &
     &            ptiestw,kplmon,omask)

!      write(io_stdo_bgc,*) 'OCPROD'

#ifdef PBGC_CK_TIMESTEP   
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after OCPROD: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask &
     &                  ,0)
!      call INVENTORY_BGC(kpie,kpje,kpke)
#endif

!      call chksumbgc(ocetra(1,1,1,isco212),kpke,'sco212')
!      call chksumbgc(ocetra(1,1,1,ialkali),kpke,'alkali')
!      call chksumbgc(ocetra(1,1,1,iphosph),kpke,'phosph')
!      call chksumbgc(ocetra(1,1,1,ioxygen),kpke,'oxygen')
!      call chksumbgc(ocetra(1,1,1,igasnit),kpke,'gasnit')
!      call chksumbgc(ocetra(1,1,1,iano3),kpke,'ano3')
!      call chksumbgc(ocetra(1,1,1,isilica),kpke,'silica')
!      call chksumbgc(ocetra(1,1,1,idoc),kpke,'doc')
!      call chksumbgc(ocetra(1,1,1,iphy),kpke,'phy')
!      call chksumbgc(ocetra(1,1,1,izoo),kpke,'zoo')
!      call chksumbgc(ocetra(1,1,1,ian2o),kpke,'ian2o')
!      call chksumbgc(ocetra(1,1,1,idms),kpke,'dms')
!      call chksumbgc(ocetra(1,1,1,iiron),kpke,'iron')
!      call chksumbgc(ocetra(1,1,1,ifdust),kpke,'dust')
!      call chksumbgc(ocetra(1,1,1,isco213),kpke,'sco213')
!      call chksumbgc(ocetra(1,1,1,isco214),kpke,'sco214')
!      call chksumbgc(ocetra(1,1,1,idet),kpke,'det')
!      call chksumbgc(ocetra(1,1,1,icalc),kpke,'calc')
!      call chksumbgc(ocetra(1,1,1,iopal),kpke,'opal')
!      call chksumbgc(ocetra(1,1,1,idet13),kpke,'det13')
!      call chksumbgc(ocetra(1,1,1,icalc13),kpke,'calc13')
!      call chksumbgc(ocetra(1,1,1,idet14),kpke,'det14')
!      call chksumbgc(ocetra(1,1,1,icalc14),kpke,'calc14')
 
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
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask &
     &                  ,0)
!      call INVENTORY_BGC(kpie,kpje,kpke)
#endif

      CALL CYANO(kpie,kpje,kpke,pddpo,omask)

!      write(io_stdo_bgc,*) 'CYANO'

#ifdef PBGC_CK_TIMESTEP   
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after CYANO: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask &
     &                  ,0)
!      call INVENTORY_BGC(kpie,kpje,kpke)
#endif

      CALL CARCHM                                                      &
     &(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,psao,ptho,psicomo,             &
     & pfu10,kplyear,kplmon,kplday,ndtdaybgc,ldtrunbgc,ptiestu,        &
     & kmonlen,omask)
!      CALL CARCHM                                                      &
!     &(kpie,kpje,kpke,pddpo,psao,ptho,psicomo,                         &
!     & pfu10,kplmon,kplday,kmonlen,omask)

!      write(io_stdo_bgc,*) 'CARCHM'

#ifdef PBGC_CK_TIMESTEP   
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after CARCHM: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask &
     &                  ,0)
!      call INVENTORY_BGC(kpie,kpje,kpke)
#endif

!      
#ifdef DIFFAT     
!      write(lp,*) 'SATM_STEP'
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
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask &
     &                  ,0)
!      call INVENTORY_BGC(kpie,kpje,kpke)
#endif	 
!_______________________________________________________________
!
!     Sediment module

      CALL POWACH(kpie,kpje,kpke,pdlxp,pdlyp,psao,omask)
!
#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after POWACH: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask &
     &                  ,0)
!      call INVENTORY_BGC(kpie,kpje,kpke)
#endif	 

      IF(KLDTDAY .EQ. 1) THEN
         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)                                          &
     &   'Sediment shifting ...'
         ENDIF

         CALL SEDSHI(kpie,kpje,omask)

      ENDIF

!     accumulate sediments
      DO l=1,nbgc 
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
      ENDDO 

!
#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after BGC: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask &
     &                  ,0)
!      call INVENTORY_BGC(kpie,kpje,kpke)
#endif	 

!--------------------------------------------------------------------
! pass tracer fields out to ocean model 
!
      dummy_tr(1:kpie,1:kpje,:,itrbgc:itrbgc+ntrbgc-1)=ocetra(:,:,:,:)

!
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


!---------------------------------------------------------------------
! write output files
!
      DO l=1,nbgc 
        nacc_bgc(l)=nacc_bgc(l)+1
        if (bgcwrt(l).gt.0.5) then
          if (GLB_INVENTORY(l).ne.0)                                    & 
     &      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
          call ncwrt_bgc(l)
          nacc_bgc(l)=0 
        endif
      ENDDO
!
!--------------------------------------------------------------------

      RETURN
      END
