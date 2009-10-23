      SUBROUTINE HAMOCC4BCM(kpie,kpje,kpke,kpbe,                       &
     &    pfswr,psicomo,ptho,psao,pddpo,pdlxp,pdlyp,ptiestu,ptiestw,   &
     &    pdpio,pfu10,patmco2,pflxco2,kplyear,kplmon,kplday,kmonlen,   &
     &    kldtmon,kldtday,omask,dummy_tr,ndtr,bgc3dwrt,bgc2dwrt,       &
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
      INTEGER :: kpie,kpje,kpke,kpbe,i,j,k,l,ndtr,ntr

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
      REAL dummy_tr(1-kpbe:kpie+kpbe,1-kpbe:kpje+kpbe,kpke,ndtr)
      REAL bgc3dwrt,bgc2dwrt
      INTEGER :: kplyear,kplmon,kplday,kmonlen,kldtmon,kldtday
      INTEGER :: days_in_yr

      REAL emissions

      IF (mnproc.eq.1) THEN
      write(io_stdo_bgc,*) 'HAMOCC'
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
      ocetra(:,:,:,:)=dummy_tr(1:kpie,1:kpje,:,:)

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
      CALL SATM_STEP(atmflx,atm)
!      write(lp,*) 'SATM_STEP'
!$OMP PARALLEL DO 
      do j=1,kpje
      do i=1,kpie
        bgcm2d(i,j,jatmco2)=bgcm2d(i,j,jatmco2)+atm(i,j,iatmco2)
        bgcm2d(i,j,jatmo2) =bgcm2d(i,j,jatmo2) +atm(i,j,iatmo2)
        bgcm2d(i,j,jatmn2) =bgcm2d(i,j,jatmn2) +atm(i,j,iatmn2)
      enddo
      enddo
!$OMP END PARALLEL DO
!      call xcstop('(hamocc4bcn)')
!             stop('(hamocc4bcn)')
#elif CCSMCOUPLED
!$OMP PARALLEL DO 
      do j=1,kpje
      do i=1,kpie
        bgcm2d(i,j,jatmco2)=bgcm2d(i,j,jatmco2)+atm(i,j,iatmco2)
      enddo
      enddo
!$OMP END PARALLEL DO
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

      IF(MOD(ldtrunbgc,ndtdaybgc) .EQ. 0 ) THEN
         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)                                          &
     &   'Sediment shifting at runstep ',ldtrunbgc,' ...'
         ENDIF

         CALL SEDSHI(kpie,kpje,omask)

      ENDIF

!js sediment
      DO k=1,ks
!$OMP PARALLEL DO 
      DO j=1,kpje
      DO i=1,kpie
         IF(omask(i,j).GT.0.5) THEN
             bgct_sed(i,j,k,jpowaic)  =                      &
     &       bgct_sed(i,j,k,jpowaic)  + powtra(i,j,k,ipowaic)
             bgct_sed(i,j,k,jpowaal)  =                      &
     &       bgct_sed(i,j,k,jpowaal)  + powtra(i,j,k,ipowaal)
             bgct_sed(i,j,k,jpowaph)  =                      &
     &       bgct_sed(i,j,k,jpowaph)  + powtra(i,j,k,ipowaph)
             bgct_sed(i,j,k,jpowaox)  =                      &
     &       bgct_sed(i,j,k,jpowaox)  + powtra(i,j,k,ipowaox)
             bgct_sed(i,j,k,jpown2)   =                      &
     &       bgct_sed(i,j,k,jpown2)   + powtra(i,j,k,ipown2)
             bgct_sed(i,j,k,jpowno3)  =                      &
     &       bgct_sed(i,j,k,jpowno3)  + powtra(i,j,k,ipowno3)
             bgct_sed(i,j,k,jpowasi)  =                      &
     &       bgct_sed(i,j,k,jpowasi)  + powtra(i,j,k,ipowasi)
             bgct_sed(i,j,k,jssso12)  =                      &
     &       bgct_sed(i,j,k,jssso12)  + sedlay(i,j,k,issso12)
             bgct_sed(i,j,k,jssssil)  =                      &
     &       bgct_sed(i,j,k,jssssil)  + sedlay(i,j,k,issssil)
             bgct_sed(i,j,k,jsssc12)  =                      &
     &       bgct_sed(i,j,k,jsssc12)  + sedlay(i,j,k,isssc12)
             bgct_sed(i,j,k,jssster)  =                      &
     &       bgct_sed(i,j,k,jssster)  + sedlay(i,j,k,issster)
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
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
      dummy_tr(1:kpie,1:kpje,:,:)=ocetra(:,:,:,:)
!
!--------------------------------------------------------------------
! Pass co2 flux. Convert unit from kmol/m^2 to kg/m^2/s.

!$OMP PARALLEL DO
      DO  j=1,kpje
      DO  i=1,kpie
        pflxco2(i,j)=-44.*atmflx(i,j,iatmco2)/dtbgc
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

!---------------------------------------------------------------------
! write output files
!
      nacc_bgc2d=nacc_bgc2d+1
      if (bgc3dwrt.gt.0.5) then
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask &
     &                  ,0)
      call ncwrt3_bgc(bgcm3d)
      call ncwrt3_sed(bgct_sed,nacc_bgc2d)
      endif
      if (bgc2dwrt.gt.0.5) then
!      do i=1,kpie
!      do j=1,kpje
!         if(omask(i,j).gt.0.5) then
!           bgcm2d(i,j,1)=float(kwrbioz(i,j))
!           bgcm2d(i,j,2)=ptiestu(i,j,kwrbioz(i,j))
!           bgcm2d(i,j,3)=ptiestw(i,j,kwrbioz(i,j)+1)
!         endif
!      enddo
!      enddo
      call ncwrt2_bgc(bgcm2d,nbgcm2d,nacc_bgc2d)
      nacc_bgc2d=0 
      endif
!      write(io_stdo_bgc,*) 'WRITE'
!
!--------------------------------------------------------------------

      RETURN
      END
