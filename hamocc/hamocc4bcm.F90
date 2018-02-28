      SUBROUTINE HAMOCC4BCM(kpie,kpje,kpke,pglat,                        &
     &    pfswr,psicomo,ptho,psao,ppao,prho,pddpo,pdlxp,pdlyp,ptiestu,   &
     &    ptiestw,pdpio,pfu10,patmco2,pflxco2,kplyear,kplmon,kplday,     &
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
!     J.Schwinger       *GFI, Bergen*    2014-05-21
!     - moved copying of tracer field to ocetra to micom2hamocc 
!       and hamocc2micom
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
!     *REAL*    *pdpio*      - inverse size of grid cell (1/depth)[1/m].
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
!      USE mo_timeser_bgc
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
      REAL    :: pdpio  (kpie,kpje,kpke)
      REAL    :: ptiestu(kpie,kpje,kpke+1)
      REAL    :: ptiestw(kpie,kpje,kpke+1)
      REAL    :: omask  (kpie,kpje)
      INTEGER :: kplyear,kplmon,kplday,kmonlen,kldtmon,kldtday
      INTEGER :: days_in_yr

      INTEGER :: i,j,k,l
      INTEGER :: ind1(kpie,kpje), ind2(kpie,kpje)
      REAL    :: wghts(kpie,kpje,ddm)
      REAL    :: emissions      
#ifdef cisonew
      REAL, DIMENSION(kpie,kpje,kpke) :: d13c,d14c,bigd14c
#endif
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
      CALL n_deposition(kpie,kpje,kpke,kplyear,kplmon,pddpo,omask)
      CALL OCPROD(kpie,kpje,kpke,ptho,pddpo,pdlxp,pdlyp,pdpio,ptiestu,  &
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

!---------------------------------------------------------------------
!     Accumulate global fields and write output files (note: should 
!     eventually be moved to own subroutine)

#ifdef cisonew
! Calculation d13C and Dd14C: Delta notation for output
      d13c(:,:,:)=0.
      d14c(:,:,:)=0.
      bigd14c(:,:,:)=0.
      do k=1,kpke
      do j=1,kpje
      do i=1,kpie
        if(omask(i,j).gt.0.5.and.pddpo(i,j,k).gt.dp_min) then
             d13c(i,j,k)=                        &
     &    (ocetra(i,j,k,isco213)/(ocetra(i,j,k,isco212)-ocetra(i,j,k,isco213)+safediv)&
     &        /re1312-1.)*1000.

             d14c(i,j,k)=                        &
     &    (ocetra(i,j,k,isco214)*c14fac/(ocetra(i,j,k,isco212)+safediv) &
     &        /re14to-1.)*1000.

             bigd14c(i,j,k)=d14c(i,j,k)-2.*(d13c(i,j,k)+25.)*(1.+d14c(i,j,k)/1000.)
        endif
      enddo
      enddo
      enddo
#endif

!     Accumulate atmosphere fields
      call accsrf(jatmco2,atm(1,1,iatmco2),omask,0)
      call accsrf(jn2ofx,atmflx(1,1,iatmn2o),omask,0)
#ifdef DIFFAT     
      call accsrf(jatmo2 ,atm(1,1,iatmo2),omask,0)
      call accsrf(jatmn2 ,atm(1,1,iatmn2),omask,0)
#endif

!     Accumulate srf diagnostics
      call accsrf(jsrfphosph,ocetra(1,1,1,iphosph),omask,0)
      call accsrf(jsrfoxygen,ocetra(1,1,1,ioxygen),omask,0)
      call accsrf(jsrfiron,ocetra(1,1,1,iiron),omask,0)
      call accsrf(jsrfano3,ocetra(1,1,1,iano3),omask,0)
      call accsrf(jsrfalkali,ocetra(1,1,1,ialkali),omask,0)
      call accsrf(jsrfsilica,ocetra(1,1,1,isilica),omask,0)
      call accsrf(jsrfdic,ocetra(1,1,1,isco212),omask,0)

!     Accumulate layer diagnostics
      call acclyr(jdp,pddpo,pddpo,0)
      call acclyr(jphyto,ocetra(1,1,1,iphy),pddpo,1)   
      call acclyr(jgrazer,ocetra(1,1,1,izoo),pddpo,1) 
      call acclyr(jphosph,ocetra(1,1,1,iphosph),pddpo,1)
      call acclyr(joxygen,ocetra(1,1,1,ioxygen),pddpo,1)
      call acclyr(jiron,ocetra(1,1,1,iiron),pddpo,1)    
      call acclyr(jano3,ocetra(1,1,1,iano3),pddpo,1)    
      call acclyr(jalkali,ocetra(1,1,1,ialkali),pddpo,1)
      call acclyr(jsilica,ocetra(1,1,1,isilica),pddpo,1)
      call acclyr(jdic,ocetra(1,1,1,isco212),pddpo,1)    
      call acclyr(jdoc,ocetra(1,1,1,idoc),pddpo,1)       
      call acclyr(jpoc,ocetra(1,1,1,idet),pddpo,1)       
      call acclyr(jcalc,ocetra(1,1,1,icalc),pddpo,1)    
      call acclyr(jopal,ocetra(1,1,1,iopal),pddpo,1)    
      call acclyr(jn2o,ocetra(1,1,1,ian2o),pddpo,1) 
      call acclyr(jco3,co3,pddpo,1)                      
      call acclyr(jph,hi,pddpo,1)
      call acclyr(jomegac,OmegaC,pddpo,1)
#ifdef natDIC
      call acclyr(jnatalkali,ocetra(1,1,1,inatalkali),pddpo,1)
      call acclyr(jnatdic,ocetra(1,1,1,inatsco212),pddpo,1)
      call acclyr(jnatcalc,ocetra(1,1,1,inatcalc),pddpo,1)
      call acclyr(jnatco3,natco3,pddpo,1)                      
      call acclyr(jnatomegac,natOmegaC,pddpo,1)
#endif
#ifdef cisonew
      call acclyr(jdic13,ocetra(1,1,1,isco213),pddpo,1)    
      call acclyr(jdic14,ocetra(1,1,1,isco214),pddpo,1)    
      call acclyr(jd13c,d13c,pddpo,1)    
      call acclyr(jpoc13,ocetra(1,1,1,idet13),pddpo,1)
      call acclyr(jdoc13,ocetra(1,1,1,idoc13),pddpo,1)
      call acclyr(jcalc13,ocetra(1,1,1,icalc13),pddpo,1)
      call acclyr(jphyto13,ocetra(1,1,1,iphy13),pddpo,1)   
      call acclyr(jgrazer13,ocetra(1,1,1,izoo13),pddpo,1)  
#endif 
#ifdef AGG
      call acclyr(jnos,ocetra(1,1,1,inos),pddpo,1)      
#endif     
#ifdef CFC
      call acclyr(jcfc11,ocetra(1,1,1,icfc11),pddpo,1)
      call acclyr(jcfc12,ocetra(1,1,1,icfc12),pddpo,1)
      call acclyr(jsf6,ocetra(1,1,1,isf6),pddpo,1)
#endif
      call acclyr(jprefo2,ocetra(1,1,1,iprefo2),pddpo,1)
      call acclyr(jprefpo4,ocetra(1,1,1,iprefpo4),pddpo,1)
      call acclyr(jprefalk,ocetra(1,1,1,iprefalk),pddpo,1)
      call acclyr(jprefdic,ocetra(1,1,1,iprefdic),pddpo,1)
      call acclyr(jdicsat,ocetra(1,1,1,idicsat),pddpo,1)

!     Accumulate level diagnostics
      IF (SUM(jlvlphyto+jlvlgrazer+jlvlphosph+jlvloxygen+jlvliron+      &
     &  jlvlano3+jlvlalkali+jlvlsilica+jlvldic+jlvldoc+jlvlpoc+jlvlcalc+&
     &  jlvlopal+jlvln2o+jlvlco3+jlvlph+jlvlomegac++jlvldic13+jlvldic14+&
     &  jlvlnos+jlvld13c+jlvlpoc13+jlvldoc13+jlvlcalc13+jlvlphyto13+    &
     &  jlvlgrazer13+jlvlnatdic+jlvlnatalkali+jlvlnatcalc+jlvlnatco3+   &
     &  jlvlnatomegac+jlvlnos+                                          &
     &  jlvlcfc11+jlvlcfc12+jlvlsf6+jlvlprefo2+jlvlprefpo4+jlvlprefalk+ &
     &  jlvlprefdic+jlvldicsat).NE.0) THEN
        DO k=1,kpke
          call bgczlv(pddpo,k,ind1,ind2,wghts)
          call acclvl(jlvlphyto,ocetra(1,1,1,iphy),k,ind1,ind2,wghts)
          call acclvl(jlvlgrazer,ocetra(1,1,1,izoo),k,ind1,ind2,wghts)
          call acclvl(jlvlphosph,ocetra(1,1,1,iphosph),k,ind1,ind2,wghts)
          call acclvl(jlvloxygen,ocetra(1,1,1,ioxygen),k,ind1,ind2,wghts)
          call acclvl(jlvliron,ocetra(1,1,1,iiron),k,ind1,ind2,wghts)
          call acclvl(jlvlano3,ocetra(1,1,1,iano3),k,ind1,ind2,wghts)
          call acclvl(jlvlalkali,ocetra(1,1,1,ialkali),k,ind1,ind2,wghts)
          call acclvl(jlvlsilica,ocetra(1,1,1,isilica),k,ind1,ind2,wghts)
          call acclvl(jlvldic,ocetra(1,1,1,isco212),k,ind1,ind2,wghts)
          call acclvl(jlvldoc,ocetra(1,1,1,idoc),k,ind1,ind2,wghts)
          call acclvl(jlvlpoc,ocetra(1,1,1,idet),k,ind1,ind2,wghts)
          call acclvl(jlvlcalc,ocetra(1,1,1,icalc),k,ind1,ind2,wghts)
          call acclvl(jlvlopal,ocetra(1,1,1,iopal),k,ind1,ind2,wghts)
          call acclvl(jlvln2o,ocetra(1,1,1,ian2o),k,ind1,ind2,wghts)          
          call acclvl(jlvlco3,co3,k,ind1,ind2,wghts)
          call acclvl(jlvlph,hi,k,ind1,ind2,wghts)
          call acclvl(jlvlomegac,OmegaC,k,ind1,ind2,wghts)
#ifdef natDIC
          call acclvl(jlvlnatdic,ocetra(1,1,1,inatsco212),k,ind1,ind2,wghts)
          call acclvl(jlvlnatalkali,ocetra(1,1,1,inatalkali),k,ind1,ind2,wghts)
          call acclvl(jlvlnatcalc,ocetra(1,1,1,inatcalc),k,ind1,ind2,wghts)
          call acclvl(jlvlnatco3,natco3,k,ind1,ind2,wghts)
          call acclvl(jlvlnatomegac,natOmegaC,k,ind1,ind2,wghts)
#endif
#ifdef cisonew
          call acclvl(jlvld13c,d13c,k,ind1,ind2,wghts)
          call acclvl(jlvldic13,ocetra(1,1,1,isco213),k,ind1,ind2,wghts)
          call acclvl(jlvldic14,ocetra(1,1,1,isco214),k,ind1,ind2,wghts)
          call acclvl(jlvlpoc13,ocetra(1,1,1,idet13),k,ind1,ind2,wghts)
          call acclvl(jlvldoc13,ocetra(1,1,1,idoc13),k,ind1,ind2,wghts)
          call acclvl(jlvlcalc13,ocetra(1,1,1,icalc13),k,ind1,ind2,wghts)
          call acclvl(jlvlphyto13,ocetra(1,1,1,iphy13),k,ind1,ind2,wghts)
          call acclvl(jlvlgrazer13,ocetra(1,1,1,izoo13),k,ind1,ind2,wghts)
#endif
#ifdef AGG
          call acclvl(jlvlnos,ocetra(1,1,1,inos),k,ind1,ind2,wghts)
#endif     
#ifdef CFC
          call acclvl(jlvlcfc11,ocetra(1,1,1,icfc11),k,ind1,ind2,wghts)
          call acclvl(jlvlcfc12,ocetra(1,1,1,icfc12),k,ind1,ind2,wghts)
          call acclvl(jlvlsf6,ocetra(1,1,1,isf6),k,ind1,ind2,wghts)
#endif
          call acclvl(jlvlprefo2,ocetra(1,1,1,iprefo2),k,ind1,ind2,wghts)
          call acclvl(jlvlprefpo4,ocetra(1,1,1,iprefpo4),k,ind1,ind2,wghts)
          call acclvl(jlvlprefalk,ocetra(1,1,1,iprefalk),k,ind1,ind2,wghts)
          call acclvl(jlvlprefdic,ocetra(1,1,1,iprefdic),k,ind1,ind2,wghts)
          call acclvl(jlvldicsat,ocetra(1,1,1,idicsat),k,ind1,ind2,wghts)
        ENDDO
      ENDIF

#ifndef sedbypass
!     Accumulate sediments
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

!     Accumulate sediment burial
      call accbur(jburssso12,burial(1,1,issso12))
      call accbur(jburssssil,burial(1,1,issssil))
      call accbur(jbursssc12,burial(1,1,isssc12))
      call accbur(jburssster,burial(1,1,issster))
#endif

#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after BGC: call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif	 

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
