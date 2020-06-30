      SUBROUTINE BELEG_VARS(kpaufr,kpie,kpje,kpke,kbnd,pddpo,prho,omask,      &
                            pglon,pglat)
!******************************************************************************
!
! BELEG_VARS - initialize bgc variables.
!
!  Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!  Modified
!  --------
!  J.Schwinger,        *NORCE Climate, Bergen*    2020-05-19
!   -split the original BELEG_BGC in two parts, BELEG_PARM and BELEG_VARS
!
!
!  Purpose
!  -------
!  - set initial values for bgc variables.
!
!
!  Parameter list:
!  ---------------
!     *INTEGER*   *kpaufr*  - 1/0 flag, 1 indicating a restart run
!     *INTEGER*   *kpie*    - 1st dimension of model grid.
!     *INTEGER*   *kpje*    - 2nd dimension of model grid.
!     *INTEGER*   *kpke*    - 3rd (vertical) dimension of model grid.
!     *INTEGER*   *kbnd*    - nb of halo grid points
!     *REAL*      *pddpo*   - size of grid cell (3rd dimension) [m].
!     *REAL*      *prho*    - density [g/cm^3].
!     *REAL*      *omask*   - ocean mask.
!     *REAL*      *pglon*   - longitude of grid cell [deg].
!     *REAL*      *pglat*   - latitude  of grid cell [deg].
!
!
!******************************************************************************
      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      USE mo_control_bgc, only: rmasks
      use mo_param1_bgc
      use mo_vgrid, only: kmle,kbo
      USE mod_xc,   only: mnproc

      implicit none      

      INTEGER, intent(in) :: kpaufr,kpie,kpje,kpke,kbnd
      REAL,    intent(in) :: pddpo(kpie,kpje,kpke)
      REAL,    intent(in) :: prho (kpie,kpje,kpke)
      REAL,    intent(in) :: omask(kpie,kpje)
      REAL,    intent(in) :: pglon(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)
      REAL,    intent(in) :: pglat(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd)

      ! local variables
      INTEGER :: i,j,k,l


#ifdef FB_BGC_OCE
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
        abs_oce(i,j,k)=1.
      ENDDO
      ENDDO
      ENDDO
#endif

!
! Initialisation of ocean tracers and sediment
!

! Initialise ocean tracers with WOA and GLODAP data. This is done even in case
! of a restart since some tracers (e.g. C-isotopes) might not be in the restart 
! file and aufr.f90 instead expects an initialised field.
      call profile_gd(kpie,kpje,kpke,kbnd,pglon,pglat,omask)

! If this is a restart run initialisation is done in aufr.F90 
      IF(kpaufr.EQ.1) RETURN

      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
        IF (omask(i,j) .GT. 0.5 ) THEN
          ! convert WOA tracers kmol/m^3 -> mol/kg; GLODAP dic and alk
          ! are already in mol/kg. We need these units here, since after 
          ! initialisation the tracer field is passed to the ocean model
          ! first where units are mol/kg.
          ocetra(i,j,k,iphosph) = ocetra(i,j,k,iphosph)/prho(i,j,k)
          ocetra(i,j,k,ioxygen) = ocetra(i,j,k,ioxygen)/prho(i,j,k)
          ocetra(i,j,k,iano3)   = ocetra(i,j,k,iano3)  /prho(i,j,k)
          ocetra(i,j,k,isilica) = ocetra(i,j,k,isilica)/prho(i,j,k)
#ifdef cisonew
          ! d13C based on Eide data is read in above (profile_gd)                        
          ! Convert to 13C using model initial (ie GLODAP) total C
          ! If restarting, this is redone with model total C from restart in aufr_bgc.F90 
          beta13=ocetra(i,j,k,isco213)/1000.+1.
          ocetra(i,j,k,isco213) = ocetra(i,j,k,isco212)*beta13*re1312/(1.+beta13*re1312)

          ! 14C is read in as small delta14C (calculated from R. Key, 2003 and Eide et al. 2017)
          ! Convert to 14C using model total C, and normalize by c14fac to prevent numerical errors
          beta14=ocetra(i,j,k,isco214)/1000.+1.
          ocetra(i,j,k,isco214) = ocetra(i,j,k,isco212)*beta14*re14to/c14fac
#endif
        ENDIF
      ENDDO
      ENDDO
      ENDDO

! Initialise remaining ocean tracers 
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
        IF(omask(i,j) .GT. 0.5) THEN
          ocetra(i,j,k,igasnit)=1.e-10
          ocetra(i,j,k,idoc)   =1.e-8
          ocetra(i,j,k,iphy)   =1.e-8 
          ocetra(i,j,k,izoo)   =1.e-8 
          ocetra(i,j,k,idet)   =1.e-8 
          ocetra(i,j,k,icalc)  =0. 
          ocetra(i,j,k,iopal)  =1.e-8 
          ocetra(i,j,k,ian2o)  =0. 
          ocetra(i,j,k,idms)   =0. 
          ocetra(i,j,k,ifdust) =0. 
          ocetra(i,j,k,iiron)  =fesoly
          ocetra(i,j,k,iprefo2)=0.
          ocetra(i,j,k,iprefpo4)=0.
          ocetra(i,j,k,iprefalk)=0.
          ocetra(i,j,k,iprefdic)=0.
          ocetra(i,j,k,idicsat)=1.e-8
          hi(i,j,k)            =1.e-8
          co3(i,j,k)           =0.
          co2star(i,j,k)       =20.e-6
#ifdef AGG
! calculate initial numbers from mass, to start with appropriate size distribution
          snow = (ocetra(i,j,k,iphy)+ocetra(i,j,k,idet))*1.e+6
          ocetra(i,j,k,inos)   = snow / cellmass / (FractDim+1.)
          ocetra(i,j,k,iadust) =0. 
#endif /*AGG*/
#ifdef CFC
          ocetra(i,j,k,icfc11)   =0.
          ocetra(i,j,k,icfc12)   =0.
          ocetra(i,j,k,isf6)     =0.
#endif
#ifdef natDIC
          nathi(i,j,k)           =1.e-8
          natco3(i,j,k)          =0.
          ocetra(i,j,k,inatcalc) =0. 
#endif
#ifdef cisonew
          rco213=ocetra(i,j,k,isco213)/(ocetra(i,j,k,isco212)+safediv)
          rco214=ocetra(i,j,k,isco214)/(ocetra(i,j,k,isco212)+safediv)
          ocetra(i,j,k,iphy13) =ocetra(i,j,k,iphy)*rco213*bifr13
          ocetra(i,j,k,iphy14) =ocetra(i,j,k,iphy)*rco214*bifr14
          ocetra(i,j,k,izoo13) =ocetra(i,j,k,izoo)*rco213*bifr13
          ocetra(i,j,k,izoo14) =ocetra(i,j,k,izoo)*rco214*bifr14
          ocetra(i,j,k,idoc13) =ocetra(i,j,k,idoc)*rco213*bifr13
          ocetra(i,j,k,idoc14) =ocetra(i,j,k,idoc)*rco214*bifr14
          ocetra(i,j,k,idet13) =ocetra(i,j,k,idet)*rco213*bifr13
          ocetra(i,j,k,idet14) =ocetra(i,j,k,idet)*rco214*bifr14
          ocetra(i,j,k,icalc13)=ocetra(i,j,k,icalc)*rco213
          ocetra(i,j,k,icalc14)=ocetra(i,j,k,icalc)*rco214
#endif
        ENDIF ! omask > 0.5
      ENDDO
      ENDDO
      ENDDO

! Initialise preformed tracers in the mixed layer; note that the 
! whole field has been initialised to zero above
      DO k=1,kmle
      DO j=1,kpje
      DO i=1,kpie
        IF(omask(i,j) .GT. 0.5) THEN
          ocetra(i,j,k,iprefo2) =ocetra(i,j,k,ioxygen)
          ocetra(i,j,k,iprefpo4)=ocetra(i,j,k,iphosph)
          ocetra(i,j,k,iprefalk)=ocetra(i,j,k,ialkali)
          ocetra(i,j,k,iprefdic)=ocetra(i,j,k,isco212)
        ENDIF
      ENDDO
      ENDDO
      ENDDO


! Initial values for sediment
#ifndef sedbypass
      DO  k=1,ks
      DO  j=1,kpje
      DO  i=1,kpie 
        IF(omask(i,j) .GT. 0.5) THEN
          powtra(i,j,k,ipowaic)=ocetra(i,j,kbo(i,j),isco212)
          powtra(i,j,k,ipowaal)=ocetra(i,j,kbo(i,j),ialkali)
          powtra(i,j,k,ipowaph)=ocetra(i,j,kbo(i,j),iphosph)
          powtra(i,j,k,ipowaox)=ocetra(i,j,kbo(i,j),ioxygen)
          powtra(i,j,k,ipown2) =0.
          powtra(i,j,k,ipowno3)=ocetra(i,j,kbo(i,j),iano3)
          powtra(i,j,k,ipowasi)=ocetra(i,j,kbo(i,j),isilica)      
          sedlay(i,j,k,issso12)=1.e-8
          sedlay(i,j,k,isssc12)=1.e-8
          sedlay(i,j,k,issster)=30.
          sedlay(i,j,k,issssil)=1.e-8
          sedhpl(i,j,k)        =hi(i,j,kbo(i,j))
#ifdef cisonew
          rco213=ocetra(i,j,kbo(i,j),isco213)/(ocetra(i,j,kbo(i,j),isco212)+safediv)
          rco214=ocetra(i,j,kbo(i,j),isco214)/(ocetra(i,j,kbo(i,j),isco212)+safediv)
          powtra(i,j,k,ipowc13)=powtra(i,j,k,ipowaic)*rco213*bifr13
	  powtra(i,j,k,ipowc14)=powtra(i,j,k,ipowaic)*rco214*bifr14
          sedlay(i,j,k,issso13)=sedlay(i,j,k,issso12)*rco213*bifr13
          sedlay(i,j,k,issso14)=sedlay(i,j,k,issso12)*rco214*bifr14
          sedlay(i,j,k,isssc13)=sedlay(i,j,k,isssc12)*rco213
          sedlay(i,j,k,isssc14)=sedlay(i,j,k,isssc12)*rco214
#endif
        ELSE
          powtra(i,j,k,ipowno3)=rmasks
          powtra(i,j,k,ipown2) =rmasks
          powtra(i,j,k,ipowaic)=rmasks
          powtra(i,j,k,ipowaal)=rmasks
          powtra(i,j,k,ipowaph)=rmasks
          powtra(i,j,k,ipowaox)=rmasks
          powtra(i,j,k,ipowasi)=rmasks
          sedlay(i,j,k,issso12)=rmasks
          sedlay(i,j,k,isssc12)=rmasks
          sedlay(i,j,k,issssil)=rmasks
          sedlay(i,j,k,issster)=rmasks
          sedlay(i,j,k,issssil)=rmasks
          sedhpl(i,j,k)        =rmasks
#ifdef cisonew
          powtra(i,j,k,ipowc13)=rmasks
	  powtra(i,j,k,ipowc14)=rmasks
          sedlay(i,j,k,issso13)=rmasks
          sedlay(i,j,k,issso14)=rmasks
          sedlay(i,j,k,isssc13)=rmasks
          sedlay(i,j,k,isssc14)=rmasks
#endif
        ENDIF
      ENDDO
      ENDDO
      ENDDO

      ! last and final sediment layer
      DO  l=1,nsedtra
      DO  j=1,kpje
      DO  i=1,kpie
         burial(i,j,l)=0.
      ENDDO
      ENDDO
      ENDDO
#endif

      return
!******************************************************************************
      end subroutine beleg_vars


