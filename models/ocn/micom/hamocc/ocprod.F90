      SUBROUTINE OCPROD(kpie,kpje,kpke,ptho,pddpo,                     &
     &               pdlxp,pdlyp,pdpio,ptiestu,ptiestw,kplmon,omask)

!$Source: /scratch/local1/m212047/patrick/SRC_MPI/src_hamocc/RCS/ocprod.f90,v $\\
!$Revision: 1.1 $\\
!$Date: 2005/01/28 08:37:45 $\\

!**********************************************************************
!
!**** *OCPROD* - .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!
!     Purpose
!     -------
!     compute biological production, settling of debris, and related biogeochemistry
!
!     Method:
!     ------
!     kchck=1 can be used to check max/min of bgc arrays on wet/dry cells.
!     Note: prosil is used only for k=1,2. It is adressed, however, for
!           k=1,4 in loop 100. To save memory,  ???
!
!     _ant fields are natural PLUS anthropogenic (not anthropogenic only!!!)
!
!     *CALL*       *OCPROD*
!
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *ptho*    - potential temperature [deg C].
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st dimension) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (2nd dimension) [m].
!     *REAL*    *pdpio*   - inverse size of grid cell (3rd dimension)[m].
!
!     Externals
!     ---------
!     .
!**********************************************************************
! nocetra is the number of all BGC element (size of ocetra(,,,l))

!      USE mo_timeser_bgc
      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      use mo_param1_bgc 

      USE mo_control_bgc
      USE mo_bgcmean
      use mod_xc , only : isp,ifp,ilp

      implicit none

      INTEGER :: kplmon,kpie,kpje,kpke
      INTEGER :: i,j,k,l,kab1
      INTEGER :: kinv,kdonor,found
      REAL :: ptho (kpie,kpje,kpke)
      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: pdpio(kpie,kpje,kpke)
      REAL :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL :: ptiestu(kpie,kpje,kpke+1)
      REAL :: ptiestw(kpie,kpje,kpke+1)
      REAL :: abs_bgc(kpie,kpje,kpke)
      REAL :: omask(kpie,kpje)
      
      REAL :: dmsp1,dmsp2,dmsp3,dmsp4,dmsp5,dmsp6
      REAL :: atten,avphy,avanut,avanfe,pho,xa,xn,ya,yn,phosy,         &
     &        volcell,avgra,grazing,avsil,graton,                      &
     &        gratpoc,grawa,bacfra,phymor,zoomor,excdoc,exud,          &
     &        export, delsil, delcar, sterph, sterzo, remin,           &
     &        docrem, opalrem, remin2o, aou,refra,pocrem,phyrem
      
      REAL :: zoothresh,phythresh
      REAL :: temfa,phofa                  ! temperature and irradiation factor for photosynthesis
      REAL :: dustinp, gutscale
      REAL :: fopa, fdet, fcal
      REAL :: absorption
      REAL :: dmsprod,dms_bac,dms_uv 
      REAL :: bdp,dtr,dp_ez 
      REAL :: detref, detrl
!#ifdef __c_isotopes
      REAL :: rem13,rem14
      REAL :: rl13, rl14
      REAL :: rocean13, rocean14, flui13, flui14
      REAL :: fcal13,fcal14
      REAL :: d14C
!ib 
      REAL, DIMENSION(kpie,kpje,kpke) :: d13C, dd14C
!ib
!#endif
#ifdef AGG
      REAL :: fphy
      REAL :: wmass(kpie,kpje,kpke)
      REAL :: wnumb(kpie,kpje,kpke)
      REAL :: aggregate(kpie,kpje,kpke)
      REAL :: dustagg(kpie,kpje,kpke)

      REAL :: avmass, avsize, avsmin1, avsmax1, avnos, anosloss     
      INTEGER :: nosin1, nosde1, nosin2, nosde2, nosin3, nosde3
      REAL :: zmornos, avm, avn, eps, e1,e2,e3,e4,es1,es3
      REAL :: TopM, TopF, snow,fshear,sagg1,sagg2,sagg4
      REAL :: sett_agg,shear_agg,effsti,dfirst,dshagg,dsett
      REAL :: checksize,nacheck,flar1,flar2,flar3
      REAL :: fTSFac,fTMFac,fTopF,fTopM,wphy,wphyup,wnosup,wnos
#endif 
!ka      REAL :: psum0,psum1,pinv0,pinv1,dpinv
!ib
      INTEGER, DIMENSION(kpie,kpje)   :: ind1,ind2
      REAL, DIMENSION(kpie,kpje,ddm)  :: wghts
      REAL, DIMENSION(kpie,kpje)      :: aux2d_dmsprod 
      REAL, DIMENSION(kpie,kpje)      :: aux2d_dms_bac
      REAL, DIMENSION(kpie,kpje)      :: aux2d_dms_uv 
      REAL, DIMENSION(kpie,kpje)      :: aux2d_export 
      REAL, DIMENSION(kpie,kpje)      :: aux2d_expoca 
      REAL, DIMENSION(kpie,kpje)      :: aux2d_exposi 
      REAL, DIMENSION(kpie,kpje,kpke) :: aux3d_phosy
 
!kma
      aux2d_dmsprod(:,:)=0. 
      aux2d_dms_bac(:,:)=0. 
      aux2d_dms_uv (:,:)=0.  
      aux2d_export (:,:)=0.  
      aux2d_expoca (:,:)=0.  
      aux2d_exposi (:,:)=0.  
      aux3d_phosy  (:,:,:)=0.

!ib
!
! Constant parameters
!
! parameter definition in BELEG_BGC.F

      dmsp6=dmspar(6)
      dmsp5=dmspar(5)
      dmsp4=dmspar(4)
      dmsp3=dmspar(3)
      dmsp2=dmspar(2)
      dmsp1=dmspar(1)

! Calculate swr absorption by water and phytoplankton

! Almost half of the SWR is absorbed in the surface layer. --> *0.6
! Implicit, it's faster then: abs_bgc(i,j,k)=abs_bgc(i,j,k-1)*exp(-atten*pddpo(i,j,k-1))
#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'beginning of OCRPOD '
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask &
     &                  ,0)
!      call INVENTORY_BGC(kpie,kpje,kpke)
#endif   
!$OMP PARALLEL DO PRIVATE(absorption,atten,kab1)
      DO j=1,kpje
      DO i=1,kpie
        IF(omask(i,j).GT.0.5) THEN
         atten=atten_w+atten_c*max(0.,ocetra(i,j,1,iphy))
         absorption=1.
!JT
         abs_bgc(i,j,1)=absorption
!         abs_bgc(i,j,1)=((absorption/atten)*                            &
!     &                  (1.-exp(-atten*MIN(pddpo(i,j,1),100.))))        &
!     &                  /MIN(pddpo(i,j,1),100.)
         kab1=1
#ifdef FB_BGC_OCE     
         abs_oce(i,j,1)=1.
#endif 	
      DO k=2,kpke

         abs_bgc(i,j,k)=0.
#ifdef FB_BGC_OCE     
         abs_oce(i,j,k)=0.
#endif 	

        IF(pddpo(i,j,k).gt.0.0) THEN
          atten=atten_w+atten_c*max(0.,ocetra(i,j,kab1,iphy))
          absorption=abs_bgc(i,j,kab1)
          abs_bgc(i,j,k)=((absorption/atten)*                           &
     &                   (1.-exp(-atten*pddpo(i,j,k))))/pddpo(i,j,k)
#ifdef FB_BGC_OCE
          abs_oce(i,j,k)=abs_oce(i,j,kab1)*absorption
          if (k.eq.2) then
          abs_oce(i,j,2)=atten_f*absorption
          endif
#endif 

          kab1=k

        ENDIF
      ENDDO

        ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

! dust flux from the atmosphere to the surface layer; dust fields are
! monthly mean values (kg/m2/year)
! dissolved iron is a fixed fraction (typically 3.5%), and immediately released

!$OMP PARALLEL DO PRIVATE(dustinp)
      do j=1,kpje
      do i=1,kpie
       if(omask(i,j).gt.0.5) then
        dustinp=dusty(i,j,kplmon)/365.*dtb*pdpio(i,j,1)
        ocetra(i,j,1,ifdust)=ocetra(i,j,1,ifdust)+dustinp 
        ocetra(i,j,1,iiron)=ocetra(i,j,1,iiron)+dustinp*perc_diron 
!ka flux of nutrients, C & alkalinity from land into surface layer
!ka constant flux calculated in beleg_bgc 
!        ocetra(i,j,1,isco212)=                                         &
!     &  ocetra(i,j,1,isco212)+sco212_sfc*dtb*pdpio(i,j,1)
!        ocetra(i,j,1,ialkali)=                                         &
!     &  ocetra(i,j,1,ialkali)+alkali_sfc*dtb*pdpio(i,j,1)
!        ocetra(i,j,1,iphosph)=                                         &
!     &  ocetra(i,j,1,iphosph)+phosph_sfc*dtb*pdpio(i,j,1)
!        ocetra(i,j,1,iano3)  =                                         &
!     &  ocetra(i,j,1,iano3)  +ano3_sfc  *dtb*pdpio(i,j,1)
!        ocetra(i,j,1,isilica)=                                         &
!     &  ocetra(i,j,1,isilica)+silica_sfc*dtb*pdpio(i,j,1)
       endif      
      enddo
      enddo
!$OMP END PARALLEL DO

!
! averaging monthly stocks and flows
! 

!ib
! accumulate layer diagnostics
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
       call acclyr(jco3,co3,pddpo,1)                      
       call acclyr(jph,hi,pddpo,1)         
#ifdef __c_isotopes
!$OMP PARALLEL DO PRIVATE(d14C) 
      DO j=1,kpje
        DO i=1,kpie
          IF(omask(i,j).GT.0.5) THEN
            DO k=1,kpke
! delta notation for 13C, 14C (js 6.3.2006)
              d13C(i,j,k)  = ((ocetra(i,j,k,isco213)/                   &
     &                       ocetra(i,j,k,isco212)) -1.) *1000.
              d14C         = ((ocetra(i,j,k,isco214)/                   &
     &                       ocetra(i,j,k,isco212)) -1.) *1000.
              dd14C(i,j,k) = d14C - 2.* (d13C(i,j,k) + 25.) *           &
     &                       (1.+ d14C *1.e-3)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
       call acclyr(jdic13,d13C,pddpo,1)                 
       call acclyr(jdic14,dd14C,pddpo,1)                
#endif     
#ifdef AGG
       call acclyr(jnos,ocetra(1,1,1,inos),pddpo,1)      
#endif     

! accumulate level diagnostics
      IF (SUM(jlvlphyto+jlvlgrazer+jlvlphosph+jlvloxygen+jlvliron+      &
     &  jlvlano3+jlvlalkali+jlvlsilica+jlvldic+jlvldoc+jlvlpoc+jlvlcalc+&
     &  jlvlco3+jlvlph+jlvldic13+jlvldic14+jlvlnos).NE.0) &
     &  THEN
        DO k=1,kpke
          call bgczlv(pddpo,k,ind1,ind2,wghts)
          call acclvl(jlvlphyto,ocetra(1,1,1,iphy),k,ind1,ind2,wghts)
          call acclvl(jlvlgrazer,ocetra(1,1,1,izoo),k,ind1,ind2,wghts)
          call acclvl(jlvlphosph,ocetra(1,1,1,iphosph),k,ind1,ind2,     &
     &      wghts)
          call acclvl(jlvloxygen,ocetra(1,1,1,ioxygen),k,ind1,ind2,     &
     &      wghts)
          call acclvl(jlvliron,ocetra(1,1,1,iiron),k,ind1,ind2,wghts)
          call acclvl(jlvlano3,ocetra(1,1,1,iano3),k,ind1,ind2,wghts)
          call acclvl(jlvlalkali,ocetra(1,1,1,ialkali),k,ind1,ind2,     &
     &      wghts)
          call acclvl(jlvlsilica,ocetra(1,1,1,isilica),k,ind1,ind2,     &
     &      wghts)
          call acclvl(jlvldic,ocetra(1,1,1,isco212),k,ind1,ind2,wghts)
          call acclvl(jlvldoc,ocetra(1,1,1,idoc),k,ind1,ind2,wghts)
          call acclvl(jlvlpoc,ocetra(1,1,1,idet),k,ind1,ind2,wghts)
          call acclvl(jlvlcalc,ocetra(1,1,1,icalc),k,ind1,ind2,wghts)
          call acclvl(jlvlopal,ocetra(1,1,1,iopal),k,ind1,ind2,wghts)
          call acclvl(jlvlco3,co3,k,ind1,ind2,wghts)
          call acclvl(jlvlph,hi,k,ind1,ind2,wghts)
#ifdef __c_isotopes
          call acclvl(jlvldic13,d13C,k,ind1,ind2,wghts)
          call acclvl(jlvldic14,dd14C,k,ind1,ind2,wghts)
#endif     
#ifdef AGG
          call acclvl(jlvlnos,ocetra(1,1,1,inos),k,ind1,ind2,wghts)
#endif     
        ENDDO
      ENDIF


!
! Sampling timeseries-1 : global inventory
!

!
! Sampling timeseries-1 : concentrations at specific positions
!
 
#ifdef AGG
!***********************************************************************
!  Have a special resetting for numbers, that fixes their conc. to one
!  depending on mass of marine snow:
!  Compartments have already been set to 0 in 
!  ADVECTION_BGC.h and OCTDIFF_BGC.h.
!  Ensure that if there is no mass, there are no particles, and 
!  that the number of particles is in the right range (this is crude, but
!  is supposed to happen only due to numerical errors such as truncation or 
!  overshoots during advection)
! (1) avnos<<avmass, such that eps = FractDim + 1: increase numbers
!     such that eps = FractDim + 1 + safe (currently set to 1.e-6 in BELEG_BGC) 
! (2) avnos>>avmass, such that  Nbar (=Mass/Nos/cellmass) <=1: decrease numbers
!     such that Nbar=1.1 (i.e. 1.1 cells per aggregate, set in BELEG_BGC) 
!************************************************************************

       DO  k=1,kpke
         DO j=1,kpje
           DO i=1,kpie
            IF(pddpo(i,j,k).GT.1.e-12.and.omask(i,j).gt.0.5) THEN
             avmass = ocetra(i,j,k,iphy) + ocetra(i,j,k,idet)
             snow = avmass*1.e+6
! look for max. and min average size = Nbar             
             if(ocetra(i,j,k,inos).gt.0.) then
                avsize=snow/ocetra(i,j,k,inos)/cellmass
                avsmin1=MIN(avsize,avsmin1)
                avsmax1=MAX(avsize,avsmax1)
             endif
! check whether the numbers had to be decreased or increased     
             if (snow*pupper.gt.ocetra(i,j,k,inos)) then
               nosin1 = nosin1 + 1
             endif
             if (snow*plower.lt.ocetra(i,j,k,inos)) then
               nosde1 = nosde1 + 1
             endif
             ocetra(i,j,k,inos) = MAX(snow*pupper,ocetra(i,j,k,inos)) 
             ocetra(i,j,k,inos) = MIN(snow*plower,ocetra(i,j,k,inos)) 
            ENDIF
           ENDDO
       ENDDO
      ENDDO
         
#endif  /*AGG*/

#define bioprod
#ifdef bioprod
!      pinv0=0.0
!      DO j=1,kpje
!      DO i=1,kpie
!        IF(omask(i,j).gt.0.0) THEN
!      DO k=1,kpke
!      dpinv=                                                      &
!     &   ocetra(i,j,k,idet)+ocetra(i,j,k,idoc)+ocetra(i,j,k,iphy) &
!     &  +ocetra(i,j,k,izoo)+ocetra(i,j,k,iphosph)                  
!      pinv0=pinv0+dpinv*pddpo(i,j,k)*pdlxp(i,j)*pdlyp(i,j)      
!      ENDDO
!        ENDIF
!      ENDDO
!      ENDDO
!
! Biological productivity reaches down to about 80-100m
!
!
! first calculate which layer lie in the euphotic zone
! (at the moment the layer has to lie completely in it)

      call calc_ez(kpie,kpje,kpke,pddpo,ptiestw)

! then let it grow... 

!$OMP PARALLEL DO                            &                
!$OMP&PRIVATE(avphy,avgra,avsil,avanut,avanfe,pho,xa,xn,phosy,  &
!$OMP&        ya,yn,grazing,graton,gratpoc,grawa,bacfra,phymor, &  
!$OMP&        zoomor,excdoc,exud,export,delsil,delcar,dmsprod,  &
!$OMP&        dms_bac,dms_uv,bdp,dtr,rocean13,rocean14,flui13,  &
!$OMP&        flui14,dp_ez,phofa,temfa,zoothresh)
      DO 1 j=1,kpje
      DO 1 i=1,kpie

      DO 100 K=1,kwrbioz(i,j)

      IF(pddpo(i,j,k).GT.1.e-12.and.omask(i,j).gt.0.5) THEN

!      psum0=                                                          &
!     &   ocetra(i,j,k,idet)+ocetra(i,j,k,idoc)+ocetra(i,j,k,iphy)     &
!     &  +ocetra(i,j,k,izoo)+ocetra(i,j,k,iphosph)  

! depth of euphotic zone
      dp_ez=100.0
      if(ptiestw(i,j,k+1).gt.dp_ez) then
        bdp=abs(ptiestw(i,j,k)-dp_ez)
        if(bdp.gt.pddpo(i,j,k)) then
          write(io_stdo_bgc,*) 'bdp gt pddpo(i,j,k) 1',i,j,k
        endif
      else
        bdp=pddpo(i,j,k)
      endif

#ifdef AGG
         avmass = ocetra(i,j,k,iphy)+ocetra(i,j,k,idet)
#endif /*AGG*/

#ifdef __c_isotopes
         rocean13=ocetra(i,j,k,isco213)/ocetra(i,j,k,isco212)                ! "ratio oceanic 13C/12C" (at least in carchm)
         rocean14=ocetra(i,j,k,isco214)/ocetra(i,j,k,isco212)                ! ratio oceanic 14C/12C
#endif

         phofa=pi_alpha*strahl(i,j)*abs_bgc(i,j,k) 
         temfa= 0.6* 1.066**ptho(i,j,k)                 
!taylor:      temfa= 0.6*(1. + 0.0639*ptho(i,j,k) *               &
!    &               (1. + 0.0639*ptho(i,j,k)/2. * (1. + 0.0639*ptho(i,j,k)/3.)))
         pho= dtb* phofa*temfa/sqrt(phofa**2 + temfa**2)     

         avphy=MAX(0.,ocetra(i,j,k,iphy))  
         avgra=MAX(0.,ocetra(i,j,k,izoo))  

!ka phytos & zoos present in top 100m only
!      if(ptiestw(i,j,k+1).gt.dp_ez) then
!        avphy=avphy*(pddpo(i,j,k)/bdp)
!        avgra=avgra*(pddpo(i,j,k)/bdp)
!      ENDIF

         avsil=MAX(0.,ocetra(i,j,k,isilica))
         avanut=MAX(0.,MIN(ocetra(i,j,k,iphosph),                      &
     &          rnoi*ocetra(i,j,k,iano3)))
         avanfe=MAX(0.,MIN(avanut,ocetra(i,j,k,iiron)/riron))
         xa=avanfe
         xn=xa/(1.+pho*avphy/(xa+bkphy)) 
         phosy=MAX(0.,xa-xn)
         xn=MAX(xn,1.e-10)
         ya=avphy+phosy
         yn=(ya+grazra*avgra*phytomi/(avphy+bkzoo))                    &
     &            /(1.+grazra*avgra/(avphy+bkzoo))
         grazing=MAX(0.,ya-yn)
         graton=epsher*(1.-zinges)*grazing
         gratpoc=(1.-epsher)*grazing ! epsher=0.8
         grawa=epsher*zinges*grazing
         bacfra=remido*ocetra(i,j,k,idoc)
         phymor=dyphy*MAX(0.,(ocetra(i,j,k,iphy)-2.*phytomi)) ! dyphy=.008*dt
!new iris
         zoothresh=MAX(0.,(ocetra(i,j,k,izoo)-2.*grami))                  ! quadratic mortality
         zoomor=spemor*zoothresh*zoothresh                                ! *10 compared to linear in tropics (tinka)
         excdoc=gammaz*zoothresh                                          ! excretion of doc by zooplankton
         exud=gammap*MAX(0.,(ocetra(i,j,k,iphy)-2.*phytomi))

         export= zoomor*(1.-ecan) + phymor + gratpoc ! ecan=.95, gratpoc= .2*grazing

! new from emr version 4.5.06
#ifdef __c_isotopes
         flui13=max(rocean13-1.,0.)        ! assumes rocean >1 in euphotic layer 
         flui14=max(rocean14-1.,0.)        ! flui means what? (flui14 never used) flux...?
#endif /*__c_isotopes*/

#ifdef AGG	 
         delsil=MIN(ropal*phosy*avsil/(avsil+bkopal),0.5*avsil) 
	 delcar=rcalc*MIN(calmax*phosy,(phosy-delsil/ropal))
#else
         delsil=MIN(ropal*export*avsil/(avsil+bkopal),0.8*avsil) 
!         delsil=MIN(ropal*export*avsil/(avsil+bkopal),0.5*avsil) 
         delcar=rcalc * export * bkopal/(avsil+bkopal)
#endif

         dmsprod = (dmsp5*delsil+dmsp4*delcar)                        &
     &            *(1.+1./(ptho(i,j,k)+dmsp1)**2)  	 
	 dms_bac = dmsp3*abs(ptho(i,j,k)+3.)*ocetra(i,j,k,idms)       &
     &             *(ocetra(i,j,k,idms)/(dmsp6+ocetra(i,j,k,idms)))	 
	 dms_uv  = dmsp2*4.*pho*ocetra(i,j,k,idms)

         dtr=bacfra-phosy+graton+ecan*zoomor 

         ocetra(i,j,k,iphosph)=                                        &
     &   ((ocetra(i,j,k,iphosph)+dtr)*bdp                              &
     &    +ocetra(i,j,k,iphosph)     *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)

         ocetra(i,j,k,iano3)=                                          &
     &   ((ocetra(i,j,k,iano3)+dtr*rnit)*bdp                           &
     &    +ocetra(i,j,k,iano3)          *(pddpo(i,j,k)-bdp))/pddpo(i,j,k) 

         ocetra(i,j,k,idet)=                                           &
     &   ((ocetra(i,j,k,idet)+export)*bdp                              &
     &    +ocetra(i,j,k,idet)        *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)

! new from emr version 4.5.06
#ifdef __c_isotopes
         ocetra(i,j,k,idet14)=                                        &
     &   ((ocetra(i,j,k,idet14)+rcar*export*bifr14)*bdp               &
     &    +ocetra(i,j,k,idet14)*(pddpo(i,j,k)-bdp))/pddpo(i,j,k) 

         ocetra(i,j,k,idet13)=                                        &
     &   ((ocetra(i,j,k,idet13)+rcar*export*bifr13)*bdp              &
     &    +ocetra(i,j,k,idet13)*(pddpo(i,j,k)-bdp))/pddpo(i,j,k)   

         ocetra(i,j,k,isco214)=                                        &
     &   ((ocetra(i,j,k,isco214)-rcar*export*bifr14)*bdp    &
     &    +ocetra(i,j,k,isco214)*(pddpo(i,j,k)-bdp))/pddpo(i,j,k)    

! 13C, 14C 'removal' by exportproduction
! are these two needed? ('delcar' is removed later) (also there is no equivalent for isco212)
! short term removal from dissolved pool "for gas exchange" (in carchm), effect is fairly small
! 12C in P-units, 13C/14C in C-units (*122)  

         ocetra(i,j,k,isco213)=                                        &
     &   ((ocetra(i,j,k,isco213) -rcar*export* bifr13 * flui13)*bdp    &
     &    +ocetra(i,j,k,isco213)*(pddpo(i,j,k)-bdp))/pddpo(i,j,k)           !bifr = biogenic fractionation (0.98)  ! *flui13
         ocetra(i,j,k,isco214)=                                        &
     &   ((ocetra(i,j,k,isco214) -rcar*export* bifr14 * flui14)*bdp    &
     &    +ocetra(i,j,k,isco214)*(pddpo(i,j,k)-bdp))/pddpo(i,j,k)        ! *rocean14 (--> flui14?) 4.5.06

         ocetra(i,j,k,idet13)=                                        &
     &   ((ocetra(i,j,k,idet13) +rcar*export* bifr13 * flui13)*bdp    &
     &    +ocetra(i,j,k,idet13)*(pddpo(i,j,k)-bdp))/pddpo(i,j,k)         ! det in P-units, det13 in C-units (*122)
         ocetra(i,j,k,idet14)=                                        &
     &   ((ocetra(i,j,k,idet14) +rcar*export* bifr14 * flui14)*bdp    &
     &    +ocetra(i,j,k,idet14)*(pddpo(i,j,k)-bdp))/pddpo(i,j,k)         ! det in P-units, det13 in C-units (*122)
#endif /*__c_isotopes*/

         ocetra(i,j,k,idms) =                                          &
     &   ((ocetra(i,j,k,idms)+dmsprod-dms_bac-dms_uv)*bdp              &
     &	  +ocetra(i,j,k,idms)        *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)  

         ocetra(i,j,k,isco212)=                                        &
     &   ((ocetra(i,j,k,isco212)-delcar+rcar*dtr)*bdp                  &
     &    +ocetra(i,j,k,isco212)     *(pddpo(i,j,k)-bdp))/pddpo(i,j,k) 

#ifdef __c_isotopes
         ocetra(i,j,k,isco213)=                                        &
     &   ((ocetra(i,j,k,isco213)-delcar*rocean13                       &
     &   +bifr13*rcar*( bacfra - phosy + graton + ecan*zoomor))*bdp     & ! bifr=0.98
     &   +ocetra(i,j,k,isco213)     *(pddpo(i,j,k)-bdp))/pddpo(i,j,k) 

          ocetra(i,j,k,isco214)=                                       &
     &    ((ocetra(i,j,k,isco214)-delcar*rocean14)*bdp                 &
     &    +ocetra(i,j,k,isco214)     *(pddpo(i,j,k)-bdp))/pddpo(i,j,k) 
! js: for efficency below line (which should in principle be there) is neglected (additional tracer field would be needed 
!     to account for radioactive decay of 14C in particles)
!     &      + bifr14*rcar*( bacfra - phosy + graton + ecan*zoomor)
#endif /*__c_isotopes*/

         ocetra(i,j,k,ialkali)=                                        &
     &   ((ocetra(i,j,k,ialkali)-2.*delcar-rnit*dtr)*bdp               &
     &    +ocetra(i,j,k,ialkali)     *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)

         ocetra(i,j,k,ioxygen)=                                        &
     &   ((ocetra(i,j,k,ioxygen)+ro2ut*(phosy-bacfra)                  &
     &                          -(graton+ecan*zoomor)*ro2ut)*bdp       &
     &    +ocetra(i,j,k,ioxygen)     *(pddpo(i,j,k)-bdp))/pddpo(i,j,k) 

        ocetra(i,j,k,iphy)=                                           &
     &   ((ocetra(i,j,k,iphy)+phosy-grazing-phymor-exud)*bdp           &
     &    +ocetra(i,j,k,iphy)     *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)

         ocetra(i,j,k,izoo)=                                           &
     &   ((ocetra(i,j,k,izoo)+grawa-excdoc-zoomor)*bdp                 &
     &    +ocetra(i,j,k,izoo)     *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)

!         ocetra(i,j,k,iphy)=                                           &
!     &     (ocetra(i,j,k,iphy)*(pddpo(i,j,k)/bdp)+phosy-grazing-phymor-exud)&
!     &    *(bdp/pddpo(i,j,k))       

!         ocetra(i,j,k,izoo)=                                           &
!     &     (ocetra(i,j,k,izoo)*(pddpo(i,j,k)/bdp)+grawa-excdoc-zoomor) &
!     &    *(bdp/pddpo(i,j,k))       

        ocetra(i,j,k,idoc)=                                           &
     &   ((ocetra(i,j,k,idoc)-bacfra+excdoc+exud)*bdp                  &
     &    +ocetra(i,j,k,idoc)       *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)  

         ocetra(i,j,k,icalc)=                                          &
     &   ((ocetra(i,j,k,icalc)+delcar)*bdp                             &
     &    +ocetra(i,j,k,icalc)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k) 

#ifdef __c_isotopes
         ocetra(i,j,k,icalc13)=                                        &
     &   ((ocetra(i,j,k,icalc13)+delcar*rocean13)*bdp                  &
     &    +ocetra(i,j,k,icalc13)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k) 

         ocetra(i,j,k,icalc14)=                                        &
     &   ((ocetra(i,j,k,icalc14)+delcar*rocean14)*bdp                  &
     &    +ocetra(i,j,k,icalc14)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k) 
#endif

         ocetra(i,j,k,isilica)=                                        &
     &  ((ocetra(i,j,k,isilica)-delsil+dremopal*ocetra(i,j,k,iopal))*bdp&
     &    +ocetra(i,j,k,isilica)    *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)  

         ocetra(i,j,k,iopal)=                                          &
     &   ((ocetra(i,j,k,iopal)+delsil-dremopal*ocetra(i,j,k,iopal))*bdp &
     &    +ocetra(i,j,k,iopal)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)  

         ocetra(i,j,k,iiron)=                                          &
     &   ((ocetra(i,j,k,iiron)+dtr*riron                               &
     &           -relaxfe*MAX(ocetra(i,j,k,iiron)-fesoly,0.))*bdp      &
     &    +ocetra(i,j,k,iiron)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)  


#ifdef AGG
!***********************************************************************
! effects of biological processes on number of particles:
! photosynthesis creates POM
! exudation deletes POM
! grazing deletes POM; but only the fraction that is not egested as 
! fecal pellets again (grawa remains in zoo, graton goes to po4)     
! none of the processes at the current time is assumed to change
! the size distribution (subject to change)
! NOTE that phosy, exud etc. are in kmol/m3! 
! Thus divide by avmass (kmol/m3)
!**********************************************************************

        if(avmass.gt.0.) then
           avnos = ocetra(i,j,k,inos) 
           anosloss = (phosy-exud-graton-grawa)*avnos/avmass
           ocetra(i,j,k,inos) = ocetra(i,j,k,inos)+anosloss
        endif  

!***********************************************************************
! dead zooplankton corpses come with their own, flat distribution
! this flow even takes place if there is neither nos nor mass
! NOTE: zoomor is in kmol/m3!! Thus multiply flow by 1.e+6
!***********************************************************************

      zmornos = zoomor * (1.-ecan) * zdis * 1.e+6
      ocetra(i,j,k,inos) = ocetra(i,j,k,inos)+zmornos

#endif /*AGG*/

!
! add up for total inventory
!
         expoor(i,j)=expoor(i,j)+bdp*export*rcar
         expoca(i,j)=expoca(i,j)+bdp*delcar
         exposi(i,j)=exposi(i,j)+bdp*delsil
!
! write output for bgcmean
!

!kma
         aux2d_dmsprod(i,j)   = aux2d_dmsprod(i,j)+dmsprod*bdp 
         aux2d_dms_bac(i,j)   = aux2d_dms_bac(i,j)+dms_bac*bdp 
         aux2d_dms_uv(i,j)    = aux2d_dms_uv (i,j)+dms_uv*bdp 
         aux2d_export(i,j)    = aux2d_export(i,j) +export*rcar*bdp 
         aux2d_expoca(i,j)    = aux2d_expoca(i,j) +delcar*bdp
         aux2d_exposi(i,j)    = aux2d_exposi(i,j) +delsil*bdp 
!kma primary production in g C m-2
         aux3d_phosy(i,j,k)   = phosy*rcar*(bdp/pddpo(i,j,k))  

!      psum1=                                                          &
!     &   ocetra(i,j,k,idet)+ocetra(i,j,k,idoc)+ocetra(i,j,k,iphy)     &
!     &  +ocetra(i,j,k,izoo)+ocetra(i,j,k,iphosph)  

!      IF(abs(psum1-psum0).gt.1.e-12) THEN
!        write(io_stdo_bgc,*) 'OCPROD EZ P',psum1,psum0
!      ENDIF
     
      ENDIF      ! pddpo(i,j,k).GT.0.5
100   CONTINUE   ! kwrbioz
1     CONTINUE   ! kpie, kpje
!$OMP END PARALLEL DO

#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'in OCRPOD after 1st bio prod'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask &
     &                  ,0)
!      call INVENTORY_BGC(kpie,kpje,kpke)
#endif 
!ib
! Accumulate 2d diagnostics 
       call accsrf(jdmsprod,aux2d_dmsprod,omask,0)    
       call accsrf(jdms_bac,aux2d_dms_bac,omask,0)  
       call accsrf(jdms_uv,aux2d_dms_uv,omask,0)     
       call accsrf(jexport,aux2d_export,omask,0)      
       call accsrf(jexpoca,aux2d_expoca,omask,0)     
       call accsrf(jexposi,aux2d_exposi,omask,0)     

! Accumulate primary production 
      call acclyr(jphosy,aux3d_phosy,pddpo,1)
      IF (SUM(jlvlphosy).NE.0) THEN
        DO k=1,kpke
          call bgczlv(pddpo,k,ind1,ind2,wghts)
          call acclvl(jlvlphosy,aux3d_phosy,k,ind1,ind2,wghts)
        ENDDO
      ENDIF
!ib


!      pinv1=0.0
!      DO j=1,kpje
!      DO i=1,kpie
!        IF(omask(i,j).gt.0.0) THEN
!      DO k=1,kpke
!      dpinv=                                                      &
!     &   ocetra(i,j,k,idet)+ocetra(i,j,k,idoc)+ocetra(i,j,k,iphy) &
!     &  +ocetra(i,j,k,izoo)+ocetra(i,j,k,iphosph)                  
!      pinv1=pinv1+dpinv*pddpo(i,j,k)*pdlxp(i,j)*pdlyp(i,j)      
!      ENDDO
!        ENDIF
!      ENDDO
!      ENDDO

!      WRITE(io_stdo_bgc,*) 'PINV BIOPROD',pinv0,pinv1
!      pinv0=pinv1

!bioprod
#endif       
!bioprod

#ifdef AGG
       DO  k=1,kpke
         DO j=1,kpje
           DO i=1,kpie
            IF(pddpo(i,j,k).gt.1.e-6.and.omask(i,j).GT.0.5) THEN
             avmass = ocetra(i,j,k,iphy) + ocetra(i,j,k,idet)
             snow = avmass*1.e+6
!  check whether the numbers had to be decreased or increased
             if (snow*pupper.gt.ocetra(i,j,k,inos)) then
               nosin2 = nosin2 + 1
             endif
             if (snow/cellmass.lt.ocetra(i,j,k,inos)) then
               nosde2 = nosde2 + 1
             endif  
            ENDIF
           ENDDO
       ENDDO
      ENDDO
#endif /*AGG*/

!      CALL contro(148)

!$OMP PARALLEL DO                                                & 
!$OMP&PRIVATE(sterph,sterzo,remin,docrem,opalrem,aou, &
!$OMP&        refra,dms_bac,detref,rem13,rem14,phythresh,  &
!$OMP&        pocrem,phyrem,bdp,dp_ez)  
      DO 201 j=1,kpje
      DO 201 i=1,kpie
         DO 20 k=kwrbioz(i,j),kpke
            IF(pddpo(i,j,k).gt.1.e-12.and.omask(i,j).GT.0.5) THEN

! depth of euphotic zone
      dp_ez=100.0
      if(ptiestw(i,j,k).lt.dp_ez.and.ptiestw(i,j,k+1).gt.dp_ez) then
        bdp=abs(ptiestw(i,j,k+1)-dp_ez)
        if(bdp.gt.pddpo(i,j,k)) then
          write(io_stdo_bgc,*) 'bdp gt pddpo(i,j,k) 2',i,j,k
          write(io_stdo_bgc,*) 'bdp gt pddpo 2',pddpo(i,j,k),bdp
          write(io_stdo_bgc,*) 'bdp gt pddpo 2',ptiestw(i,j,k),ptiestw(i,j,k+1)
        endif
      else
        bdp=pddpo(i,j,k)
      endif

#ifdef AGG
            avmass=ocetra(i,j,k,iphy)+ocetra(i,j,k,idet)
#endif /*AGG*/	    
            phythresh=MAX(0.,ocetra(i,j,k,iphy)-phytomi)
            sterph=0.5*dphymor* phythresh                               !phytoplankton to detritus
            sterzo=dzoomor*MAX(0.,ocetra(i,j,k,izoo)-grami)   

       	    ocetra(i,j,k,iphy)=((ocetra(i,j,k,iphy)-sterph)*bdp         & 
     &     +ocetra(i,j,k,iphy)         *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
       	    ocetra(i,j,k,izoo)=((ocetra(i,j,k,izoo)-sterzo)*bdp         & 
     &     +ocetra(i,j,k,izoo)         *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)

            IF(ocetra(i,j,k,ioxygen).gt.5.e-8) THEN
               pocrem=MIN(drempoc*ocetra(i,j,k,idet),                  &
     &                   0.33*ocetra(i,j,k,ioxygen)/ro2ut)
               docrem=MIN(dremdoc*ocetra(i,j,k,idoc),                 &
     &                   0.33*ocetra(i,j,k,ioxygen)/ro2ut)
               if (docrem.lt.0.) then
                write(*,*)'negative doc remineralization',docrem
                docrem=0.;
               endif
               phyrem=MIN(0.5*dphymor*phythresh,                       &
     &                   0.33*ocetra(i,j,k,ioxygen)/ro2ut)
               detref=pocrem/(ocetra(i,j,k,idet)+1.e-20)                      ! 'detritus remineralized fraction' (?)
#ifdef __c_isotopes
               rem13=detref*ocetra(i,j,k,idet13)                             ! remineralization of poc13
               rem14=detref*ocetra(i,j,k,idet14)                             !                     poc14
#endif
            else
               pocrem=0.
               docrem=0.
               phyrem=0.
#ifdef __c_isotopes
               rem13 =0.
               rem14 =0.
#endif
	       docrem=0.
            endif 
	    
            ocetra(i,j,k,idet)=                                        &
     &    ((ocetra(i,j,k,idet)-pocrem+sterph+sterzo)*bdp               &
     &     +ocetra(i,j,k,idet)        *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
#ifdef __c_isotopes
            ocetra(i,j,k,idet13)=                                      &
     &    ((ocetra(i,j,k,idet13)+rcar*bifr13*(sterph+sterzo)-rem13)*bdp&
     &     +ocetra(i,j,k,idet13)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)

            ocetra(i,j,k,idet14)=                                      &
     &    ((ocetra(i,j,k,idet14)-rem14)*bdp                            &
     &     +ocetra(i,j,k,idet14)       *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
#endif
            ocetra(i,j,k,idoc)=                                        &
     &    ((ocetra(i,j,k,idoc)-docrem)*bdp                             &
     &     +ocetra(i,j,k,idoc)        *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
            ocetra(i,j,k,iphy)=                                        &
     &    ((ocetra(i,j,k,iphy)-phyrem)*bdp                             &
     &     +ocetra(i,j,k,iphy)        *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)

            remin=pocrem+docrem+phyrem

            ocetra(i,j,k,iphosph)=                                     &
     &    ((ocetra(i,j,k,iphosph)+remin)*bdp                           &
     &     +ocetra(i,j,k,iphosph)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
            ocetra(i,j,k,iano3)=                                       &
     &    ((ocetra(i,j,k,iano3)+remin*rnit)*bdp                        &
     &     +ocetra(i,j,k,iano3)        *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
            ocetra(i,j,k,isco212)=                                     &
     &    ((ocetra(i,j,k,isco212)+rcar*remin)*bdp                      &
     &     +ocetra(i,j,k,isco212)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
            ocetra(i,j,k,ialkali)=                                     &
     &    ((ocetra(i,j,k,ialkali)-rnit*remin)*bdp                      &
     &     +ocetra(i,j,k,ialkali)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
            ocetra(i,j,k,ioxygen)=                                     &
     &    ((ocetra(i,j,k,ioxygen)-ro2ut*remin)*bdp                     &
     &     +ocetra(i,j,k,ioxygen)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
            ocetra(i,j,k,iiron)=                                       &
     &    ((ocetra(i,j,k,iiron)+remin*riron                            &
     &        -relaxfe*MAX(ocetra(i,j,k,iiron)-fesoly,0.))*bdp         &
     &     +ocetra(i,j,k,iiron)        *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)

#ifdef __c_isotopes
            ocetra(i,j,k,isco213)=                                     &
     &    ((ocetra(i,j,k,isco213)+rcar*docrem*bifr13 + rem13)*bdp      &
     &     +ocetra(i,j,k,isco213)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
          ! rem13 not *rcar as idet13 is in C-units

!           ocetra(i,j,k,isco214)= ocetra(i,j,k,isco214)                          ! no biogenic component for 14C
#endif
!***********************************************************************
! as ragueneau (2000) notes, Si(OH)4sat is about 1000 umol, but
! Si(OH)4 varies only between 0-100 umol
! so the expression dremopal*(Si(OH)4sat-Si(OH)4) would change the 
! rate only from 0 to 100%     
!***********************************************************************
            opalrem=dremopal*0.1*(ptho(i,j,k)+3.)*ocetra(i,j,k,iopal)
            ocetra(i,j,k,iopal)=                                       &
     &    ((ocetra(i,j,k,iopal)-opalrem)*bdp                           &
     &     +ocetra(i,j,k,iopal)        *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
            ocetra(i,j,k,isilica)=                                     &
     &    ((ocetra(i,j,k,isilica)+opalrem)*bdp                         &
     &     +ocetra(i,j,k,isilica)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)

!***********************************************************************
!           There is about 1.e4 O2 on 1 N2O molecule (Broeker&Peng)
!           refra : Tim Rixton, private communication
!***********************************************************************
            aou=satoxy(i,j,k)-ocetra(i,j,k,ioxygen)
            refra=1.+3.*(0.5+sign(0.5,aou-1.97e-4))
            ocetra(i,j,k,ian2o)=                                       &
     &    ((ocetra(i,j,k,ian2o)+remin*1.e-4*ro2ut*refra)*bdp           &
     &     +ocetra(i,j,k,ian2o)        *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
            ocetra(i,j,k,igasnit)=                                     &
     &    ((ocetra(i,j,k,igasnit)-remin*1.e-4*ro2ut*refra)*bdp         &
     &     +ocetra(i,j,k,igasnit)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
            ocetra(i,j,k,ioxygen)=                                     &
     &    ((ocetra(i,j,k,ioxygen)-remin*1.e-4*ro2ut*refra*0.5)*bdp     &
     &     +ocetra(i,j,k,ioxygen)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)

!careful, pho not defined this far down
!            ocetra(i,j,k,idms)=ocetra(i,j,k,idms)                      &
!     &                        -dmsp2*8.*pho*ocetra(i,j,k,idms)         &
!     &           -dmsp3*abs(ptho(i,j,k)+3.)*ocetra(i,j,k,idms)

	 dms_bac = dmsp3*abs(ptho(i,j,k)+3.)*ocetra(i,j,k,idms)	 
	 
         ocetra(i,j,k,idms)=                                        &
     & ((ocetra(i,j,k,idms)-dms_bac)*bdp                            &
     &  +ocetra(i,j,k,idms)         *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
	 
#ifdef AGG
!***********************************************************************
! loss of snow numbers due to remineralization of poc
! gain of snow numbers due to zooplankton mortality
! NOTE that remin is in kmol/m3. Thus divide by avmass (kmol/m3)
!***********************************************************************
           if(avmass.gt.0.) then  
              avnos = ocetra(i,j,k,inos)
              ocetra(i,j,k,inos) = ocetra(i,j,k,inos)                  & 
     &                           - remin * avnos/avmass
           endif
!***********************************************************************
! dead zooplankton corpses come with their own, flat distribution
! this flow even takes place if there is neither nos nor mass
! NOTE: zoomor is in kmol/m3!! Thus multiply flow by 1.e+6
!***********************************************************************
           zmornos = sterzo * zdis * 1.e+6
           ocetra(i,j,k,inos) = ocetra(i,j,k,inos) + zmornos 
#endif /*AGG*/

            ENDIF
 20      CONTINUE
 201   CONTINUE     
!$OMP END PARALLEL DO

#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'in OCRPOD after poc remin'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask &
     &                  ,0)
!      call INVENTORY_BGC(kpie,kpje,kpke)
#endif 
!$OMP PARALLEL DO                                                   &
!$OMP&PRIVATE(remin,remin2o,rem13,rem14,detref,detrl,rl13,rl14,bdp,dp_ez) 
       DO 30 j=1,kpje
       DO 30 i=1,kpie
         DO 30 k=kwrbioz(i,j),kpke
         IF(omask(i,j).GT.0.5) THEN
         IF(ocetra(i,j,k,ioxygen).LT.5.e-7.and.pddpo(i,j,k).gt.1.e-12) THEN

! depth of euphotic zone
!JT
      dp_ez=100.0
!      dp_ez=90.0
      if(ptiestw(i,j,k).lt.dp_ez.and.ptiestw(i,j,k+1).gt.dp_ez) then
        bdp=abs(ptiestw(i,j,k+1)-dp_ez)
        if(bdp.gt.pddpo(i,j,k)) then
          write(io_stdo_bgc,*) 'bdp gt pddpo(i,j,k) 3',i,j,k
        endif
      else
        bdp=pddpo(i,j,k)
      endif

#ifdef AGG
               avmass = ocetra(i,j,k,iphy)+ocetra(i,j,k,idet)
#endif /*AGG*/	    
	    
               remin=0.05*drempoc*MIN(ocetra(i,j,k,idet),               &
     &                           0.5*ocetra(i,j,k,iano3)/rnit23)

               detref=remin/(ocetra(i,j,k,idet)+1.e-60)                        ! P-units
#ifdef __c_isotopes
               rem13=detref*ocetra(i,j,k,idet13)                               ! C-units
               rem14=detref*ocetra(i,j,k,idet14)                               ! C-units
#endif

               remin2o=dremn2o*MIN(ocetra(i,j,k,idet),                    &
     &	                        0.003*ocetra(i,j,k,ian2o)/(2*ro2ut))
               detrl=remin2o/(ocetra(i,j,k,idet)+1.e-60)                       ! detrl?
#ifdef __c_isotopes
               rl13=detrl*ocetra(i,j,k,idet13)                                 ! C-units
               rl14=detrl*ocetra(i,j,k,idet14)                                 ! C-units
#endif

            ocetra(i,j,k,ialkali)=                                     &
     &    ((ocetra(i,j,k,ialkali)-rnit*(remin+remin2o))*bdp            &
     &     +ocetra(i,j,k,ialkali)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
            ocetra(i,j,k,isco212)=                                     &
     &    ((ocetra(i,j,k,isco212)+rcar*(remin+remin2o))*bdp            &
     &     +ocetra(i,j,k,isco212)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)

! proxies 13C, 14C 
#ifdef __c_isotopes
            ocetra(i,j,k,isco213)=                                     &
     &    ((ocetra(i,j,k,isco213)+(rem13+rl13))*bdp                    &
     &     +ocetra(i,j,k,isco213)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
!    &                               +rcar* (rem13+rl13  ) ! changed 3.5.2006
            ocetra(i,j,k,isco214)=                                     &
     &    ((ocetra(i,j,k,isco214)+(rem14+rl14))*bdp                    &
     &     +ocetra(i,j,k,isco214)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
#endif

            ocetra(i,j,k,idet)=                                        &
     &    ((ocetra(i,j,k,idet)-(remin+remin2o))*bdp                    &
     &     +ocetra(i,j,k,idet)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
! proxies
#ifdef __c_isotopes
            ocetra(i,j,k,idet13)=                                     &
     &    ((ocetra(i,j,k,idet13)-(rem13+rl13))*bdp                    &
     &     +ocetra(i,j,k,idet13)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
            ocetra(i,j,k,idet14)=                                     &
     &    ((ocetra(i,j,k,idet14)-(rem14+rl14))*bdp                    &
     &     +ocetra(i,j,k,idet14)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
#endif

            ocetra(i,j,k,iphosph)=                                     &
     &    ((ocetra(i,j,k,iphosph)+(remin + remin2o))*bdp               &
     &     +ocetra(i,j,k,iphosph)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
            ocetra(i,j,k,iano3)=                                       &
     &    ((ocetra(i,j,k,iano3)-rnit23*remin+rnit*(remin + remin2o))*bdp&
     &     +ocetra(i,j,k,iano3)        *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
            ocetra(i,j,k,igasnit)=                                     &
     &    ((ocetra(i,j,k,igasnit)+rnit13*remin + 2*ro2ut*remin2o)*bdp  &
     &     +ocetra(i,j,k,igasnit)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
            ocetra(i,j,k,iiron)=                                       &
     &    ((ocetra(i,j,k,iiron)+riron*(remin + remin2o))*bdp           &
     &     +ocetra(i,j,k,iiron)        *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
            ocetra(i,j,k,ian2o)=                                       &
     &    ((ocetra(i,j,k,ian2o)-2*ro2ut*remin2o)*bdp                   &
     &     +ocetra(i,j,k,ian2o)        *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)

#ifdef AGG
!***********************************************************************
! loss of snow numbers due to remineralization of poc
! NOTE that remin is in kmol/m3. Thus divide by avmass (kmol/m3)
!***********************************************************************
           if(avmass.gt.0.) then  
              avnos = ocetra(i,j,k,inos)
              ocetra(i,j,k,inos) = ocetra(i,j,k,inos)         &
     &                - (remin+remin2o)*avnos/avmass
           endif
#endif /*AGG*/

         ENDIF
         ENDIF
30    CONTINUE
!$OMP END PARALLEL DO
#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'in OCRPOD after remin n2o'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask &
     &                  ,0)
!      call INVENTORY_BGC(kpie,kpje,kpke)
#endif 
!sulphate reduction   ! introduced 11.5.2007 to improve poc-remineralisation in the 
!                       oxygen minimum zone in the subsurface equatorial Pacific
!                       assumption of endless pool of SO4 (typical concentration are on the order of mmol/l)
!      js 02072007:    for other runs than current millenium (cosmos-setup) experiments this seems
!                      to cause trouble as phosphate concentrations are too high at the depth of the oxygen
!                      minimum in the equatorial pacific/atlantic
!                      does it make sense to check for oxygen and nitrate deficit?

!$OMP PARALLEL DO                                                   &
!$OMP&PRIVATE(remin,rem13,rem14,detref,bdp,dp_ez) 
      DO 301 j=1,kpje
      DO 301 i=1,kpie
        DO 301 k=kwrbioz(i,j),kpke
            IF(omask(i,j).gt.0.5.and.pddpo(i,j,k).gt.1.e-12) then  
            IF(ocetra(i,j,k,ioxygen).lt.3.e-6.and.ocetra(i,j,k,iano3).lt.3.e-6) THEN

! depth of euphotic zone
!JT
      dp_ez=100.0
!JT      dp_ez=90.0
      if(ptiestw(i,j,k).lt.dp_ez.and.ptiestw(i,j,k+1).gt.dp_ez) then
        bdp=abs(ptiestw(i,j,k+1)-dp_ez)
        if(bdp.gt.pddpo(i,j,k)) then
          write(io_stdo_bgc,*) 'bdp gt pddpo(i,j,k) 4',i,j,k
        endif
      else
        bdp=pddpo(i,j,k)
      endif

#ifdef AGG
               avmass = ocetra(i,j,k,iphy)+ocetra(i,j,k,idet)
#endif /*AGG*/

               remin=dremsul*ocetra(i,j,k,idet)
               detref=dremsul
#ifdef __c_isotopes
               rem13=detref*ocetra(i,j,k,idet13)
               rem14=detref*ocetra(i,j,k,idet14)
#endif /*__c_isotopes*/

               ocetra(i,j,k,idet)=                                     &
     &       ((ocetra(i,j,k,idet)-remin)*bdp                           &
     &        +ocetra(i,j,k,idet)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
               ocetra(i,j,k,ialkali)=                                  &
     &       ((ocetra(i,j,k,ialkali)-rnit*remin)*bdp                   &
     &        +ocetra(i,j,k,ialkali)   *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
               ocetra(i,j,k,isco212)=                                  &
     &       ((ocetra(i,j,k,isco212)+rcar*remin)*bdp                   &
     &        +ocetra(i,j,k,isco212)   *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
#ifdef __c_isotopes
               ocetra(i,j,k,isco213)=                                  &
     &       ((ocetra(i,j,k,isco213)+rem13)*bdp                        &
     &        +ocetra(i,j,k,isco213)   *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
               ocetra(i,j,k,isco214)=                                  &
     &       ((ocetra(i,j,k,isco214)+rem14)*bdp                        &
     &        +ocetra(i,j,k,isco214)   *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
               ocetra(i,j,k,idet13)=                                   &
     &       ((ocetra(i,j,k,idet13)-rem13)*bdp                         &
     &        +ocetra(i,j,k,idet13)    *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
               ocetra(i,j,k,idet14)=                                   &
     &       ((ocetra(i,j,k,idet14)-rem14)*bdp                         &
     &        +ocetra(i,j,k,idet14)    *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
#endif /*__c_isotopes*/
               ocetra(i,j,k,iphosph)=                                  &
     &       ((ocetra(i,j,k,iphosph)+remin)*bdp                        &
     &        +ocetra(i,j,k,iphosph)   *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
       ! remin from sulphate reduction k=kwrbioz,ke
               ocetra(i,j,k,iano3)=                                    &
     &       ((ocetra(i,j,k,iano3)+rnit*remin)*bdp                   &
     &        +ocetra(i,j,k,iano3)     *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
               ocetra(i,j,k,iiron)=                                    &
     &       ((ocetra(i,j,k,iiron)+riron*remin)*bdp                    &
     &        +ocetra(i,j,k,iiron)    *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
   
#ifdef AGG
!***********************************************************************
! loss of snow numbers due to remineralization of poc
! NOTE that remin is in kmol/m3. Thus divide by avmass (kmol/m3)
!***********************************************************************
            if(avmass.gt.0.) then
               avnos = ocetra(i,j,k,inos)
               ocetra(i,j,k,inos) = ocetra(i,j,k,inos)         &
      &                - (remin)*avnos/avmass
            endif
#endif /*AGG*/

            ENDIF
            ENDIF
301    CONTINUE
!$OMP END PARALLEL DO
!    end sulphate reduction
#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'in OCRPOD after sulphate reduction '
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask &
     &                  ,0)
!      call INVENTORY_BGC(kpie,kpje,kpke)
#endif 
#ifdef AGG
       DO  k=1,kpke
         DO j=1,kpje
           DO i=1,kpie
            IF(pddpo(i,j,k).gt.1.e-6.and.omask(i,j).GT.0.5) THEN
             avmass = ocetra(i,j,k,iphy) + ocetra(i,j,k,idet)
             snow = avmass*1.e+6
! check whether the numbers had to be decreased or increased
             if (snow*pupper.gt.ocetra(i,j,k,inos)) then
               nosin3 = nosin3 + 1
             endif
             if (snow/cellmass.lt.ocetra(i,j,k,inos)) then
               nosde3 = nosde3 + 1
             endif
            ENDIF
           ENDDO
       ENDDO
      ENDDO
#endif /*AGG*/

#ifdef AGG

! **********************AGGREGATION************************************
! General:
! Sinking speed, size distribution and aggregation are calculated 
! as in Kriest and Evans, 2000.
! I assume that opal and calcium carbonate sink at the same speed as P (mass).
!
! Sinking speed and aggregation: I assume that if there is no phosphorous mass,
! the sinking speed is the maximal sinking speed of aggregates. I further
! assume that then there are no particles, and that the rate of aggregation
! is 0. This scheme removes no P in the absence of P, but still opal and/or
! calcium carbonate.
! This could or should be changed, because silica as well as carbonate
! shell will add to the aggregate mass, and should be considered.
! Puh. Does anyone know functional relationships between
! size and Si or CaCO3? Perhaps on a later version, I have to
! take the relationship bewteen weight and size?
!
! 1. Size distribution and resulting loss of marine snow aggregates due to aggregation 
! (aggregate(i,j,k)) and sinking speed of mass and numbers (wmass(i,j,k)
! and wnumb(i,j,k) are calculated in a loop over 2-kpke. 
!
! 2. The depth of the first layer may change due to ice drift, etc.
! This puts a restriction onto the maximum sinking speed.
! I currently set the max. size for size dependent sinking onto
! one appropriate for this depth.
!
! 3. The fluxes out of the last layer are calculated from sinking speed
! and mass concentration, and form the boundary condition for the sediment.
!
! 4. The fluxes in layer kpke->2 are calculated from sinking speed and mass
! concentration, sinking speed, aggregation and number concentration. 

! 5. The fluxes in layer 1 are calculated from sinking speed and mass
! concentration, sinking speed, aggregation and number concentration.  (??)
!************************************************************************

      do k=2,kpke
      do i=1,kpie
      do j=1,kpje
         IF(pddpo(i,j,k).gt.1.e-6.and.omask(i,j).GT.0.5) THEN
          avm = ocetra(i,j,k,iphy)+ocetra(i,j,k,idet)
          if(avm.gt.0.) then
          snow = avm*1.e+6
          avn = ocetra(i,j,k,inos)
          eps = ((1.+ FractDim)*snow-avn*cellmass) /                   &
     &           (snow-avn*cellmass)

! prevent epsilon from becoming exactly one of the values which are 
! needed for the division (guide from??js)
          if (abs(eps-3.).lt.1.e-15) eps=3.+ vsmall
          if (abs(eps-4.).lt.1.e-15) eps=4.+ vsmall
          if (abs(eps-3.-SinkExp).lt.1.e-15) eps=3.+SinkExp+vsmall
          if (abs(eps-1.-SinkExp-FractDim).lt.1.e-15)                  &
     &        eps=1.+SinkExp+FractDim+vsmall

          e1 = 1. - eps
          e2 = 2. - eps
          e3 = 3. - eps
          e4 = 4. - eps
          es1 = e1 + SinkExp
          es3 = e3 + SinkExp
          TopF = (alar1/alow1)**e1
          TopM = TopF*TMFac

! SINKING SPEED FOR THIS LAYER
          wmass(i,j,k) = cellsink * ( (FractDim+e1)/ (FractDim+es1)    &
     &         +TopM*TSFac*SinkExp/ (FractDim+es1))
          wnumb(i,j,k) = cellsink * (e1/es1+TopF*TSFac*SinkExp/es1)

! AGGREGATION

! As a first step, assume that shear in the upper 4 layers is high and 
! zero below. Subject to change. js: this is for 20 layer version.
!                                    include 40 layer version
!                                  should be replaced by check for detph. done 29072005js
!         if (k.lt.5.and.ke.eq.20.or.k.lt.10.and.ke.eq.40) then
          if (ptiestu(i,j,k).le.100.) then
            fshear = fsh
          else
            fshear = 0.
          endif     


! shear kernel:
      sagg1 = (TopF-1.)*(TopF*alar3-alow3)*e1/e4                       &
     &   + 3.*(TopF*alar1-alow1)*(TopF*alar2-alow2)*e1*e1/(e2*e3)
      sagg2 = TopF*(                                                   &
     &    (alar3+3.*(alar2*alow1*e1/e2+alar1*alow2*e1/e3)+alow3*e1/e4) &
     &   - TopF*alar3*(1.+3*(       e1/e2+       e1/e3)+     e1/e4))
      sagg4 = TopF*TopF*4.*alar3
      shear_agg = (sagg1+sagg2+sagg4)*fshear

! settlement kernel:
      sagg1 = (TopF * TopF * alar2 * TSFac - alow2)                    &
     &   * SinkExp / (es3 * e3 * (es3 + e1))                           &
     &   + alow2 * ((1. - TopF * TSFac) / (e3 * es1)                   &
     &   - (1. - TopF) / (es3*e1))
      sagg2 = TopF * e1 * (TSFac * ( alow2 - TopF * alar2) / e3        &
     &   - (alow2 - TopF * alar2 * TSFac) / es3)
      sett_agg =  (e1*e1*sagg1+sagg2)*fse
      
      effsti=Stick*(ocetra(i,j,k,iopal)*1.e+6/ropal)/                  &
     &  ((ocetra(i,j,k,iopal)*1.e+6/ropal)+snow)

      aggregate(i,j,k) = (shear_agg+sett_agg)*effsti*avn*avn

! dust aggregation:
! shear kernel:
      dfirst=dustd3+3.*dustd2*alar1+3.*dustd1*alar2+alar3
      dshagg=e1*fsh*(dfirst*TopF/e1-(                                  &
     &  (TopF-1.)/e1*dustd3+3.*(TopF*alar1-alow1)/e2*dustd2            &
     &   +3.*(TopF*alar2-alow2)/e3*dustd1+(TopF*alar3-alow3)/e4))

! settlement kernel:
      dsett=fse*dustd2*((e1+SinkExp*TopF*TSFac)/es1-dustsink/cellsink)
      
      dustagg(i,j,k) = effsti*avn*ocetra(i,j,k,ifdust)                 &
     &                *(dshagg+dsett)

      else
        wmass(i,j,k)=TSFac*cellsink
        wnumb(i,j,k)=0.
        aggregate(i,j,k)=0.
        dustagg(i,j,k)=0.
        ocetra(i,j,k,inos)=0.
      endif
        endif
      enddo
      enddo
      enddo

! EVALUATE SINKING RATE AND AGGREGATION FOR FIRST LAYER, WHICH MAY BE
! SMALLER THAN INITIALLY SET BECAUSE OF EVAPORATION, ICE ETC.

      DO j=1,kpje
      DO i=1,kpie
         if(omask(i,j).gt.0.5) then

!ik evaluate safe length scale for size dependent sinking and
!ik aggregation, and the resulting sinking rate and aggregation rate.
!ik zo may reduce the first layer depth to values that are small and
!ik may cause the sinking length to exceed that depth.
!ik to be safe, for this upper layer set the upper size such that
!ik loss due to sinking is at max the whole inventory of this box.
!ik aggregation will be calculated accordingly.

         checksize = (pddpo(i,j,1)/cellsink)**(1./SinkExp)*alow1
         if(alar1.gt.checksize) then
           nacheck=nacheck+1
         endif
         flar1 = MIN(alar1,checksize)
         flar2 = flar1 * flar1
         flar3 = flar2 * flar1
         fTSFac = (flar1/alow1)**SinkExp
         fTMFac = (flar1/alow1)**FractDim

! SIZE DITRIBUTION
         avm = ocetra(i,j,1,iphy)+ocetra(i,j,1,idet)
         if(avm.gt.0.) then
         snow = avm*1.e+6
         avn = ocetra(i,j,1,inos)
         eps = ((1.+ FractDim)*snow-avn*cellmass) /                    &
     &           (snow-avn*cellmass)
     
         if (abs(eps-3.).lt.1.e-15) eps=3.+ vsmall
         if (abs(eps-4.).lt.1.e-15) eps=4.+ vsmall
         if (abs(eps-3.-SinkExp).lt.1.e-15) eps=3.+SinkExp+vsmall
         if (abs(eps-1.-SinkExp-FractDim).lt.1.e-15)                   &
     &        eps=1.+SinkExp+FractDim+vsmall

         e1 = 1. - eps
         e2 = 2. - eps
         e3 = 3. - eps
         e4 = 4. - eps
         es1 = e1 + SinkExp
         es3 = e3 + SinkExp

         fTopF = (flar1/alow1)**e1
         fTopM = fTopF*fTMFac

! SINKING SPEEDS
         wmass(i,j,1) = cellsink * ( (FractDim+e1)/ (FractDim+es1)     &
     &          +fTopM*fTSFac*SinkExp/ (FractDim+es1))
         wnumb(i,j,1) = cellsink * (e1/es1+fTopF*fTSFac*SinkExp/es1)

! AGGREGATION
      sagg1 = (fTopF-1.)*(fTopF*flar3-alow3)*e1/e4                     &
     &   + 3.*(fTopF*flar1-alow1)*(fTopF*flar2-alow2)*e1*e1/(e2*e3)
      sagg2 = fTopF*(                                                  &
     &    (flar3+3.*(flar2*alow1*e1/e2+flar1*alow2*e1/e3)+alow3*e1/e4) &
     &   - fTopF*flar3*(1.+3*(       e1/e2+       e1/e3)+     e1/e4))
      sagg4 = fTopF*fTopF*4.*flar3
      shear_agg = (sagg1+sagg2+sagg4)*fsh

      sagg1 = (fTopF * fTopF * flar2 * fTSFac - alow2)                 &
     &   * SinkExp / (es3 * e3 * (es3 + e1))                           &
     &   + alow2 * ((1. - fTopF * fTSFac) / (e3 * es1)                 &
     &   - (1. - fTopF) / (es3*e1))
      sagg2 = fTopF * e1 * (fTSFac * ( alow2 - fTopF * flar2) / e3     &
     &   - (alow2 - fTopF * flar2 * fTSFac) / es3)
      sett_agg =  (e1*e1*sagg1+sagg2)*fse


      effsti=Stick*(ocetra(i,j,1,iopal)*1.e+6/ropal)/                  &
     &  ((ocetra(i,j,1,iopal)*1.e+6/ropal)+snow)

      aggregate(i,j,1) = (shear_agg+sett_agg)*effsti*avn*avn

! dust aggregation:
! shear kernel:
      dfirst=dustd3+3.*dustd2*flar1+3.*dustd1*flar2+flar3
      dshagg=e1*fsh*(dfirst*fTopF/e1-(                                 &
     &  (fTopF-1.)/e1*dustd3+3.*(fTopF*flar1-alow1)/e2*dustd2          &
     &   +3.*(fTopF*flar2-alow2)/e3*dustd1+(fTopF*flar3-alow3)/e4))

! settlement kernel:
      dsett=fse*dustd2*((e1+SinkExp*fTopF*fTSFac)/es1-dustsink/cellsink)
      
      dustagg(i,j,1) = effsti*avn*ocetra(i,j,1,ifdust)                 &
     &                *(dshagg+dsett)

      else
        wmass(i,j,1)=fTSFac*cellsink
        wnumb(i,j,1)=0.
        aggregate(i,j,1)=0.
        dustagg(i,j,1)=0.
        ocetra(i,j,1,inos)=0.
      endif
      endif
      enddo
      enddo

! EVALUATE SINKING RATE AND AGGREGATION FOR LAST LAYER, WHICH MAY BE
! SMALLER THAN THE MINIMUM LAYER DEPTH

      DO j=1,kpje
      DO i=1,kpie
         if(omask(i,j).gt.0.5) then
         if(alar1max(i,j).lt.alar1) then

!ik take safe length scale for size dependent sinking and
!ik aggregation, and the resulting sinking rate and aggregation rate.

         flar1 = alar1max(i,j)
         flar2 = flar1 * flar1
         flar3 = flar2 * flar1
         fTSFac = TSFmax(i,j)
         fTMFac = TMFmax(i,j)

! SIZE DITRIBUTION
         avm = ocetra(i,j,kbo(i,j),iphy)+ocetra(i,j,kbo(i,j),idet)
         if(avm.gt.0.) then
         snow = avm*1.e+6
         avn = ocetra(i,j,kbo(i,j),inos)
         eps = ((1.+ FractDim)*snow-avn*cellmass) /                    &
     &           (snow-avn*cellmass)
     
         if (abs(eps-3.).lt.1.e-15) eps=3.+ vsmall
         if (abs(eps-4.).lt.1.e-15) eps=4.+ vsmall
         if (abs(eps-3.-SinkExp).lt.1.e-15) eps=3.+SinkExp+vsmall
         if (abs(eps-1.-SinkExp-FractDim).lt.1.e-15)                   &
     &        eps=1.+SinkExp+FractDim+vsmall

         e1 = 1. - eps
         e2 = 2. - eps
         e3 = 3. - eps
         e4 = 4. - eps
         es1 = e1 + SinkExp
         es3 = e3 + SinkExp

         fTopF = (flar1/alow1)**e1
         fTopM = fTopF*fTMFac

! SINKING SPEEDS
         wmass(i,j,kbo(i,j)) = cellsink *                              &
     &        ( (FractDim+e1)/ (FractDim+es1)                          &
     &          +fTopM*fTSFac*SinkExp/ (FractDim+es1))
         wnumb(i,j,kbo(i,j)) = cellsink *                              &
     &          (e1/es1+fTopF*fTSFac*SinkExp/es1)

! AGGREGATION
      sagg1 = (fTopF-1.)*(fTopF*flar3-alow3)*e1/e4                     &
     &   + 3.*(fTopF*flar1-alow1)*(fTopF*flar2-alow2)*e1*e1/(e2*e3)
      sagg2 = fTopF*(                                                  &
     &    (flar3+3.*(flar2*alow1*e1/e2+flar1*alow2*e1/e3)+alow3*e1/e4) &
     &   - fTopF*flar3*(1.+3*(       e1/e2+       e1/e3)+     e1/e4))
      sagg4 = fTopF*fTopF*4.*flar3
      shear_agg = (sagg1+sagg2+sagg4)*fsh

      sagg1 = (fTopF * fTopF * flar2 * fTSFac - alow2)                 &
     &   * SinkExp / (es3 * e3 * (es3 + e1))                           &
     &   + alow2 * ((1. - fTopF * fTSFac) / (e3 * es1)                 &
     &   - (1. - fTopF) / (es3*e1))
      sagg2 = fTopF * e1 * (fTSFac * ( alow2 - fTopF * flar2) / e3     &
     &   - (alow2 - fTopF * flar2 * fTSFac) / es3)
      sett_agg =  (e1*e1*sagg1+sagg2)*fse


      effsti=Stick*(ocetra(i,j,kbo(i,j),iopal)*1.e+6/ropal)/           &
     &  ((ocetra(i,j,kbo(i,j),iopal)*1.e+6/ropal)+snow)

      aggregate(i,j,kbo(i,j)) = (shear_agg+sett_agg)*effsti*avn*avn

! dust aggregation:
! shear kernel:
      dfirst=dustd3+3.*dustd2*flar1+3.*dustd1*flar2+flar3
      dshagg=e1*fsh*(dfirst*fTopF/e1-(                                 &
     &  (fTopF-1.)/e1*dustd3+3.*(fTopF*flar1-alow1)/e2*dustd2          &
     &   +3.*(fTopF*flar2-alow2)/e3*dustd1+(fTopF*flar3-alow3)/e4))

! settlement kernel:
      dsett=fse*dustd2*((e1+SinkExp*fTopF*fTSFac)/es1-dustsink/cellsink)
      
      dustagg(i,j,kbo(i,j)) = effsti*avn*ocetra(i,j,kbo(i,j),ifdust)   &
     &                *(dshagg+dsett)

      else
        wmass(i,j,kbo(i,j))=fTSFac*cellsink
        wnumb(i,j,kbo(i,j))=0.
        aggregate(i,j,kbo(i,j))=0.
        dustagg(i,j,kbo(i,j))=0.
        ocetra(i,j,kbo(i,j),inos)=0.
      endif ! avm
      endif ! alar1max
      endif ! pddpo
      enddo
      enddo
      
!
! Sampling timeseries-1 : sedimentation at specific positions
!
!      DO l=1,nts
!         i = its1(l)-p_ioff
!         j = jts1(l)-p_joff
!         IF(i<=1 .OR. i>=kpie .OR. j<=1 .OR. j>=kpje) CYCLE

! first depth 

!         if(k1ts1(l).gt.0) then
!             wphy = wmass(i,j,k1ts1(l))
	     
!             fphy=wphy*ocetra(i,j,k1ts1(l),iphy)  
!             fopa=wphy*ocetra(i,j,k1ts1(l),iopal)  
!             fdet=wphy*ocetra(i,j,k1ts1(l),idet)  
!             fcal=wphy*ocetra(i,j,k1ts1(l),icalc) 
!             ts1(its1fdet,l+1,lts1) = ts1(its1fdet,l+1,lts1) + fphy + fdet
!             ts1(its1fopa,l+1,lts1) = ts1(its1fopa,l+1,lts1) + fopa
!             ts1(its1fcal,l+1,lts1) = ts1(its1fcal,l+1,lts1) + fcal
!          else
!             ts1(its1fdet,l+1,lts1) = ts1(its1fdet,l+1,lts1) -9999.
!             ts1(its1fopa,l+1,lts1) = ts1(its1fopa,l+1,lts1) -9999.
!             ts1(its1fcal,l+1,lts1) = ts1(its1fcal,l+1,lts1) -9999.
!          endif          

! second depth 

!         if(k2ts1(l).gt.0) then
!             wphy = wmass(i,j,k2ts1(l))
	     
!             fphy=wphy*ocetra(i,j,k2ts1(l),iphy)  
!             fopa=wphy*ocetra(i,j,k2ts1(l),iopal)  
!             fdet=wphy*ocetra(i,j,k2ts1(l),idet)  
!             fcal=wphy*ocetra(i,j,k2ts1(l),icalc) 
!             ts1(its2fdet,l+1,lts1) = ts1(its2fdet,l+1,lts1) + fphy + fdet
!             ts1(its2fopa,l+1,lts1) = ts1(its2fopa,l+1,lts1) + fopa
!             ts1(its2fcal,l+1,lts1) = ts1(its2fcal,l+1,lts1) + fcal
!          else
!             ts1(its2fdet,l+1,lts1) = ts1(its2fdet,l+1,lts1) -9999.
!             ts1(its2fopa,l+1,lts1) = ts1(its2fopa,l+1,lts1) -9999.
!             ts1(its2fcal,l+1,lts1) = ts1(its2fcal,l+1,lts1) -9999.
!          endif
! third depth 
!         if(k3ts1(l).gt.0) then
!             wphy = wmass(i,j,k3ts1(l))
	     
!             fphy=wphy*ocetra(i,j,k3ts1(l),iphy)  
!             fopa=wphy*ocetra(i,j,k3ts1(l),iopal)  
!             fdet=wphy*ocetra(i,j,k3ts1(l),idet)  
!             fcal=wphy*ocetra(i,j,k2ts1(l),icalc) 
!             ts1(its3fdet,l+1,lts1) = ts1(its3fdet,l+1,lts1) + fphy + fdet
!             ts1(its3fopa,l+1,lts1) = ts1(its3fopa,l+1,lts1) + fopa
!             ts1(its3fcal,l+1,lts1) = ts1(its3fcal,l+1,lts1) + fcal
!          else
!             ts1(its3fdet,l+1,lts1) = ts1(its3fdet,l+1,lts1) -9999.
!             ts1(its3fopa,l+1,lts1) = ts1(its3fopa,l+1,lts1) -9999.
!             ts1(its3fcal,l+1,lts1) = ts1(its3fcal,l+1,lts1) -9999.
!         endif
!      ENDDO
!
! prepare output for bgcmean files
!
!      DO j=1,kpje
!      DO i=1,kpie
!         if(pddpo(i,j,kbo(i,j)).gt.0.5) then
! fluxes at 90 m
!         bgcm2d(i,j,jcoex90) = bgcm2d(i,j,jcoex90) +     &
!     &   (ocetra(i,j,n90depth,iphy)+ocetra(i,j,8,idet))*wmass(i,j,8)
!         bgcm2d(i,j,jopex90) = bgcm2d(i,j,jopex90) +     &
!     &   ocetra(i,j,n90depth,iopal)*wmass(i,j,8)
!         bgcm2d(i,j,jcaex90) = bgcm2d(i,j,jcaex90) +     &
!     &   ocetra(i,j,n90depth,icalc)*wmass(i,j,8)
! fluxes at about 1000 m
!         bgcm2d(i,j,jcoex1000) = bgcm2d(i,j,jcoex1000) + &
!     &   (ocetra(i,j,n1000depth,iphy)+ocetra(i,j,23,idet))*wmass(i,j,23)
!         bgcm2d(i,j,jopex1000) = bgcm2d(i,j,jopex1000) + &
!     &   ocetra(i,j,n1000depth,iopal)*wmass(i,j,23)
!         bgcm2d(i,j,jcaex1000) = bgcm2d(i,j,jcaex1000) + &
!     &   ocetra(i,j,n1000depth,icalc)*wmass(i,j,23)
!! fluxes at about 1950 m
!         bgcm2d(i,j,jcoex2000) = bgcm2d(i,j,jcoex2000) + &
!     &   (ocetra(i,j,n2000depth,iphy)+ocetra(i,j,29,idet))*wmass(i,j,15)
!         bgcm2d(i,j,jopex2000) = bgcm2d(i,j,jopex2000) + &
!     &   ocetra(i,j,n2000depth,iopal)*wmass(i,j,29)
!         bgcm2d(i,j,jcaex2000) = bgcm2d(i,j,jcaex2000) + &
!     &   ocetra(i,j,n2000depth,icalc)*wmass(i,j,29)
!         endif
!      ENDDO
!      ENDDO

!IK COMPUTE FLUXES FOR BOUNDARY CONDITION/BOTTOM LAYER

!js fluxes to sediment
	
      DO 36 j=1,kpje
      DO 36 i=1,kpie
         if(pddpo(i,j,kbo(i,j)).gt.1.e-6.and.omask(i,j).gt.0.5) then
         wphy = wmass(i,j,kbo(i,j))

         prorca(i,j) = ocetra(i,j,kbo(i,j),iphy)  *wphy                &
     &               + ocetra(i,j,kbo(i,j),idet)  *wphy
         prcaca(i,j) = ocetra(i,j,kbo(i,j),icalc) *wphy
         silpro(i,j) = ocetra(i,j,kbo(i,j),iopal) *wphy
         produs(i,j) = ocetra(i,j,kbo(i,j),ifdust)*dustsink            &
     &               + ocetra(i,j,kbo(i,j),iadust)*wphy    

!
! prepare output for bgcmean files 
!
!         bgct2d(i,j,jprorca)=bgct2d(i,j,jprorca) +prorca(i,j)
!         bgct2d(i,j,jprcaca)=bgct2d(i,j,jprcaca) +prcaca(i,j)
!         bgct2d(i,j,jsilpro)=bgct2d(i,j,jsilpro) +silpro(i,j)
!         bgct2d(i,j,jprodus)=bgct2d(i,j,jprodus) +produs(i,j)

        endif
36      CONTINUE

! COMPUTE FLUXES FOR LAYERS 2 TO kpke
      DO 2 K=kpke,2,-1
      DO 34 j=1,kpje
      DO 34 i=1,kpie
         if(omask(i,j).gt.0.5.and.pddpo(i,j,k).gt.1.e-6) then

         if(pddpo(i,j,k-1).gt.1.e-6) then
           kdonor=k-1
         else
           found=0
           do kinv=k-1,1,-1
             if (pddpo(i,j,kinv).gt.1.e-6.and.found.eq.0) then
               kdonor=kinv
               found=1
             endif 
           enddo
         endif 

! SINKING SPEED FOR UPPER LAYER
          wphyup = wmass(i,j,kdonor)    ! settling velocity of mass
          wnosup = wnumb(i,j,kdonor)    ! settling velocity of number of marine snow aggregates

! SINKING SPEED FOR ACTUAL LAYER
          wphy = wmass(i,j,k)
          wnos = wnumb(i,j,k)

! SUM FLUXES (compute new concentrations)
        ocetra(i,j,k,iphy) =ocetra(i,j,k,iphy) +                       &
     &      (ocetra(i,j,kdonor,iphy)*wphyup-ocetra(i,j,k,iphy)*wphy)   &
     &       *pdpio(i,j,k)
        ocetra(i,j,k,idet) =ocetra(i,j,k,idet) +                       &
     &      (ocetra(i,j,kdonor,idet)*wphyup-ocetra(i,j,k,idet)*wphy)   &
     &       *pdpio(i,j,k)
        ocetra(i,j,k,icalc) =ocetra(i,j,k,icalc) +                     &
     &      (ocetra(i,j,kdonor,icalc)*wphyup-ocetra(i,j,k,icalc)*wphy) &
     &       *pdpio(i,j,k)
        ocetra(i,j,k,iopal) =ocetra(i,j,k,iopal) +                     &
     &      (ocetra(i,j,kdonor,iopal)*wphyup-ocetra(i,j,k,iopal)*wphy) &
     &       *pdpio(i,j,k)
        ocetra(i,j,k,inos) =ocetra(i,j,k,inos) - aggregate(i,j,k) +    &
     &      (ocetra(i,j,kdonor,inos)*wnosup-ocetra(i,j,k,inos)*wnos)   &
     &       *pdpio(i,j,k)
! sinking of free dust and loss due to attachment to aggregated dust
        ocetra(i,j,k,ifdust) =ocetra(i,j,k,ifdust) - dustagg(i,j,k) +  &
     &      (ocetra(i,j,kdonor,ifdust)-ocetra(i,j,k,ifdust))*dustsink  &
     &       *pdpio(i,j,k)
! sinking of aggregated dust and gain due to attachment of free dust
        ocetra(i,j,k,iadust) =ocetra(i,j,k,iadust) + dustagg(i,j,k) +  &
     &      (ocetra(i,j,kdonor,iadust)*wphyup-ocetra(i,j,k,iadust)*wphy)&
     &       *pdpio(i,j,k)

      endif
34    CONTINUE
2     CONTINUE



!IK  EVALUATE FLUXES FOR FIRST LAYER

      DO 35 j=1,kpje
      DO 35 i=1,kpie
         if(omask(i,j).gt.0.5) then
! special depth for first layer:
         wphy = wmass(i,j,1)
         wnos = wnumb(i,j,1)

! SUM FLUXES
        ocetra(i,j,1,iphy) =ocetra(i,j,1,iphy)                          &
     &          -ocetra(i,j,1,iphy)*wphy*pdpio(i,j,1) 
        ocetra(i,j,1,idet) =ocetra(i,j,1,idet)                          &
     &          -ocetra(i,j,1,idet)*wphy*pdpio(i,j,1)
        ocetra(i,j,1,icalc) =ocetra(i,j,1,icalc)                        &
     &          -ocetra(i,j,1,icalc)*wphy*pdpio(i,j,1)
        ocetra(i,j,1,iopal) =ocetra(i,j,1,iopal)                        &
     &         -ocetra(i,j,1,iopal)*wphy*pdpio(i,j,1)
        ocetra(i,j,1,inos) =ocetra(i,j,1,inos) - aggregate(i,j,1)       &
     &         -ocetra(i,j,1,inos)*wnos*pdpio(i,j,1)
! sinking of free dust and loss due to attachment to aggregate dust
        ocetra(i,j,1,ifdust) =ocetra(i,j,1,ifdust) - dustagg(i,j,1)     &
     &         -ocetra(i,j,1,ifdust)*dustsink*pdpio(i,j,1)
! sinking of aggregated dust and gain due to attachment to aggregate dust
        ocetra(i,j,1,iadust) =ocetra(i,j,1,iadust) + dustagg(i,j,1)     &
     &         -ocetra(i,j,1,iadust)*wphy*pdpio(i,j,1)
          endif
35    CONTINUE 


#endif /*AGG*/

#ifndef AGG
#define sed_now
#ifdef sed_now
!
! implicit method:
! C(k,T+dt)=C(k,T) + (w*dt/ddpo(k))*(C(k-1,T+1)-C(k,T+1))
! --> 	    
! C(k,T+dt)=(ddpo(k)*C(k,T)+w*dt*C(k-1,T+dt))/(ddpo(k)+w*dt)
! sedimentation=w*dt*C(ks,T+dt)
!

       k=1
!$OMP PARALLEL DO
       DO j=1,kpje
       DO i=1,kpie
          IF(omask(i,j).gt.0.5.and.pddpo(i,j,k).GT.0.0) THEN

             ocetra(i,j,k,idet) =(ocetra(i,j,k,idet)*pddpo(i,j,k))/    &
     &                           (pddpo(i,j,k)+wpoc)
!proxies
#ifdef __c_isotopes
         ocetra(i,j,k,idet13) =(ocetra(i,j,k,idet13)*pddpo(i,j,k))/    &
     &                           (pddpo(i,j,k)+wpoc)
         ocetra(i,j,k,idet14) =(ocetra(i,j,k,idet14)*pddpo(i,j,k))/    &
     &                           (pddpo(i,j,k)+wpoc)
#endif

             ocetra(i,j,k,icalc)=(ocetra(i,j,k,icalc)*pddpo(i,j,k))/   &
     &                           (pddpo(i,j,k)+wcal)
!proxies
#ifdef __c_isotopes
         ocetra(i,j,k,icalc13)=(ocetra(i,j,k,icalc13)*pddpo(i,j,k))/   &
     &                           (pddpo(i,j,k)+wcal)
         ocetra(i,j,k,icalc14)=(ocetra(i,j,k,icalc14)*pddpo(i,j,k))/   &
     &                           (pddpo(i,j,k)+wcal)
#endif
 
             ocetra(i,j,k,iopal)=(ocetra(i,j,k,iopal)*pddpo(i,j,k))/   &
     &                           (pddpo(i,j,k)+wopal)    
 
             ocetra(i,j,k,ifdust)=(ocetra(i,j,k,ifdust)*pddpo(i,j,k))/ &
     &                            (pddpo(i,j,k)+wdust)    
     
          ENDIF
      enddo
      enddo
!$OMP END PARALLEL DO
    
!$OMP PARALLEL DO PRIVATE(kdonor)
      DO 12 j=1,kpje
      DO 12 i=1,kpie
        kdonor=1
        DO 10 k=2,kbo(i,j)
          IF(omask(i,j).gt.0.5) THEN
          IF(pddpo(i,j,k).GT.0.0) THEN

             ocetra(i,j,k,idet)=(ocetra(i,j,k  ,idet)*pddpo(i,j,k)     &
     &	                        +ocetra(i,j,kdonor,idet)*wpoc)/        &
     &                          (pddpo(i,j,k)+wpoc)
#ifdef __c_isotopes
             ocetra(i,j,k,idet13)=(ocetra(i,j,k  ,idet13)*pddpo(i,j,k) &
     &                          +ocetra(i,j,kdonor,idet13)*wpoc)/         &
     &                          (pddpo(i,j,k)+wpoc)
             ocetra(i,j,k,idet14)=(ocetra(i,j,k  ,idet14)*pddpo(i,j,k) &
     &                          +ocetra(i,j,kdonor,idet14)*wpoc)/         &
     &                          (pddpo(i,j,k)+wpoc)
#endif

             ocetra(i,j,k,icalc)=(ocetra(i,j,k  ,icalc)*pddpo(i,j,k)   &
     &	                         +ocetra(i,j,kdonor,icalc)*wcal)/      &
     &                           (pddpo(i,j,k)+wcal)
#ifdef __c_isotopes
         ocetra(i,j,k,icalc13)=(ocetra(i,j,k  ,icalc13)*pddpo(i,j,k)   &
     &                         +ocetra(i,j,kdonor,icalc13)*wcal)/         &
     &                           (pddpo(i,j,k)+wcal)
         ocetra(i,j,k,icalc14)=(ocetra(i,j,k  ,icalc14)*pddpo(i,j,k)   &
     &                         +ocetra(i,j,kdonor,icalc14)*wcal)/         &
     &                         (pddpo(i,j,k)+wcal)
#endif
 
             ocetra(i,j,k,iopal)=(ocetra(i,j,k  ,iopal)*pddpo(i,j,k)   &
     &	                         +ocetra(i,j,kdonor,iopal)*wopal)/     &
     &                           (pddpo(i,j,k)+wopal)        
 
             ocetra(i,j,k,ifdust)=(ocetra(i,j,k ,ifdust)*pddpo(i,j,k)  &
     &	                          +ocetra(i,j,kdonor,ifdust)*wdust)/   &
     &                            (pddpo(i,j,k)+wdust)        

          kdonor=k
          ELSE
             ocetra(i,j,k,idet)    =ocetra(i,j,kdonor,idet)
#ifdef __c_isotopes
            ocetra(i,j,k,idet13)=ocetra(i,j,kdonor,idet13)
            ocetra(i,j,k,idet14)=ocetra(i,j,kdonor,idet14)
#endif
             ocetra(i,j,k,icalc)   =ocetra(i,j,kdonor,icalc)
#ifdef __c_isotopes
            ocetra(i,j,k,icalc13)=ocetra(i,j,kdonor,icalc13)
            ocetra(i,j,k,icalc14)=ocetra(i,j,kdonor,icalc14)
#endif
             ocetra(i,j,k,iopal)   =ocetra(i,j,kdonor,iopal)
             ocetra(i,j,k,ifdust)  =ocetra(i,j,kdonor,ifdust)
          ENDIF
          ENDIF
10      CONTINUE
12    CONTINUE
!$OMP END PARALLEL DO
#ifdef PBGC_CK_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'in OCRPOD after sinking poc '
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask &
     &                  ,0)
!      call INVENTORY_BGC(kpie,kpje,kpke)
#endif 
!$OMP PARALLEL DO
      DO 33 j=1,kpje
      DO 33 i=1,kpie
        IF(omask(i,j).gt.0.5) THEN
         prorca(i,j)=ocetra(i,j,kbo(i,j),idet ) *wpoc
         prcaca(i,j)=ocetra(i,j,kbo(i,j),icalc) *wcal
#ifdef __c_isotopes
         pror13(i,j)=ocetra(i,j,kbo(i,j),idet13 )*wpoc
         prca13(i,j)=ocetra(i,j,kbo(i,j),icalc13)*wcal
         pror14(i,j)=ocetra(i,j,kbo(i,j),idet14 )*wpoc
         prca14(i,j)=ocetra(i,j,kbo(i,j),icalc14)*wcal
#endif
         silpro(i,j)=ocetra(i,j,kbo(i,j),iopal) *wopal
         produs(i,j)=ocetra(i,j,kbo(i,j),ifdust)*wdust

!
! Paddy: write output for bgcmean       
!
!         bgct2d(i,j,jprorca)=bgct2d(i,j,jprorca) + prorca(i,j)
!         bgct2d(i,j,jprcaca)=bgct2d(i,j,jprcaca) + prcaca(i,j)
!         bgct2d(i,j,jsilpro)=bgct2d(i,j,jsilpro) + silpro(i,j)
!         bgct2d(i,j,jprodus)=bgct2d(i,j,jprodus) + produs(i,j)

       ENDIF
33    CONTINUE
!$OMP END PARALLEL DO

#endif
#endif /* not AGG*/

!      pinv1=0.0
!      DO j=1,kpje
!      DO i=1,kpie
!        IF(omask(i,j).gt.0.0) THEN
!      DO k=1,kpke
!      dpinv=                                                      &
!     &   ocetra(i,j,k,idet)+ocetra(i,j,k,idoc)+ocetra(i,j,k,iphy) &
!     &  +ocetra(i,j,k,izoo)+ocetra(i,j,k,iphosph)                  
!      pinv1=pinv1+dpinv*pddpo(i,j,k)*pdlxp(i,j)*pdlyp(i,j)      
!      ENDDO
!      pinv1=pinv1+prorca(i,j)*pdlxp(i,j)*pdlyp(i,j)
!        ENDIF
!      ENDDO
!      ENDDO
!
!      WRITE(io_stdo_bgc,*) 'PINV sink',pinv0,pinv1

!
! Check maximum/minimum values on wet/dry cells.
!
!      IF( kchck .EQ. 1) CALL CHCK_BGC(io_stdo_bgc,icyclibgc,           &
!     &'Check values of ocean tracer at exit from SBR OCPROD :',        &
!     & kpie,kpje,kpke,pddpo)

      RETURN
      END
