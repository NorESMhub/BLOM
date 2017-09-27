      SUBROUTINE OCPROD(kpie,kpje,kpke,ptho,pddpo,                     &
     &                  pdlxp,pdlyp,pdpio,ptiestu,ptiestw,kplmon,omask)
!**********************************************************************
!
!**** *OCPROD* - .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,             *MPI-MaD, HH*    2010-04-01
!     J.Schwinger,           *GFI, UiB*       2013-04-22
!      - Corrected bug in light penetration formulation
!      - Cautious code clean-up
!     J.Tjiputra,            *UNI-RESEARCH*   2015-11-25
!      - Implemented natural DIC/ALK/CALC
!     I.Kriest,              *GEOMAR*         11.08.2016
!      - Modified stoichiometry for denitrification (affects NO3, N2, Alk)
!     J.Schwinger,           *UNI-RESEARCH*   2017-08-30
!      - Removed split of the layer that only partly falls into the 
!        euphotic zone. Loops are now calculated over 
!        (1) layers that are completely or partly in the euphotoc zone
!        (2) layers that do not lie within the euphotic zone.
!      - Moved the accumulation of global fields for output to routine
!        hamocc4bgc. The accumulation of local fields has been moved to
!        the end of this routine.
!
!     Purpose
!     -------
!     compute biological production, settling of debris, and related 
!     biogeochemistry
!
!
!
!**** Parameter list:
!     ---------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *ptho*    - potential temperature [deg C].
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st dimension) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (2nd dimension) [m].
!     *REAL*    *pdpio*   - inverse size of grid cell (3rd dimension)[m].
!     *REAL*    *ptiestu* - depth of layer centres
!     *REAL*    *ptiestw* - depth of layer interfaces (upper boundary)
!     *INTEGER* *kplmon*  - number of current month
!     *REAL*    *omask*   - land/ocean mask
!
!**********************************************************************

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
      REAL :: ptho (kpie,kpje,kpke)
      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: pdpio(kpie,kpje,kpke)
      REAL :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL :: ptiestu(kpie,kpje,kpke+1)
      REAL :: ptiestw(kpie,kpje,kpke+1)
      REAL :: abs_bgc(kpie,kpje,kpke)
      REAL :: omask(kpie,kpje)
      
      INTEGER :: i,j,k,l
      INTEGER :: is,kdonor
      INTEGER, PARAMETER :: nsinkmax=12
      REAL :: tco(nsinkmax),tcn(nsinkmax),q(nsinkmax)
      REAL :: dmsp1,dmsp2,dmsp3,dmsp4,dmsp5,dmsp6,dms_gamma,dms_ph
      REAL :: atten,avphy,avanut,avanfe,pho,xa,xn,ya,yn,phosy,          &
     &        avgra,grazing,avsil,avdic,graton,                         &
     &        gratpoc,grawa,bacfra,phymor,zoomor,excdoc,exud,           &
     &        export, delsil, delcar, sterph, sterzo, remin,            &
     &        docrem, opalrem, remin2o, aou,refra,pocrem,phyrem
      
      REAL :: zoothresh,phythresh
      REAL :: temfa,phofa                  ! temperature and irradiation factor for photosynthesis
      REAL :: dustinp
      REAL :: absorption
      REAL :: dmsprod,dms_bac,dms_uv 
      REAL :: dtr,dz
      REAL :: wpocd,wcald,wopald,dagg
#ifdef AGG
      REAL :: wmass(kpie,kpje,kpke)
      REAL :: wnumb(kpie,kpje,kpke)
      REAL :: aggregate(kpie,kpje,kpke)
      REAL :: dustagg(kpie,kpje,kpke)
      REAL :: avmass, avnos, anosloss     
      REAL :: zmornos, eps, e1,e2,e3,e4,es1,es3
      REAL :: TopM,TopF, snow,fshear,sagg1,sagg2,sagg4
      REAL :: sett_agg,shear_agg,effsti,dfirst,dshagg,dsett
      REAL :: wnos,wnosd
#endif 
      INTEGER, DIMENSION(kpie,kpje)   :: ind1,ind2
      REAL, DIMENSION(kpie,kpje,ddm)  :: wghts
      REAL, DIMENSION(kpie,kpje)      :: aux2d_dmsprod 
      REAL, DIMENSION(kpie,kpje)      :: aux2d_dms_bac
      REAL, DIMENSION(kpie,kpje)      :: aux2d_dms_uv 
      REAL, DIMENSION(kpie,kpje)      :: aux2d_export 
      REAL, DIMENSION(kpie,kpje)      :: aux2d_expoca 
      REAL, DIMENSION(kpie,kpje)      :: aux2d_exposi 
      REAL, DIMENSION(kpie,kpje)      :: aux2d_carflx0100 
      REAL, DIMENSION(kpie,kpje)      :: aux2d_carflx0500 
      REAL, DIMENSION(kpie,kpje)      :: aux2d_carflx1000
      REAL, DIMENSION(kpie,kpje)      :: aux2d_carflx2000 
      REAL, DIMENSION(kpie,kpje)      :: aux2d_carflx4000 
      REAL, DIMENSION(kpie,kpje)      :: aux2d_bsiflx0100
      REAL, DIMENSION(kpie,kpje)      :: aux2d_bsiflx0500
      REAL, DIMENSION(kpie,kpje)      :: aux2d_bsiflx1000
      REAL, DIMENSION(kpie,kpje)      :: aux2d_bsiflx2000
      REAL, DIMENSION(kpie,kpje)      :: aux2d_bsiflx4000
      REAL, DIMENSION(kpie,kpje)      :: aux2d_calflx0100
      REAL, DIMENSION(kpie,kpje)      :: aux2d_calflx0500
      REAL, DIMENSION(kpie,kpje)      :: aux2d_calflx1000
      REAL, DIMENSION(kpie,kpje)      :: aux2d_calflx2000
      REAL, DIMENSION(kpie,kpje)      :: aux2d_calflx4000
      REAL, DIMENSION(kpie,kpje)      :: aux2d_phosy 
      REAL, DIMENSION(kpie,kpje)      :: aux2d_dnit 
      REAL, DIMENSION(kpie,kpje,kpke) :: aux3d_phosy
#ifdef AGG
      REAL, DIMENSION(kpie,kpje,kpke) :: aux3d_eps
      REAL, DIMENSION(kpie,kpje,kpke) :: aux3d_asize
#endif 

      aux2d_dmsprod   (:,:)=0. 
      aux2d_dms_bac   (:,:)=0. 
      aux2d_dms_uv    (:,:)=0.  
      aux2d_export    (:,:)=0.  
      aux2d_expoca    (:,:)=0.  
      aux2d_exposi    (:,:)=0.  
      aux2d_carflx0100(:,:)=0.
      aux2d_carflx0500(:,:)=0.
      aux2d_carflx1000(:,:)=0.
      aux2d_carflx2000(:,:)=0.
      aux2d_carflx4000(:,:)=0.
      aux2d_bsiflx0100(:,:)=0.
      aux2d_bsiflx0500(:,:)=0.
      aux2d_bsiflx1000(:,:)=0.
      aux2d_bsiflx2000(:,:)=0.
      aux2d_bsiflx4000(:,:)=0.
      aux2d_calflx0100(:,:)=0.
      aux2d_calflx0500(:,:)=0.
      aux2d_calflx1000(:,:)=0.
      aux2d_calflx2000(:,:)=0.
      aux2d_calflx4000(:,:)=0.
      aux2d_phosy     (:,:)=0.  
      aux2d_dnit      (:,:)=0.  
      aux3d_phosy   (:,:,:)=0.
#ifdef AGG
      aux3d_eps     (:,:,:)=0.
      aux3d_asize   (:,:,:)=0.
#endif 

! Constant parameters
!
! parameter definition in BELEG_BGC.F

      dmsp6=dmspar(6)
      dmsp5=dmspar(5)
      dmsp4=dmspar(4)
      dmsp3=dmspar(3)
      dmsp2=dmspar(2)
      dmsp1=dmspar(1)
      dms_gamma=0.87


#ifdef PBGC_OCNP_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'beginning of OCRPOD '
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif   

! Calculate bottommost layer in the euphotic zone (kwrbioz)

      call calc_idepth(kpie,kpje,kpke,pddpo,ptiestw)


! Calculate swr absorption by water and phytoplankton

      abs_bgc(:,:,:)=0.
#ifdef FB_BGC_OCE     
      abs_oce(:,:,:)=0.
      abs_oce(:,:,1)=1.
#endif

!$OMP PARALLEL DO PRIVATE(absorption,atten,dz)
      DO j=1,kpje
      DO i=1,kpie

        IF(omask(i,j).GT.0.5) THEN

          absorption=1.

          vloop: DO k=1,kwrbioz(i,j)

          IF(pddpo(i,j,k).gt.0.0) THEN
          
          dz = pddpo(i,j,k)

          ! Average light intensity in layer k
          atten=atten_w+atten_c*max(0.,ocetra(i,j,k,iphy))
          abs_bgc(i,j,k)=((absorption/atten)*(1.-exp(-atten*dz)))/dz

#ifdef FB_BGC_OCE
          abs_oce(i,j,k)=abs_oce(i,j,k)*absorption
          if (k.eq.2) then
          abs_oce(i,j,2)=atten_f*absorption
          endif
#endif 

          ! Radiation intensity I_0 at the top of next layer
          absorption=absorption*exp(-atten*dz)

          ENDIF
          ENDDO vloop

        ENDIF ! omask GT 0.5

      ENDDO
      ENDDO
!$OMP END PARALLEL DO


! dust flux from the atmosphere to the surface layer; dust fields are
! monthly mean values (kg/m2/month - assume 30 days per month here)
! dissolved iron is a fixed fraction (typically 3.5%), and immediately released

!$OMP PARALLEL DO PRIVATE(dustinp)
      do j=1,kpje
      do i=1,kpie
        if(omask(i,j).gt.0.5) then
          dustinp=dusty(i,j,kplmon)/30.*dtb*pdpio(i,j,1)
          ocetra(i,j,1,ifdust)=ocetra(i,j,1,ifdust)+dustinp 
          ocetra(i,j,1,iiron)=ocetra(i,j,1,iiron)+dustinp*perc_diron 
        endif      
      enddo
      enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO                                                       &                
!$OMP&PRIVATE(avphy,avgra,avsil,avanut,avanfe,pho,xa,xn,phosy,          &
!$OMP&        ya,yn,grazing,graton,gratpoc,grawa,bacfra,phymor,         &  
!$OMP&        zoomor,excdoc,exud,export,delsil,delcar,dmsprod,          &
!$OMP&        dms_bac,dms_uv,dtr,phofa,temfa,zoothresh,dms_ph,dz)
      DO 1 j=1,kpje
      DO 1 i=1,kpie

      DO 100 K=1,kwrbioz(i,j)

      IF(pddpo(i,j,k).GT.dp_min.and.omask(i,j).gt.0.5) THEN


#ifdef AGG
        avmass = ocetra(i,j,k,iphy)+ocetra(i,j,k,idet)
#endif /*AGG*/

        phofa=pi_alpha*strahl(i,j)*abs_bgc(i,j,k) 
        temfa= 0.6* 1.066**ptho(i,j,k)                 
!taylor:     temfa= 0.6*(1. + 0.0639*ptho(i,j,k) *               &
!    &              (1. + 0.0639*ptho(i,j,k)/2. * (1. + 0.0639*ptho(i,j,k)/3.)))
        pho= dtb* phofa*temfa/sqrt(phofa**2 + temfa**2)     

        avphy=MAX(phytomi,ocetra(i,j,k,iphy))                   ! 'available' phytoplankton
        avgra=MAX(grami,ocetra(i,j,k,izoo))                     ! 'available' zooplankton
        avsil=MAX(0.,ocetra(i,j,k,isilica))
        avdic=MAX(0.,ocetra(i,j,k,isco212))
        avanut=MAX(0.,MIN(ocetra(i,j,k,iphosph),                        &
     &         rnoi*ocetra(i,j,k,iano3)))
        avanfe=MAX(0.,MIN(avanut,ocetra(i,j,k,iiron)/riron))
        xa=avanfe
        xn=xa/(1.+pho*avphy/(xa+bkphy)) 
        phosy=MAX(0.,xa-xn)
        phosy=MERGE(avdic/rcar,phosy,avdic.le.rcar*phosy)       ! limit phosy by available DIC
        ya=avphy+phosy
        yn=(ya+grazra*avgra*phytomi/(avphy+bkzoo))                      &
     &           /(1.+grazra*avgra/(avphy+bkzoo))
        grazing=MAX(0.,ya-yn)
        graton=epsher*(1.-zinges)*grazing
        gratpoc=(1.-epsher)*grazing
        grawa=epsher*zinges*grazing
        bacfra=remido*ocetra(i,j,k,idoc)

        phythresh=MAX(0.,(ocetra(i,j,k,iphy)-2.*phytomi))
        zoothresh=MAX(0.,(ocetra(i,j,k,izoo)-2.*grami))
        phymor=dyphy*phythresh
        exud=gammap*phythresh
        zoomor=spemor*zoothresh*zoothresh           ! *10 compared to linear in tropics (tinka)
        excdoc=gammaz*zoothresh                     ! excretion of doc by zooplankton  
        export= zoomor*(1.-ecan) + phymor + gratpoc ! ecan=.95, gratpoc= .2*grazing

#ifdef AGG	 
        delsil=MIN(ropal*phosy*avsil/(avsil+bkopal),0.5*avsil) 
	delcar=rcalc*MIN(calmax*phosy,(phosy-delsil/ropal))
#else
        delsil=MIN(ropal*export*avsil/(avsil+bkopal),0.5*avsil) 
        delcar=rcalc * export * bkopal/(avsil+bkopal)
#endif
!        dms_ph  = 1+(-log10(hi(i,j,1))-pi_ph(i,j,kplmon))*dms_gamma
        dms_ph  = 1. 
        dmsprod = (dmsp5*delsil+dmsp4*delcar)                           &
     &           *(1.+1./(ptho(i,j,k)+dmsp1)**2)*dms_ph        
        dms_bac = dmsp3*dtb*abs(ptho(i,j,k)+3.)*ocetra(i,j,k,idms)      &
     &             *(ocetra(i,j,k,idms)/(dmsp6+ocetra(i,j,k,idms)))     
        dms_uv  = dmsp2*dtb*phofa/pi_alpha*ocetra(i,j,k,idms)

        dtr=bacfra-phosy+graton+ecan*zoomor 

        ocetra(i,j,k,iphosph)=ocetra(i,j,k,iphosph)+dtr
        ocetra(i,j,k,iano3)=ocetra(i,j,k,iano3)+dtr*rnit
        ocetra(i,j,k,idet)=ocetra(i,j,k,idet)+export
        ocetra(i,j,k,idms)=ocetra(i,j,k,idms)+dmsprod-dms_bac-dms_uv
        ocetra(i,j,k,isco212)=ocetra(i,j,k,isco212)-delcar+rcar*dtr
        ocetra(i,j,k,ialkali)=ocetra(i,j,k,ialkali)-2.*delcar-(rnit+1)*dtr
        ocetra(i,j,k,ioxygen)=ocetra(i,j,k,ioxygen)-dtr*ro2ut
        ocetra(i,j,k,iphy)=ocetra(i,j,k,iphy)+phosy-grazing-phymor-exud
        ocetra(i,j,k,izoo)=ocetra(i,j,k,izoo)+grawa-excdoc-zoomor
        ocetra(i,j,k,idoc)=ocetra(i,j,k,idoc)-bacfra+excdoc+exud
        ocetra(i,j,k,icalc)=ocetra(i,j,k,icalc)+delcar
#ifdef natDIC
        ocetra(i,j,k,inatsco212)=ocetra(i,j,k,inatsco212)-delcar+rcar*dtr
        ocetra(i,j,k,inatalkali)=ocetra(i,j,k,inatalkali)-2.*delcar-(rnit+1)*dtr
        ocetra(i,j,k,inatcalc)=ocetra(i,j,k,inatcalc)+delcar
#endif
        ocetra(i,j,k,isilica)=ocetra(i,j,k,isilica)-delsil+dremopal*ocetra(i,j,k,iopal)
        ocetra(i,j,k,iopal)=ocetra(i,j,k,iopal)+delsil-dremopal*ocetra(i,j,k,iopal)
        ocetra(i,j,k,iiron)=ocetra(i,j,k,iiron)+dtr*riron              &
     &         -relaxfe*MAX(ocetra(i,j,k,iiron)-fesoly,0.)


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


! add up for total inventory and output
        dz = pddpo(i,j,k)

        expoor(i,j)=expoor(i,j)+dz*export*rcar
        expoca(i,j)=expoca(i,j)+dz*delcar
        exposi(i,j)=exposi(i,j)+dz*delsil

        aux2d_dmsprod(i,j)   = aux2d_dmsprod(i,j)+dmsprod*dz 
        aux2d_dms_bac(i,j)   = aux2d_dms_bac(i,j)+dms_bac*dz 
        aux2d_dms_uv(i,j)    = aux2d_dms_uv (i,j)+dms_uv*dz 
        aux2d_export(i,j)    = aux2d_export(i,j) +export*rcar*dz 
        aux2d_expoca(i,j)    = aux2d_expoca(i,j) +delcar*dz
        aux2d_exposi(i,j)    = aux2d_exposi(i,j) +delsil*dz 
        aux2d_phosy(i,j)     = aux2d_phosy(i,j)  +phosy*rcar*dz  ! primary production in kmol C m-2
        aux3d_phosy(i,j,k)   = phosy*rcar                        ! primary production in kmol C m-3

     
      ENDIF      ! pddpo(i,j,k).GT.dp_min
100   CONTINUE   ! kwrbioz
1     CONTINUE   ! kpie, kpje
!$OMP END PARALLEL DO

#ifdef PBGC_OCNP_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'in OCRPOD after 1st bio prod'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif 



!$OMP PARALLEL DO                                                       & 
!$OMP&PRIVATE(phythresh,zoothresh,sterph,sterzo,remin,opalrem,aou,      &
!$OMP&        refra,dms_bac,pocrem,docrem,phyrem,dz)  
      DO 201 j=1,kpje
      DO 201 i=1,kpie
         DO 20 k=kwrbioz(i,j)+1,kpke
            IF(pddpo(i,j,k).gt.dp_min.and.omask(i,j).GT.0.5) THEN

#ifdef AGG
            avmass=ocetra(i,j,k,iphy)+ocetra(i,j,k,idet)
#endif /*AGG*/	    
            phythresh=MAX(0.,(ocetra(i,j,k,iphy)-2.*phytomi))
            zoothresh=MAX(0.,(ocetra(i,j,k,izoo)-2.*grami))             
            sterph=0.5*dphymor*phythresh                                ! phytoplankton to detritus
            sterzo=dzoomor*zoothresh*zoothresh                          ! quadratic mortality

       	    ocetra(i,j,k,iphy)=ocetra(i,j,k,iphy)-sterph
       	    ocetra(i,j,k,izoo)=ocetra(i,j,k,izoo)-sterzo

            IF(ocetra(i,j,k,ioxygen).gt.5.e-8) THEN
               pocrem=MIN(drempoc*ocetra(i,j,k,idet),                   &
     &                   0.33*ocetra(i,j,k,ioxygen)/ro2ut)
               docrem=MIN(dremdoc*ocetra(i,j,k,idoc),                   &
     &                   0.33*ocetra(i,j,k,ioxygen)/ro2ut)
               phyrem=MIN(0.5*dphymor*phythresh,                        &
     &                   0.33*ocetra(i,j,k,ioxygen)/ro2ut)
            else
               pocrem=0.
               docrem=0.
               phyrem=0.
            endif 
	    
            ocetra(i,j,k,idet)=ocetra(i,j,k,idet)-pocrem+sterph+sterzo
            ocetra(i,j,k,idoc)=ocetra(i,j,k,idoc)-docrem
            ocetra(i,j,k,iphy)=ocetra(i,j,k,iphy)-phyrem

            remin=pocrem+docrem+phyrem

            ocetra(i,j,k,iphosph)=ocetra(i,j,k,iphosph)+remin
            ocetra(i,j,k,iano3)=ocetra(i,j,k,iano3)+remin*rnit
            ocetra(i,j,k,isco212)=ocetra(i,j,k,isco212)+rcar*remin
            ocetra(i,j,k,ialkali)=ocetra(i,j,k,ialkali)-(rnit+1)*remin
            ocetra(i,j,k,ioxygen)=ocetra(i,j,k,ioxygen)-ro2ut*remin
            ocetra(i,j,k,iiron)=ocetra(i,j,k,iiron)+remin*riron         & 
     &             -relaxfe*MAX(ocetra(i,j,k,iiron)-fesoly,0.)
#ifdef natDIC
            ocetra(i,j,k,inatsco212)=ocetra(i,j,k,inatsco212)+rcar*remin
            ocetra(i,j,k,inatalkali)=ocetra(i,j,k,inatalkali)-(rnit+1)*remin
#endif
!***********************************************************************
! as ragueneau (2000) notes, Si(OH)4sat is about 1000 umol, but
! Si(OH)4 varies only between 0-100 umol
! so the expression dremopal*(Si(OH)4sat-Si(OH)4) would change the 
! rate only from 0 to 100%     
!***********************************************************************
            opalrem=dremopal*0.1*(ptho(i,j,k)+3.)*ocetra(i,j,k,iopal)
            ocetra(i,j,k,iopal)=ocetra(i,j,k,iopal)-opalrem
            ocetra(i,j,k,isilica)=ocetra(i,j,k,isilica)+opalrem

!***********************************************************************
!           There is about 1.e4 O2 on 1 N2O molecule (Broeker&Peng)
!           refra : Tim Rixton, private communication
!***********************************************************************
            aou=satoxy(i,j,k)-ocetra(i,j,k,ioxygen)
            refra=1.+3.*(0.5+sign(0.5,aou-1.97e-4))
            dms_bac = dmsp3*dtb*abs(ptho(i,j,k)+3.)*ocetra(i,j,k,idms)   &
     &             *(ocetra(i,j,k,idms)/(dmsp6+ocetra(i,j,k,idms)))     
            ocetra(i,j,k,ian2o)=ocetra(i,j,k,ian2o)+remin*1.e-4*ro2ut*refra
            ocetra(i,j,k,igasnit)=ocetra(i,j,k,igasnit)-remin*1.e-4*ro2ut*refra
            ocetra(i,j,k,ioxygen)=ocetra(i,j,k,ioxygen)-remin*1.e-4*ro2ut*refra*0.5
            ocetra(i,j,k,idms)=ocetra(i,j,k,idms)-dms_bac

            dz = pddpo(i,j,k)
            aux2d_dms_bac(i,j)= aux2d_dms_bac(i,j)+dms_bac*dz

#ifdef AGG
!***********************************************************************
! loss of snow numbers due to remineralization of poc
! gain of snow numbers due to zooplankton mortality
! NOTE that remin is in kmol/m3. Thus divide by avmass (kmol/m3)
!***********************************************************************
           if(avmass.gt.0.) then  
              avnos = ocetra(i,j,k,inos)
              ocetra(i,j,k,inos) = ocetra(i,j,k,inos)-remin*avnos/avmass
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

#ifdef PBGC_OCNP_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'in OCRPOD after poc remin'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif


!$OMP PARALLEL DO                                                       &
!$OMP&PRIVATE(remin,remin2o,dz) 
       DO 30 j=1,kpje
       DO 30 i=1,kpie
         DO 30 k=kwrbioz(i,j)+1,kpke
         IF(omask(i,j).GT.0.5) THEN
         IF(ocetra(i,j,k,ioxygen).LT.5.e-7.and.pddpo(i,j,k).gt.dp_min) THEN


#ifdef AGG
           avmass = ocetra(i,j,k,iphy)+ocetra(i,j,k,idet)
#endif /*AGG*/	    
	    
           remin=0.05*drempoc*MIN(ocetra(i,j,k,idet),                   &
     &                        0.5*ocetra(i,j,k,iano3)/rdnit1)
           remin2o=dremn2o*MIN(ocetra(i,j,k,idet),                      &
     &	                 0.003*ocetra(i,j,k,ian2o)/rdn2o1)

           ocetra(i,j,k,ialkali)=ocetra(i,j,k,ialkali)+(rdnit1-1)*remin-remin2o
           ocetra(i,j,k,isco212)=ocetra(i,j,k,isco212)+rcar*(remin+remin2o)
#ifdef natDIC
           ocetra(i,j,k,inatalkali)=ocetra(i,j,k,inatalkali)+(rdnit1-1)*remin-remin2o
           ocetra(i,j,k,inatsco212)=ocetra(i,j,k,inatsco212)+rcar*(remin+remin2o)
#endif
           ocetra(i,j,k,idet)=ocetra(i,j,k,idet)-(remin+remin2o)
           ocetra(i,j,k,iphosph)=ocetra(i,j,k,iphosph)+(remin+remin2o)
           ocetra(i,j,k,iano3)=ocetra(i,j,k,iano3)-rdnit1*remin
           ocetra(i,j,k,igasnit)=ocetra(i,j,k,igasnit)+rdnit2*remin+rdn2o2*remin2o
           ocetra(i,j,k,iiron)=ocetra(i,j,k,iiron)+riron*(remin+remin2o)
           ocetra(i,j,k,ian2o)=ocetra(i,j,k,ian2o)-rdn2o1*remin2o

! nitrate loss through denitrification in kmol N m-2
           dz = pddpo(i,j,k)
           aux2d_dnit(i,j) = aux2d_dnit(i,j) + rdnit0*remin*dz 


#ifdef AGG
!***********************************************************************
! loss of snow numbers due to remineralization of poc
! NOTE that remin is in kmol/m3. Thus divide by avmass (kmol/m3)
!***********************************************************************
           if(avmass.gt.0.) then  
              avnos = ocetra(i,j,k,inos)
              ocetra(i,j,k,inos)=ocetra(i,j,k,inos)-(remin+remin2o)*avnos/avmass
           endif
#endif /*AGG*/

         ENDIF
         ENDIF
30    CONTINUE
!$OMP END PARALLEL DO


#ifdef PBGC_OCNP_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'in OCRPOD after remin n2o'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif 


!sulphate reduction   ! introduced 11.5.2007 to improve poc-remineralisation in the 
!                       oxygen minimum zone in the subsurface equatorial Pacific
!                       assumption of endless pool of SO4 (typical concentration are on the order of mmol/l)
!      js 02072007:    for other runs than current millenium (cosmos-setup) experiments this seems
!                      to cause trouble as phosphate concentrations are too high at the depth of the oxygen
!                      minimum in the equatorial pacific/atlantic
!                      does it make sense to check for oxygen and nitrate deficit?

!$OMP PARALLEL DO                                                      &
!$OMP&PRIVATE(remin) 
      DO 301 j=1,kpje
      DO 301 i=1,kpie
        DO 301 k=kwrbioz(i,j)+1,kpke
            IF(omask(i,j).gt.0.5.and.pddpo(i,j,k).gt.dp_min) then  
            IF(ocetra(i,j,k,ioxygen).lt.3.e-6.and.ocetra(i,j,k,iano3).lt.3.e-6) THEN

#ifdef AGG
               avmass = ocetra(i,j,k,iphy)+ocetra(i,j,k,idet)
#endif /*AGG*/
               remin=dremsul*ocetra(i,j,k,idet)

               ocetra(i,j,k,idet)=ocetra(i,j,k,idet)-remin
               ocetra(i,j,k,ialkali)=ocetra(i,j,k,ialkali)-(rnit+1)*remin
               ocetra(i,j,k,isco212)=ocetra(i,j,k,isco212)+rcar*remin
#ifdef natDIC
               ocetra(i,j,k,inatalkali)=ocetra(i,j,k,inatalkali)-(rnit+1)*remin
               ocetra(i,j,k,inatsco212)=ocetra(i,j,k,inatsco212)+rcar*remin
#endif
               ocetra(i,j,k,iphosph)=ocetra(i,j,k,iphosph)+remin
               ocetra(i,j,k,iano3)=ocetra(i,j,k,iano3)+rnit*remin
               ocetra(i,j,k,iiron)=ocetra(i,j,k,iiron)+riron*remin
   
#ifdef AGG
!***********************************************************************
! loss of snow numbers due to remineralization of poc
! NOTE that remin is in kmol/m3. Thus divide by avmass (kmol/m3)
!***********************************************************************
            if(avmass.gt.0.) then
               avnos = ocetra(i,j,k,inos)
               ocetra(i,j,k,inos) = ocetra(i,j,k,inos)-remin*avnos/avmass
            endif
#endif /*AGG*/

            ENDIF
            ENDIF
301    CONTINUE
!$OMP END PARALLEL DO
!    end sulphate reduction

#ifdef PBGC_OCNP_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'in OCRPOD after sulphate reduction '
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif 


#ifdef AGG

!**********************AGGREGATION***************************************
! General:
! Sinking speed, size distribution and aggregation are calculated 
! as in Kriest and Evans, 2000. I assume that opal and calcium carbonate
! sink at the same speed as P (mass).
!
! Sinking speed and aggregation: I assume that if there is no phosphorous mass,
! the sinking speed is the minimum sinking speed of aggregates. I further
! assume that then there are no particles, and that the rate of aggregation
! is 0. This scheme removes no P in the absence of P, but still opal and/or
! calcium carbonate.
! This could or should be changed, because silica as well as carbonate
! shell will add to the aggregate mass, and should be considered.
! Puh. Does anyone know functional relationships between
! size and Si or CaCO3? Perhaps on a later version, I have to
! take the relationship bewteen weight and size?
!
! Size distribution and resulting loss of marine snow aggregates due to  
! aggregation (aggregate(i,j,k)) and sinking speed of mass and numbers 
! (wmass(i,j,k) and wnumb(i,j,k) are calculated in a loop over 2-kpke. 
!
!************************************************************************

      wmass(:,:,:)     = 0.0
      wnumb(:,:,:)     = 0.0
      aggregate(:,:,:) = 0.0
      dustagg(:,:,:)   = 0.0

      do k=1,kpke
      do j=1,kpje
      do i=1,kpie

        IF(pddpo(i,j,k).gt.dp_min.and.omask(i,j).GT.0.5) THEN

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
          avmass = ocetra(i,j,k,iphy)+ocetra(i,j,k,idet)
          snow  = avmass*1.e+6

          if(avmass.gt.0.) then

! Set minimum particle number to nmldmin in the mixed layer. This is to prevent
! very small values of nos (and asscociated high sinking speed if there is mass)
! in high latitudes during winter        
          if ( k .le. kmle ) then  
            ocetra(i,j,k,inos) = MAX(nmldmin,ocetra(i,j,k,inos)) 
          endif

          ocetra(i,j,k,inos) = MAX(snow*pupper,ocetra(i,j,k,inos)) 
          ocetra(i,j,k,inos) = MIN(snow*plower,ocetra(i,j,k,inos)) 

          avnos = ocetra(i,j,k,inos)
          eps   = ((1.+ FractDim)*snow-avnos*cellmass)/(snow-avnos*cellmass)

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

! As a first step, assume that shear in the mixed layer is high and 
! zero below. 
          if ( k .le. kmle ) then  
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

      aggregate(i,j,k) = (shear_agg+sett_agg)*effsti*avnos*avnos

! dust aggregation:
! shear kernel:
      dfirst=dustd3+3.*dustd2*alar1+3.*dustd1*alar2+alar3
      dshagg=e1*fsh*(dfirst*TopF/e1-(                                  &
     &  (TopF-1.)/e1*dustd3+3.*(TopF*alar1-alow1)/e2*dustd2            &
     &   +3.*(TopF*alar2-alow2)/e3*dustd1+(TopF*alar3-alow3)/e4))

! settlement kernel:
      dsett=fse*dustd2*((e1+SinkExp*TopF*TSFac)/es1-dustsink/cellsink)
      
      dustagg(i,j,k) = effsti*avnos*ocetra(i,j,k,ifdust)               &
     &                *(dshagg+dsett)

        aux3d_eps(i,j,k)   = eps
        aux3d_asize(i,j,k) = snow/avnos/cellmass
  
      else 

        wmass(i,j,k)=cellsink
        wnumb(i,j,k)=0.
        aggregate(i,j,k)=0.
        dustagg(i,j,k)=0.
        ocetra(i,j,k,inos)=0.

        aux3d_eps(i,j,k)   = 1.
        aux3d_asize(i,j,k) = 0.

      endif ! avmass.gt.0

      endif ! pddpo>dp_min .and. omask>0.5
      enddo ! i=1,kpie
      enddo ! j=1,kpje
      enddo ! k=1,kpke

#endif /*AGG*/


!
! implicit method for sinking of particles:
! C(k,T+dt)=C(k,T) + (w*dt/ddpo(k))*(C(k-1,T+1)-C(k,T+1))
! --> 	    
! C(k,T+dt)=(ddpo(k)*C(k,T)+w*dt*C(k-1,T+dt))/(ddpo(k)+w*dt)
! sedimentation=w*dt*C(ks,T+dt)
!
!$OMP PARALLEL DO PRIVATE(kdonor,wpoc,wpocd,wcal,wcald,wopal,wopald,wnos,wnosd,dagg)
      DO j=1,kpje
      DO i=1,kpie

        tco(:) = 0.0 
        tcn(:) = 0.0

        IF(omask(i,j).gt.0.5) THEN

          kdonor=1
          DO k=1,kpke

          ! Sum up total column inventory before sinking scheme
          tco( 1) = tco( 1) + ocetra(i,j,k,idet  )*pddpo(i,j,k) 
          tco( 2) = tco( 2) + ocetra(i,j,k,icalc )*pddpo(i,j,k) 
#ifdef natDIC
          tco( 3) = tco( 3) + ocetra(i,j,k,inatcalc)*pddpo(i,j,k) 
#endif
          tco( 4) = tco( 4) + ocetra(i,j,k,iopal )*pddpo(i,j,k) 
          tco( 5) = tco( 5) + ocetra(i,j,k,ifdust)*pddpo(i,j,k) 
#if defined(AGG)
          tco( 6) = tco( 6) + ocetra(i,j,k,iphy  )*pddpo(i,j,k) 
          tco( 7) = tco( 7) + ocetra(i,j,k,inos  )*pddpo(i,j,k) 
          tco( 8) = tco( 8) + ocetra(i,j,k,iadust)*pddpo(i,j,k) 
#endif
#ifdef __c_isotopes
          tco( 9) = tco( 9) + ocetra(i,j,k,idet13 )*pddpo(i,j,k) 
          tco(10) = tco(10) + ocetra(i,j,k,idet14 )*pddpo(i,j,k) 
          tco(11) = tco(11) + ocetra(i,j,k,icalc13)*pddpo(i,j,k) 
          tco(12) = tco(12) + ocetra(i,j,k,icalc14)*pddpo(i,j,k) 
#endif

          IF(pddpo(i,j,k).GT.dp_min_sink) THEN

#if defined(AGG)
            wpoc   = wmass(i,j,k)
            wpocd  = wmass(i,j,kdonor)
            wcal   = wmass(i,j,k)
            wcald  = wmass(i,j,kdonor)
            wopal  = wmass(i,j,k)
            wopald = wmass(i,j,kdonor)
            wnos   = wnumb(i,j,k)
            wnosd  = wnumb(i,j,kdonor)
            wdust  = dustsink
            dagg   = dustagg(i,j,k)
#elif defined(WLIN)
            wpoc   = min(wmin+wlin*ptiestu(i,j,k),     wmax)
            wpocd  = min(wmin+wlin*ptiestu(i,j,kdonor),wmax)
            wcald  = wcal
            wopald = wopal
            dagg   = 0.0
#else
            wpocd  = wpoc
            wcald  = wcal
            wopald = wopal
            dagg   = 0.0
#endif

            IF( k==1 ) THEN
              wpocd  = 0.0
              wcald  = 0.0
              wopald = 0.0
#if defined(AGG)
              wnosd  = 0.0
#elif defined(WLIN)
              wpoc   = wmin
#endif
            ENDIF

            ocetra(i,j,k,idet)  =(ocetra(i,j,k     ,idet)*pddpo(i,j,k)    &
     &	                         +ocetra(i,j,kdonor,idet)*wpocd)/         &
     &                           (pddpo(i,j,k)+wpoc)
            ocetra(i,j,k,icalc) =(ocetra(i,j,k     ,icalc)*pddpo(i,j,k)   &
     &	                         +ocetra(i,j,kdonor,icalc)*wcald)/        &
     &                           (pddpo(i,j,k)+wcal)
#ifdef natDIC
            ocetra(i,j,k,inatcalc) =(ocetra(i,j,k,inatcalc)*pddpo(i,j,k)  &
     &	                         +ocetra(i,j,kdonor,inatcalc)*wcald)/     &
     &                           (pddpo(i,j,k)+wcal)
#endif
            ocetra(i,j,k,iopal) =(ocetra(i,j,k     ,iopal)*pddpo(i,j,k)   &
     &	                         +ocetra(i,j,kdonor,iopal)*wopald)/       &
     &                           (pddpo(i,j,k)+wopal)        
            ocetra(i,j,k,ifdust)=(ocetra(i,j,k     ,ifdust)*pddpo(i,j,k)  &
     &	                         +ocetra(i,j,kdonor,ifdust)*wdust)/       &
     &                           (pddpo(i,j,k)+wdust) - dagg        
#ifdef AGG
            ocetra(i,j,k,iphy)  =(ocetra(i,j,k     ,iphy)*pddpo(i,j,k)    &
     &	                         +ocetra(i,j,kdonor,iphy)*wpocd)/         &
     &                           (pddpo(i,j,k)+wpoc)
            ocetra(i,j,k,inos)  =(ocetra(i,j,k     ,inos)*pddpo(i,j,k)    &
     &	                         +ocetra(i,j,kdonor,inos)*wnosd)/         &
     &                           (pddpo(i,j,k)+wnos) - aggregate(i,j,k)
            ocetra(i,j,k,iadust)=(ocetra(i,j,k     ,iadust)*pddpo(i,j,k)  &
     &	                         +ocetra(i,j,kdonor,iadust)*wpocd)/       &
     &                           (pddpo(i,j,k)+wpoc)  + dagg 
#endif
#ifdef __c_isotopes
            ocetra(i,j,k,idet13)=(ocetra(i,j,k     ,idet13)*pddpo(i,j,k)  &
     &                           +ocetra(i,j,kdonor,idet13)*wpocd)/       &
     &                           (pddpo(i,j,k)+wpoc)
            ocetra(i,j,k,idet14)=(ocetra(i,j,k     ,idet14)*pddpo(i,j,k)  &
     &                           +ocetra(i,j,kdonor,idet14)*wpocd)/       &
     &                           (pddpo(i,j,k)+wpoc)
            ocetra(i,j,k,icalc13)=(ocetra(i,j,k     ,icalc13)*pddpo(i,j,k)&
     &                            +ocetra(i,j,kdonor,icalc13)*wcald)/     &
     &                            (pddpo(i,j,k)+wcal)
            ocetra(i,j,k,icalc14)=(ocetra(i,j,k     ,icalc14)*pddpo(i,j,k)&
     &                            +ocetra(i,j,kdonor,icalc14)*wcald)/     &
     &                            (pddpo(i,j,k)+wcal)
#endif
            kdonor=k

          ELSE IF( pddpo(i,j,k).GT.dp_min ) THEN

            ocetra(i,j,k,idet)   = ocetra(i,j,kdonor,idet)
            ocetra(i,j,k,icalc)  = ocetra(i,j,kdonor,icalc)
#ifdef natDIC
            ocetra(i,j,k,inatcalc)= ocetra(i,j,kdonor,inatcalc)
#endif
            ocetra(i,j,k,iopal)  = ocetra(i,j,kdonor,iopal)
            ocetra(i,j,k,ifdust) = ocetra(i,j,kdonor,ifdust)
#ifdef AGG
            ocetra(i,j,k,iphy)   = ocetra(i,j,kdonor,iphy)
            ocetra(i,j,k,inos)   = ocetra(i,j,kdonor,inos)
            ocetra(i,j,k,iadust) = ocetra(i,j,kdonor,iadust)
#endif
#ifdef __c_isotopes
            ocetra(i,j,k,idet13) = ocetra(i,j,kdonor,idet13)
            ocetra(i,j,k,idet14) = ocetra(i,j,kdonor,idet14)
            ocetra(i,j,k,icalc13)= ocetra(i,j,kdonor,icalc13)
            ocetra(i,j,k,icalc14)= ocetra(i,j,kdonor,icalc14)
#endif

          ENDIF  ! pddpo.GT.dp_min_sink

          ! Sum up total column inventory after sinking scheme
          ! flux to sediment added after kpke-loop
          IF( pddpo(i,j,k).GT.dp_min ) THEN
            tcn( 1) = tcn( 1) + ocetra(i,j,k,idet  )*pddpo(i,j,k) 
            tcn( 2) = tcn( 2) + ocetra(i,j,k,icalc )*pddpo(i,j,k) 
#ifdef natDIC
            tcn( 3) = tcn( 3) + ocetra(i,j,k,inatcalc)*pddpo(i,j,k) 
#endif
            tcn( 4) = tcn( 4) + ocetra(i,j,k,iopal )*pddpo(i,j,k) 
            tcn( 5) = tcn( 5) + ocetra(i,j,k,ifdust)*pddpo(i,j,k) 
#if defined(AGG)
            tcn( 6) = tcn( 6) + ocetra(i,j,k,iphy  )*pddpo(i,j,k) 
            tcn( 7) = tcn( 7) + ocetra(i,j,k,inos  )*pddpo(i,j,k) 
            tcn( 8) = tcn( 8) + ocetra(i,j,k,iadust)*pddpo(i,j,k) 
#endif
#ifdef __c_isotopes
            tcn( 9) = tcn( 9) + ocetra(i,j,k,idet13 )*pddpo(i,j,k) 
            tcn(10) = tcn(10) + ocetra(i,j,k,idet14 )*pddpo(i,j,k) 
            tcn(11) = tcn(11) + ocetra(i,j,k,icalc13)*pddpo(i,j,k) 
            tcn(12) = tcn(12) + ocetra(i,j,k,icalc14)*pddpo(i,j,k) 
#endif
          ENDIF

          ENDDO  ! loop k=1,kpke


          ! Add fluxes to sediment to new total column inventory
          tcn( 1) = tcn( 1) + ocetra(i,j,kdonor,idet  )*wpoc
          tcn( 2) = tcn( 2) + ocetra(i,j,kdonor,icalc )*wcal
#ifdef natDIC
          tcn( 3) = tcn( 3) + ocetra(i,j,kdonor,inatcalc)*wcal
#endif
          tcn( 4) = tcn( 4) + ocetra(i,j,kdonor,iopal )*wopal
          tcn( 5) = tcn( 5) + ocetra(i,j,kdonor,ifdust)*wdust
#if defined(AGG)
          tcn( 6) = tcn( 6) + ocetra(i,j,kdonor,iphy  )*wpoc
          tcn( 7) = tcn( 7) + ocetra(i,j,kdonor,inos  )*wnos
          tcn( 8) = tcn( 8) + ocetra(i,j,kdonor,iadust)*wpoc 
#endif
#ifdef __c_isotopes
          tcn( 9) = tcn( 9) + ocetra(i,j,kdonor,idet13 )*wpoc
          tcn(10) = tcn(10) + ocetra(i,j,kdonor,idet14 )*wpoc
          tcn(11) = tcn(11) + ocetra(i,j,kdonor,icalc13)*wcal
          tcn(12) = tcn(12) + ocetra(i,j,kdonor,icalc14)*wcal
#endif

          ! Do columnwise multiplicative mass conservation correction
          q(:)=1.0
          DO is=1,nsinkmax
            IF( tco(is)>1.e-12 .and. tcn(is)>1.e-12 ) q(is)=tco(is)/tcn(is)
          ENDDO
          DO k=1,kpke
            IF( pddpo(i,j,k).GT.dp_min ) THEN
              ocetra(i,j,k,idet  )=ocetra(i,j,k,idet  )*q(1)
              ocetra(i,j,k,icalc )=ocetra(i,j,k,icalc )*q(2)
#ifdef natDIC
              ocetra(i,j,k,inatcalc)=ocetra(i,j,k,inatcalc)*q(3)
#endif
              ocetra(i,j,k,iopal )=ocetra(i,j,k,iopal )*q(4)
              ocetra(i,j,k,ifdust)=ocetra(i,j,k,ifdust)*q(5)
#if defined(AGG)
              ocetra(i,j,k,iphy  )=ocetra(i,j,k,iphy  )*q(6)
              ocetra(i,j,k,inos  )=ocetra(i,j,k,inos  )*q(7)
              ocetra(i,j,k,iadust)=ocetra(i,j,k,iadust)*q(8)
#endif
#ifdef __c_isotopes
              ocetra(i,j,k,idet13 )=ocetra(i,j,k,idet13 )*q(9)
              ocetra(i,j,k,idet14 )=ocetra(i,j,k,idet14 )*q(10)
              ocetra(i,j,k,icalc13)=ocetra(i,j,k,icalc13)*q(11)
              ocetra(i,j,k,icalc14)=ocetra(i,j,k,icalc14)*q(12)
#endif
            ENDIF
          ENDDO

! Fluxes to sediment, layers thinner than dp_min_sink are ignored.
! Note that kdonor=kbo(i,j) by definition since kbo is the lowermost
! layer thicker than dp_min_sink.
#if defined(AGG)
          prorca(i,j)=ocetra(i,j,kdonor,iphy  )*wpoc                    &
     &              + ocetra(i,j,kdonor,idet  )*wpoc
          prcaca(i,j)=ocetra(i,j,kdonor,icalc )*wcal
          silpro(i,j)=ocetra(i,j,kdonor,iopal )*wopal
          produs(i,j)=ocetra(i,j,kdonor,ifdust)*wdust                   &
     &              + ocetra(i,j,kdonor,iadust)*wpoc    
#else
          prorca(i,j)=ocetra(i,j,kdonor,idet  )*wpoc
          prcaca(i,j)=ocetra(i,j,kdonor,icalc )*wcal
          silpro(i,j)=ocetra(i,j,kdonor,iopal )*wopal
          produs(i,j)=ocetra(i,j,kdonor,ifdust)*wdust
#endif
#ifdef __c_isotopes
          pror13(i,j)=ocetra(i,j,kdonor,idet13 )*wpoc
          prca13(i,j)=ocetra(i,j,kdonor,icalc13)*wcal
          pror14(i,j)=ocetra(i,j,kdonor,idet14 )*wpoc
          prca14(i,j)=ocetra(i,j,kdonor,icalc14)*wcal
#endif


        ENDIF  ! omask.gt.0.5
      ENDDO    ! loop i=1,kpie
      ENDDO    ! loop j=1,kpje
!$OMP END PARALLEL DO




! Calculate mass sinking flux for carbon, opal and calcium carbonate
! through the 100 m, 500 m, 1000 m, and 2000 m depth surfaces. These 
! fluxes are intentionally calculated using values at the NEW timelevel
! to be fully consistent with the implicit sinking scheme
      IF( domassfluxes ) THEN

!$OMP PARALLEL DO PRIVATE(wpoc,wcal,wopal)
      DO j=1,kpje
      DO i=1,kpie
        IF(omask(i,j).gt.0.5) THEN

          ! 100 m
          k = k0100(i,j)
          if(k>0) then
#if defined(AGG)
            wpoc  = wmass(i,j,k)
            wcal  = wmass(i,j,k)
            wopal = wmass(i,j,k)
#elif defined(WLIN)
            wpoc  = min(wmin+wlin*ptiestu(i,j,k), wmax)
#endif

#if defined(AGG)
            aux2d_carflx0100(i,j) = (ocetra(i,j,k,idet)+ocetra(i,j,k,iphy))*rcar*wpoc
#else
            aux2d_carflx0100(i,j) = ocetra(i,j,k,idet)*rcar*wpoc
#endif
            aux2d_bsiflx0100(i,j) = ocetra(i,j,k,iopal)*wopal
            aux2d_calflx0100(i,j) = ocetra(i,j,k,icalc)*wcal
          endif

          ! 500 m
          k = k0500(i,j)
          if(k>0) then
#if defined(AGG)
            wpoc  = wmass(i,j,k)
            wcal  = wmass(i,j,k)
            wopal = wmass(i,j,k)
#elif defined(WLIN)
            wpoc  = min(wmin+wlin*ptiestu(i,j,k), wmax)
#endif

#if defined(AGG)
            aux2d_carflx0500(i,j) = (ocetra(i,j,k,idet)+ocetra(i,j,k,iphy))*rcar*wpoc
#else
            aux2d_carflx0500(i,j) = ocetra(i,j,k,idet)*rcar*wpoc
#endif
            aux2d_bsiflx0500(i,j) = ocetra(i,j,k,iopal)*wopal
            aux2d_calflx0500(i,j) = ocetra(i,j,k,icalc)*wcal
          endif

          ! 1000 m
          k = k1000(i,j)
          if(k>0) then
#if defined(AGG)
            wpoc  = wmass(i,j,k)
            wcal  = wmass(i,j,k)
            wopal = wmass(i,j,k)
#elif defined(WLIN)
            wpoc  = min(wmin+wlin*ptiestu(i,j,k), wmax)
#endif

#if defined(AGG)
            aux2d_carflx1000(i,j) = (ocetra(i,j,k,idet)+ocetra(i,j,k,iphy))*rcar*wpoc
#else
            aux2d_carflx1000(i,j) = ocetra(i,j,k,idet)*rcar*wpoc
#endif
            aux2d_bsiflx1000(i,j) = ocetra(i,j,k,iopal)*wopal
            aux2d_calflx1000(i,j) = ocetra(i,j,k,icalc)*wcal
          endif

          ! 2000 m
          k = k2000(i,j)
          if(k>0) then
#if defined(AGG)
            wpoc  = wmass(i,j,k)
            wcal  = wmass(i,j,k)
            wopal = wmass(i,j,k)
#elif defined(WLIN)
            wpoc  = min(wmin+wlin*ptiestu(i,j,k), wmax)
#endif

#if defined(AGG)
            aux2d_carflx2000(i,j) = (ocetra(i,j,k,idet)+ocetra(i,j,k,iphy))*rcar*wpoc
#else
            aux2d_carflx2000(i,j) = ocetra(i,j,k,idet)*rcar*wpoc
#endif
            aux2d_bsiflx2000(i,j) = ocetra(i,j,k,iopal)*wopal
            aux2d_calflx2000(i,j) = ocetra(i,j,k,icalc)*wcal
          endif

          ! 4000 m
          k = k4000(i,j)
          if(k>0) then
#if defined(AGG)
            wpoc  = wmass(i,j,k)
            wcal  = wmass(i,j,k)
            wopal = wmass(i,j,k)
#elif defined(WLIN)
            wpoc  = min(wmin+wlin*ptiestu(i,j,k), wmax)
#endif

#if defined(AGG)
            aux2d_carflx4000(i,j) = (ocetra(i,j,k,idet)+ocetra(i,j,k,iphy))*rcar*wpoc
#else
            aux2d_carflx4000(i,j) = ocetra(i,j,k,idet)*rcar*wpoc
#endif
            aux2d_bsiflx4000(i,j) = ocetra(i,j,k,iopal)*wopal
            aux2d_calflx4000(i,j) = ocetra(i,j,k,icalc)*wcal
          endif

        ENDIF ! omask > 0.5
      enddo
      enddo
!$OMP END PARALLEL DO

      ENDIF ! domassfluxes


!     Accumulate local 2d diagnostics 
      call accsrf(jdmsprod,aux2d_dmsprod,omask,0)    
      call accsrf(jdms_uv,aux2d_dms_uv,omask,0)     
      call accsrf(jdms_bac,aux2d_dms_bac,omask,0) 
      call accsrf(jexport,aux2d_export,omask,0)      
      call accsrf(jexpoca,aux2d_expoca,omask,0)     
      call accsrf(jexposi,aux2d_exposi,omask,0)     
      call accsrf(jintphosy,aux2d_phosy,omask,0)     
      call accsrf(jintdnit,aux2d_dnit,omask,0) 

!     Accumulate local layer diagnostics
      call acclyr(jphosy,aux3d_phosy,pddpo,1)
#ifdef AGG
      call acclyr(jwphy, wmass/dtb,  pddpo,1)
      call acclyr(jwnos, wnumb/dtb,  pddpo,1)
      call acclyr(jeps,  aux3d_eps,  pddpo,1)
      call acclyr(jasize,aux3d_asize,pddpo,1)
#endif     

!     Accumulate local level diagnostics
      IF (SUM(jlvlphosy+jlvlwphy+jlvlwnos+jlvleps+jlvlasize).NE.0) THEN
        DO k=1,kpke
          call bgczlv(pddpo,k,ind1,ind2,wghts)
          call acclvl(jlvlphosy,aux3d_phosy,k,ind1,ind2,wghts)
#ifdef AGG
          call acclvl(jlvlwphy, wmass/dtb,  k,ind1,ind2,wghts)
          call acclvl(jlvlwnos, wnumb/dtb,  k,ind1,ind2,wghts)
          call acclvl(jlvleps,  aux3d_eps,  k,ind1,ind2,wghts)
          call acclvl(jlvlasize,aux3d_asize,k,ind1,ind2,wghts)
#endif     
        ENDDO
      ENDIF

!     Accumulate the diagnostic mass sinking field 
      IF( domassfluxes ) THEN
        call accsrf(jcarflx0100,aux2d_carflx0100,omask,0)    
        call accsrf(jbsiflx0100,aux2d_bsiflx0100,omask,0)    
        call accsrf(jcalflx0100,aux2d_calflx0100,omask,0)    
        call accsrf(jcarflx0500,aux2d_carflx0500,omask,0)    
        call accsrf(jbsiflx0500,aux2d_bsiflx0500,omask,0)    
        call accsrf(jcalflx0500,aux2d_calflx0500,omask,0)    
        call accsrf(jcarflx1000,aux2d_carflx1000,omask,0)    
        call accsrf(jbsiflx1000,aux2d_bsiflx1000,omask,0)    
        call accsrf(jcalflx1000,aux2d_calflx1000,omask,0)    
        call accsrf(jcarflx2000,aux2d_carflx2000,omask,0)    
        call accsrf(jbsiflx2000,aux2d_bsiflx2000,omask,0)    
        call accsrf(jcalflx2000,aux2d_calflx2000,omask,0)    
        call accsrf(jcarflx4000,aux2d_carflx4000,omask,0)    
        call accsrf(jbsiflx4000,aux2d_bsiflx4000,omask,0)    
        call accsrf(jcalflx4000,aux2d_calflx4000,omask,0)    
        call accsrf(jcarflx_bot,prorca*rcar,     omask,0)    
        call accsrf(jbsiflx_bot,silpro,          omask,0)    
        call accsrf(jcalflx_bot,prcaca,          omask,0)    
      ENDIF ! domassfluxes


#ifdef PBGC_OCNP_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'in OCRPOD after sinking poc '
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif 


      RETURN
      END
