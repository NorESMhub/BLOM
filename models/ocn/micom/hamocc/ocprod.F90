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
!     S.Legutke,             *MPI-MaD, HH*    10.04.01
!     J.Schwinger            *GFI, UiB*       2013-04-22
!      - Corrected bug in light penetration formulation
!      - Cautious code clean-up
!
!
!     Purpose
!     -------
!     compute biological production, settling of debris, and related 
!     biogeochemistry
!
!     Note: 
!     _ant fields are natural PLUS anthropogenic (not anthropogenic only!!!)
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
      INTEGER :: i,j,k,l
      INTEGER :: kinv,kdonor
      REAL :: ptho (kpie,kpje,kpke)
      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: pdpio(kpie,kpje,kpke)
      REAL :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL :: ptiestu(kpie,kpje,kpke+1)
      REAL :: ptiestw(kpie,kpje,kpke+1)
      REAL :: abs_bgc(kpie,kpje,kpke)
      REAL :: omask(kpie,kpje)
      
      REAL :: dmsp1,dmsp2,dmsp3,dmsp4,dmsp5,dmsp6,dms_gamma,dms_ph
      REAL :: atten,avphy,avanut,avanfe,pho,xa,xn,ya,yn,phosy,         &
     &        avgra,grazing,avsil,graton,                              &
     &        gratpoc,grawa,bacfra,phymor,zoomor,excdoc,exud,          &
     &        export, delsil, delcar, sterph, sterzo, remin,           &
     &        docrem, opalrem, remin2o, aou,refra,pocrem,phyrem
      
      REAL :: zoothresh,phythresh
      REAL :: temfa,phofa                  ! temperature and irradiation factor for photosynthesis
      REAL :: dustinp
      REAL :: absorption
      REAL :: dmsprod,dms_bac,dms_uv 
      REAL :: bdp,dtr 
      REAL :: detref, detrl
      REAL :: dz_light
      REAL :: wpocd,wcald,wopald,dagg
#ifdef __c_isotopes
      REAL :: rem13,rem14
      REAL :: rl13, rl14
      REAL :: rocean13, rocean14, flui13, flui14, phyto_prev
      REAL :: aDIC13, aDIC14
      REAL, DIMENSION(kpie,kpje,kpke) :: d13C, d14C, dd14C
#endif
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

!$OMP PARALLEL DO PRIVATE(absorption,atten,dz_light)
      DO j=1,kpje
      DO i=1,kpie

        IF(omask(i,j).GT.0.5) THEN

          absorption=1.

          vloop: DO k=1,kwrbioz(i,j)

          IF(pddpo(i,j,k).gt.0.0) THEN

          if( ptiestw(i,j,k+1) >= dp_ez .and. ptiestw(i,j,k) < dp_ez ) then
            dz_light = dp_ez-ptiestw(i,j,k)
          else
            dz_light = pddpo(i,j,k)
          end if

          ! Average light intensity in layer k
          atten=atten_w+atten_c*max(0.,ocetra(i,j,k,iphy))
          abs_bgc(i,j,k)=((absorption/atten)*                           &
     &                   (1.-exp(-atten*dz_light)))/dz_light

#ifdef FB_BGC_OCE
          abs_oce(i,j,k)=abs_oce(i,j,k)*absorption
          if (k.eq.2) then
          abs_oce(i,j,2)=atten_f*absorption
          endif
#endif 

          ! Radiation intensity I_0 at the top of next layer
          absorption=absorption*exp(-atten*dz_light)

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


! accumulate srf diagnostics
      call accsrf(jsrfphosph,ocetra(1,1,1,iphosph),omask,0)
      call accsrf(jsrfoxygen,ocetra(1,1,1,ioxygen),omask,0)
      call accsrf(jsrfiron,ocetra(1,1,1,iiron),omask,0)
      call accsrf(jsrfano3,ocetra(1,1,1,iano3),omask,0)
      call accsrf(jsrfalkali,ocetra(1,1,1,ialkali),omask,0)
      call accsrf(jsrfsilica,ocetra(1,1,1,isilica),omask,0)
      call accsrf(jsrfdic,ocetra(1,1,1,isco212),omask,0)

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
! Delta notation delta 13C & capital delta 14C (dd14C) [promille]

!$OMP PARALLEL DO PRIVATE(d14C) 
      DO j=1,kpje
      DO i=1,kpie
        IF(omask(i,j).GT.0.5) THEN
          DO k=1,kpke
          
            IF(pddpo(i,j,k).GT.dp_min)
            ! 1. Calculate absolute model values for all isotope parameters using 
            ! the calibration factor [kmol/m3]
            aDIC13  = ocetra(i,j,k,isco213) * factor_13c
            aDIC14  = ocetra(i,j,k,isco214) * factor_14c

            ! 2. Calculate the calibrated absolute values for the ocean in [promille]
            ! Realize that model isco212 actually is 12C + 13C; therefore subtraction 12C - 13C
            d13C(i,j,k)  = ((aDIC13/((ocetra(i,j,k,isco212)-aDIC13)+1.e-25)/PDB)-1.)*1000.

            ! For 14C two steps are needed: first to d14C [promille], then to dd14C [promille]
            d14C(i,j,k)  = ((aDIC14(i,j,k)/ocetra(i,j,k,isco212)   / ref14c) -1.) * 1000.
            dd14C(i,j,k) = d14C(i,j,k) - 2.*(d13C(i,j,k) + 25.) * (1.+ d14C(i,j,k)/1000.)
            ENDIF

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



#define bioprod
#ifdef bioprod


!$OMP PARALLEL DO                                               &                
!$OMP&PRIVATE(avphy,avgra,avsil,avanut,avanfe,pho,xa,xn,phosy,  &
!$OMP&        ya,yn,grazing,graton,gratpoc,grawa,bacfra,phymor, &  
!$OMP&        zoomor,excdoc,exud,export,delsil,delcar,dmsprod,  &
!$OMP&        dms_bac,dms_uv,bdp,dtr,rocean13,rocean14,flui13,  &
!$OMP&        flui14,dp_ez,phofa,temfa,zoothresh,dms_ph)
      DO 1 j=1,kpje
      DO 1 i=1,kpie

      DO 100 K=1,kwrbioz(i,j)

      IF(pddpo(i,j,k).GT.dp_min.and.omask(i,j).gt.0.5) THEN


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

         flui13=max(rocean13-1.,0.)        				     ! assumes rocean >1 in euphotic layer 
         flui14=max(rocean14-1.,0.)        				     ! 
#endif

         phofa=pi_alpha*strahl(i,j)*abs_bgc(i,j,k) 
         temfa= 0.6* 1.066**ptho(i,j,k)                 
!taylor:      temfa= 0.6*(1. + 0.0639*ptho(i,j,k) *               &
!    &               (1. + 0.0639*ptho(i,j,k)/2. * (1. + 0.0639*ptho(i,j,k)/3.)))
         pho= dtb* phofa*temfa/sqrt(phofa**2 + temfa**2)     

         avphy=MAX(0.,ocetra(i,j,k,iphy))  
         avgra=MAX(0.,ocetra(i,j,k,izoo))  
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

#ifdef AGG	 
         delsil=MIN(ropal*phosy*avsil/(avsil+bkopal),0.5*avsil) 
	 delcar=rcalc*MIN(calmax*phosy,(phosy-delsil/ropal))
#else
         delsil=MIN(ropal*export*avsil/(avsil+bkopal),0.5*avsil) 
         delcar=rcalc * export * bkopal/(avsil+bkopal)
#endif
!         dms_ph  = 1+(-log10(hi(i,j,1))-pi_ph(i,j,kplmon))*dms_gamma
         dms_ph  = 1. 
         dmsprod = (dmsp5*delsil+dmsp4*delcar)                        &
     &            *(1.+1./(ptho(i,j,k)+dmsp1)**2)*dms_ph        
         dms_bac = dmsp3*dtb*abs(ptho(i,j,k)+3.)*ocetra(i,j,k,idms)   &
     &             *(ocetra(i,j,k,idms)/(dmsp6+ocetra(i,j,k,idms)))     
         dms_uv  = dmsp2*dtb*phofa/pi_alpha*ocetra(i,j,k,idms)

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
     &   ((ocetra(i,j,k,idet14)+rcar*export*bifr14)*bdp               & ! rcar: because idet in P units and idet13 in C units
     &    +ocetra(i,j,k,idet14)*(pddpo(i,j,k)-bdp))/pddpo(i,j,k) 

         ocetra(i,j,k,idet13)=                                        &
     &   ((ocetra(i,j,k,idet13)+rcar*export*bifr13)*bdp              &
     &    +ocetra(i,j,k,idet13)*(pddpo(i,j,k)-bdp))/pddpo(i,j,k)   

         ocetra(i,j,k,isco214)=                                        &
     &   ((ocetra(i,j,k,isco214)-rcar*export*bifr14)*bdp    &
     &    +ocetra(i,j,k,isco214)*(pddpo(i,j,k)-bdp))/pddpo(i,j,k)   

         ocetra(i,j,k,isco213)=                                        &
     &   ((ocetra(i,j,k,isco213)-rcar*export*bifr13)*bdp    &
     &    +ocetra(i,j,k,isco213)*(pddpo(i,j,k)-bdp))/pddpo(i,j,k)  

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
! Update DIC isotope fields and include fractionation during CaCO3 and OM formation
         ocetra(i,j,k,isco213)=                                        &
     &   ((ocetra(i,j,k,isco213)-delcar*rocean13*frac_caco313 +        &
     &    rcar*dtr)*bdp +ocetra(i,j,k,isco213)  *(pddpo(i,j,k)-bdp))/pddpo(i,j,k) 

         ocetra(i,j,k,isco214)=                                        &
     &   ((ocetra(i,j,k,isco214)-delcar*rocean14*frac_caco314 +        &
     &    rcar*dtr)*bdp+ocetra(i,j,k,isco214)   *(pddpo(i,j,k)-bdp))/pddpo(i,j,k) 
#endif /*__c_isotopes/

!#ifdef __c_isotopes
!         ocetra(i,j,k,isco213)=                                        &
!     &   ((ocetra(i,j,k,isco213)-delcar*rocean13                       &
!     &   +bifr13*rcar*( bacfra - phosy + graton + ecan*zoomor))*bdp     & ! bifr=0.98
!     &   +ocetra(i,j,k,isco213)     *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)

!          ocetra(i,j,k,isco214)=                                       &
!     &    ((ocetra(i,j,k,isco214)-delcar*rocean14)*bdp                 &
!     &    +ocetra(i,j,k,isco214)     *(pddpo(i,j,k)-bdp))/pddpo(i,j,k) 
! js: for efficency below line (which should in principle be there) is neglected (additional tracer field would be needed 
!     to account for radioactive decay of 14C in particles)
!     &      + bifr14*rcar*( bacfra - phosy + graton + ecan*zoomor)
!#endif /*__c_isotopes*/

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

         ocetra(i,j,k,idoc)=                                           &
     &   ((ocetra(i,j,k,idoc)-bacfra+excdoc+exud)*bdp                  &
     &    +ocetra(i,j,k,idoc)       *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)  

         ocetra(i,j,k,icalc)=                                          &
     &   ((ocetra(i,j,k,icalc)+delcar)*bdp                             &
     &    +ocetra(i,j,k,icalc)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k) 

#ifdef __c_isotopes
! Enrichment of CaCO3 resrevoir in 13C&14C isotopes due to CaCO3 formation: fractionation
         ocetra(i,j,k,icalc13)=                                        &
     &   ((ocetra(i,j,k,icalc13)+(delcar*rocean13*(ABS(frac_caco313-1.)+1.)))*bdp   &
     &    +ocetra(i,j,k,icalc13)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k) 

         ocetra(i,j,k,icalc14)=                                        &
     &   ((ocetra(i,j,k,icalc14)+(delcar*rocean14*(ABS(frac_caco314-1.)+1.)))*bdp   &
     &    +ocetra(i,j,k,icalc14)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)

         phyto_growth(i,j,k) = ((ocetra(i,j,k,iphy)+phosy)/ocetra(i,j,k,iphy))/dtb ! Growth rate phytoplankton [1/d]
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
         aux2d_dmsprod(i,j)   = aux2d_dmsprod(i,j)+dmsprod*bdp 
         aux2d_dms_bac(i,j)   = aux2d_dms_bac(i,j)+dms_bac*bdp 
         aux2d_dms_uv(i,j)    = aux2d_dms_uv (i,j)+dms_uv*bdp 
         aux2d_export(i,j)    = aux2d_export(i,j) +export*rcar*bdp 
         aux2d_expoca(i,j)    = aux2d_expoca(i,j) +delcar*bdp
         aux2d_exposi(i,j)    = aux2d_exposi(i,j) +delsil*bdp 
         aux2d_phosy(i,j)     = aux2d_phosy(i,j)  +phosy*rcar*bdp ! primary production in kmol C m-2
         aux3d_phosy(i,j,k)   = phosy*rcar*(bdp/pddpo(i,j,k))     ! primary production in kmol C m-3

     
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

! Accumulate 2d diagnostics 
       call accsrf(jdmsprod,aux2d_dmsprod,omask,0)    
       call accsrf(jdms_uv,aux2d_dms_uv,omask,0)     
       call accsrf(jexport,aux2d_export,omask,0)      
       call accsrf(jexpoca,aux2d_expoca,omask,0)     
       call accsrf(jexposi,aux2d_exposi,omask,0)     
       call accsrf(jintphosy,aux2d_phosy,omask,0)     

! Accumulate primary production 
      call acclyr(jphosy,aux3d_phosy,pddpo,1)
      IF (SUM(jlvlphosy).NE.0) THEN
        DO k=1,kpke
          call bgczlv(pddpo,k,ind1,ind2,wghts)
          call acclvl(jlvlphosy,aux3d_phosy,k,ind1,ind2,wghts)
        ENDDO
      ENDIF

#endif /*bioprod*/



!$OMP PARALLEL DO                                                & 
!$OMP&PRIVATE(sterph,sterzo,remin,docrem,opalrem,aou, &
!$OMP&        refra,dms_bac,detref,rem13,rem14,phythresh,  &
!$OMP&        pocrem,phyrem,bdp,dp_ez)  
      DO 201 j=1,kpje
      DO 201 i=1,kpie
         DO 20 k=kwrbioz(i,j),kpke
            IF(pddpo(i,j,k).gt.dp_min.and.omask(i,j).GT.0.5) THEN

      if(ptiestw(i,j,k).lt.dp_ez.and.ptiestw(i,j,k+1).gt.dp_ez) then
        bdp=abs(ptiestw(i,j,k+1)-dp_ez)
        if(bdp.gt.pddpo(i,j,k)) then
          write(io_stdo_bgc,*) 'bdp gt pddpo(i,j,k) 2',i,j,k
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
               docrem=MIN(dremdoc*ocetra(i,j,k,idoc),                  &
     &                   0.33*ocetra(i,j,k,ioxygen)/ro2ut)
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
            dms_bac = dmsp3*dtb*abs(ptho(i,j,k)+3.)*ocetra(i,j,k,idms)   &
     &             *(ocetra(i,j,k,idms)/(dmsp6+ocetra(i,j,k,idms)))     
            ocetra(i,j,k,ian2o)=                                       &
     &    ((ocetra(i,j,k,ian2o)+remin*1.e-4*ro2ut*refra)*bdp           &
     &     +ocetra(i,j,k,ian2o)        *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
            ocetra(i,j,k,igasnit)=                                     &
     &    ((ocetra(i,j,k,igasnit)-remin*1.e-4*ro2ut*refra)*bdp         &
     &     +ocetra(i,j,k,igasnit)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
            ocetra(i,j,k,ioxygen)=                                     &
     &    ((ocetra(i,j,k,ioxygen)-remin*1.e-4*ro2ut*refra*0.5)*bdp     &
     &     +ocetra(i,j,k,ioxygen)      *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)	 
            ocetra(i,j,k,idms)=                                        &
     &    ((ocetra(i,j,k,idms)-dms_bac)*bdp                            &
     &     +ocetra(i,j,k,idms)         *(pddpo(i,j,k)-bdp))/pddpo(i,j,k)
            aux2d_dms_bac(i,j)   = aux2d_dms_bac(i,j)+dms_bac*bdp
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
      call accsrf(jdms_bac,aux2d_dms_bac,omask,0) 

#ifdef PBGC_OCNP_TIMESTEP 
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'in OCRPOD after poc remin'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif 
!$OMP PARALLEL DO                                                   &
!$OMP&PRIVATE(remin,remin2o,rem13,rem14,detref,detrl,rl13,rl14,bdp,dp_ez) 
       DO 30 j=1,kpje
       DO 30 i=1,kpie
         DO 30 k=kwrbioz(i,j),kpke
         IF(omask(i,j).GT.0.5) THEN
         IF(ocetra(i,j,k,ioxygen).LT.5.e-7.and.pddpo(i,j,k).gt.dp_min) THEN

           if(ptiestw(i,j,k).lt.dp_ez.and.ptiestw(i,j,k+1).gt.dp_ez) then
             bdp=abs(ptiestw(i,j,k+1)-dp_ez)
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

!$OMP PARALLEL DO                                                   &
!$OMP&PRIVATE(remin,rem13,rem14,detref,bdp,dp_ez) 
      DO 301 j=1,kpje
      DO 301 i=1,kpie
        DO 301 k=kwrbioz(i,j),kpke
            IF(omask(i,j).gt.0.5.and.pddpo(i,j,k).gt.dp_min) then  
            IF(ocetra(i,j,k,ioxygen).lt.3.e-6.and.ocetra(i,j,k,iano3).lt.3.e-6) THEN

              if(ptiestw(i,j,k).lt.dp_ez.and.ptiestw(i,j,k+1).gt.dp_ez) then
                bdp=abs(ptiestw(i,j,k+1)-dp_ez)
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


      ! Accumulate diagnostics for agtgregate sinking
      call acclyr(jwphy, wmass/dtb,  pddpo,1)
      call acclyr(jwnos, wnumb/dtb,  pddpo,1)
      call acclyr(jeps,  aux3d_eps,  pddpo,1)
      call acclyr(jasize,aux3d_asize,pddpo,1)
      IF (SUM(jlvlwphy+jlvlwnos+jlvleps+jlvlasize).NE.0) THEN
        DO k=1,kpke
          call bgczlv(pddpo,k,ind1,ind2,wghts)
          call acclvl(jlvlwphy, wmass/dtb,  k,ind1,ind2,wghts)
          call acclvl(jlvlwnos, wnumb/dtb,  k,ind1,ind2,wghts)
          call acclvl(jlvleps,  aux3d_eps,  k,ind1,ind2,wghts)
          call acclvl(jlvlasize,aux3d_asize,k,ind1,ind2,wghts)
        ENDDO
      ENDIF


#endif /*AGG*/


!
! implicit method for sinking of particles:
! C(k,T+dt)=C(k,T) + (w*dt/ddpo(k))*(C(k-1,T+1)-C(k,T+1))
! --> 	    
! C(k,T+dt)=(ddpo(k)*C(k,T)+w*dt*C(k-1,T+dt))/(ddpo(k)+w*dt)
! sedimentation=w*dt*C(ks,T+dt)
!
      k=1
!$OMP PARALLEL DO PRIVATE(wpoc,wpocd,wcal,wcald,wopal,wopald,wnos,wnosd,dagg)
      DO j=1,kpje
      DO i=1,kpie
        IF(omask(i,j).gt.0.5) THEN

#if defined(AGG)
          wpoc  = wmass(i,j,k)
          wcal  = wmass(i,j,k)
          wopal = wmass(i,j,k)
          wnos  = wnumb(i,j,k)
          wdust = dustsink
          dagg  = dustagg(i,j,k)
#elif defined(WLIN)
          wpoc  = wmin
          dagg  = 0.0
#else
          dagg  = 0.0
#endif

          ocetra(i,j,k,idet)  =(ocetra(i,j,k,idet)*pddpo(i,j,k))/        &
     &                         (pddpo(i,j,k)+wpoc)
          ocetra(i,j,k,icalc) =(ocetra(i,j,k,icalc)*pddpo(i,j,k))/       &
     &                         (pddpo(i,j,k)+wcal)
          ocetra(i,j,k,iopal) =(ocetra(i,j,k,iopal)*pddpo(i,j,k))/       &
     &                         (pddpo(i,j,k)+wopal)    
          ocetra(i,j,k,ifdust)=(ocetra(i,j,k,ifdust)*pddpo(i,j,k))/      &
     &                         (pddpo(i,j,k)+wdust) - dagg 
#ifdef AGG
          ocetra(i,j,k,iphy)  =(ocetra(i,j,k,iphy)*pddpo(i,j,k))/        &
     &                         (pddpo(i,j,k)+wpoc)
          ocetra(i,j,k,inos)  =(ocetra(i,j,k,inos)*pddpo(i,j,k))/        &
     &                         (pddpo(i,j,k)+wnos)  - aggregate(i,j,k)
          ocetra(i,j,k,iadust)=(ocetra(i,j,k,iadust)*pddpo(i,j,k))/      &
     &                         (pddpo(i,j,k)+wpoc)  + dagg   
#endif
#ifdef __c_isotopes
          ocetra(i,j,k,idet13) =(ocetra(i,j,k,idet13)*pddpo(i,j,k))/     &
     &                          (pddpo(i,j,k)+wpoc)
          ocetra(i,j,k,idet14) =(ocetra(i,j,k,idet14)*pddpo(i,j,k))/     &
     &                          (pddpo(i,j,k)+wpoc)
          ocetra(i,j,k,icalc13)=(ocetra(i,j,k,icalc13)*pddpo(i,j,k))/    &
     &                          (pddpo(i,j,k)+wcal)
          ocetra(i,j,k,icalc14)=(ocetra(i,j,k,icalc14)*pddpo(i,j,k))/    &
     &                          (pddpo(i,j,k)+wcal)
#endif     
        ENDIF
      enddo
      enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(kdonor,wpoc,wpocd,wcal,wcald,wopal,wopald,wnos,wnosd,dagg)
      DO j=1,kpje
      DO i=1,kpie
        IF(omask(i,j).gt.0.5) THEN
          kdonor=1
          DO k=2,kpke

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

            ocetra(i,j,k,idet)  =(ocetra(i,j,k     ,idet)*pddpo(i,j,k)    &
     &	                         +ocetra(i,j,kdonor,idet)*wpocd)/         &
     &                           (pddpo(i,j,k)+wpoc)
            ocetra(i,j,k,icalc) =(ocetra(i,j,k     ,icalc)*pddpo(i,j,k)   &
     &	                         +ocetra(i,j,kdonor,icalc)*wcald)/        &
     &                           (pddpo(i,j,k)+wcal)
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
          ENDDO  ! k=2,kpke

! Fluxes to sediment, layers thinner than dp_min_sink are ignored.
! Note that kdonor=kbo(i,j) by definition since kbo is the lowermost
! layer thicker than dp_min_sink.
#if defined(AGG)
          wpoc  = wmass(i,j,kdonor)
          wcal  = wmass(i,j,kdonor)
          wopal = wmass(i,j,kdonor)
          prorca(i,j) = ocetra(i,j,kdonor,iphy)  *wpoc                  &
     &                + ocetra(i,j,kdonor,idet)  *wpoc
          prcaca(i,j) = ocetra(i,j,kdonor,icalc) *wcal
          silpro(i,j) = ocetra(i,j,kdonor,iopal) *wopal
          produs(i,j) = ocetra(i,j,kdonor,ifdust)*dustsink              &
     &                + ocetra(i,j,kdonor,iadust)*wpoc    
#elif defined(WLIN)
          wpoc  = min(wmin+wlin*ptiestu(i,j,kdonor),wmax)
          prorca(i,j)=ocetra(i,j,kdonor,idet ) *wpoc
          prcaca(i,j)=ocetra(i,j,kdonor,icalc) *wcal
          silpro(i,j)=ocetra(i,j,kdonor,iopal) *wopal
          produs(i,j)=ocetra(i,j,kdonor,ifdust)*wdust
#else
          prorca(i,j)=ocetra(i,j,kdonor,idet ) *wpoc
          prcaca(i,j)=ocetra(i,j,kdonor,icalc) *wcal
          silpro(i,j)=ocetra(i,j,kdonor,iopal) *wopal
          produs(i,j)=ocetra(i,j,kdonor,ifdust)*wdust
#endif
#ifdef __c_isotopes
          pror13(i,j)=ocetra(i,j,kdonor,idet13 )*wpoc
          prca13(i,j)=ocetra(i,j,kdonor,icalc13)*wcal
          pror14(i,j)=ocetra(i,j,kdonor,idet14 )*wpoc
          prca14(i,j)=ocetra(i,j,kdonor,icalc14)*wcal
#endif

        ENDIF  ! omask.gt.0.5
      ENDDO
      ENDDO


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

      ! Accumulate the diagnostic mass sinking field 
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
