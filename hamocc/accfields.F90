! Copyright (C) 2020  J. Schwinger, A. Moree
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


      SUBROUTINE ACCFIELDS(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask)       
!**********************************************************************
!
!**** *ACCFIELDS* - .
!
!     J.Schwinger,    *UNI-RESEARCH*    2018-03-22
!
!     Modified
!     --------
!
!     Purpose
!     -------
!     Accumulate fields for time-avaraged output and write output
!
!
!
!**** Parameter list:
!     ---------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st dimension) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (2nd dimension) [m].
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!     *REAL*    *omask*   - land/ocean mask
!
!**********************************************************************
      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      USE mo_bgcmean
      USE mo_control_bgc
      use mo_param1_bgc 
      use mo_vgrid, only: dp_min
      use mod_xc

      implicit none
      INTEGER :: kpie,kpje,kpke
      REAL    :: pdlxp(kpie,kpje)
      REAL    :: pdlyp(kpie,kpje)
      REAL    :: pddpo(kpie,kpje,kpke)
      REAL    :: omask(kpie,kpje)

! Local variables
      INTEGER :: i,j,k,l
      INTEGER :: ind1(kpie,kpje),ind2(kpie,kpje)
      REAL    :: wghts(kpie,kpje,ddm)

#ifdef cisonew
      REAL    :: di12c
      REAL    :: d13c(kpie,kpje,kpke)
      REAL    :: d14c(kpie,kpje,kpke)
      REAL    :: bigd14c(kpie,kpje,kpke)


! Calculation d13C, d14C and Dd14C: Delta notation for output
      d13c(:,:,:)=0.
      d14c(:,:,:)=0.
      bigd14c(:,:,:)=0.
      do k=1,kpke
      do j=1,kpje
      do i=1,kpie
        if(omask(i,j).gt.0.5.and.pddpo(i,j,k).gt.dp_min) then

        di12c=max(ocetra(i,j,k,isco212)-ocetra(i,j,k,isco213),0.)
        d13c(i,j,k)=(ocetra(i,j,k,isco213)/(di12c+safediv)/re1312-1.)*1000.
        d14c(i,j,k)=(ocetra(i,j,k,isco214)*c14fac/(ocetra(i,j,k,isco212)+safediv)/re14to-1.)*1000.
        bigd14c(i,j,k)=d14c(i,j,k)-2.*(d13c(i,j,k)+25.)*(1.+d14c(i,j,k)/1000.)

        endif
      enddo
      enddo
      enddo
#endif


! Accumulated fluxes for inventory.F90. Note that these are currently not written to restart!
! Division by 2 is to account for leap-frog timestepping (but this is not exact)
      do j=1,kpje
      do i=1,kpie
        if(omask(i,j).gt.0.5) then
 
        bgct2d(i,j,jco2flux) = bgct2d(i,j,jco2flux) + atmflx(i,j,iatmco2)/2.0
        bgct2d(i,j,jo2flux)  = bgct2d(i,j,jo2flux)  + atmflx(i,j,iatmo2)/2.0
        bgct2d(i,j,jn2flux)  = bgct2d(i,j,jn2flux)  + atmflx(i,j,iatmn2)/2.0
        bgct2d(i,j,jn2oflux) = bgct2d(i,j,jn2oflux) + atmflx(i,j,iatmn2o)/2.0
    
        endif
      enddo
      enddo


! Accumulate atmosphere fields and fluxes
      call accsrf(jatmco2,atm(1,1,iatmco2),omask,0)
#if defined(BOXATM)
      call accsrf(jatmo2 ,atm(1,1,iatmo2),omask,0)
      call accsrf(jatmn2 ,atm(1,1,iatmn2),omask,0)
#endif
      call accsrf(joxflux,atmflx(1,1,iatmo2),omask,0)
      call accsrf(jniflux,atmflx(1,1,iatmn2),omask,0)
      call accsrf(jn2ofx,atmflx(1,1,iatmn2o),omask,0)
      call accsrf(jdmsflux,atmflx(1,1,iatmdms),omask,0)
#ifdef CFC
      call accsrf(jcfc11fx,atmflx(1,1,iatmf11),omask,0)
      call accsrf(jcfc12fx,atmflx(1,1,iatmf12),omask,0)
      call accsrf(jsf6fx,atmflx(1,1,iatmsf6),omask,0)
#endif
#ifdef natDIC
      call accsrf(jnatco2fx,atmflx(1,1,iatmnco2),omask,0)
#endif
#ifdef cisonew
      call accsrf(jatmc13,atm(1,1,iatmc13),omask,0)
      call accsrf(jatmc14,atm(1,1,iatmc14),omask,0)
#endif

      ! Save up and downward fluxes for CO2 seperately
      call accsrf(jco2fxd,co2fxd,omask,0)
      call accsrf(jco2fxu,co2fxu,omask,0)
#ifdef cisonew
      call accsrf(jco213fxd,co213fxd,omask,0)
      call accsrf(jco213fxu,co213fxu,omask,0)
      call accsrf(jco214fxd,co214fxd,omask,0)
      call accsrf(jco214fxu,co214fxu,omask,0)
#endif

! Accumulate 2d diagnostics
      call accsrf(jpco2,pco2d,omask,0)
      call accsrf(jkwco2,kwco2sol,omask,0)
      call accsrf(jsrfphosph,ocetra(1,1,1,iphosph),omask,0)
      call accsrf(jsrfoxygen,ocetra(1,1,1,ioxygen),omask,0)
      call accsrf(jsrfiron,ocetra(1,1,1,iiron),omask,0)
      call accsrf(jsrfano3,ocetra(1,1,1,iano3),omask,0)
      call accsrf(jsrfalkali,ocetra(1,1,1,ialkali),omask,0)
      call accsrf(jsrfsilica,ocetra(1,1,1,isilica),omask,0)
      call accsrf(jsrfdic,ocetra(1,1,1,isco212),omask,0)
      call accsrf(jsrfphyto,ocetra(1,1,1,iphy),omask,0)
      call accsrf(jdms,ocetra(1,1,1,idms),omask,0)
      call accsrf(jexport,expoor,omask,0)      
      call accsrf(jexpoca,expoca,omask,0)     
      call accsrf(jexposi,exposi,omask,0)     
      call accsrf(jdmsprod,intdmsprod,omask,0)    
      call accsrf(jdms_uv,intdms_uv,omask,0)     
      call accsrf(jdms_bac,intdms_bac,omask,0) 
      call accsrf(jintphosy,intphosy,omask,0)     
      call accsrf(jintdnit,intdnit,omask,0)
      call accsrf(jintnfix,intnfix,omask,0)
#ifdef natDIC
      call accsrf(jsrfnatdic,ocetra(1,1,1,inatsco212),omask,0)
      call accsrf(jsrfnatalk,ocetra(1,1,1,inatalkali),omask,0)
      call accsrf(jnatpco2,natpco2d,omask,0)
#endif

! Accumulate the diagnostic mass sinking field 
      IF( domassfluxes ) THEN
        call accsrf(jcarflx0100,carflx0100,omask,0)    
        call accsrf(jbsiflx0100,bsiflx0100,omask,0)    
        call accsrf(jcalflx0100,calflx0100,omask,0)    
        call accsrf(jcarflx0500,carflx0500,omask,0)    
        call accsrf(jbsiflx0500,bsiflx0500,omask,0)    
        call accsrf(jcalflx0500,calflx0500,omask,0)    
        call accsrf(jcarflx1000,carflx1000,omask,0)    
        call accsrf(jbsiflx1000,bsiflx1000,omask,0)    
        call accsrf(jcalflx1000,calflx1000,omask,0)    
        call accsrf(jcarflx2000,carflx2000,omask,0)    
        call accsrf(jbsiflx2000,bsiflx2000,omask,0)    
        call accsrf(jcalflx2000,calflx2000,omask,0)    
        call accsrf(jcarflx4000,carflx4000,omask,0)    
        call accsrf(jbsiflx4000,bsiflx4000,omask,0)    
        call accsrf(jcalflx4000,calflx4000,omask,0)    
        call accsrf(jcarflx_bot,carflx_bot,omask,0)    
        call accsrf(jbsiflx_bot,bsiflx_bot,omask,0)    
        call accsrf(jcalflx_bot,calflx_bot,omask,0)    
      ENDIF

! Accumulate layer diagnostics
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
      call acclyr(jomegaa,OmegaA,pddpo,1)
      call acclyr(jomegac,OmegaC,pddpo,1)
      call acclyr(jphosy,phosy3d,pddpo,1)
      call acclyr(jo2sat,satoxy,pddpo,1) 
      call acclyr(jprefo2,ocetra(1,1,1,iprefo2),pddpo,1)
      call acclyr(jprefpo4,ocetra(1,1,1,iprefpo4),pddpo,1)
      call acclyr(jprefalk,ocetra(1,1,1,iprefalk),pddpo,1)
      call acclyr(jprefdic,ocetra(1,1,1,iprefdic),pddpo,1)
      call acclyr(jdicsat,ocetra(1,1,1,idicsat),pddpo,1)
#ifdef natDIC
      call acclyr(jnatalkali,ocetra(1,1,1,inatalkali),pddpo,1)
      call acclyr(jnatdic,ocetra(1,1,1,inatsco212),pddpo,1)
      call acclyr(jnatcalc,ocetra(1,1,1,inatcalc),pddpo,1)
      call acclyr(jnatco3,natco3,pddpo,1)                      
      call acclyr(jnatph,nathi,pddpo,1)
      call acclyr(jnatomegaa,natOmegaA,pddpo,1)
      call acclyr(jnatomegac,natOmegaC,pddpo,1)
#endif
#ifdef cisonew
      call acclyr(jdic13,ocetra(1,1,1,isco213),pddpo,1)    
      call acclyr(jdic14,ocetra(1,1,1,isco214),pddpo,1)    
      call acclyr(jd13c,d13c,pddpo,1)    
      call acclyr(jd14c,d14c,pddpo,1)    
      call acclyr(jbigd14c,bigd14c,pddpo,1)    
      call acclyr(jpoc13,ocetra(1,1,1,idet13),pddpo,1)
      call acclyr(jdoc13,ocetra(1,1,1,idoc13),pddpo,1)
      call acclyr(jcalc13,ocetra(1,1,1,icalc13),pddpo,1)
      call acclyr(jphyto13,ocetra(1,1,1,iphy13),pddpo,1)   
      call acclyr(jgrazer13,ocetra(1,1,1,izoo13),pddpo,1)  
#endif 
#ifdef AGG
      call acclyr(jnos,ocetra(1,1,1,inos),pddpo,1)      
      call acclyr(jwphy, wmass/dtb,pddpo,1)
      call acclyr(jwnos, wnumb/dtb,pddpo,1)
      call acclyr(jeps,  eps3d,    pddpo,1)
      call acclyr(jasize,asize3d,  pddpo,1)
#endif     
#ifdef CFC
      call acclyr(jcfc11,ocetra(1,1,1,icfc11),pddpo,1)
      call acclyr(jcfc12,ocetra(1,1,1,icfc12),pddpo,1)
      call acclyr(jsf6,ocetra(1,1,1,isf6),pddpo,1)
#endif


! Accumulate level diagnostics
      IF (SUM(jlvlphyto+jlvlgrazer+jlvlphosph+jlvloxygen+jlvliron+      &
     &  jlvlano3+jlvlalkali+jlvlsilica+jlvldic+jlvldoc+jlvlpoc+jlvlcalc+&
     &  jlvlopal+jlvln2o+jlvlco3+jlvlph+jlvlomegaa+jlvlomegac+jlvlphosy+&
     &  jlvlo2sat+jlvlprefo2+jlvlprefpo4+jlvlprefalk+jlvlprefdic+       &
     &  jlvldicsat+jlvlnatdic+jlvlnatalkali+jlvlnatcalc+jlvlnatco3+     &
     &  jlvlnatomegaa+jlvlnatomegac+jlvldic13+jlvldic14+jlvld13c+       &
     &  jlvld14c+jlvlbigd14c+jlvlpoc13+jlvldoc13+jlvlcalc13+jlvlphyto13+&
     &  jlvlgrazer13+jlvlnos+jlvlwphy+jlvlwnos+jlvleps+jlvlasize+       &
     &  jlvlcfc11+jlvlcfc12+jlvlsf6).NE.0) THEN
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
          call acclvl(jlvlomegaa,OmegaA,k,ind1,ind2,wghts)
          call acclvl(jlvlomegac,OmegaC,k,ind1,ind2,wghts)
          call acclvl(jlvlphosy,phosy3d,k,ind1,ind2,wghts)
          call acclvl(jlvlo2sat,satoxy,k,ind1,ind2,wghts)          
          call acclvl(jlvlprefo2,ocetra(1,1,1,iprefo2),k,ind1,ind2,wghts)
          call acclvl(jlvlprefpo4,ocetra(1,1,1,iprefpo4),k,ind1,ind2,wghts)
          call acclvl(jlvlprefalk,ocetra(1,1,1,iprefalk),k,ind1,ind2,wghts)
          call acclvl(jlvlprefdic,ocetra(1,1,1,iprefdic),k,ind1,ind2,wghts)
          call acclvl(jlvldicsat,ocetra(1,1,1,idicsat),k,ind1,ind2,wghts)
#ifdef natDIC
          call acclvl(jlvlnatdic,ocetra(1,1,1,inatsco212),k,ind1,ind2,wghts)
          call acclvl(jlvlnatalkali,ocetra(1,1,1,inatalkali),k,ind1,ind2,wghts)
          call acclvl(jlvlnatcalc,ocetra(1,1,1,inatcalc),k,ind1,ind2,wghts)
          call acclvl(jlvlnatco3,natco3,k,ind1,ind2,wghts)
          call acclvl(jlvlnatph,nathi,k,ind1,ind2,wghts)
          call acclvl(jlvlnatomegaa,natOmegaA,k,ind1,ind2,wghts)
          call acclvl(jlvlnatomegac,natOmegaC,k,ind1,ind2,wghts)
#endif
#ifdef cisonew
          call acclvl(jlvld13c,d13c,k,ind1,ind2,wghts)
          call acclvl(jlvld14c,d14c,k,ind1,ind2,wghts)
          call acclvl(jlvlbigd14c,bigd14c,k,ind1,ind2,wghts)
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
          call acclvl(jlvlwphy, wmass/dtb,k,ind1,ind2,wghts)
          call acclvl(jlvlwnos, wnumb/dtb,k,ind1,ind2,wghts)
          call acclvl(jlvleps,  eps3d,    k,ind1,ind2,wghts)
          call acclvl(jlvlasize,asize3d,  k,ind1,ind2,wghts)
#endif     
#ifdef CFC
          call acclvl(jlvlcfc11,ocetra(1,1,1,icfc11),k,ind1,ind2,wghts)
          call acclvl(jlvlcfc12,ocetra(1,1,1,icfc12),k,ind1,ind2,wghts)
          call acclvl(jlvlsf6,ocetra(1,1,1,isf6),k,ind1,ind2,wghts)
#endif
        ENDDO
      ENDIF


#ifndef sedbypass
! Accumulate sediments
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

! Accumulate sediment burial
      call accbur(jburssso12,burial(1,1,issso12))
      call accbur(jburssssil,burial(1,1,issssil))
      call accbur(jbursssc12,burial(1,1,isssc12))
      call accbur(jburssster,burial(1,1,issster))
#endif


! Write output if requested
      DO l=1,nbgc 
        nacc_bgc(l)=nacc_bgc(l)+1
        if (bgcwrt(l).gt.0.5) then
          if (GLB_INVENTORY(l).ne.0)                                    & 
     &      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
          call ncwrt_bgc(l)
          nacc_bgc(l)=0 
        endif
      ENDDO

     RETURN
     END
