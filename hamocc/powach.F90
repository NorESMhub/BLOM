! Copyright (C) 2001  Ernst Maier-Reimer, S. Legutke
! Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger
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


      SUBROUTINE POWACH(kpie,kpje,kpke,kbnd,prho,omask,psao)
!******************************************************************************
!
!**** *POWACH* - .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!
!     Purpose
!     -------
!     .
!
!     Method
!     -------
!     .
!
!**   Interface.
!     ----------
!
!     *CALL*       *POWACH*
!
!     *COMMON*     *PARAM1_BGC.h* - declaration of ocean/sediment tracer.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *INTEGER* *kbnd*    - nb of halo grid points
!     *REAL*    *prho*    - seawater density [g/cm^3].
!     *REAL*    *psao*    - salinity [psu].
!     *REAL*    *omask*   - land/ocean mask
!
!     Externals
!     ---------
!     none.
!
!******************************************************************************
      USE mo_carbch
      USE mo_chemcon, only: calcon
      USE mo_sedmnt
      USE mo_biomod
      USE mo_control_bgc
      use mo_param1_bgc 
      use mo_vgrid, only: kbo,bolay

      implicit none

      INTEGER, intent(in) :: kpie,kpje,kpke,kbnd
      REAL,    intent(in) :: prho(kpie,kpje,kpke)
      REAL,    intent(in) :: omask(kpie,kpje)
      REAL,    intent(in) :: psao(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd,kpke)

      ! Local variables
      INTEGER :: i,j,k,l
      REAL :: sedb1(kpie,0:ks),sediso(kpie,0:ks)
      REAL :: solrat(kpie,ks),powcar(kpie,ks)
      REAL :: aerob(kpie,ks),anaerob(kpie,ks)
#ifdef cisonew
      REAL :: aerob13(kpie,ks),anaerob13(kpie,ks)
      REAL :: aerob14(kpie,ks),anaerob14(kpie,ks)
#endif
      REAL :: disso, dissot, undsa, silsat, posol 
      REAL :: umfa,denit,saln,rrho,alk,c,sit,pt
      REAL :: K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p
      REAL :: ah1,ac,cu,cb,cc,satlev
      REAL :: ratc13,ratc14,rato13,rato14,poso13,poso14
      ! number of iterations for carchm_solve
      INTEGER,PARAMETER :: niter=5
! *****************************************************************
! accelerated sediment
! needed for boundary layer vertilation in fast sediment routine      

      REAL :: bolven(kpie)

! A LOOP OVER J
! RJ: This loop must go from 1 to kpje in the parallel version,
!     otherways we had to do a boundary exchange


!$OMP PARALLEL DO                                                &                
!$OMP&PRIVATE(sedb1,sediso,solrat,powcar,aerob,anaerob,          &
!$OMP&        disso,dissot,undsa,silsat,posol,                   &
!$OMP&        umfa,denit,saln,rrho,alk,c,sit,pt,                 &
!$OMP&        K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,                &
!$OMP&        ah1,ac,cu,cb,cc,satlev,bolven,                     &     
!$OMP&        ratc13,ratc14,rato13,rato14,poso13,poso14)
      DO 8888 j=1,kpje

      DO 1189 k=1,ks
      DO 1189 i=1,kpie
         solrat(i,k) =0.
         powcar(i,k) =0.
         anaerob(i,k)=0.
         aerob(i,k)  =0.
#ifdef cisonew
         anaerob13(i,k)=0.
         aerob13(i,k)  =0.	 
         anaerob14(i,k)=0.
         aerob14(i,k)  =0.	 
#endif

1189  CONTINUE

! calculate bottom ventilation rate for scaling of sediment-water exchange
      do 1170 i=1,kpie      
      bolven(i) = 1.
1170  continue

      DO 1171 k=0,ks
      DO 1171 i=1,kpie
        sedb1(i,k)=0.
        sediso(i,k)=0.
1171  CONTINUE


! CALCULATE SILICATE-OPAL CYCLE AND SIMULTANEOUS SILICATE DIFFUSION      
!******************************************************************

! Dissolution rate constant of opal (disso) [1/(kmol Si(OH)4/m3)*1/sec]

!      disso=1.e-8 
      disso=1.e-6 ! test vom 03.03.04 half live sil ca. 20.000 yr 
      dissot=disso*dtbgc

! Silicate saturation concentration is 1 mol/m3

      silsat=0.001

! Evaluate boundary conditions for sediment-water column exchange.
! Current undersaturation of bottom water: sedb(i,0) and 
! Approximation for new solid sediment, as from sedimentation flux: solrat(i,1)

      DO 3 i=1,kpie
         IF(omask(i,j).GT.0.5) THEN
!ka         IF(bolay(i,j).GT.0.) THEN
            undsa=silsat-powtra(i,j,1,ipowasi)
            sedb1(i,0)=bolay(i,j)*(silsat-ocetra(i,j,kbo(i,j),isilica)) &
     &                 *bolven(i)       
            solrat(i,1)=                                                &
     &      (sedlay(i,j,1,issssil)+silpro(i,j)/(porsol(1)*seddw(1)))    &
     &      *dissot/(1.+dissot*undsa)*porsol(1)/porwat(1)
         ENDIF
3     CONTINUE


! Evaluate sediment undersaturation and degradation.
! Current undersaturation in pore water: sedb(i,k) and
! Approximation for new solid sediment, as from degradation: solrat(i,k) 

      DO 2 k=1,ks
      DO 2 i=1,kpie
         IF(omask(i,j).GT.0.5) THEN
!ka         IF(bolay(i,j).GT.0.) THEN
            undsa=silsat-powtra(i,j,k,ipowasi)
            sedb1(i,k)=seddw(k)*porwat(k)*(silsat-powtra(i,j,k,ipowasi))
            IF(k.GT.1)solrat(i,k)=sedlay(i,j,k,issssil)                 &
     &                 *dissot/(1.+dissot*undsa)*porsol(k)/porwat(k)
         ENDIF
2     CONTINUE
     
! Solve for new undersaturation sediso, from current undersaturation sedb1,
! and first guess of new solid sediment solrat.     

      CALL powadi(j,kpie,kpje,solrat,sedb1,sediso,bolven,omask)

! Update water column silicate, and store the flux for budget.
! Add sedimentation to first layer.

      DO 4 i=1,kpie
         IF(omask(i,j).GT.0.5) THEN
!ka         IF(bolay(i,j).GT.0.) THEN
         sedfluxo(i,j,ipowasi)=sedfluxo(i,j,ipowasi) +                &
     &  (silsat-sediso(i,0)-ocetra(i,j,kbo(i,j),isilica))*bolay(i,j)

            ocetra(i,j,kbo(i,j),isilica)=silsat-sediso(i,0)
            sedlay(i,j,1,issssil)=                                    &
     &        sedlay(i,j,1,issssil)+silpro(i,j)/(porsol(1)*seddw(1))
         ENDIF
4     CONTINUE


! Calculate updated degradation rate from updated undersaturation.
! Calculate new solid sediment.
! Update pore water concentration from new undersaturation.

      DO 5 k=1,ks
      DO 5 i=1,kpie
         IF(omask(i,j).GT.0.5) THEN
!ka         IF(bolay(i,j).GT.0.) THEN
            umfa=porsol(k)/porwat(k)
            solrat(i,k)=sedlay(i,j,k,issssil)                        &
     &                  *dissot/(1.+dissot*sediso(i,k))
            posol=sediso(i,k)*solrat(i,k)
            sedlay(i,j,k,issssil)=                                   &
     &          sedlay(i,j,k,issssil)-posol
            powtra(i,j,k,ipowasi)=silsat-sediso(i,k)
         ENDIF
5     CONTINUE

     
! CALCULATE OXYGEN-POC CYCLE AND SIMULTANEOUS OXYGEN DIFFUSION      
!*************************************************************

! Degradation rate constant of POP (disso) [1/(kmol O2/m3)*1/sec]

      disso=0.01/86400.  !  disso=3.e-5 was quite high
      dissot=disso*dtbgc

! This scheme is not based on undersaturation, but on O2 itself

! Evaluate boundary conditions for sediment-water column exchange.
! Current concentration of bottom water: sedb(i,0) and 
! Approximation for new solid sediment, as from sedimentation flux: solrat(i,1)

      DO 13 i=1,kpie
         IF(omask(i,j).GT.0.5) THEN
!ka         IF(bolay(i,j).GT.0.) THEN
            undsa=powtra(i,j,1,ipowaox)
            sedb1(i,0)=bolay(i,j)*ocetra(i,j,kbo(i,j),ioxygen)         &
     &                 *bolven(i)       
            solrat(i,1)=                                               &
     &       (sedlay(i,j,1,issso12)+prorca(i,j)/(porsol(1)*seddw(1)))  &
     &          *ro2ut*dissot/(1.+dissot*undsa)*porsol(1)/porwat(1)
         ENDIF
13    CONTINUE

! Evaluate sediment concentration and degradation.
! Current concentration in pore water: sedb(i,k) and
! Approximation for new solid sediment, as from degradation: solrat(i,k) 

      DO 12 k=1,ks
      DO 12 i=1,kpie
!ka         IF(bolay(i,j).GT.0.) THEN
         IF(omask(i,j).GT.0.5) THEN
            undsa=powtra(i,j,k,ipowaox)
            sedb1(i,k)=seddw(k)*porwat(k)*powtra(i,j,k,ipowaox)
            IF(k.GT.1)solrat(i,k)=sedlay(i,j,k,issso12)               &
     &         *ro2ut*dissot/(1.+dissot*undsa)*porsol(k)/porwat(k)
         ENDIF
12     CONTINUE

! Solve for new O2 concentration sediso, from current concentration sedb1,
! and first guess of new solid sediment solrat.     

      CALL powadi(j,kpie,kpje,solrat,sedb1,sediso,bolven,omask)

! Update water column oxygen, and store the flux for budget (opwflux).
! Add sedimentation to first layer.

      DO 14 i=1,kpie
!ka         IF(bolay(i,j).GT.0.) THEN
         IF(omask(i,j).GT.0.5) THEN
            ocetra(i,j,kbo(i,j),ioxygen)=sediso(i,0)
            sedlay(i,j,1,issso12)                                     &
     &      =sedlay(i,j,1,issso12)+prorca(i,j)/(porsol(1)*seddw(1))
         prorca(i,j)=0.
#ifdef cisonew
            sedlay(i,j,1,issso13)                                     &
     &      =sedlay(i,j,1,issso13)+pror13(i,j)/(porsol(1)*seddw(1))
            sedlay(i,j,1,issso14)                                     &
     &      =sedlay(i,j,1,issso14)+pror14(i,j)/(porsol(1)*seddw(1))
         pror13(i,j)=0.
         pror14(i,j)=0.
#endif
         ENDIF
14    CONTINUE


! Calculate updated degradation rate from updated concentration.
! Calculate new solid sediment.
! Update pore water concentration.
! Store flux in array aerob, for later computation of DIC and alkalinity.
      DO 15 k=1,ks
      DO 15 i=1,kpie
!         IF(bolay(i,j).GT.0.) THEN
         IF(omask(i,j).GT.0.5) THEN
            umfa=porsol(k)/porwat(k)
            solrat(i,k)=sedlay(i,j,k,issso12)                         &
     &                 *dissot/(1.+dissot*sediso(i,k))
            posol=sediso(i,k)*solrat(i,k)
            aerob(i,k)=posol*umfa !this has P units: kmol P/m3 of pore water
#ifdef cisonew
            rato13=sedlay(i,j,k,issso13)/(sedlay(i,j,k,issso12)+safediv)
            rato14=sedlay(i,j,k,issso14)/(sedlay(i,j,k,issso12)+safediv)
            poso13=posol*rato13
            poso14=posol*rato14
            aerob13(i,k)=poso13*umfa !this has P units: kmol P/m3 of pore water
            aerob14(i,k)=poso14*umfa !this has P units: kmol P/m3 of pore water
#endif
            sedlay(i,j,k,issso12)=sedlay(i,j,k,issso12)-posol
            powtra(i,j,k,ipowaph)=powtra(i,j,k,ipowaph)+posol*umfa
            powtra(i,j,k,ipowno3)=powtra(i,j,k,ipowno3)+posol*rnit*umfa
            powtra(i,j,k,ipowaox)=sediso(i,k)
#ifdef cisonew
            sedlay(i,j,k,issso13)=sedlay(i,j,k,issso13)-poso13
            sedlay(i,j,k,issso14)=sedlay(i,j,k,issso14)-poso14
            ! is this correct? no correspondance in the lines above
            powtra(i,j,k,ipowc13)=powtra(i,j,k,ipowc13)+poso13*umfa
            powtra(i,j,k,ipowc14)=powtra(i,j,k,ipowc14)+poso14*umfa
#endif
         ENDIF
15    CONTINUE     

! CALCULATE NITRATE REDUCTION UNDRE ANAEROBIC CONDITIONS EXPLICITELY
!*******************************************************************

! Denitrification rate constant of POP (disso) [1/sec]
! Store flux in array anaerob, for later computation of DIC and alkalinity.
      
!ik      denit = 1.e-6*dtbgc
      denit = 0.01/86400. *dtbgc
      DO 124 k=1,ks
      DO 124 i=1,kpie
!ka         if(bolay(i,j).gt.0.) then
         IF(omask(i,j).GT.0.5) THEN
         IF(powtra(i,j,k,ipowaox).LT.1.e-6) THEN
            posol=denit*MIN(0.5*powtra(i,j,k,ipowno3)/114.,          &
     &                          sedlay(i,j,k,issso12))
            umfa=porsol(k)/porwat(k)
            anaerob(i,k)=posol*umfa !this has P units: kmol P/m3 of pore water
#ifdef cisonew
            rato13=sedlay(i,j,k,issso13)/(sedlay(i,j,k,issso12)+safediv)
            rato14=sedlay(i,j,k,issso14)/(sedlay(i,j,k,issso12)+safediv)
            poso13=posol*rato13
            poso14=posol*rato14
            anaerob13(i,k)=poso13*umfa !this has P units: kmol P/m3 of pore water
            anaerob14(i,k)=poso14*umfa !this has P units: kmol P/m3 of pore water
#endif
            sedlay(i,j,k,issso12)=sedlay(i,j,k,issso12)-posol
            powtra(i,j,k,ipowaph)=powtra(i,j,k,ipowaph)+posol*umfa
            powtra(i,j,k,ipowno3)=powtra(i,j,k,ipowno3)-98.*posol*umfa
            powtra(i,j,k,ipown2)=powtra(i,j,k,ipown2)+57.*posol*umfa
#ifdef cisonew
            sedlay(i,j,k,issso13)=sedlay(i,j,k,issso13)-poso13
            sedlay(i,j,k,issso14)=sedlay(i,j,k,issso14)-poso14
            ! is this correct? no corresponance in the lines above
            powtra(i,j,k,ipowc13)=powtra(i,j,k,ipowc13)+poso13*umfa
            powtra(i,j,k,ipowc14)=powtra(i,j,k,ipowc14)+poso14*umfa
#endif
         endif
         ENDIF
124   CONTINUE


!    sulphate reduction in sediments   
     DO 125 k=1,ks
     DO 125 i=1,kpie
        IF(omask(i,j).GT.0.5) THEN
        IF(powtra(i,j,k,ipowaox).lt.3.e-6.and.powtra(i,j,k,ipowno3).lt.3.e-6) THEN
           posol=denit* sedlay(i,j,k,issso12)                                    ! remineralization of poc
           umfa=porsol(k)/porwat(k)
           anaerob(i,k)=anaerob(i,k)+posol*umfa !this has P units: kmol P/m3 of pore water  
                                                !this overwrites anaerob from denitrification. added =anaerob+..., works
#ifdef cisonew
           rato13=sedlay(i,j,k,issso13)/(sedlay(i,j,k,issso12)+safediv)
           rato14=sedlay(i,j,k,issso14)/(sedlay(i,j,k,issso12)+safediv)
           poso13=posol*rato13
           poso14=posol*rato14
           anaerob13(i,k)=anaerob13(i,k)+poso13*umfa !this has P units: kmol P/m3 of pore water  
           anaerob14(i,k)=anaerob13(i,k)+poso14*umfa !this has P units: kmol P/m3 of pore water
#endif
           sedlay(i,j,k,issso12)=sedlay(i,j,k,issso12)-posol
           powtra(i,j,k,ipowaph)=powtra(i,j,k,ipowaph)+posol*umfa
           powtra(i,j,k,ipowno3)=powtra(i,j,k,ipowno3)+posol*umfa*rno3
#ifdef cisonew
           sedlay(i,j,k,issso13)=sedlay(i,j,k,issso13)-poso13
           sedlay(i,j,k,issso14)=sedlay(i,j,k,issso14)-poso14
#endif
        endif
         ENDIF
 125   CONTINUE     ! end sulphate reduction


! CALCULATE CaCO3-CO3 CYCLE AND SIMULTANEOUS CO3-UNDERSATURATION DIFFUSION      
!*************************************************************************


! COMPUTE NEW POWCAR=CARBONATE ION CONCENTRATION IN THE SEDIMENT
! FROM CHANGED ALKALINITY (NITRATE PRODUCTION DURING REMINERALISATION)
! AND DIC GAIN. ITERATE 5 TIMES. THIS CHANGES PH (SEDHPL) OF SEDIMENT.

      DO 1 K=1,KS
      DO 1 i=1,kpie
!ka         IF(bolay(i,j).GT.0.) THEN
         IF(omask(i,j).GT.0.5) THEN
            saln= min(40.,max( 0.,psao(i,j,kbo(i,j))))
            rrho= prho(i,j,kbo(i,j))
            alk = (powtra(i,j,k,ipowaal)-(anaerob(i,k)+aerob(i,k))*16.)  / rrho
            c   = (powtra(i,j,k,ipowaic)+(anaerob(i,k)+aerob(i,k))*122.) / rrho
            sit =  powtra(i,j,k,ipowasi) / rrho
            pt  =  powtra(i,j,k,ipowaph) / rrho
            ah1 = sedhpl(i,j,k)
            K1  = keqb( 1,i,j)
            K2  = keqb( 2,i,j)
            Kb  = keqb( 3,i,j)
            Kw  = keqb( 4,i,j)
            Ks1 = keqb( 5,i,j)
            Kf  = keqb( 6,i,j)
            Ksi = keqb( 7,i,j)
            K1p = keqb( 8,i,j)
            K2p = keqb( 9,i,j)
            K3p = keqb(10,i,j)

            CALL CARCHM_SOLVE(saln,c,alk,sit,pt,                  &
                              K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p, &
                              ah1,ac,niter)

            cu = ( 2. * c - ac ) / ( 2. + K1 / ah1 )
            cb = K1 * cu / ah1
            cc = K2 * cb / ah1
            sedhpl(i,j,k) = max(1.e-20,ah1)
            powcar(i,k)   = cc * rrho 

         ENDIF
1     CONTINUE


! Dissolution rate constant of CaCO3 (disso) [1/(kmol CO3--/m3)*1/sec]
      disso=1.e-7
      dissot=disso*dtbgc

! Evaluate boundary conditions for sediment-water column exchange.
! Current undersaturation of bottom water: sedb(i,0) and 
! Approximation for new solid sediment, as from sedimentation flux: solrat(i,1)

! CO3 saturation concentration is aksp/calcon as in CARCHM 
! (calcon defined in BELEG_BGC with 1.03e-2; 1/calcon =~ 97.) 

      DO 23 i=1,kpie
!ka         IF(bolay(i,j).GT.0.) THEN
         IF(omask(i,j).GT.0.5) THEN
            satlev=keqb(11,i,j)/calcon+2.e-5
            undsa=MAX(satlev-powcar(i,1),0.)
            sedb1(i,0)=bolay(i,j)*(satlev-co3(i,j,kbo(i,j)))             &
     &                 *bolven(i)       
            solrat(i,1)=                                                 &
     &         (sedlay(i,j,1,isssc12)+prcaca(i,j)/(porsol(1)*seddw(1)))  &
     &         *dissot/(1.+dissot*undsa)*porsol(1)/porwat(1)
         ENDIF
23     CONTINUE

! Evaluate sediment undersaturation and degradation.
! Current undersaturation in pore water: sedb(i,k) and
! Approximation for new solid sediment, as from degradation: solrat(i,k) 

      DO 22 k=1,ks
      DO 22 i=1,kpie
!ka         IF(bolay(i,j).GT.0.) THEN
         IF(omask(i,j).GT.0.5) THEN
            undsa=MAX(keqb(11,i,j)/calcon-powcar(i,k),0.)
            sedb1(i,k)=seddw(k)*porwat(k)*undsa
            IF(k.GT.1)solrat(i,k)=sedlay(i,j,k,isssc12)                 &
     &          *dissot/(1.+dissot*undsa)*porsol(k)/porwat(k)
            IF(undsa.LE.0.) solrat(i,k)=0.
         ENDIF
22     CONTINUE

! Solve for new undersaturation sediso, from current undersaturation sedb1,
! and first guess of new solid sediment solrat.     

      CALL powadi(j,kpie,kpje,solrat,sedb1,sediso,bolven,omask)
   
! There is no exchange between water and sediment with respect to co3 so far.
! Add sedimentation to first layer.
      DO 24 i=1,kpie
         IF(omask(i,j).GT.0.5) THEN
!ka         IF(bolay(i,j).GT.0.) THEN
            sedlay(i,j,1,isssc12)=                                     &
     &      sedlay(i,j,1,isssc12)+prcaca(i,j)/(porsol(1)*seddw(1))
         prcaca(i,j)=0.
#ifdef cisonew
            sedlay(i,j,1,isssc13)=                                     &
     &      sedlay(i,j,1,isssc13)+prca13(i,j)/(porsol(1)*seddw(1))
            sedlay(i,j,1,isssc14)=                                     &
     &      sedlay(i,j,1,isssc14)+prca14(i,j)/(porsol(1)*seddw(1))
         prca13(i,j)=0.
         prca14(i,j)=0.
#endif
         ENDIF
24    CONTINUE

! Calculate updated degradation rate from updated undersaturation.
! Calculate new solid sediment.
! No update of powcar pore water concentration from new undersaturation so far.
! Instead, only update DIC, and, of course, alkalinity.
! This also includes gains from aerobic and anaerobic decomposition.

      DO 25 k=1,ks
      DO 25 i=1,kpie
         IF(omask(i,j).GT.0.5) THEN
!ka         IF(bolay(i,j).GT.0.) THEN
           umfa=porsol(k)/porwat(k)
           solrat(i,k)=sedlay(i,j,k,isssc12)                           &
     &                 *dissot/(1.+dissot*sediso(i,k))
           posol=sediso(i,k)*solrat(i,k)
#ifdef cisonew
           ratc13=sedlay(i,j,k,isssc13)/(sedlay(i,j,k,isssc12)+safediv)
           ratc14=sedlay(i,j,k,isssc14)/(sedlay(i,j,k,isssc12)+safediv)
           poso13=posol*ratc13
           poso14=posol*ratc14
#endif
           sedlay(i,j,k,isssc12)=sedlay(i,j,k,isssc12)-posol
           powtra(i,j,k,ipowaic)=powtra(i,j,k,ipowaic)                 &
     &        +posol*umfa+(aerob(i,k)+anaerob(i,k))*122.
           powtra(i,j,k,ipowaal)=powtra(i,j,k,ipowaal)                 &
     &        +2.*posol*umfa-16.*(aerob(i,k)+anaerob(i,k))
#ifdef cisonew
           sedlay(i,j,k,isssc13)=sedlay(i,j,k,isssc13)-poso13
           sedlay(i,j,k,isssc14)=sedlay(i,j,k,isssc14)-poso14
           powtra(i,j,k,ipowc13)=powtra(i,j,k,ipowc13)+poso13*umfa   &
     &          +(aerob13(i,k)+anaerob13(i,k))*122.
           powtra(i,j,k,ipowc14)=powtra(i,j,k,ipowc14)+poso14*umfa   &
     &          +(aerob14(i,k)+anaerob14(i,k))*122.
#endif
         ENDIF
25    CONTINUE     


8888  CONTINUE
!$OMP END PARALLEL DO

      CALL DIPOWA(kpie,kpje,kpke,omask)


!ik add clay sedimentation onto sediment
!ik this is currently assumed to depend on total and corg sedimentation:
!ik f(POC) [kg C] / f(total) [kg] = 0.05
!ik thus it is 
!$OMP PARALLEL DO  
       do j=1,kpje
       do i=1,kpie
       sedlay(i,j,1,issster) = sedlay(i,j,1,issster)                 &
     &                       + produs(i,j)/(porsol(1)*seddw(1))     
       enddo
       enddo
!$OMP END PARALLEL DO
       

!$OMP PARALLEL DO  
       DO 91 j=1,kpje
       DO 91 i=1,kpie
         silpro(i,j)=0.
         prorca(i,j)=0.
         prcaca(i,j)=0.
#ifdef cisonew
         pror13(i,j)=0.
         pror14(i,j)=0.
         prca13(i,j)=0.
         prca14(i,j)=0.
#endif
         produs(i,j)=0.
91     CONTINUE
!$OMP END PARALLEL DO

      RETURN
      END
