! Copyright (C) 2002  P. Wetzel
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


      SUBROUTINE INVENTORY_BGC(kpie,kpje,kpke,dlxp,dlyp,ddpo,WETO     &
     &                        ,volchck)
!*******************************************************************
!
!**** *INVENTORY_BGC* - calculate the BGC inventory.
!
!     P.Wetzel,              *MPI-Met, HH*    29.07.02
!
!     Modified
!     --------
!     
!     Purpose
!     -------
!     - calculate the BGC inventory.
!
!     Method
!     -------
!     -
!
!**   Interface.
!     ----------
!
!     *CALL*       *INVENTORY_BGC*
!
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************
    
      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      USE mo_control_bgc
      USE mo_bgcmean
      USE mo_param1_bgc 
      use mo_vgrid, only: dp_min
      USE mod_xc
      
      implicit none

      INTEGER :: kpie,kpje,kpke,i,j,k,l,volchck

      REAL :: ddpo(kpie,kpje,kpke)
      REAL :: dlxp(kpie,kpje),dlyp(kpie,kpje)
      REAL :: WETO(kpie,kpje)
      REAL :: ztmp1(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy)
      REAL :: ztmp2(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy)

      REAL :: zpowtrato(npowtra)
      REAL :: zsedlayto(nsedtra),zburial(nsedtra)
      REAL :: zocetrato(nocetra)
      REAL :: zalkali
      REAL :: ztotvol,ztotarea,vol,zsedhplto
      REAL :: zhito,zco3to,sum,zprorca,zprcaca,zsilpro
      REAL :: zatmco2,zatmo2,zatmn2
      REAL :: co2flux,so2flux,sn2flux,sn2oflux
      REAL :: totalcarbon,totalphos,totalsil,totalnitr,totaloxy
      REAL :: ppm2con, co2atm

! aqueous sediment tracer
!----------------------------------------------------------------------
#ifdef sedbypass
      zpowtrato(:)=0.0
      zsedlayto(:)=0.0
      zburial(:)=0.0
      zsedhplto=0.0
#else
      ztmp1(:,:)=0.0
      DO k=1,ks
      DO j=1,kpje
      DO i=1,kpie
        ztmp1(i,j)=ztmp1(i,j)+WETO(I,J)*seddw(k)                    &
     &             *dlxp(i,j)*dlyp(i,j)*porwat(k)
      ENDDO
      ENDDO
      ENDDO

      CALL xcsum(ztotvol,ztmp1,ips)

      DO l=1,npowtra

         ztmp1(:,:)=0.0
         DO k=1,ks
         DO j=1,kpje
         DO i=1,kpie
            vol    = seddw(k)*dlxp(i,j)*dlyp(i,j)*porwat(k)
            ztmp1(i,j)=ztmp1(i,j)+WETO(I,J)*powtra(i,j,k,l)*vol
         ENDDO
         ENDDO
         ENDDO

         CALL xcsum(zpowtrato(l),ztmp1,ips)

      ENDDO


      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*)'Global inventory of aqueous sediment tracer'
      WRITE(io_stdo_bgc,*)'-------------------------------------------'
      WRITE(io_stdo_bgc,*) '       total[kmol]    concentration[mol/L]'
      DO l=1,npowtra
      WRITE(io_stdo_bgc,*)'No. ',l,' ', &
     &       zpowtrato(l),'  ',zpowtrato(l)/ztotvol,'  ',ztotvol
      ENDDO
      WRITE(io_stdo_bgc,*) ' '
      ENDIF

! non aqueous sediment tracer
!----------------------------------------------------------------------
      DO l=1,nsedtra

         ztmp1(:,:)=0.0
         DO k=1,ks
         DO j=1,kpje
         DO i=1,kpie
            vol = porsol(k)*seddw(k)*dlxp(i,j)*dlyp(i,j)
            ztmp1(i,j)=ztmp1(i,j)+WETO(I,J)*sedlay(i,j,k,l)*vol
         ENDDO
         ENDDO
         ENDDO

         CALL xcsum(zsedlayto(l),ztmp1,ips)

      ENDDO

      DO l=1,nsedtra

         ztmp1(:,:)=0.0
         DO j=1,kpje
         DO i=1,kpie
            ztmp1(i,j)=burial(i,j,l)*WETO(I,J)*dlxp(i,j)*dlyp(i,j)
         ENDDO
         ENDDO

         CALL xcsum(zburial(l),ztmp1,ips)

      ENDDO

      ztmp1(:,:)=0.0
      DO k=1,ks
      DO j=1,kpje
      DO i=1,kpie
        vol = porsol(k)*seddw(k)*dlxp(i,j)*dlyp(i,j)
        ztmp1(i,j)=ztmp1(i,j)+WETO(I,J)*sedhpl(i,j,k)*vol
      ENDDO
      ENDDO
      ENDDO

      CALL xcsum(zsedhplto,ztmp1,ips)




      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*)                                            &
     &     'Global inventory of solid sediment constituents'
      WRITE(io_stdo_bgc,*)                                            &
     &     '----------------------------------------------------'
      WRITE(io_stdo_bgc,*) '        [kmol]'

      DO l=1,nsedtra
        WRITE(io_stdo_bgc,*) 'Sediment No. ',l,' ',                   &
     &                        zsedlayto(l)
        WRITE(io_stdo_bgc,*) 'Burial No. ',l,' ',                     &
     &                        zburial(l)
      ENDDO
      WRITE(io_stdo_bgc,*) 'hpl ',                                    &
     &       zsedhplto   
      WRITE(io_stdo_bgc,*) ' '   
      ENDIF

#endif

!  oceanic tracers
!----------------------------------------------------------------------

      ztmp1(:,:)=0.0
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
         IF(ddpo(i,j,k).gt.dp_min) THEN
            ztmp1(i,j)=ztmp1(i,j)+WETO(I,J)*dlxp(i,j)*dlyp(i,j)     &
     &                 *DDPO(i,j,k)
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      CALL xcsum(ztotvol,ztmp1,ips)

      DO l=1,nocetra

         ztmp1(:,:)=0.0
         DO k=1,kpke
         DO j=1,kpje
         DO i=1,kpie
           IF(ddpo(i,j,k).gt.dp_min) THEN
             vol = dlxp(i,j)*dlyp(i,j)*DDPO(i,j,k)
             ztmp1(i,j)=ztmp1(i,j)+WETO(I,J)*ocetra(i,j,k,l)*vol
!             if (ocetra(i,j,k,l).lt.0.0) then
!      WRITE(io_stdo_bgc,*) 'ocetra -ve', l,ocetra(i,j,k,l)
!             endif
           ENDIF
         ENDDO
         ENDDO
         ENDDO

         CALL xcsum(zocetrato(l),ztmp1,ips)

      ENDDO
    
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global inventory of advected ocean tracers'
      WRITE(io_stdo_bgc,*) '------------------------------------------'
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) '       total[kmol]  concentration[kmol/m^3]'
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'ztotvol',ztotvol
      DO l=1,nocetra
      WRITE(io_stdo_bgc,*) 'No. ',l,    &
     &         zocetrato(l),zocetrato(l)/ztotvol
      ENDDO
      ENDIF

! additional ocean tracer
!----------------------------------------------------------------------

      ztmp1(:,:)=0.0
      ztmp2(:,:)=0.0
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
         IF(ddpo(i,j,k).gt.dp_min) THEN
         vol = dlxp(i,j)*dlyp(i,j)*DDPO(i,j,k)
         ztmp1(i,j) = ztmp1(i,j) + WETO(I,J)*hi(i,j,k) *vol
         ztmp2(i,j) = ztmp2(i,j) + WETO(I,J)*co3(i,j,k)*vol
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      CALL xcsum(zhito ,ztmp1,ips)
      CALL xcsum(zco3to,ztmp2,ips)
 
!      IF (mnproc.eq.1) THEN
!      WRITE(io_stdo_bgc,*) ' '
!      WRITE(io_stdo_bgc,*) 'Glob. inventory of additional ocean tracer'
!      WRITE(io_stdo_bgc,*) '------------------------------------------'
!      WRITE(io_stdo_bgc,*) '      total[kmol]  concentration[kmol/m^3]'
!      WRITE(io_stdo_bgc,*) ' '
!
!      WRITE(io_stdo_bgc,*) ' hi',            &
!     &       zhito,zhito/ztotvol
!      WRITE(io_stdo_bgc,*) ' co3',           &
!     &       zco3to,zco3to/ztotvol
!      WRITE(io_stdo_bgc,*) ' '
!      ENDIF


! alkalinity of the first layer
!-------------------------------------------------------------------- 

      k=1
      ztmp1(:,:)=0.0
      ztmp2(:,:)=0.0
      DO j=1,kpje
      DO i=1,kpie
         ztmp1(i,j) = WETO(I,J)*dlxp(i,j)*dlyp(i,j)*DDPO(i,j,k)    
         ztmp2(i,j) = ocetra(i,j,k,ialkali)*ztmp1(i,j)
      ENDDO
      ENDDO

      CALL xcsum(ztotvol,ztmp1,ips)
      CALL xcsum(zalkali,ztmp2,ips)
    
!      IF (mnproc.eq.1) THEN
!      WRITE(io_stdo_bgc,*) ' '
!      WRITE(io_stdo_bgc,*) 'Global inventory of first layer alkalinity'
!      WRITE(io_stdo_bgc,*) '------------------------------------------'
!      WRITE(io_stdo_bgc,*) ' '
!      WRITE(io_stdo_bgc,*) '       total[kmol]  concentration[kmol/m^3]'
!      WRITE(io_stdo_bgc,*) ' '
!      WRITE(io_stdo_bgc,*) zalkali,zalkali/ztotvol
!      ENDIF
     


! atmosphere flux and atmospheric CO2
!-------------------------------------------------------------------- 

!      IF (mnproc.eq.1) THEN
!      WRITE(io_stdo_bgc,*) ' '
!      WRITE(io_stdo_bgc,*) 'Global fluxes into atmosphere'
!      WRITE(io_stdo_bgc,*) '-----------------------------'
!      WRITE(io_stdo_bgc,*) '        [kmol]'
!      ENDIF

      co2flux  =0.
      so2flux  =0.
      sn2flux  =0.
      sn2oflux =0.
      zatmco2  =0.
      zatmo2   =0.
      zatmn2   =0.
      ztotarea =0.
      ppm2con=0.35e-3
      do i=1,kpie
      do j=1,kpje
        co2flux =co2flux +bgct2d(i,j,jco2flux)*dlxp(i,j)*dlyp(i,j)
        so2flux =so2flux +bgct2d(i,j,jo2flux) *dlxp(i,j)*dlyp(i,j)
        sn2flux =sn2flux +bgct2d(i,j,jn2flux) *dlxp(i,j)*dlyp(i,j)
        sn2oflux=sn2oflux+bgct2d(i,j,jn2oflux)*dlxp(i,j)*dlyp(i,j)
        ztotarea = ztotarea + dlxp(i,j)*dlyp(i,j)
        zatmco2 =zatmco2 + atm(i,j,iatmco2)*dlxp(i,j)*dlyp(i,j)
#if defined(BOXATM)
        zatmo2= zatmo2  + atm(i,j,iatmo2) *dlxp(i,j)*dlyp(i,j)
        zatmn2= zatmn2  + atm(i,j,iatmn2) *dlxp(i,j)*dlyp(i,j)	
#endif
      enddo
      enddo

      ztmp1(:,:)=0.0
      DO j=1,kpje
      DO i=1,kpie
        ztmp1(i,j) = bgct2d(i,j,jco2flux)*dlxp(i,j)*dlyp(i,j)
      ENDDO
      ENDDO

      CALL xcsum(co2flux,ztmp1,ips)

      ztmp1(:,:)=0.0
      DO j=1,kpje
      DO i=1,kpie
        ztmp1(i,j) = bgct2d(i,j,jo2flux)*dlxp(i,j)*dlyp(i,j)
      ENDDO
      ENDDO

      CALL xcsum(so2flux,ztmp1,ips)

      ztmp1(:,:)=0.0
      DO j=1,kpje
      DO i=1,kpie
        ztmp1(i,j) = bgct2d(i,j,jn2flux)*dlxp(i,j)*dlyp(i,j)
      ENDDO
      ENDDO

      CALL xcsum(sn2flux,ztmp1,ips)

      ztmp1(:,:)=0.0
      DO j=1,kpje
      DO i=1,kpie
        ztmp1(i,j) = bgct2d(i,j,jn2oflux)*dlxp(i,j)*dlyp(i,j)
      ENDDO
      ENDDO

      CALL xcsum(sn2oflux,ztmp1,ips)

      ztmp1(:,:)=0.0
      DO j=1,kpje
      DO i=1,kpie
        ztmp1(i,j) = dlxp(i,j)*dlyp(i,j)
      ENDDO
      ENDDO

      CALL xcsum(ztotarea,ztmp1,ips)

      ztmp1(:,:)=0.0
      DO j=1,kpje
      DO i=1,kpie
        ztmp1(i,j) = atm(i,j,iatmco2)*dlxp(i,j)*dlyp(i,j)
      ENDDO
      ENDDO

      CALL xcsum(zatmco2,ztmp1,ips)

#if defined(BOXATM)
      ztmp1(:,:)=0.0
      DO j=1,kpje
      DO i=1,kpie
        ztmp1(i,j) = atm(i,j,iatmo2)*dlxp(i,j)*dlyp(i,j)
      ENDDO
      ENDDO

      CALL xcsum(zatmo2,ztmp1,ips)

      ztmp1(:,:)=0.0
      DO j=1,kpje
      DO i=1,kpie
        ztmp1(i,j) = atm(i,j,iatmn2)*dlxp(i,j)*dlyp(i,j)
      ENDDO
      ENDDO

      CALL xcsum(zatmn2,ztmp1,ips)
#endif

!      IF (mnproc.eq.1) THEN
!      WRITE(io_stdo_bgc,*) ' '
!      WRITE(io_stdo_bgc,*) 'CO2Flux  :',co2flux
!      WRITE(io_stdo_bgc,*) 'O2 Flux  :',so2flux
!      WRITE(io_stdo_bgc,*) 'N2 Flux  :',sn2flux
!      WRITE(io_stdo_bgc,*) 'N2O Flux :',sn2oflux
!      WRITE(io_stdo_bgc,*) ' '
#if defined(BOXATM)	      
!      WRITE(io_stdo_bgc,*) 'global atm. CO2[ppm] / kmol: ',          &
!     &                               zatmco2/ztotarea,zatmco2*ppm2con       
!      WRITE(io_stdo_bgc,*) 'global atm. O2[ppm] / kmol : ',          &
!     &                               zatmo2/ztotarea,zatmo2*ppm2con 
!      WRITE(io_stdo_bgc,*) 'global atm. N2[ppm] / kmol : ',          &
!     &                               zatmn2/ztotarea,zatmn2*ppm2con 
!      ENDIF
     
#endif

! Complete sum of inventory in between bgc.f90

      ztmp1(:,:)=0.0
      DO j=1,kpje
      DO i=1,kpie
        ztmp1(i,j) = prorca(i,j)*WETO(I,J)*dlxp(i,j)*dlyp(i,j)
      ENDDO
      ENDDO

      CALL xcsum(zprorca,ztmp1,ips)

      ztmp1(:,:)=0.0
      DO j=1,kpje
      DO i=1,kpie
        ztmp1(i,j) = prcaca(i,j)*WETO(I,J)*dlxp(i,j)*dlyp(i,j)
      ENDDO
      ENDDO

      CALL xcsum(zprcaca,ztmp1,ips)

      ztmp1(:,:)=0.0
      DO j=1,kpje
      DO i=1,kpie
        ztmp1(i,j) = silpro(i,j)*WETO(I,J)*dlxp(i,j)*dlyp(i,j)
      ENDDO
      ENDDO

      CALL xcsum(zsilpro,ztmp1,ips)

!      IF (mnproc.eq.1) THEN
!      WRITE(io_stdo_bgc,*) ' '  
!      WRITE(io_stdo_bgc,*) 'Should be zero at the end: '      
!      WRITE(io_stdo_bgc,*) 'prorca, prcaca, silpro  ',                &
!     &       zprorca, zprcaca, zsilpro  
!      WRITE(io_stdo_bgc,*) ' '  
!      ENDIF


! ppm2con: atmospheric weight: ~10000kg/m^2, avrg. ~29 g/mol
! --> 350 kmol/m^2 --> 1ppm ~ 0.35e-3 kmol/m^2

! Sum of inventory
!----------------------------------------------------------------------
! Units in P have a C:P Ratio of 122:1  
  
!      totalcarbon=                                                    &
!     & (zocetrato(idet)+zocetrato(idoc)+zocetrato(iphy)               &
!     & +zocetrato(izoo))*rcar+zocetrato(isco212)+zocetrato(icalc)      

      totalcarbon=                                                    &
     & (zocetrato(idet)+zocetrato(idoc)+zocetrato(iphy)               &
     & +zocetrato(izoo))*rcar+zocetrato(isco212)+zocetrato(icalc)     & 
     & +zpowtrato(ipowaic)+zsedlayto(isssc12)+zsedlayto(issso12)*rcar &
     & +zburial(isssc12)+zburial(issso12)*rcar+zprorca*rcar+zprcaca   & 
     & +zatmco2*ppm2con  

      totalnitr=                                                      &
     &   (zocetrato(idet)+zocetrato(idoc)+zocetrato(iphy)             &
     &  +zocetrato(izoo))*rnit+zocetrato(iano3)+zocetrato(igasnit)*2  &
     &  +zpowtrato(ipowno3)+zpowtrato(ipown2)*2                       &
     &  +zsedlayto(issso12)*rnit+zburial(issso12)*rnit                &
     &  +zocetrato(ian2o)*2                                           &
     &  +zprorca*rnit                                                 &
#if defined(BOXATM)
    &  +zatmn2*ppm2con*2                        
#else
     & +sn2flux*2+sn2oflux*2
#endif     

      totalphos=                                                      &
     &   zocetrato(idet)+zocetrato(idoc)+zocetrato(iphy)              &
     &  +zocetrato(izoo)+zocetrato(iphosph)                           & 
     &  +zpowtrato(ipowaph)+zsedlayto(issso12)                        &
     &  +zburial(issso12)                                             &
     &  +zprorca

      totalsil=                                                       &
     &   zocetrato(isilica)+zocetrato(iopal)                          & 
     &  +zpowtrato(ipowasi)+zsedlayto(issssil)+zburial(issssil)       &
     &  +zsilpro

      totaloxy=                                                       &
     &  (zocetrato(idet)+zocetrato(idoc)+zocetrato(iphy)              &
     &  +zocetrato(izoo))*(-24.)+zocetrato(ioxygen)                   &
     &  +zocetrato(iphosph)*2 +zocetrato(isco212)+zocetrato(icalc)    &
     &  +zocetrato(iano3)*1.5+zocetrato(ian2o)*0.5                    & 
     &  +zsedlayto(issso12)*(-24.) + zsedlayto(isssc12)               &
!     &  +zburial(issso12)*(-24.)   +   zburial(isssc12)               &
     &  +zpowtrato(ipowno3)*1.5+zpowtrato(ipowaic)                    &
     &  +zpowtrato(ipowaox)+zpowtrato(ipowaph)*2                      &
     &  +zprorca*(-24.)+zprcaca                                       & 
#if defined(BOXATM)
     &  +zatmo2*ppm2con+zatmco2*ppm2con
#else
     & +so2flux+sn2oflux*0.5+co2flux
#endif     

      IF (mnproc.eq.1) THEN
!      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global total[kmol] of carbon   : ',       &
     & totalcarbon
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global total[kmol] of phosph.  : ',       &
     & totalphos
      WRITE(io_stdo_bgc,*) ' '	    
      WRITE(io_stdo_bgc,*) 'Global total[kmol] of silicate : ',       &
     & totalsil
      WRITE(io_stdo_bgc,*) ' '	    
      WRITE(io_stdo_bgc,*) 'Global total[kmol] of nitrogen.  : ',     &
     & totalnitr
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global total[kmol] of oxygen.  : ',     &
     & totaloxy
      ENDIF
  
! Write sediment fluxes

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global fluxes into and out of the sediment'
      WRITE(io_stdo_bgc,*) '------------------------------------------'
      WRITE(io_stdo_bgc,*) '        [kmol]'
      ENDIF

      ztmp1(:,:)=0.0
      DO j=1,kpje
      DO i=1,kpie
        ztmp1(i,j) = bgct2d(i,j,jprorca)*WETO(I,J)*dlxp(i,j)*dlyp(i,j)
      ENDDO
      ENDDO

      CALL xcsum(zprorca,ztmp1,ips)

      ztmp1(:,:)=0.0
      DO j=1,kpje
      DO i=1,kpie
        ztmp1(i,j) = bgct2d(i,j,jprcaca)*WETO(I,J)*dlxp(i,j)*dlyp(i,j)
      ENDDO
      ENDDO

      CALL xcsum(zprcaca,ztmp1,ips)

      ztmp1(:,:)=0.0
      DO j=1,kpje
      DO i=1,kpie
        ztmp1(i,j) = bgct2d(i,j,jsilpro)*WETO(I,J)*dlxp(i,j)*dlyp(i,j)
      ENDDO
      ENDDO

      CALL xcsum(zsilpro,ztmp1,ips)

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*) ' '      
      WRITE(io_stdo_bgc,*) 'Detritus, Calcium Carbonate, Silicate  ', &
     &       zprorca, zprcaca, zsilpro  
      WRITE(io_stdo_bgc,*) ' '      
      ENDIF
      
      DO l=1,npowtra

        ztmp1(:,:)=0.0
        DO j=1,kpje
        DO i=1,kpie
          ztmp1(i,j) = sedfluxo(i,j,l)*dlxp(i,j)*dlyp(i,j)
        ENDDO
        ENDDO

        CALL xcsum(sum,ztmp1,ips)
      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*) 'No. ',l,' ',sum
      ENDIF
      ENDDO

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global total export production'
      WRITE(io_stdo_bgc,*) '------------------------------'
      WRITE(io_stdo_bgc,*) '        [kmol]'
      ENDIF

      ztmp1(:,:)=0.0
      DO j=1,kpje
      DO i=1,kpie
        ztmp1(i,j) = expoor(i,j)*dlxp(i,j)*dlyp(i,j)
      ENDDO
      ENDDO

      CALL xcsum(sum,ztmp1,ips)

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*) 'carbon   : ',sum
      ENDIF

      ztmp1(:,:)=0.0
      DO j=1,kpje
      DO i=1,kpie
        ztmp1(i,j) = expoca(i,j)*dlxp(i,j)*dlyp(i,j)
      ENDDO
      ENDDO

      CALL xcsum(sum,ztmp1,ips)

     IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*) 'carbonate: ',sum
      ENDIF

      ztmp1(:,:)=0.0
      DO j=1,kpje
      DO i=1,kpie
        ztmp1(i,j) = exposi(i,j)*dlxp(i,j)*dlyp(i,j)
      ENDDO
      ENDDO

      CALL xcsum(sum,ztmp1,ips)

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*) 'silitat  : ',sum
      WRITE(io_stdo_bgc,*) ' '
      ENDIF

      RETURN
      END
