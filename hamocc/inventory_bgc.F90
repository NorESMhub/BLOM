! Copyright (C) 2002  P. Wetzel
! Copyright (C) 2022  K. Assmann, J. Tjiputra, J. Schwinger, T. Torsvik
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


SUBROUTINE INVENTORY_BGC(kpie,kpje,kpke,dlxp,dlyp,ddpo,omask,iogrp)
!*******************************************************************
!
!**** *INVENTORY_BGC* - calculate the BGC inventory.
!
!     P.Wetzel,              *MPI-Met, HH*    29.07.02
!
!     Modified
!     --------
!     T. Torsvik             *UiB*            22.02.22
!        Include option for writing inventory to netCDF file.
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

  use mo_carbch
  use mo_sedmnt
  use mo_biomod
  use mo_control_bgc
  use mo_bgcmean
  use mo_param1_bgc
  use mo_vgrid,       only: dp_min
  use mod_xc,         only: ips, mnproc, nbdy, xcsum
  use mo_riverinpt,   only: rivinflx,irdin,irdip,irsi,iralk,irdoc,irdet,nriv

  implicit none

  integer, intent(in) :: kpie,kpje,kpke
  integer, intent(in) :: iogrp
  real,    intent(in) :: dlxp(kpie,kpje)
  real,    intent(in) :: dlyp(kpie,kpje)
  real,    intent(in) :: ddpo(kpie,kpje,kpke)
  real,    intent(in) :: omask(kpie,kpje)

  ! Local variables
  integer :: i,j,k,l

  real :: ztmp1(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy)
  real :: ztmp2(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy)
  real :: vol

  ! ppm2con: atmospheric weight: ~10000kg/m^2, avrg. ~29 g/mol
  ! --> 350 kmol/m^2 --> 1ppm ~ 0.35e-3 kmol/m^2
  real, parameter :: ppm2con = 0.35e-3

  !=== Variables for global sums
  real :: ztotvol                  ! Total ocean volume
  real :: ztotarea                 ! Total sea surface area
  !--- aqueous sediment tracer
  real :: zsedtotvol               ! Total pore water volume
  real :: zpowtrato(npowtra)       ! Sum : Pore water tracers
  !--- non aqueous sediment tracer
  real :: zsedhplto                ! Total sediment accumulated hydrogen ions
  real :: zsedlayto(nsedtra)       ! Sum : Sediment layer tracers
  real :: zburial(nsedtra)         ! Sum : Sediment burial tracers
  !--- oceanic tracers
  real :: zocetratot(nocetra)      ! Sum : Ocean tracers
  real :: zocetratoc(nocetra)      ! Mean concentration of ocean racers
  !--- additional ocean tracer
  real :: zhito                    ! Total hydrogen ion tracer
  real :: zco3to                   ! Total dissolved carbonate (CO3) tracer
  !--- alkalinity of the first layer
  real :: zvoltop                  ! Total volume of top ocean layer
  real :: zalkali                  ! Total alkalinity of top ocean layer
  !--- river fluxes
  real :: srivflux(nriv)           ! sum of riverfluxes
  !--- atmosphere flux and atmospheric CO2
  real :: sndepflux                ! sum of N dep fluxes
  real :: zatmco2,zatmo2,zatmn2
  real :: co2flux,so2flux,sn2flux,sn2oflux
  real :: zprorca,zprcaca,zsilpro
  !--- total tracer budgets
  real :: totalcarbon,totalphos,totalsil,totalnitr,totaloxy
  !--- sediment fluxes
  real :: sum_zprorca
  real :: sum_zprcaca
  real :: sum_zsilpro
  real :: sum_sedfluxo(npowtra)
  !--- export production
  real :: sum_expoor
  real :: sum_expoca
  real :: sum_exposi

!=== aqueous sediment tracer
!----------------------------------------------------------------------
#ifdef sedbypass
  zsedtotvol = 0.0
  zpowtrato(:)=0.0
  zsedlayto(:)=0.0
  zburial(:)=0.0
  zsedhplto=0.0
#else
  ztmp1(:,:)=0.0
  DO k=1,ks
  DO j=1,kpje
  DO i=1,kpie
     ztmp1(i,j) = ztmp1(i,j) + omask(i,j)*seddw(k)                             &
          &       *dlxp(i,j)*dlyp(i,j)*porwat(k)
  ENDDO
  ENDDO
  ENDDO

  CALL xcsum(zsedtotvol,ztmp1,ips)

  DO l=1,npowtra
     ztmp1(:,:)=0.0
     DO k=1,ks
     DO j=1,kpje
     DO i=1,kpie
        vol    = seddw(k)*dlxp(i,j)*dlyp(i,j)*porwat(k)
        ztmp1(i,j)= ztmp1(i,j) + omask(i,j)*powtra(i,j,k,l)*vol
     ENDDO
     ENDDO
     ENDDO

     CALL xcsum(zpowtrato(l),ztmp1,ips)
  ENDDO

!=== non aqueous sediment tracer
!----------------------------------------------------------------------
  zburial = sum2d_array(burial, nsedtra)

  DO l=1,nsedtra
     ztmp1(:,:)=0.0
     DO k=1,ks
     DO j=1,kpje
     DO i=1,kpie
        vol = porsol(k)*seddw(k)*dlxp(i,j)*dlyp(i,j)
        ztmp1(i,j) = ztmp1(i,j) + omask(i,j)*sedlay(i,j,k,l)*vol
     ENDDO
     ENDDO
     ENDDO

     CALL xcsum(zsedlayto(l),ztmp1,ips)
  ENDDO

  ztmp1(:,:)=0.0
  DO k=1,ks
  DO j=1,kpje
  DO i=1,kpie
     vol = porsol(k)*seddw(k)*dlxp(i,j)*dlyp(i,j)
     ztmp1(i,j) = ztmp1(i,j) + omask(i,j)*sedhpl(i,j,k)*vol
  ENDDO
  ENDDO
  ENDDO

  CALL xcsum(zsedhplto,ztmp1,ips)
#endif

!=== oceanic tracers
!----------------------------------------------------------------------
  ztotvol    = 0.
  zocetratot = 0.
  zocetratoc = 0.

  ztmp1(:,:)=0.0
  DO k=1,kpke
  DO j=1,kpje
  DO i=1,kpie
     IF(ddpo(i,j,k).gt.dp_min) THEN
        ztmp1(i,j) = ztmp1(i,j) + omask(i,j)*dlxp(i,j)*dlyp(i,j)*ddpo(i,j,k)
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
           vol = dlxp(i,j)*dlyp(i,j)*ddpo(i,j,k)
           ztmp1(i,j) = ztmp1(i,j) + omask(i,j)*ocetra(i,j,k,l)*vol
!             if (ocetra(i,j,k,l).lt.0.0) then
!      WRITE(io_stdo_bgc,*) 'ocetra -ve', l,ocetra(i,j,k,l)
!             endif
        ENDIF
     ENDDO
     ENDDO
     ENDDO

     CALL xcsum(zocetratot(l),ztmp1,ips)
     zocetratoc(l) = zocetratot(l)/ztotvol
  ENDDO

!=== additional ocean tracer
!----------------------------------------------------------------------
  zhito  = 0.
  zco3to = 0.

  ztmp1(:,:)=0.0
  ztmp2(:,:)=0.0
  DO k=1,kpke
  DO j=1,kpje
  DO i=1,kpie
     IF(ddpo(i,j,k).gt.dp_min) THEN
        vol = dlxp(i,j)*dlyp(i,j)*ddpo(i,j,k)
        ztmp1(i,j) = ztmp1(i,j) + omask(i,j)*hi(i,j,k) *vol
        ztmp2(i,j) = ztmp2(i,j) + omask(i,j)*co3(i,j,k)*vol
     ENDIF
  ENDDO
  ENDDO
  ENDDO

  CALL xcsum(zhito ,ztmp1,ips)
  CALL xcsum(zco3to,ztmp2,ips)

!=== alkalinity of the first layer
!--------------------------------------------------------------------
  zvoltop = 0.
  zalkali = 0.

  k=1
  ztmp1(:,:)=0.0
  ztmp2(:,:)=0.0
  DO j=1,kpje
  DO i=1,kpie
     ztmp1(i,j) = omask(i,j)*dlxp(i,j)*dlyp(i,j)*ddpo(i,j,k)
     ztmp2(i,j) = ocetra(i,j,k,ialkali)*ztmp1(i,j)
  ENDDO
  ENDDO

  CALL xcsum(zvoltop,ztmp1,ips)
  CALL xcsum(zalkali,ztmp2,ips)

!=== atmosphere flux and atmospheric CO2
!--------------------------------------------------------------------
  ztotarea =0.
  co2flux  =0.
  so2flux  =0.
  sn2flux  =0.
  sn2oflux =0.
  sndepflux=0.
  srivflux =0.
  zatmco2  =0.
  zatmo2   =0.
  zatmn2   =0.

  ztmp1(:,:)=0.0
  DO j=1,kpje
     DO i=1,kpie
        ztmp1(i,j) = dlxp(i,j)*dlyp(i,j)
     ENDDO
  ENDDO
  CALL xcsum(ztotarea,ztmp1,ips)

  co2flux = sum2d(bgct2d(:,:,jco2flux))
  so2flux = sum2d(bgct2d(:,:,jo2flux))
  sn2flux = sum2d(bgct2d(:,:,jn2flux))
  sn2oflux = sum2d(bgct2d(:,:,jn2oflux))
  zatmco2 = sum2d(atm(:,:,iatmco2))
#if defined(BOXATM)
  zatmo2 = sum2d(atm(:,:,iatmo2))
  zatmn2 = sum2d(atm(:,:,iatmn2))
#endif
  if(do_ndep) then
     sndepflux = sum2d(bgct2d(:,:,jndep))
  endif

!--- Complete sum of inventory in between bgc.f90
  zprorca = sum2d(prorca)
  zprcaca = sum2d(prcaca)
  zsilpro = sum2d(silpro)

!------------------------riverine fluxes
  IF(do_rivinpt)THEN
     DO l=1,nriv
        ztmp1(:,:)=0.0
        DO j=1,kpje
        DO i=1,kpie
           ztmp1(i,j) = bgct2d(i,j,jirdin+l-1)*dlxp(i,j)*dlyp(i,j)
        ENDDO
        ENDDO
        CALL xcsum(srivflux(l),ztmp1,ips)
     ENDDO
  ENDIF

!=== Sum of inventory
!----------------------------------------------------------------------
! Units in P have a C:P Ratio of 122:1

!      totalcarbon=                                                    &
!     & (zocetratot(idet)+zocetratot(idoc)+zocetratot(iphy)               &
!     & +zocetratot(izoo))*rcar+zocetratot(isco212)+zocetratot(icalc)

  totalcarbon=                                                                 &
       & (zocetratot(idet)+zocetratot(idoc)+zocetratot(iphy)                   &
       & +zocetratot(izoo))*rcar+zocetratot(isco212)+zocetratot(icalc)         &
       & +zpowtrato(ipowaic)+zsedlayto(isssc12)+zsedlayto(issso12)*rcar        &
       & +zburial(isssc12)+zburial(issso12)*rcar                               &
#if defined(BOXATM)
       & +zatmco2*ppm2con
#else
       & +co2flux
#endif

  totalnitr=                                                                   &
       &   (zocetratot(idet)+zocetratot(idoc)+zocetratot(iphy)                 &
       &  +zocetratot(izoo))*rnit+zocetratot(iano3)+zocetratot(igasnit)*2      &
       &  +zpowtrato(ipowno3)+zpowtrato(ipown2)*2                              &
       &  +zsedlayto(issso12)*rnit+zburial(issso12)*rnit                       &
       &  +zocetratot(ian2o)*2                                                 &
       &  - sndepflux                                                          &
#if defined(BOXATM)
       &  +zatmn2*ppm2con*2
#else
       &  +sn2flux*2+sn2oflux*2
#endif

  totalphos=                                                                   &
       &   zocetratot(idet)+zocetratot(idoc)+zocetratot(iphy)                  &
       &  +zocetratot(izoo)+zocetratot(iphosph)                                &
       &  +zpowtrato(ipowaph)+zsedlayto(issso12)                               &
       &  +zburial(issso12)

  totalsil=                                                                    &
       &   zocetratot(isilica)+zocetratot(iopal)                               &
       &  +zpowtrato(ipowasi)+zsedlayto(issssil)+zburial(issssil)

  totaloxy=                                                                    &
       &  (zocetratot(idet)+zocetratot(idoc)+zocetratot(iphy)                  &
       &  +zocetratot(izoo))*(-24.)+zocetratot(ioxygen)                        &
       &  +zocetratot(iphosph)*2 +zocetratot(isco212)+zocetratot(icalc)        &
       &  +zocetratot(iano3)*1.5+zocetratot(ian2o)*0.5                         &
       &  +zsedlayto(issso12)*(-24.) + zsedlayto(isssc12)                      &
!       &  +zburial(issso12)*(-24.)   +   zburial(isssc12)                      &
       &  +zpowtrato(ipowno3)*1.5+zpowtrato(ipowaic)                           &
       &  +zpowtrato(ipowaox)+zpowtrato(ipowaph)*2                             &
       &  - sndepflux*1.5                                                      &
#if defined(BOXATM)
       &  +zatmo2*ppm2con+zatmco2*ppm2con
#else
       &  +so2flux+sn2oflux*0.5+co2flux
#endif

  IF (do_rivinpt) THEN
     totalcarbon = totalcarbon - (srivflux(irdoc)+srivflux(irdet))*rcar        &
          &             -(srivflux(iralk)+srivflux(irdin)+srivflux(irdip))     ! =sco212
     totalnitr = totalnitr - (srivflux(irdoc)+srivflux(irdet))*rnit            &
          &             - srivflux(irdin)
     totalphos = totalphos                                                     &
          &             -(srivflux(irdoc)+srivflux(irdet)+srivflux(irdip))
     totalsil = totalsil - srivflux(irsi)
     totaloxy = totaloxy - (srivflux(irdoc)+srivflux(irdet))*(-24.)            &
          &             - srivflux(irdin)*1.5 - srivflux(irdip)*2.             &
          &             - (srivflux(iralk)+srivflux(irdin)+srivflux(irdip))    ! =sco212
  ENDIF

!=== Compute sediment fluxes
!----------------------------------------------------------------------
  sum_zprorca = sum2d(bgct2d(:,:,jprorca))
  sum_zprcaca = sum2d(bgct2d(:,:,jprcaca))
  sum_zsilpro = sum2d(bgct2d(:,:,jsilpro))

  sum_sedfluxo = sum2d_array(sedfluxo, npowtra)

  sum_expoor = sum2d(expoor)
  sum_expoca = sum2d(expoca)
  sum_exposi = sum2d(exposi)

!=== Write output to netCDF file or stdout
!----------------------------------------------------------------------
  if (mnproc == 1) then
     if (GLB_INVENTORY(iogrp) == 2) then
        call write_netcdf(iogrp)
     else
        call write_stdout
     end if
  end if

  return

contains


function sum2d(var2d) result(total)
!**********************************************************************
!**** Sum 2D scalar fields
!**********************************************************************
  implicit none
  real, dimension(kpie,kpje), intent(in) :: var2d
  real :: total

  ! Local variables
  integer :: i,j
  !--- input to xcsum require halo indices
  real, dimension(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy) :: ztmp

  ztmp(:,:)=0.0
  do j=1,kpje
     do i=1,kpie
        ztmp(i,j) = var2d(i,j)*dlxp(i,j)*dlyp(i,j)*omask(i,j)
     enddo
  enddo
  call xcsum(total,ztmp,ips)
end function sum2d


function sum2d_array(var3d, narr) result(total)
!**********************************************************************
!**** Sum 2D array fields
!**********************************************************************
  implicit none
  integer, intent(in) :: narr
  real, dimension(kpie,kpje,narr), intent(in) :: var3d
  real, dimension(narr) :: total

  ! Local variables
  integer :: i,j,k
  !--- input to xcsum require halo indices
  real, dimension(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy) :: ztmp

  ztmp(:,:)=0.0
  do k=1,narr
     do j=1,kpje
        do i=1,kpie
           ztmp(i,j) = var3d(i,j,k)*dlxp(i,j)*dlyp(i,j)*omask(i,j)
        enddo
     enddo
     call xcsum(total(k),ztmp,ips)
  end do
end function sum2d_array


subroutine write_stdout
!**********************************************************************
!**** Write inventory to log file.
!**********************************************************************
  implicit none

  integer :: l

#ifndef sedbypass
  !=== aqueous sediment tracer
  !------------------------------------------------------------------
  WRITE(io_stdo_bgc,*) ' '
  WRITE(io_stdo_bgc,*)'Global inventory of aqueous sediment tracer'
  WRITE(io_stdo_bgc,*)'-------------------------------------------'
  WRITE(io_stdo_bgc,*) '       total[kmol]    concentration[mol/L]'
  DO l=1,npowtra
     WRITE(io_stdo_bgc,*)'No. ',l,' ',zpowtrato(l),                            &
          &     '  ',zpowtrato(l)/zsedtotvol,'  ',zsedtotvol
  ENDDO
  WRITE(io_stdo_bgc,*) ' '

  !=== non aqueous sediment tracer
  !------------------------------------------------------------------
  WRITE(io_stdo_bgc,*) ' '
  WRITE(io_stdo_bgc,*)                                                         &
       &     'Global inventory of solid sediment constituents'
  WRITE(io_stdo_bgc,*)                                                         &
       &     '----------------------------------------------------'
  WRITE(io_stdo_bgc,*) '        [kmol]'

  DO l=1,nsedtra
     WRITE(io_stdo_bgc,*) 'Sediment No. ',l,' ', zsedlayto(l)
     WRITE(io_stdo_bgc,*) 'Burial No. ',l,' ', zburial(l)
  ENDDO
  WRITE(io_stdo_bgc,*) 'hpl ', zsedhplto
  WRITE(io_stdo_bgc,*) ' '
#endif

  !=== oceanic tracers
  !------------------------------------------------------------------
  WRITE(io_stdo_bgc,*) ' '
  WRITE(io_stdo_bgc,*) 'Global inventory of advected ocean tracers'
  WRITE(io_stdo_bgc,*) '------------------------------------------'
  WRITE(io_stdo_bgc,*) ' '
  WRITE(io_stdo_bgc,*) '       total[kmol]  concentration[kmol/m^3]'
  WRITE(io_stdo_bgc,*) ' '
  WRITE(io_stdo_bgc,*) 'ztotvol',ztotvol
  DO l=1,nocetra
     WRITE(io_stdo_bgc,*) 'No. ',l, zocetratot(l), zocetratoc(l)
  ENDDO

  !=== additional ocean tracer
  !------------------------------------------------------------------
  ! WRITE(io_stdo_bgc,*) ' '
  ! WRITE(io_stdo_bgc,*) 'Glob. inventory of additional ocean tracer'
  ! WRITE(io_stdo_bgc,*) '------------------------------------------'
  ! WRITE(io_stdo_bgc,*) '      total[kmol]  concentration[kmol/m^3]'
  ! WRITE(io_stdo_bgc,*) ' '
  ! WRITE(io_stdo_bgc,*) ' hi', zhito, zhito/ztotvol
  ! WRITE(io_stdo_bgc,*) ' co3', zco3to, zco3to/ztotvol
  ! WRITE(io_stdo_bgc,*) ' '

  !=== alkalinity of the first layer
  !------------------------------------------------------------------
  ! WRITE(io_stdo_bgc,*) ' '
  ! WRITE(io_stdo_bgc,*) 'Global inventory of first layer alkalinity'
  ! WRITE(io_stdo_bgc,*) '------------------------------------------'
  ! WRITE(io_stdo_bgc,*) ' '
  ! WRITE(io_stdo_bgc,*) '       total[kmol]  concentration[kmol/m^3]'
  ! WRITE(io_stdo_bgc,*) ' '
  ! WRITE(io_stdo_bgc,*) zalkali, zalkali/zvoltop

  !=== atmosphere flux and atmospheric CO2
  !------------------------------------------------------------------
  ! WRITE(io_stdo_bgc,*) ' '
  ! WRITE(io_stdo_bgc,*) 'Global fluxes into atmosphere'
  ! WRITE(io_stdo_bgc,*) '-----------------------------'
  ! WRITE(io_stdo_bgc,*) '        [kmol]'
  ! WRITE(io_stdo_bgc,*) ' '
  ! WRITE(io_stdo_bgc,*) 'CO2Flux  :',co2flux
  ! WRITE(io_stdo_bgc,*) 'O2 Flux  :',so2flux
  ! WRITE(io_stdo_bgc,*) 'N2 Flux  :',sn2flux
  ! WRITE(io_stdo_bgc,*) 'N2O Flux :',sn2oflux
  ! WRITE(io_stdo_bgc,*) ' '
#if defined(BOXATM)
  ! WRITE(io_stdo_bgc,*) 'global atm. CO2[ppm] / kmol: ',                 &
  !      &               zatmco2/ztotarea,zatmco2*ppm2con
  ! WRITE(io_stdo_bgc,*) 'global atm. O2[ppm] / kmol : ',                 &
  !      &               zatmo2/ztotarea,zatmo2*ppm2con
  ! WRITE(io_stdo_bgc,*) 'global atm. N2[ppm] / kmol : ',                 &
  !      &               zatmn2/ztotarea,zatmn2*ppm2con
#endif
  ! WRITE(io_stdo_bgc,*) ' '
  ! WRITE(io_stdo_bgc,*) 'Should be zero at the end: '
  ! WRITE(io_stdo_bgc,*) 'prorca, prcaca, silpro  ',                      &
  !      &               zprorca, zprcaca, zsilpro
  ! WRITE(io_stdo_bgc,*) ' '

  IF(do_ndep) WRITE(io_stdo_bgc,*) 'NdepFlux :',sndepflux

  ! riverine fluxes
  !------------------------------------------------------------------
  IF(do_rivinpt)THEN
     WRITE(io_stdo_bgc,*) 'Riverine fluxes:'
     DO l=1,nriv
        WRITE(io_stdo_bgc,*) 'No. ',l,srivflux(l)
     ENDDO
  ENDIF

  !=== Sum of inventory
  !------------------------------------------------------------------
  ! Units in P have a C:P Ratio of 122:1
  WRITE(io_stdo_bgc,*) 'Global total[kmol] of carbon   : ', totalcarbon
  WRITE(io_stdo_bgc,*) ' '
  WRITE(io_stdo_bgc,*) 'Global total[kmol] of phosph.  : ', totalphos
  WRITE(io_stdo_bgc,*) ' '
  WRITE(io_stdo_bgc,*) 'Global total[kmol] of silicate : ', totalsil
  WRITE(io_stdo_bgc,*) ' '
  WRITE(io_stdo_bgc,*) 'Global total[kmol] of nitrogen. : ', totalnitr
  WRITE(io_stdo_bgc,*) ' '
  WRITE(io_stdo_bgc,*) 'Global total[kmol] of oxygen.  : ', totaloxy

  !=== Write sediment fluxes
  !------------------------------------------------------------------
  WRITE(io_stdo_bgc,*) ' '
  WRITE(io_stdo_bgc,*) 'Global fluxes into and out of the sediment'
  WRITE(io_stdo_bgc,*) '------------------------------------------'
  WRITE(io_stdo_bgc,*) '        [kmol]'
  WRITE(io_stdo_bgc,*) ' '
  WRITE(io_stdo_bgc,*) 'Detritus, Calcium Carbonate, Silicate  ',              &
       &               sum_zprorca, sum_zprcaca, sum_zsilpro
  WRITE(io_stdo_bgc,*) ' '
  DO l=1,npowtra
     WRITE(io_stdo_bgc,*) 'No. ',l,' ',sum_sedfluxo(l)
  ENDDO
  WRITE(io_stdo_bgc,*) ' '
  WRITE(io_stdo_bgc,*) 'Global total export production'
  WRITE(io_stdo_bgc,*) '------------------------------'
  WRITE(io_stdo_bgc,*) '        [kmol]'
  WRITE(io_stdo_bgc,*) 'carbon   : ',sum_expoor
  WRITE(io_stdo_bgc,*) 'carbonate: ',sum_expoca
  WRITE(io_stdo_bgc,*) 'silicate : ',sum_exposi
  WRITE(io_stdo_bgc,*) ' '

end subroutine write_stdout


subroutine write_netcdf(iogrp)
!**********************************************************************
!**** Write inventory to netCDF file.
!**********************************************************************
  use netcdf
  use mod_types,    only: r8
  use mod_config,   only: expcnf, runid, inst_suffix
  use mod_time,     only: date0, time0, date, time, nstep

  implicit none

  integer, intent(in) :: iogrp

  !=== Save filename and counter variables
  !--- netCDF output file names
  character(len=256), dimension(nbgcmax), save :: fname_inv
  integer, dimension(nbgcmax), save :: ncrec = 0

  !=== Local variables
  character(len=:), allocatable :: prefix, sep1, sep2
  character(len=20) :: tstamp
  character(len=30) :: timeunits
  integer :: l
  real(r8) :: datenum

  !=== Variables for netcdf
  integer :: ncid, ncvarid, ncstat
  integer :: wrstart(1)
  !--- dimension id
  integer :: time_dimid
  integer :: nocetra_dimid
  ! integer :: zpowtra_dimids(2)     ! aqueous sediment dimensions
  ! integer :: zsedtra_dimids(2)     ! solid sediments dimensions
  ! integer :: zocetra_dimids(2)     ! ocean tracer dimensions
  ! integer :: zocetra_wrstart(2)
  ! integer :: zocetra_count(2)
  !--- variable id
  integer :: time_varid
#ifndef sedbypass
  !--- aqueous sediment tracers
  integer :: zsedtotvol_varid      ! Total sediment volume
  integer :: zpowtratot_varid      ! Total aqueous sediment tracer [kmol]
  integer :: zpowtratoc_varid      ! Sediment tracer concentration [kmol/L]
  !--- non-aqueous sediment tracers
  integer :: zsedlayto_varid
  integer :: zburial_varid
  integer :: zsedhplto_varid
#endif
  !--- oceanic tracers
  !--- Write total sum zt_<variable>_varid, and mean concentration zc_<variable>_varid
  integer :: ztotvol_varid         ! Total ocean volume
  ! integer :: zocetratot_varid      ! Total ocean tracer [kmol]
  ! integer :: zocetratoc_varid      ! Ocean tracer concentration [kmol/m^3]
  integer :: zt_sco212_varid,    zc_sco212_varid      ! Dissolved CO2
  integer :: zt_alkali_varid,    zc_alkali_varid      ! Alkalinity
  integer :: zt_phosph_varid,    zc_phosph_varid      ! Dissolved phosphate
  integer :: zt_oxygen_varid,    zc_oxygen_varid      ! Dissolved oxygen
  integer :: zt_gasnit_varid,    zc_gasnit_varid      ! Gaseous nitrogen (N2)
  integer :: zt_ano3_varid,      zc_ano3_varid        ! Dissolved nitrate
  integer :: zt_silica_varid,    zc_silica_varid      ! Silicid acid (Si(OH)4)
  integer :: zt_doc_varid,       zc_doc_varid         ! Dissolved organic carbon
  integer :: zt_poc_varid,       zc_poc_varid         ! Particulate organic carbon
  integer :: zt_phyto_varid,     zc_phyto_varid       ! Phytoplankton concentration
  integer :: zt_grazer_varid,    zc_grazer_varid      ! Zooplankton concentration
  integer :: zt_calciu_varid,    zc_calciu_varid      ! Calcium carbonate
  integer :: zt_opal_varid,      zc_opal_varid        ! Biogenic silica
  integer :: zt_n2o_varid,       zc_n2o_varid         ! Laughing gas (N2O)
  integer :: zt_dms_varid,       zc_dms_varid         ! DiMethylSulfide
  integer :: zt_fdust_varid,     zc_fdust_varid       ! Non-aggregated dust
  integer :: zt_iron_varid,      zc_iron_varid        ! Dissolved iron
  integer :: zt_prefo2_varid,    zc_prefo2_varid      ! Preformed oxygen
  integer :: zt_prefpo4_varid,   zc_prefpo4_varid     ! Preformed phosphate
  integer :: zt_prefalk_varid,   zc_prefalk_varid     ! Preformed alkalinity
  integer :: zt_prefdic_varid,   zc_prefdic_varid     ! Preformed DIC
  integer :: zt_dicsat_varid,    zc_dicsat_varid      ! Saturated DIC
#ifdef cisonew
  integer :: zt_sco213_varid,    zc_sco213_varid      ! Dissolved CO2-C13
  integer :: zt_sco214_varid,    zc_sco214_varid      ! Dissolved CO2-C14
  integer :: zt_doc13_varid,     zc_doc13_varid       ! Dissolved organic carbon-C13
  integer :: zt_doc14_varid,     zc_doc14_varid       ! Dissolved organic carbon-C14
  integer :: zt_poc13_varid,     zc_poc13_varid       ! Particulate organic carbon-C13
  integer :: zt_poc14_varid,     zc_poc14_varid       ! Particulate organic carbon-C14
  integer :: zt_phyto13_varid,   zc_phyto13_varid     ! Phytoplankton concentration-C13
  integer :: zt_phyto14_varid,   zc_phyto14_varid     ! Phytoplankton concentration-C14
  integer :: zt_grazer13_varid,  zc_grazer13_varid    ! Zooplankton concentration-C13
  integer :: zt_grazer14_varid,  zc_grazer14_varid    ! Zooplankton concentration-C14
  integer :: zt_calciu13_varid,  zc_calciu13_varid    ! Calcium carbonate-C13
  integer :: zt_calciu14_varid,  zc_calciu14_varid    ! Calcium carbonate-C14
#endif
#ifdef AGG
  integer :: zt_snos_varid,      zc_snos_varid        ! Marine snow aggregates per g sea water
  integer :: zt_adust_varid,     zc_adust_varid       ! Aggregated dust
#endif
#ifdef CFC
  integer :: zt_cfc11_varid,     zc_cfc11_varid       ! CFC-11 : Trichlorofluoromethane
  integer :: zt_cfc12_varid,     zc_cfc12_varid       ! CFC-12 : Dichlorodifluoromethane
  integer :: zt_sf6_varid,       zc_sf6_varid         ! SF6    : Sulfur hexafluoride
#endif
#ifdef natDIC
  integer :: zt_natsco212_varid, zc_natsco212_varid   ! Natural dissolved CO2
  integer :: zt_natalkali_varid, zc_natalkali_varid   ! Natural alkalinity
  integer :: zt_natcalciu_varid, zc_natcalciu_varid   ! Natural calcium carbonate
#endif
#ifdef BROMO
  integer :: zt_bromo_varid,     zc_bromo_varid       ! Bromoform
#endif
  !--- sum of inventory
  integer :: totcarb_varid, totphos_varid, totsili_varid, totnitr_varid
  integer :: totoxyg_varid
  !--- sediment fluxes
  integer :: sum_zprorca_varid, sum_zprcaca_varid, sum_zsilpro_varid
  integer :: sum_sedfluxo_varid
  integer :: sum_expoor_varid, sum_expoca_varid, sum_exposi_varid


  !=== Create new or open existing netCDF file
  if (ncrec(iogrp) == 0) then
     !--- file name : fname_inv(iogrp)
     if (expcnf.eq.'cesm') then
        prefix=trim(runid)//'.blom'//trim(inst_suffix)
        sep1='.'
        sep2='-'
     else
        prefix=trim(runid)
        sep1='_'
        sep2='.'
     endif
     write(tstamp,'(i4.4,a1,i2.2,a1,i2.2)')                                    &
          &    date%year,sep2,date%month,sep2,date%day
     fname_inv(iogrp) = prefix//sep1//trim(glb_fnametag(iogrp))//sep1//        &
          &    'i'//sep1//trim(tstamp)//'.nc'

     !--- create a new netCDF file
     write(io_stdo_bgc,*) 'Create BGC inventory file : ',trim(fname_inv(iogrp))
     call nccheck( NF90_CREATE(trim(fname_inv(iogrp)), NF90_CLOBBER, ncid) )

     !--- set time information
     timeunits=' '
     write(timeunits,'(a11,i4.4,a1,i2.2,a1,i2.2,a6)')                          &
          &   'days since ',date0%year,'-',date0%month,'-',date0%day,' 00:00'

     !--- Define global attributes
     call nccheck( NF90_PUT_ATT(ncid, NF90_GLOBAL, 'title',                    &
          &    'Global inventory for marine bgc') )
     call nccheck( NF90_PUT_ATT(ncid, NF90_GLOBAL, 'history',                  &
          &    'Global inventory for marine bgc') )
     call nccheck( NF90_PUT_ATT(ncid, NF90_GLOBAL, 'date', timeunits) )

     !--- Define dimensions
     call nccheck( NF90_DEF_DIM(ncid, 'nocetra', nocetra, nocetra_dimid) )
     call nccheck( NF90_DEF_DIM(ncid, 'time', NF90_UNLIMITED, time_dimid) )

     !--- Dimensions for arrays.
     !--- The unlimited "time" dimension must come last in the list of dimensions.
     ! zocetra_dimids = (/ nocetra_dimid, time_dimid /)

     !--- Define variables : time
     call nccheck( NF90_DEF_VAR(ncid, 'time', NF90_DOUBLE, time_dimid,         &
          &    time_varid) )
     call nccheck( NF90_PUT_ATT(ncid, time_varid, 'units', 'days') )

     !--- Define variables : oceanic tracers
     call nccheck( NF90_DEF_VAR(ncid, 'ztotvol', NF90_DOUBLE, time_dimid,      &
          &    ztotvol_varid) )
     call nccheck( NF90_PUT_ATT(ncid, ztotvol_varid, 'long_name',              &
          &    'Total ocean volume') )
     call nccheck( NF90_PUT_ATT(ncid, ztotvol_varid, 'units', 'm^3') )

     ! call nccheck( NF90_DEF_VAR(ncid, 'zocetratot', NF90_DOUBLE,               &
     !      &    zocetra_dimids, zocetratot_varid) )
     ! call nccheck( NF90_PUT_ATT(ncid, zocetratot_varid, 'long_name',           &
     !      &    'Total ocean tracer') )
     ! call nccheck( NF90_PUT_ATT(ncid, zocetratot_varid, 'units', 'kmol') )

     ! call nccheck( NF90_DEF_VAR(ncid, 'zocetratoc', NF90_DOUBLE,               &
     !      &    zocetra_dimids, zocetratoc_varid) )
     ! call nccheck( NF90_PUT_ATT(ncid, zocetratoc_varid, 'long_name',           &
     !      &    'Ocean tracer concentration') )
     ! call nccheck( NF90_PUT_ATT(ncid, zocetratoc_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_sco212', NF90_DOUBLE,                &
          &    time_dimid, zt_sco212_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_sco212_varid, 'long_name',            &
          &    'Total dissolved CO2 tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_sco212_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_sco212', NF90_DOUBLE,                &
          &    time_dimid, zc_sco212_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_sco212_varid, 'long_name',            &
          &    'Mean dissolved CO2 concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_sco212_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_alkali', NF90_DOUBLE,                &
          &    time_dimid, zt_alkali_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_alkali_varid, 'long_name',            &
          &    'Total alkalinity tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_alkali_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_alkali', NF90_DOUBLE,                &
          &    time_dimid, zc_alkali_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_alkali_varid, 'long_name',            &
          &    'Mean alkalinity concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_alkali_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_phosph', NF90_DOUBLE,                &
          &    time_dimid, zt_phosph_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_phosph_varid, 'long_name',            &
          &    'Total dissolved phosphate tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_phosph_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_phosph', NF90_DOUBLE,                &
          &    time_dimid, zc_phosph_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_phosph_varid, 'long_name',            &
          &    'Mean dissolved phosphate concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_phosph_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_oxygen', NF90_DOUBLE,                &
          &    time_dimid, zt_oxygen_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_oxygen_varid, 'long_name',            &
          &    'Total dissolved oxygen tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_oxygen_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_oxygen', NF90_DOUBLE,                &
          &    time_dimid, zc_oxygen_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_oxygen_varid, 'long_name',            &
          &    'Mean dissolved oxygen concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_oxygen_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_gasnit', NF90_DOUBLE,                &
          &    time_dimid, zt_gasnit_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_gasnit_varid, 'long_name',            &
          &    'Total gaseous nitrogen (N2) tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_gasnit_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_gasnit', NF90_DOUBLE,                &
          &    time_dimid, zc_gasnit_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_gasnit_varid, 'long_name',            &
          &    'Mean gaseous nitrogen (N2) concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_gasnit_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_ano3', NF90_DOUBLE,                  &
          &    time_dimid, zt_ano3_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_ano3_varid, 'long_name',              &
          &    'Total dissolved nitrate tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_ano3_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_ano3', NF90_DOUBLE,                  &
          &    time_dimid, zc_ano3_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_ano3_varid, 'long_name',              &
          &    'Mean dissolved nitrate concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_ano3_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_silica', NF90_DOUBLE,                &
          &    time_dimid, zt_silica_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_silica_varid, 'long_name',            &
          &    'Total silicid acid (Si(OH)4) tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_silica_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_silica', NF90_DOUBLE,                &
          &    time_dimid, zc_silica_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_silica_varid, 'long_name',            &
          &    'Mean silicid acid (Si(OH)4) concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_silica_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_doc', NF90_DOUBLE,                   &
          &    time_dimid, zt_doc_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_doc_varid, 'long_name',               &
          &    'Total dissolved organic carbon tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_doc_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_doc', NF90_DOUBLE,                   &
          &    time_dimid, zc_doc_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_doc_varid, 'long_name',               &
          &    'Mean dissolved organic carbon concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_doc_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_poc', NF90_DOUBLE,                   &
          &    time_dimid, zt_poc_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_poc_varid, 'long_name',               &
          &    'Total particulate organic carbon tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_poc_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_poc', NF90_DOUBLE,                   &
          &    time_dimid, zc_poc_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_poc_varid, 'long_name',               &
          &    'Mean particulate organic carbon concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_poc_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_phyto', NF90_DOUBLE,                 &
          &    time_dimid, zt_phyto_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_phyto_varid, 'long_name',             &
          &    'Total phytoplankton tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_phyto_varid, 'units', 'kmolP') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_phyto', NF90_DOUBLE,                 &
          &    time_dimid, zc_phyto_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_phyto_varid, 'long_name',             &
          &    'Mean phytoplankton concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_phyto_varid, 'units', 'kmolP/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_grazer', NF90_DOUBLE,                &
          &    time_dimid, zt_grazer_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_grazer_varid, 'long_name',            &
          &    'Total zooplankton tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_grazer_varid, 'units', 'kmolP') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_grazer', NF90_DOUBLE,                &
          &    time_dimid, zc_grazer_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_grazer_varid, 'long_name',            &
          &    'Mean zooplankton concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_grazer_varid, 'units', 'kmolP/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_calciu', NF90_DOUBLE,                &
          &    time_dimid, zt_calciu_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_calciu_varid, 'long_name',            &
          &    'Total calcium carbonate tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_calciu_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_calciu', NF90_DOUBLE,                &
          &    time_dimid, zc_calciu_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_calciu_varid, 'long_name',            &
          &    'Mean calcium carbonate concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_calciu_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_opal', NF90_DOUBLE,                  &
          &    time_dimid, zt_opal_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_opal_varid, 'long_name',              &
          &    'Total biogenic silica tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_opal_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_opal', NF90_DOUBLE,                  &
          &    time_dimid, zc_opal_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_opal_varid, 'long_name',              &
          &    'Mean biogenic silica concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_opal_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_n2o', NF90_DOUBLE,                   &
          &    time_dimid, zt_n2o_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_n2o_varid, 'long_name',               &
          &    'Total laughing gas (N2O) tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_n2o_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_n2o', NF90_DOUBLE,                   &
          &    time_dimid, zc_n2o_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_n2o_varid, 'long_name',               &
          &    'Mean laughing gas (N2O) concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_n2o_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_dms', NF90_DOUBLE,                   &
          &    time_dimid, zt_dms_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_dms_varid, 'long_name',               &
          &    'Total DiMethylSulfide tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_dms_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_dms', NF90_DOUBLE,                   &
          &    time_dimid, zc_dms_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_dms_varid, 'long_name',               &
          &    'Mean DiMethylSulfide concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_dms_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_fdust', NF90_DOUBLE,                 &
          &    time_dimid, zt_fdust_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_fdust_varid, 'long_name',             &
          &    'Total non-aggregated dust tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_fdust_varid, 'units', 'Mg') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_fdust', NF90_DOUBLE,                 &
          &    time_dimid, zc_fdust_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_fdust_varid, 'long_name',             &
          &    'Mean non-aggregate dust concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_fdust_varid, 'units', 'Mg/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_iron', NF90_DOUBLE,                  &
          &    time_dimid, zt_iron_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_iron_varid, 'long_name',              &
          &    'Total dissolved iron tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_iron_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_iron', NF90_DOUBLE,                  &
          &    time_dimid, zc_iron_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_iron_varid, 'long_name',              &
          &    'Mean dissolved iron concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_iron_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_prefo2', NF90_DOUBLE,                &
          &    time_dimid, zt_prefo2_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_prefo2_varid, 'long_name',            &
          &    'Total preformed oxygen tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_prefo2_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_prefo2', NF90_DOUBLE,                &
          &    time_dimid, zc_prefo2_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_prefo2_varid, 'long_name',            &
          &    'Mean preformed oxygen concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_prefo2_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_prefpo4', NF90_DOUBLE,               &
          &    time_dimid, zt_prefpo4_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_prefpo4_varid, 'long_name',           &
          &    'Total preformed phosphate tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_prefpo4_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_prefpo4', NF90_DOUBLE,               &
          &    time_dimid, zc_prefpo4_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_prefpo4_varid, 'long_name',           &
          &    'Mean preformed phosphate concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_prefpo4_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_prefalk', NF90_DOUBLE,               &
          &    time_dimid, zt_prefalk_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_prefalk_varid, 'long_name',           &
          &    'Total preformed alkalinity tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_prefalk_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_prefalk', NF90_DOUBLE,               &
          &    time_dimid, zc_prefalk_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_prefalk_varid, 'long_name',           &
          &    'Mean preformed alkalinity concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_prefalk_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_prefdic', NF90_DOUBLE,               &
          &    time_dimid, zt_prefdic_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_prefdic_varid, 'long_name',           &
          &    'Total preformed DIC tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_prefdic_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_prefdic', NF90_DOUBLE,               &
          &    time_dimid, zc_prefdic_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_prefdic_varid, 'long_name',           &
          &    'Mean preformed DIC concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_prefdic_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_dicsat', NF90_DOUBLE,                &
          &    time_dimid, zt_dicsat_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_dicsat_varid, 'long_name',            &
          &    'Total saturated DIC tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_dicsat_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_dicsat', NF90_DOUBLE,                &
          &    time_dimid, zc_dicsat_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_dicsat_varid, 'long_name',            &
          &    'Mean saturated DIC concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_dicsat_varid, 'units', 'kmol/m^3') )

#ifdef cisonew
     call nccheck( NF90_DEF_VAR(ncid, 'zt_sco213', NF90_DOUBLE,                &
          &    time_dimid, zt_sco213_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_sco213_varid, 'long_name',            &
          &    'Total dissolved CO2-C13 tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_sco213_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_sco213', NF90_DOUBLE,                &
          &    time_dimid, zc_sco213_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_sco213_varid, 'long_name',            &
          &    'Mean dissolved CO2-C13 concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_sco213_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_sco214', NF90_DOUBLE,                &
          &    time_dimid, zt_sco214_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_sco214_varid, 'long_name',            &
          &    'Total dissolved CO2-C14 tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_sco214_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_sco214', NF90_DOUBLE,                &
          &    time_dimid, zc_sco214_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_sco214_varid, 'long_name',            &
          &    'Mean dissolved CO2-C14 concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_sco214_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_doc13', NF90_DOUBLE,                 &
          &    time_dimid, zt_doc13_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_doc13_varid, 'long_name',             &
          &    'Total dissolved organic carbon-C13 tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_doc13_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_doc13', NF90_DOUBLE,                 &
          &    time_dimid, zc_doc13_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_doc13_varid, 'long_name',             &
          &    'Mean dissolved organic carbon-C13 concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_doc13_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_doc14', NF90_DOUBLE,                 &
          &    time_dimid, zt_doc14_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_doc14_varid, 'long_name',             &
          &    'Total dissolved organic carbon-C14 tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_doc14_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_doc14', NF90_DOUBLE,                 &
          &    time_dimid, zc_doc14_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_doc14_varid, 'long_name',             &
          &    'Mean dissolved organic carbon-C14 concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_doc14_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_poc13', NF90_DOUBLE,                 &
          &    time_dimid, zt_poc13_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_poc13_varid, 'long_name',             &
          &    'Total particulate organic carbon-C13 tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_poc13_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_poc13', NF90_DOUBLE,                 &
          &    time_dimid, zc_poc13_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_poc13_varid, 'long_name',             &
          &    'Mean particulate organic carbon-C13 concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_poc13_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_poc14', NF90_DOUBLE,                 &
          &    time_dimid, zt_poc14_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_poc14_varid, 'long_name',             &
          &    'Total particulate organic carbon-C14 tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_poc14_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_poc14', NF90_DOUBLE,                 &
          &    time_dimid, zc_poc14_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_poc14_varid, 'long_name',             &
          &    'Mean particulate organic carbon-C14 concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_poc14_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_phyto13', NF90_DOUBLE,               &
          &    time_dimid, zt_phyto13_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_phyto13_varid, 'long_name',           &
          &    'Total phytoplankton-C13 tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_phyto13_varid, 'units', 'kmolP') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_phyto13', NF90_DOUBLE,               &
          &    time_dimid, zc_phyto13_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_phyto13_varid, 'long_name',           &
          &    'Mean phytoplankton-C13 concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_phyto13_varid, 'units', 'kmolP/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_phyto14', NF90_DOUBLE,               &
          &    time_dimid, zt_phyto14_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_phyto14_varid, 'long_name',           &
          &    'Total phytoplankton-C14 tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_phyto14_varid, 'units', 'kmolP') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_phyto14', NF90_DOUBLE,               &
          &    time_dimid, zc_phyto14_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_phyto14_varid, 'long_name',           &
          &    'Mean phytoplankton-C14 concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_phyto14_varid, 'units', 'kmolP/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_grazer13', NF90_DOUBLE,              &
          &    time_dimid, zt_grazer13_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_grazer13_varid, 'long_name',          &
          &    'Total zooplankton-C13 tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_grazer13_varid, 'units', 'kmolP') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_grazer13', NF90_DOUBLE,              &
          &    time_dimid, zc_grazer13_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_grazer13_varid, 'long_name',          &
          &    'Mean zooplankton-C13 concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_grazer13_varid, 'units',              &
          &    'kmolP/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_grazer14', NF90_DOUBLE,              &
          &    time_dimid, zt_grazer14_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_grazer14_varid, 'long_name',          &
          &    'Total zooplankton-C14 tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_grazer14_varid, 'units', 'kmolP') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_grazer14', NF90_DOUBLE,              &
          &    time_dimid, zc_grazer14_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_grazer14_varid, 'long_name',          &
          &    'Mean zooplankton-C14 concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_grazer14_varid, 'units',              &
          &    'kmolP/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_calciu13', NF90_DOUBLE,              &
          &    time_dimid, zt_calciu13_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_calciu13_varid, 'long_name',          &
          &    'Total calcium carbonate-C13 tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_calciu13_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_calciu13', NF90_DOUBLE,              &
          &    time_dimid, zc_calciu13_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_calciu13_varid, 'long_name',          &
          &    'Mean calcium carbonate-C13 concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_calciu13_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_calciu14', NF90_DOUBLE,              &
          &    time_dimid, zt_calciu14_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_calciu14_varid, 'long_name',          &
          &    'Total calcium carbonate-C14 tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_calciu14_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_calciu14', NF90_DOUBLE,              &
          &    time_dimid, zc_calciu14_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_calciu14_varid, 'long_name',          &
          &    'Mean calcium carbonate-C14 concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_calciu14_varid, 'units', 'kmol/m^3') )
#endif

#ifdef AGG
     call nccheck( NF90_DEF_VAR(ncid, 'zt_snos', NF90_DOUBLE,                  &
          &    time_dimid, zt_snos_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_snos_varid, 'long_name',              &
          &    'Total marine snow aggrerates tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_snos_varid, 'units', '---') )           ! What is the unit?

     call nccheck( NF90_DEF_VAR(ncid, 'zc_snos', NF90_DOUBLE,                  &
          &    time_dimid, zc_snos_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_snos_varid, 'long_name',              &
          &    'Mean marine snow aggregates concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_snos_varid, 'units', '---/m^3') )       ! What is the unit?

     call nccheck( NF90_DEF_VAR(ncid, 'zt_adust', NF90_DOUBLE,                 &
          &    time_dimid, zt_adust_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_adust_varid, 'long_name',             &
          &    'Total aggregated dust tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_adust_varid, 'units', '---') )          ! What is the unit?

     call nccheck( NF90_DEF_VAR(ncid, 'zc_adust', NF90_DOUBLE,                 &
          &    time_dimid, zc_adust_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_adust_varid, 'long_name',             &
          &    'Mean aggregated dust concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_adust_varid, 'units', '---/m^3') )      ! What is the unit?
#endif

#ifdef CFC
     call nccheck( NF90_DEF_VAR(ncid, 'zt_cfc11', NF90_DOUBLE,                 &
          &    time_dimid, zt_cfc11_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_cfc11_varid, 'long_name',             &
          &    'Total CFC-11 tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_cfc11_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_cfc11', NF90_DOUBLE,                 &
          &    time_dimid, zc_cfc11_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_cfc11_varid, 'long_name',             &
          &    'Mean CFC-11 concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_cfc11_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_cfc12', NF90_DOUBLE,                 &
          &    time_dimid, zt_cfc12_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_cfc12_varid, 'long_name',             &
          &    'Total CFC-12 tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_cfc12_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_cfc12', NF90_DOUBLE,                 &
          &    time_dimid, zc_cfc12_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_cfc12_varid, 'long_name',             &
          &    'Mean CFC-12 concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_cfc12_varid, 'units', 'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_sf6', NF90_DOUBLE,                   &
          &    time_dimid, zt_sf6_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_sf6_varid, 'long_name',               &
          &    'Total SF6 tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_sf6_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_sf6', NF90_DOUBLE,                   &
          &    time_dimid, zc_sf6_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_sf6_varid, 'long_name',               &
          &    'Mean SF6 concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_sf6_varid, 'units', 'kmol/m^3') )
#endif

#ifdef natDIC
     call nccheck( NF90_DEF_VAR(ncid, 'zt_natsco212', NF90_DOUBLE,             &
          &    time_dimid, zt_natsco212_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_natsco212_varid, 'long_name',         &
          &    'Total natural dissolved CO2 tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_natsco212_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_natsco212', NF90_DOUBLE,             &
          &    time_dimid, zc_natsco212_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_natsco212_varid, 'long_name',         &
          &    'Mean natural dissolved CO2 concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_natsco212_varid, 'units',             &
          &    'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_natalkali', NF90_DOUBLE,             &
          &    time_dimid, zt_natalkali_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_natalkali_varid, 'long_name',         &
          &    'Total natural alkalinity tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_natalkali_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_natalkali', NF90_DOUBLE,             &
          &    time_dimid, zc_natalkali_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_natalkali_varid, 'long_name',         &
          &    'Mean natural alkalinity concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_natalkali_varid, 'units',             &
          &    'kmol/m^3') )

     call nccheck( NF90_DEF_VAR(ncid, 'zt_natcalciu', NF90_DOUBLE,             &
          &    time_dimid, zt_natcalciu_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_natcalciu_varid, 'long_name',         &
          &    'Total natural calcium carbonate tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_natcalciu_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_natcalciu', NF90_DOUBLE,             &
          &    time_dimid, zc_natcalciu_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_natcalciu_varid, 'long_name',         &
          &    'Mean natural calcium carbonate concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_natcalciu_varid, 'units',             &
          &    'kmol/m^3') )
#endif

#ifdef BROMO
     call nccheck( NF90_DEF_VAR(ncid, 'zt_bromo', NF90_DOUBLE,                 &
          &    time_dimid, zt_bromo_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zt_bromo_varid, 'long_name',             &
          &    'Total bromoform tracer') )
     call nccheck( NF90_PUT_ATT(ncid, zt_bromo_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'zc_bromo', NF90_DOUBLE,                 &
          &    time_dimid, zc_bromo_varid) )
     call nccheck( NF90_PUT_ATT(ncid, zc_bromo_varid, 'long_name',             &
          &    'Mean bromoform concentration') )
     call nccheck( NF90_PUT_ATT(ncid, zc_bromo_varid, 'units', 'kmol/m^3') )
#endif

     !--- Define variables : sum of inventory
     call nccheck( NF90_DEF_VAR(ncid, 'totcarb', NF90_DOUBLE, time_dimid,      &
          &    totcarb_varid) )
     call nccheck( NF90_PUT_ATT(ncid, totcarb_varid, 'long_name',              &
          &    'Global total of carbon') )
     call nccheck( NF90_PUT_ATT(ncid, totcarb_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'totphos', NF90_DOUBLE, time_dimid,      &
          &    totphos_varid) )
     call nccheck( NF90_PUT_ATT(ncid, totphos_varid, 'long_name',              &
          &    'Global total of phosphorous') )
     call nccheck( NF90_PUT_ATT(ncid, totphos_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'totsili', NF90_DOUBLE, time_dimid,      &
          &    totsili_varid) )
     call nccheck( NF90_PUT_ATT(ncid, totsili_varid, 'long_name',              &
          &    'Global total of silicate') )
     call nccheck( NF90_PUT_ATT(ncid, totsili_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'totnitr', NF90_DOUBLE, time_dimid,      &
          &    totnitr_varid) )
     call nccheck( NF90_PUT_ATT(ncid, totnitr_varid, 'long_name',              &
          &    'Global total of nitrogen') )
     call nccheck( NF90_PUT_ATT(ncid, totnitr_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'totoxyg', NF90_DOUBLE, time_dimid,      &
          &    totoxyg_varid) )
     call nccheck( NF90_PUT_ATT(ncid, totoxyg_varid, 'long_name',              &
          &    'Global total of oxygen') )
     call nccheck( NF90_PUT_ATT(ncid, totoxyg_varid, 'units', 'kmol') )

     !--- Define variables : sediment fluxes
     call nccheck( NF90_DEF_VAR(ncid, 'sum_zprorca', NF90_DOUBLE,              &
          &    time_dimid, sum_zprorca_varid) )
     call nccheck( NF90_PUT_ATT(ncid, sum_zprorca_varid, 'long_name',          &
          &    'Global flux of detritus into and out of sediments') )
     call nccheck( NF90_PUT_ATT(ncid, sum_zprorca_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'sum_zprcaca', NF90_DOUBLE,              &
          &    time_dimid, sum_zprcaca_varid) )
     call nccheck( NF90_PUT_ATT(ncid, sum_zprcaca_varid, 'long_name',          &
          &    'Global flux of calcium carbonate into and out of sediments') )
     call nccheck( NF90_PUT_ATT(ncid, sum_zprcaca_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'sum_zsilpro', NF90_DOUBLE,              &
          &    time_dimid, sum_zsilpro_varid) )
     call nccheck( NF90_PUT_ATT(ncid, sum_zsilpro_varid, 'long_name',          &
          &    'Global flux of silicate into and out of sediments') )
     call nccheck( NF90_PUT_ATT(ncid, sum_zsilpro_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'sum_expoor', NF90_DOUBLE,               &
          &    time_dimid, sum_expoor_varid) )
     call nccheck( NF90_PUT_ATT(ncid, sum_expoor_varid, 'long_name',           &
          &    'Global total export production of carbon') )
     call nccheck( NF90_PUT_ATT(ncid, sum_expoor_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'sum_expoca', NF90_DOUBLE,               &
          &    time_dimid, sum_expoca_varid) )
     call nccheck( NF90_PUT_ATT(ncid, sum_expoca_varid, 'long_name',           &
          &    'Global total export production of carbonate') )
     call nccheck( NF90_PUT_ATT(ncid, sum_expoca_varid, 'units', 'kmol') )

     call nccheck( NF90_DEF_VAR(ncid, 'sum_exposi', NF90_DOUBLE,               &
          &    time_dimid, sum_exposi_varid) )
     call nccheck( NF90_PUT_ATT(ncid, sum_exposi_varid, 'long_name',           &
          &    'Global total export production of silicate') )
     call nccheck( NF90_PUT_ATT(ncid, sum_exposi_varid, 'units', 'kmol') )

     !--- End define mode.
     call nccheck( NF90_ENDDEF(ncid) )

  else
     !=== Open existing netCDF file
     write(io_stdo_bgc,*) 'Write BGC inventory to file : ',                    &
          &    trim(fname_inv(iogrp))
     call nccheck( NF90_OPEN(trim(fname_inv(iogrp)), NF90_WRITE, ncid) )
     !--- Inquire dimid
     call nccheck( NF90_INQ_DIMID(ncid, "time", time_dimid) )
     call nccheck( NF90_INQ_DIMID(ncid, "nocetra", nocetra_dimid) )
     !--- Inquire varid : time
     call nccheck( NF90_INQ_VARID(ncid, "time", time_varid) )
     !--- Inquire varid : ocean tracers
     call nccheck( NF90_INQ_VARID(ncid, "ztotvol", ztotvol_varid) )
     ! call nccheck( NF90_INQ_VARID(ncid, "zocetratot", zocetratot_varid) )
     ! call nccheck( NF90_INQ_VARID(ncid, "zocetratoc", zocetratoc_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_sco212", zt_sco212_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_sco212", zc_sco212_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_alkali", zt_alkali_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_alkali", zc_alkali_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_phosph", zt_phosph_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_phosph", zc_phosph_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_oxygen", zt_oxygen_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_oxygen", zc_oxygen_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_gasnit", zt_gasnit_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_gasnit", zc_gasnit_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_ano3", zt_ano3_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_ano3", zc_ano3_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_silica", zt_silica_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_silica", zc_silica_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_doc", zt_doc_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_doc", zc_doc_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_poc", zt_poc_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_poc", zc_poc_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_phyto", zt_phyto_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_phyto", zc_phyto_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_grazer", zt_grazer_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_grazer", zc_grazer_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_calciu", zt_calciu_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_calciu", zc_calciu_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_opal", zt_opal_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_opal", zc_opal_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_n2o", zt_n2o_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_n2o", zc_n2o_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_dms", zt_dms_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_dms", zc_dms_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_fdust", zt_fdust_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_fdust", zc_fdust_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_iron", zt_iron_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_iron", zc_iron_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_prefo2", zt_prefo2_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_prefo2", zc_prefo2_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_prefpo4", zt_prefpo4_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_prefpo4", zc_prefpo4_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_prefalk", zt_prefalk_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_prefalk", zc_prefalk_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_prefdic", zt_prefdic_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_prefdic", zc_prefdic_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_dicsat", zt_dicsat_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_dicsat", zc_dicsat_varid) )
#ifdef cisonew
     call nccheck( NF90_INQ_VARID(ncid, "zt_sco213", zt_sco213_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_sco213", zc_sco213_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_sco214", zt_sco214_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_sco214", zc_sco214_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_doc13", zt_doc13_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_doc13", zc_doc13_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_doc14", zt_doc14_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_doc14", zc_doc14_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_poc13", zt_poc13_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_poc13", zc_poc13_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_poc14", zt_poc14_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_poc14", zc_poc14_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_phyto13", zt_phyto13_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_phyto13", zc_phyto13_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_phyto14", zt_phyto14_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_phyto14", zc_phyto14_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_grazer13", zt_grazer13_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_grazer13", zc_grazer13_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_grazer14", zt_grazer14_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_grazer14", zc_grazer14_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_calciu13", zt_calciu13_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_calciu13", zc_calciu13_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_calciu14", zt_calciu14_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_calciu14", zc_calciu14_varid) )
#endif
#ifdef AGG
     call nccheck( NF90_INQ_VARID(ncid, "zt_snos", zt_snos_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_snos", zc_snos_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_adust", zt_adust_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_adust", zc_adust_varid) )
#endif
#ifdef CFC
     call nccheck( NF90_INQ_VARID(ncid, "zt_cfc11", zt_cfc11_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_cfc11", zc_cfc11_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_cfc12", zt_cfc12_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_cfc12", zc_cfc12_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_sf6", zt_sf6_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_sf6", zc_sf6_varid) )
#endif
#ifdef natDIC
     call nccheck( NF90_INQ_VARID(ncid, "zt_natsco212", zt_natsco212_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_natsco212", zc_natsco212_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_natalkali", zt_natalkali_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_natalkali", zc_natalkali_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zt_natcalciu", zt_natcalciu_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_natcalciu", zc_natcalciu_varid) )
#endif
#ifdef BROMO
     call nccheck( NF90_INQ_VARID(ncid, "zt_bromo", zt_bromo_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "zc_bromo", zc_bromo_varid) )
#endif
     !--- Inquire varid : sum of inventory
     call nccheck( NF90_INQ_VARID(ncid, "totcarb", totcarb_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "totphos", totphos_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "totsili", totsili_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "totnitr", totnitr_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "totoxyg", totoxyg_varid) )
     !--- Inquire varid : sediment fluxes
     call nccheck( NF90_INQ_VARID(ncid, "sum_zprorca", sum_zprorca_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "sum_zprcaca", sum_zprcaca_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "sum_zsilpro", sum_zsilpro_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "sum_expoor", sum_expoor_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "sum_expoca", sum_expoca_varid) )
     call nccheck( NF90_INQ_VARID(ncid, "sum_exposi", sum_exposi_varid) )
  end if

  !--- Increment record by 1, reset start and count arrays
  ncrec(iogrp) = ncrec(iogrp) + 1
  wrstart = (/ ncrec(iogrp) /)
  ! zocetra_wrstart = (/ 1, ncrec(iogrp) /)
  ! zocetra_count = (/ nocetra, 1 /)

  !=== Write output data to netCDF file
  !--- Write data : time
  datenum = time - time0
  call nccheck( NF90_PUT_VAR(ncid, time_varid, datenum, start = wrstart) )
  !--- Write data : ocean tracers
  call nccheck( NF90_PUT_VAR(ncid, ztotvol_varid, ztotvol, start = wrstart) )
  ! call nccheck( NF90_PUT_VAR(ncid, zocetratot_varid, zocetratot,               &
  !      &    start = zocetra_wrstart, count = zocetra_count) )
  ! call nccheck( NF90_PUT_VAR(ncid, zocetratoc_varid, zocetratoc,               &
  !      &    start = zocetra_wrstart, count = zocetra_count) )
  call nccheck( NF90_PUT_VAR(ncid, zt_sco212_varid,                            &
       &    zocetratot(isco212), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_sco212_varid,                            &
       &    zocetratoc(isco212), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_alkali_varid,                            &
       &    zocetratot(ialkali), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_alkali_varid,                            &
       &    zocetratoc(ialkali), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_phosph_varid,                            &
       &    zocetratot(iphosph), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_phosph_varid,                            &
       &    zocetratoc(iphosph), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_oxygen_varid,                            &
       &    zocetratot(ioxygen), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_oxygen_varid,                            &
       &    zocetratoc(ioxygen), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_gasnit_varid,                            &
       &    zocetratot(igasnit), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_gasnit_varid,                            &
       &    zocetratoc(igasnit), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_ano3_varid,                              &
       &    zocetratot(iano3), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_ano3_varid,                              &
       &    zocetratoc(iano3), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_silica_varid,                            &
       &    zocetratot(isilica), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_silica_varid,                            &
       &    zocetratoc(isilica), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_doc_varid,                               &
       &    zocetratot(idoc), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_doc_varid,                               &
       &    zocetratoc(idoc), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_poc_varid,                               &
       &    zocetratot(idet), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_poc_varid,                               &
       &    zocetratoc(idet), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_phyto_varid,                             &
       &    zocetratot(iphy), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_phyto_varid,                             &
       &    zocetratoc(iphy), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_grazer_varid,                            &
       &    zocetratot(izoo), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_grazer_varid,                            &
       &    zocetratoc(izoo), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_calciu_varid,                            &
       &    zocetratot(icalc), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_calciu_varid,                            &
       &    zocetratoc(icalc), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_opal_varid,                              &
       &    zocetratot(iopal), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_opal_varid,                              &
       &    zocetratoc(iopal), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_n2o_varid,                               &
       &    zocetratot(ian2o), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_n2o_varid,                               &
       &    zocetratoc(ian2o), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_dms_varid,                               &
       &    zocetratot(idms), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_dms_varid,                               &
       &    zocetratoc(idms), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_fdust_varid,                             &
       &    zocetratot(ifdust), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_fdust_varid,                             &
       &    zocetratoc(ifdust), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_iron_varid,                              &
       &    zocetratot(iiron), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_iron_varid,                              &
       &    zocetratoc(iiron), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_prefo2_varid,                            &
       &    zocetratot(iprefo2), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_prefo2_varid,                            &
       &    zocetratoc(iprefo2), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_prefpo4_varid,                           &
       &    zocetratot(iprefpo4), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_prefpo4_varid,                           &
       &    zocetratoc(iprefpo4), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_prefalk_varid,                           &
       &    zocetratot(iprefalk), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_prefalk_varid,                           &
       &    zocetratoc(iprefalk), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_prefdic_varid,                           &
       &    zocetratot(iprefdic), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_prefdic_varid,                           &
       &    zocetratoc(iprefdic), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_dicsat_varid,                            &
       &    zocetratot(idicsat), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_dicsat_varid,                            &
       &    zocetratoc(idicsat), start = wrstart) )
#ifdef cisonew
  call nccheck( NF90_PUT_VAR(ncid, zt_sco213_varid,                            &
       &    zocetratot(isco213), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_sco213_varid,                            &
       &    zocetratoc(isco213), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_sco214_varid,                            &
       &    zocetratot(isco214), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_sco214_varid,                            &
       &    zocetratoc(isco214), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_doc13_varid,                             &
       &    zocetratot(idoc13), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_doc13_varid,                             &
       &    zocetratoc(idoc13), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_doc14_varid,                             &
       &    zocetratot(idoc14), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_doc14_varid,                             &
       &    zocetratoc(idoc14), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_poc13_varid,                             &
       &    zocetratot(idet13), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_poc13_varid,                             &
       &    zocetratoc(idet13), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_poc14_varid,                             &
       &    zocetratot(idet14), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_poc14_varid,                             &
       &    zocetratoc(idet14), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_phyto13_varid,                           &
       &    zocetratot(iphy13), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_phyto13_varid,                           &
       &    zocetratoc(iphy13), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_phyto14_varid,                           &
       &    zocetratot(iphy14), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_phyto14_varid,                           &
       &    zocetratoc(iphy14), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_grazer13_varid,                          &
       &    zocetratot(izoo13), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_grazer13_varid,                          &
       &    zocetratoc(izoo13), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_grazer14_varid,                          &
       &    zocetratot(izoo14), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_grazer14_varid,                          &
       &    zocetratoc(izoo14), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_calciu13_varid,                          &
       &    zocetratot(icalc13), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_calciu13_varid,                          &
       &    zocetratoc(icalc13), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_calciu14_varid,                          &
       &    zocetratot(icalc14), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_calciu14_varid,                          &
       &    zocetratoc(icalc14), start = wrstart) )
#endif
#ifdef AGG
  call nccheck( NF90_PUT_VAR(ncid, zt_snos_varid,                              &
       &    zocetratot(inos), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_snos_varid,                              &
       &    zocetratoc(inos), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_adust_varid,                             &
       &    zocetratot(iadust), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_adust_varid,                             &
       &    zocetratoc(iadust), start = wrstart) )
#endif
#ifdef CFC
  call nccheck( NF90_PUT_VAR(ncid, zt_cfc11_varid,                             &
       &    zocetratot(icfc11), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_cfc11_varid,                             &
       &    zocetratoc(icfc11), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_cfc12_varid,                             &
       &    zocetratot(icfc12), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_cfc12_varid,                             &
       &    zocetratoc(icfc12), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_sf6_varid,                               &
       &    zocetratot(isf6), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_sf6_varid,                               &
       &    zocetratoc(isf6), start = wrstart) )
#endif
#ifdef natDIC
  call nccheck( NF90_PUT_VAR(ncid, zt_natsco212_varid,                         &
       &    zocetratot(inatsco212), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_natsco212_varid,                         &
       &    zocetratoc(inatsco212), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_natalkali_varid,                         &
       &    zocetratot(inatalkali), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_natalkali_varid,                         &
       &    zocetratoc(inatalkali), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zt_natcalciu_varid,                         &
       &    zocetratot(inatcalc), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_natcalciu_varid,                         &
       &    zocetratoc(inatcalc), start = wrstart) )
#endif
#ifdef BROMO
  call nccheck( NF90_PUT_VAR(ncid, zt_bromo_varid,                             &
       &    zocetratot(ibromo), start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, zc_bromo_varid,                             &
       &    zocetratoc(ibromo), start = wrstart) )
#endif
  !--- Write data : sum of inventory
  call nccheck( NF90_PUT_VAR(ncid, totcarb_varid, totalcarbon,                 &
       &    start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, totphos_varid, totalphos,                   &
       &    start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, totsili_varid, totalsil,                    &
       &    start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, totnitr_varid, totalnitr,                   &
       &    start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, totoxyg_varid, totaloxy,                    &
       &    start = wrstart) )
  !--- Write data : fluxes into and out of sediments
  call nccheck( NF90_PUT_VAR(ncid, sum_zprorca_varid, sum_zprorca,             &
       &    start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, sum_zprcaca_varid, sum_zprcaca,             &
       &    start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, sum_zsilpro_varid, sum_zsilpro,             &
       &    start = wrstart) )
  !--- Write data : global total export production
  call nccheck( NF90_PUT_VAR(ncid, sum_expoor_varid, sum_expoor,               &
       &    start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, sum_expoca_varid, sum_expoca,               &
       &    start = wrstart) )
  call nccheck( NF90_PUT_VAR(ncid, sum_exposi_varid, sum_exposi,               &
       &    start = wrstart) )

  !--- Close netCDF file
  call nccheck( NF90_CLOSE(ncid) )

end subroutine write_netcdf


subroutine nccheck(status)
  use netcdf
  implicit none

  integer, intent(in) :: status

  if (status /= nf90_noerr) then
     call xchalt('(inventory_bgc: Problem with netCDF)')
     stop '(inventory_bgc: Problem with netCDF)'
  end if
end subroutine nccheck


END SUBROUTINE INVENTORY_BGC
