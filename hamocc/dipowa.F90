! Copyright (C) 2001  Ernst Maier-Reimer, S. Legutke
! Copyright (C) 2020  J. Schwinger
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


SUBROUTINE dipowa(kpie,kpje,kpke,omask,lspin)
!**********************************************************************
!
!**** *DIPOWA* - 'diffusion of pore water'
!      vertical diffusion of sediment pore water tracers
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!     - all npowtra-1 properties are diffused in 1 go.
!     js: not mass conserving check c13/powtra/ocetra
!
!     Purpose
!     -------
!     calculate vertical diffusion of sediment pore water properties
!     and diffusive flux through the ocean/sediment interface.
!     integration.
!
!     Method
!     -------
!     implicit formulation;
!     constant diffusion coefficient : 1.e-9 set in ini_sedmnt in mo_sedmnt
!     diffusion coefficient : zcoefsu/zcoeflo for upper/lower
!     sediment layer boundary.
!
!**   Interface.
!     ----------
!
!     *CALL*       *DIPOWA*
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************

  use mo_carbch,     only: ocetra, sedfluxo
  use mo_sedmnt,     only: powtra,porwat,porwah,seddw,zcoefsu,zcoeflo
  use mo_param1_bgc, only: ks,npowtra,map_por2octra
  use mo_vgrid,      only: kbo,bolay
  ! cisonew
  use mo_param1_bgc, only: ipowc13,ipowc14,isco213,isco214
  ! natDIC
  use mo_param1_bgc, only: ialkali,inatalkali,inatsco212,isco212
  use mo_control_bgc, only: use_natDIC

  implicit none

  integer, intent(in) :: kpie, kpje, kpke
  real,    intent(in) :: omask(kpie,kpje)
  logical, intent(in) :: lspin

  ! Local variables
  integer :: i,j,k,l,iv
  integer :: iv_oc                     ! index of ocetra in powtra loop

  real :: sedb1(kpie,0:ks,npowtra)     ! ????
  real :: tredsy(kpie,0:kpke,3)        ! redsy for 'reduced system'?
  real :: aprior                       ! start value of oceanic tracer in bottom layer


!$OMP PARALLEL DO                            &
!$OMP&PRIVATE(i,k,iv,l,tredsy,sedb1,aprior,iv_oc)
  j_loop: do j=1,kpje

  k = 0
  do i = 1,kpie
     tredsy(i,k,1) = zcoefsu(i,j,k)
     tredsy(i,k,3) = zcoeflo(i,j,k)
     tredsy(i,k,2) = bolay(i,j) - tredsy(i,k,1) - tredsy(i,k,3)
     !                  dz(kbo) - diff upper    - diff lower
  enddo

  k = 0
  do iv = 1,npowtra      ! loop over pore water tracers
     iv_oc = map_por2octra(iv)
     do i = 1,kpie
        sedb1(i,k,iv) = 0.
        if (omask(i,j) > 0.5) then
           sedb1(i,k,iv) = ocetra(i,j,kbo(i,j),iv_oc) * bolay(i,j)
           !               tracer_concentration(kbo)  * dz(kbo)
        endif
     enddo
  enddo

  do k = 1,ks
     do i = 1,kpie
        tredsy(i,k,1) = zcoefsu(i,j,k)
        tredsy(i,k,3) = zcoeflo(i,j,k)
        tredsy(i,k,2) = seddw(k)*porwat(i,j,k) -tredsy(i,k,1) -tredsy(i,k,3)
     enddo
  enddo

  do iv = 1,npowtra
     do k = 1,ks
        do i = 1,kpie
           ! tracer_concentration(k[1:ks]) * porewater fraction(k) * dz(k)
           sedb1(i,k,iv) = powtra(i,j,k,iv) * porwat(i,j,k) * seddw(k)
        enddo
     enddo
  enddo

  do k = 1,ks
     do i = 1,kpie
        if (omask(i,j) > 0.5) then
           ! this overwrites tredsy(k=0) for k=1
           tredsy(i,k-1,1) = tredsy(i,k,1) / tredsy(i,k-1,2)
           !                 diff upper    / conc (k-1)
           tredsy(i,k,2)   = tredsy(i,k,2)                                     &
                &  - tredsy(i,k-1,3) * tredsy(i,k,1) / tredsy(i,k-1,2)
           ! concentration - diff lower * diff upper / conc(k-1)
        endif
     enddo
  enddo

! diffusion from above
  do iv = 1,npowtra
     do k = 1,ks
        do i = 1,kpie
           sedb1(i,k,iv) = sedb1(i,k,iv) - tredsy(i,k-1,1) * sedb1(i,k-1,iv)
        enddo
     enddo
  enddo

! sediment bottom layer
  k = ks
  do iv = 1,npowtra
     do i = 1,kpie
        if (omask(i,j) > 0.5) then
           powtra(i,j,k,iv) = sedb1(i,k,iv) / tredsy(i,k,2)
        endif
     enddo
  enddo

! sediment column
  do iv = 1,npowtra
     do k = 1,ks-1
        l = ks-k
        do i = 1,kpie
           if (omask(i,j) > 0.5) then
              powtra(i,j,l,iv) = ( sedb1(i,l,iv)                               &
                   &  - tredsy(i,l,3) * powtra(i,j,l+1,iv) ) / tredsy(i,l,2)
           endif
        enddo
     enddo
  enddo

  if(.not. lspin) then
! sediment ocean interface
  do iv = 1, npowtra
     iv_oc = map_por2octra(iv)
     do i = 1,kpie
        l = 0
        if (omask(i,j) > 0.5) then
           aprior = ocetra(i,j,kbo(i,j),iv_oc)
           ocetra(i,j,kbo(i,j),iv_oc) =                                        &
                &  ( sedb1(i,l,iv) - tredsy(i,l,3) * powtra(i,j,l+1,iv) )      &
                &  / tredsy(i,l,2)

           ! diffusive fluxes (positive downward)
           sedfluxo(i,j,iv) = sedfluxo(i,j,iv)                                 &
                &  -(ocetra(i,j,kbo(i,j),iv_oc) - aprior)* bolay(i,j)
           if (use_natDIC) then
              ! workaround as long as natDIC is not implemented throughout the sediment MODULE 
              if (iv_oc==isco212) ocetra(i,j,kbo(i,j),inatsco212) =            &
                   &  ocetra(i,j,kbo(i,j),inatsco212) +                        &
                   &  ocetra(i,j,kbo(i,j),isco212) - aprior
              if (iv_oc==ialkali) ocetra(i,j,kbo(i,j),inatalkali) =            &
                   &  ocetra(i,j,kbo(i,j),inatalkali) +                        &
                   &  ocetra(i,j,kbo(i,j),ialkali) - aprior
           endif
        endif
     enddo
  enddo

  endif ! .not. lspin
  
  enddo j_loop

END SUBROUTINE dipowa
