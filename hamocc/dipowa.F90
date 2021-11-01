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


subroutine dipowa(kpie,kpje,kpke,omask)
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
!     constant diffusion coefficient : 1.e-9 set in BODENSED.
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

  use mo_carbch
  use mo_sedmnt
  use mo_biomod
  use mo_param1_bgc
  use mo_control_bgc
  use mo_vgrid, only: kbo,bolay

  implicit none

  integer, intent(in) :: kpie, kpje, kpke
  real, dimension(kpie,kpje), intent(in) :: omask

  ! Local variables
  integer :: i,j,k,l,iv
  integer :: iv_oc                     ! index of ocetra in powtra loop

  real :: sedb1(kpie,0:ks,npowtra)     ! ????
  real :: zcoefsu(0:ks),zcoeflo(0:ks)  ! diffusion coefficients (upper/lower)
  real :: tredsy(kpie,0:kpke,3)        ! redsy for 'reduced system'?
  real :: aprior                       ! start value of oceanic tracer in bottom layer

!ik accelerated sediment
!ik needed for boundary layer ventilation in fast sediment routine
  real :: bolven(kpie)                 ! bottom layer ventilation rate

  zcoefsu(0) = 0.0
  do k = 1,ks
     ! sediment diffusion coefficient * 1/dz * fraction of pore water at half depths
     zcoefsu(k  ) = -sedict * seddzi(k) * porwah(k)
     zcoeflo(k-1) = -sedict * seddzi(k) * porwah(k)    ! why the same ?
  enddo
  zcoeflo(ks) = 0.0                    ! diffusion coefficient for bottom sediment layer

!$OMP PARALLEL DO                            &
!$OMP&PRIVATE(i,k,iv,l,bolven,tredsy,sedb1,aprior,iv_oc)
  j_loop: do j=1,kpje

! calculate bottom ventilation rate for scaling of sediment-water exchange
  do i = 1,kpie
     bolven(i) = 1.
  enddo

  k = 0
  do i = 1,kpie
     tredsy(i,k,1) = zcoefsu(k)
     tredsy(i,k,3) = zcoeflo(k)
     tredsy(i,k,2) = bolven(i)*bolay(i,j) - tredsy(i,k,1) - tredsy(i,k,3)
     !                            dz(kbo) - diff upper    - diff lower
  enddo

  k = 0
  do iv = 1,npowtra      ! loop over pore water tracers
     iv_oc = iv
#ifdef cisonew
     if (iv == ipowc13) iv_oc = isco213
     if (iv == ipowc14) iv_oc = isco214
#endif
     do i = 1,kpie
        sedb1(i,k,iv) = 0.
        if (omask(i,j) > 0.5) then
           sedb1(i,k,iv) = ocetra(i,j,kbo(i,j),iv_oc) * bolay(i,j)*bolven(i)
           !               tracer_concentration(kbo)  * dz(kbo)
        endif
     enddo
  enddo

  do k = 1,ks
     do i = 1,kpie
        tredsy(i,k,1) = zcoefsu(k)
        tredsy(i,k,3) = zcoeflo(k)
        tredsy(i,k,2) = seddw(k)*porwat(k) -tredsy(i,k,1) -tredsy(i,k,3)
     enddo
  enddo

  do iv = 1,npowtra
     do k = 1,ks
        do i = 1,kpie
           ! tracer_concentration(k[1:ks]) * porewater fraction(k) * dz(k)
           sedb1(i,k,iv) = powtra(i,j,k,iv) * porwat(k) * seddw(k)
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

!  call maschk(kpie,kpje,kpke,23)
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
!  call maschk(kpie,kpje,kpke,24)

! sediment ocean interface
!
! CAUTION - the following assumes same indecees for ocetra and powtra
!           test npowa_base 071106
!           check mo_param1_bgc.f90 for consistency
  do iv = 1, npowtra
     iv_oc = iv
#ifdef cisonew
     if (iv == ipowc13) iv_oc=isco213
     if (iv == ipowc14) iv_oc=isco214
#endif
     do i = 1,kpie
        l = 0
        if (omask(i,j) > 0.5) then
           aprior = ocetra(i,j,kbo(i,j),iv_oc)
           ocetra(i,j,kbo(i,j),iv_oc) =                                        &
                &  ( sedb1(i,l,iv) - tredsy(i,l,3) * powtra(i,j,l+1,iv) )      &
                &  / tredsy(i,l,2)

           ! used in inventory_bgc/maschk (diagnostics)
           sedfluxo(i,j,iv) = sedfluxo(i,j,iv)                                 &
                &  + ocetra(i,j,kbo(i,j),iv) - aprior
#ifdef natDIC
           if (iv==isco212) ocetra(i,j,kbo(i,j),inatsco212) =                  &
                &  ocetra(i,j,kbo(i,j),inatsco212) +                           &
                &  ocetra(i,j,kbo(i,j),iv) - aprior
           if (iv==ialkali) ocetra(i,j,kbo(i,j),inatalkali) =                  &
                &  ocetra(i,j,kbo(i,j),inatalkali) +                           &
                &  ocetra(i,j,kbo(i,j),iv) - aprior
#endif
        endif
     enddo
  enddo
!  call maschk(kpie,kpje,kpke,25)

  enddo j_loop

end subroutine dipowa
