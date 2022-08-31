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


subroutine bodensed(kpie,kpje,kpke,pddpo)
!**********************************************************************
!
!**** *BODENSED* - .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!
!     Purpose
!     -------
!     set up of sediment layer.
!
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!
!**********************************************************************
  
  use mo_sedmnt,      only: calcwei,calfa,clafa,claydens,calcdens,opaldens,opalwei,oplfa,orgdens,orgfa,seddzi,porwat,porwah,       &
                          & porsol,dzs,seddw,sedict,solfu,orgwei,zcoefsu,zcoeflo,disso_sil,silsat,disso_poc,sed_denit,disso_caco3
  use mo_control_bgc, only: dtbgc,io_stdo_bgc
  use mo_param1_bgc,  only: ks
  use mod_xc,         only: mnproc

  implicit none

  integer, intent(in) :: kpie,kpje,kpke
  real,    intent(in) :: pddpo(kpie,kpje,kpke)

  ! Local variables
  integer             :: i,j,k
  real                :: sumsed

  dzs(1) = 0.001
  dzs(2) = 0.003
  dzs(3) = 0.005
  dzs(4) = 0.007
  dzs(5) = 0.009
  dzs(6) = 0.011
  dzs(7) = 0.013
  dzs(8) = 0.015
  dzs(9) = 0.017
  dzs(10) = 0.019
  dzs(11) = 0.021
  dzs(12) = 0.023
  dzs(13) = 0.025

  if (mnproc == 1) then
     write(io_stdo_bgc,*)  ' '
     write(io_stdo_bgc,*)  'Sediment layer thickness [m] : '
     write(io_stdo_bgc,'(5F9.3)') dzs
     write(io_stdo_bgc,*)  ' '
  endif

  ! this initialization can be done later via reading a porosity map
  porwat(:,:,1) = 0.85
  porwat(:,:,2) = 0.83
  porwat(:,:,3) = 0.8
  porwat(:,:,4) = 0.79
  porwat(:,:,5) = 0.77
  porwat(:,:,6) = 0.75
  porwat(:,:,7) = 0.73
  porwat(:,:,8) = 0.7
  porwat(:,:,9) = 0.68
  porwat(:,:,10) = 0.66
  porwat(:,:,11) = 0.64
  porwat(:,:,12) = 0.62

  if (mnproc == 1) then
     write(io_stdo_bgc,*)  'Pore water in sediment initialized'
  endif

  seddzi(1) = 500.
  do k = 1, ks
     seddzi(k+1) = 1. / dzs(k+1)
     seddw(k) = 0.5 * (dzs(k) + dzs(k+1))
     do i = 1, kpie
     do j = 1, kpje
        porsol(i,j,k) = 1. - porwat(i,j,k)
        if(k >= 2) porwah(i,j,k) = 0.5 * (porwat(i,j,k) + porwat(i,j,k-1))
        if(k == 1) porwah(i,j,k) = 0.5 * (1. + porwat(i,j,1))
     enddo
     enddo
  enddo
  
  sedict = 1.e-9 * dtbgc ! Moecular diffusion coefficient
  ! Dissolution rate constant of opal (disso) [1/(kmol Si(OH)4/m3)*1/sec]
  ! THIS NEEDS TO BE CHANGED TO disso=3.e-8! THIS IS ONLY KEPT FOR THE MOMENT
  ! FOR BACKWARDS COMPATIBILITY
  !disso=3.e-8  ! (2011-01-04) EMR
  !disso=1.e-6 ! test vom 03.03.04 half live sil ca. 20.000 yr 
  disso_sil = 1.e-6*dtbgc
  ! Silicate saturation concentration is 1 mol/m3
  silsat    = 0.001

  ! Degradation rate constant of POP (disso) [1/(kmol O2/m3)*1/sec]
  disso_poc = 0.01 / 86400. * dtbgc  !  disso=3.e-5 was quite high

  ! Denitrification rate constant of POP (disso) [1/sec]
  sed_denit =  0.01/86400. * dtbgc !ik      denit = 1.e-6*dtbgc

  ! Dissolution rate constant of CaCO3 (disso) [1/(kmol CO3--/m3)*1/sec]
  disso_caco3 = 1.e-7 * dtbgc

! ******************************************************************
! densities etc. for SEDIMENT SHIFTING

! define weight of calcium carbonate, opal, and poc [kg/kmol]
  calcwei = 100.           ! 40+12+3*16 kg/kmol C
  opalwei = 60.            ! 28 + 2*16  kg/kmol Si
  orgwei  = 30.            ! from 12 kg/kmol * 2.5 POC[kg]/DW[kg]
                           ! after Alldredge, 1998:
                           ! POC(g)/DW(g) = 0.4 of diatom marine snow, size 1mm3

! define densities of opal, caco3, poc [kg/m3]
  calcdens = 2600.
  opaldens = 2200.
  orgdens  = 1000.
  claydens = 2600.         !quartz

! define volumes occupied by solid constituents [m3/kmol]
  calfa = calcwei / calcdens
  oplfa = opalwei / opaldens
  orgfa = orgwei / orgdens
  clafa = 1. / claydens    !clay is calculated in kg/m3

! determine total solid sediment volume
  solfu = 0.
  do i = 1, kpie
  do j = 1, kpje
  do k = 1, ks
     solfu(i,j) = solfu(i,j) + seddw(k) * porsol(i,j,k)
  enddo
  enddo
  enddo

! Initialize porosity-dependent diffusion coefficients of sediment
  zcoefsu(:,:,0) = 0.0
  do i = 1, kpie
  do j = 1, kpje
  do k = 1,ks
     ! sediment diffusion coefficient * 1/dz * fraction of pore water at half depths
     zcoefsu(i,j,k  ) = -sedict * seddzi(k) * porwah(i,j,k)
     zcoeflo(i,j,k-1) = -sedict * seddzi(k) * porwah(i,j,k)    ! why the same ?
  enddo
  enddo
  enddo
  zcoeflo(:,:,ks) = 0.0                    ! diffusion coefficient for bottom sediment layer

  
end subroutine bodensed
