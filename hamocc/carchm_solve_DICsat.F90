! Copyright (C) 2020  J. Tjiputra, J. Schwinger
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


SUBROUTINE carchm_solve_DICsat(saln,pco2,ta,sit,pt,                   &
                        Kh,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,        &
                        tc_sat,niter)
!**********************************************************************
!
!**** *CARCHM_SOLVE_DICsat* - .
!
!     J. Tjiputra,    *BCCR, Bergen*    25.01.17
!
!     Modified
!     --------
!
!     Purpose
!     -------
!     Solve DICsat from TALK and pCO2.
!
!     Method
!     -------
!
!
!**** Parameter list:
!     ---------------
!     *REAL*    *saln*    - salinity [psu].
!     *REAL*    *pco2*    - partial pressure of CO2 [ppm].
!     *REAL*    *ta*      - total alkalinity [eq/kg].
!     *REAL*    *sit*     - silicate concentration [mol/kg].
!     *REAL*    *pt*      - phosphate concentration [mol/kg].
!     *REAL*    *Kh*      - equilibrium constant K0  = [H2CO3]/pCO2.
!     *REAL*    *K1*      - equilibrium constant K1  = [H][HCO3]/[H2CO3].
!     *REAL*    *K2*      - equilibrium constant K2  = [H][CO3]/[HCO3].
!     *REAL*    *Kb*      - equilibrium constant Kb  = [H][BO2]/[HBO2].
!     *REAL*    *Kw*      - equilibrium constant Kw  = [H][OH].
!     *REAL*    *Ks1*     - equilibrium constant Ks1 = [H][SO4]/[HSO4].
!     *REAL*    *Kf*      - equilibrium constant Kf  = [H][F]/[HF].
!     *REAL*    *Ksi*     - equilibrium constant Ksi = [H][SiO(OH)3]/[Si(OH)4].
!     *REAL*    *K1p*     - equilibrium constant K1p = [H][H2PO4]/[H3PO4].
!     *REAL*    *K2p*     - equilibrium constant K2p = [H][HPO4]/[H2PO4].
!     *REAL*    *K3p*     - equilibrium constant K3p = [H][PO4]/[HPO4].
!     *REAL*    *tc_sat*  - saturated total DIC concentration [mol/kg].
!     *INTEGER* *niter*   - maximum number of iteration
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************
USE mo_chemcon, only: bor1,bor2,salchl

IMPLICIT NONE
REAL,    INTENT(IN)    :: saln,pco2,ta,sit,pt
REAL,    INTENT(IN)    :: Kh,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p
REAL,    INTENT(OUT)   :: tc_sat
INTEGER, INTENT(IN)    :: niter
  
! Parameters to set accuracy of iteration 
REAL,    PARAMETER     :: eps=5.e-5

! Local varibles
INTEGER                :: jit
REAL                   :: s,scl,borat,sti,ft
REAL                   :: hso4,hf,hsi,hpo4,ab,aw,ah2o,ah2,erel
REAL                   :: dic_h2co3,dic_hco3,dic_co3,ah1,ac



! Calculate concentrations for borate, sulfate, and fluoride; see Dickson, A.G.,
! Sabine, C.L. and Christian, J.R. (Eds.) 2007. Guide to best practices 
! for ocean CO2 measurements. PICES Special Publication 3, chapter 5 p. 10
s = MAX(25.,saln)
scl = s * salchl
borat = bor1 * scl * bor2      ! Uppstrom (1974)
sti = 0.14 * scl / 96.062      ! Morris & Riley (1966)
ft = 0.000067 * scl / 18.9984  ! Riley (1965)
ah1=1.e-8
dic_h2co3 = Kh * pco2 * 1e-6 

iflag: DO jit = 1,niter
   hso4 = sti / ( 1. + Ks1 / ( ah1 / ( 1. + sti / Ks1 ) ) )
   hf   = 1. / ( 1. + Kf / ah1 )
   hsi  = 1./ ( 1. + ah1 / Ksi )
   hpo4 = ( K1p * K2p * ( ah1 + 2. * K3p ) - ah1**3 ) /    & 
          ( ah1**3 + K1p * ah1**2 + K1p * K2p * ah1 + K1p * K2p * K3p )
   ab   = borat / ( 1. + ah1 / Kb )
   aw   = Kw / ah1 - ah1 / ( 1. + sti / Ks1 )
   ac   = ta + hso4 - sit * hsi - ab - aw + ft * hf - pt * hpo4
   ah2o = SQRT((K1*dic_h2co3)**2 + 4.*ac*2.*K1*k2*dic_h2co3) 
   ah2  = (K1*dic_h2co3 + ah2o)/(2.*ac)
   erel = ( ah2 - ah1 ) / ah2
   if (abs( erel ).ge.eps) then
      ah1 = ah2
   else
      exit iflag
   endif
ENDDO iflag

dic_hco3  = Kh * K1 *      pco2 * 1e-6 / ah1
dic_co3   = Kh * K1 * K2 * pco2 * 1e-6 / ah1**2
tc_sat    = dic_h2co3 + dic_hco3 + dic_co3 

END SUBROUTINE carchm_solve_DICsat




