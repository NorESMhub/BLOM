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


SUBROUTINE CARCHM_SOLVE(saln,tc,ta,sit,pt,                            &
                        K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,           &
                        ah1,ac,niter)
!**********************************************************************
!
!**** *CARCHM_SOLVE* - .
!
!     J. Schwinger,    *BCCR, Bergen*    09.02.16
!
!     Modified
!     --------
!
!     Purpose
!     -------
!     Solve carbon chemistry.
!
!     Method
!     -------
!
!
!**** Parameter list:
!     ---------------
!     *REAL*    *saln*    - salinity [psu].
!     *REAL*    *tc*      - total DIC concentraion [mol/kg].
!     *REAL*    *ta*      - total alkalinity [eq/kg].
!     *REAL*    *sit*     - silicate concentration [mol/kg].
!     *REAL*    *pt*      - phosphate concentration [mol/kg].
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
!     *REAL*    *ah1*     - hydrogen ion concentration.
!     *REAL*    *ac*      - carbonate alkalinity.
!     *INTEGER* *niter*   - maximum number of iteration
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************
USE mo_chemcon, only: bor1,bor2,salchl

IMPLICIT NONE
REAL,    INTENT(IN)    :: saln,tc,ta,sit,pt
REAL,    INTENT(IN)    :: K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p
REAL,    INTENT(INOUT) :: ah1
REAL,    INTENT(OUT)   :: ac
INTEGER, INTENT(IN)    :: niter
  
! Parameters to set accuracy of iteration 
REAL,    PARAMETER     :: eps=5.e-5

! Local varibles
INTEGER                :: jit
REAL                   :: s,scl,borat,sti,ft
REAL                   :: hso4,hf,hsi,hpo4,ab,aw,ah2o,ah2,erel



! Calculate concentrations for borate, sulfate, and fluoride; see Dickson, A.G.,
! Sabine, C.L. and Christian, J.R. (Eds.) 2007. Guide to best practices 
! for ocean CO2 measurements. PICES Special Publication 3, chapter 5 p. 10
s = MAX(25.,saln)
scl = s * salchl
borat = bor1 * scl * bor2      ! Uppstrom (1974)
sti = 0.14 * scl / 96.062      ! Morris & Riley (1966)
ft = 0.000067 * scl / 18.9984  ! Riley (1965)


iflag: DO jit = 1,niter
   hso4 = sti / ( 1. + Ks1 / ( ah1 / ( 1. + sti / Ks1 ) ) )
   hf   = 1. / ( 1. + Kf / ah1 )
   hsi  = 1./ ( 1. + ah1 / Ksi )
   hpo4 = ( K1p * K2p * ( ah1 + 2. * K3p ) - ah1**3 ) /    & 
          ( ah1**3 + K1p * ah1**2 + K1p * K2p * ah1 + K1p * K2p * K3p )
   ab   = borat / ( 1. + ah1 / Kb )
   aw   = Kw / ah1 - ah1 / ( 1. + sti / Ks1 )
   ac   = ta + hso4 - sit * hsi - ab - aw + ft * hf - pt * hpo4
   ah2o = SQRT( ( tc - ac )**2 + 4. * ( ac * K2 / K1 ) * ( 2. * tc - ac ) )
   ah2  = 0.5 * K1 / ac *( ( tc - ac ) + ah2o )
   erel = ( ah2 - ah1 ) / ah2
   if (abs( erel ).ge.eps) then
      ah1 = ah2
   else
      exit iflag
   endif
ENDDO iflag

END SUBROUTINE CARCHM_SOLVE

