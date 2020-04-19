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


SUBROUTINE CARCHM_KEQUI(temp,saln,prb,Kh,Khd,K1,K2,Kb,Kw,Ks1,Kf,Ksi,           &
                        K1p,K2p,K3p,Kspc,Kspa)
!*******************************************************************************
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
!     Calculate equilibrium constant for the carbonate system
!
!     Method
!     -------
!     J.Schwinger,      *Uni Research, Bergen*   2018-04-12
!     - added output Khd (CO2 solubility w.r.t. dry air) and
!       Kspa
!
!
!**** Parameter list:
!     ---------------
!
!     *REAL*    *temp*    - potential temperature [degr C].
!     *REAL*    *saln*    - salinity [psu].
!     *REAL*    *prb*     - pressure [bar].
!     *REAL*    *Kh*      - equilibrium constant Kh  =  [CO2]/pCO2, moist air.
!     *REAL*    *Khd*     - equilibrium constant Kh  =  [CO2]/pCO2, dry air.
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
!     *REAL*    *Kspc*    - equilibrium constant Kspc= [Ca2+]T [CO3]T.
!     *REAL*    *Kspa*    - equilibrium constant Kspa= [Ca2+]T [CO3]T.
!
!     Externals
!     ---------
!     none.
!
!*******************************************************************************
USE mo_chemcon, only: tzero,rgas,bor1,bor2,salchl,  &
                      ac1,ac2,ac3,ac4,bc1,bc2,bc3,  &
                      ad1,ad2,ad3,bd1,bd2,bd3,      &
                      a0,a1,a2,b0,b1,b2

IMPLICIT NONE
REAL,    INTENT(IN)    :: temp,saln,prb
REAL,    INTENT(OUT)   :: Kh,Khd,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,Kspc,Kspa

! Local varibles
INTEGER                :: js
REAL                   :: tk,tk100,invtk,dlogtk
REAL                   :: s,is,is2,sqrtis,s15,s2,sqrts,scl
REAL                   :: nKhwe74,deltav,deltak,zprb,zprb2
REAL                   :: lnkpok0(11)

s = MAX(25.,saln)
tk = temp + tzero
tk100 = tk/100.0
invtk = 1.0 / tk
dlogtk = log(tk)
is = 19.924 * s / ( 1000. - 1.005 * s )
is2 = is * is
sqrtis = SQRT(is)
s15    = s**1.5
s2     = s * s
sqrts  = SQRT(s)
scl    = s * salchl

! Kh = [CO2]/ p CO2
! Weiss (1974), refitted for moist air Weiss and Price (1980) [mol/kg/atm]
nKhwe74 = ac1+ac2/tk100+ac3*log(tk100)+ac4*tk100**2+s*(bc1+bc2*tk100+bc3*tk100**2)
Kh      = exp( nKhwe74 )
! Khd = [CO2]/ p CO2
! Weiss (1974) for dry air [mol/kg/atm]
nKhwe74 = ad1+ad2/tk100+ad3*log(tk100)+s*(bd1+bd2*tk100+bd3*tk100**2)
Khd     = exp( nKhwe74 )
! K1 = [H][HCO3]/[H2CO3]   ; K2 = [H][CO3]/[HCO3]
! Millero p.664 (1995) using Mehrbach et al. data on seawater scale
K1 = 10**( -1.0 * ( 3670.7 * invtk - 62.008 + 9.7944 * dlogtk - 0.0118 * s + 0.000116 * s2 ) )
K2 = 10**( -1.0 * ( 1394.7 * invtk + 4.777 - 0.0184 * s + 0.000118 * s2 ) )
! Kb = [H][BO2]/[HBO2] !
! Millero p.669 (1995) using DATA from Dickson (1990)
Kb = exp( ( -8966.90 - 2890.53  * sqrts - 77.942  * s + 1.728 * s15 - 0.0996 * s2 ) * invtk +    &
          ( 148.0248 + 137.1942 * sqrts + 1.62142 * s ) +                                        &
          ( -24.4344 - 25.085   * sqrts - 0.2474  * s ) * dlogtk + 0.053105 * sqrts * tk )
! K1p = [H][H2PO4]/[H3PO4] ; K2p = [H][HPO4]/[H2PO4] ; K3p = [H][PO4]/[HPO4]
! DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
K1p = exp( -4576.752 * invtk + 115.525 - 18.453 * dlogtk + ( -106.736 * invtk + 0.69171 ) *      &
           sqrts + ( -0.65643 * invtk - 0.01844 ) * s )
K2p = exp( -8814.715 * invtk + 172.0883 - 27.927 * dlogtk + ( -160.340 * invtk + 1.3566 ) *      &
           sqrts + ( 0.37335 * invtk - 0.05778 ) *s );
K3p = exp( -3070.75 * invtk - 18.141 + ( 17.27039 * invtk + 2.81197 ) * sqrts + ( -44.99486 *    &
           invtk - 0.09984 ) * s );
! Ksi = [H][SiO(OH)3]/[Si(OH)4]
! Millero p.671 (1995) using data from Yao and Millero (1995)
Ksi = exp( -8904.2 * invtk + 117.385 - 19.334 * dlogtk + ( -458.79 * invtk + 3.5913 ) * sqrtis   & 
       + ( 188.74 * invtk - 1.5998) * is + ( -12.1652 * invtk + 0.07871) * is2 +                 &
           log(1.0-0.001005*s))
! Kw = [H][OH] 
! Millero p.670 (1995) using composite data
Kw = exp( -13847.26 * invtk + 148.9652 - 23.6521 * dlogtk + ( 118.67 * invtk - 5.977 + 1.0495 *  &
          dlogtk ) * sqrts - 0.01615 * s)
! Ks = [H][SO4]/[HSO4]
! Dickson (1990, J. chem. Thermodynamics 22, 113)
Ks1 = exp( -4276.1 * invtk + 141.328 - 23.093 * dlogtk + ( -13856. * invtk + 324.57 - 47.986 *   &
           dlogtk ) * sqrtis + ( 35474. * invtk - 771.54 + 114.723 * dlogtk ) * is - 2698. *     &
           invtk * is**1.5 + 1776. * invtk * is2 + log(1.0 - 0.001005 * s ) )
! Kf = [H][F]/[HF]
! Dickson and Riley (1979) -- change pH scale to total
Kf = exp( 1590.2 * invtk - 12.641 + 1.525 * sqrtis + log( 1.0 - 0.001005 * s ) + log( 1.0 + (    &
          0.1400 / 96.062 ) * scl / Ks1 ) )
! Kspc (calcite)
! apparent solubility product of calcite : Kspc = [Ca2+]T [CO32-]T
! where $[]_T$ refers to the equilibrium total (free + complexed) ion concentration.
!          Mucci 1983 mol/kg-soln
Kspc = 10**( -171.9065 - 0.077993 * tk + 2839.319 / tk + 71.595 * log10( tk ) + ( - 0.77712 +    &
             0.0028426 * tk + 178.34 / tk ) * sqrts - 0.07711 * s + 0.0041249 * s15 );
! Kspa (aragonite)
! apparent solubility product of aragonite : Kspa = [Ca2+]T [CO32-]T
! where $[]_T$ refers to the equilibrium total (free + complexed) ion concentration.
!          Mucci 1983 mol/kg-soln
Kspa = 10**( -171.945 - 0.077993 * tk + 2903.293 / tk  + 71.595 * log10( tk ) + ( -0.068393 +    &
             0.0017276 * tk + 88.135 / tk ) * sqrts - 0.10018 * s + 0.0059415 * s15 );


!---------------------- Pressure effect on Ks (Millero, 95) --------------------
! index: K1 1, K2 2, Kb 3, Kw 4, Ks 5, Kf 6, Kspc 7, Kspa 8, K1p 9, K2p 10, K3p 11
DO js = 1,11
   deltav      = a0(js) + a1(js) * temp + a2(js) * temp * temp
   deltak      = b0(js) + b1(js) * temp + b2(js) * temp * temp
   zprb        = prb / ( rgas * tk )
   zprb2       = prb * zprb
   lnkpok0(js) = - ( deltav * zprb + 0.5 * deltak * zprb2 )
ENDDO

K1   = K1   * exp( lnkpok0(1)  )
K2   = K2   * exp( lnkpok0(2)  )
Kb   = Kb   * exp( lnkpok0(3)  )
Kw   = Kw   * exp( lnkpok0(4)  )
Ks1  = Ks1  * exp( lnkpok0(5)  )
Kf   = Kf   * exp( lnkpok0(6)  )
Kspc = Kspc * exp( lnkpok0(7)  )
Kspa = Kspa * exp( lnkpok0(8)  )
K1p  = K1p  * exp( lnkpok0(9)  )
K2p  = K2p  * exp( lnkpok0(10) )
K3p  = K3p  * exp( lnkpok0(11) )


END SUBROUTINE CARCHM_KEQUI
