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


module mo_chemcon

  !*************************************************************************************************
  ! Parameter definitions for chemical formulas (previously defined in subroutine chemcon)
  !
  ! J. Schwinger,    *UiB-GfI, Bergen*    2013-08-21
  !
  ! Modified
  ! J.Schwinger,      *Uni Research, Bergen*   2018-04-12
  ! - added constants for Kh CO2 w.r.t. dry air (Weiss, 1974)
  !*************************************************************************************************

  use mo_kind, only: rp

  implicit none
  public

  !     -----------------------------------------------------------------
  !*    BORON CONCENTRATION IN SEA WATER IN G/KG PER O/OO CL
  !     (RILEY AND SKIRROW, 1965, P.250)
  !
  real, parameter :: BOR1=0.000232_rp

  !     -----------------------------------------------------------------
  !*    INVERSE OF ATOMIC WEIGHT OF BORON [G**-1]
  !     (USED TO CONVERT SPECIFIC TOTAL BORAT INTO CONCENTRATIONS)
  !
  real, parameter :: BOR2=1._rp/10.811_rp


  !     -----------------------------------------------------------------
  !*    CONVERSION FACTOR SALINITY -> CHLORINITY
  !     (AFTER WOOSTER ET AL., 1969)
  !
  real, parameter :: SALCHL=1._rp/1.80655_rp
  real, parameter :: rrrcl=salchl*1.025_rp*bor1*bor2

  !     -----------------------------------------------------------------
  !*    ZERO DEG CENTIGRADE AT KELVIN SCALE
  !
  real, parameter :: tzero=273.15_rp

  !     -----------------------------------------------------------------
  !*    SET MEAN TOTAL [CA++] IN SEAWATER (MOLES/KG) (SEE BROECKER
  !     A. PENG, 1982, P. 26; [CA++](MOLES/KG)=1.028E-2*(S/35.); Value
  !     taken from Sarmiento and Gruber, 2006, p. 365
  !
  real, parameter :: CALCON=0.01028_rp

  !     -----------------------------------------------------------------
  !*    INVERS OF NORMAL MOLAL VOLUME OF AN IDEAL GAS [mol/ml] at 0C
  !
  real, parameter :: OXYCO=1._rp/22414.4_rp

  !     -----------------------------------------------------------------
  !*    VOLUMETRIC SOLUBILITY CONSTANTS FOR O2 IN ML/L from moist air at
  !     one atm total pressure. Table 2 in WEISS, R. F. (1970) THE
  !     SOLUBILITY OF NITROGEN OXYGEN AND ARGON IN WATER AND SEAWATER.
  !     DEEP-SEA RESEARCH, VOL. 17, 721-735.
  !
  real, parameter :: OX0=-173.4292_rp
  real, parameter :: OX1=249.6339_rp
  real, parameter :: OX2=143.3483_rp
  real, parameter :: OX3=-21.8492_rp
  real, parameter :: OX4=-0.033096_rp
  real, parameter :: OX5=0.014259_rp
  real, parameter :: OX6=-0.0017_rp

  !     -----------------------------------------------------------------
  !*    VOLUMETRIC SOLUBILITY CONSTANTS FOR N2 IN ML/L from moist air at
  !     one atm total pressure. Table 2 in WEISS, R. F. (1970) THE
  !     SOLUBILITY OF NITROGEN OXYGEN AND ARGON IN WATER AND SEAWATER.
  !     DEEP-SEA RESEARCH, VOL. 17, 721-735.
  !
  real, parameter :: AN0=-172.4965_rp
  real, parameter :: AN1=248.4262_rp
  real, parameter :: AN2=143.0738_rp
  real, parameter :: AN3=-21.7120_rp
  real, parameter :: AN4=-0.049781_rp
  real, parameter :: AN5=0.025018_rp
  real, parameter :: AN6=-0.0034861_rp

  !     -----------------------------------------------------------------
  !     Constants for CO2 solubility in mol/kg/atm from moist
  !     air at one atm total pressure. Table 6 in WEISS, R.F.,
  !     NITROUS OXIDE SOLUBILITY IN WATER AND SEAWATER,
  !     Marine Chemistry, 8, 347-359, 1980
  !
  real, parameter :: ac1= -162.8301_rp
  real, parameter :: ac2=  218.2968_rp
  real, parameter :: ac3=   90.9241_rp
  real, parameter :: ac4=   -1.47696_rp
  real, parameter :: bc1=    0.025695_rp
  real, parameter :: bc2=   -0.025225_rp
  real, parameter :: bc3=    0.0049867_rp

  !     -----------------------------------------------------------------
  !     Constants for CO2 solubility in mol/kg/atm for dry
  !     air at one atm total pressure. Table 1 in WEISS, R.F.,
  !     CARBON DIOXIDE IN WATER AND SEAWATER: THE SOLUBILITY OF
  !     A NON - IDEAL GAS, Marine Chemistry, 2, 203-215, 1974
  !
  real, parameter :: ad1= -60.2409_rp
  real, parameter :: ad2=  93.4517_rp
  real, parameter :: ad3=  23.3585_rp
  real, parameter :: bd1=   0.023517_rp
  real, parameter :: bd2=  -0.023656_rp
  real, parameter :: bd3=   0.0047036_rp

  !     -----------------------------------------------------------------
  !     Constants for laughing gas solubility in mol/l/atm from moist
  !     air at one atm total pressure. Table 2 in WEISS, R.F.,
  !     NITROUS OXIDE SOLUBILITY IN WATER AND SEAWATER,
  !     Marine Chemistry, 8, 347-359, 1980
  !
  real, parameter :: al1= -165.8806_rp
  real, parameter :: al2=  222.8743_rp
  real, parameter :: al3=   92.0792_rp
  real, parameter :: al4=   -1.48425_rp
  real, parameter :: bl1=   -0.056235_rp
  real, parameter :: bl2=    0.031619_rp
  real, parameter :: bl3=   -0.0048472_rp

  !     -----------------------------------------------------------------
  !     NH3/NH4
  !     Tsilingiris 2008
  !     moist air dynamic viscosity parameters
  real, parameter :: SV0_air =  1.715747771e-5_rp
  real, parameter :: SV1_air =  4.722402075e-8_rp
  real, parameter :: SV2_air = -3.663027156e-10_rp
  real, parameter :: SV3_air =  1.873236686e-12_rp
  real, parameter :: SV4_air = -8.050218737e-14_rp

  ! moist air density parameters
  real, parameter :: SD0_air =  1.293393662_rp
  real, parameter :: SD1_air = -5.538444326e-3_rp
  real, parameter :: SD2_air =  3.860201577e-5_rp
  real, parameter :: SD3_air = -5.2536065e-7_rp

  ! diffusion of NH3 in water and air
  real, parameter :: Va_air  = 20.1_rp  ! Johnson 2010
  real, parameter :: Ma_air  = 28.97_rp ! Johnson 2010
  real, parameter :: Mb_nh3  = 17.03_rp ! Johnson 2010, Tang 2014
  real, parameter :: Vb_nh3  = 20.7_rp  ! Johnson 2010
  real, parameter :: M_nh3   =  (1._rp/Ma_air + 1._rp/Mb_nh3)**0.5_rp &
                                / (Va_air**(1._rp/3._rp)+Vb_nh3**(1._rp/3._rp))**2
  real, parameter :: kappa   = 0.4_rp   ! von Karman constant

  !     -----------------------------------------------------------------
  !     Molecular weights
  real, parameter :: mw_nitrogen = 14.00674_rp  ! [g/mol N]   nitrogen mol-weight as defined by CAM
  real, parameter :: mw_nh3      = 17.028940_rp ! [g/mol NH3] ammonia mol-weight as defined by CAM
  real, parameter :: mw_n2o      = 44.012880_rp ! [g/mol N2O] nitrous oxide mol-weight as defined by CAM
  real, parameter :: mw_fe       = 55.85_rp     ! [g/mol Fe]  iron mol-weight

  !     -----------------------------------------------------------------
  !     Constants needed for pressure correction of equilibrium constants
  !     F. Millero, Thermodynamics of the carbon dioxide system in the oceans,
  !     Geochimica et Cosmochimica Acta, Vol. 59, No. 4, pp. 661-677, 1995
  real, dimension(11) :: a0, a1, a2, b0, b1, b2
  data a0 /-25.5_rp, -15.82_rp, -29.48_rp, -25.60_rp, -18.03_rp, -9.78_rp, -48.76_rp, &
           -46._rp, -14.51_rp, -23.12_rp, -26.57_rp/
  data a1 /0.1271_rp, -0.0219_rp, 0.1622_rp, 0.2324_rp, 0.0466_rp, -0.0090_rp,     &
           0.5304_rp, 0.5304_rp, 0.1211_rp, 0.1758_rp, 0.2020_rp/
  data a2 /0.0_rp, 0.0_rp, 2.608e-3_rp, -3.6246e-3_rp, 0.316e-3_rp,             &
          -0.942e-3_rp, 0.0_rp, 0.0_rp, -0.321e-3_rp, -2.647e-3_rp, -3.042e-3_rp/
  data b0 /-3.08e-3_rp, 1.13e-3_rp, -2.84e-3_rp, -5.13e-3_rp, -4.53e-3_rp,      &
           -3.91e-3_rp, -11.76e-3_rp, -11.76e-3_rp, -2.67e-3_rp, -5.15e-3_rp,   &
           -4.08e-3_rp/
  data b1 /0.0877e-3_rp, -0.1475e-3_rp, 0.0_rp, 0.0794e-3_rp, 0.09e-3_rp,       &
           0.054e-3_rp, 0.3692e-3_rp, 0.3692e-3_rp, 0.0427e-3_rp,            &
           0.09e-3_rp, 0.0714e-3_rp/
  data b2 /0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp/

  !     -----------------------------------------------------------------
  !     Gas constant, value as used by Millero (1995)
  real, parameter :: rgas = 83.131_rp

end module mo_chemcon
