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

  implicit none
  public

  !     -----------------------------------------------------------------
  !*    BORON CONCENTRATION IN SEA WATER IN G/KG PER O/OO CL
  !     (RILEY AND SKIRROW, 1965, P.250)
  !
  real, parameter :: BOR1=0.000232

  !     -----------------------------------------------------------------
  !*    INVERSE OF ATOMIC WEIGHT OF BORON [G**-1]
  !     (USED TO CONVERT SPECIFIC TOTAL BORAT INTO CONCENTRATIONS)
  !
  real, parameter :: BOR2=1./10.811


  !     -----------------------------------------------------------------
  !*    CONVERSION FACTOR SALINITY -> CHLORINITY
  !     (AFTER WOOSTER ET AL., 1969)
  !
  real, parameter :: SALCHL=1./1.80655
  real, parameter :: rrrcl=salchl*1.025*bor1*bor2

  !     -----------------------------------------------------------------
  !*    ZERO DEG CENTIGRADE AT KELVIN SCALE
  !
  real, parameter :: tzero=273.15

  !     -----------------------------------------------------------------
  !*    SET MEAN TOTAL [CA++] IN SEAWATER (MOLES/KG) (SEE BROECKER
  !     A. PENG, 1982, P. 26; [CA++](MOLES/KG)=1.028E-2*(S/35.); Value
  !     taken from Sarmiento and Gruber, 2006, p. 365
  !
  real, parameter :: CALCON=0.01028

  !     -----------------------------------------------------------------
  !*    INVERS OF NORMAL MOLAL VOLUME OF AN IDEAL GAS [mol/ml] at 0C
  !
  real, parameter :: OXYCO=1./22414.4

  !     -----------------------------------------------------------------
  !*    VOLUMETRIC SOLUBILITY CONSTANTS FOR O2 IN ML/L from moist air at
  !     one atm total pressure. Table 2 in WEISS, R. F. (1970) THE
  !     SOLUBILITY OF NITROGEN OXYGEN AND ARGON IN WATER AND SEAWATER.
  !     DEEP-SEA RESEARCH, VOL. 17, 721-735.
  !
  real, parameter :: OX0=-173.4292
  real, parameter :: OX1=249.6339
  real, parameter :: OX2=143.3483
  real, parameter :: OX3=-21.8492
  real, parameter :: OX4=-0.033096
  real, parameter :: OX5=0.014259
  real, parameter :: OX6=-0.0017

  !     -----------------------------------------------------------------
  !*    VOLUMETRIC SOLUBILITY CONSTANTS FOR N2 IN ML/L from moist air at
  !     one atm total pressure. Table 2 in WEISS, R. F. (1970) THE
  !     SOLUBILITY OF NITROGEN OXYGEN AND ARGON IN WATER AND SEAWATER.
  !     DEEP-SEA RESEARCH, VOL. 17, 721-735.
  !
  real, parameter :: AN0=-172.4965
  real, parameter :: AN1=248.4262
  real, parameter :: AN2=143.0738
  real, parameter :: AN3=-21.7120
  real, parameter :: AN4=-0.049781
  real, parameter :: AN5=0.025018
  real, parameter :: AN6=-0.0034861

  !     -----------------------------------------------------------------
  !     Constants for CO2 solubility in mol/kg/atm from moist
  !     air at one atm total pressure. Table 6 in WEISS, R.F.,
  !     NITROUS OXIDE SOLUBILITY IN WATER AND SEAWATER,
  !     Marine Chemistry, 8, 347-359, 1980
  !
  real, parameter :: ac1= -162.8301
  real, parameter :: ac2=  218.2968
  real, parameter :: ac3=   90.9241
  real, parameter :: ac4=   -1.47696
  real, parameter :: bc1=    0.025695
  real, parameter :: bc2=   -0.025225
  real, parameter :: bc3=    0.0049867

  !     -----------------------------------------------------------------
  !     Constants for CO2 solubility in mol/kg/atm for dry
  !     air at one atm total pressure. Table 1 in WEISS, R.F.,
  !     CARBON DIOXIDE IN WATER AND SEAWATER: THE SOLUBILITY OF
  !     A NON - IDEAL GAS, Marine Chemistry, 2, 203-215, 1974
  !
  real, parameter :: ad1= -60.2409
  real, parameter :: ad2=  93.4517
  real, parameter :: ad3=  23.3585
  real, parameter :: bd1=   0.023517
  real, parameter :: bd2=  -0.023656
  real, parameter :: bd3=   0.0047036

  !     -----------------------------------------------------------------
  !     Constants for laughing gas solubility in mol/l/atm from moist
  !     air at one atm total pressure. Table 2 in WEISS, R.F.,
  !     NITROUS OXIDE SOLUBILITY IN WATER AND SEAWATER,
  !     Marine Chemistry, 8, 347-359, 1980
  !
  real, parameter :: al1= -165.8806
  real, parameter :: al2=  222.8743
  real, parameter :: al3=   92.0792
  real, parameter :: al4=   -1.48425
  real, parameter :: bl1=   -0.056235
  real, parameter :: bl2=    0.031619
  real, parameter :: bl3=   -0.0048472

#ifdef extNcycle
      ! Tsilingiris 2008
      ! moist air dynamic viscosity parameters
      real, parameter :: SV0_air =  1.715747771e-5
      real, parameter :: SV1_air =  4.722402075e-8
      real, parameter :: SV2_air = -3.663027156e-10
      real, parameter :: SV3_air =  1.873236686e-12
      real, parameter :: SV4_air = -8.050218737e-14

      ! moist air density parameters
      real, parameter :: SD0_air =  1.293393662
      real, parameter :: SD1_air = -5.538444326e-3
      real, parameter :: SD2_air =  3.860201577e-5
      real, parameter :: SD3_air = -5.2536065e-7

      ! diffusion of NH3 in water and air 
      real, parameter :: Va_air  = 20.1  ! Johnson 2010
      real, parameter :: Ma_air  = 28.97 ! Johnson 2010
      real, parameter :: Mb_nh3  = 17.03 ! Johnson 2010, Tang 2014
      real, parameter :: Vb_nh3  = 20.7  ! Johnson 2010
      real, parameter :: M_nh3   =  (1./Ma_air + 1./Mb_nh3)**0.5 / (Va_air**(1./3.)+Vb_nh3**(1./3.))**2.
      real, parameter :: kappa   = 0.4   ! von Karman constant

      real, parameter :: mw_nitrogen = 14.00674  ! [g/mol N]   nitrogen mol-weight as defined by CAM
      real, parameter :: mw_nh3      = 17.028940 ! [g/mol NH3] ammonia mol-weight as defined by CAM
      real, parameter :: mw_n2o      = 44.012880 ! [g/mol N2O] nitrous oxide mol-weight as defined by CAM
#endif

  !     -----------------------------------------------------------------
  !     Constants needed for pressure correction of equilibrium constants
  !     F. Millero, Thermodynamics of the carbon dioxide system in the oceans,
  !     Geochimica et Cosmochimica Acta, Vol. 59, No. 4, pp. 661-677, 1995
  real, dimension(11) :: a0, a1, a2, b0, b1, b2
  data a0 /-25.5, -15.82, -29.48, -25.60, -18.03, -9.78, -48.76, &
           -46., -14.51, -23.12, -26.57/
  data a1 /0.1271, -0.0219, 0.1622, 0.2324, 0.0466, -0.0090,     &
           0.5304, 0.5304, 0.1211, 0.1758, 0.2020/
  data a2 /0.0, 0.0, 2.608e-3, -3.6246e-3, 0.316e-3,             &
          -0.942e-3, 0.0, 0.0, -0.321e-3, -2.647e-3, -3.042e-3/
  data b0 /-3.08e-3, 1.13e-3, -2.84e-3, -5.13e-3, -4.53e-3,      &
           -3.91e-3, -11.76e-3, -11.76e-3, -2.67e-3, -5.15e-3,   &
           -4.08e-3/
  data b1 /0.0877e-3, -0.1475e-3, 0.0, 0.0794e-3, 0.09e-3,       &
           0.054e-3, 0.3692e-3, 0.3692e-3, 0.0427e-3,            &
           0.09e-3, 0.0714e-3/
  data b2 /0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

  !     -----------------------------------------------------------------
  !     Gas constant, value as used by Millero (1995)
  real, parameter :: rgas = 83.131

end module mo_chemcon
