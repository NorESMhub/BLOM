      MODULE mo_chemcon

!**********************************************************************
!
!**** *MODULE mo_chemcon* - Parameter definitions for chemical formulas
!
!     J. Schwinger,    *UiB-GfI, Bergen*    2013-08-21
!
!     Modified
!     --------
!     
!     
!     Purpose
!     -------
!     - declare chemical parameters previously defined in
!       subroutine chemcon
!
!**********************************************************************


      implicit none


!      real, parameter :: ZERO=0.
!      real, parameter :: TENM7=10.**(-7.0)
!      real, parameter :: SMICR=1.E-6
!      real, parameter :: THOUSI=1./1000.
!      real, parameter :: PERC=0.01
!      real, parameter :: FOURTH=0.25
!      real, parameter :: THIRD=1./3.
!      real, parameter :: HALF=0.5
!      real, parameter :: ONE=1.
!      real, parameter :: TWO=2.
!      real, parameter :: TEN=10.


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
!      Constants for CO2 solubility in mol/kg/atm from moist 
!      air at one atm total pressure. Table 6 in WEISS, R.F.,
!      NITROUS OXIDE SOLUBILITY IN WATER AND SEAWATER, 
!      Marine Chemistry, 8, 347-359, 1980


       real, parameter :: ac1= -162.8301
       real, parameter :: ac2=  218.2968
       real, parameter :: ac3=   90.9241
       real, parameter :: ac4=   -1.47696
       real, parameter :: bc1=    0.025695
       real, parameter :: bc2=   -0.025225
       real, parameter :: bc3=    0.0049867



!     -----------------------------------------------------------------
!      Constants for laughing gas solubility in mol/l/atm from moist 
!      air at one atm total pressure. Table 2 in WEISS, R.F.,
!      NITROUS OXIDE SOLUBILITY IN WATER AND SEAWATER, 
!      Marine Chemistry, 8, 347-359, 1980


       real, parameter :: al1= -165.8806
       real, parameter :: al2=  222.8743
       real, parameter :: al3=   92.0792
       real, parameter :: al4=   -1.48425
       real, parameter :: bl1=   -0.056235
       real, parameter :: bl2=    0.031619
       real, parameter :: bl3=   -0.0048472


!      Atmospheric mixing ratio of N2O around 1980 300 ppb
       real, parameter :: atn2o=3.e-7
       

       END MODULE mo_chemcon
