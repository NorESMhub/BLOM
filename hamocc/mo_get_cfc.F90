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

module mo_get_cfc

  implicit none
  private

  public :: get_cfc

contains

  subroutine get_cfc(kplyear,atm_cfc11_nh,atm_cfc12_nh,atm_sf6_nh,                                 &
                     atm_cfc11_sh,atm_cfc12_sh,atm_sf6_sh)

    !***********************************************************************************************
    !     Jerry Tjiputra         *BCCR*          05.12_rp.2012
    !***********************************************************************************************

    use mo_kind,        only: rp
    use mo_control_bgc, only: io_stdo_bgc
    use mod_xc,         only: mnproc

    ! Arguments
    integer, intent(in)  :: kplyear
    real,    intent(out) :: atm_cfc11_nh,atm_cfc12_nh,atm_sf6_nh
    real,    intent(out) :: atm_cfc11_sh,atm_cfc12_sh,atm_sf6_sh

    ! Local variables
    integer :: i
    integer :: yr_dat(105)
    integer :: start_yr
    real    :: cfc_11_nh(105),cfc_12_nh(105),sf_6_nh(105)
    real    :: cfc_11_sh(105),cfc_12_sh(105),sf_6_sh(105)
    integer, save :: kplyear_old = 0

    ! **********************************************************************************************
    ! Data from EMil Jeansson (Bullister, 2008; Walker et al. 2000; Maiss and Brenninkmeijer (1998)
    ! First (last) data represents year 1910.5 (2014.5), Units are all in [ppt]
    data cfc_11_nh                                                                         &
         & /    0.0_rp,                                                                    &
         &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
         &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
         &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
         &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.1_rp,  0.1_rp,  0.1_rp,  0.2_rp,  0.4_rp,  0.7_rp,  &
         &     1.01_rp, 1.51_rp, 2.21_rp, 3.02_rp, 4.12_rp, 5.33_rp, 6.83_rp, 8.14_rp, 9.45_rp,11.06_rp,  &
         &    13.27_rp,16.18_rp,19.60_rp,23.72_rp,28.44_rp,33.67_rp,39.40_rp,46.03_rp,53.77_rp,62.41_rp,  &
         &    72.06_rp,  82.71_rp,  94.87_rp, 108.34_rp, 121.41_rp,                        &
         &   133.97_rp, 145.93_rp, 156.58_rp, 168.34_rp, 176.68_rp,                        &
         &   184.32_rp, 191.46_rp, 199.30_rp, 208.04_rp, 217.99_rp,                        &
         &   229.35_rp, 241.61_rp, 252.86_rp, 259.30_rp, 265.83_rp,                        &
         &   268.24_rp, 268.14_rp, 269.55_rp, 269.65_rp, 268.34_rp,                        &
         &   266.93_rp, 265.73_rp, 264.52_rp, 263.12_rp, 261.71_rp,                        &
         &   260.00_rp, 258.19_rp, 256.18_rp, 253.97_rp, 251.96_rp,                        &
         &   249.55_rp, 247.54_rp, 245.63_rp, 243.61_rp, 241.33_rp,                        &
         &   239.41_rp, 236.60_rp, 235.08_rp, 233.55_rp/

    data cfc_11_sh                                                                         &
         & /    0.0_rp,                                                                    &
         &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
         &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
         &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
         &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.1_rp,  0.1_rp,  0.1_rp,  0.2_rp,  0.4_rp,  &
         &      0.7_rp, 1.01_rp, 1.51_rp, 2.21_rp, 3.02_rp, 4.02_rp, 5.23_rp, 6.53_rp, 7.84_rp, 9.15_rp,  &
         &    10.85_rp,13.07_rp,15.78_rp,19.20_rp,23.12_rp,27.64_rp,32.66_rp,38.29_rp,44.82_rp,52.26_rp,  &
         &    60.70_rp,  69.95_rp,  80.40_rp,  92.16_rp, 104.72_rp,                        &
         &   117.09_rp, 129.35_rp, 140.80_rp, 148.74_rp, 159.30_rp,                        &
         &   167.84_rp, 176.08_rp, 184.52_rp, 192.46_rp, 202.01_rp,                        &
         &   211.36_rp, 222.21_rp, 233.27_rp, 242.11_rp, 251.06_rp,                        &
         &   256.68_rp, 260.80_rp, 262.51_rp, 263.72_rp, 263.22_rp,                        &
         &   262.91_rp, 262.01_rp, 261.01_rp, 259.90_rp, 258.29_rp,                        &
         &   256.98_rp, 255.08_rp, 253.27_rp, 251.36_rp, 249.15_rp,                        &
         &   247.34_rp, 245.03_rp, 243.12_rp, 241.07_rp, 239.19_rp,                        &
         &   236.92_rp, 234.60_rp, 233.29_rp, 231.97_rp/

    data cfc_12_nh                                                                         &
         & /    0.0_rp,                                                                    &
         &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
         &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
         &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.1_rp,  0.1_rp,  0.2_rp,  0.3_rp,  0.4_rp,  &
         &      0.5_rp,  0.7_rp,  0.9_rp,  1.2_rp,  1.7_rp,  2.3_rp,  3.4_rp,  4.8_rp,  6.1_rp,  7.6_rp,  &
         &      9.2_rp, 11.0_rp, 12.8_rp, 15.0_rp, 17.4_rp, 20.2_rp, 23.4_rp, 26.8_rp, 30.5_rp, 35.0_rp,  &
         &     40.0_rp, 45.8_rp, 52.5_rp, 60.4_rp, 69.3_rp, 79.2_rp, 90.3_rp,102.8_rp,116.8_rp,132.00_rp, &
         &   148.40_rp, 166.10_rp, 185.80_rp, 207.10_rp, 228.20_rp,                        &
         &   248.10_rp, 266.90_rp, 284.30_rp, 306.10_rp, 323.20_rp,                        &
         &   339.60_rp, 353.40_rp, 369.00_rp, 385.70_rp, 403.40_rp,                        &
         &   424.30_rp, 444.00_rp, 465.40_rp, 483.60_rp, 497.70_rp,                        &
         &   506.00_rp, 516.30_rp, 523.20_rp, 528.50_rp, 533.40_rp,                        &
         &   537.30_rp, 540.10_rp, 542.90_rp, 544.40_rp, 545.90_rp,                        &
         &   546.50_rp, 546.70_rp, 546.70_rp, 545.70_rp, 544.90_rp,                        &
         &   543.10_rp, 541.10_rp, 538.60_rp, 536.11_rp, 533.30_rp,                        &
         &   530.67_rp, 527.16_rp, 525.26_rp, 523.36_rp/

    data cfc_12_sh                                                                         &
         & /    0.0_rp,                                                                    &
         &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
         &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
         &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.1_rp,  0.1_rp,  0.2_rp,  0.3_rp,  &
         &      0.4_rp,  0.5_rp,  0.7_rp,  0.9_rp,  1.2_rp,  1.7_rp,  2.4_rp,  3.4_rp,  4.7_rp,  6.0_rp,  &
         &      7.4_rp,  9.0_rp, 10.7_rp, 12.6_rp, 14.7_rp, 17.1_rp, 19.9_rp, 23.0_rp, 26.3_rp, 30.1_rp,  &
         &     34.4_rp, 39.4_rp, 45.1_rp, 51.8_rp, 59.5_rp, 68.2_rp, 77.9_rp, 88.8_rp,101.1_rp,114.7_rp,  &
         &    129.6_rp,145.7_rp,163.3_rp,182.5_rp,202.9_rp,223.2_rp,242.7_rp,261.2_rp,273.5_rp,292.3_rp,  &
         &    308.8_rp,325.5_rp,342.6_rp,359.4_rp,378.2_rp,396.5_rp,416.3_rp,435.8_rp,454.4_rp,472.7_rp,  &
         &    487.3_rp,498.3_rp,507.0_rp,514.8_rp,521.0_rp,526.5_rp,530.8_rp,534.3_rp,537.2_rp,539.0_rp,  &
         &    540.6_rp,  541.3_rp,  541.6_rp,  541.5_rp,  540.7_rp,                        &
         &    539.8_rp,  538.1_rp,  536.2_rp, 533.53_rp, 530.94_rp,                        &
         &   528.47_rp, 525.88_rp, 523.48_rp, 521.08_rp/

    data sf_6_nh                                                                           &
         & /   0.00_rp,                                                                    &
         &     0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp,  &
         &     0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp,  &
         &     0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp,  &
         &     0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp,  &
         &    0.000_rp, 0.000_rp, 0.042_rp, 0.043_rp, 0.043_rp,                            &
         &    0.044_rp, 0.046_rp, 0.048_rp, 0.051_rp, 0.055_rp,                            &
         &    0.061_rp, 0.068_rp, 0.078_rp, 0.091_rp, 0.109_rp,                            &
         &    0.131_rp, 0.155_rp, 0.181_rp, 0.207_rp, 0.235_rp,                            &
         &    0.266_rp, 0.301_rp, 0.341_rp, 0.386_rp, 0.438_rp,                            &
         &    0.501_rp, 0.579_rp, 0.665_rp, 0.766_rp, 0.887_rp,                            &
         &    1.011_rp, 1.141_rp, 1.273_rp, 1.409_rp, 1.562_rp,                            &
         &    1.722_rp, 1.892_rp, 2.063_rp, 2.237_rp, 2.427_rp,                            &
         &    2.640_rp, 2.868_rp, 3.104_rp, 3.350_rp, 3.600_rp,                            &
         &    3.861_rp, 4.080_rp, 4.262_rp, 4.485_rp, 4.690_rp,                            &
         &    4.909_rp, 5.135_rp, 5.360_rp, 5.580_rp, 5.795_rp,                            &
         &    6.034_rp, 6.324_rp, 6.613_rp, 6.876_rp, 7.191_rp,                            &
         &    7.439_rp, 7.715_rp, 8.066_rp, 8.417_rp/

    data sf_6_sh                                                                           &
         & /   0.00_rp,                                                                    &
         &     0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp,  &
         &     0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp,  &
         &     0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp,  &
         &     0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp,  &
         &     0.000_rp, 0.000_rp, 0.039_rp, 0.039_rp, 0.040_rp,                           &
         &     0.041_rp, 0.042_rp, 0.044_rp, 0.047_rp, 0.051_rp,                           &
         &     0.056_rp, 0.062_rp, 0.071_rp, 0.084_rp, 0.100_rp,                           &
         &     0.120_rp, 0.142_rp, 0.166_rp, 0.190_rp, 0.215_rp,                           &
         &     0.243_rp, 0.276_rp, 0.312_rp, 0.354_rp, 0.401_rp,                           &
         &     0.459_rp, 0.531_rp, 0.610_rp, 0.703_rp, 0.813_rp,                           &
         &     0.927_rp, 1.046_rp, 1.167_rp, 1.292_rp, 1.432_rp,                           &
         &     1.579_rp, 1.735_rp, 1.892_rp, 2.051_rp, 2.225_rp,                           &
         &     2.420_rp, 2.629_rp, 2.846_rp, 3.071_rp, 3.300_rp,                           &
         &     3.560_rp, 3.824_rp, 4.026_rp, 4.262_rp, 4.471_rp,                           &
         &     4.657_rp, 4.887_rp, 5.081_rp, 5.305_rp, 5.513_rp,                           &
         &     5.749_rp, 6.028_rp, 6.286_rp, 6.576_rp, 6.856_rp,                           &
         &     7.159_rp, 7.424_rp, 7.754_rp, 8.084_rp/

    start_yr=1910
    do i=1,105
      yr_dat(i)=start_yr+i-1
    enddo

    ! ******************************************************************
    !if (kplyear.lt.start_yr) then
    atm_cfc11_nh=0.0_rp
    atm_cfc11_sh=0.0_rp
    atm_cfc12_nh=0.0_rp
    atm_cfc12_sh=0.0_rp
    atm_sf6_nh=0.0_rp
    atm_sf6_sh=0.0_rp

    do i=1,105
      if (kplyear.eq.yr_dat(i)) then
        atm_cfc11_nh=cfc_11_nh(i)
        atm_cfc11_sh=cfc_11_sh(i)
        atm_cfc12_nh=cfc_12_nh(i)
        atm_cfc12_sh=cfc_12_sh(i)
        atm_sf6_nh=sf_6_nh(i)
        atm_sf6_sh=sf_6_sh(i)
      endif
    enddo

    if (mnproc==1 .and. kplyear > kplyear_old) then
      write(io_stdo_bgc,*) 'ATM NH CFC11, CFC12, SF6=',             &
           &    kplyear,atm_cfc11_nh,atm_cfc12_nh,atm_sf6_nh
      write(io_stdo_bgc,*) 'ATM SH CFC11, CFC12, SF6=',             &
           &    kplyear,atm_cfc11_sh,atm_cfc12_sh,atm_sf6_sh
      kplyear_old = kplyear
    endif

  end subroutine get_cfc

end module mo_get_cfc
