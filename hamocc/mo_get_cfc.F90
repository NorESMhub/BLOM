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
    !     Jerry Tjiputra         *BCCR*          05.12.2012
    !***********************************************************************************************

    use mo_kind,        only: rp
    use mo_control_bgc, only: io_stdo_bgc
    use mod_xc,         only: mnproc

    ! Arguments
    integer, intent(in)  :: kplyear
    real(rp),intent(out) :: atm_cfc11_nh,atm_cfc12_nh,atm_sf6_nh
    real(rp),intent(out) :: atm_cfc11_sh,atm_cfc12_sh,atm_sf6_sh

    ! Local variables
    integer, parameter   :: start_yr = 1910   ! first year of data
    integer, parameter   :: nyears   = 113    ! nb of years in data
    integer              :: i
    integer              :: yr_dat(nyears)
    real(rp)             :: cfc_11_nh(nyears),cfc_12_nh(nyears),sf_6_nh(nyears)
    real(rp)             :: cfc_11_sh(nyears),cfc_12_sh(nyears),sf_6_sh(nyears)
    integer, save        :: kplyear_old = 0

    ! **********************************************************************************************
    ! Data from taken from Bullister, John L.; Warner, Mark J. (2017).
    ! Atmospheric Histories (1765-2022) for CFC-11, CFC-12, CFC-113, CCl4, SF6 and
    ! N2O. NOAA National Centers for Environmental Information. Dataset. 
    ! https://doi.org/10.3334/cdiac/otg.cfc_atm_hist_2015.
    ! First (last) data represents year 1910.5 (2022.5), Units are all in [ppt]
    ! **********************************************************************************************
    data cfc_11_nh                                                                                    &
     & /    0.0_rp,                                                                                   &
     &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
     &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
     &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
     &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.1_rp,  0.1_rp,  0.1_rp,  0.2_rp,  0.4_rp,  0.7_rp,  &
     &     1.01_rp, 1.51_rp, 2.21_rp, 3.02_rp, 4.12_rp, 5.33_rp, 6.83_rp, 8.14_rp, 9.45_rp,11.06_rp,  &
     &    13.27_rp,16.18_rp,19.60_rp,23.72_rp,28.44_rp,33.67_rp,39.40_rp,46.03_rp,53.77_rp,62.41_rp,  &
     &    72.06_rp,  82.71_rp,  94.87_rp, 108.34_rp, 121.41_rp,                                       &
     &   133.97_rp, 145.93_rp, 156.58_rp, 168.34_rp, 176.68_rp,                                       &
     &   184.32_rp, 191.46_rp, 199.30_rp, 208.04_rp, 217.99_rp,                                       &
     &   229.35_rp, 241.61_rp, 252.86_rp, 259.30_rp, 265.83_rp,                                       &
     &   268.24_rp, 268.14_rp, 269.55_rp, 269.65_rp, 268.34_rp,                                       &
     &   266.93_rp, 265.73_rp, 264.52_rp, 263.12_rp, 261.71_rp,                                       &
     &   260.00_rp, 258.19_rp, 256.18_rp, 253.97_rp, 251.96_rp,                                       &
     &   249.55_rp, 247.54_rp, 245.63_rp, 243.78_rp, 241.94_rp,                                       &
     &   239.89_rp, 237.83_rp, 236.36_rp, 235.38_rp, 234.14_rp,                                       &
     &   233.18_rp, 232.30_rp, 231.27_rp, 229.33_rp, 226.58_rp,                                       &
     &   224.01_rp, 221.71_rp/

    data cfc_11_sh                                                                                    &
     & /    0.0_rp,                                                                                   &
     &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
     &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
     &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
     &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.1_rp,  0.1_rp,  0.1_rp,  0.2_rp,  0.4_rp,  &
     &      0.7_rp, 1.01_rp, 1.51_rp, 2.21_rp, 3.02_rp, 4.02_rp, 5.23_rp, 6.53_rp, 7.84_rp, 9.15_rp,  &
     &    10.85_rp,13.07_rp,15.78_rp,19.20_rp,23.12_rp,27.64_rp,32.66_rp,38.29_rp,44.82_rp,52.26_rp,  &
     &    60.70_rp,  69.95_rp,  80.40_rp,  92.16_rp, 104.72_rp,                                       &
     &   117.09_rp, 129.35_rp, 140.80_rp, 148.74_rp, 159.30_rp,                                       &
     &   167.84_rp, 176.08_rp, 184.52_rp, 192.46_rp, 202.01_rp,                                       &
     &   211.36_rp, 222.21_rp, 233.27_rp, 242.11_rp, 251.06_rp,                                       &
     &   256.68_rp, 260.80_rp, 262.51_rp, 263.72_rp, 263.22_rp,                                       &
     &   262.91_rp, 262.01_rp, 261.01_rp, 259.90_rp, 258.29_rp,                                       &
     &   256.98_rp, 255.08_rp, 253.27_rp, 251.36_rp, 249.15_rp,                                       &
     &   247.34_rp, 245.03_rp, 243.12_rp, 241.29_rp, 239.46_rp,                                       &
     &   237.36_rp, 235.40_rp, 233.52_rp, 231.81_rp, 230.65_rp,                                       &
     &   229.64_rp, 228.58_rp, 228.08_rp, 226.51_rp, 224.24_rp,                                       &
     &   222.19_rp, 220.08_rp/

    data cfc_12_nh                                                                                    &
     & /    0.0_rp,                                                                                   &
     &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
     &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
     &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.1_rp,  0.1_rp,  0.2_rp,  0.3_rp,  0.4_rp,  &
     &      0.5_rp,  0.7_rp,  0.9_rp,  1.2_rp,  1.7_rp,  2.3_rp,  3.4_rp,  4.8_rp,  6.1_rp,  7.6_rp,  &
     &      9.2_rp, 11.0_rp, 12.8_rp, 15.0_rp, 17.4_rp, 20.2_rp, 23.4_rp, 26.8_rp, 30.5_rp, 35.0_rp,  &
     &     40.0_rp, 45.8_rp, 52.5_rp, 60.4_rp, 69.3_rp, 79.2_rp, 90.3_rp,102.8_rp,116.8_rp,132.00_rp, &
     &   148.40_rp, 166.10_rp, 185.80_rp, 207.10_rp, 228.20_rp,                                       &
     &   248.10_rp, 266.90_rp, 284.30_rp, 306.10_rp, 323.20_rp,                                       &
     &   339.60_rp, 353.40_rp, 369.00_rp, 385.70_rp, 403.40_rp,                                       &
     &   424.30_rp, 444.00_rp, 465.40_rp, 483.60_rp, 497.70_rp,                                       &
     &   506.00_rp, 516.30_rp, 523.20_rp, 528.50_rp, 533.40_rp,                                       &
     &   537.30_rp, 540.10_rp, 542.90_rp, 544.40_rp, 545.90_rp,                                       &
     &   546.50_rp, 546.70_rp, 546.70_rp, 545.70_rp, 544.90_rp,                                       &
     &   543.10_rp, 541.10_rp, 538.60_rp, 536.02_rp, 532.46_rp,                                       &
     &   529.72_rp, 527.00_rp, 523.48_rp, 520.77_rp, 517.93_rp,                                       &
     &   514.76_rp, 511.08_rp, 508.06_rp, 503.82_rp, 499.81_rp,                                       &
     &   495.61_rp, 491.93_rp/

    data cfc_12_sh                                                                                    &
     & /    0.0_rp,                                                                                   &
     &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
     &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  &
     &      0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.1_rp,  0.1_rp,  0.2_rp,  0.3_rp,  &
     &      0.4_rp,  0.5_rp,  0.7_rp,  0.9_rp,  1.2_rp,  1.7_rp,  2.4_rp,  3.4_rp,  4.7_rp,  6.0_rp,  &
     &      7.4_rp,  9.0_rp, 10.7_rp, 12.6_rp, 14.7_rp, 17.1_rp, 19.9_rp, 23.0_rp, 26.3_rp, 30.1_rp,  &
     &     34.4_rp, 39.4_rp, 45.1_rp, 51.8_rp, 59.5_rp, 68.2_rp, 77.9_rp, 88.8_rp,101.1_rp,114.7_rp,  &
     &    129.6_rp,145.7_rp,163.3_rp,182.5_rp,202.9_rp,223.2_rp,242.7_rp,261.2_rp,273.5_rp,292.3_rp,  &
     &    308.8_rp,325.5_rp,342.6_rp,359.4_rp,378.2_rp,396.5_rp,416.3_rp,435.8_rp,454.4_rp,472.7_rp,  &
     &    487.3_rp,498.3_rp,507.0_rp,514.8_rp,521.0_rp,526.5_rp,530.8_rp,534.3_rp,537.2_rp,539.0_rp,  &
     &    540.6_rp,  541.3_rp,  541.6_rp,  541.5_rp,  540.7_rp,                                       &
     &    539.8_rp,  538.1_rp,  536.2_rp, 533.71_rp, 530.20_rp,                                       &
     &   527.58_rp, 524.98_rp, 521.72_rp, 518.62_rp, 515.57_rp,                                       &
     &   512.94_rp, 509.63_rp, 507.01_rp, 502.69_rp, 498.39_rp,                                       &
     &   494.81_rp, 491.46_rp/

    data sf_6_nh                                                                                      &
     & /   0.00_rp,                                                                                   &
     &     0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp,  &
     &     0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp,  &
     &     0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp,  &
     &     0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp,  &
     &     0.00_rp, 0.00_rp, 0.04_rp, 0.04_rp, 0.04_rp, 0.04_rp, 0.05_rp, 0.05_rp, 0.05_rp, 0.05_rp,  &
     &     0.06_rp, 0.07_rp, 0.08_rp, 0.09_rp, 0.11_rp, 0.13_rp, 0.15_rp, 0.18_rp, 0.21_rp, 0.23_rp,  &
     &     0.26_rp, 0.30_rp, 0.34_rp, 0.38_rp, 0.44_rp, 0.50_rp, 0.58_rp, 0.66_rp, 0.76_rp, 0.88_rp,  &
     &     1.00_rp, 1.13_rp, 1.27_rp, 1.40_rp, 1.55_rp, 1.71_rp, 1.88_rp, 2.05_rp, 2.22_rp, 2.41_rp,  &
     &     2.62_rp, 2.85_rp, 3.09_rp, 3.33_rp, 3.55_rp, 3.82_rp, 4.05_rp, 4.26_rp, 4.46_rp, 4.65_rp,  &
     &     4.86_rp, 5.11_rp, 5.35_rp, 5.58_rp, 5.80_rp, 6.04_rp, 6.32_rp, 6.62_rp, 6.91_rp, 7.21_rp,  &
     &     7.47_rp, 7.76_rp, 8.09_rp, 8.43_rp, 8.76_rp, 9.10_rp, 9.44_rp, 9.78_rp,10.12_rp,10.44_rp,  &
     &    10.84_rp,11.23_rp/

    data sf_6_sh                                                                                      &
     & /   0.00_rp,                                                                                   &
     &     0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp,  &
     &     0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp,  &
     &     0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp,  &
     &     0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp, 0.00_rp,  &
     &     0.00_rp, 0.00_rp, 0.04_rp, 0.04_rp, 0.04_rp, 0.04_rp, 0.04_rp, 0.04_rp, 0.05_rp, 0.05_rp,  &
     &     0.06_rp, 0.06_rp, 0.07_rp, 0.08_rp, 0.10_rp, 0.12_rp, 0.14_rp, 0.17_rp, 0.19_rp, 0.22_rp,  &
     &     0.24_rp, 0.28_rp, 0.31_rp, 0.35_rp, 0.40_rp, 0.46_rp, 0.53_rp, 0.61_rp, 0.70_rp, 0.81_rp,  &
     &     0.93_rp, 1.04_rp, 1.17_rp, 1.29_rp, 1.43_rp, 1.58_rp, 1.73_rp, 1.89_rp, 2.05_rp, 2.22_rp,  &
     &     2.42_rp, 2.63_rp, 2.84_rp, 3.07_rp, 3.30_rp, 3.55_rp, 3.81_rp, 4.01_rp, 4.25_rp, 4.46_rp,  &
     &     4.65_rp, 4.86_rp, 5.08_rp, 5.31_rp, 5.53_rp, 5.76_rp, 6.03_rp, 6.30_rp, 6.59_rp, 6.87_rp,  &
     &     7.15_rp, 7.43_rp, 7.73_rp, 8.06_rp, 8.40_rp, 8.72_rp, 9.06_rp, 9.42_rp, 9.76_rp,10.08_rp,  &
     &    10.44_rp,10.82_rp/

    do i=1,nyears
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

    do i=1,nyears
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
