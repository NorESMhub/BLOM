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
    data cfc_11_nh                                                          &
         & /    0.0,                                                        &
         &      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  &
         &      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  &
         &      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  &
         &      0.0,  0.0,  0.0,  0.0,  0.1,  0.1,  0.1,  0.2,  0.4,  0.7,  &
         &     1.01, 1.51, 2.21, 3.02, 4.12, 5.33, 6.83, 8.14, 9.45,11.06,  &
         &    13.27,16.18,19.60,23.72,28.44,33.67,39.40,46.03,53.77,62.41,  &
         &    72.06,  82.71,  94.87, 108.34, 121.41,                        &
         &   133.97, 145.93, 156.58, 168.34, 176.68,                        &
         &   184.32, 191.46, 199.30, 208.04, 217.99,                        &
         &   229.35, 241.61, 252.86, 259.30, 265.83,                        &
         &   268.24, 268.14, 269.55, 269.65, 268.34,                        &
         &   266.93, 265.73, 264.52, 263.12, 261.71,                        &
         &   260.00, 258.19, 256.18, 253.97, 251.96,                        &
         &   249.55, 247.54, 245.63, 243.61, 241.33,                        &
         &   239.41, 236.60, 235.08, 233.55/

    data cfc_11_sh                                                          &
         & /    0.0,                                                        &
         &      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  &
         &      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  &
         &      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  &
         &      0.0,  0.0,  0.0,  0.0,  0.0,  0.1,  0.1,  0.1,  0.2,  0.4,  &
         &      0.7, 1.01, 1.51, 2.21, 3.02, 4.02, 5.23, 6.53, 7.84, 9.15,  &
         &    10.85,13.07,15.78,19.20,23.12,27.64,32.66,38.29,44.82,52.26,  &
         &    60.70,  69.95,  80.40,  92.16, 104.72,                        &
         &   117.09, 129.35, 140.80, 148.74, 159.30,                        &
         &   167.84, 176.08, 184.52, 192.46, 202.01,                        &
         &   211.36, 222.21, 233.27, 242.11, 251.06,                        &
         &   256.68, 260.80, 262.51, 263.72, 263.22,                        &
         &   262.91, 262.01, 261.01, 259.90, 258.29,                        &
         &   256.98, 255.08, 253.27, 251.36, 249.15,                        &
         &   247.34, 245.03, 243.12, 241.07, 239.19,                        &
         &   236.92, 234.60, 233.29, 231.97/

    data cfc_12_nh                                                          &
         & /    0.0,                                                        &
         &      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  &
         &      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  &
         &      0.0,  0.0,  0.0,  0.0,  0.0,  0.1,  0.1,  0.2,  0.3,  0.4,  &
         &      0.5,  0.7,  0.9,  1.2,  1.7,  2.3,  3.4,  4.8,  6.1,  7.6,  &
         &      9.2, 11.0, 12.8, 15.0, 17.4, 20.2, 23.4, 26.8, 30.5, 35.0,  &
         &     40.0, 45.8, 52.5, 60.4, 69.3, 79.2, 90.3,102.8,116.8,132.00, &
         &   148.40, 166.10, 185.80, 207.10, 228.20,                        &
         &   248.10, 266.90, 284.30, 306.10, 323.20,                        &
         &   339.60, 353.40, 369.00, 385.70, 403.40,                        &
         &   424.30, 444.00, 465.40, 483.60, 497.70,                        &
         &   506.00, 516.30, 523.20, 528.50, 533.40,                        &
         &   537.30, 540.10, 542.90, 544.40, 545.90,                        &
         &   546.50, 546.70, 546.70, 545.70, 544.90,                        &
         &   543.10, 541.10, 538.60, 536.11, 533.30,                        &
         &   530.67, 527.16, 525.26, 523.36/

    data cfc_12_sh                                                          &
         & /    0.0,                                                        &
         &      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  &
         &      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  &
         &      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.1,  0.1,  0.2,  0.3,  &
         &      0.4,  0.5,  0.7,  0.9,  1.2,  1.7,  2.4,  3.4,  4.7,  6.0,  &
         &      7.4,  9.0, 10.7, 12.6, 14.7, 17.1, 19.9, 23.0, 26.3, 30.1,  &
         &     34.4, 39.4, 45.1, 51.8, 59.5, 68.2, 77.9, 88.8,101.1,114.7,  &
         &    129.6,145.7,163.3,182.5,202.9,223.2,242.7,261.2,273.5,292.3,  &
         &    308.8,325.5,342.6,359.4,378.2,396.5,416.3,435.8,454.4,472.7,  &
         &    487.3,498.3,507.0,514.8,521.0,526.5,530.8,534.3,537.2,539.0,  &
         &    540.6,  541.3,  541.6,  541.5,  540.7,                        &
         &    539.8,  538.1,  536.2, 533.53, 530.94,                        &
         &   528.47, 525.88, 523.48, 521.08/

    data sf_6_nh                                                            &
         & /   0.00,                                                        &
         &     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,  &
         &     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,  &
         &     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,  &
         &     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,  &
         &    0.000, 0.000, 0.042, 0.043, 0.043,                            &
         &    0.044, 0.046, 0.048, 0.051, 0.055,                            &
         &    0.061, 0.068, 0.078, 0.091, 0.109,                            &
         &    0.131, 0.155, 0.181, 0.207, 0.235,                            &
         &    0.266, 0.301, 0.341, 0.386, 0.438,                            &
         &    0.501, 0.579, 0.665, 0.766, 0.887,                            &
         &    1.011, 1.141, 1.273, 1.409, 1.562,                            &
         &    1.722, 1.892, 2.063, 2.237, 2.427,                            &
         &    2.640, 2.868, 3.104, 3.350, 3.600,                            &
         &    3.861, 4.080, 4.262, 4.485, 4.690,                            &
         &    4.909, 5.135, 5.360, 5.580, 5.795,                            &
         &    6.034, 6.324, 6.613, 6.876, 7.191,                            &
         &    7.439, 7.715, 8.066, 8.417/

    data sf_6_sh                                                            &
         & /   0.00,                                                        &
         &     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,  &
         &     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,  &
         &     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,  &
         &     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,  &
         &     0.000, 0.000, 0.039, 0.039, 0.040,                           &
         &     0.041, 0.042, 0.044, 0.047, 0.051,                           &
         &     0.056, 0.062, 0.071, 0.084, 0.100,                           &
         &     0.120, 0.142, 0.166, 0.190, 0.215,                           &
         &     0.243, 0.276, 0.312, 0.354, 0.401,                           &
         &     0.459, 0.531, 0.610, 0.703, 0.813,                           &
         &     0.927, 1.046, 1.167, 1.292, 1.432,                           &
         &     1.579, 1.735, 1.892, 2.051, 2.225,                           &
         &     2.420, 2.629, 2.846, 3.071, 3.300,                           &
         &     3.560, 3.824, 4.026, 4.262, 4.471,                           &
         &     4.657, 4.887, 5.081, 5.305, 5.513,                           &
         &     5.749, 6.028, 6.286, 6.576, 6.856,                           &
         &     7.159, 7.424, 7.754, 8.084/

    start_yr=1910
    do i=1,105
      yr_dat(i)=start_yr+i-1
    enddo

    ! ******************************************************************
    !if (kplyear.lt.start_yr) then
    atm_cfc11_nh=0.0
    atm_cfc11_sh=0.0
    atm_cfc12_nh=0.0
    atm_cfc12_sh=0.0
    atm_sf6_nh=0.0
    atm_sf6_sh=0.0

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
