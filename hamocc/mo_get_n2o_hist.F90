! Copyright (C) 2020  J. Tjiputra, J. Schwinger, j maerz
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

module mo_get_n2o_hist

  implicit none
  private

  public get_n2o_hist

contains

  subroutine get_n2o_hist(kplyear)
    use mod_xc,         only: mnproc
    use mo_kind,        only: rp
    use mo_control_bgc, only: io_stdo_bgc
    use mo_param1_bgc,  only: iatmn2o
    use mo_carbch,      only: atm
    use mo_control_bgc, only: do_n2o_hist

    ! Arguments
    integer, intent(in) :: kplyear

    ! Local variables
    integer, parameter   :: start_yr = 1750   ! first year of data
    integer, parameter   :: nyears   = 266    ! nb of years in data
    real(rp)             :: atm_n2o_conc(nyears)
    integer, save        :: kplyear_old = 0

    ! N2O atmospheric concentration
    ! Meinshausen et al. 2017: Historical greenhouse gas concentrations for climate modelling (CMIP6)
    ! generated from shared LBC_1750-2015_CMIP6_GlobAnnAvg_c180926.nc file.
    ! First (last) data represents year 1750.5 (2015.5), units are [ppb]

    data atm_n2o_conc /                                                                            &
              273864.96875_rp, 273894.06250_rp, 273927.09375_rp, 273972.81250_rp, 274016.18750_rp, &
              274051.12500_rp, 274079.00000_rp, 274096.03125_rp, 274129.00000_rp, 274165.09375_rp, &
              274193.09375_rp, 274222.00000_rp, 274249.90625_rp, 274301.96875_rp, 274339.93750_rp, &
              274371.96875_rp, 274391.90625_rp, 274426.00000_rp, 274446.81250_rp, 274468.03125_rp, &
              274495.06250_rp, 274513.03125_rp, 274530.93750_rp, 274547.09375_rp, 274565.03125_rp, &
              274578.96875_rp, 274593.93750_rp, 274606.06250_rp, 274617.87500_rp, 274629.00000_rp, &
              274639.21875_rp, 274646.93750_rp, 274654.00000_rp, 274657.93750_rp, 274663.93750_rp, &
              274669.18750_rp, 274676.03125_rp, 274678.06250_rp, 274677.00000_rp, 274678.06250_rp, &
              274677.00000_rp, 274673.96875_rp, 274673.96875_rp, 274659.00000_rp, 274648.03125_rp, &
              274631.03125_rp, 274615.15625_rp, 274601.00000_rp, 274584.18750_rp, 274570.03125_rp, &
              274547.09375_rp, 274531.96875_rp, 274523.00000_rp, 274505.06250_rp, 274485.96875_rp, &
              274465.96875_rp, 274446.81250_rp, 274428.12500_rp, 274411.12500_rp, 274383.18750_rp, &
              274360.06250_rp, 274334.00000_rp, 274304.06250_rp, 274267.12500_rp, 274230.81250_rp, &
              274181.00000_rp, 274143.06250_rp, 274107.12500_rp, 274070.96875_rp, 274050.09375_rp, &
              274007.09375_rp, 273970.00000_rp, 273933.18750_rp, 273894.96875_rp, 273846.96875_rp, &
              273786.06250_rp, 273728.06250_rp, 273677.03125_rp, 273640.87500_rp, 273602.15625_rp, &
              273564.00000_rp, 273520.09375_rp, 273444.03125_rp, 273372.09375_rp, 273214.09375_rp, &
              273078.96875_rp, 272953.96875_rp, 272873.18750_rp, 272866.09375_rp, 272857.96875_rp, &
              272799.03125_rp, 272689.18750_rp, 272582.00000_rp, 272606.96875_rp, 272647.93750_rp, &
              272640.09375_rp, 272730.00000_rp, 272812.96875_rp, 272884.96875_rp, 272953.93750_rp, &
              273021.18750_rp, 273093.96875_rp, 273167.90625_rp, 273264.06250_rp, 273363.00000_rp, &
              273470.00000_rp, 273578.09375_rp, 273675.06250_rp, 273756.03125_rp, 273894.96875_rp, &
              274055.96875_rp, 274240.00000_rp, 274416.12500_rp, 274570.96875_rp, 274719.93750_rp, &
              274880.15625_rp, 275048.21875_rp, 275212.90625_rp, 275390.96875_rp, 275561.18750_rp, &
              275719.93750_rp, 275902.93750_rp, 276078.09375_rp, 276249.93750_rp, 276416.09375_rp, &
              276584.87500_rp, 276734.93750_rp, 276864.09375_rp, 277002.00000_rp, 277133.00000_rp, &
              277265.09375_rp, 277373.87500_rp, 277485.96875_rp, 277590.96875_rp, 277695.06250_rp, &
              277794.93750_rp, 277885.00000_rp, 277994.96875_rp, 278079.84375_rp, 278190.87500_rp, &
              278273.09375_rp, 278346.75000_rp, 278442.09375_rp, 278554.09375_rp, 278689.00000_rp, &
              278831.00000_rp, 278938.90625_rp, 279050.96875_rp, 279163.06250_rp, 279308.12500_rp, &
              279454.06250_rp, 279613.15625_rp, 279860.68750_rp, 280155.96875_rp, 280431.81250_rp, &
              280704.87500_rp, 280980.00000_rp, 281276.03125_rp, 281611.12500_rp, 281949.90625_rp, &
              282314.18750_rp, 282720.96875_rp, 283019.12500_rp, 283361.96875_rp, 283715.87500_rp, &
              284047.00000_rp, 284311.93750_rp, 284614.96875_rp, 284805.00000_rp, 284850.93750_rp, &
              284928.84375_rp, 285038.78125_rp, 285170.00000_rp, 285467.00000_rp, 285605.06250_rp, &
              285652.03125_rp, 285692.09375_rp, 285740.03125_rp, 285832.84375_rp, 285891.03125_rp, &
              285938.12500_rp, 286124.03125_rp, 286221.87500_rp, 286371.18750_rp, 286466.93750_rp, &
              286587.03125_rp, 286747.06250_rp, 286951.00000_rp, 287190.96875_rp, 287387.12500_rp, &
              287619.06250_rp, 287864.12500_rp, 288137.96875_rp, 288781.03125_rp, 289000.09375_rp, &
              289227.06250_rp, 289427.12500_rp, 289510.96875_rp, 289556.21875_rp, 289598.06250_rp, &
              289738.90625_rp, 289860.09375_rp, 290024.93750_rp, 290334.15625_rp, 290547.81250_rp, &
              290844.03125_rp, 291186.96875_rp, 291511.90625_rp, 291772.03125_rp, 291987.00000_rp, &
              292282.96875_rp, 292601.96875_rp, 292945.06250_rp, 293327.03125_rp, 293685.18750_rp, &
              294045.09375_rp, 294453.18750_rp, 294859.93750_rp, 295269.03125_rp, 295681.00000_rp, &
              296098.25000_rp, 296521.90625_rp, 296954.96875_rp, 297399.06250_rp, 297854.93750_rp, &
              298326.12500_rp, 298813.96875_rp, 299318.84375_rp, 299844.93750_rp, 300392.90625_rp, &
              300965.15625_rp, 301562.03125_rp, 302186.96875_rp, 302841.96875_rp, 303528.12500_rp, &
              304247.15625_rp, 305001.96875_rp, 305792.87500_rp, 306623.90625_rp, 307831.12500_rp, &
              308683.00000_rp, 309233.15625_rp, 309725.03125_rp, 310099.09375_rp, 310807.96875_rp, &
              311278.90625_rp, 312298.12500_rp, 313182.93750_rp, 313906.93750_rp, 314708.84375_rp, &
              315759.18750_rp, 316493.00000_rp, 317100.81250_rp, 317729.90625_rp, 318357.00000_rp, &
              319130.03125_rp, 319933.15625_rp, 320645.96875_rp, 321574.84375_rp, 322274.96875_rp, &
              323141.00000_rp, 324159.12500_rp, 325004.78125_rp, 325918.87500_rp, 326987.93750_rp, &
              326987.93750_rp /

    if (do_n2o_hist) then
      if ((kplyear >= start_yr) .and. (kplyear < (start_yr + nyears))) then
        atm(:,:,iatmn2o) = atm_n2o_conc(kplyear - start_yr + 1)

        if (mnproc==1 .and. kplyear > kplyear_old) then
          write(io_stdo_bgc,*) 'ATM N2O CONC=',kplyear,atm_n2o_conc(kplyear - start_yr + 1)
          kplyear_old = kplyear
        endif
      endif
    else
      ! Remain with atm(:,:,iatmn2o) = atm_n2o as initialized in mo_ini_fields
      return
    endif


  end subroutine get_n2o_hist

end module mo_get_n2o_hist

