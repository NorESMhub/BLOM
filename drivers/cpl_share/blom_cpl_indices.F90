! ------------------------------------------------------------------------------
! Copyright (C) 2011-2020 Mats Bentsen, Alok Kumar Gupta, Jerry Tjiputra,
!                         Ingo Bethke
!
! This file is part of BLOM.
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
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

module blom_cpl_indices
 
  use seq_flds_mod, only : seq_flds_x2o_fields, seq_flds_o2x_fields 
!  use seq_flds_mod, only : seq_flds_i2o_per_cat, ice_ncat 
!  use mcog, only : mcog_ncols, lmcog_flds_sent
  use mct_mod

  implicit none

  SAVE
  public                               ! By default make data private

  ! ocn -> drv

  integer :: index_o2x_So_t      
  integer :: index_o2x_So_u
  integer :: index_o2x_So_v
  integer :: index_o2x_So_s
  integer :: index_o2x_So_dhdx
  integer :: index_o2x_So_dhdy
  ! QL, 150526, to wav, boundary layer depth
  integer :: index_o2x_So_bldepth
  integer :: index_o2x_Fioo_q
  integer :: index_o2x_Faoo_fco2_ocn
  integer :: index_o2x_Faoo_fdms_ocn

  ! drv -> ocn

  integer :: index_x2o_Si_ifrac        ! fractional ice wrt ocean
  integer :: index_x2o_So_duu10n       ! 10m wind speed squared           (m^2/s^2)
  integer :: index_x2o_Sa_pslv         ! sea-level pressure               (Pa)
  integer :: index_x2o_Sa_co2prog      ! bottom atm level prognostic CO2
  integer :: index_x2o_Sa_co2diag      ! bottom atm level diagnostic CO2
  integer :: index_x2o_Faxa_nhx        ! nitrogen deposition (nhx) flux from atm (kgNm2/sec)
  integer :: index_x2o_Faxa_noy        ! nitrogen deposition (noy) flux from atm (kgNm2/sec)

  ! QL, 150526, from wav
  integer :: index_x2o_Sw_lamult       ! wave model langmuir multiplier
  integer :: index_x2o_Sw_ustokes      ! surface Stokes drift, x-component
  integer :: index_x2o_Sw_vstokes      ! surface Stokes drift, y-component
  integer :: index_x2o_Foxx_taux       ! zonal wind stress (taux)         (W/m2   )
  integer :: index_x2o_Foxx_tauy       ! meridonal wind stress (tauy)     (W/m2   )
  integer :: index_x2o_Foxx_swnet      ! net short-wave heat flux         (W/m2   )
  integer :: index_x2o_Foxx_sen        ! sensible heat flux               (W/m2   )
  integer :: index_x2o_Foxx_lat        
  integer :: index_x2o_Foxx_lwup       ! longwave radiation (up)          (W/m2   )
  integer :: index_x2o_Faxa_lwdn       ! longwave radiation (down)        (W/m2   )
  integer :: index_x2o_Fioi_melth      ! heat flux from snow & ice melt   (W/m2   )
  integer :: index_x2o_Fioi_meltw      ! snow melt flux                   (kg/m2/s)
  integer :: index_x2o_Fioi_bcpho      ! flux: Black Carbon hydrophobic release from sea ice component
  integer :: index_x2o_Fioi_bcphi      ! flux: Black Carbon hydrophilic release from sea ice component
  integer :: index_x2o_Fioi_flxdst     ! flux: dust release from sea ice component
  integer :: index_x2o_Fioi_salt       ! salt                             (kg(salt)/m2/s)
  integer :: index_x2o_Foxx_evap       ! evaporation flux                 (kg/m2/s)
  integer :: index_x2o_Faxa_prec         
  integer :: index_x2o_Faxa_snow       ! water flux due to snow           (kg/m2/s)
  integer :: index_x2o_Faxa_rain       ! water flux due to rain           (kg/m2/s)
  integer :: index_x2o_Faxa_bcphidry   ! flux: Black   Carbon hydrophilic dry deposition
  integer :: index_x2o_Faxa_bcphodry   ! flux: Black   Carbon hydrophobic dry deposition
  integer :: index_x2o_Faxa_bcphiwet   ! flux: Black   Carbon hydrophilic wet deposition
  integer :: index_x2o_Faxa_ocphidry   ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_x2o_Faxa_ocphodry   ! flux: Organic Carbon hydrophobic dry deposition
  integer :: index_x2o_Faxa_ocphiwet   ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_x2o_Faxa_dstwet1    ! flux: Size 1 dust -- wet deposition
  integer :: index_x2o_Faxa_dstwet2    ! flux: Size 2 dust -- wet deposition
  integer :: index_x2o_Faxa_dstwet3    ! flux: Size 3 dust -- wet deposition
  integer :: index_x2o_Faxa_dstwet4    ! flux: Size 4 dust -- wet deposition
  integer :: index_x2o_Faxa_dstdry1    ! flux: Size 1 dust -- dry deposition
  integer :: index_x2o_Faxa_dstdry2    ! flux: Size 2 dust -- dry deposition
  integer :: index_x2o_Faxa_dstdry3    ! flux: Size 3 dust -- dry deposition
  integer :: index_x2o_Faxa_dstdry4    ! flux: Size 4 dust -- dry deposition
  integer :: index_x2o_Foxx_rofl       ! river runoff flux                (kg/m2/s)
  integer :: index_x2o_Foxx_rofi       ! ice runoff flux                  (kg/m2/s)

! optional per thickness category fields

!  integer, dimension(:), allocatable :: index_x2o_frac_col       ! fraction of ocean cell, per column
!  integer, dimension(:), allocatable :: index_x2o_fracr_col      ! fraction of ocean cell used in radiation computations, per column
!  integer, dimension(:), allocatable :: index_x2o_qsw_fracr_col  ! qsw * fracr, per column

contains

  subroutine blom_cpl_indices_set( )

    type(mct_aVect) :: o2x      ! temporary
    type(mct_aVect) :: x2o      ! temporary

    integer          :: ncat  ! thickness category index
    character(len=2) :: cncat ! character version of ncat
    integer          :: ncol  ! column index

    ! Determine attribute vector indices

    ! create temporary attribute vectors
    call mct_aVect_init(x2o, rList=seq_flds_x2o_fields, lsize=1)
    call mct_aVect_init(o2x, rList=seq_flds_o2x_fields, lsize=1)

    index_o2x_So_t          = mct_avect_indexra(o2x,'So_t')
    index_o2x_So_u          = mct_avect_indexra(o2x,'So_u')
    index_o2x_So_v          = mct_avect_indexra(o2x,'So_v')
    index_o2x_So_s          = mct_avect_indexra(o2x,'So_s')
    index_o2x_So_dhdx       = mct_avect_indexra(o2x,'So_dhdx')
    index_o2x_So_dhdy       = mct_avect_indexra(o2x,'So_dhdy')
    ! QL, 150526, to wav, boundary layer depth
    index_o2x_So_bldepth    = mct_avect_indexra(o2x,'So_bldepth')
    index_o2x_Fioo_q        = mct_avect_indexra(o2x,'Fioo_q')
    index_o2x_Faoo_fco2_ocn = mct_avect_indexra(o2x,'Faoo_fco2_ocn',perrWith='quiet')
    index_o2x_Faoo_fdms_ocn = mct_avect_indexra(o2x,'Faoo_fdms_ocn',perrWith='quiet')
    index_x2o_Si_ifrac      = mct_avect_indexra(x2o,'Si_ifrac')
    index_x2o_Sa_pslv       = mct_avect_indexra(x2o,'Sa_pslv')
    index_x2o_So_duu10n     = mct_avect_indexra(x2o,'So_duu10n')
    ! QL, 150526, from wav
    index_x2o_Sw_lamult     = mct_avect_indexra(x2o,'Sw_lamult')
    index_x2o_Sw_ustokes    = mct_avect_indexra(x2o,'Sw_ustokes')
    index_x2o_Sw_vstokes    = mct_avect_indexra(x2o,'Sw_vstokes')

    index_x2o_Foxx_tauy     = mct_avect_indexra(x2o,'Foxx_tauy')
    index_x2o_Foxx_taux     = mct_avect_indexra(x2o,'Foxx_taux')
    index_x2o_Foxx_swnet    = mct_avect_indexra(x2o,'Foxx_swnet')
    index_x2o_Foxx_lat      = mct_avect_indexra(x2o,'Foxx_lat')
    index_x2o_Foxx_sen      = mct_avect_indexra(x2o,'Foxx_sen')
    index_x2o_Foxx_lwup     = mct_avect_indexra(x2o,'Foxx_lwup')
    index_x2o_Faxa_lwdn     = mct_avect_indexra(x2o,'Faxa_lwdn')
    index_x2o_Fioi_melth    = mct_avect_indexra(x2o,'Fioi_melth')   
    index_x2o_Fioi_meltw    = mct_avect_indexra(x2o,'Fioi_meltw')
    index_x2o_Fioi_salt     = mct_avect_indexra(x2o,'Fioi_salt')   
    index_x2o_Fioi_bcpho    = mct_avect_indexra(x2o,'Fioi_bcpho')
    index_x2o_Fioi_bcphi    = mct_avect_indexra(x2o,'Fioi_bcphi')
    index_x2o_Fioi_flxdst   = mct_avect_indexra(x2o,'Fioi_flxdst')
    index_x2o_Faxa_prec     = mct_avect_indexra(x2o,'Faxa_prec')   
    index_x2o_Faxa_snow     = mct_avect_indexra(x2o,'Faxa_snow')   
    index_x2o_Faxa_rain     = mct_avect_indexra(x2o,'Faxa_rain')   
    index_x2o_Foxx_evap     = mct_avect_indexra(x2o,'Foxx_evap')
    index_x2o_Foxx_rofl     = mct_avect_indexra(x2o,'Foxx_rofl')
    index_x2o_Foxx_rofi     = mct_avect_indexra(x2o,'Foxx_rofi')
    index_x2o_Faxa_bcphidry = mct_avect_indexra(x2o,'Faxa_bcphidry')
    index_x2o_Faxa_bcphodry = mct_avect_indexra(x2o,'Faxa_bcphodry')
    index_x2o_Faxa_bcphiwet = mct_avect_indexra(x2o,'Faxa_bcphiwet')
    index_x2o_Faxa_ocphidry = mct_avect_indexra(x2o,'Faxa_ocphidry')
    index_x2o_Faxa_ocphodry = mct_avect_indexra(x2o,'Faxa_ocphodry')
    index_x2o_Faxa_ocphiwet = mct_avect_indexra(x2o,'Faxa_ocphiwet')
    index_x2o_Faxa_dstdry1  = mct_avect_indexra(x2o,'Faxa_dstdry1')
    index_x2o_Faxa_dstdry2  = mct_avect_indexra(x2o,'Faxa_dstdry2')
    index_x2o_Faxa_dstdry3  = mct_avect_indexra(x2o,'Faxa_dstdry3')
    index_x2o_Faxa_dstdry4  = mct_avect_indexra(x2o,'Faxa_dstdry4')
    index_x2o_Faxa_dstwet1  = mct_avect_indexra(x2o,'Faxa_dstwet1')
    index_x2o_Faxa_dstwet2  = mct_avect_indexra(x2o,'Faxa_dstwet2')
    index_x2o_Faxa_dstwet3  = mct_avect_indexra(x2o,'Faxa_dstwet3')
    index_x2o_Faxa_dstwet4  = mct_avect_indexra(x2o,'Faxa_dstwet4')
    index_x2o_Sa_co2prog    = mct_avect_indexra(x2o,'Sa_co2prog',perrWith='quiet')
    index_x2o_Sa_co2diag    = mct_avect_indexra(x2o,'Sa_co2diag',perrWith='quiet')
    index_x2o_Faxa_nhx      = mct_avect_indexra(x2o,'Faxa_nhx',perrWith='quiet')
    index_x2o_Faxa_noy      = mct_avect_indexra(x2o,'Faxa_noy',perrWith='quiet')

    ! optional per thickness category fields

    ! convert cpl indices to mcog column indices
    ! this implementation only handles columns due to ice thickness categories

!    lmcog_flds_sent = seq_flds_i2o_per_cat

!    if (seq_flds_i2o_per_cat) then
!      mcog_ncols = ice_ncat+1
!      allocate(index_x2o_frac_col(mcog_ncols))
!      allocate(index_x2o_fracr_col(mcog_ncols))
!      allocate(index_x2o_qsw_fracr_col(mcog_ncols))

!      ncol = 1
!      index_x2o_frac_col(ncol)        = mct_avect_indexra(x2o,'Sf_afrac')
!      index_x2o_fracr_col(ncol)       = mct_avect_indexra(x2o,'Sf_afracr')
!      index_x2o_qsw_fracr_col(ncol)   = mct_avect_indexra(x2o,'Foxx_swnet_afracr')

!      do ncat = 1, ice_ncat
!        write(cncat,'(i2.2)') ncat
!        ncol = ncat+1
!        index_x2o_frac_col(ncol)      = mct_avect_indexra(x2o,'Si_ifrac_'//cncat)
!        index_x2o_fracr_col(ncol)     = index_x2o_frac_col(ncol)
!        index_x2o_qsw_fracr_col(ncol) = mct_avect_indexra(x2o,'PFioi_swpen_ifrac_'//cncat)
!      enddo
!    else
!      mcog_ncols = 1
!    endif

    call mct_aVect_clean(x2o)
    call mct_aVect_clean(o2x)

  end subroutine blom_cpl_indices_set

end module blom_cpl_indices
