module micom_cpl_indices
  
  use seq_flds_mod
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
  integer :: index_o2x_Fioo_q
  integer :: index_o2x_Faoo_fco2_ocn
  integer :: index_o2x_Faoo_fdms_ocn

  ! drv -> ocn

  integer :: index_x2o_Si_ifrac        ! fractional ice wrt ocean
  integer :: index_x2o_So_duu10n       ! 10m wind speed squared           (m^2/s^2)
  integer :: index_x2o_Sa_pslv         ! sea-level pressure               (Pa)
  integer :: index_x2o_Sa_co2prog      ! bottom atm level prognostic CO2
  integer :: index_x2o_Sa_co2diag      ! bottom atm level diagnostic CO2
  integer :: index_x2o_Foxx_taux       ! zonal wind stress (taux)         (W/m2   )
  integer :: index_x2o_Foxx_tauy       ! meridonal wind stress (tauy)     (W/m2   )
  integer :: index_x2o_Foxx_swnet      ! net short-wave heat flux         (W/m2   )
  integer :: index_x2o_Foxx_sen        ! sensible heat flux               (W/m2   )
  integer :: index_x2o_Foxx_lat        
  integer :: index_x2o_Foxx_lwdn       ! longwave radiation (down)        (W/m2   )
  integer :: index_x2o_Foxx_lwup       ! longwave radiation (up)          (W/m2   )
  integer :: index_x2o_Foxx_melth      ! heat flux from snow & ice melt   (W/m2   )
  integer :: index_x2o_Foxx_salt       ! salt                             (kg(salt)/m2/s)
  integer :: index_x2o_Foxx_prec         
  integer :: index_x2o_Foxx_snow       ! water flux due to snow           (kg/m2/s)
  integer :: index_x2o_Foxx_rain       ! water flux due to rain           (kg/m2/s)
  integer :: index_x2o_Foxx_evap       ! evaporation flux                 (kg/m2/s)
  integer :: index_x2o_Foxx_meltw      ! snow melt flux                   (kg/m2/s)
  integer :: index_x2o_Forr_roff       ! river runoff flux                (kg/m2/s)
  integer :: index_x2o_Forr_ioff       ! ice runoff flux                  (kg/m2/s)
  integer :: index_x2o_Foxx_bcphidry   ! flux: Black   Carbon hydrophilic dry deposition
  integer :: index_x2o_Foxx_bcphodry   ! flux: Black   Carbon hydrophobic dry deposition
  integer :: index_x2o_Foxx_bcphiwet   ! flux: Black   Carbon hydrophilic wet deposition
  integer :: index_x2o_Foxx_ocphidry   ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_x2o_Foxx_ocphodry   ! flux: Organic Carbon hydrophobic dry deposition
  integer :: index_x2o_Foxx_ocphiwet   ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_x2o_Foxx_dstwet1    ! flux: Size 1 dust -- wet deposition
  integer :: index_x2o_Foxx_dstwet2    ! flux: Size 2 dust -- wet deposition
  integer :: index_x2o_Foxx_dstwet3    ! flux: Size 3 dust -- wet deposition
  integer :: index_x2o_Foxx_dstwet4    ! flux: Size 4 dust -- wet deposition
  integer :: index_x2o_Foxx_dstdry1    ! flux: Size 1 dust -- dry deposition
  integer :: index_x2o_Foxx_dstdry2    ! flux: Size 2 dust -- dry deposition
  integer :: index_x2o_Foxx_dstdry3    ! flux: Size 3 dust -- dry deposition
  integer :: index_x2o_Foxx_dstdry4    ! flux: Size 4 dust -- dry deposition

contains

  subroutine micom_cpl_indices_set( )

    type(mct_aVect) :: o2x      ! temporary
    type(mct_aVect) :: x2o      ! temporary

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
    index_o2x_Fioo_q        = mct_avect_indexra(o2x,'Fioo_q')

    !index_o2x_Faoo_fco2_ocn= mct_avect_indexra(o2x,'Faoo_fco2_ocn',perrWith='quiet')
    index_o2x_Faoo_fco2_ocn = mct_avect_indexra(o2x,'Faoo_fco2'    ,perrWith='quiet')

    index_o2x_Faoo_fdms_ocn = mct_avect_indexra(o2x,'Faoo_fdms_ocn',perrWith='quiet')

    index_x2o_Si_ifrac      = mct_avect_indexra(x2o,'Si_ifrac')
    index_x2o_Sa_pslv       = mct_avect_indexra(x2o,'Sa_pslv')

    !index_x2o_So_duu10n    = mct_avect_indexra(x2o,'So_duu10n')
    index_x2o_So_duu10n     = mct_avect_indexra(x2o,'Sx_duu10n')

    index_x2o_Foxx_tauy     = mct_avect_indexra(x2o,'Foxx_tauy')
    index_x2o_Foxx_taux     = mct_avect_indexra(x2o,'Foxx_taux')
    index_x2o_Foxx_swnet    = mct_avect_indexra(x2o,'Foxx_swnet')
    index_x2o_Foxx_lat      = mct_avect_indexra(x2o,'Foxx_lat')
    index_x2o_Foxx_sen      = mct_avect_indexra(x2o,'Foxx_sen')
    index_x2o_Foxx_lwdn     = mct_avect_indexra(x2o,'Foxx_lwdn')
    index_x2o_Foxx_lwup     = mct_avect_indexra(x2o,'Foxx_lwup')
    index_x2o_Foxx_melth    = mct_avect_indexra(x2o,'Foxx_melth')   
    index_x2o_Foxx_salt     = mct_avect_indexra(x2o,'Foxx_salt')   
    index_x2o_Foxx_prec     = mct_avect_indexra(x2o,'Foxx_prec')   
    index_x2o_Foxx_snow     = mct_avect_indexra(x2o,'Foxx_snow')   
    index_x2o_Foxx_rain     = mct_avect_indexra(x2o,'Foxx_rain')   
    index_x2o_Foxx_evap     = mct_avect_indexra(x2o,'Foxx_evap')
    index_x2o_Foxx_meltw    = mct_avect_indexra(x2o,'Foxx_meltw')
    index_x2o_Forr_roff     = mct_avect_indexra(x2o,'Forr_roff')
    index_x2o_Forr_ioff     = mct_avect_indexra(x2o,'Forr_ioff')
    index_x2o_Foxx_bcphidry = mct_avect_indexra(x2o,'Foxx_bcphidry')
    index_x2o_Foxx_bcphodry = mct_avect_indexra(x2o,'Foxx_bcphodry')
    index_x2o_Foxx_bcphiwet = mct_avect_indexra(x2o,'Foxx_bcphiwet')
    index_x2o_Foxx_ocphidry = mct_avect_indexra(x2o,'Foxx_ocphidry')
    index_x2o_Foxx_ocphodry = mct_avect_indexra(x2o,'Foxx_ocphodry')
    index_x2o_Foxx_ocphiwet = mct_avect_indexra(x2o,'Foxx_ocphiwet')
    index_x2o_Foxx_dstdry1  = mct_avect_indexra(x2o,'Foxx_dstdry1')
    index_x2o_Foxx_dstdry2  = mct_avect_indexra(x2o,'Foxx_dstdry2')
    index_x2o_Foxx_dstdry3  = mct_avect_indexra(x2o,'Foxx_dstdry3')
    index_x2o_Foxx_dstdry4  = mct_avect_indexra(x2o,'Foxx_dstdry4')
    index_x2o_Foxx_dstwet1  = mct_avect_indexra(x2o,'Foxx_dstwet1')
    index_x2o_Foxx_dstwet2  = mct_avect_indexra(x2o,'Foxx_dstwet2')
    index_x2o_Foxx_dstwet3  = mct_avect_indexra(x2o,'Foxx_dstwet3')
    index_x2o_Foxx_dstwet4  = mct_avect_indexra(x2o,'Foxx_dstwet4')
    index_x2o_Sa_co2prog    = mct_avect_indexra(x2o,'Sa_co2prog',perrWith='quiet')
    index_x2o_Sa_co2diag    = mct_avect_indexra(x2o,'Sa_co2diag',perrWith='quiet')

    call mct_aVect_clean(x2o)
    call mct_aVect_clean(o2x)

  end subroutine micom_cpl_indices_set

end module micom_cpl_indices
