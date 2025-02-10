#!/bin/bash
# printf "some data for the file\nAnd a new line" >> fileName
module load NCO/5.1.9-iomkl-2022a
module load CDO/2.0.6-gompi-2022a

basename="ndeposition"
version=$(date '+%Y%m%d')
declare -a TPS=("pi" "HIST" "SSP119" "SSP126" "SSP245" "SSP370" "SSP434" "SSP460" "SSP534os" "SSP585")
declare -a GRIDS=("tnx2v1" "tnx1v4" "tnx0.5v1" "tnx0.25v4" "tnx0.125v4)

#declare -a TPS=("SSP585")
#declare -a GRIDS=("tnx0.5")

for GRID in "${GRIDS[@]}"
do
  echo "$GRID"
  if [ ${GRID} == 'tnx2v1' ]; then
     GRIDFILE="grid_${GRID}_20130206.nc"
  elif [ ${GRID} == 'tnx1v4' ]; then
     GRIDFILE="grid_${GRID}_20170622.nc"
  elif [ ${GRID} == 'tnx0.5v1' ]; then
     GRIDFILE="grid_tnx0.5v1_20240702.nc"
  elif [ ${GRID} == 'tnx0.25v4' ]; then
     GRIDFILE="grid_tnx0.25v4_20170622.nc"
  elif [ ${GRID} == 'tnx0.125v4' ]; then
     GRIDFILE="grid_tnx0.125v4_20221013.nc"
  fi

  #-----------------------------------------------------------------generate gridfile:
  # get the respective variables:
  ncks -A -v plon,plat,parea,pclon,pclat $GRIDFILE gridtmp.nc
  # rename plon and plat
  ncrename -v plon,lon -v plat,lat gridtmp.nc 
  # add lon, lat variables tas coordinates to parea
  ncatted -a coordinates,parea,c,c,"lon lat" gridtmp.nc
  # rename attribute name corners to bounds (cf-conventions)
  ncrename -a lon@corners,bounds -a lat@corners,bounds gridtmp.nc
    
  # re-order coordinates to enable cdo to read/interpret them 
  ncpdq -O --rdr=y,x,nv gridtmp.nc gridfile_${GRID}.nc  
    
  # check the structure of the gridfile
  cdo verifygrid gridfile_${GRID}.nc
    
    
  # alternative: way of generating weights file
  #   cdo griddes gridfile.nc > grid.txt
  #   cdo gencon,grid.txt ndeposition-kg.nc weights.nc

  for TP in "${TPS[@]}"
  do
    echo "  $TP"    
    if [ ${TP} == 'pi' ]; then
      YNAME="1850_CMIP6"
      VCHECK=
      DRYNHx=drynhx_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-185012.nc
      DRYNOy=drynoy_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-185012.nc
      WETNHx=wetnhx_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-185012.nc
      WETNOy=wetnoy_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-185012.nc
      STARTYEAR=1850
      ENDYEAR=1850
    elif [ ${TP} == 'HIST' ]; then
      YNAME="185001-201412"
      # historical time
      DRYNHx=drynhx_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-201412.nc
      DRYNOy=drynoy_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-201412.nc
      WETNHx=wetnhx_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-201412.nc
      WETNOy=wetnoy_input4MIPs_surfaceFluxes_CMIP_NCAR-CCMI-2-0_gn_185001-201412.nc
      STARTYEAR=1850
      ENDYEAR=2014
    elif [ ${TP} == 'SSP119' ]; then
      YNAME="201501-210012-ssp119"
      DRYNHx=drynhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp119-1-0_gn_201501-210012.nc
      DRYNOy=drynoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp119-1-0_gn_201501-210012.nc
      WETNHx=wetnhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp119-1-0_gn_201501-210012.nc
      WETNOy=wetnoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp119-1-0_gn_201501-210012.nc
      STARTYEAR=2015
      ENDYEAR=2100
    elif [ ${TP} == 'SSP126' ]; then
      YNAME="201501-210012-ssp126"
      DRYNHx=drynhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp126-1-0_gn_201501-209912.nc
      DRYNOy=drynoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp126-1-0_gn_201501-209912.nc
      WETNHx=wetnhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp126-1-0_gn_201501-209912.nc
      WETNOy=wetnoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp126-2-0_gn_201501-209912.nc
      STARTYEAR=2015
      ENDYEAR=2100
    elif [ ${TP} == 'SSP245' ]; then
      YNAME="201501-210012-ssp245"
      DRYNHx=drynhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp245-1-0_gn_201501-209912.nc
      DRYNOy=drynoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp245-1-0_gn_201501-209912.nc
      WETNHx=wetnhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp245-1-0_gn_201501-209912.nc
      WETNOy=wetnoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp245-2-0_gn_201501-209912.nc
      STARTYEAR=2015
      ENDYEAR=2100
    elif [ ${TP} == 'SSP370' ]; then
      YNAME="201501-210012-ssp370"
      DRYNHx=drynhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp370-1-0_gn_201501-209912.nc
      DRYNOy=drynoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp370-1-0_gn_201501-209912.nc
      WETNHx=wetnhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp370-1-0_gn_201501-209912.nc
      WETNOy=wetnoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp370-2-0_gn_201501-209912.nc
      STARTYEAR=2015
      ENDYEAR=2100
    elif [ ${TP} == 'SSP434' ]; then
      YNAME="201501-210012-ssp434"
      DRYNHx=drynhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp434-1-0_gn_201501-210012.nc
      DRYNOy=drynoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp434-1-0_gn_201501-210012.nc
      WETNHx=wetnhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp434-1-0_gn_201501-210012.nc
      WETNOy=wetnoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp434-1-0_gn_201501-210012.nc
      STARTYEAR=2015
      ENDYEAR=2100
    elif [ ${TP} == 'SSP460' ]; then
      YNAME="201501-210012-ssp460"
      DRYNHx=drynhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp460-1-0_gn_201501-210012.nc
      DRYNOy=drynoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp460-1-0_gn_201501-210012.nc
      WETNHx=wetnhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp460-1-0_gn_201501-210012.nc
      WETNOy=wetnoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp460-1-0_gn_201501-210012.nc
      STARTYEAR=2015
      ENDYEAR=2100
    elif [ ${TP} == 'SSP534os' ]; then
      YNAME="201501-210012-ssp534os"
      DRYNHx=drynhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp534os-1-0_gn_201501-210012.nc
      DRYNOy=drynoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp534os-1-0_gn_201501-210012.nc
      WETNHx=wetnhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp534os-1-0_gn_201501-210012.nc
      WETNOy=wetnoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp534os-1-0_gn_201501-210012.nc
      STARTYEAR=2015
      ENDYEAR=2100
    elif [ ${TP} == 'SSP585' ]; then
      YNAME="201501-210012-ssp585"
      DRYNHx=drynhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp585-1-0_gn_201501-209912.nc
      DRYNOy=drynoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp585-1-0_gn_201501-209912.nc
      WETNHx=wetnhx_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp585-1-0_gn_201501-209912.nc
      WETNOy=wetnoy_input4MIPs_surfaceFluxes_ScenarioMIP_NCAR-CCMI-ssp585-2-0_gn_201501-209912.nc
      STARTYEAR=2015
      ENDYEAR=2100
    fi
    OUTFILE=${basename}_${YNAME}_${GRID}_${version}.nc
    CHECKTRY=ndep_${YNAME}_${GRID}_
    CHECKFILE=`ls ./${CHECKTRY}*.nc`
    echo "  ${GRIDFILE}"
    echo "  ${DRYNHx}"
    echo "  ${STARTYEAR}"  
    echo "  ${ENDYEAR}"
    echo "  ${OUTFILE}"  
    echo "  CHECK: ${CHECKFILE}"  


    # pi-control input4MIPs files downloaded from esgf: https://esgf-node.llnl.gov/search/input4mips/
    
    
    # dry nhx & noy: kg/m2/s
    # wet nhx & noy: kg/m2/s
    
    # Constant for unit conversion:
    # iHAMOCC units for the input file: kmol N m-2 yr-1
    # Possible options for unit conversions:
    # - considering no leap years, but mol-weight 14.01
    #   s -> y: 365*86400 s/yr
    #   kg -> kmol: 1/14.01
    #   c = 2250963.597
    # - considering leap years (probably would be the most correct one)
    #   s -> yr: 365.25 * 86400
    #   kg -> kmol: 1/14.01
    #   c: 2252505.3533190577
    # - considering leap years, but take mol-weight of 14:
    #   s -> yr: 365.25 * 86400
    #   kg -> kmol: 1/14.
    #   c: 2254114.285714286
    # - not considering leap years and take mol-weight of 14:
    #   s -> yr: 365. * 86400
    #   kg -> kmol: 1/14.
    #   c: 2252571.4285714286
    # The last option is closest to the currently provided deposition file for iHAMOCC
    
    year=365.
    day=86400.
    molweightN=14. # CAM: 14.00674
    
    
    #-----------------------------------------------------------------------------------
    # Combine the deposition files into one file, where NHx and NOy deposition is separated 
    # drynhx: kg/m2/s
    # wetnhx: kg/m2/s
    txt=`cdo -showattribute,drynhx@original_name  $DRYNHx` 
    attdrynhx=${txt##*=}
    txt=`cdo -showattribute,wetnhx@original_name  $WETNHx` 
    attwetnhx=${txt##*=}
    nhxorgname="${attdrynhx//\"} and ${attwetnhx//\"}"
    echo "-------> set attribute for nhxdep: ${nhxorgname}"
    cdo -setattribute,nhxdep@standard_name="tendency_of_atmosphere_mass_content_of_nhx_expressed_as_nitrogen_due_to_dry_and_wet_deposition" \
        -setattribute,nhxdep@long_name="dry_and_wet_nhx_deposition" \
        -setattribute,nhxdep@original_name="${nhxorgname}" \
        -setname,nhxdep -add -selvar,drynhx $DRYNHx -selvar,wetnhx $WETNHx nhxdep-temp.nc
    
    # drynoy: kg/m2/s
    # wetnoy: kg/m2/s
    txt=`cdo -showattribute,drynoy@original_name  $DRYNOy` 
    attdrynoy=${txt##*=}
    txt=`cdo -showattribute,wetnoy@original_name  $WETNOy` 
    attwetnoy=${txt##*=}
    noyorgname="${attdrynoy//\"} and ${attwetnoy//\"}"
    echo "-------> set attribute for noydep: ${noyorgname}"
    cdo -setattribute,noydep@standard_name="tendency_of_atmosphere_mass_content_of_noy_expressed_as_nitrogen_due_to_dry_and_wet_deposition" \
        -setattribute,noydep@long_name="dry_and_wet_noy_deposition" \
        -setattribute,noydep@original_name="${noyorgname}" \
        -setname,noydep -add -selvar,drynoy $DRYNOy -selvar,wetnoy $WETNOy noydep-temp.nc
    
    txt=`cdo -showattribute,nhxdep@original_name nhxdep-temp.nc`
    attnhx=${txt##*=}
    txt=`cdo -showattribute,noydep@original_name noydep-temp.nc`
    attnoy=${txt##*=}
    ndeporgname="sum of ${attnhx//\"} and ${attnoy//\"}"
    echo "-------> set attribute for ndep: ${ndeporgname}"
    cdo -setattribute,ndep@standard_name="tendency_of_atmosphere_mass_content_of_nitrogen_expressed_as_nitrogen_due_to_dry_and_wet_deposition" \
        -setattribute,ndep@long_name="dry_and_wet_sum_noy_and_nhx_deposition" \
        -setattribute,ndep@original_name="${ndeporgname}" \
        -setname,ndep -add -selvar,noydep  noydep-temp.nc -selvar,nhxdep nhxdep-temp.nc ndep-temp.nc
    
    # merge the nitrogen deposition files into one:
    cdo merge nhxdep-temp.nc noydep-temp.nc ndep-temp.nc ndeposition-kg_${GRID}.nc

    if [ ! -f  weights_${GRID}.nc ]; then  
      # generate the weights for conservative remapping
      cdo gencon,gridfile_${GRID}.nc ndeposition-kg_${GRID}.nc weights_${GRID}.nc
    fi

 
    #----------------------------------------------- Conservative remapping to BLOM grid and unit conversio:
    c=`python -c "print(${year}*${day}/${molweightN})"`
    echo "----------> Unit conversion factor (kg m-2 s-1 to kmol m-2 yr-1): $c"
    
    if [[ $DRYNHx =~ "2099" ]];then 
      # setting some few missing values to zero (land areas anyway) and perform the unit conversion
      cdo -setattribute,ndep@units="kmol N m-2 yr-1"   \
          -setattribute,noydep@units="kmol N m-2 yr-1" \
          -setattribute,nhxdep@units="kmol N m-2 yr-1" \
          -mulc,$c -setmisstoc,0.                      \
          -remap,gridfile_${GRID}.nc,weights_${GRID}.nc  ndeposition-kg_${GRID}.nc tmpfull.nc
      # alternative: direct conversion without writing the weights file
      # -remapcon,gridfile_${GRID}.nc  ndeposition-kg.nc ndeposition_${GRID}.nc
   
      cdo -selyear,"2099" tmpfull.nc tmp2099.nc
      cdo -setyear,"2100" tmp2099.nc tmp2100.nc
      rm tmp2099.nc
      cdo mergetime tmpfull.nc tmp2100.nc ${OUTFILE}
      rm tmp2100.nc
      rm tmpfull.nc
    else 
      # setting some few missing values to zero (land areas anyway) and perform the unit conversion
      cdo -setattribute,ndep@units="kmol N m-2 yr-1"   \
          -setattribute,noydep@units="kmol N m-2 yr-1" \
          -setattribute,nhxdep@units="kmol N m-2 yr-1" \
          -mulc,$c -setmisstoc,0.                      \
          -remap,gridfile_${GRID}.nc,weights_${GRID}.nc  ndeposition-kg_${GRID}.nc ${OUTFILE}
    fi

 
    # add global attributes start- and endyear to the file (required by iHAMOCC) 
    ncatted --glb_att_add startyear=${STARTYEAR} ${OUTFILE} 
    ncatted --glb_att_add endyear=${ENDYEAR} ${OUTFILE}
    # modify from char to integer values
    ncatted -O -h -a startyear,global,m,i,${STARTYEAR} ${OUTFILE} 
    ncatted -O -h -a endyear,global,m,i,${ENDYEAR} ${OUTFILE} 
    
    if [ ${TP} == 'HIST' ]; then
      cdo -selyear,2000 ${OUTFILE} ${basename}_2000_CMIP6_${GRID}_${version}.nc
    fi 
    
    if [ -n "${CHECKFILE}" ]; then
      #----------------------------------------------- Check wrt old input file
      echo "----------> Check generated N-deposition files wrt former N-deposition"
      cdo -setattribute,err@units="unitless" \
          -setattribute,err@long_name="error_ratio_ndep_new_div_old_minus_one" \
          -setname,err -subc,1 -div -add -selvar,noydep ${OUTFILE} -selvar,nhxdep ${OUTFILE}  ${CHECKFILE} ${basename}_${YNAME}_${GRID}_${version}_check_ratio_error.nc
    
      cdo -setattribute,perc_err@units="%" \
          -setattribute,perc_err@long_name="percent_error_wet_and_dry_depostion_wrt_former_ndep" \
          -setname,perc_err -mulc,100 -div -sub -add -selvar,noydep ${OUTFILE} -selvar,nhxdep ${OUTFILE}  ${CHECKFILE}  ${CHECKFILE} ${basename}_${YNAME}_${GRID}_${version}_check_percent_error.nc
    fi

    # Clean-up
    rm nhxdep-temp.nc
    rm noydep-temp.nc
    rm ndep-temp.nc
    rm ndeposition-kg_${GRID}.nc
    CHECKFIL=""
  done
  rm gridtmp.nc
  rm gridfile_${GRID}.nc
  rm weights_${GRID}.nc
done
exit 

